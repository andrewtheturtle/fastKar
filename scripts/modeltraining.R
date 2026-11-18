#
# modeltraining.R: train every parameter of the fastKar Hi-C forward model from GM12878 (Rao 2014, 4DN .hic, hg38, chr1-22 + X).
#   distance spline (point density, for the pixel-integrated kernel)  -> data/scalefree_spline.rda (splineobj$distance_spline)
#   diagonal fit (log diagonal density vs log bin area)               -> data/scalefree_spline.rda (splineobj$diag_spline)
#   interchrom_density D (genome-average trans density per copy pair)  -> data/interchrom_density.rda
#   depth_calibration (off-diagonal contacts per kept copy-bp per unit depth, for estimate.depthratio()) -> data/depth_calibration.rda
#   kinf model JSON (inst/python/kinf/data/fastkar_model.json)        <- scripts/export_kinf_model.R
#   theta (NB dispersion; fit and reported; the package uses theta = 2, user decision 2026-11-13)
# Training tables and fits are written to outdir; ~/projects/karyotyping/scripts/model_figure.R makes the Extended Data figure from them.
# Kept-bin recipe (user decisions 2026-11-13; originally the iraj_blog.R sections of 2026-11-13 "retrain the Hi-C model ... on kept bins only"
# and the ones after it):
#   mask at 10 kb: a bin is dropped if >= half covered by the hg38 mask (reference mask only; the low-coverage filter was dropped 2026-11-17);
#   coarse bins keep only their kept 10 kb sub-bins: counts from kept x kept 10 kb pixels, area = kept bp_I x kept bp_J;
#   depth (calibration, 2026-11-17) = valid-pair coverage of the kept bins: read ends in kept bins (row sums + diagonal again) x read length
#   / kept bp; read length 100.78 bp per end (SRA bases / spots / 2 over the 27 source experiments of the .hic); depth = 2 x haploid;
#   D = kept trans counts / (kept trans area x depth x 2) (= / haploid x 4 copy pairs);
#   distance spline from the genome-wide pooled 10 kb densities: target = density / depth - 2D (the kernel adds the trans background to
#   every cis pixel), up to 1e8, then 1/d to 3e8 (tail knots at multiples of 10 kb, needed by kinf's knot lookup table), deconvolved over
#   10 kb pixels to a point density (forward_pixel_cpp);
#   diagonal: mean diagonal count of fully kept coarse bins per resolution, linear fit in log area;
#   theta: direct NB ML (MASS::theta.ml diverges here) on kept pixels (both bins >= half kept), diagonal excluded, pixel integration and
#   balanced trans on; 1 Mb all pixels, 100 kb a seeded 1% sample.
# Run from anywhere (~30 min, 6 cores): Rscript ~/git/fastKar/scripts/modeltraining.R
#
library(data.table)
pkg = '~/git/fastKar'
suppressMessages(devtools::load_all(pkg, quiet = TRUE))
hic.file = '/gpfs/commons/groups/imielinski_lab/data/PoreC/Rao2014/4DNFI1UEG1HD.hic'
mask.file = '/gpfs/commons/groups/imielinski_lab/home/sclarke/git/scFlow/Flow/files/masks/hg38/hg38_mask_zc_stash_nochr.rds'
outdir = '~/projects/karyotyping/data/model_training_refmask/'; dir.create(outdir, showWarnings = FALSE)
read_length = 1090641403300 / 5411108475 / 2  # bp per read end: SRA bases / spots (pairs), 33 runs of the 27 source experiments
fine = 1e4
resolutions = c(1e4,2e4,5e4,1e5,2e5,5e5,1e6,2e6,5e6)
maxval = 1e8  # spline knots up to here, 1/d beyond (as the kernel's maxval)
ncore = 6
chs = c(1:22,'X'); cl = strawr::readHicChroms(hic.file); chrlen = setNames(as.numeric(cl$length[match(chs,cl$name)]), chs)
pairs = CJ(a = seq_along(chs), b = seq_along(chs))[a <= b]

# ---- 10 kb bins and mask
bins = gr2dt(gr.tile(si2gr(chrlen), fine))[, .(seqnames = as.character(seqnames), start, end, width)]
bins[, k := seq_len(.N), by = seqnames][, gid := .I]
off = bins[, .(o = min(gid) - 1L), by = seqnames]; off = setNames(off$o, off$seqnames)
mask38 = reduce(readRDS(mask.file))
gr = dt2gr(bins); ov = findOverlaps(gr, mask38)
mbp = tapply(width(pintersect(gr[queryHits(ov)], mask38[subjectHits(ov)])), queryHits(ov), sum)
bins[, masked := 0][as.integer(names(mbp)), masked := as.numeric(mbp)]

# pass 1: raw total coverage per 10 kb bin (cis + trans; diagonal counted once)
covs = mclapply(seq_len(nrow(pairs)), function(r){
	a = chs[pairs$a[r]]; b = chs[pairs$b[r]]
	x = as.data.table(strawr::straw('NONE', hic.file, a, b, 'BP', fine))
	x[, `:=`(i = off[a] + as.integer(x/fine) + 1L, j = off[b] + as.integer(y/fine) + 1L)]
	rbind(x[, .(c = sum(counts), d = sum(counts[i == j])), by = .(gid = i)], x[i != j, .(c = sum(counts), d = 0), by = .(gid = j)])
}, mc.cores = ncore, mc.preschedule = FALSE)
covt = rbindlist(covs)[, .(c = sum(c), d = sum(d)), by = gid]; rm(covs)
bins[, `:=`(cov = 0, dcount = 0)][covt$gid, `:=`(cov = covt$c, dcount = covt$d)]
bins[, keep := masked < width/2]
bins[, kbp := keep * width]
cat('10 kb bins kept:', sum(bins$keep), '/', nrow(bins), ' kept bp fraction:', sum(bins$kbp)/sum(bins$width), '\n')
# depth calibration: valid-pair coverage of the kept bins (read ends: row sum, diagonal counted once there, + the diagonal again)
depth = bins[keep == TRUE, sum(cov + dcount)] * read_length / sum(bins$kbp)
cat('depth (valid-pair coverage of kept bins)', depth, ' old formula (sum x 300 / 3e9) was 387.01\n')

# autocorrelation sum_I a_I a_{I+k}, k = 0..n-1 (kept area per offset)
acf_sum = function(a){ n = length(a); p = 2^ceiling(log2(2*n)); f = fft(c(a, rep(0, p - n)))
	pmax(Re(fft(Mod(f)^2, inverse = TRUE))[1:n]/p, 0) }

# pass 2: cis, per chromosome and resolution: kept counts and kept area per offset k (dist = k x res), diagonal per coarse bin
# (the chromosome is `ch`, not `a`: inside ak[...] a name `a` would resolve to its kept-bp column)
cis = mclapply(chs, function(ch){
	x = as.data.table(strawr::straw('NONE', hic.file, ch, ch, 'BP', fine))
	x[, `:=`(i = as.integer(x/fine) + 1L, j = as.integer(y/fine) + 1L)]
	bb = bins[seqnames == ch]
	x = x[bb$keep[i] & bb$keep[j]]
	res = lapply(resolutions, function(r){
		f = as.integer(r/fine)
		I = (x$i - 1L) %/% f; J = (x$j - 1L) %/% f
		cnt = data.table(k = J - I, c = x$counts)[, .(counts = sum(c)), by = k]
		ak = bb[, .(a = sum(kbp), full = sum(width)), by = .(I = (k - 1L) %/% f)]
		ar = data.table(k = seq_len(nrow(ak)) - 1L, area = acf_sum(ak$a))
		o = merge(ar, cnt, by = 'k', all.x = TRUE)[is.na(counts), counts := 0][, `:=`(res = r, chr = ch)]
		dg = data.table(I = I[I == J], c = x$counts[I == J])[, .(c = sum(c)), by = I]
		ak = merge(ak, dg, by = 'I', all.x = TRUE)[is.na(c), c := 0][, `:=`(res = r, chr = ch)]
		list(off = o, diag = ak)
	})
	list(off = rbindlist(lapply(res, `[[`, 'off')), diag = rbindlist(lapply(res, `[[`, 'diag')))
}, mc.cores = ncore, mc.preschedule = FALSE)
offk = rbindlist(lapply(cis, `[[`, 'off')); diagb = rbindlist(lapply(cis, `[[`, 'diag')); rm(cis)

# pass 3: trans, kept x kept 10 kb pixels; area = kept bp_a x kept bp_b
trans = rbindlist(mclapply(which(pairs$a < pairs$b), function(r){
	a = chs[pairs$a[r]]; b = chs[pairs$b[r]]
	x = as.data.table(strawr::straw('NONE', hic.file, a, b, 'BP', fine))
	ka = bins[seqnames == a]$keep; kb = bins[seqnames == b]$keep
	data.table(a = a, b = b, counts = sum(x$counts[ka[as.integer(x$x/fine) + 1L] & kb[as.integer(x$y/fine) + 1L]]), counts_all = sum(x$counts),
		area = sum(bins[seqnames == a]$kbp) * sum(bins[seqnames == b]$kbp))
}, mc.cores = ncore, mc.preschedule = FALSE))

# pooled genome-wide distance table (dist > 0), all resolutions; value = mean counts per full res x res pixel, dens = counts / kept area
distk = offk[k > 0, .(counts = sum(counts), area = sum(area)), by = .(res, k)][area > 0][, `:=`(dist = k * res, dens = counts/area)][, value := dens * res^2]
# diagonal table: fully kept coarse bins (kept bp = full width = res), mean count vs area res^2
diagt = diagb[a == full & full == res, .(value = mean(c), n = .N), by = res][, widthprod := res^2]
saveRDS(list(bins = bins, offk = offk, diagb = diagb, trans = trans, distk = distk, diagt = diagt, depth = depth), paste0(outdir, 'kept_training_tables.rds'))

# ---- parameters
# D: kept trans counts per kept area, per copy pair (x depth = 2 x haploid, x 2 more = 4 diploid copy pairs)
interchrom_density = sum(trans$counts) / (sum(trans$area) * depth * 2)

# depth calibration for estimate.depthratio(): off-diagonal kept x kept 10 kb contacts (10 kb self-pixels excluded) per kept copy-bp
# (2 copies) per unit depth, so estimate.depthratio() returns `depth` on these data
depth_calibration = (sum(offk[res == fine & k > 0]$counts) + sum(trans$counts)) / (depth * 2 * sum(bins$kbp))

# diagonal: log(diagonal count / depth) linear in log(area)
dfit = lm(log(value/depth) ~ log(widthprod), data = diagt)
spline_diag = eval(substitute(function(x_new){ a + b * x_new }, list(a = unname(coef(dfit)[1]), b = unname(coef(dfit)[2]))))

# distance: point density rho(s) such that the 10 kb x 10 kb pixel average reproduces the target density / depth - 2D at offsets k x 10 kb.
# Damped (half-step) fixed point: the k = 1 pixel is dominated by the power-law extension below 10 kb, whose level and slope follow the
# first knots; full steps oscillate.
sd10 = distk[res == fine & counts > 0 & dist <= maxval][order(dist)]
target = sd10$dens/depth - 2*interchrom_density
stopifnot(all(target > 0))
xs = log(sd10$dist)
xt = log(unique(round(10^seq(log10(maxval), log10(3e8), length.out = 30)[-1] / fine) * fine))  # 1/d tail knots on the 10 kb grid
mk = function(y) splinefun(c(xs, xt), c(y, tail(y, 1) - (xt - log(maxval))))
y = log(target)
for (it in 1:100){
	f = mk(y); z = environment(f)$z; lmin = z$x[1]
	pav = forward_pixel_cpp(rep(fine, length(xs)), rep(fine, length(xs)), sd10$dist, z, lmin, f(lmin, deriv = 1), maxval, exp(f(log(maxval)))) / fine^2
	err = log(target/pav); y = y + 0.5*err
	if (max(abs(err)) < 1e-6) break
}
cat('deconvolution: iterations', it, ' max |log(target / pixel average)|', max(abs(err)), '\n')
spline_dist = mk(y)

splineobj = list(diag_spline = spline_diag, distance_spline = spline_dist)
cat('depth', depth, ' D', interchrom_density, ' diagonal slope', coef(dfit)[2], ' depth calibration', depth_calibration, '\n')
saveRDS(list(splineobj = splineobj, interchrom_density = interchrom_density, depth = depth, depth_calibration = depth_calibration,
	spline_note = 'point density rho(s) (pixel-integrated kernel): 10 kb kept-bin density / depth - 2D deconvolved over 10 kb pixels, up to 1e8; 1/d to 3e8 (tail knots at multiples of 10 kb)'),
	paste0(outdir, 'kept_model.rds'))

# ---- install as package data, export the kinf model JSON, reload
save(splineobj, file = file.path(pkg, 'data/scalefree_spline.rda'))
save(interchrom_density, file = file.path(pkg, 'data/interchrom_density.rda'))
save(depth_calibration, file = file.path(pkg, 'data/depth_calibration.rda'))
local({ owd = setwd(path.expand(pkg)); on.exit(setwd(owd)); source('scripts/export_kinf_model.R') })
suppressMessages(devtools::load_all(pkg, quiet = TRUE))

# ---- theta: NB ML on GM12878 pixels under the installed model (pixel integration always on; balanced trans on)
# expected cis = simulate_diploid_map(depth = calibrated depth / 2) (haploid depth, as forward_simulate's diploid component);
# trans = depth/2 x D x w_i w_j x B_i B_j, B = reference_trans_weights(); both x kept fractions f_i f_j. Observed = kept x kept 10 kb counts.
options(fastKar.conserve = TRUE)
cm = conserve_model(); cat('conserve K', cm$K, ' S_ref', cm$S_ref, '\n')
theta_px = function(rr, frac){
	f = as.integer(rr/fine)
	cb = bins[, .(a = sum(kbp), w = sum(width), start = min(start), end = max(end)), by = .(seqnames, I = (k - 1L) %/% f)][, frac := a/w]
	out = mclapply(seq_len(nrow(pairs)), function(r){
		set.seed(r)
		A = chs[pairs$a[r]]; B = chs[pairs$b[r]]
		x = as.data.table(strawr::straw('NONE', hic.file, A, B, 'BP', fine))
		ka = bins[seqnames == A]$keep; kb = bins[seqnames == B]$keep
		x = x[ka[as.integer(x/fine) + 1L] & kb[as.integer(y/fine) + 1L]]
		obs = x[, .(obs = sum(counts)), by = .(i = as.integer(x/fine) %/% f + 1L, j = as.integer(y/fine) %/% f + 1L)]
		ta = cb[seqnames == A]; tb = cb[seqnames == B]
		if (A == B){
			sim = simulate_diploid_map(ta[, .(seqnames, start, end, width = w)], depth = depth/2, gm.out = FALSE)$dat
			m = sim[i < j, .(i, j, mu = value * ta$frac[i] * ta$frac[j])]
		} else {
			Ba = reference_trans_weights(ta[, .(seqnames, start, end, width = w)]); Bb = reference_trans_weights(tb[, .(seqnames, start, end, width = w)])
			m = CJ(i = seq_len(nrow(ta)), j = seq_len(nrow(tb)))[, mu := depth/2 * interchrom_density * ta$w[i] * tb$w[j] * Ba[i] * Bb[j] * ta$frac[i] * tb$frac[j]]
		}
		m = m[ta$frac[i] >= 0.5 & tb$frac[j] >= 0.5]
		m = merge(m, obs, by = c('i', 'j'), all.x = TRUE)[is.na(obs), obs := 0]
		m = m[sample(.N, ceiling(frac * .N))]
		m[, `:=`(cis = A == B, k = if (A == B) j - i else NA_integer_, chrA = A, chrB = B)]
	}, mc.cores = ncore, mc.preschedule = FALSE)
	err = sapply(out, inherits, 'try-error'); if (any(err)) stop(out[[which(err)[1]]])
	rbindlist(out)[, res := rr]
}
nb_theta = function(y, mu) exp(optimize(function(lt) -sum(dnbinom(y, size = exp(lt), mu = mu, log = TRUE)), c(-7, 7))$minimum)
px = rbind(theta_px(1e6, 1), theta_px(1e5, 0.01))
fits = px[, .(theta_cis = nb_theta(obs[cis], mu[cis]), theta_trans = nb_theta(obs[!cis], mu[!cis]), theta_all = nb_theta(obs, mu),
	n_cis = sum(cis), n_trans = sum(!cis), oe_cis = sum(obs[cis])/sum(mu[cis]), oe_trans = sum(obs[!cis])/sum(mu[!cis])), by = res]
saveRDS(list(px = px, fits = fits, K = cm$K, S_ref = cm$S_ref), paste0(outdir, 'theta_kept_conserve.rds'))
print(fits)
