# Contact-conserving ("balanced") trans term of the Hi-C forward model.
#
# Principle: every bp of every copy has the same total number of contacts K. Its cis contacts are the distance kernel
# integrated along its own walk, F(x) + F(L - x) at position x of a linear walk of length L (2 F(L/2) on a circle),
# F(d) = integral of distance_density from 10 kb to d. Trans takes the rest: per copy-bp b = (K - cis) / (D S), and the
# trans pixel is width_i width_j D B_i B_j, B_i = the sum of b over the copies on tile i (replaces CN_i in the current
# term width_i width_j D CN_i CN_j). Since trans is rank one, S has a closed form: S^2 = sum over all copies (K - cis) w / D.
# S runs over the whole genome: the walks' copies, each node's remaining graph copies (node CN - copies the walks use,
# placed at their reference positions) and, outside the graph's nodes, 2 reference copies of chr1-22 + X.
# K is fixed (diploid calibration on hg38 chr1-22 + X): K = D 2G + mean reference cis, so a diploid reference genome has
# bp-weighted mean b = 1 (its trans equals the current model's) and S_ref = 2G. A fixed K makes the genome-wide total
# contacts linear in ploidy, as estimate.depthratio() assumes.
# On with options(fastKar.conserve = TRUE) (or conserve = TRUE in forward_simulate(), simulate_walks(), simulate_diploid_map()).
# Derivation and checks: ~/contexts/fastKar/docs/conservation_model_handoff.md.

#' hg38 lengths of chr1-22 and X (from the GM12878 .hic the kernel was trained on): calibration genome of K
conserve_chrlen_hg38 = c(`1` = 248956422, `2` = 242193529, `3` = 198295559, `4` = 190214555, `5` = 181538259, `6` = 170805979,
                         `7` = 159345973, `8` = 145138636, `9` = 138394717, `10` = 133797422, `11` = 135086622, `12` = 133275309,
                         `13` = 114364328, `14` = 107043718, `15` = 101991189, `16` = 90338345, `17` = 83257441, `18` = 80373285,
                         `19` = 58617616, `20` = 64444167, `21` = 46709983, `22` = 50818468, `X` = 156040895)

.conserve_cache = new.env(parent = emptyenv())

#' Constants of the conserved trans term: D, F, G (= integral of F), K, S_ref. Built once per interchrom_density.
#'
#' F and G are the Hi-C kernel's own cumulative integrals (forward_FG_cpp(), the G of the pixel-integrated kernel), restarted
#' at 10 kb: contacts below 10 kb are left out of F, since they are about the same for every bp, so K absorbs them.
#' @return list(D, z, F, G, K, S_ref, chrlen)
conserve_model = function(){
    D = fastKar::interchrom_density
    z = environment(fastKar::splineobj$distance_spline)$z  # cache keyed on the spline too (K and S_ref depend on it)
    m = .conserve_cache$model
    if (!is.null(m) && identical(m$D, D) && identical(m$z, z)) return(m)
    ka = forward_kernel_args(0); x0 = 1e4
    FG = function(s) forward_FG_cpp(as.numeric(s), ka$z, ka$lmin, ka$slope, ka$maxval, ka$maxval_density)
    FG0 = FG(x0)
    m = list(D = D, z = z,
             F = function(d) FG(pmax(d, x0))$F - FG0$F,                                     # int_{10 kb}^{d} rho
             G = function(x){ x = pmax(x, x0); FG(x)$G - FG0$G - FG0$F * (x - x0) },         # int_{10 kb}^{x} F
             chrlen = conserve_chrlen_hg38)
    L = m$chrlen; G2 = 2 * sum(L); refcis = .conserve_refint(m, 0, L, L)
    m$K = D * G2 + 2 * sum(refcis) / G2
    m$S_ref = sqrt(2 * sum(m$K * L - refcis) / D)
    .conserve_cache$model = m
    m
}

#' Integral of the reference cis contacts per bp over [a, b] (0-based continuous coordinates) of a chromosome of length L
.conserve_refint = function(m, a, b, L) m$G(b) - m$G(a) + m$G(L - a) - m$G(L - b)

#' Chromosome lengths for reference positions: given seqlengths where known, else the calibration (hg38) length of the
#' same chromosome ('chr' prefix ignored), else the largest end seen (ends)
.conserve_lengths = function(chr, seqlengths = NULL, ends = NULL){
    chr = as.character(chr)
    L = if (is.null(seqlengths)) rep(NA_real_, length(chr)) else as.numeric(seqlengths[chr])
    miss = is.na(L)
    L[miss] = conserve_chrlen_hg38[sub('^chr', '', chr[miss])]
    miss = is.na(L)
    if (any(miss) && !is.null(ends)) L[miss] = stats::ave(as.numeric(ends), chr, FUN = max)[miss]
    L
}

#' Per-tile trans weights B of a set of walks under the conserved trans term
#'
#' @param walks gWalk or list(graph, snode.id, circular) as in simulate_walks(); every walk counts (also those off the target)
#' @param tiled.target,widthdt from prep_for_sim()
#' @return list(B = per-tile weight, replacing tiled.target$cn in the trans term; S; outside_bp = genome bp outside the graph's nodes)
trans_weights = function(walks, tiled.target, widthdt){
    if (is.null(walks$graph)) stop('conserve = TRUE needs the genome graph: pass a gWalk or list(graph, snode.id, circular)')
    m = conserve_model(); K = m$K; D = m$D
    nd = data.table::as.data.table(walks$graph$nodes$dt)[, .(node.id, chr = as.character(seqnames), start, end, cn = as.numeric(cn))]
    nd[is.na(cn), cn := 0]
    sl = GenomeInfoDb::seqlengths(walks$graph$gr)
    nd[, L := .conserve_lengths(chr, sl, end)]
    sn = walks$snode.id
    wcn = walks$dt$cn; if (is.null(wcn)) wcn = rep(1, length(sn))
    circ = walks$circular
    # every row (tile or untiled node part) of every walk, its position along the walk and its cis integral I
    rows = walk_rows(sn, tiled.target, widthdt)
    r = data.table(w = rep(seq_along(sn), diff(rows$walk_ptr)), orig = rows$orig, width = rows$width)
    r[, x2 := cumsum(width), by = w][, `:=`(x1 = x2 - width, Lw = sum(width)), by = w][, `:=`(cn = wcn[w], circ = circ[w])]
    r[, I := ifelse(circ, width * 2 * m$F(Lw / 2), m$G(x2) - m$G(x1) + m$G(Lw - x1) - m$G(Lw - x2))]
    # remaining graph copies of each node, at reference positions
    used = data.table(node.id = abs(unlist(sn)), u = rep(wcn, base::lengths(sn)))[, .(used = sum(u)), by = node.id]
    nd = merge(nd, used, by = 'node.id', all.x = TRUE)[is.na(used), used := 0][, res := cn - used]
    if (any(nd$res < 0)){
        warning(sum(nd$res < 0), ' nodes are used by the walks more often than their graph CN; their remaining copies are set to 0')
        nd[res < 0, res := 0]
    }
    nd[, Iref := .conserve_refint(m, start - 1, end, L)]
    og = conserve_outside(walks$graph)
    S = sqrt((r[, sum(cn * (K * width - I))] + nd[, sum(res * (K * (end - start + 1) - Iref))] + og$Q_out) / D)
    # per tile: the walks' copies + the node's remaining copies
    tt = tiled.target
    num = numeric(nrow(tt))
    bw = r[!is.na(orig), .(v = sum(cn * (K * width - I))), by = orig]
    num[bw$orig] = bw$v
    res = nd$res[match(tt$node.id, nd$node.id)]; res[is.na(res)] = 0
    ttL = .conserve_lengths(tt$seqnames, sl, tt$end)
    num = num + res * (K * tt$width - .conserve_refint(m, tt$start - 1, tt$end, ttL))
    list(B = num / (tt$width * D * S), S = S, outside_bp = og$outside_bp)
}

#' Contribution of the genome outside a graph's nodes to D S^2 of the conserved trans term
#'
#' Chromosomes 1-22 + X at the graph's own lengths (else the calibration lengths), minus the graph's nodes, as 2 reference
#' copies: Q_out = sum 2 (K w - I_ref). Used by trans_weights() and written by prep_data_for_py() for kinf.
#' @param gg gGraph
#' @return list(Q_out, outside_bp)
conserve_outside = function(gg){
    m = conserve_model(); K = m$K
    nd = data.table::as.data.table(gg$nodes$dt)[, .(chr = as.character(seqnames), start, end)]
    sl = GenomeInfoDb::seqlengths(gg$gr)
    key = names(m$chrlen); gchr = unique(nd$chr)
    style = ifelse(paste0('chr', key) %in% gchr & !(key %in% gchr), paste0('chr', key), key)
    Lfill = .conserve_lengths(style, sl)
    out = data.table::as.data.table(GenomicRanges::setdiff(GenomicRanges::GRanges(style, IRanges::IRanges(1, Lfill)),
                                                           GenomicRanges::GRanges(nd$chr, IRanges::IRanges(nd$start, nd$end))))
    out[, L := Lfill[match(as.character(seqnames), style)]]
    list(Q_out = if (nrow(out)) out[, sum(2 * (K * width - .conserve_refint(m, start - 1, end, L)))] else 0,
         outside_bp = sum(out$width))
}

#' Per-tile trans weights of a diploid reference genome (2 copies at reference positions, calibration S_ref)
#'
#' @param tiled.target data.table with seqnames, start, end, width
#' @param seqlengths optional named chromosome lengths (else the hg38 calibration lengths; else the largest tile end)
reference_trans_weights = function(tiled.target, seqlengths = NULL){
    m = conserve_model()
    L = .conserve_lengths(tiled.target$seqnames, seqlengths, tiled.target$end)
    2 * (m$K * tiled.target$width - .conserve_refint(m, tiled.target$start - 1, tiled.target$end, L)) /
        (tiled.target$width * m$D * m$S_ref)
}
