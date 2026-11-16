## Unit tests of the C++ forward-simulation kernel (src/forward_kernel.cpp), oracle-independent:
## the spline port against R's splinefun, the lift semantics against a dense t(T) S T, and the
## pixel id formula against make_template_dat_cpp().

kernel_args = function(model = 0) forward_kernel_args(model)

rel_err = function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))

test_that('C++ distance density matches distance_density() on a dense log grid, knots and joins', {
    ka = kernel_args()
    x = ka$z$x
    eps = 1e-9
    lgrid = c(seq(log(1), log(1e9), length.out = 20001),
              x, x - eps, x + eps, x * (1 - 1e-15), x * (1 + 1e-15),
              ka$lmin + c(-1e-6, -1e-12, 0, 1e-12, 1e-6),
              log(1e8) + c(-1e-6, 0, 1e-6),
              log(c(1e4 - 1, 1e4, 1e4 + 1, 1e8 - 1, 1e8, 1e8 + 1, 2.49e8, 5e8)))
    d = exp(lgrid)
    d = d[d >= 1 & d <= 1e9]
    cpp = forward_density_cpp(d, ka$z, ka$lmin, ka$slope)
    r = distance_density(d)
    expect_true(all(is.finite(cpp)))
    expect_lte(rel_err(cpp, r), 1e-13)
    ## the spline is read from the package data, not hard-coded
    expect_identical(ka$z, environment(fastKar::splineobj$distance_spline)$z)
    expect_identical(ka$slope, fastKar::splineobj$distance_spline(ka$lmin, deriv = 1))
})

## dense reference for one hand-built walk: simulate_map() on the local rows, symmetrize,
## lift with a dense 0/1 matrix, upper triangle in template order
dense_lift = function(orig, width, circular = FALSE, model = 0, n_ref) {
    walk.dt = data.table(orig.id = orig, compartment = 'A', width = width, tile.id = seq_along(orig))
    walk.dt[, end := cumsum(width)][, start := end - width + 1][, mid := (start + end) / 2]
    loc = simulate_map(walk.dt, circular, model = model)
    n = length(orig)
    L = matrix(0, n, n)
    L[cbind(loc$i, loc$j)] = loc$value
    S = L + t(L)
    diag(S) = diag(L)
    Tm = matrix(0, n, n_ref)
    ok = !is.na(orig)
    Tm[cbind(which(ok), orig[ok])] = 1
    ## base::`%*%`: gUtils (attached via Depends) masks %*% with a GRanges operator
    R = base::`%*%`(base::`%*%`(t(Tm), S), Tm)
    tmpl = make_template_dat_cpp(data.frame(width = rep(1, n_ref)))
    R[cbind(tmpl$i, tmpl$j)]
}

kernel_walks = function(orig_list, width_list, circular, cn, n_ref, model = 0) {
    ka = kernel_args(model)
    ptr = as.integer(c(0, cumsum(base::lengths(orig_list))))
    forward_walks_cpp(ptr, as.integer(unlist(orig_list)), as.numeric(unlist(width_list)), FALSE,
                      as.logical(circular), as.numeric(cn), as.integer(n_ref), ka$z, ka$lmin, ka$slope,
                      ka$maxval, ka$maxval_density, ka$model, ka$lr_denom, ka$diagf)
}

lift_cases = list(
    revisit_same_tile = list(orig = c(1, 2, 2, 3), width = c(1e5, 2, 2, 1e5), circ = FALSE, n_ref = 3),
    revisit_far = list(orig = c(1, 2, 3, 1, 2), width = c(5e4, 5e4, 5e4, 5e4, 5e4), circ = FALSE, n_ref = 3),
    minus_strand_order = list(orig = c(4, 3, 2, 1), width = c(1e4, 2e4, 3e4, 4e4), circ = FALSE, n_ref = 4),
    circular_wrap = list(orig = c(1, 2, 3, 4, 5), width = c(1e6, 1e4, 1e4, 1e4, 1e6), circ = TRUE, n_ref = 5),
    mixed_target = list(orig = c(NA, 1, NA, 2, 1, NA), width = c(3e6, 1e5, 2e8, 1e5, 1e5, 10), circ = FALSE, n_ref = 2),
    mixed_circular = list(orig = c(2, NA, 1, NA), width = c(1e5, 5e7, 1e5, 1e3), circ = TRUE, n_ref = 3),
    tail_beyond_1e8 = list(orig = c(1, NA, 2), width = c(1e4, 3e8, 1e4), circ = FALSE, n_ref = 2)
)

test_that('kernel lift reproduces the dense t(T) S T on hand-built walks (Hi-C and long-read)', {
    for (nm in names(lift_cases)) {
        cs = lift_cases[[nm]]
        for (model in c(0, 1e5)) {
            inf = paste(nm, model)
            expected = dense_lift(cs$orig, cs$width, cs$circ, model, cs$n_ref)
            got = kernel_walks(list(cs$orig), list(cs$width), cs$circ, 1, cs$n_ref, model)
            expect_length(got, cs$n_ref * (cs$n_ref + 1) / 2)
            expect_equal(got, expected, tolerance = 1e-13, info = inf)
            ## walk CN multiplies
            got3 = kernel_walks(list(cs$orig), list(cs$width), cs$circ, 3, cs$n_ref, model)
            expect_equal(got3, 3 * expected, tolerance = 1e-13, info = paste(inf, 'cn'))
        }
    }
})

test_that('revisits mapping to one tile add 2 L_ab to the diagonal pixel', {
    ## two rows of the same tile: local (1,1), (2,2) and 2 x (1,2) all land on reference (1,1)
    orig = c(1, 1); width = c(1e4, 1e4)
    walk.dt = data.table(orig.id = orig, compartment = 'A', width = width, tile.id = 1:2)
    walk.dt[, end := cumsum(width)][, start := end - width + 1][, mid := (start + end) / 2]
    loc = simulate_map(walk.dt)
    want = sum(loc[i == j]$value) + 2 * loc[i == 1 & j == 2]$value
    expect_equal(kernel_walks(list(orig), list(width), FALSE, 1, 1), want, tolerance = 1e-14)
})

## reference pixel integral: int rho(s) T(s) ds, T = overlap of tiles wa, wb at separation s (trapezoid around d),
## composite Simpson in log s on each smooth piece (1e-3 bp lower limit for touching tiles: T ~ s there)
pixel_quadrature = function(wa, wb, d, ka) {
    rho = function(s) {
        out = distance_density(pmin(s, ka$maxval))
        hi = s > ka$maxval
        out[hi] = ka$maxval_density * ka$maxval / s[hi]
        out
    }
    A = (wa + wb) / 2; B = abs(wa - wb) / 2
    Tw = function(s) pmax(0, pmin(min(wa, wb), A - abs(s - d)))
    br = sort(unique(c(max(d - A, 1e-3), d - B, d + B, d + A, exp(ka$lmin), ka$maxval)))
    br = br[br >= max(d - A, 1e-3) & br <= d + A]
    tot = 0
    for (k in seq_len(length(br) - 1)) {
        u = seq(log(br[k]), log(br[k + 1]), length.out = 20001); h = diff(u[1:2])
        s = exp(u); y = rho(s) * Tw(s) * s
        tot = tot + h / 3 * (y[1] + y[20001] + 4 * sum(y[seq(2, 20000, 2)]) + 2 * sum(y[seq(3, 19999, 2)]))
    }
    tot
}

test_that('pixel-integrated kernel matches quadrature; beyond 20 widths it is the midpoint value', {
    ka = kernel_args()
    cases = list(c(1e4, 1e4), c(1e6, 1e6), c(3.7e4, 1e6), c(5e6, 2e5), c(5e7, 5e7))
    for (w in cases) {
        wa = w[1]; wb = w[2]; wm = max(w)
        d = (wa + wb) / 2 * c(1, 1.5, 3)            # touching and near
        d = c(d, wm * c(10, 19.9, 20.1, 50))        # both sides of the exact / midpoint switch (20 widths)
        got = forward_pixel_cpp(rep(wa, length(d)), rep(wb, length(d)), d, ka$z, ka$lmin, ka$slope, ka$maxval, ka$maxval_density)
        ref = sapply(d, function(x) pixel_quadrature(wa, wb, x, ka))
        exact = d < 20 * wm
        expect_lte(max(abs(got[exact] / ref[exact] - 1)), 1e-3, label = paste('exact', wa, wb))
        expect_lte(max(abs(got[!exact] / ref[!exact] - 1)), 1e-2, label = paste('midpoint zone', wa, wb))
        mid = ifelse(d > ka$maxval, ka$maxval_density * ka$maxval / d, distance_density(pmin(d, ka$maxval))) * wa * wb
        expect_equal(got[!exact], mid[!exact], tolerance = 1e-13)
    }
})

test_that('kernel sums several walks with their CNs and is additive over walks', {
    a = lift_cases$revisit_far; b = lift_cases$minus_strand_order
    n_ref = 4
    one = kernel_walks(list(a$orig), list(a$width), FALSE, 2, n_ref)
    two = kernel_walks(list(b$orig), list(b$width), FALSE, 5, n_ref)
    both = kernel_walks(list(a$orig, b$orig), list(a$width, b$width), c(FALSE, FALSE), c(2, 5), n_ref)
    expect_equal(both, one + two, tolerance = 1e-14)
})

test_that('pixel id formula and assembled tables match make_template_dat_cpp for n = 1..200', {
    for (n in 1:200) {
        tm = make_template_dat_cpp(data.frame(width = seq_len(n) * 7))
        expect_identical(pixel_id(tm$i, tm$j, n), tm$id)
        as = forward_assemble_cpp(as.numeric(seq_len(n) * 7), rep(2, n), numeric(nrow(tm)), 1, 0, FALSE)
        expect_identical(as$i, tm$i)
        expect_identical(as$j, tm$j)
        expect_identical(as$id, tm$id)
        ## the kernel's diagonal pixel (p, p) lands on the template row with i == j == p
        if (n <= 20) {
            for (p in seq_len(n)) {
                got = kernel_walks(list(p), list(1e4), FALSE, 1, n)
                expect_identical(which(got != 0), tm$id[tm$i == p & tm$j == p])
            }
        }
    }
})

test_that('kernel reproduces the reference cumsum NA / integer overflow behaviour', {
    ka = kernel_args()
    ## NA width: rows from there on have NA coordinates and contribute nothing
    got = kernel_walks(list(c(1, NA, 2)), list(c(1e4, NA, 1e4)), FALSE, 1, 2)
    expect_equal(got, c(dense_lift(1, 1e4, FALSE, 0, 1), 0, 0), tolerance = 1e-13)
    ## integer widths overflowing .Machine$integer.max: the reference cumsum gives NA (and warns)
    w = c(1e4, 2.2e9, 1e4)
    expect_warning(got <- forward_walks_cpp(c(0L, 3L), c(1L, NA, 2L), w, TRUE, FALSE, 1, 2L, ka$z, ka$lmin,
                                            ka$slope, ka$maxval, ka$maxval_density, 0, 1, ka$diagf), 'overflow')
    expect_equal(got[2:3], c(0, 0))
    expect_gt(got[1], 0)
})
