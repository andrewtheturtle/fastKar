## Contact-conserving trans term (R/conserve.R): every copy-bp has the same total contacts.

## per-tile row sums of an upper-triangle pixel table on n tiles (each off-diagonal pixel counts for both tiles)
row_totals = function(d, n) {
    off = d$i != d$j
    tabulate_w = function(k, w) { v = numeric(n); s = tapply(w, k, sum); v[as.integer(names(s))] = s; v }
    tabulate_w(d$i, d$value) + tabulate_w(d$j[off], d$value[off])
}

test_that('conserve: total contacts per copy-bp are the same in every simulated bin (linear copies and a circle)', {
    ## chr1 of 10 Mb as nodes A-B-C; two linear copies A-B-C and a 2 Mb circle on B at CN 5 (node CN = copies used)
    nodes = GenomicRanges::GRanges('1', IRanges::IRanges(c(1, 4e6 + 1, 6e6 + 1), c(4e6, 6e6, 1e7)), seqlengths = c('1' = 1e7))
    edges = data.table(n1 = c(1, 2, 2), n1.side = 'right', n2 = c(2, 3, 2), n2.side = 'left')
    gg = gGnome::gG(nodes = nodes, edges = edges)
    gg$nodes$mark(cn = c(2, 7, 2))
    gw = gGnome::gW(graph = gg, snode.id = list(c(1, 2, 3), 2), circular = c(FALSE, TRUE))
    gw$set(cn = c(2, 5))

    pp = prep_for_sim(gw, gg$footprint, 5e5)
    tt = pp$tiled.target; n = nrow(tt)
    sim = function(interchr) simulate_walks(gw, tt, pp$widthdt, if.interchr = interchr, gm.out = FALSE, conserve = TRUE)$dat
    in.map = row_totals(sim(TRUE), n)                      # cis + trans to the copies inside the map
    ## trans to the rest of the genome (2 reference copies of chr2-22 + X), not in the map: D w_i B_i (S - sum_j w_j B_j)
    tw = trans_weights(gw, tt, pp$widthdt)
    expect_equal(tw$outside_bp, sum(conserve_chrlen_hg38[-1]))
    outside = conserve_model()$D * tt$width * tw$B * (tw$S - sum(tt$width * tw$B))
    per.bp = (in.map + outside) / (tt$cn * tt$width)
    ## the same in every bin; contacts below 10 kb are outside the conserved budget, so the bins at the walk ends
    ## (which miss part of them) sit ~4e-4 low; every other bin, circle included, is exact (the kernel and conserve_model()
    ## share F and G)
    rel = per.bp / stats::median(per.bp) - 1
    ends = c(1, n)
    expect_lt(max(abs(rel[ends])), 1e-3)
    expect_lt(max(abs(rel[-ends])), 1e-10)
})
