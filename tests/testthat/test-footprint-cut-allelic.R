## identifiability_context(graph_mode = 'footprint') cuts the graph to the buffered footprint. On an allelic
## graph the h1/h2 copies of a locus share coordinates; the cut must keep them as separate nodes with their own
## CN and keep allele-specific edges (gGraph$disjoin() merges same-coordinate nodes unless collapse = FALSE).
## Toy graph, synthetic coordinates: 3 loci x 2 haplotypes on chromosome 1, one h1 -> h2 junction.

toy_allelic_graph = function() {
    loc = GenomicRanges::GRanges('1', IRanges::IRanges(c(1, 100001, 200001), c(100000, 200000, 300000)),
                                 seqlengths = c('1' = 300000))
    nd = c(loc, loc)
    nd$haplotype = rep(c('h1', 'h2'), each = 3)
    nd$cn = c(2, 1, 1, 1, 1, 2)              # node ids: h1 = 1:3, h2 = 4:6
    ed = data.table::data.table(n1 = c(1, 2, 4, 5, 1), n1.side = 'right', n2 = c(2, 3, 5, 6, 6), n2.side = 'left',
                                cn = 1, type = c('REF', 'REF', 'REF', 'REF', 'ALT'))
    gGnome::gG(nodes = nd, edges = ed)
}

cut_tables = function(g) {
    n = g$nodes$dt
    n$key = paste(n$seqnames, n$start, n$end, n$haplotype)
    e = g$edges$dt
    e$key = paste(n$key[match(e$n1, n$node.id)], e$n1.side, n$key[match(e$n2, n$node.id)], e$n2.side, e$type)
    list(nodes = n, edges = e)
}

test_that("footprint cut keeps h1/h2 nodes and allele-specific edges on an allelic graph", {
    skip_if_no_loosefix()
    g = toy_allelic_graph()
    fp = GenomicRanges::GRanges('1', IRanges::IRanges(90001, 210000))   # buffered: 1:66001-234000
    ctx = suppressWarnings(opt('identifiability_context')(g, footprint = fp, graph_mode = 'footprint'))
    tb = cut_tables(ctx$graph)
    expect_equal(nrow(tb$nodes), 6)
    expect_setequal(tb$nodes$haplotype, c('h1', 'h2'))
    ## per-haplotype CN of each locus piece is unchanged
    want = c('1 66001 100000 h1' = 2, '1 66001 100000 h2' = 1, '1 100001 200000 h1' = 1,
             '1 100001 200000 h2' = 1, '1 200001 234000 h1' = 1, '1 200001 234000 h2' = 2)
    expect_equal(setNames(tb$nodes$cn, tb$nodes$key)[names(want)], want)
    ## the h1 -> h2 junction still joins the h1 copy of locus 1 to the h2 copy of locus 3
    alt = tb$edges[type == 'ALT']
    expect_equal(alt$key, '1 66001 100000 h1 right 1 200001 234000 h2 left ALT')
    expect_equal(alt$cn, 1)
    ## REF edges stay within a haplotype, with CN 1 each (not merged into CN 2)
    ref = tb$edges[type == 'REF']
    expect_equal(nrow(ref), 4)
    expect_true(all(ref$cn == 1))
})

test_that("footprint cut of a graph without overlapping nodes or parallel edges is unchanged by collapse = FALSE", {
    skip_if_no_loosefix()
    g0 = toy_allelic_graph()
    g = gGnome::gG(nodes = g0$nodes$gr[1:3, 'cn'],
                   edges = data.table::data.table(n1 = 1:2, n1.side = 'right', n2 = 2:3, n2.side = 'left', cn = 1, type = 'REF'))
    fp = GenomicRanges::GRanges('1', IRanges::IRanges(90001, 210000))
    ctx = suppressWarnings(opt('identifiability_context')(g, footprint = fp, graph_mode = 'footprint'))
    old = suppressWarnings(opt('loosefix')(g$copy$disjoin(ctx$target_region) %&% ctx$target_region))
    expect_identical(gUtils::gr.string(ctx$graph$nodes$gr), gUtils::gr.string(old$nodes$gr))
    expect_equal(ctx$graph$nodes$dt$cn, old$nodes$dt$cn)
    expect_equal(ctx$graph$edges$dt[, paste(n1, n1.side, n2, n2.side, cn)],
                 old$edges$dt[, paste(n1, n1.side, n2, n2.side, cn)])
})
