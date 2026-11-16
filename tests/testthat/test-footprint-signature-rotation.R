# .local_footprint_signature(): a closed cycle wholly inside the footprint must get the same
# signature whatever its rotation, orientation or h1/h2 labelling once a locus map is given.
# All node ids below are synthetic.

test_that(".canonical_cycle is invariant to rotation and orientation", {
  v <- c(3L, -5L, 7L, 2L)
  canon <- .canonical_cycle(v)
  for (i in 0:3) {
    r <- v[(seq_along(v) + i - 1L) %% 4L + 1L]
    expect_identical(.canonical_cycle(r), canon)
    expect_identical(.canonical_cycle(-rev(r)), canon)
  }
  expect_true(nzchar(canon))
  expect_false(identical(.canonical_cycle(c(3L, -5L, 7L, 9L)), canon))
})

test_that("closed in-footprint cycle: relabelled h1/h2 copies share a signature under a locus map", {
  ids <- 1:8                       # nodes 5..8 are the second copies of loci 1..4
  ctx <- list(
    graph = list(nodes = list(dt = data.frame(node.id = ids, width = 100))),
    use.nodes = ids, edit_widths = NULL,
    split_ends = list(left = setNames(rep(TRUE, 8), ids), right = setNames(rep(TRUE, 8), ids)))
  locus <- c(1:4, 1:4)
  w1 <- list(snode.id = list(c(1L, 2L, 3L, 4L)), circular = TRUE)
  w2 <- list(snode.id = list(c(1L, 6L, 3L, 8L)), circular = TRUE)   # same loci, other copies
  w3 <- list(snode.id = list(c(5L, 2L, 7L, 4L)), circular = TRUE)   # starts at a different node
  s <- function(w, loc = locus) .local_footprint_signature(w, ctx, loc)
  expect_identical(s(w1), s(w2))
  expect_identical(s(w1), s(w3))
  # without a locus map nothing changes: the copies stay distinct
  expect_false(identical(s(w1, NULL), s(w2, NULL)))
  # a different cycle of loci stays distinct
  expect_false(identical(s(w1), s(list(snode.id = list(c(1L, 3L, 2L, 4L)), circular = TRUE))))
})
