# .local_footprint_signature(): a closed cycle wholly inside the footprint must get the same
# signature whatever its rotation, orientation or h1/h2 labelling (the last once a locus map
# is given), and a different one from a linear fragment with the same nodes.
# All node ids below are synthetic.

sig_ctx <- function(ids) list(
  graph = list(nodes = list(dt = data.frame(node.id = ids, width = 100))),
  use.nodes = ids, edit_widths = NULL,
  split_ends = list(left = setNames(rep(TRUE, length(ids)), ids),
                    right = setNames(rep(TRUE, length(ids)), ids)))

test_that("closed in-footprint cycle: signature is invariant to rotation and orientation", {
  ctx <- sig_ctx(1:9)
  v <- c(3L, -5L, 7L, 2L)
  s <- function(x, loc = NULL) .local_footprint_signature(list(snode.id = list(x), circular = TRUE), ctx, loc)
  for (loc in list(NULL, 1:9)) {
    canon <- s(v, loc)
    for (i in 0:3) {
      r <- v[(seq_along(v) + i - 1L) %% 4L + 1L]
      expect_identical(s(r, loc), canon)
      expect_identical(s(-rev(r), loc), canon)
    }
    expect_false(identical(s(c(3L, -5L, 7L, 9L), loc), canon))
  }
})

test_that("closed in-footprint cycle and linear fragment with the same nodes differ", {
  ctx <- sig_ctx(1:4)
  cyc <- list(snode.id = list(1:4), circular = TRUE)
  lin <- list(snode.id = list(1:4), circular = FALSE)
  expect_false(identical(.local_footprint_signature(cyc, ctx), .local_footprint_signature(lin, ctx)))
  expect_false(identical(.local_footprint_signature(cyc, ctx, 1:4), .local_footprint_signature(lin, ctx, 1:4)))
})

test_that("closed in-footprint cycle: relabelled h1/h2 copies share a signature under a locus map", {
  ids <- 1:8                       # nodes 5..8 are the second copies of loci 1..4
  ctx <- sig_ctx(ids)
  locus <- c(1:4, 1:4)
  w1 <- list(snode.id = list(c(1L, 2L, 3L, 4L)), circular = TRUE)
  w2 <- list(snode.id = list(c(1L, 6L, 3L, 8L)), circular = TRUE)   # same loci, other copies
  w3 <- list(snode.id = list(c(5L, 2L, 7L, 4L)), circular = TRUE)   # starts at a different node
  s <- function(w, loc = locus) .local_footprint_signature(w, ctx, loc)
  expect_identical(s(w1), s(w2))
  expect_identical(s(w1), s(w3))
  # without a locus map the copies stay distinct
  expect_false(identical(s(w1, NULL), s(w2, NULL)))
  # a different cycle of loci stays distinct
  expect_false(identical(s(w1), s(list(snode.id = list(c(1L, 3L, 2L, 4L)), circular = TRUE))))
})
