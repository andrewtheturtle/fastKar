# .walks_to_wiring() and .map_truth_walks(): a truth karyotype can be re-created on another
# copy of the same graph (synthetic toy graph only).

toy_truth_setup <- function() {
  set.seed(7)
  gg <- load_fixture("toy_bfb")$gg
  wiring <- gg.to.wiring(gg)
  ie <- wiring$internal.edges
  le <- wiring$loose.ends
  ie[, right := if (.N > 1) sample(right) else right, by = n]
  list(gg = gg, ie = ie, le = le, walks = traverse_graph_cpp(ie, le))
}

hash_w <- function(w) hash_karyotype_cpp(w$snode.id, w$circular)

test_that(".walks_to_wiring reproduces a walk set from a differently wired graph", {
  skip_if_no_loosefix()
  s <- toy_truth_setup()
  other <- data.table::copy(s$ie)
  set.seed(99)
  other[, right := if (.N > 1) sample(right) else right, by = n]
  rewired <- opt(".walks_to_wiring")(s$walks, other)
  expect_equal(nrow(rewired), nrow(s$ie))
  expect_equal(hash_w(traverse_graph_cpp(rewired, s$le)), hash_w(s$walks))
})

test_that(".map_truth_walks is the identity between two copies of one graph", {
  skip_if_no_loosefix()
  s <- toy_truth_setup()
  mapped <- opt(".map_truth_walks")(s$walks, s$gg, s$gg)
  expect_equal(mapped$snode.id, s$walks$snode.id)
  expect_equal(mapped$circular, s$walks$circular)
})

test_that(".walks_to_wiring errors when the graph lacks the needed edge copies", {
  skip_if_no_loosefix()
  s <- toy_truth_setup()
  bad <- s$walks
  bad$snode.id[[1]] <- rev(bad$snode.id[[1]])[-1]
  expect_error(opt(".walks_to_wiring")(bad, s$ie))
})
