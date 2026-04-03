# tests/testthat/test-network-utils.R
# Unit tests for DPI pruning and KO matrix bridging

# --- prune_network_dpi ---

test_that("prune_network_dpi removes indirect edges", {
  # Triangle: A-B=0.9, B-C=0.8, A-C=0.1 (weak indirect edge)
  W <- matrix(c(0,   0.9, 0.1,
                0.9, 0,   0.8,
                0.1, 0.8, 0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A","B","C"), c("A","B","C")))

  result <- prune_network_dpi(W, eps = 0, verbose = FALSE)

  # A-C should be removed (weakest in triangle)
  expect_equal(result["A", "C"], 0)
  expect_equal(result["C", "A"], 0)

  # A-B and B-C should survive
  expect_gt(result["A", "B"], 0)
  expect_gt(result["B", "C"], 0)
})

test_that("prune_network_dpi preserves edges in chain", {
  # Chain: A -> B -> C (no triangle, nothing to prune)
  W <- matrix(c(0,   0.8, 0,
                0,   0,   0.6,
                0,   0,   0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A","B","C"), c("A","B","C")))

  result <- prune_network_dpi(W, eps = 0, verbose = FALSE)

  # All original edges should survive
  expect_equal(result["A","B"], 0.8)
  expect_equal(result["B","C"], 0.6)
})

test_that("prune_network_dpi rejects non-matrix input", {
  expect_error(prune_network_dpi(data.frame(a=1)), "matrix")
})

test_that("prune_network_dpi_fast down-weights indirect edges", {
  W <- matrix(c(0,   0.8, 0.3,
                0,   0,   0.9,
                0,   0,   0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A","B","C"), c("A","B","C")))

  result <- prune_network_dpi_fast(W, dpi_threshold = 0.5)

  # A-C is indirect (A->B->C is stronger), should be down-weighted
  expect_lt(result["A", "C"], W["A", "C"])
  # A-B and B-C should remain unchanged (no stronger path bypasses them)
  expect_equal(result["A", "B"], W["A", "B"])
  expect_equal(result["B", "C"], W["B", "C"])
})

# --- build_ko_matrix ---

test_that("build_ko_matrix shift method works", {
  WT <- matrix(c(10, 5, 3,
                 8,  4, 2),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("S1","S2"), c("G1","G2","G3")))

  ko_vec <- c(G1 = 0, G2 = 4, G3 = 2)  # G1 knocked out
  baseline <- c(G1 = 9, G2 = 4.5, G3 = 2.5)

  result <- build_ko_matrix(WT, ko_vec, baseline, method = "shift")

  # G1 should be reduced by delta (0 - 9 = -9 shift)
  expect_equal(result["S1", "G1"], max(0, 10 - 9))
  expect_equal(result["S2", "G1"], max(0, 8 - 9))  # clamped to 0

  # G2 unchanged (ko - baseline = 4 - 4.5 = -0.5 shift)
  expect_equal(result["S1", "G2"], 5 - 0.5)
})

test_that("build_ko_matrix replace method works", {
  WT <- matrix(c(10, 5,
                 8,  4),
               nrow = 2, byrow = TRUE,
               dimnames = list(c("S1","S2"), c("G1","G2")))

  ko_vec <- c(G1 = 0, G2 = 3)
  baseline <- c(G1 = 9, G2 = 4.5)

  result <- build_ko_matrix(WT, ko_vec, baseline, method = "replace")

  # All rows should be identical to ko_vec
  expect_equal(result["S1", "G1"], 0)
  expect_equal(result["S2", "G1"], 0)
  expect_equal(result["S1", "G2"], 3)
  expect_equal(result["S2", "G2"], 3)
})

test_that("build_ko_matrix rejects invalid method", {
  WT <- matrix(1, nrow = 2, ncol = 2, dimnames = list(c("S1","S2"), c("G1","G2")))
  ko_vec <- c(G1 = 0, G2 = 1)
  baseline <- c(G1 = 1, G2 = 1)

  expect_error(build_ko_matrix(WT, ko_vec, baseline, method = "invalid"), "shift.*replace")
})

test_that("build_ko_matrix rejects dimension mismatch", {
  WT <- matrix(1, nrow = 2, ncol = 3, dimnames = list(c("S1","S2"), c("G1","G2","G3")))
  ko_vec <- c(G1 = 0, G2 = 1)  # Only 2 genes vs 3
  baseline <- c(G1 = 1, G2 = 1)

  expect_error(build_ko_matrix(WT, ko_vec, baseline), "must match|cols")
})

# --- refine_network_aracne (wrapper) ---

test_that("refine_network_aracne delegates to prune_network_dpi", {
  W <- matrix(c(0,   0.9, 0.1,
                0.9, 0,   0.8,
                0.1, 0.8, 0),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("A","B","C"), c("A","B","C")))

  result1 <- refine_network_aracne(W, eps = 0)
  result2 <- prune_network_dpi(W, eps = 0, verbose = FALSE)

  expect_identical(result1, result2)
})
