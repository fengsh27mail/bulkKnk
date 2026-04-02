# tests/testthat/test-run_virtual_knockout.R
# Unit tests for the core RWR virtual knockout engine

test_that("run_virtual_knockout basic knockout works", {
  # Build a small 5-gene network: A -> B -> C, A -> D, D -> E
  genes <- c("A", "B", "C", "D", "E")
  W <- matrix(0, nrow = 5, ncol = 5, dimnames = list(genes, genes))
  W["A", "B"] <- 0.5
  W["B", "C"] <- 0.5
  W["A", "D"] <- 0.3
  W["D", "E"] <- 0.4

  E <- c(A = 10, B = 5, C = 3, D = 4, E = 2)

  # Knock out gene A
  result <- run_virtual_knockout(
    target_genes = "A",
    E_init = E,
    adj_matrix = W,
    steps = 50,
    tol = 1e-10
  )

  # A should be 0
  expect_equal(result["A"], 0)

  # B and D should decrease since A was their only regulator
  expect_lt(result["B"], E["B"])
  expect_lt(result["D"], E["D"])

  # C should decrease downstream of B
  expect_lt(result["C"], E["C"])
})

test_that("run_virtual_knockout with NULL target returns baseline", {
  genes <- c("X", "Y", "Z")
  W <- matrix(c(0, 0.5, 0.3,
                0, 0,   0.2,
                0, 0,   0),
              nrow = 3, byrow = TRUE,
              dimnames = list(genes, genes))
  E <- c(X = 5, Y = 3, Z = 1)

  # No knockout — should converge to some steady state
  result <- run_virtual_knockout(
    target_genes = NULL,
    E_init = E,
    adj_matrix = W,
    steps = 100,
    tol = 1e-10
  )

  expect_type(result, "double")
  expect_length(result, 3)
  expect_true(all(names(result) == genes))
})

test_that("run_virtual_knockout rejects invalid inputs", {
  W <- matrix(0.5, nrow = 3, ncol = 3, dimnames = list(c("A","B","C"), c("A","B","C")))
  E <- c(A = 1, B = 2, C = 3)

  # Wrong length E_init
  expect_error(run_virtual_knockout(E_init = c(1,2), adj_matrix = W),
               "must match")

  # Non-numeric E_init
  expect_error(run_virtual_knockout(E_init = c("a","b","c"), adj_matrix = W),
               "numeric")

  # alpha out of range
  expect_error(run_virtual_knockout(E_init = E, adj_matrix = W, alpha = 0),
               "alpha")
  expect_error(run_virtual_knockout(E_init = E, adj_matrix = W, alpha = 1),
               "alpha")

  # Invalid target gene
  expect_error(run_virtual_knockout(target_genes = "NOT_HERE", E_init = E, adj_matrix = W),
               "not found")
})

test_that("run_virtual_knockout convergence detection works", {
  W <- matrix(0.1, nrow = 3, ncol = 3,
              dimnames = list(c("A","B","C"), c("A","B","C")))
  diag(W) <- 0
  E <- c(A = 5, B = 3, C = 1)

  # With generous tol, should converge quickly
  result <- run_virtual_knockout(
    target_genes = "A",
    E_init = E,
    adj_matrix = W,
    steps = 100,
    tol = 0.01
  )

  expect_type(result, "double")
  expect_equal(result["A"], 0)
})

test_that("run_virtual_knockout multi-target knockout works", {
  genes <- c("A", "B", "C", "D")
  W <- matrix(0.2, nrow = 4, ncol = 4, dimnames = list(genes, genes))
  diag(W) <- 0
  E <- c(A = 10, B = 8, C = 6, D = 4)

  result <- run_virtual_knockout(
    target_genes = c("A", "B"),
    E_init = E,
    adj_matrix = W,
    steps = 50
  )

  expect_equal(result["A"], 0)
  expect_equal(result["B"], 0)
})
