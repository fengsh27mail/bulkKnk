# tests/testthat/test-ai-predictor.R
# Unit tests for ssGSEA scoring and AI phenotype prediction

test_that(".ssgsea_score returns correct dimensions", {
  # Small expression matrix: 10 genes x 3 samples
  set.seed(42)
  genes <- paste0("GENE", 1:10)
  samples <- c("S1", "S2", "S3")
  expr_mat <- matrix(rnorm(30), nrow = 10, ncol = 3,
                     dimnames = list(genes, samples))

  gene_sets <- list(
    set1 = c("GENE1", "GENE2", "GENE3", "GENE4"),
    set2 = c("GENE7", "GENE8", "GENE9", "GENE10")
  )

  scores <- bulkKnk:::.ssgsea_score(expr_mat, gene_sets)

  expect_true(is.matrix(scores))
  expect_equal(nrow(scores), 2)   # 2 gene sets
  expect_equal(ncol(scores), 3)   # 3 samples
  expect_true(all(rownames(scores) == c("set1", "set2")))
})

test_that(".ssgsea_score skips small gene sets", {
  genes <- paste0("GENE", 1:10)
  expr_mat <- matrix(rnorm(20), nrow = 10, ncol = 2,
                     dimnames = list(genes, c("S1","S2")))

  gene_sets <- list(
    big_set = c("GENE1", "GENE2", "GENE3", "GENE4"),
    tiny_set = c("GENE1")  # Only 1 gene, below threshold of 3
  )

  scores <- bulkKnk:::.ssgsea_score(expr_mat, gene_sets)

  # tiny_set should be NA
  expect_true(is.na(scores["tiny_set", "S1"]))
  # big_set should have a value
  expect_false(is.na(scores["big_set", "S1"]))
})

test_that(".get_builtin_signatures returns expected structure", {
  sigs <- bulkKnk:::.get_builtin_signatures()

  expect_type(sigs, "list")
  expect_true(is.list(sigs))
  expect_true(length(sigs) >= 20)  # At least 20 signature categories

  # Check key signatures exist
  expect_true("immune_CD8T" %in% names(sigs))
  expect_true("hallmark_Proliferation" %in% names(sigs))
  expect_true("drug_EGFR_pathway" %in% names(sigs))

  # Each should be a character vector
  for (nm in names(sigs)) {
    expect_type(sigs[[nm]], "character")
    expect_true(length(sigs[[nm]]) >= 3)
  }
})
