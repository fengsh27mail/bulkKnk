#' Execute network-wide virtual gene knockout (core RWR perturbation engine)
#'
#' Simulates gene knockout using a Random Walk with Restart (RWR) algorithm.
#' The target gene's expression is set to 0 and all incoming/outgoing edges are
#' severed, then iterative propagation computes the new network steady state.
#'
#' @param target_genes Character vector or NULL. Gene(s) to knock out.
#'   Set to NULL for a no-intervention baseline run.
#' @param E_init Numeric vector. Baseline expression for all genes (typically group mean).
#' @param adj_matrix Matrix. Normalized directed regulatory network (row sums <= 1).
#' @param alpha Numeric. Restart probability for RWR, default 0.7.
#' @param steps Integer. Maximum iterations, default 30.
#' @param tol Numeric. Convergence tolerance. Stops early if max|delta| < tol, default 1e-6.
#' @param check_convergence Logical. Print per-step convergence error, default FALSE.
#'
#' @return Numeric vector of post-intervention steady-state expression values.
#' @export
run_virtual_knockout <- function(target_genes = NULL,
                                 E_init,
                                 adj_matrix,
                                 alpha = 0.7,
                                 steps = 30,
                                 tol = 1e-6,
                                 check_convergence = FALSE) {

  # --- Input validation ---
  if (!is.numeric(E_init)) stop("E_init must be a numeric vector.")
  if (!is.matrix(adj_matrix)) stop("adj_matrix must be a matrix.")
  if (any(is.na(E_init))) stop("E_init must not contain NA values.")
  if (any(is.na(adj_matrix))) stop("adj_matrix must not contain NA values.")

  if (length(E_init) != nrow(adj_matrix)) {
    stop(sprintf("E_init length (%d) must match adj_matrix rows (%d).",
                 length(E_init), nrow(adj_matrix)))
  }
  if (!is.null(rownames(adj_matrix)) && !is.null(names(E_init))) {
    if (!all(names(E_init) %in% rownames(adj_matrix)) ||
        !all(rownames(adj_matrix) %in% names(E_init))) {
      stop("Names of E_init and adj_matrix must match.")
    }
  }
  if (alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1).")
  if (steps < 1) stop("steps must be >= 1.")
  if (tol < 0) stop("tol must be > 0.")

  E_anchor  <- E_init
  E_current <- E_init
  W_current <- adj_matrix

  if (!is.null(target_genes)) {
    valid_targets <- intersect(target_genes, names(E_init))

    if (length(valid_targets) == 0) {
      stop("None of the target genes are found in the network. Check gene names.")
    }

    if (length(valid_targets) < length(target_genes)) {
      warning(sprintf(
        "Only %d/%d target genes found in network. Missing: %s",
        length(valid_targets), length(target_genes),
        paste(setdiff(target_genes, valid_targets), collapse = ", ")
      ))
    }

    message(sprintf(">>> [RWR Engine] Virtual knockout: %d target(s)", length(valid_targets)))

    # Sever target gene: zero out expression and all edges
    E_anchor[valid_targets] <- 0
    E_current[valid_targets] <- 0
    W_current[, valid_targets] <- 0   # Cut input edges to target
    W_current[valid_targets, ] <- 0   # Cut output edges from target
  }

  # --- RWR iteration with convergence check ---
  converged <- FALSE
  for (i in seq_len(steps)) {
    E_new <- alpha * (W_current %*% E_current) + (1 - alpha) * E_anchor
    E_new <- as.vector(E_new)
    names(E_new) <- names(E_init)

    # Enforce non-negative expression
    E_new[E_new < 0] <- 0

    max_err <- max(abs(E_new - E_current))

    if (check_convergence) {
      message(sprintf("   Step %02d/%d | Max error: %.2e", i, steps, max_err))
    }

    if (max_err < tol) {
      message(sprintf(">>> [RWR Engine] Converged at step %d (error %.2e < tol %.0e)", i, max_err, tol))
      converged <- TRUE
      break
    }
    E_current <- E_new
  }

  if (!converged) {
    message(sprintf(">>> [RWR Engine] Did not fully converge after %d steps (final error: %.2e)", steps, max_err))
  }

  return(E_current)
}


#' Batch virtual knockout for a sample cohort
#'
#' Automatically extracts target phenotype samples, runs RWR knockout on group mean,
#' then broadcasts the perturbation delta back to individual samples.
#'
#' @param expr_matrix Matrix. Expression matrix (rows=samples, cols=genes).
#' @param pheno_vec Vector. Sample phenotype labels, length must equal nrow(expr_matrix).
#' @param target_pheno String. Phenotype to perturb (e.g. "Tumor").
#' @param target_genes Character vector. Gene(s) to knock out.
#' @param adj_matrix Matrix. Normalized regulatory network.
#' @param method String. Bridge strategy: "shift" (recommended) or "replace".
#'
#' @return List with WT and KO matrices.
#' @export
batch_virtual_knockout <- function(expr_matrix, pheno_vec, target_pheno,
                                   target_genes, adj_matrix, method = "shift") {

  if (!is.matrix(expr_matrix)) stop("expr_matrix must be a matrix.")
  if (length(pheno_vec) != nrow(expr_matrix)) {
    stop("pheno_vec length must match nrow(expr_matrix).")
  }
  if (!target_pheno %in% pheno_vec) {
    stop(sprintf("Phenotype '%s' not found in pheno_vec.", target_pheno))
  }

  target_idx <- which(pheno_vec == target_pheno)
  WT_mat <- expr_matrix[target_idx, , drop = FALSE]
  message(sprintf(">>> [Batch KO] Selected %d [%s] samples for perturbation",
                  length(target_idx), target_pheno))

  # Compute group baseline
  E_baseline <- colMeans(WT_mat)

  # Run RWR knockout on group mean
  ko_vector <- run_virtual_knockout(
    target_genes = target_genes,
    E_init       = E_baseline,
    adj_matrix   = adj_matrix
  )

  # Broadcast delta to individual samples
  KO_mat <- build_ko_matrix(
    WT_matrix   = WT_mat,
    ko_vector   = ko_vector,
    WT_baseline = E_baseline,
    method      = method
  )

  message(">>> [Batch KO] Knockout matrix built.")

  return(list(WT = WT_mat, KO = KO_mat))
}
