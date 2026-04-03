# =============================================================================
#  network_utils.R
#  网络工具集：① ARACNE 数据处理不等式后处理（净化 GENIE3 网络）
#              ② KO 向量→样本矩阵桥接函数（连接 Module 2 与 Module 4）
# =============================================================================

# -----------------------------------------------------------------------------
#' 用数据处理不等式 (DPI) 净化 GENIE3 调控网络 (ARACNE-inspired pruning)
#'
#' @param weight_matrix 矩阵。GENIE3 输出的有向权重邻接矩阵（行=调控者，列=被调控者）。
#' @param eps 数值。DPI 容忍度，默认 0。
#' @param verbose 逻辑值。是否打印进度，默认 TRUE。
#'
#' @return 经 DPI 净化后的稀疏权重矩阵。
#' @export
prune_network_dpi <- function(weight_matrix, eps = 0, verbose = TRUE) {
  if (!is.matrix(weight_matrix)) stop("weight_matrix 必须是标准矩阵。")
  genes    <- rownames(weight_matrix)
  n        <- length(genes)
  pruned   <- weight_matrix
  removed  <- 0

  if (verbose) message(sprintf(">> [ARACNE-DPI] 开始净化，共 %d 个基因，分析 %d 个三角形...", n, choose(n, 3)))

  for (i in seq_len(n - 2)) {
    for (j in seq(i + 1, n - 1)) {
      for (k in seq(j + 1, n)) {
        w_ij <- max(pruned[i, j], pruned[j, i])
        w_jk <- max(pruned[j, k], pruned[k, j])
        w_ik <- max(pruned[i, k], pruned[k, i])

        if (w_ij == 0 || w_jk == 0 || w_ik == 0) next
        min_w <- min(w_ij, w_jk, w_ik)

        if (min_w == w_ij && (min_w + eps) < min(w_jk, w_ik)) {
          pruned[i, j] <- 0; pruned[j, i] <- 0; removed <- removed + 1
        } else if (min_w == w_jk && (min_w + eps) < min(w_ij, w_ik)) {
          pruned[j, k] <- 0; pruned[k, j] <- 0; removed <- removed + 1
        } else if (min_w == w_ik && (min_w + eps) < min(w_ij, w_jk)) {
          pruned[i, k] <- 0; pruned[k, i] <- 0; removed <- removed + 1
        }
      }
    }
  }
  if (verbose) message(sprintf(">> [ARACNE-DPI] 完成！剩余有效边: %d 条", sum(pruned > 0)))
  return(pruned)
}

# -----------------------------------------------------------------------------
#' 快速版 ARACNE-DPI（基于矩阵运算，适合大型网络）
#'
#' @param weight_matrix 矩阵。GENIE3 输出的有向权重邻接矩阵。
#' @param dpi_threshold 数值。间接调控降权比例，默认 0.5。
#'
#' @return 净化后的权重矩阵。
#' @export
prune_network_dpi_fast <- function(weight_matrix, dpi_threshold = 0.5) {
  if (!is.matrix(weight_matrix)) stop("weight_matrix 必须是标准矩阵。")
  n      <- nrow(weight_matrix)
  pruned <- weight_matrix

  message(">> [ARACNE-Fast] 正在进行向量化 DPI 净化...")
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j || weight_matrix[i, j] == 0) next
      w_ij <- weight_matrix[i, j]
      intermediaries <- which(weight_matrix[i, ] > w_ij & weight_matrix[, j] > w_ij)
      if (length(intermediaries) > 0) {
        pruned[i, j] <- w_ij * dpi_threshold
      }
    }
  }
  return(pruned)
}

# -----------------------------------------------------------------------------
#' 将 KO 稳态向量扩展为样本矩阵（Module 2 → Module 4 桥接函数）
#'
#' @param WT_matrix 矩阵。原始野生型表达矩阵。
#' @param ko_vector 数值型向量。run_virtual_knockout() 的返回值。
#' @param WT_baseline 数值型向量。计算 ko_vector 时用的基线。
#' @param method 字符串。"shift"（默认）或 "replace"。
#'
#' @return 与 WT_matrix 维度相同的 KO 表达矩阵。
#' @export
build_ko_matrix <- function(WT_matrix, ko_vector, WT_baseline, method = "shift") {

  if (!is.matrix(WT_matrix)) stop("WT_matrix must be a matrix.")
  if (ncol(WT_matrix) != length(ko_vector))
    stop(sprintf("WT_matrix cols (%d) must match ko_vector length (%d).",
         ncol(WT_matrix), length(ko_vector)))

  # Align by gene names
  common_genes <- intersect(colnames(WT_matrix), names(ko_vector))
  if (length(common_genes) == 0)
    stop("No common gene names between WT_matrix colnames and ko_vector names.")

  if (length(common_genes) < ncol(WT_matrix) * 0.9)
    warning(sprintf(
      "Gene name overlap only %.1f%%. Check data consistency.",
      length(common_genes) / ncol(WT_matrix) * 100
    ))

  # Align by WT_matrix column order
  ko_vector   <- ko_vector[colnames(WT_matrix)]
  WT_baseline <- WT_baseline[colnames(WT_matrix)]

  # Handle NAs from misaligned names
  na_ko <- is.na(ko_vector)
  if (any(na_ko)) {
    warning(sprintf("%d genes have NA ko_vector values, filling with WT baseline.", sum(na_ko)))
    ko_vector[na_ko] <- WT_baseline[na_ko]
  }

  if (method == "shift") {
    delta <- ko_vector - WT_baseline
    KO_matrix <- sweep(WT_matrix, 2, delta, FUN = "+")
    KO_matrix[KO_matrix < 0] <- 0
  } else if (method == "replace") {
    KO_matrix <- matrix(rep(ko_vector, nrow(WT_matrix)),
                        nrow = nrow(WT_matrix), byrow = TRUE,
                        dimnames = list(rownames(WT_matrix), names(ko_vector)))
  } else {
    stop("method must be 'shift' or 'replace'.")
  }
  return(KO_matrix)
}
