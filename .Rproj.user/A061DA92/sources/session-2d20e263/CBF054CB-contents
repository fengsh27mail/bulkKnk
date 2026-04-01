#' 利用 ARACNE 逻辑净化调控网络 (Module 1C)
#'
#' 该函数应用数据处理不等式 (DPI) 过滤网络中的间接调控边。
#' 对于每一组基因三元组 (A, B, C)，它会比较三条边的权重，并移除最弱的一条。
#'
#' @param adj_matrix 矩阵。由 infer_causal_network 生成的基因调控权重矩阵。
#' @param eps 数值。DPI 容忍阈值 (0-1)，默认 0.1。较大的 eps 会保留更多边。
#'
#' @return 经过 DPI 过滤后的净化矩阵。
#' @export
refine_network_aracne <- function(adj_matrix, eps = 0.1) {

  message(">>> [ARACNE] 正在应用 DPI 不等式执行网络净化...")

  refined_mat <- adj_matrix
  genes <- rownames(adj_matrix)
  n <- length(genes)

  # 为了计算效率，仅处理非零权重的边
  # 遍历所有可能的三元组 (i, j, k)
  for (i in 1:(n - 2)) {
    for (j in (i + 1):(n - 1)) {
      # 检查边 (i, j) 是否存在
      if (adj_matrix[i, j] == 0) next

      for (k in (j + 1):n) {
        # 获取三元组的三条边权重
        w_ij <- adj_matrix[i, j]
        w_jk <- adj_matrix[j, k]
        w_ik <- adj_matrix[i, k]

        # 如果三条边不构成闭环（存在零边），则跳过
        if (w_jk == 0 || w_ik == 0) next

        # 识别三元组中最弱的边
        min_w <- min(w_ij, w_jk, w_ik)

        # 应用 DPI 规则：剔除权重远小于其他两者的那条边
        if (min_w == w_ik && w_ik < min(w_ij, w_jk) * (1 - eps)) {
          refined_mat[i, k] <- refined_mat[k, i] <- 0
        } else if (min_w == w_ij && w_ij < min(w_jk, w_ik) * (1 - eps)) {
          refined_mat[i, j] <- refined_mat[j, i] <- 0
        } else if (min_w == w_jk && w_jk < min(w_ij, w_ik) * (1 - eps)) {
          refined_mat[j, k] <- refined_mat[k, j] <- 0
        }
      }
    }
  }

  removed_count <- sum(adj_matrix > 0) - sum(refined_mat > 0)
  message(paste(">>> [ARACNE] 净化完成！共剔除", removed_count, "条间接调控边。"))

  return(refined_mat)
}
