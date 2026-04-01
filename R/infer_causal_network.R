#' 推断因果调控网络并进行动态去噪
#'
#' 本函数利用 GENIE3 算法推断基因间的有向调控关系，并自动基于数据分布特征
#' （分位数法或均值标准差法）动态设定硬阈值，过滤背景噪音，输出干净的邻接矩阵。
#'
#' @param expr_matrix 表达量矩阵。格式要求：行是样本 (Samples)，列是基因 (Genes)。
#' @param threshold_method 动态阈值策略，可选 "top_quantile" (默认) 或 "mean_sd"。
#' @param q_val 当使用 "top_quantile" 时，设定的分位数阈值，默认 0.95 (即 Top 5%)。
#' @param sd_multiplier 当使用 "mean_sd" 时，标准差的倍数，默认 3。
#' @param n_cores 并行计算的线程数，默认 1 (如果是庞大网络建议设置更大)。
#'
#' @return 返回一个 List，包含清洗后的权重矩阵、阈值信息和保留的连线数。
#' @export
#'
#' @examples
#' # 假设你已经有了 simulated_data
#' # network_res <- infer_causal_network(simulated_data, threshold_method = "top_quantile")
infer_causal_network <- function(expr_matrix,
                                 threshold_method = "top_quantile",
                                 q_val = 0.95,
                                 sd_multiplier = 3,
                                 n_cores = 1) {

  # 1. 输入检查与格式转换 (GENIE3 要求行是基因，列是样本)
  if(is.null(colnames(expr_matrix))) stop("表达矩阵必须包含基因名 (colnames)！")
  message(">> 正在启动 GENIE3 核心推断算法，请稍候...")

  # 强制转置以符合 GENIE3 输入要求
  genie3_input <- t(as.matrix(expr_matrix))

  # 2. 运行 GENIE3
  # 注意：在描述文件 DESCRIPTION 里声明了 Imports: GENIE3 后，
  # 这里建议使用 GENIE3:: 前缀调用，更为严谨。
  raw_weight_matrix <- GENIE3::GENIE3(genie3_input, nCores = n_cores)

  # 3. 提取有效权重并计算动态阈值
  all_weights <- as.vector(raw_weight_matrix)
  valid_weights <- all_weights[all_weights > 0]

  # 处理如果没有有效权重的情况
  if(length(valid_weights) == 0) stop("GENIE3 未能推断出任何有效连线。")

  if (threshold_method == "top_quantile") {
    final_threshold <- quantile(valid_weights, q_val)
    message(paste(">> 采用 Top", round((1-q_val)*100, 2), "% 策略，计算得出动态阈值:", round(final_threshold, 4)))
  } else if (threshold_method == "mean_sd") {
    final_threshold <- mean(valid_weights) + sd_multiplier * sd(valid_weights)
    message(paste(">> 采用 Mean +", sd_multiplier, "SD 策略，计算得出动态阈值:", round(final_threshold, 4)))
  } else {
    stop("threshold_method 参数错误！请选择 'top_quantile' 或 'mean_sd'。")
  }

  # 4. 矩阵清洗去噪
  clean_weight_matrix <- raw_weight_matrix
  clean_weight_matrix[clean_weight_matrix < final_threshold] <- 0

  # 5. 统计网络拓扑信息
  total_edges <- length(valid_weights)
  retained_edges <- sum(clean_weight_matrix > 0)
  message(paste(">> 网络修剪完成！保留了", retained_edges, "条有效连线，剔除了", total_edges - retained_edges, "条噪音连线。"))

  # 返回结果对象
  return(list(
    clean_matrix = clean_weight_matrix,
    threshold_used = final_threshold,
    edge_stats = c(Total = total_edges, Retained = retained_edges)
  ))
}
