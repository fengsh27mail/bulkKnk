#' 智能特征筛选：保送核心签名基因或自适应高变基因 (Module 0.1)
#'
#' 在进入网络推断前，自动确保 AI 预测模块所需的关键生物学签名基因被保留。
#' 若检测到输入数据与内置人类签名极度不匹配（如非人类物种、小样本靶向测序），
#' 将自动降级为纯数据驱动的高变基因 (HVG) 模式，确保跨物种与跨平台的绝对通用性。
#'
#' @param expr_matrix 矩阵。原始表达矩阵 (行=样本, 列=基因)。
#' @param n_features 整数。最终保留的总基因数，默认 500。
#' @param fallback_threshold 整数。当匹配到的签名基因少于此阈值时，触发通用模式，默认 10。
#'
#' @return 筛选后的表达矩阵。
#' @export
select_advanced_features <- function(expr_matrix, n_features = 500, fallback_threshold = 10) {

  # 1. 安全获取 predict_ko_phenotype 中定义的所有内置签名基因
  sig_genes <- tryCatch({
    unique(unlist(.get_builtin_signatures()))
  }, error = function(e) {
    # 【修复】明确区分"函数未加载"和"其他错误"，避免静默误导用户
    if (grepl("could not find function", conditionMessage(e))) {
      stop("找不到内置签名函数 .get_builtin_signatures()！\n请确认已按顺序加载：source('R/ai_predictor.R')")
    }
    character(0)
  })

  # 2. 识别用户数据中存在的签名基因
  present_sigs <- intersect(sig_genes, colnames(expr_matrix))

  # --- 核心通用化改造逻辑 ---
  if (length(present_sigs) < fallback_threshold) {
    message(sprintf(">>> [通用模式] 仅匹配到 %d 个内置签名基因，低于阈值 (%d)。",
                    length(present_sigs), fallback_threshold))
    message(">>> [通用模式] 推测为非人类物种或靶向测序数据，已自动切换为'高变基因 (HVG) 驱动模式'。")
    message(">>> 提示：在运行 predict_ko_phenotype 时，强烈建议使用 custom_signatures 参数传入您的物种专属基因集。")

    # 清空 present_sigs，彻底交由数据的自然方差 (SD) 来决定保留哪些基因
    present_sigs <- character(0)
  } else {
    message(paste(">>> [核心保送] 完美匹配人类背景，已自动锁定", length(present_sigs), "个关键生物学签名基因。"))
  }
  # --------------------------

  # 3. 计算所有基因的变异度 (标准差 SD)
  gene_sds <- apply(expr_matrix, 2, stats::sd)

  # 4. 排除已保送的基因，从剩余基因中按 SD 降序补齐
  remaining_genes <- setdiff(colnames(expr_matrix), present_sigs)
  needed_count <- n_features - length(present_sigs)

  if (needed_count > 0) {
    other_sds <- gene_sds[remaining_genes]
    top_hvgs <- names(sort(other_sds, decreasing = TRUE))[1:min(needed_count, length(other_sds))]
    final_genes <- c(present_sigs, top_hvgs)
  } else {
    final_genes <- present_sigs[1:n_features]
  }

  message(paste(">>> [特征筛选] 矩阵构建完成，总特征数:", length(final_genes)))

  return(expr_matrix[, final_genes, drop = FALSE])
}
