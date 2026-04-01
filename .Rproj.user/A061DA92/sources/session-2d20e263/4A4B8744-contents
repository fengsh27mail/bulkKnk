# ==============================================================================
# import_xena.R
# UCSC Xena 大队列极速接驳舱 (Module 0.0)
# ==============================================================================

#' 极速构建 bulkKnk 标准化分析队列 (适配 UCSC Xena)
#'
#' 本函数旨在提供“一键式”的数据加载体验。它会极速读取 UCSC Xena 下载的
#' 表达矩阵和生存数据，自动完成 Ensembl -> SYMBOL 的转换、去除死寂基因、
#' 矩阵转置，并将生存数据与表达矩阵进行绝对严格的样本对齐。
#'
#' @param expr_file 字符串。UCSC Xena 表达矩阵文件路径 (如 .tsv.gz)。
#' @param surv_file 字符串。UCSC Xena 生存数据文件路径 (如 survival.tsv)。
#' @param convert_id 逻辑值。是否自动将 Ensembl ID 转换为 SYMBOL，默认 TRUE。
#' @param species 字符串。物种，默认 "human"。用于 ID 转换。
#'
#' @return 一个 bulkKnk_cohort (S3) 对象，包含:
#'   \itemize{
#'     \item \code{expr_matrix}: 清洗完毕的表达矩阵 [行=样本, 列=基因]
#'     \item \code{clinical_data}: 完美对齐的生存数据框 (包含 OS_time, OS_status)
#'     \item \code{pheno_vector}: 默认生成的全 Tumor 表型向量
#'   }
#' @export
import_xena_cohort <- function(expr_file,
                               surv_file,
                               convert_id = TRUE,
                               species = "human") {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop(">>> 请先安装极速读取引擎: install.packages('data.table')")

  message("========================================================")
  message("🛸 bulkKnk Xena 接驳舱启动: 正在构建标准化队列...")
  message("========================================================")

  # ── 1. 极速读取表达矩阵 ───────────────────────────────────────────────
  message(">>> [1/4] 正在加载表达矩阵 (data.table 极速模式)...")
  expr_dt <- data.table::fread(expr_file, data.table = FALSE)

  # Xena 数据第一列通常是基因 ID
  gene_ids <- expr_dt[, 1]
  expr_matrix <- as.matrix(expr_dt[, -1])
  rownames(expr_matrix) <- gene_ids

  # ── 2. 智能 ID 转换 (Ensembl -> SYMBOL) ───────────────────────────────
  if (convert_id) {
    message(">>> [2/4] 正在清洗基因名并转换为 SYMBOL...")
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("需安装 clusterProfiler。")

    orgDb_name <- switch(tolower(species), "human" = "org.Hs.eg.db", "mouse" = "org.Mm.eg.db")
    if (!requireNamespace(orgDb_name, quietly = TRUE)) stop(sprintf("需安装 %s。", orgDb_name))

    # 暴力切除小数点版本号 (ENSG00000000003.14 -> ENSG00000000003)
    clean_ids <- gsub("\\..*$", "", rownames(expr_matrix))
    rownames(expr_matrix) <- clean_ids

    # 转换 ID
    suppressMessages({
      gene_trans <- clusterProfiler::bitr(clean_ids, fromType = "ENSEMBL",
                                          toType = "SYMBOL", OrgDb = get(orgDb_name))
    })

    # 去重并应用到矩阵
    gene_trans <- gene_trans[!duplicated(gene_trans$SYMBOL), ]
    expr_matrix <- expr_matrix[gene_trans$ENSEMBL, ]
    rownames(expr_matrix) <- gene_trans$SYMBOL
  }

  # ── 3. 矩阵翻转与死寂基因清理 ─────────────────────────────────────────
  message(">>> [3/4] 正在执行矩阵翻转与死寂基因过滤...")
  expr_matrix <- t(expr_matrix) # 翻转为 [行=样本, 列=基因]
  gene_vars <- apply(expr_matrix, 2, var)
  expr_matrix <- expr_matrix[, gene_vars > 0, drop = FALSE]

  # ── 4. 临床数据无缝对齐 ───────────────────────────────────────────────
  message(">>> [4/4] 正在解析临床生存数据并执行严格样本对齐...")
  surv_df <- read.table(surv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # 寻找完全匹配的样本交集
  common_samples <- intersect(rownames(expr_matrix), surv_df$sample)
  if (length(common_samples) < 10) {
    stop("🚨 严重错误: 表达矩阵与生存数据的样本名匹配数量不足 10 个，请检查文件是否匹配！")
  }

  # 严格切片与对齐
  expr_matrix <- expr_matrix[common_samples, ]
  surv_df <- surv_df[match(common_samples, surv_df$sample), ]

  # 构建符合 bulkKnk 层4 要求的标准化 data.frame
  clinical_data <- data.frame(
    sample_id = surv_df$sample,
    OS_time   = surv_df$OS.time / 30, # Xena 默认为天，除以 30 转为月
    OS_status = surv_df$OS,           # 1=死亡, 0=存活
    stringsAsFactors = FALSE
  )

  # 默认生成全 Tumor 表型（Xena 主队列通常是全肿瘤）
  pheno_vector <- rep("Tumor", nrow(expr_matrix))

  message(sprintf(">>> [接驳完成] 队列构建成功！有效样本: %d, 有效基因: %d",
                  nrow(expr_matrix), ncol(expr_matrix)))

  # 封装为 S3 对象返回
  res <- list(
    expr_matrix   = expr_matrix,
    clinical_data = clinical_data,
    pheno_vector  = pheno_vector,
    target_pheno  = "Tumor"
  )
  class(res) <- "bulkKnk_cohort"
  return(res)
}
