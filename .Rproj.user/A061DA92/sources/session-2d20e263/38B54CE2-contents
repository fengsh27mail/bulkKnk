# 全局变量声明，防止 R CMD check 报 note
utils::globalVariables(c("log2FC", "pvalue", "Status", "Log2FC", "WT", "moduleColors", "Gene", "KO", "Log2FC_Virtual"))

#' 绘制真实与虚拟对照的双重火山图 (Dual Volcano Plot)
#'
#' @param real_deg_df 数据框或 NULL。包含真实差异分析结果，列名须含 "Gene", "log2FC", "pvalue"。若为 NULL，则仅绘制虚拟敲除图。
#' @param WT_matrix 矩阵或向量。野生型表达矩阵 (行=样本, 列=基因) 或平均表达向量。
#' @param KO_matrix 矩阵或向量。虚拟敲除后稳态矩阵或平均表达向量。
#' @param target_genes 字符向量。被敲除的靶点名称。
#' @param p_cutoff 数值。真实火山图的 P 值阈值，默认 0.05。
#' @param fc_cutoff 数值。真实火山图的 Log2FC 阈值，默认 1.5。
#'
#' @return 返回一个由 patchwork 拼接好的 ggplot2 对比图对象。
#' @import ggplot2
#' @importFrom patchwork plot_layout
#' @export
plot_dual_volcano <- function(real_deg_df = NULL, WT_matrix, KO_matrix, target_genes, p_cutoff = 0.05, fc_cutoff = 1.5) {

  if (!requireNamespace("patchwork", quietly = TRUE)) stop("请先安装 patchwork 包: install.packages('patchwork')")

  # 1. 智能处理输入维度 (支持向量或矩阵)
  WT_expr <- if(is.matrix(WT_matrix) || is.data.frame(WT_matrix)) colMeans(WT_matrix) else WT_matrix
  KO_expr <- if(is.matrix(KO_matrix) || is.data.frame(KO_matrix)) colMeans(KO_matrix) else KO_matrix

  if(length(WT_expr) != length(KO_expr)) stop("WT 和 KO 的基因数量不一致！")

  # 2. 图 B: 虚拟敲除火山图 (散点图: 变化倍数 vs 基线表达)
  virt_df <- data.frame(Gene = names(WT_expr), WT = WT_expr, KO = KO_expr)
  # 添加伪计数避免 log2(0)
  virt_df$Log2FC_Virtual <- log2((virt_df$KO + 0.1) / (virt_df$WT + 0.1))

  virt_df$Status <- "Stable"
  virt_df$Status[virt_df$Log2FC_Virtual < -0.1] <- "Downregulated"
  virt_df$Status[virt_df$Log2FC_Virtual > 0.1]  <- "Upregulated"
  virt_df$Status[virt_df$Gene %in% target_genes] <- "Target(s)"

  p_virt <- ggplot(virt_df, aes(x = Log2FC_Virtual, y = WT, color = Status)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Stable" = "grey80",
                                  "Downregulated" = "#2980b9",
                                  "Upregulated" = "#c0392b",
                                  "Target(s)" = "#8e44ad")) +
    theme_minimal(base_size = 12) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = "Virtual Knockout Effect",
         subtitle = paste("Perturbed Targets:", paste(target_genes, collapse = ", ")),
         x = "Virtual Log2FC (KO vs WT)",
         y = "Baseline Expression (WT)")

  # 3. 灵活处理：如果没有真实 DEG 数据，只返回虚拟图
  if (is.null(real_deg_df)) {
    return(p_virt)
  }

  # 4. 图 A: 真实临床差异基因
  real_deg_df$Status <- "Not Sig"
  real_deg_df$Status[real_deg_df$log2FC >= fc_cutoff & real_deg_df$pvalue <= p_cutoff] <- "Up"
  real_deg_df$Status[real_deg_df$log2FC <= -fc_cutoff & real_deg_df$pvalue <= p_cutoff] <- "Down"

  p_real <- ggplot(real_deg_df, aes(x = log2FC, y = -log10(pvalue), color = Status)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "#c0392b", "Down" = "#2980b9", "Not Sig" = "grey80")) +
    theme_minimal(base_size = 12) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    labs(title = "Real Clinical DEGs", x = "Real Log2FC", y = "-Log10(P-value)") +
    theme(legend.position = "none")

  # 终极拼图
  return(p_real + p_virt + patchwork::plot_layout(widths = c(1, 1.2)))
}

#' 运行虚拟敲除后的 GO 富集分析 (全物种自适应版)
#'
#' @param virtual_degs 字符向量。虚拟敲除后显著下调/上调的基因名 (通常为 SYMBOL 格式)。
#' @param species 字符串。快捷物种选项："human", "mouse", "rat"。
#' @param custom_orgDb 字符串。若物种不在快捷选项内，可传入物种专用的注释包名 (例如 "org.At.tair.db" 拟南芥)。
#'
#' @return clusterProfiler 的 enrichResult 对象
#' @export
run_vk_enrichment <- function(virtual_degs, species = "human", custom_orgDb = NULL, database = "GO") {

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) stop("请安装 clusterProfiler 包。")

  # 1. 确定使用的本地注释数据库 (用于 ID 转换)
  if (!is.null(custom_orgDb)) {
    org_db_name <- custom_orgDb
  } else {
    org_db_name <- switch(tolower(species),
                          "human" = "org.Hs.eg.db",
                          "mouse" = "org.Mm.eg.db",
                          "rat"   = "org.Rn.eg.db",
                          stop("快捷物种仅支持 human, mouse, rat。其它物种请使用 custom_orgDb 参数。"))
  }

  if (!requireNamespace(org_db_name, quietly = TRUE)) {
    stop(sprintf("请先安装注释包: BiocManager::install('%s')", org_db_name))
  }

  # 2. 基因 ID 转换 (SYMBOL -> ENTREZID)
  gene_ids <- tryCatch({
    clusterProfiler::bitr(virtual_degs, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db_name)
  }, error = function(e) {
    stop("基因 ID 转换失败。请确认输入的基因名格式是否为 SYMBOL。")
  })

  if(nrow(gene_ids) == 0) stop("成功转换的基因 ID 数量为 0，停止富集分析。")
  message(sprintf(">> [Enrichment] 成功转换 %d 个基因，准备进行 %s 分析...", nrow(gene_ids), toupper(database)))

  # 3. 双引擎分流：GO 或 KEGG
  if (toupper(database) == "GO") {
    res <- clusterProfiler::enrichGO(
      gene          = gene_ids$ENTREZID,
      OrgDb         = org_db_name,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE  # 自动将 ENTREZID 转换回可读的 SYMBOL
    )

  } else if (toupper(database) == "KEGG") {
    # KEGG 需要特定的物种缩写代号
    kegg_org <- switch(tolower(species),
                       "human" = "hsa",
                       "mouse" = "mmu",
                       "rat"   = "rno",
                       "hsa") # 默认给 hsa

    message(paste(">> [Enrichment] 正在联网获取最新 KEGG 数据，物种代号:", kegg_org))
    res <- clusterProfiler::enrichKEGG(
      gene          = gene_ids$ENTREZID,
      organism      = kegg_org,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )

    # 将 KEGG 结果中的 ENTREZID 转换回 SYMBOL，方便用户画图和阅读
    if (!is.null(res) && nrow(res) > 0) {
      res <- clusterProfiler::setReadable(res, OrgDb = org_db_name, keyType="ENTREZID")
    }

  } else {
    stop("database 参数错误！仅支持 'GO' 或 'KEGG'。")
  }

  # 4. 打印报告
  if (is.null(res) || nrow(res) == 0) {
    message(sprintf(">> [Enrichment] 提示：未能富集到显著的 %s 通路。", toupper(database)))
  } else {
    message(sprintf(">> [Enrichment] 完成！发现 %d 条显著富集的 %s 通路。", nrow(res), toupper(database)))
  }

  return(res)
}
