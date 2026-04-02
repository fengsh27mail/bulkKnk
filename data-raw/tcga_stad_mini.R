# =========================================================================
# 脚本：生成 bulkKnk 的内置真实测试数据 (TCGA-STAD subset)
# =========================================================================
library(TCGAbiolinks)
library(SummarizedExperiment)

message(">>> 1. 正在向 GDC 提交 TCGA-STAD (胃癌) RNA-seq 查询请求...")
query <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

all_barcodes <- getResults(query, cols = c("cases"))
tumor_cases <- head(all_barcodes[grepl("-01A", all_barcodes)], 50)
normal_cases <- head(all_barcodes[grepl("-11A", all_barcodes)], 10)

query_mini <- GDCquery(
  project = "TCGA-STAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = c(tumor_cases, normal_cases)
)

message(">>> 2. 正在下载表达数据 (这可能需要几分钟，请保持网络畅通)...")
GDCdownload(query_mini)

message(">>> 3. 正在组装表达矩阵与临床信息...")
tcga_data <- GDCprepare(query_mini)

tpm_matrix <- assay(tcga_data, "tpm_unstrand")
gene_info <- rowData(tcga_data)
rownames(tpm_matrix) <- gene_info$gene_name

clinical_info <- colData(tcga_data)

message(">>> 4. 数据降维与浓缩 (提纯核心网络基因)...")
tpm_matrix <- tpm_matrix[rowSums(tpm_matrix > 1) > (ncol(tpm_matrix) * 0.2), ]

mads <- apply(tpm_matrix, 1, mad)
top_var_genes <- names(sort(mads, decreasing = TRUE)[1:500])

immune_genes <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CXCL9", "CXCL10", "STAT3", "MYC")
final_genes <- unique(c(immune_genes, top_var_genes))
final_genes <- intersect(final_genes, rownames(tpm_matrix))

mini_matrix <- tpm_matrix[final_genes, ]
tcga_stad_expr <- t(mini_matrix)

tcga_stad_pheno <- data.frame(
  SampleID = rownames(tcga_stad_expr),
  SampleType = ifelse(grepl("-11A", rownames(tcga_stad_expr)), "Normal", "Tumor"),
  stringsAsFactors = FALSE
)

if ("tumor_stage" %in% colnames(clinical_info)) {
  tcga_stad_pheno$Stage <- clinical_info$tumor_stage
} else if ("ajcc_pathologic_stage" %in% colnames(clinical_info)) {
  tcga_stad_pheno$Stage <- clinical_info$ajcc_pathologic_stage
} else {
  tcga_stad_pheno$Stage <- "Unknown"
}

if ("vital_status" %in% colnames(clinical_info)) {
  tcga_stad_pheno$VitalStatus <- clinical_info$vital_status
} else {
  tcga_stad_pheno$VitalStatus <- "Unknown"
}

rownames(tcga_stad_pheno) <- tcga_stad_pheno$SampleID

tcga_stad_mini <- list(
  expr_matrix = tcga_stad_expr,
  pheno_data = tcga_stad_pheno
)

message(">>> 5. 正在将真实的 TCGA 数据打包进 bulkKnk 引擎中...")
usethis::use_data(tcga_stad_mini, overwrite = TRUE)

message(">>> 完美！测试数据集构建完毕。")
