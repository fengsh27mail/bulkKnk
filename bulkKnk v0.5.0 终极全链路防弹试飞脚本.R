# ==============================================================================
# bulkKnk v0.5.0 终极全链路防弹试飞脚本 (run_test.R)
# ==============================================================================

# 0. 加载你的包和数据
devtools::load_all()
data("tcga_stad_mini")

pheno_vector <- tcga_stad_mini$pheno_data$SampleType
target_group <- "Tumor"

message("\n========================================================")
message("🚀 启动 bulkKnk 全链路试飞")
message("========================================================\n")

# ------------------------------------------------------------------------------
# Module 0 & 1: 特征提取与靶点发现
# ------------------------------------------------------------------------------
# 1. 提炼 500 个核心特征
sub_expr <- select_advanced_features(tcga_stad_mini$expr_matrix, n_features = 500)

# 2. WGCNA 寻找与肿瘤最相关的基因模块
trait_num <- as.numeric(as.factor(pheno_vector))
wgcna_res <- identify_hub_modules(sub_expr, trait_vec = trait_num)

# 3. 提取最多 100 个核心基因防止网络推断卡死
top_100_hubs <- head(wgcna_res$best_genes, 100)
top_target <- top_100_hubs[1] # 取第一个作为我们的打击靶点！
message(sprintf(">>> [决策] 我们即将敲除的超级靶点是: %s", top_target))

# ------------------------------------------------------------------------------
# Module 2: 构建因果网络
# ------------------------------------------------------------------------------
# 4. 提取这 100 个基因的矩阵 (强制保持矩阵格式)
focus_mat <- as.matrix(sub_expr[, top_100_hubs, drop = FALSE])

# 5. GENIE3 推断 (你包里的 infer_causal_network 内部会自动转置，这里直接喂)
net_raw <- infer_causal_network(focus_mat)

# 6. DPI 净化
net_clean <- prune_network_dpi_fast(net_raw$clean_matrix, dpi_threshold = 0.5)

# 🚨 终极防弹装甲：强制为网络矩阵贴回基因名，防止 GENIE3 吞名字！
rownames(net_clean) <- colnames(focus_mat)
colnames(net_clean) <- colnames(focus_mat)

# 7. 网络归一化
W_norm <- sweep(t(net_clean), 2, colSums(t(net_clean)) + 1e-9, FUN="/")

# ------------------------------------------------------------------------------
# Module 3: 虚拟敲除模拟
# ------------------------------------------------------------------------------
# 8. 执行批量干预
ko_results <- batch_virtual_knockout(
  expr_matrix  = focus_mat,
  pheno_vec    = pheno_vector,
  target_pheno = target_group,
  target_genes = top_target,
  adj_matrix   = W_norm,
  method       = "shift"
)

# ------------------------------------------------------------------------------
# Module 4: AI 表型逆转预测
# ------------------------------------------------------------------------------
# 9. 设定基线标签 (高风险 vs 低风险)
baseline_labels <- ifelse(ko_results$WT[, top_target] > median(ko_results$WT[, top_target]),
                          "HighRisk", "LowRisk")

# 10. 运行终极预测！(🚨注意：AI预测要求行=基因，列=样本，所以这里必须用 t() 转置)
report <- predict_ko_phenotype(
  WT_matrix    = t(ko_results$WT),
  KO_matrix    = t(ko_results$KO),
  pheno_labels = baseline_labels,
  target_pheno = "HighRisk",
  ko_gene      = top_target
)

# ------------------------------------------------------------------------------
# 终点：出图
# ------------------------------------------------------------------------------
p <- plot_ko_delta(report, color_by = "direction")
print(p)




# ==============================================================================
# 附加环节：下游高级分析 (火山图 & GO 富集)
# 依赖模块: downstream_eval.R
# ==============================================================================

# 1. 快速计算真实的临床差异 (Tumor vs Normal) 作为对照
message("   > 正在快速计算真实临床队列的 T 检验 P-value 与 Log2FC...")
is_tumor <- pheno_vector == "Tumor"
expr_tumor <- sub_expr[is_tumor, , drop = FALSE]
expr_normal <- sub_expr[!is_tumor, , drop = FALSE]

mean_tumor <- colMeans(expr_tumor)
mean_normal <- colMeans(expr_normal)

# 算真实 Log2FC
real_log2fc <- log2((mean_tumor + 0.1) / (mean_normal + 0.1))

# 批量算真实 P-value
real_pvals <- sapply(seq_len(ncol(sub_expr)), function(i) {
  # 加上 tryCatch 防止某些基因方差为 0 导致 t.test 报错
  tryCatch(t.test(expr_tumor[, i], expr_normal[, i])$p.value, error = function(e) 1)
})

# 组装出你的函数需要的 real_deg_df
my_real_deg <- data.frame(
  Gene   = colnames(sub_expr),
  log2FC = real_log2fc,
  pvalue = real_pvals
)

# 2. 召唤完全体：双重火山图！
p2_dual <- plot_dual_volcano(
  real_deg_df  = my_real_deg,   # 这次我们把真实的差异数据传进去了！
  WT_matrix    = ko_results$WT,
  KO_matrix    = ko_results$KO,
  target_genes = top_target
)

print(p2_dual)


message("\n--- Step 8: 运行级联反应 GO 富集分析 ---")
# 1. 提取 WT 和 KO 的群体平均基线
wt_mean <- colMeans(ko_results$WT)
ko_mean <- colMeans(ko_results$KO)

# 2. 计算虚拟 Log2FC 变化倍数 (加 0.1 伪计数防止 log2(0))
log2fc <- log2((ko_mean + 0.1) / (wt_mean + 0.1))

# 3. 提取发生显著级联变化的基因（比如 Log2FC 绝对值 > 0.05 的节点）
# 注意：我们要把靶点自己剔除，只看它引起了哪些别的基因改变
virtual_degs <- names(log2fc)[abs(log2fc) > 0.05]
virtual_degs <- setdiff(virtual_degs, top_target)

if (length(virtual_degs) > 0) {
  message(sprintf(">>> [Enrichment] 发现 %d 个受级联影响发生偏转的基因，启动 GO 富集...", length(virtual_degs)))

  # 调用你写好的全物种自适应富集函数！(这里演示用 human)
  go_res <- run_vk_enrichment(virtual_degs = virtual_degs, species = "human")

  # 如果富集出了显著通路，直接画一张高大上的气泡图！
  if (!is.null(go_res) && nrow(go_res) > 0) {
    if (requireNamespace("clusterProfiler", quietly = TRUE)) {
      p3 <- clusterProfiler::dotplot(go_res, showCategory = 10,
                                     title = paste("Pathways Shifted after", top_target, "KO"))
      print(p3)
    }
  }
} else {
  message(">>> [Enrichment] 提示：该靶点敲除未在局部网络内引起广泛级联反应，跳过富集分析。")
}


# 在提取出 virtual_degs 后...

# 1. 跑 GO 并画图
go_res <- run_vk_enrichment(virtual_degs, species = "human", database = "GO")
if (!is.null(go_res) && nrow(go_res) > 0) {
  print(clusterProfiler::dotplot(go_res, title = "GO: Biological Processes"))
}

# 2. 跑 KEGG 并画图
kegg_res <- run_vk_enrichment(virtual_degs, species = "human", database = "KEGG")
if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
  print(clusterProfiler::dotplot(kegg_res, title = "KEGG: Signaling Pathways"))
}
