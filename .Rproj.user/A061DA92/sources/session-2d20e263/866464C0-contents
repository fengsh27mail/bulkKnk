# bulkKnk (v0.5.0) 🧬

**A Universal, Data-Driven Virtual Perturbation & Phenotype Reversal Platform for Bulk RNA-seq**

`bulkKnk` 是一个面向**全物种、全疾病模型**的系统生物学工具箱。
它将加权基因共表达网络 (WGCNA)、因果网络推断 (GENIE3 + ARACNE) 与机器学习 (Random Forest) 深度融合，
旨在从常规的 Bulk RNA-seq 数据中：
1. **数据驱动地**发现核心枢纽基因 (Hub Genes)。
2. **在硅 (in silico) 模拟**对这些靶点进行基因敲除/干预后的全网转录组稳态漂移。
3. **多维预测**干预操作带来的表型逆转概率与通路活性改变。

无论是人类肿瘤微环境重塑、小鼠发育阻滞挽救，还是植物抗逆性通路的调控，
`bulkKnk` 都能提供从“靶点发现”到“疗效预演”的端到端解决方案。

---

## 🌟 核心特性 (Core Features)

- **🌍 物种与疾病无关 (Species & Disease Agnostic)**：完全基于数据驱动的无监督/半监督算法，支持自定义任意物种的生物学通路签名 (Signatures)。
- **🛡️ 智能特征锁定 (Smart Feature Selection)**：自动识别并“保送”关键通路基因，辅以高变基因 (HVGs) 补齐，确保信号捕捉不遗漏。
- **🕸️ 高精度因果网络 (ARACNE-Pruned Network)**：利用数据处理不等式 (DPI) 净化 GENIE3 产生的间接调控噪音，构建精准的有向物理网络。
- **💥 虚拟扰动引擎 (Virtual RWR Engine)**：通过改进的带重启随机游走算法，模拟单靶点或多靶点联合物理敲除后的级联反应。
- **🤖 AI 预言家 (AI Phenotype Predictor)**：内置随机森林分类器，量化评估敲除后的样本是否成功向“健康/敏感/低风险”表型发生逆转。

---

## 📦 安装 (Installation)

你可以通过 R 控制台安装开发版本：

```R
```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("fengsh27mail/bulkKnk") # 请替换为实际仓库地址

#快速开始 (Quick Start)
#以下展示如何使用 bulkKnk 跑通一个标准的**“表型干预预演”**工作流。

#(注：示例使用内置的 tcga_stad_mini 数据集，但该流程完美适用于任何包含 expr_matrix 和表型分组的转录组数据。)

library(bulkKnk)#加载包
data("tcga_stad_mini")#加载示例数据

# 0. 数据准备：定义你的“目标干预群” (例如: Tumor, Mutant, Treated)
target_group <- "Tumor" 
pheno_vector <- tcga_stad_mini$pheno_data$SampleType

# 1. 智能特征构建 (自动保留关键通路基因 + 高变基因)
sub_expr <- select_advanced_features(tcga_stad_mini$expr_matrix, n_features = 500)

# 2. 靶点发现：识别与目标表型最相关的 WGCNA 模块
trait_num <- as.numeric(as.factor(pheno_vector))
wgcna_res <- identify_hub_modules(sub_expr, trait_vec = trait_num)
focus_genes <- wgcna_res$best_genes
top_target <- wgcna_res$hub_genes[1] # 自动提取网络核心度最高的基因作为敲除靶点

# 3. 构建并净化局部因果网络
net_raw <- infer_causal_network(sub_expr[, focus_genes])
net_clean <- prune_network_dpi_fast(net_raw$clean_matrix) # DPI 净化去噪
W_norm <- sweep(t(net_clean), 2, colSums(t(net_clean)) + 1e-9, FUN="/")

# 4. 执行群体级虚拟敲除 (仅对目标群体施加干预)
ko_results <- batch_virtual_knockout(
  expr_matrix  = sub_expr[, focus_genes], 
  pheno_vec    = pheno_vector,
  target_pheno = target_group, 
  target_genes = top_target,  # 敲除 WGCNA 算出的 Top Hub
  adj_matrix   = W_norm
)

# 5. AI 表型逆转评估
# 使用核心基因的中位数，将目标群体划分为 "High_Risk" 和 "Low_Risk" 基线
baseline_labels <- ifelse(ko_results$WT[, top_target] > median(ko_results$WT[, top_target]), 
                          "High_Risk", "Low_Risk")

# 支持传入用户自定义的物种专属基因集 (可选)
my_custom_pathways <- list(
  My_Metabolism = c("GeneA", "GeneB", "GeneC"),
  My_Stress_Response = c("GeneX", "GeneY", "GeneZ")
)

report <- predict_ko_phenotype(
  WT_matrix         = t(ko_results$WT), 
  KO_matrix         = t(ko_results$KO), 
  pheno_labels      = baseline_labels, 
  target_pheno      = "High_Risk",
  ko_gene           = top_target,
  custom_signatures = my_custom_pathways # 非人类物种强烈建议传入此参数
)

# 6. 可视化干预带来的全网特征改变
plot_ko_delta(report, color_by = "direction")


#进阶用法：如何适配非人类物种数据？
#bulkKnk 的网络推断和敲除引擎是完全数学化的，无需任何修改。
#对于步骤 5 中的 predict_ko_phenotype，
#由于内置签名库 (ssGSEA Signatures) 主要是人类 HGNC 基因符号，
#如果是小鼠、斑马鱼或拟南芥数据，请务必通过 custom_signatures 参数传入你的自定义通路列表
#（例如从 KEGG/GO 导出的目标物种基因集），AI 引擎将自动使用你的自定义通路进行机器学习分类和预测。
#引用 (Citation)
#如果你在科学研究中使用了 bulkKnk，请引用：

#Bioinformatics Night-Shift Hero (shuofeng MUST). (2026). bulkKnk: A Data-Driven Virtual Perturbation Platform for Universal Bulk RNA-seq. 
GitHub：fengsh27mail
