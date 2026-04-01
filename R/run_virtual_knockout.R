#' 执行全网络虚拟基因或通路干预 (核心网络扰动引擎)
#'
#' @param target_genes 字符向量。支持单个或多个靶点 (如 c("GeneA", "GeneB"))。
#' @param E_init 数值型向量。全网基因的基线表达量 (通常为目标群体的均值)。
#' @param adj_matrix 矩阵。经过清洗和归一化的有向调控网络。
#' @param alpha 数值。网络信息保留系数 (restart probability)，默认 0.7。
#' @param steps 整数。迭代寻找稳态的次数，默认 30 次。
#' @param check_convergence 逻辑值。是否打印每步收敛误差，默认 FALSE。
#'
#' @return 返回一个数值型向量，代表虚拟干预后的全网稳态绝对表达量。
#' @export
run_virtual_knockout <- function(target_genes = NULL,
                                 E_init,
                                 adj_matrix,
                                 alpha = 0.7,
                                 steps = 30,
                                 check_convergence = FALSE) {

  if (!is.matrix(adj_matrix)) stop("adj_matrix 必须是一个标准的矩阵。")
  if (length(E_init) != nrow(adj_matrix)) stop("E_init 的长度必须与 adj_matrix 的行数相等！")

  E_anchor <- E_init   # 锚点：RWR 重启时始终回归此稳态
  E_current <- E_init
  W_current <- adj_matrix

  if (!is.null(target_genes)) {
    valid_targets <- intersect(target_genes, names(E_init))

    if (length(valid_targets) == 0) {
      stop("你输入的靶基因都不在这个网络中！请检查基因名。")
    }

    message(sprintf(">>> [核心引擎] 执行拓扑级虚拟干预，涉及 %d 个有效靶点...", length(valid_targets)))

    # 物理阻断：将靶基因的初始表达和网络输入/输出边全部清零
    E_anchor[valid_targets] <- 0
    E_current[valid_targets] <- 0
    W_current[, valid_targets] <- 0   # 切断"其他基因→靶点"的输入边
    W_current[valid_targets, ] <- 0   # 切断"靶点→其他基因"的输出边
  }

  # RWR 迭代求稳态
  for (i in seq_len(steps)) {
    E_new <- alpha * (W_current %*% E_current) + (1 - alpha) * E_anchor
    E_new <- as.vector(E_new)
    names(E_new) <- names(E_init)

    if (check_convergence) {
      conv_err <- max(abs(E_new - E_current))
      message(sprintf("   迭代 %02d/%d | 最大收敛误差: %.2e", i, steps, conv_err))
    }
    E_current <- E_new
  }

  return(E_current)
}


#' 队列级批量虚拟干预引擎 (Batch Virtual Perturbation)
#'
#' 自动提取目标群体，结合 run_virtual_knockout 与 build_ko_matrix，
#' 高效生成用于 AI 预测模块的多样本干预矩阵。
#'
#' @param expr_matrix 矩阵。原始表达矩阵 (行=样本, 列=基因)。
#' @param pheno_vec 向量。样本表型标签，长度须等于 nrow(expr_matrix)。
#' @param target_pheno 字符串。指定要进行干预的样本类型 (例如 "Tumor", "Mutant", "Disease")。
#' @param target_genes 字符向量。干预靶点名称。
#' @param adj_matrix 矩阵。由 Module 1 生成并净化的网络。
#' @param method 字符串。桥接策略，"shift" (保留个体差异，推荐) 或 "replace" (绝对替换)。
#'
#' @return 包含 WT 矩阵和对应 KO 矩阵的列表。
#' @export
batch_virtual_knockout <- function(expr_matrix, pheno_vec, target_pheno, target_genes, adj_matrix, method = "shift") {

  # 1. 提取目标群体的 WT 矩阵
  target_idx <- which(pheno_vec == target_pheno)
  if(length(target_idx) == 0) stop(sprintf("未在 pheno_vec 中找到表型为 '%s' 的样本！", target_pheno))

  WT_mat <- expr_matrix[target_idx, , drop = FALSE]
  message(sprintf(">>> [Batch Perturbation] 锁定 %d 个 [%s] 样本进行模拟干预...", length(target_idx), target_pheno))

  # 2. 计算目标群体的平均基线 (E_init)
  E_baseline <- colMeans(WT_mat)

  # 3. 运行核心 RWR 算法，获取干预后的群体平均稳态向量
  ko_vector <- run_virtual_knockout(
    target_genes = target_genes,
    E_init       = E_baseline,
    adj_matrix   = adj_matrix
  )

  # 4. 利用桥接函数，将平均偏移量广播回所有个体样本 (依托 network_utils.R)
  KO_mat <- build_ko_matrix(
    WT_matrix   = WT_mat,
    ko_vector   = ko_vector,
    WT_baseline = E_baseline,
    method      = method
  )

  message(">>> [Batch Perturbation] 队列级干预矩阵构建完成。")

  return(list(
    WT = WT_mat,
    KO = KO_mat
  ))
}
