# ==============================================================================
# bulkKnk v0.6.0 — Full Workflow Demo with TCGA-STAD Data
# ==============================================================================
#
# This script runs the complete bulkKnk pipeline end-to-end:
#   Module 0: Feature selection (preserve signature genes + HVGs)
#   Module 1: WGCNA hub gene discovery
#   Module 2: GENIE3 causal network inference + DPI pruning
#   Module 3: RWR virtual knockout simulation
#   Module 4: AI phenotype prediction (ssGSEA + RF + survival)
#   Downstream: Volcano plots, GO/KEGG enrichment
#
# Requirements: bulkKnk installed or devtools::load_all()
# ==============================================================================

# --- 0. Bootstrap --------------------------------------------------------------
cat("\n")
cat("==============================================================\n")
cat("  bulkKnk v0.6.0 — Full Workflow Demo (TCGA-STAD)\n")
cat("==============================================================\n\n")

# Load package (always use load_all to pick up latest source changes)
cat(">>> Loading bulkKnk (devtools::load_all)...\n")
devtools::load_all()

# Load built-in demo data
data("tcga_stad_mini")

expr_mat  <- tcga_stad_mini$expr_matrix   # rows = samples, cols = genes
pheno_vec <- tcga_stad_mini$pheno_data$SampleType
target_pheno <- "Tumor"

n_samples <- nrow(expr_mat)
n_genes   <- ncol(expr_mat)
n_tumor   <- sum(pheno_vec == "Tumor")
n_normal  <- sum(pheno_vec == "Normal")

cat(sprintf(">>> Data loaded: %d samples (%d Tumor, %d Normal) x %d genes\n",
            n_samples, n_tumor, n_normal, n_genes))

# Create output directory for figures
out_dir <- file.path(tempdir(), "bulkKnk_demo_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf(">>> Output directory: %s\n", out_dir))

# Check optional dependencies
has_survival   <- requireNamespace("survival",   quietly = TRUE)
has_survminer  <- requireNamespace("survminer",  quietly = TRUE)
has_ggplot2    <- requireNamespace("ggplot2",    quietly = TRUE)
has_patchwork  <- requireNamespace("patchwork",  quietly = TRUE)
has_cp         <- requireNamespace("clusterProfiler", quietly = TRUE)

cat(sprintf(">>> Optional deps: survival=%s, survminer=%s, ggplot2=%s, patchwork=%s, clusterProfiler=%s\n",
            has_survival, has_survminer, has_ggplot2, has_patchwork, has_cp))


# ==============================================================================
# Module 0: Feature Selection
# ==============================================================================
cat("\n--- Module 0: Feature Selection -----------------------------------\n")

sub_expr <- select_advanced_features(expr_mat, n_features = 500)
cat(sprintf("    Result: %d samples x %d genes\n", nrow(sub_expr), ncol(sub_expr)))


# ==============================================================================
# Module 1: WGCNA Hub Gene Discovery
# ==============================================================================
cat("\n--- Module 1: WGCNA Hub Gene Discovery ---------------------------\n")

# Convert phenotype to numeric for WGCNA correlation
trait_num <- as.numeric(as.factor(pheno_vec))

wgcna_res <- tryCatch({
  identify_hub_modules(sub_expr, trait_vec = trait_num)
}, error = function(e) {
  cat(sprintf("    WGCNA failed: %s\n", e$message))
  cat("    Falling back to variance-based gene selection...\n")
  NULL
})

if (!is.null(wgcna_res)) {
  hub_genes  <- head(wgcna_res$best_genes, 100)
  ko_target  <- hub_genes[1]
  cat(sprintf("    Best module: %s with %d genes\n",
              wgcna_res$best_module, length(wgcna_res$best_genes)))
  cat(sprintf("    Selected KO target: %s (rank #1 hub)\n", ko_target))
} else {
  # Fallback: use top-variance gene as target
  gene_vars <- apply(sub_expr, 2, var)
  hub_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]
  ko_target <- hub_genes[1]
  cat(sprintf("    Fallback KO target: %s (top-variance gene)\n", ko_target))
}


# ==============================================================================
# Module 2: Causal Network Inference
# ==============================================================================
cat("\n--- Module 2: GENIE3 + DPI Pruning -------------------------------\n")

focus_mat <- as.matrix(sub_expr[, hub_genes, drop = FALSE])
cat(sprintf("    Focus matrix: %d samples x %d genes\n", nrow(focus_mat), ncol(focus_mat)))

# GENIE3 inference (may take ~30s for 100 genes)
cat("    Running GENIE3 (this may take a moment)...\n")
net_res <- tryCatch({
  infer_causal_network(focus_mat, n_cores = 1)
}, error = function(e) {
  cat(sprintf("    GENIE3 failed: %s\n", e$message))
  NULL
})

if (is.null(net_res)) {
  stop("Network inference failed. Cannot continue pipeline.")
}

cat(sprintf("    GENIE3 done: %d edges retained (threshold=%.4f)\n",
            net_res$edge_stats["Retained"], net_res$threshold_used))

# DPI pruning
net_clean <- prune_network_dpi_fast(net_res$clean_matrix, dpi_threshold = 0.5)

# Re-attach gene names (safety net in case GENIE3 drops them)
rownames(net_clean) <- colnames(focus_mat)
colnames(net_clean) <- colnames(focus_mat)

# Normalize: column-normalize the transpose so row sums <= 1
W_norm <- t(net_clean)
col_sums <- colSums(W_norm) + 1e-9
W_norm <- sweep(W_norm, 2, col_sums, FUN = "/")

cat(sprintf("    Network matrix: %d x %d, density: %.3f\n",
            nrow(W_norm), ncol(W_norm), mean(W_norm > 0)))


# ==============================================================================
# Module 3: Virtual Knockout Simulation
# ==============================================================================
cat("\n--- Module 3: RWR Virtual Knockout -------------------------------\n")

ko_results <- batch_virtual_knockout(
  expr_matrix  = focus_mat,
  pheno_vec    = pheno_vec,
  target_pheno = target_pheno,
  target_genes = ko_target,
  adj_matrix   = W_norm
)

cat(sprintf("    WT matrix: %d x %d\n", nrow(ko_results$WT), ncol(ko_results$WT)))
cat(sprintf("    KO matrix: %d x %d\n", nrow(ko_results$KO), ncol(ko_results$KO)))

# Compute delta expression for the target
wt_mean_target  <- mean(ko_results$WT[, ko_target])
ko_mean_target  <- mean(ko_results$KO[, ko_target])
cat(sprintf("    Target [%s]: WT mean=%.2f -> KO mean=%.2f (delta=%.2f)\n",
            ko_target, wt_mean_target, ko_mean_target, ko_mean_target - wt_mean_target))


# ==============================================================================
# Module 4: AI Phenotype Prediction (Layers 1-3)
# ==============================================================================
cat("\n--- Module 4: AI Phenotype Prediction ----------------------------\n")

# Baseline labels: split tumor samples by median target expression
baseline_labels <- ifelse(
  ko_results$WT[, ko_target] > median(ko_results$WT[, ko_target]),
  "HighRisk", "LowRisk"
)

# predict_ko_phenotype requires: rows=genes, cols=samples -> transpose
# NOTE: Some signature scores may be NA because the 100-gene focus set
# may not contain all genes from the builtin signatures (e.g., CD8A, MKI67).
# This is expected with small gene sets. For production use, increase the
# network size or provide custom_signatures.
report <- predict_ko_phenotype(
  WT_matrix    = t(ko_results$WT),
  KO_matrix    = t(ko_results$KO),
  pheno_labels = baseline_labels,
  target_pheno = "HighRisk",
  ko_gene      = ko_target
)

cat(sprintf("    OOB accuracy: %.1f%%\n", report$oob_accuracy * 100))
cat(sprintf("    Reversal rate: %.1f%%\n", report$reversal_rate * 100))
cat(sprintf("    Benefit probability: %.1f%%\n",
            ifelse(is.na(report$benefit_probability), 0, report$benefit_probability) * 100))


# --- Figure 1: KO Delta Signature Plot ----------------------------------------
if (has_ggplot2) {
  cat("\n    Generating Figure 1: KO Delta Signature Plot...\n")
  p1 <- tryCatch({
    plot_ko_delta(report, color_by = "category")
  }, error = function(e) {
    cat(sprintf("    Warning: %s\n", e$message))
    NULL
  })
  if (!is.null(p1)) {
    fig1_path <- file.path(out_dir, "fig1_ko_delta_signatures.png")
    ggplot2::ggsave(fig1_path, p1, width = 10, height = 6, dpi = 150)
    cat(sprintf("    Saved: %s\n", fig1_path))
  }
}


# ==============================================================================
# Module 4 — Layer 4: Survival Analysis (with simulated clinical data)
# ==============================================================================
cat("\n--- Module 4 — Layer 4: Survival Prognosis -----------------------\n")

if (has_survival) {
  # The demo data lacks OS_time/OS_status, so we simulate realistic clinical data
  cat("    Simulating clinical survival data for Layer 4 demo...\n")
  set.seed(2026)
  n_tumor_samples <- nrow(ko_results$WT)

  # Simulate OS times: log-normal distribution (~24 months median)
  sim_os_time <- round(rlnorm(n_tumor_samples, log(24), 0.6), 1)
  # Simulate events: ~40% death rate
  sim_os_status <- rbinom(n_tumor_samples, 1, 0.4)

  sim_clinical <- data.frame(
    sample_id = rownames(ko_results$WT),
    OS_time   = sim_os_time,
    OS_status = sim_os_status,
    stringsAsFactors = FALSE
  )

  # Re-run prediction with clinical data for Layer 4
  report_full <- tryCatch({
    predict_ko_phenotype(
      WT_matrix    = t(ko_results$WT),
      KO_matrix    = t(ko_results$KO),
      pheno_labels = baseline_labels,
      target_pheno = "HighRisk",
      ko_gene      = ko_target,
      clinical_data = sim_clinical,
      time_unit    = "month"
    )
  }, error = function(e) {
    cat(sprintf("    Layer 4 failed: %s\n", e$message))
    report  # Fall back to the 3-layer report
  })

  # --- Figure 2: KM Survival Curves -------------------------------------------
  if (has_ggplot2 && has_patchwork &&
      !is.null(report_full$prognosis) &&
      !is.null(report_full$prognosis$km_data_wt)) {
    cat("\n    Generating Figure 2: KM Survival Curves (WT vs KO)...\n")
    fig2_path <- file.path(out_dir, "fig2_km_survival.png")
    p2 <- tryCatch({
      plot_km_comparison(report_full)
    }, error = function(e) {
      cat(sprintf("    Warning: plot_km_comparison failed: %s\n", e$message))
      NULL
    })
    if (!is.null(p2)) {
      saved_fig2 <- tryCatch({
        ggplot2::ggsave(fig2_path, p2, width = 14, height = 7, dpi = 150)
        file.exists(fig2_path)
      }, error = function(e) {
        cat(sprintf("    ggsave failed (%s), trying png fallback...\n", e$message))
        try(dev.off(), silent = TRUE)
        tryCatch({
          png(fig2_path, width = 2100, height = 1050, res = 150)
          print(p2)
          dev.off()
          file.exists(fig2_path)
        }, error = function(e2) {
          cat(sprintf("    Warning: Could not save KM plot: %s\n", e2$message))
          FALSE
        })
      })
      if (saved_fig2) cat(sprintf("    Saved: %s\n", fig2_path))
    }
  }

  # --- Figure 3: Risk Score Distribution --------------------------------------
  if (has_ggplot2 &&
      !is.null(report_full$prognosis) &&
      !is.null(report_full$prognosis$risk_score_wt)) {
    cat("    Generating Figure 3: Risk Score Distribution...\n")
    p3 <- tryCatch({
      plot_risk_score_shift(report_full)
    }, error = function(e) {
      cat(sprintf("    Warning: %s\n", e$message))
      NULL
    })
    if (!is.null(p3)) {
      fig3_path <- file.path(out_dir, "fig3_risk_score_shift.png")
      ggplot2::ggsave(fig3_path, p3, width = 8, height = 5, dpi = 150)
      cat(sprintf("    Saved: %s\n", fig3_path))
    }
  }

  # --- Figure 4: Cox Forest Plot ----------------------------------------------
  if (has_ggplot2 &&
      !is.null(report_full$prognosis) &&
      !is.null(report_full$prognosis$cox_results)) {
    cat("    Generating Figure 4: Cox Forest Plot...\n")
    p4 <- tryCatch({
      plot_cox_forest(report_full)
    }, error = function(e) {
      cat(sprintf("    Warning: %s\n", e$message))
      NULL
    })
    if (!is.null(p4)) {
      fig4_path <- file.path(out_dir, "fig4_cox_forest.png")
      ggplot2::ggsave(fig4_path, p4, width = 9, height = 6, dpi = 150)
      cat(sprintf("    Saved: %s\n", fig4_path))
    }
  }

} else {
  cat("    Skipped (survival package not available)\n")
  report_full <- report
}


# ==============================================================================
# Downstream: Dual Volcano + Enrichment
# ==============================================================================
cat("\n--- Downstream: Dual Volcano Plot --------------------------------\n")

# Compute real clinical DEGs (Tumor vs Normal) for the volcano
is_tumor   <- pheno_vec == "Tumor"
expr_tumor <- sub_expr[is_tumor, , drop = FALSE]
expr_norm  <- sub_expr[!is_tumor, , drop = FALSE]

mean_t <- colMeans(expr_tumor)
mean_n <- colMeans(expr_norm)

real_log2fc <- log2((mean_t + 0.1) / (mean_n + 0.1))

real_pvals <- sapply(seq_len(ncol(sub_expr)), function(i) {
  tryCatch(t.test(expr_tumor[, i], expr_norm[, i])$p.value, error = function(e) 1)
})

real_deg_df <- data.frame(
  Gene   = colnames(sub_expr),
  log2FC = real_log2fc,
  pvalue = real_pvals,
  stringsAsFactors = FALSE
)

if (has_ggplot2 && has_patchwork) {
  cat("    Generating Figure 5: Dual Volcano Plot...\n")
  p5 <- tryCatch({
    plot_dual_volcano(
      real_deg_df  = real_deg_df,
      WT_matrix    = ko_results$WT,
      KO_matrix    = ko_results$KO,
      target_genes = ko_target
    )
  }, error = function(e) {
    cat(sprintf("    Warning: %s\n", e$message))
    NULL
  })
  if (!is.null(p5)) {
    fig5_path <- file.path(out_dir, "fig5_dual_volcano.png")
    ggplot2::ggsave(fig5_path, p5, width = 14, height = 6, dpi = 150)
    cat(sprintf("    Saved: %s\n", fig5_path))
  }
}


# --- Enrichment Analysis -------------------------------------------------------
cat("\n--- Downstream: GO/KEGG Enrichment -------------------------------\n")

wt_mean <- colMeans(ko_results$WT)
ko_mean <- colMeans(ko_results$KO)
log2fc_virtual <- log2((ko_mean + 0.1) / (wt_mean + 0.1))

virtual_degs <- names(log2fc_virtual)[abs(log2fc_virtual) > 0.05]
virtual_degs <- setdiff(virtual_degs, ko_target)  # Exclude target itself

cat(sprintf("    Found %d cascade-affected genes (|log2FC| > 0.05)\n", length(virtual_degs)))

if (length(virtual_degs) >= 5 && has_cp) {
  # GO enrichment
  cat("    Running GO enrichment...\n")
  go_res <- tryCatch({
    run_vk_enrichment(virtual_degs, species = "human", database = "GO")
  }, error = function(e) {
    cat(sprintf("    GO enrichment failed: %s\n", e$message))
    NULL
  })

  if (!is.null(go_res) && nrow(go_res) > 0 && has_ggplot2) {
    cat("    Generating Figure 6: GO Enrichment Dotplot...\n")
    p6 <- tryCatch({
      clusterProfiler::dotplot(go_res, showCategory = 10,
                               title = paste("GO: BP after", ko_target, "KO"))
    }, error = function(e) {
      cat(sprintf("    Warning: %s\n", e$message))
      NULL
    })
    if (!is.null(p6)) {
      fig6_path <- file.path(out_dir, "fig6_go_enrichment.png")
      ggplot2::ggsave(fig6_path, p6, width = 10, height = 6, dpi = 150)
      cat(sprintf("    Saved: %s\n", fig6_path))
    }
  }

  # KEGG enrichment
  cat("    Running KEGG enrichment...\n")
  kegg_res <- tryCatch({
    run_vk_enrichment(virtual_degs, species = "human", database = "KEGG")
  }, error = function(e) {
    cat(sprintf("    KEGG enrichment failed: %s\n", e$message))
    NULL
  })

  if (!is.null(kegg_res) && nrow(kegg_res) > 0 && has_ggplot2) {
    cat("    Generating Figure 7: KEGG Enrichment Dotplot...\n")
    p7 <- tryCatch({
      clusterProfiler::dotplot(kegg_res, showCategory = 10,
                               title = paste("KEGG after", ko_target, "KO"))
    }, error = function(e) {
      cat(sprintf("    Warning: %s\n", e$message))
      NULL
    })
    if (!is.null(p7)) {
      fig7_path <- file.path(out_dir, "fig7_kegg_enrichment.png")
      ggplot2::ggsave(fig7_path, p7, width = 10, height = 6, dpi = 150)
      cat(sprintf("    Saved: %s\n", fig7_path))
    }
  }
} else if (length(virtual_degs) < 5) {
  cat("    Too few cascade-affected genes for enrichment analysis.\n")
} else {
  cat("    Skipped (clusterProfiler not available)\n")
}


# ==============================================================================
# Final Report
# ==============================================================================
cat("\n================================================================\n")
cat("  FINAL REPORT\n")
cat("================================================================\n")
cat(sprintf("  KO target gene    : %s\n", ko_target))
cat(sprintf("  Samples analyzed  : %d Tumor (perturbed)\n", n_tumor))
cat(sprintf("  Network genes     : %d\n", length(hub_genes)))
cat(sprintf("  RF OOB accuracy   : %.1f%%\n", report$oob_accuracy * 100))
cat(sprintf("  Phenotype reversal: %.1f%%\n", report$reversal_rate * 100))
cat(sprintf("  Benefit prob      : %.1f%%\n",
            ifelse(is.na(report$benefit_probability), 0,
                   report$benefit_probability) * 100))

if (!is.null(report_full$prognosis) && !is.null(report_full$prognosis$logrank_p)) {
  prog <- report_full$prognosis
  cat(sprintf("  log-rank p        : %.4f\n", prog$logrank_p))
  cat(sprintf("  Rescue rate       : %.1f%%\n", prog$rescue_rate * 100))
  if (!is.na(prog$predicted_os_gain_months)) {
    cat(sprintf("  Predicted OS gain : +%.1f months\n", prog$predicted_os_gain_months))
  }
}

cat("\n  Key immune signature shifts:\n")
top_immune <- head(report$immune_summary, 3)
for (i in seq_len(nrow(top_immune))) {
  cat(sprintf("    %-30s %s (delta=%.4f)\n",
              top_immune$signature[i], top_immune$direction[i], top_immune$delta_score[i]))
}

cat("\n  Key hallmark shifts:\n")
top_hall <- head(report$malignancy_summary, 3)
for (i in seq_len(nrow(top_hall))) {
  cat(sprintf("    %-30s %s (delta=%.4f)\n",
              top_hall$signature[i], top_hall$direction[i], top_hall$delta_score[i]))
}

cat(sprintf("\n  Figures saved to: %s\n", out_dir))

# List generated figures
figs <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
for (f in figs) cat(sprintf("    - %s\n", basename(f)))

cat("\n  AI Interpretation:\n")
cat(paste0("  ", report$interpretation, "\n"))

cat("\n================================================================\n")
cat("  Pipeline complete.\n")
cat("================================================================\n")
