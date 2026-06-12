#####################################################
### Class-imbalance-adjusted Performance Analysis ###
#####################################################

# ==================================================
# 1. PERFORMANCE STATISTICS OBTAINED BY CHANCE ALONE
# ==================================================

# --- Dataset parameters ---

# TiAb comparisons (human reference standard)
N_tiab        <- 793
N_inc_human   <- 61
N_exc_human   <- N_tiab - N_inc_human    # 732

# TiAb comparison (Catchii reference standard)
N_inc_catchii <- 58
N_exc_catchii <- N_tiab - N_inc_catchii  # 735

# Full-text comparison (human reference standard)
N_ft          <- 61
N_inc_ft      <- 41
N_exc_ft      <- 20

# ============================================================
# FUNCTION: compute null model metrics
# majority_class = "exclude" (predict all out)
#                = "include" (predict all in)
# ============================================================

null_model <- function(N_total, N_positive, N_negative,
                       majority_class = "exclude") {
  
  if (majority_class == "exclude") {
    TP <- 0
    TN <- N_negative
    FP <- 0
    FN <- N_positive
  } else {
    TP <- N_positive
    TN <- 0
    FP <- N_negative
    FN <- 0
  }
  
  sensitivity <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  specificity <- ifelse((TN + FP) > 0, TN / (TN + FP), 0)
  PPV         <- ifelse((TP + FP) > 0, TP / (TP + FP), NA)
  NPV         <- ifelse((TN + FN) > 0, TN / (TN + FN), NA)
  accuracy    <- (TP + TN) / N_total
  F1          <- ifelse((2*TP + FP + FN) > 0,
                        2*TP / (2*TP + FP + FN), 0)
  
  p_obs       <- (TP + TN) / N_total
  p_exp       <- ((TP + FP)/N_total) * ((TP + FN)/N_total) +
    ((TN + FN)/N_total) * ((TN + FP)/N_total)
  kappa       <- ifelse((1 - p_exp) > 0,
                        (p_obs - p_exp) / (1 - p_exp), 0)
  
  data.frame(
    Metric     = c("Sensitivity", "Specificity", "PPV",
                   "NPV", "Kappa", "Accuracy", "F1 Score"),
    Null_Value = c(
      round(sensitivity, 2),
      round(specificity, 2),
      ifelse(is.na(PPV), "Undefined", as.character(round(PPV, 2))),
      ifelse(is.na(NPV), "Undefined", as.character(round(NPV, 2))),
      round(kappa, 2),
      round(accuracy, 2),
      round(F1, 2)
    )
  )
}

# ============================================================
# COMPUTE NULL MODEL FOR EACH COMPARISON
# ============================================================

null_cat_h_tiab  <- null_model(N_tiab, N_inc_human,
                               N_exc_human,
                               majority_class = "exclude")

null_ll_h_tiab   <- null_model(N_tiab, N_inc_human,
                               N_exc_human,
                               majority_class = "exclude")

null_ll_cat_tiab <- null_model(N_tiab, N_inc_catchii,
                               N_exc_catchii,
                               majority_class = "exclude")

null_ll_h_ft     <- null_model(N_ft, N_inc_ft,
                               N_exc_ft,
                               majority_class = "exclude")

# ============================================================
# PRINT TABLES
# ============================================================

cat("Null Model: Catchii vs. Human - TiAb\n")
print(null_cat_h_tiab, row.names = FALSE)

cat("\nNull Model: Loon Lens vs. Human - TiAb\n")
print(null_ll_h_tiab, row.names = FALSE)

cat("\nNull Model: Loon Lens vs. Catchii - TiAb\n")
print(null_ll_cat_tiab, row.names = FALSE)

cat("\nNull Model: Loon Lens vs. Human - Full Text\n")
print(null_ll_h_ft, row.names = FALSE)

# ============================================================
# OPTIONAL: combine all into one summary table
# ============================================================

null_cat_h_tiab$Comparison  <- "Catchii vs. Human (TiAb)"
null_ll_h_tiab$Comparison   <- "Loon Lens vs. Human (TiAb)"
null_ll_cat_tiab$Comparison <- "Loon Lens vs. Catchii (TiAb)"
null_ll_h_ft$Comparison     <- "Loon Lens vs. Human (FT)"

null_summary <- rbind(null_cat_h_tiab,
                      null_ll_h_tiab,
                      null_ll_cat_tiab,
                      null_ll_h_ft)

null_summary <- null_summary[, c("Comparison", "Metric",
                                 "Null_Value")]

cat("\nSummary: All Null Model Values\n")
print(null_summary, row.names = FALSE)

# _____________________________________________________________________

# ============================================================
# 2. CLASS-IMBALANCE-ADJUSTED PERFORMANCE ANALYSIS
# ============================================================
# Positive event = "In"
# Null model = excludes all articles (majority class = "Out")
# R = 1000 bootstrap resamples
# ============================================================

library(boot)

set.seed(123)


run_analysis <- function(data, ref_col, ai_col,
                         comparison_label) {

  cat("\n============================================================\n")
  cat(comparison_label, "\n")
  cat("============================================================\n")

  fns   <- make_boot_fns(ref_col, ai_col)
  nulls <- null_metrics(data, ref_col)

  metrics <- list(
    list(name = "Sensitivity",  fn = fns$sens,
         null = nulls$sensitivity),
    list(name = "Specificity",  fn = fns$spec,
         null = nulls$specificity),
    list(name = "PPV",          fn = fns$ppv,
         null = nulls$PPV),
    list(name = "NPV",          fn = fns$npv,
         null = nulls$NPV),
    list(name = "Kappa",        fn = fns$kappa,
         null = nulls$kappa),
    list(name = "Consistency",  fn = fns$consistency,
         null = nulls$consistency),
    list(name = "F1 Score",     fn = fns$f1,
         null = nulls$F1)
  )

  results <- data.frame(
    Metric        = character(),
    Estimate      = numeric(),
    CI_Lower      = character(),
    CI_Upper      = character(),
    Null_Value    = character(),
    Exceeds_Null  = character(),
    stringsAsFactors = FALSE
  )

  for (m in metrics) {

    res <- boot(data, m$fn, R = 1000)
    est <- round(m$fn(data, 1:nrow(data)), 2)

    # Handle cases where boot.ci cannot be computed
    ci_result <- tryCatch(
      {
        ci  <- boot.ci(res, type = "perc")
        lo  <- round(ci$percent[4], 2)
        hi  <- round(ci$percent[5], 2)
        list(lo = lo, hi = hi, note = "")
      },
      error = function(e) {
        list(lo = "N/A", hi = "N/A",
             note = " [CI not estimable: all bootstrap values equal]")
      },
      warning = function(w) {
        list(lo = "N/A", hi = "N/A",
             note = " [CI not estimable: all bootstrap values equal]")
      }
    )

    lo   <- ci_result$lo
    hi   <- ci_result$hi
    note <- ci_result$note

    null_val <- m$null

    # Check if lower CI bound exceeds null value
    if (is.character(null_val) && null_val == "Undefined") {
      exceeds <- "N/A (PPV undefined for null)"
    } else if (lo == "N/A") {
      exceeds <- paste0("N/A", note)
    } else {
      exceeds <- ifelse(as.numeric(lo) > as.numeric(null_val),
                        "YES", "NO")
    }

    results <- rbind(results, data.frame(
      Metric       = m$name,
      Estimate     = est,
      CI_Lower     = as.character(lo),
      CI_Upper     = as.character(hi),
      Null_Value   = as.character(null_val),
      Exceeds_Null = exceeds,
      stringsAsFactors = FALSE
    ))
  }

  print(results, row.names = FALSE)
  cat("\nNote: 'Exceeds_Null' = YES indicates the AI tool performs\n")
  cat("statistically significantly better than majority-class guessing.\n")

  invisible(results)
}

# ------------------------------------------------------------
# Re-run all four comparisons with updated function
# ------------------------------------------------------------

results_cat_hum_tiab <- run_analysis(
  data             = cat_hum_tiab,
  ref_col          = "human_include",
  ai_col           = "catchii_include",
  comparison_label = "Catchii vs. Human --- TiAb"
)

results_ll_hum_tiab <- run_analysis(
  data             = ll_hum_tiab,
  ref_col          = "human_include",
  ai_col           = "loon_include",
  comparison_label = "Loon Lens vs. Human --- TiAb"
)

results_ll_cat_tiab <- run_analysis(
  data             = ll_cat_tiab,
  ref_col          = "catchii_include",
  ai_col           = "loon_include",
  comparison_label = "Loon Lens vs. Catchii --- TiAb"
)

results_ll_hum_ft <- run_analysis(
  data             = ll_hum_ft,
  ref_col          = "human_include",
  ai_col           = "loon_include",
  comparison_label = "Loon Lens vs. Human --- Full Text"
)