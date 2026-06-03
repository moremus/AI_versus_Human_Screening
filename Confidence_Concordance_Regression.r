# ============================================================
#  Confidence_Concordance_Regression
#  AI vs. Human Screening in Systematic Reviews
# ============================================================
#
#  Required packages:
#    install.packages(c("dplyr", "gtsummary", "gt"))
# ============================================================
 
 
library(dplyr)
library(gtsummary)
library(gt)


# ── 1. Regression - Title and Abstract Screening ─────────────────────────────────────────────────
df$Concordance <- factor(df$Concordance, levels = c("Agree", "Disagree"))

df <- df %>%
  mutate(Confidence2 = case_when(
    Confidence == "Very High" ~ "Very High",
    Confidence == "High"      ~ "High",
    TRUE                      ~ "Medium/Low"
  ),
  Confidence2 = factor(Confidence2, 
                       levels = c("Medium/Low", "High", "Very High")))

tiab <- glm(Concordance ~ Confidence2, family = "binomial", data = df)
summary(tiab)

t2 <-
  tbl_regression(tiab, exp=T,
                 label=list(Confidence2 = "Confidence")) |> modify_header(label ~ "**Variable**") |>
  bold_labels() |>
  as_gt() |>
  gt::tab_header(title = md("**Table 2:** Association between Loon Lens screening confidence and disagreement between artificial intelligence and human screening – title and abstract screening"))
 
 
# ── 2. Regression - Full-text Screening ─────────────────────────────────────────────────
ds$Concordance <- factor(ds$Concordance, levels = c("Agree", "Disagree"))

ds$Confidence <- factor(ds$Confidence, levels = c("Low", "Medium", "High", "Very High"))

ft <- glm(Concordance ~ Confidence, family = "binomial", data = ds)
summary(ft)

t3 <-
  tbl_regression(ft, exp=T,
                 label=list(Confidence="Confidence")) |> modify_header(label ~ "**Variable**") |>
  bold_labels() |>
  as_gt() |>
  gt::tab_header(title = md("**Table 3:** Association between Loon Lens screening confidence and disagreement between artificial intelligence and human screening – full-text screening"))
  
  ### END OF PROGRAM ###