##################################################################################
# Title: Title and Abstract Screening Comparison
# Description: All analyses for the title and abstract comparison of Loon Lens
# 			   versus human screening.
# Project: AI versus Human Screening
# Note: The same code can be utilized to conduct the other screening comparisons.
##################################################################################

# ------------------------------- #
# 1. Load Libraries
# ------------------------------- #

library(Hmisc)
library(revtools)
library(tidyverse)
library(common)
library(dplyr)
library(gtsummary)
library(gt)
library(questionr)
library(yardstick)
library(rsample)
library(ggplot2)
library(fmsb)
library(boot)
library(table1)

# ------------------------------- #
# 2. Import Dataset
# ------------------------------- #

doc <- read.csv("C:/Loon_Lens_Title_Abstract_Screening_Results.csv")

# ------------------------------- #
# 3. Format Variables
# ------------------------------- #

# Create Reference Standard Variable
doc$Reference <- factor(doc$Human_Decision, levels = c(1,0))
label(doc$Reference) <- "Human Screening"

# Create Concordance Variable
doc$Concordance <- ifelse(doc$Human_Decision != doc$Loon_Lens_Decision, 1,0)
doc$Concordance <- factor(doc$Concordance, levels = c(1,0), labels = c("Disagree","Agree"))

# ------------------------------- #
# 4. Tables
# ------------------------------- #

# Descriptive Statistics
doc |>
  tbl_summary(include = c(Human_Decision, Loon_Lens_Decision, 
                          Concordance)) |>
  modify_header(label ~ "**Variable**") |>
  as_gt() |>
  gt::tab_header(
    title=md("**Descriptive Statistics**")
  )

# Confusion Table
doc |>
  tbl_cross(row = Loon_Lens_Decision, col = Human_Decision) |>
  as_gt() |>
  gt::tab_header(
    title=md("**Confusion Table**")
  )

## Put the Table in a Dataframe (for comparative statistics)
Human<-c("Included", "Included", "Excluded", "Excluded")
Loon.Lens<-c("Included", "Excluded", "Included", "Excluded")
count<-c(40, 21, 4, 728)

descr1<-as.data.frame(cbind(Human, Loon.Lens, count))
descr1$Human<-factor(descr1$Human, levels=c("Included", "Excluded"))
descr1$Loon.Lens<-factor(descr1$Loon.Lens, levels=c("Included", "Excluded"))
descr1$count<-as.numeric(descr1$count)
descr2<-descr1[rep(seq_len(nrow(descr1)), descr1$count), c("Human", "Loon.Lens")]

# ------------------------------- #
# 5. Comparative Statistics
# ------------------------------- #

# Point Estimates and 95% Bootstrap CI for SENS, SPEC, PPV, NPV, Accuracy, Kappa, F1
set.seed(123)

human_include<-c("In", "In", "Out", "Out") # Create confusion table
ai_include<-c("In", "Out", "In", "Out")
count<-c(40, 21, 4, 728)
descr3<-as.data.frame(cbind(human_include, ai_include, count)) # Create a dataframe for the confusion table

dataset <- data.frame(										# Expand the confusion table into a full length dataset
  human_include = rep(descr3$human_include, descr3$count),	 
  ai_include = rep(descr3$ai_include, descr3$count))

# SENS bootstrapping
boot_sens <- function(data, indices) {
  d  <- data[indices, ]
  TP <- sum(d$ai_include == "In"  & d$human_include == "In")
  FN <- sum(d$ai_include == "Out" & d$human_include == "In")
  TP / (TP + FN)
}
results_sens <- boot(dataset, boot_sens, R = 1000)
ci_sens      <- boot.ci(results_sens, type = "perc")
cat("Sensitivity\n")
cat("  Estimate:", boot_sens(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_sens$percent[4], ",", ci_sens$percent[5], "]\n\n")

# SPEC bootstrapping
boot_spec <- function(data, indices) {
  d  <- data[indices, ]
  TN <- sum(d$ai_include == "Out" & d$human_include == "Out")
  FP <- sum(d$ai_include == "In"  & d$human_include == "Out")
  TN / (TN + FP)
}
results_spec <- boot(dataset, boot_spec, R = 1000)
ci_spec      <- boot.ci(results_spec, type = "perc")
cat("Specificity\n")
cat("  Estimate:", boot_spec(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_spec$percent[4], ",", ci_spec$percent[5], "]\n\n")

# PPV bootstrapping
boot_ppv <- function(data, indices) {
  d  <- data[indices, ]
  TP <- sum(d$ai_include == "In" & d$human_include == "In")
  FP <- sum(d$ai_include == "In" & d$human_include == "Out")
  TP / (TP + FP)
}
results_ppv <- boot(dataset, boot_ppv, R = 1000)
ci_ppv      <- boot.ci(results_ppv, type = "perc")
cat("PPV\n")
cat("  Estimate:", boot_ppv(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_ppv$percent[4], ",", ci_ppv$percent[5], "]\n\n")

# NPV bootstrapping
boot_npv <- function(data, indices) {
  d  <- data[indices, ]
  TN <- sum(d$ai_include == "Out" & d$human_include == "Out")
  FN <- sum(d$ai_include == "Out" & d$human_include == "In")
  TN / (TN + FN)
}
results_npv <- boot(dataset, boot_npv, R = 1000)
ci_npv      <- boot.ci(results_npv, type = "perc")
cat("NPV\n")
cat("  Estimate:", boot_npv(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_npv$percent[4], ",", ci_npv$percent[5], "]\n\n")

# KAPPA bootstrapping
boot_kappa <- function(data, indices) {
  d  <- data[indices, ]
  n  <- nrow(d)
  TP <- sum(d$ai_include == "In"  & d$human_include == "In")
  TN <- sum(d$ai_include == "Out" & d$human_include == "Out")
  FP <- sum(d$ai_include == "In"  & d$human_include == "Out")
  FN <- sum(d$ai_include == "Out" & d$human_include == "In")
  p_observed <- (TP + TN) / n
  p_expected <- ((TP + FP) / n) * ((TP + FN) / n) +
    ((TN + FN) / n) * ((TN + FP) / n)
  (p_observed - p_expected) / (1 - p_expected)
}
results_kappa <- boot(dataset, boot_kappa, R = 1000)
ci_kappa      <- boot.ci(results_kappa, type = "perc")
cat("Kappa\n")
cat("  Estimate:", boot_kappa(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_kappa$percent[4], ",", ci_kappa$percent[5], "]\n\n")

# ACCURACY bootstrapping
boot_accuracy <- function(data, indices) {
  d  <- data[indices, ]
  TP <- sum(d$ai_include == "In"  & d$human_include == "In")
  TN <- sum(d$ai_include == "Out" & d$human_include == "Out")
  (TP + TN) / nrow(d)
}
results_accuracy <- boot(dataset, boot_accuracy, R = 1000)
ci_accuracy      <- boot.ci(results_accuracy, type = "perc")
cat("Accuracy\n")
cat("  Estimate:", boot_accuracy(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_accuracy$percent[4], ",", ci_accuracy$percent[5], "]\n\n")

# F1 bootstrapping
boot_f1 <- function(data, indices) {
  d  <- data[indices, ]
  TP <- sum(d$ai_include == "In"  & d$human_include == "In")
  FP <- sum(d$ai_include == "In"  & d$human_include == "Out")
  FN <- sum(d$ai_include == "Out" & d$human_include == "In")
  precision <- TP / (TP + FP)
  recall    <- TP / (TP + FN)
  2 * (precision * recall) / (precision + recall)
}
results_f1 <- boot(dataset, boot_f1, R = 1000)
ci_f1      <- boot.ci(results_f1, type = "perc")
cat("F1 Score\n")
cat("  Estimate:", boot_f1(dataset, 1:nrow(dataset)), "\n")
cat("  95% CI:  [", ci_f1$percent[4], ",", ci_f1$percent[5], "]\n\n")

# ------------------------------- #
# 6. Tabularized Statistics
# ------------------------------- #

results_table <- data.frame(
  Metric   = c("Sensitivity", "Specificity", "PPV", "NPV", 
               "Kappa", "Accuracy", "F1 Score"),
  Estimate = c(
    boot_sens(dataset,     1:nrow(dataset)),
    boot_spec(dataset,     1:nrow(dataset)),
    boot_ppv(dataset,      1:nrow(dataset)),
    boot_npv(dataset,      1:nrow(dataset)),
    boot_kappa(dataset,    1:nrow(dataset)),
    boot_accuracy(dataset, 1:nrow(dataset)),
    boot_f1(dataset,       1:nrow(dataset))
  ),
  CI_Lower = c(
    ci_sens$percent[4],
    ci_spec$percent[4],
    ci_ppv$percent[4],
    ci_npv$percent[4],
    ci_kappa$percent[4],
    ci_accuracy$percent[4],
    ci_f1$percent[4]
  ),
  CI_Upper = c(
    ci_sens$percent[5],
    ci_spec$percent[5],
    ci_ppv$percent[5],
    ci_npv$percent[5],
    ci_kappa$percent[5],
    ci_accuracy$percent[5],
    ci_f1$percent[5]
  )
)

# Round all numeric columns to 2 decimal places
results_table[, 2:4] <- round(results_table[, 2:4], 2)

# Print the table
print(results_table, row.names = FALSE)

# ------------------------------- #
# 7. Radar Chart
# ------------------------------- #

# Enter chart data
chart_dat <- data.frame(
  row.names = c("Data"),
  Sensitivity = c(66),
  Specificity = c(100),
  PPV = c(91),
  NPV = c(97),
  Concordance = c(97)
)

max_min <- data.frame(
  Sensitivity = c(100, 0), Specificity = c(100, 0),
  PPV = c(100, 0),NPV = c(100, 0), 
  Concordance = c(100, 0)
)
rownames(max_min) <- c("Max", "Min")

rad <- rbind(max_min, chart_dat)

# Generate and customize the chart
radarchart(
  rad, axistype = 1, vlcex = 0.7,
  pcol = "lightgrey", pfcol = (alpha("lightgrey", 0.5)), plwd = 2.5, plty = 1,
  caxislabcol = "black",
  title = "Figure 1. Summary of Performance Metrics - Loon Lens versus Human Screening (TiAb)"
)
  
# ------------------------------- #
# 8. Concordance Analysis
# ------------------------------- # 

# Create the Outcome Variable
doc$Confidence <- factor(doc$Loon_Lens_Confidence_Level, levels = c("VERY_HIGH", "HIGH", "MEDIUM", "LOW"))

doc$out <- ifelse(doc$Human_Decision == doc$Loon_Lens_Decision, 0, 1)
doc$out <- factor(doc$out, levels = c(1,0), labels=c("Disagree", "Agree"))
doc$out <- relevel(doc$out, ref = "Agree")
labels(doc) <- list(out = "Concordance")

chisq.test(doc$Confidence, doc$out)

table1(~ Confidence | out, doc,
       overall=F,
       caption="<b>Loon Lens Confidence Level by Screening Concordance - Title and Abstract Screening</b>",
       footnote="<b>Chi-square test:</b> p < 0.0001")
  
# End of Program #