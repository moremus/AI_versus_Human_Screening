###############################################################################
# Title: Full Text Screening Comparison
# Description: All analyses for the full-text comparison of Loon Lens
# 			   versus human screening.
# Project: AI versus Human Screening
# Note1: The numbers provided in the code below are for information purposes
#		 only and do not necessarily represent actual data.
# Note2: The same code was utlized to conduct the title and abstract screening
#		 comparisons.
###############################################################################

# ------------------------------- #
# 1. Load Libraries
# ------------------------------- #

library(Hmisc)
library(revtools)
library(tidyverse)
library(common)
library(dplyr)
library(ASySD)
library(gtsummary)
library(gt)
library(questionr)
library(yardstick)
library(rsample)
library(ggplot2)
library(fmsb)
library(boot)

# ------------------------------- #
# 2. Import Dataset
# ------------------------------- #

doc <- read.csv("C:/full_text_document.csv")

# ------------------------------- #
# 3. Format Variables
# ------------------------------- #

# Create Gold variable
doc$Gold <- factor(doc$Human, levels = c(1,0))
label(doc$Gold) <- "Human Screening"

# Create concordance variable
doc$Concordance <- ifelse(doc$Human != doc$Loon.Lens, 1,0)
doc$Concordance <- factor(doc$Concordance, levels = c(1,0), labels = c("Disagree","Agree"))

# ------------------------------- #
# 4. Tables
# ------------------------------- #

# Table 1 
doc |>
  tbl_summary(include = c(Human, Loon.Lens, 
                          Concordance)) |>
  modify_header(label ~ "**Variable**") |>
  as_gt() |>
  gt::tab_header(
    title=md("**Table 1.** Descriptive Statistics")
  )

# Table 2 - Confusion Table
doc |>
  tbl_cross(row = Loon.Lens, col = Human) |>
  as_gt() |>
  gt::tab_header(
    title=md("**Table 2.** Confusion Table")
  )

## Put the Table in a Dataframe (needed to get the comparative statistics)
Human<-c("Included", "Included", "Excluded", "Excluded")
Loon.Lens<-c("Included", "Excluded", "Included", "Excluded")
count<-c(41, 2, 11, 8)

descr1<-as.data.frame(cbind(Human, Loon.Lens, count))
descr1$Human<-factor(descr1$Human, levels=c("Included", "Excluded"))
descr1$Loon.Lens<-factor(descr1$Loon.Lens, levels=c("Included", "Excluded"))
descr1$count<-as.numeric(descr1$count)
descr2<-wtd.table(descr1$Loon.Lens, descr1$Human, weights=descr1$count)

# ------------------------------- #
# 5. Comparative Statistics
# ------------------------------- #

# SENS, SPEC, PPV, NPV, Accuracy, Kappa, F1
sens(descr2, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")
spec(descr2, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")
ppv(descr2, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")
npv(descr2, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")
accuracy(descr2, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")
kap(descr2, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")
f_meas(descr2, truth = Human, estimate = Loon.Lens, beta = 1, estimator = "binary", event_level = "first")

doc$Human <- factor(ifelse(doc$Human == 1, "Included", "Excluded"),levels = c("Included", "Excluded")) # Human Variable for Bootstrapping
doc$Loon.Lens <- factor(ifelse(doc$Loon.Lens == 1, "Included", "Excluded"),levels = c("Included", "Excluded")) # Loon.Lens variable for Bootstrapping


# 95% Bootsrap CI for SENS, SPEC, PPV, NPV, Accuracy, Kappa
set.seed(1)
bootstrap_data <- bootstraps(doc, times = 1000)

bootstrap_sens <- function(splits) {
  x <- analysis(splits)
  sens(x, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")$.estimate
}

bootstrap_spec <- function(splits) {
  x <- analysis(splits)
  spec(x, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")$.estimate
}

bootstrap_ppv <- function(splits) {
  x <- analysis(splits)
  ppv(x, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")$.estimate
}

bootstrap_npv <- function(splits) {
  x <- analysis(splits)
  npv(x, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")$.estimate
}

bootstrap_acc <- function(splits) {
  x <- analysis(splits)
  accuracy(x, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")$.estimate
}

bootstrap_kap <- function(splits) {
  x <- analysis(splits)
  kap(x, truth = Human, estimate = Loon.Lens, estimator = "binary", event_level = "first")$.estimate
}
bootstrap_f1 <- function(splits) {
  x <- analysis(splits)
  f_meas(x, truth = Human, estimate = Loon.Lens, beta = 1, estimator = "binary", event_level = "first")$.estimate
}

bootstrap_data <- bootstrap_data %>% 
  mutate(
    sens = map_dbl(splits, bootstrap_sens), 
    spec = map_dbl(splits, bootstrap_spec),
    ppv = map_dbl(splits, bootstrap_ppv),
    npv = map_dbl(splits, bootstrap_npv),
    acc = map_dbl(splits, bootstrap_acc),
    kap = map_dbl(splits, bootstrap_kap)
  )

quantile(bootstrap_data$sens, probs = c(0.025, 0.975))
quantile(bootstrap_data$spec, probs = c(0.025, 0.975))
quantile(bootstrap_data$ppv, probs = c(0.025, 0.975))
quantile(bootstrap_data$npv, probs = c(0.025, 0.975))
quantile(bootstrap_data$acc, probs = c(0.025, 0.975))
quantile(bootstrap_data$kap, probs = c(0.025, 0.975))

# FI bootstrapping 
human_include<-c("In", "In", "Out", "Out") # Create confusion table
ai_include<-c("In", "Out", "In", "Out")
count<-c(40, 2, 12, 8)

descr1<-as.data.frame(cbind(human_include, ai_include, count)) # Create a dataframe for the confusion table

dataset <- data.frame(										# Expand the confusion table into a full length dataset
  human_include = rep(descr1$human_include, descr1$count),	 
  ai_include = rep(descr1$ai_include, descr1$count)
)

# F1 - Run function to get the F1 score
f1_score <- function(data) {
  
  TP <- sum(data$ai_include == "In" & data$human_include == "In") #In the function, the gold standard goes second
  FP <- sum(data$ai_include == "In" & data$human_include == "Out")
  FN <- sum(data$ai_include == "Out" & data$human_include == "In")
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  2 * (precision * recall) / (precision + recall)
}

f1_score(dataset)


# F1 (95% CI) - Run function to get the confidence interval
boot_f1 <- function(data, indices) {
  
  d <- data[indices, ]
  
  TP <- sum(d$ai_include == "In" & d$human_include == "In") #In the function, the gold standard goes second
  FP <- sum(d$ai_include == "In" & d$human_include == "Out")
  FN <- sum(d$ai_include == "Out" & d$human_include == "In")
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  2 * (precision * recall) / (precision + recall)
}

results <- boot(dataset, boot_f1, R = 2000) # R is the number of resamples to calculate the bootstrap CI; change 2000 to whatever you are doing for the other CIs

boot.ci(results, type = "perc")

# ------------------------------- #
# 6. Tabularized Statistics
# ------------------------------- #

# Table 3 - Confusion Table Statistics and Kappa
Statistic<-c("Sensitivity", "Specificity", 
             "Positive Predictive Value",
             "Negative Predictive Value",
             "Accuracy", "Kappa", "F1")
Value<-c(0.95, 0.42, 0.79, 0.80, 0.79, 0.43, 0.86)
CI<-c("0.76, 1.0", "0.19, 0.65", "0.67, 0.89",
      "0.5, 1.0", "0.69, 0.88",
      "0.18, 0.66", "0.80, 0.98")

stats<-as.data.frame(cbind(Statistic, Value,
                           CI))
stats$Value<-as.numeric(stats$Value) # Decimal places can only be set with numeric/integer variables

labels(stats)<-list(
  Statistic="Variable",
  Value="Point Estimate",
  CI="95% Confidence Interval"
) #common package

t3<-
  stats |> gt() |> tab_header(
    title = md("**Table 3.** Confusion Table Statistics and Kappa")) |>
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_column_labels()) |>
  tab_footnote(
    footnote = "Not available.",
    locations = cells_body(
      columns = CI,
      rows = CI == "NA"
    ), placement = "right"
  ) |>
  fmt_number(columns = Value, decimals = 2, use_seps = F)

# ------------------------------- #
# 7. Radar Chart
# ------------------------------- #

# Enter chart data
chart_dat <- data.frame(
  row.names = c("Data"),
  Sensitivity = c(95),
  Specificity = c(42),
  PPV = c(79),
  NPV = c(80),
  Accuracy = c(79)
)

max_min <- data.frame(
  Sensitivity = c(100, 0), Specificity = c(100, 0),
  PPV = c(100, 0),NPV = c(100, 0), 
  Accuracy = c(100, 0)
)
rownames(max_min) <- c("Max", "Min")

rad <- rbind(max_min, chart_dat)

# Generate and customize the chart
radarchart(
  rad, axistype = 1, vlcex = 0.7,
  pcol = "lightgrey", pfcol = (alpha("lightgrey", 0.5)), plwd = 2.5, plty = 1,
  caxislabcol = "black",
  title = "Figure 1. Summary of Performance Metrics"
)

glcol = "black", cglty = 1, cglwd = 1,
  
# ------------------------------- #
# 8. Concordance Analysis
# ------------------------------- # 

# Import and prepare the datatset
loon_data <- read.csv("C:/All_FT_export.csv")

loon_data$Confidence <- factor(loon_data$Confidence, levels = c("VERY_HIGH", "HIGH", "MEDIUM", "LOW"))

# Create the Outcome Variable
loon_data$out <- ifelse(loon_data$Human==loon_data$Decision, 0, 1)
loon_data$out <- factor(loon_data$out, levels=c(1,0), 
                        labels=c("Disagree", "Agree"))
loon_data$out <- relevel(loon_data$out, ref="Agree")
labels(loon_data) <- list(out="Concordance")

chisq.test(loon_data$Confidence, loon_data$out)

table1(~ Confidence | out, loon_data,
       overall=F,
       caption="<b>Table 4:</b> FT Loon Lense Confidence Level by Human versus AI Concordance",
       footnote="<b>Chi-square test:</b> p < 0.0001")

ggplot(loon_data, aes(x=Confidence, fill=out)) +
  geom_bar(position="dodge") +
  scale_fill_manual(values=c("#104862", "#008b8b")) + 
  xlab("Confidence") + 
  ylab("Articles (n)") + 
  labs(title="**Figure 2:** Loon Lens Screening Concordance with Humans by Confidence Level at Full Text",
       subtitle= "Human versus AI Concordance") +
  scale_y_continuous(limits=c(0, 100), 
                     breaks=seq(0, 100, 20)) +
  theme_classic() +
  theme(plot.title.position="plot")
  
# End of Program #