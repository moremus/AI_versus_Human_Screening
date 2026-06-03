# ============================================================
#  Forest Plots – Figure 3 (Title & Abstract Screening)
#  AI vs. Human Screening in Systematic Reviews
#
#  Layout: 3 plots stacked in one column
#          Each plot has its own y-axis (full stat names)
#          Shared x-axis at bottom only
#          X-axis: 0.25 to 1.00 by 0.25
#          Labels (est + 95% CI) above each whisker
#          Boxes (point estimates) half-size
#          Font: Arial throughout
# ============================================================
#
#  Required packages:
#    install.packages(c("ggplot2", "patchwork", "extrafont"))
#    extrafont::font_import()   # run once
#    extrafont::loadfonts()
# ============================================================

library(ggplot2)
library(patchwork)

tryCatch({
  library(extrafont)
  loadfonts(quiet = TRUE)
}, error = function(e) message("extrafont not available – falling back to sans"))

arial <- "Arial"   # change to "sans" if Arial unavailable

# ── 1. DATA ─────────────────────────────────────────────────
# Full names, ordered top-to-bottom on y-axis
stat_levels <- c(
  "F1 Score",
  "Kappa",
  "Concordance",
  "Negative Predictive Value",
  "Positive Predictive Value",
  "Specificity",
  "Sensitivity"
)

make_df <- function(est, lo, hi) {
  data.frame(
    stat  = factor(stat_levels, levels = stat_levels),
    est   = est,
    lo    = lo,
    hi    = hi,
    label = sprintf("%.2f (%.2f, %.2f)", est, lo, hi)
  )
}

df_a <- make_df(
  est = c(0.50, 0.46, 0.93, 0.96, 0.52, 0.96, 0.49),
  lo  = c(0.39, 0.34, 0.91, 0.94, 0.38, 0.95, 0.37),
  hi  = c(0.61, 0.59, 0.94, 0.97, 0.66, 0.98, 0.63)
)

df_b <- make_df(
  est = c(0.76, 0.75, 0.97, 0.97, 0.91, 1.00, 0.66),
  lo  = c(0.66, 0.65, 0.96, 0.96, 0.81, 0.99, 0.54),
  hi  = c(0.84, 0.84, 0.98, 0.98, 0.98, 1.00, 0.77)
)

df_c <- make_df(
  est = c(0.47, 0.44, 0.93, 0.96, 0.55, 0.97, 0.41),
  lo  = c(0.34, 0.30, 0.91, 0.94, 0.39, 0.96, 0.29),
  hi  = c(0.59, 0.56, 0.95, 0.97, 0.69, 0.98, 0.56)
)

# ── 2. SHARED AXIS SETTINGS ─────────────────────────────────
x_min    <- 0.25
x_max    <- 1.05
x_breaks <- seq(0.25, 1.00, by = 0.25)
x_labels <- c("0.25", "0.50", "0.75", "1.00")

# ── 3. THEME ────────────────────────────────────────────────
forest_theme <- theme_classic(base_size = 11, base_family = arial) +
  theme(
    # y-axis: labels visible, no extra title
    axis.line.y   = element_line(colour = "black"),
    axis.ticks.y  = element_blank(),
    axis.title.y  = element_blank(),
    axis.text.y   = element_text(colour = "black", size = 10,
                                 hjust = 1, family = arial),
    # x-axis: suppress on top two panels; shown only on bottom panel
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.line.x   = element_blank(),
    # gridlines
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.4,
                                      linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_blank(),
    # tight bottom margin so panels stack close together
    plot.margin = margin(t = 6, r = 10, b = 2, l = 4),
    # tag (A / B / C)
    plot.tag          = element_text(face = "bold", size = 13,
                                     colour = "black", family = arial),
    plot.tag.position = c(0.01, 0.97),
    # subtitle used as panel title
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5,
                              colour = "black", family = arial)
  )

# Theme variant for the bottom panel: x-axis shown
forest_theme_bottom <- forest_theme +
  theme(
    axis.text.x  = element_text(colour = "black", size = 9, family = arial),
    axis.ticks.x = element_line(colour = "black"),
    axis.line.x  = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black", size = 10, family = arial,
                                margin = margin(t = 6)),
    plot.margin  = margin(t = 6, r = 10, b = 10, l = 4)
  )

# ── 4. PLOT FACTORY ─────────────────────────────────────────
make_forest <- function(df, title, tag, bottom = FALSE) {

  thm <- if (bottom) forest_theme_bottom else forest_theme

  ggplot(df, aes(x = est, y = stat)) +

    geom_vline(xintercept = 1.0, colour = "grey60",
               linewidth = 0.5, linetype = "solid") +

    geom_errorbarh(aes(xmin = lo, xmax = hi),
                   height = 0.20, colour = "black", linewidth = 0.55) +

    # half-size filled square
    geom_point(shape = 15, size = 1.6, colour = "black") +

    # label centred on point estimate, placed just below the whisker
    geom_text(
      aes(x = est, label = label),
      hjust    = 0.5,
      vjust    = -1.0,
      size     = 2.8,
      colour   = "black",
      family   = arial,
      fontface = "plain"
    ) +

    scale_x_continuous(
      limits = c(x_min, x_max),
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = expansion(add = 1.1)) +

    labs(
      title  = title,
      tag    = tag,
      x      = if (bottom) "Statistic (95% Confidence Interval)" else NULL
    ) +
    thm
}

# ── 5. BUILD PANELS ─────────────────────────────────────────
pA <- make_forest(df_a, "Catchii vs. Human",     tag = "A", bottom = FALSE)
pB <- make_forest(df_b, "Loon Lens vs. Human",   tag = "B", bottom = FALSE)
pC <- make_forest(df_c, "Loon Lens vs. Catchii", tag = "C", bottom = TRUE)

# ── 6. STACK WITH PATCHWORK ─────────────────────────────────
# plot_layout(axes = "collect") aligns all x-axes to the same scale;
# the x-axis labels appear only on the bottom panel (axis.text.x blank on A & B)
combined <- pA / pB / pC +
  plot_layout(heights = c(1, 1, 1))

# ── 7. SAVE ─────────────────────────────────────────────────
ggsave("Figure_3.png", combined,
       width = 8.5, height = 11, dpi = 900)

### END OF PROGRAM ###