# For consistent visualization colors, labsethemes

# Plotting options
legend_labels <- c(
  "training" = "Training",
  "test" = "Test",
  "validation" = "Validation",
  "holdout" = "Holdout"
)

# Color-blind friendly palette
legend_colors <- c(
  "training" = "#332288",
  "test" = "#DDCC77",
  "validation" = "#117733",
  "holdout" = "#CC6677"
)

# Color blind friendly palette for ground truth clones
clone_colors <- c(
  "CloneA" = "#E69F00",
  "CloneE" = "#D55E00",
  "WT_parental" = "#009E73"
)

clone_labels <- c(
  "CloneA" = "Clone A",
  "CloneE" = "Clone E",
  "WT_parental" = "WT Parental"
)

# To denote real from shuffled results
linetype_default <- c(
  False = "solid",
  True = "dotted"
)

linetype_labels <- c(
  False = "False",
  True = "True"
)

# For a consistent custom theme
custom_theme <- ggplot2::theme_bw() + ggplot2::theme(
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 4),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    strip.text = element_text(size = 5),
    strip.background = element_rect(size = 0.5, colour = "black", fill = "#fdfff4"),
    legend.key.size = unit(2, "mm"),
    legend.key.width = unit(5, "mm"),
    legend.margin = margin(l = -2),
    panel.grid.major = element_line(size = 0.25),
    panel.grid.minor = element_line(size = 0.25)
)
