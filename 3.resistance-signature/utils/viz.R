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

# To denote resistance vs. sensitive
resistance_status_colors <- c(
    "resistant" = "#E66100",
    "sensitive" = "#5D3A9B"
)

resistance_status_labels <- c(
    "resistant" = "Resistant",
    "sensitive" = "Sensitive"
)

# For a consistent custom theme
custom_theme <- ggplot2::theme_bw() + ggplot2::theme(
    legend.title = ggplot2::element_text(size = 5),
    legend.text = ggplot2::element_text(size = 4),
    axis.title = ggplot2::element_text(size = 6),
    axis.text = ggplot2::element_text(size = 5),
    strip.text = ggplot2::element_text(size = 5),
    strip.background = ggplot2::element_rect(
      size = 0.5, colour = "black", fill = "#fdfff4"
    ),
    legend.key.size = ggplot2::unit(2, "mm"),
    legend.key.width = ggplot2::unit(5, "mm"),
    legend.margin = ggplot2::margin(l = -2),
    panel.grid.major = ggplot2::element_line(size = 0.25),
    panel.grid.minor = ggplot2::element_line(size = 0.25)
)
