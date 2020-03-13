# Functions to facilitate visualization
#
# Usage:
# Import Only

# Load Plotting Function for ttest Volcano Plots
ttest_volcano <- function(df,
                          x_string,
                          title,
                          yintercept,
                          repel_logic,
                          ggrepel_label_size,
                          title_text_size,
                          axis_text_size,
                          axis_title_size,
  
                          ymax = 10) {
  # Plot the results of a t-test on various samples for common features
  #
  # Arguments:
  # df - the dataframe storing the t-test results
  # x_string - string indicating what variable to plot on x axis
  # title - string indicating the title of the plot
  # yintercept - an alpha corrected value to plot a red dotted line
  # repel_logic - which features to highlight and name
  # ggrepel_label_size - int of the size of the feature names
  # title_text_size - int of the size of the title text
  # axis_text_size - int of the size of the text on the ggplot axes
  # axis_title_size - int of the size of the titles on the ggplot axes
  # ymax - int indicating the maximum height of the y axis
  #
  # Output:
  # The ggplot2 object for downstream saving
  ttest_gg <- ggplot(df,
                     aes_string(x = x_string,
                                y = "neglog10p")) +
    geom_point(alpha = 0.5,
               size = 0.8,
               color = ifelse(repel_logic, "red", "grey50")) +
    geom_hline(yintercept = yintercept,
               color = "red",
               linetype = "dashed") +
    xlab("t Statistic") +
    ylab("-log10 P") +
    ggtitle(title) +
    ylim(c(0, ymax)) +
    geom_text_repel(data = subset(df, repel_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    size = ggrepel_label_size,
                    segment.size = 0.1,
                    segment.alpha = 0.8,
                    force = 20,
                    aes_string(x = x_string,
                               y = "neglog10p",
                               label = "feature")) +
    theme_bw() +
    theme(plot.title = element_text(size = title_text_size),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size))

  return(ttest_gg)
}


save_figure <- function(main_figure,
                        file_base,
                        extensions = c(".png", ".pdf", ".svg"),
                        height = 6,
                        width = 8) {
  # Save figure given extensions
  #
  # Arguments:
  # main_figure - the cowplot or ggplot object
  # file_base - the name of the file without extensions
  # extensions - character vector of file extensions to save
  # height - height of plot
  # width - width of plot
  #
  # Output:
  # Will save plots to file

  for (extension in extensions) {
    ggsave(main_figure,
           filename = paste0(file_base, extension),
           height = height,
           width = width)
  }
}
