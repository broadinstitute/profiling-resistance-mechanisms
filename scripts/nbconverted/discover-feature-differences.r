
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))

doses <- c(0.7, 7)

axis_title_size <- 12
axis_text_size <- 10
strip_text_size <- 8
ggrepel_label_size <- 2.4

# Set column types for reading in data
batch_cols = readr::cols(
    .default = readr::col_double(),
    Metadata_Plate = readr::col_character(),
    Metadata_Well = readr::col_character(),
    Metadata_Assay_Plate_Barcode = readr::col_character(),
    Metadata_Plate_Map_Name = readr::col_character(),
    Metadata_Batch_Number = readr::col_integer(),
    Metadata_well_position = readr::col_character(),
    Metadata_CellLine = readr::col_character()
)

file <- file.path("data", "merged_intersected_variable_selected.csv")
data_df <- readr::read_csv(file, col_types = batch_cols) 

print(dim(data_df))
head(data_df, 3)

tstat <- c()
pval <- c()
all_features <- c()
for (feature in colnames(data_df)) {
    if (!grepl("Metadata_", feature)) {
        doseA <- data_df %>%
            dplyr::filter(Metadata_Dosage == doses[1],
                          Metadata_CellLine != "WT") %>%
            dplyr::pull(!!feature)
        
        doseB <- data_df %>%
            dplyr::filter(Metadata_Dosage == doses[2],
                          Metadata_CellLine != "WT") %>%
            dplyr::pull(!!feature)
        
        result <- t.test(doseA, doseB, var.equal = FALSE)
        
        all_features <- c(all_features, feature)
        tstat <- c(tstat, as.numeric(paste(result$statistic)))
        pval <- c(pval, result$p.value)
    }
}

result_df <- as.data.frame(cbind(tstat, pval))
result_df$neglog10p <- -log10(result_df$pval)
result_df$feature <- all_features
result_df <- result_df %>% dplyr::arrange(desc(neglog10p))

dim(result_df)
head(result_df, 3)

alpha_correction <- -log10(0.05 / dim(result_df)[1])
repel_logic <- result_df$neglog10p > alpha_correction * 1.5

ttest_gg <- ggplot(result_df, aes(x = tstat, y = neglog10p)) +
    geom_point(alpha = 0.5,
               color = ifelse(repel_logic, "red", "grey50")) +
    geom_hline(yintercept = alpha_correction,
               color = 'red',
               linetype = 'dashed') +
    xlab("t Statistic") +
    ylab("-log10 P") +
    geom_label_repel(data = subset(result_df, repel_logic),
                     arrow = arrow(length = unit(0.01, "npc")),
                     size = ggrepel_label_size,
                     segment.size = 0.1,
                     segment.alpha = 0.8,
                     force = 20,
                     aes(x = tstat, y = neglog10p, label = feature)) +
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size))

ttest_gg

# The top feature is something to do with nuclear area
top_feature <- paste(result_df$feature[1])

# top_feature <- paste(result_df$feature[3])

append_batch <- function(string) paste("Batch:", string)
append_dose <- function(string) paste("Dose:", string)

distrib_gg <- ggplot(data_df, aes_string(x = `top_feature`)) +
    geom_density(aes(fill = Metadata_CellLine),
                 alpha = 0.6) +
    geom_rug(aes(color = Metadata_CellLine),
             alpha = 0.8,
             size = 0.5) +
    facet_grid(Metadata_Dosage ~ Metadata_Batch_Number,
               scales = "free_y",
               labeller = labeller(Metadata_Batch_Number = as_labeller(append_batch),
                                   Metadata_Dosage = as_labeller(append_dose))) +
    scale_color_manual(name = "Cell Line",
                      labels = c("CloneA" = "Clone A",
                                 "CloneE" = "Clone E",
                                 "WT" = "Wild type"),
                      values = c("CloneA" = "#1b9e77",
                                 "CloneE" = "#d95f02",
                                 "WT" = "#7570b3")) +
    scale_fill_manual(name = "Cell Line",
                      labels = c("CloneA" = "Clone A",
                                 "CloneE" = "Clone E",
                                 "WT" = "Wild type"),
                      values = c("CloneA" = "#1b9e77",
                                 "CloneE" = "#d95f02",
                                 "WT" = "#7570b3")) +
    theme_bw() +
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          strip.text = element_text(size = strip_text_size),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

distrib_gg

table(data_df$Metadata_CellLine, data_df$Metadata_Dosage, data_df$Metadata_Batch_Number)

main_plot <- (
    cowplot::plot_grid(
        ttest_gg,
        distrib_gg,
        labels = c("a", "b"),
        ncol = 2,
        nrow = 1
    )
)

main_plot

file_base <- file.path("figures", "dosage_feature_figure")
for (extension in c('.png', '.pdf')) {
    ggsave(main_plot, filename = paste0(file_base, extension), height = 5, width = 10)
}
