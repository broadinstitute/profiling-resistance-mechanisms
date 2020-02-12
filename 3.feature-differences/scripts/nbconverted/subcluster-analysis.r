suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(cowplot))

set.seed(0112)

con <- curl::curl("https://raw.githubusercontent.com/broadinstitute/cytominer_scripts/master/write_gct.R")
source(con)

util_file <- file.path("..", "scripts", "processing_utils.R")
source(util_file)

viz_file <- file.path("..", "scripts", "visualization_utils.R")
source(viz_file)

file <- file.path("..", "data", "2019_06_25_Batch3_merged_intersected_variable_selected.csv")
data_df <- load_data(file)

print(dim(data_df))
head(data_df, 3)

feature_df <- data_df %>% dplyr::select(-starts_with("Metadata_"))

cor_mat <- cor(t(feature_df), method = "pearson")
dist_mat <- as.dist(cor_mat)
kmeans_result <- kmeans(dist_mat, centers = 4)

data_df <- data_df %>%
    dplyr::mutate(Metadata_cluster_groups = paste(kmeans_result$cluster))

table(kmeans_result$cluster)

gct_file <- file.path("results", "kmeans_batch3_group.gct")
write_gct(data_df, gct_file)

tstat_group <- c()
pval_group <- c()
subgroup_vector <- c()
feature_vector <- c()
shuffle_vector <- c()

for (subgroup in unique(data_df$Metadata_cluster_groups)) {
    
    data_subset_df <- data_df %>%
        dplyr::filter(Metadata_cluster_groups == !!subgroup)
    
    for (feature in colnames(data_subset_df)) {

        if (!grepl("Metadata_", feature)) {
            
            for (shuffle in c("no_shuffle", "shuffle")) {
                
                if (shuffle == "shuffle") {
                    data_subset_copy_df <- data_subset_df
                    data_subset_copy_df$Metadata_Plate <- sample(data_subset_copy_df$Metadata_Plate)
                } else {
                    data_subset_copy_df <- data_subset_df
                }
                # Perform Cell Line Experiment at 0.7uM Dose
                resistant_clones <- data_subset_copy_df %>%
                    dplyr::filter(Metadata_Plate == "MutClones") %>%
                    dplyr::pull(!!feature)

                wt_cells <- data_subset_copy_df %>%
                    dplyr::filter(Metadata_Plate == "WTClones") %>%
                    dplyr::pull(!!feature)

                result <- t.test(resistant_clones, wt_cells, var.equal = FALSE)

                tstat_group <- c(tstat_group, as.numeric(paste(result$statistic)))
                pval_group <- c(pval_group, result$p.value)
                subgroup_vector <- c(subgroup_vector, subgroup)
                feature_vector <- c(feature_vector, feature)
                shuffle_vector <- c(shuffle_vector, shuffle)
            }
            
        }
    } 
}

result_df <- as.data.frame(cbind(tstat_group, pval_group, subgroup_vector, feature_vector, shuffle_vector))

result_df$tstat_group <- as.numeric(paste(result_df$tstat_group))
result_df$pval_group <- as.numeric(paste(result_df$pval_group))
result_df$neglog10p <- -log10(result_df$pval_group)

result_file <- file.path("results", "subgroup_ttest_results.tsv")
result_df %>% readr::write_tsv(result_file)

head(result_df, 3)

yintercept = -log10(0.05 / (dim(result_df)[1] / 2))
yintercept

repel_logic <- result_df$neglog10p > yintercept * 1.4

subcluster_colors <- c("1" = "#d3aa5a",
                       "2" = "#1d6dab",
                       "3" = "#a9da7f", 
                       "4" = "#2e9627")

volcano_gg <- ggplot(result_df,
       aes(x = tstat_group,
           y = neglog10p)) +
    geom_point(aes(color = subgroup_vector),
               alpha = 0.8,
               size = .8) +
    geom_hline(yintercept = yintercept,
               color = "red",
               linetype = "dashed") +
    xlab("t Statistic") +
    ylab("-log10 P") +
    facet_wrap(~shuffle_vector) +
    scale_color_manual(name = "subgroup",
                       values = subcluster_colors,
                       labels = c("1" = "1", 
                                  "2" = "2",
                                  "3" = "3",
                                  "4" = "4")) +
    geom_text_repel(data = subset(result_df, repel_logic),
                    arrow = arrow(length = unit(0.01, "npc")),
                    size = 2,
                    segment.size = 0.1,
                    segment.alpha = 0.8,
                    force = 20,
                    aes(x = tstat_group,
                        y = neglog10p,
                        label = feature_vector)) +
    theme_bw() +
    theme(strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4"))

file_base <- file.path("figures", "subcluster_volcano_plot")
save_figure(volcano_gg, file_base, height = 6, width = 12)

volcano_gg

result_df$mutant_high <- as.integer(result_df$tstat_group > 0)

sorted_result_df <- result_df %>%
    dplyr::group_by(feature_vector, shuffle_vector) %>%
    dplyr::mutate(feature_sum = sum(neglog10p),
                  same_direction_tstat = sum(mutant_high)) %>%
    dplyr::arrange(desc(feature_sum)) %>%
    dplyr::ungroup()

sorted_result_df$feature_vector <- factor(sorted_result_df$feature_vector,
                                          levels = unique(sorted_result_df$feature_vector))


head(sorted_result_df, 10)

all_feature_gg <- ggplot(sorted_result_df,
       aes(x = feature_vector,
           y = neglog10p)) +
    geom_bar(aes(fill = subgroup_vector), stat="identity") +
    theme_bw() + 
    facet_wrap(~shuffle_vector, ncol = 1) +
    theme(axis.text.x = element_text(size = 5, angle = 90),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    scale_fill_manual(name = "subgroup",
                       values = subcluster_colors,
                       labels = c("1" = "1", 
                                  "2" = "2",
                                  "3" = "3",
                                  "4" = "4"))

file_base <- file.path("figures", "subcluster_all_feature_barchart")
save_figure(all_feature_gg, file_base, height = 10, width = 12)

all_feature_gg

all_tstat_feature_gg <- ggplot(sorted_result_df,
       aes(x = feature_vector,
           y = tstat_group)) +
    geom_bar(aes(fill = subgroup_vector), stat="identity") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 5, angle = 90),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    facet_wrap(~shuffle_vector, ncol = 1) +
    scale_fill_manual(name = "subgroup",
                       values = subcluster_colors,
                       labels = c("1" = "1", 
                                  "2" = "2",
                                  "3" = "3",
                                  "4" = "4")) +
    geom_hline(yintercept = 0,
                   linetype = "dashed",
                   lwd = 0.5)

file_base <- file.path("figures", "subcluster_all_feature_barchart_tstat")
save_figure(all_tstat_feature_gg, file_base, height = 6, width = 12)

all_tstat_feature_gg

top_result_df <- sorted_result_df %>%
    dplyr::filter(feature_sum > 9)

top_feature_gg <- ggplot(top_result_df,
       aes(x = feature_vector,
           y = neglog10p)) +
    geom_bar(aes(fill = subgroup_vector), stat="identity") +
    theme_bw() + 
    facet_wrap(~shuffle_vector, nrow = 1) +
    theme(axis.text.y = element_text(size = 10, angle = 0),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    scale_fill_manual(name = "subgroup",
                       values = subcluster_colors,
                       labels = c("1" = "1", 
                                  "2" = "2",
                                  "3" = "3",
                                  "4" = "4")) +
    coord_flip()

file_base <- file.path("figures", "subcluster_top_feature_barchart")
save_figure(top_feature_gg, file_base, height = 6, width = 9)

top_feature_gg

tstat_feature_gg <- ggplot(top_result_df,
       aes(x = feature_vector,
           y = tstat_group)) +
    geom_bar(aes(fill = subgroup_vector), stat="identity") +
    theme_bw() + 
    facet_wrap(~shuffle_vector, nrow = 1) +
    theme(axis.text.y = element_text(size = 10, angle = 0),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    scale_fill_manual(name = "subgroup",
                      values = subcluster_colors,
                      labels = c("1" = "1", 
                                  "2" = "2",
                                  "3" = "3",
                                  "4" = "4")) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               lwd = 0.5) +
    coord_flip()

file_base <- file.path("figures", "subcluster_top_feature_barchart_tstat")
save_figure(tstat_feature_gg, file_base, height = 6, width = 9)

tstat_feature_gg

same_direction_df <- sorted_result_df %>%
    dplyr::filter(same_direction_tstat %in% c(0, 4))

tstat_feature_same_direction_gg <- ggplot(same_direction_df,
       aes(x = feature_vector,
           y = tstat_group)) +
    geom_bar(aes(fill = subgroup_vector), stat="identity") +
    theme_bw() + 
    facet_wrap(~shuffle_vector, nrow = 1) +
    theme(axis.text.y = element_text(size = 10, angle = 0),
          strip.background = element_rect(colour = "black",
                                          fill = "#fdfff4")) +
    scale_fill_manual(name = "subgroup",
                       values = subcluster_colors,
                       labels = c("1" = "1", 
                                  "2" = "2",
                                  "3" = "3",
                                  "4" = "4")) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               lwd = 0.5) +
    coord_flip()

file_base <- file.path("figures", "subcluster_top_feature_barchart_tstat_samedirection")
save_figure(tstat_feature_same_direction_gg, file_base, height = 6, width = 9)

tstat_feature_same_direction_gg
