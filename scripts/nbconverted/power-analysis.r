suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pwr))

get_effect_size <- function(group_a, group_b) {
    length_a <- length(group_a) - 1
    length_b <- length(group_b) - 1
    mean_difference <- abs(mean(group_a) - mean(group_b))
    
    common_std <- length_a * var(group_a) + length_b * var(group_b)
    common_std <- common_std / (length_a + length_b)
    effect_size <- mean_difference / sqrt(common_std)
    return(effect_size)
}

# Load Profiles
batch <- "2019_06_25_Batch3"
file <- file.path("data", paste0(batch, "_merged_intersected_variable_selected.csv"))

full_df <- readr::read_csv(file, col_types = readr::cols())
print(dim(full_df))
head(full_df, 3)

# Extract out mutant vs. wildtype clones
mut_df <- full_df %>% dplyr::filter(Metadata_Assay_Plate_Barcode == "MutClones")
wt_df <- full_df %>% dplyr::filter(Metadata_Assay_Plate_Barcode == "WTClones")

# Get all cp features
cp_features <- colnames(full_df %>% dplyr::select(-starts_with("Metadata_")))

testing_ns <- c(10, 25, 50, 75, 100, 250, 500, 750, 1000, 2000, 5000, 10000)
alpha <- 0.5 / length(cp_features)
current_n <- dim(full_df)[1] / 2

all_effect_size_results <- list()
idx <- 0
for (feature in cp_features) {
    group_a <- wt_df %>% dplyr::pull(!!feature)
    group_b <- mut_df %>% dplyr::pull(!!feature)
    feature_effect_size <- get_effect_size(group_a, group_b)
    
    for (example_n in testing_ns) {
        power_result <- pwr.t.test(n = example_n, d = feature_effect_size, sig.level = alpha)
        obs_power <- power_result$power
        result <- c(feature, example_n, obs_power)
        all_effect_size_results[[paste(idx)]] <- result
        idx <- idx + 1
    }

}

full_results_df <- all_effect_size_results %>%
    dplyr::bind_cols() %>%
    t() %>%
    tibble::as_tibble(.name_repair = "minimal")

colnames(full_results_df) <- c("feature", "n", "power")

full_results_df$n <- as.numeric(paste(full_results_df$n))
full_results_df$power <- as.numeric(paste(full_results_df$power))

print(dim(full_results_df))
head(full_results_df)

summary_df <- full_results_df %>% dplyr::group_by(n) %>%
    dplyr::mutate(mean_power = mean(power),
                  stderr_power = sd(power) / sqrt(length(power))) %>%
    dplyr::distinct(n, mean_power, stderr_power)
summary_df

rect <- data.frame(xmin=250, xmax=1000, ymin=-Inf, ymax=Inf)

power_gg <- ggplot(full_results_df) +
    geom_line(aes(x = n, y = power, group = feature),
              alpha = 0.02) +
    geom_point(aes(x = n, y = power, group = feature),
              alpha = 0.02) +
    geom_errorbar(data = summary_df,
                  color = "red",
                  width = 500,
                  aes(x = n,
                      ymin = mean_power - stderr_power,
                      ymax = mean_power + stderr_power)) +
    geom_line(data = summary_df,
              color = "red",
              aes(x = n, y = mean_power)) +
    geom_point(data = summary_df,
               color = "red",
               size = 2,
               shape = 21,
               fill = "white",
               aes(x = n, y = mean_power)) +
     geom_rect(data = rect,
               inherit.aes=FALSE,
               aes(xmin = xmin,
                   xmax = xmax,
                   ymin = ymin,
                   ymax = ymax),
               color="transparent",
               fill="green", alpha=0.3) +
    geom_vline(xintercept = current_n,
               linetype = "dashed",
               color = "blue") +
    geom_vline(xintercept = 250,
               linetype = "dashed",
               color = "green") +
    geom_vline(xintercept = 1000,
               linetype = "dashed",
               color = "green") +
    xlab("n =") +
    ylab("power") +
    theme_bw()

output_file <- file.path("figures", "power_analysis.png")
ggsave(output_file, dpi = 500, height = 6, width = 6)

power_gg
