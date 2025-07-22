# Load required packages with suppressed startup messages
suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  
  # CRAN packages
  library(ggpubr)
  library(struct)
  library(reshape2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
  library(openxlsx)  
  library(janitor)   
  library(rstatix)   
})

#---------------------------------------------------------------
# Load sample metadata (common to all datasets)
sample_meta <- read.csv("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/data/sample_meta.csv")

#---------------------------------------------------------------
# Create DatasetExperiment objects for each data type

create_DE <- function(data, variable_meta, name) {
  DatasetExperiment(
    data = t(data),
    sample_meta = sample_meta,
    variable_meta = as.data.frame(variable_meta, row.names = colnames(t(data))),
    description = 'Metabolomics Testing',
    name = name
  )
}

DE_Hilic_Pos <- create_DE(data_hp, variable_meta_hp, "Hilic ESI+")
DE_Hilic_Neg <- create_DE(data_hn, variable_meta_hn, "Hilic ESI-")
DE_C18_Pos    <- create_DE(data_rpp, variable_meta_rpp, "C18 ESI+")
DE_C18_Neg    <- create_DE(data_rpn, variable_meta_rpn, "C18 ESI-")

#---------------------------------------------------------------
# Convert relevant sample_meta columns to factors for all datasets
factor_vars <- c("batch", "type", "class", "extraction_type")
datasets <- list(DE_Hilic_Pos, DE_Hilic_Neg, DE_C18_Pos, DE_C18_Neg)

for (ds in datasets) {
  for (var in factor_vars) {
    ds$sample_meta[[var]] <- factor(ds$sample_meta[[var]])
  }
}

#---------------------------------------------------------------
# Count initial features in each dataset (number of columns in data matrix)
count_pre <- sapply(datasets, function(ds) ncol(ds$data))
names(count_pre) <- c("Hilic_Pos", "Hilic_Neg", "C18_Pos", "C18_Neg")

#---------------------------------------------------------------
# Setup data cleaning filters
blk_filter <- blank_filter(
  fold_change = 15,
  blank_label = "Blank",
  qc_label = "QC",
  factor_name = "type"
)

perc_features <- mv_feature_filter(
  threshold = 60,
  method = "across",
  factor_name = "type"
)

qc_features <- rsd_filter(
  rsd_threshold = 30,
  qc_label = "QC",
  factor_name = "type"
)

#---------------------------------------------------------------
# Function to apply filters sequentially to a DatasetExperiment
apply_filters <- function(DE) {
  DE_filtered <- model_apply(blk_filter, DE) |> predicted()
  DE_filtered <- model_apply(perc_features, DE_filtered) |> predicted()
  DE_filtered <- model_apply(qc_features, DE_filtered) |> predicted()
  DE_filtered
}

# Apply filters to all datasets
Hilic_Pos <- apply_filters(DE_Hilic_Pos)
Hilic_Neg <- apply_filters(DE_Hilic_Neg)
C18_Pos   <- apply_filters(DE_C18_Pos)
C18_Neg   <- apply_filters(DE_C18_Neg)

#---------------------------------------------------------------
# Function to calculate %RSD (relative standard deviation)
calc_rsd <- function(x) {
  sd(x) * 100 / mean(x)
}

#---------------------------------------------------------------
# Function to process dataset for RSD calculation and plotting
process_rsd_data <- function(DE_filtered, title) {
  data <- DE_filtered$data
  data$class <- DE_filtered$sample_meta$class
  
  # Gather data to long format for groupwise RSD calculation
  result <- data %>%
    gather(key = "Variable", value = "Value", -class) %>%
    group_by(class, Variable) %>%
    summarise(RSD = calc_rsd(Value), .groups = "drop")
  
  # Density plot of RSD by class
  density_plot <- ggplot(result, aes(x = RSD, fill = class, color = class)) +
    geom_density(alpha = 0.1, na.rm = TRUE) +
    theme_bw(20) +
    labs(x = "RSD (%)", title = title) +
    scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
    scale_color_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
    facet_grid(~ class)
  
  # Boxplot of RSD by class
  box_plot <- ggplot(result, aes(x = class, y = RSD, fill = class)) +
    geom_boxplot(alpha = 0.5) +
    theme_bw(20) +
    labs(y = "RSD (%)", title = title) +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = rel(0.8), angle = 45, hjust = 1, vjust = 1.1)
    ) +
    scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black"))
  
  # Reshape RSD results for exporting
  rsd_table <- dcast(result, class ~ Variable, value.var = "RSD")
  rsd_table_t <- as.data.frame(t(rsd_table))
  rsd_table_t <- rsd_table_t %>% janitor::row_to_names(row_number = 1)
  rsd_table_t <- sapply(rsd_table_t, as.numeric)
  rsd_table_t <- as_tibble(rsd_table_t, .name_repair = "unique")
  
  list(
    density_plot = density_plot,
    box_plot = box_plot,
    rsd_data = rsd_table_t
  )
}

#---------------------------------------------------------------
# Process each dataset
hp_results <- process_rsd_data(Hilic_Pos, "HILIC ESI+")
hn_results <- process_rsd_data(Hilic_Neg, "HILIC ESI-")
rpp_results <- process_rsd_data(C18_Pos, "C18 ESI+")
rpn_results <- process_rsd_data(C18_Neg, "C18 ESI-")

#---------------------------------------------------------------
# Summary stats example (adjust variables accordingly)
hp_summary <- hp_results$rsd_data %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")

hn_summary <- hn_results$rsd_data %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")

rpp_summary <- rpp_results$rsd_data %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")

rpn_summary <- rpn_results$rsd_data %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")

#---------------------------------------------------------------
# Combine density plots vertically and save
(rp <- (hp_results$density_plot / hn_results$density_plot / rpp_results$density_plot / rpn_results$density_plot) +
   plot_annotation(tag_levels = 'A') +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom"))

ggsave(
  filename = "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/rsd_density.pdf",
  plot = rp,
  width = 10,
  height = 15
)

#---------------------------------------------------------------
# Combine box plots and additional TIC plots (replace tic_* with your existing plots)
box_plots <- (hp_results$box_plot | hn_results$box_plot) / (rpp_results$box_plot | rpn_results$box_plot)

# Assuming you have objects tic_Hilic_Pos_pre, tic_Hilic_Neg_pre, etc., and post_tic_* defined similarly
combined_plot <- box_plots /
  (tic_Hilic_Pos_pre | tic_Hilic_Neg_pre | tic_C18_Pos_pre | tic_C18_Neg_pre) /
  (post_tic_Hilic_Pos | post_tic_Hilic_Neg | post_tic_C18_Pos | post_tic_C18_Neg) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 10))

ggsave(
  filename = "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/tic_processed.pdf",
  plot = combined_plot,
  width = 12,
  height = 10.5
)

#---------------------------------------------------------------
# Export RSD results to Excel file with multiple sheets
data_frames <- list(
  "RSD HILIC ESI+" = hp_results$rsd_data,
  "RSD HILIC ESI-" = hn_results$rsd_data,
  "RSD C18 ESI+" = rpp_results$rsd_data,
  "RSD C18 ESI-" = rpn_results$rsd_data
)

write.xlsx(
  data_frames,
  file = "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/Supplementary Table S2.xlsx"
)

#---------------------------------------------------------------
# Optional: plot boxplots grouped in a 2x2 grid for comparison
(rsd_plot_grid <- (hp_results$box_plot | hn_results$box_plot) / (rpp_results$box_plot | rpn_results$box_plot) +
   plot_annotation(tag_levels = 'A') +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom", legend.text = element_text(size = 10)))
