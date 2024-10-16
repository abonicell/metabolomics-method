suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  
  # CRAN libraries
  library(ggpubr)
  library(struct)
  library(reshape2)
  library(dplyr)
  library(patchwork)
  library(tidyr)
})

## ---------------------------------------------------------------
# metadata are the same for all datasets
sample_meta <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/data/sample_meta.csv"
)

# create DatasetExperiment ESI+
DE_Hilic_Pos <- DatasetExperiment(
  data = t(data_hp),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_hp, row.names = colnames(t(data_hp))),
  description = 'Metabolomics Testing',
  name = "Hilic ESI+"
)

# create DatasetExperiment ESI-
DE_Hilic_Neg <- DatasetExperiment(
  data = t(data_hn),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_hn, row.names = colnames(t(data_hn))),
  description = 'Metabolomics Testing',
  name = "Hilic ESI+"
)

# create DatasetExperiment C18+
DE_C18_Pos <- DatasetExperiment(
  data = t(data_rpp),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_rpp, row.names = colnames(t(data_rpp))),
  description = 'Metabolomics Testing',
  name = "Hilic ESI+"
)

# create DatasetExperiment C18-
DE_C18_Neg <- DatasetExperiment(
  data = t(data_rpn),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_rpn, row.names = colnames(t(data_rpn))),
  description = 'Metabolomics Testing',
  name = "Hilic ESI+"
)

# convert to factors
DE_Hilic_Pos$sample_meta$batch = factor(DE_Hilic_Pos$sample_meta$batch)
DE_Hilic_Pos$sample_meta$type = factor(DE_Hilic_Pos$sample_meta$type)
DE_Hilic_Pos$sample_meta$class = factor(DE_Hilic_Pos$sample_meta$class)
DE_Hilic_Pos$sample_meta$extraction_type = factor(DE_Hilic_Pos$sample_meta$extraction_type)

DE_Hilic_Neg$sample_meta$batch = factor(DE_Hilic_Neg$sample_meta$batch)
DE_Hilic_Neg$sample_meta$type = factor(DE_Hilic_Neg$sample_meta$type)
DE_Hilic_Neg$sample_meta$class = factor(DE_Hilic_Neg$sample_meta$class)
DE_Hilic_Neg$sample_meta$extraction_type = factor(DE_Hilic_Neg$sample_meta$extraction_type)

DE_C18_Pos$sample_meta$batch = factor(DE_C18_Pos$sample_meta$batch)
DE_C18_Pos$sample_meta$type = factor(DE_C18_Pos$sample_meta$type)
DE_C18_Pos$sample_meta$class = factor(DE_C18_Pos$sample_meta$class)
DE_C18_Pos$sample_meta$extraction_type = factor(DE_C18_Pos$sample_meta$extraction_type)

DE_C18_Neg$sample_meta$batch = factor(DE_C18_Neg$sample_meta$batch)
DE_C18_Neg$sample_meta$type = factor(DE_C18_Neg$sample_meta$type)
DE_C18_Neg$sample_meta$class = factor(DE_C18_Neg$sample_meta$class)
DE_C18_Neg$sample_meta$extraction_type = factor(DE_C18_Neg$sample_meta$extraction_type)

hp_pre_count <- ncol(DE_Hilic_Pos$data)
hn_pre_count <- ncol(DE_Hilic_Neg$data)
rpp_pre_count <- ncol(DE_C18_Pos$data)
rpn_pre_count <- ncol(DE_C18_Neg$data)

count_pre <- cbind(hp_pre_count, hn_pre_count, rpp_pre_count, rpn_pre_count)

## ---------------------------------------------------------------
# matrix processing - set up cleaning steps
# blank filter
blk_filter <- blank_filter(
  fold_change = 15,
  blank_label = "Blank",
  qc_label = "QC",
  factor_name = 'type'
)

# missing feature and % value
perc_features <-
  mv_feature_filter(threshold = 60,
                    method = 'across',
                    factor_name = 'type')

# rsd qc by feature
qc_features <-
  rsd_filter(rsd_threshold = 30,
             qc_label = 'QC',
             factor_name = 'type')

##-----------------------------------------------------------------
# apply blank filter
blk_Hilic_Pos <- model_apply(blk_filter, DE_Hilic_Pos)
blk_Hilic_Pos <- predicted(blk_Hilic_Pos)

blk_Hilic_Neg <- model_apply(blk_filter, DE_Hilic_Neg)
blk_Hilic_Neg <- predicted(blk_Hilic_Neg)

blk_C18_Pos <- model_apply(blk_filter, DE_C18_Pos)
blk_C18_Pos <- predicted(blk_C18_Pos)

blk_C18_Neg <- model_apply(blk_filter, DE_C18_Neg)
blk_C18_Neg <- predicted(blk_C18_Neg)

# apply missing data filter
blk_perc_Hilic_Pos <- model_apply(perc_features, blk_Hilic_Pos)
blk_perc_Hilic_Pos <- predicted(blk_perc_Hilic_Pos)

blk_perc_Hilic_Neg <- model_apply(perc_features, blk_Hilic_Neg)
blk_perc_Hilic_Neg <- predicted(blk_perc_Hilic_Neg)

blk_perc_C18_Pos <- model_apply(perc_features, blk_C18_Pos)
blk_perc_C18_Pos <- predicted(blk_perc_C18_Pos)

blk_perc_C18_Neg <- model_apply(perc_features, blk_C18_Neg)
blk_perc_C18_Neg <- predicted(blk_perc_C18_Neg)

# apply rsd qc by feature filter
blk_perc_qc_Hilic_Pos <- model_apply(qc_features, blk_perc_Hilic_Pos)
Hilic_Pos <- predicted(blk_perc_qc_Hilic_Pos)

blk_perc_qc_Hilic_Neg <- model_apply(qc_features, blk_perc_Hilic_Neg)
Hilic_Neg <- predicted(blk_perc_qc_Hilic_Neg)

blk_perc_qc_C18_Pos <- model_apply(qc_features, blk_perc_C18_Pos)
C18_Pos <- predicted(blk_perc_qc_C18_Pos)

blk_perc_qc_C18_Neg <- model_apply(qc_features, blk_perc_C18_Neg)
C18_Neg <- predicted(blk_perc_qc_C18_Neg)

Hilic_Pos
Hilic_Neg
C18_Pos
C18_Neg

##------------------------------------------------------------------------------
# Function to calculate %RSD
calc_rsd <- function(x) {
  sd_value <- sd(x)
  mean_value <- mean(x)
  rsd <- (sd_value * 100 / mean_value)
  return(rsd)
}

##------------------------------------------------------------------------------
# Hilic pos
hp_dat <- Hilic_Pos$data
hp_dat$class <- Hilic_Pos$sample_meta$class
hp_dat_melt <- melt(hp_dat)

# Using dplyr and tidyr to calculate %RSD across groups for multiple columns
result_hp <- hp_dat %>%
  gather(key = "Variable", value = "Value", -class) %>%
  group_by(class, Variable) %>%
  summarise(RSD = rsd(Value)) %>%
  ungroup()

rsd_hp_density <- ggplot(result_hp, aes(x = RSD, fill = class, color = class), na.rm = TRUE) +
  geom_density(alpha = 0.1) + theme_bw(20) +
  labs(x = "RSD (%)", title = "HILIC ESI+") +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  scale_color_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  facet_grid(~ class)

rsd_hp_box <- ggplot(result_hp, aes(x = class, y = RSD, fill = class)) +
  geom_boxplot(alpha = 0.5) + theme_bw(20) +
  labs(y = "RSD (%)", title = "HILIC ESI+") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  ) +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black"))

rsd_results_hp <- dcast(result_hp, class ~ Variable, value.var = "RSD")
rsd_results_hp_final <- dcast(result_hp, class ~ Variable, value.var = "RSD")
rsd_results_hp <- as.data.frame(t(rsd_results_hp))
rsd_results_hp <- rsd_results_hp %>% janitor::row_to_names(row_number = 1)
rsd_results_hp <- sapply(rsd_results_hp, as.numeric)
rsd_results_hp <- as_tibble(rsd_results_hp, .name_repair = "unique")



##------------------------------------------------------------------------------
# hilic neg
hn_dat <- Hilic_Neg$data
hn_dat$class <- Hilic_Neg$sample_meta$class
hn_dat_melt <- melt(hn_dat)

# Using dplyr and tidyr to calculate %RSD across groups for multiple columns
result_hn <- hn_dat %>%
  gather(key = "Variable", value = "Value", -class) %>%
  group_by(class, Variable) %>%
  summarise(RSD = rsd(Value)) %>%
  ungroup()

rsd_hn_density <- ggplot(result_hn, aes(x = RSD, fill = class, color = class), na.rm = TRUE) +
  geom_density(alpha = 0.1) + theme_bw(20) +
  labs(x = "RSD (%)", title = "HILIC ESI-") +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  scale_color_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  facet_grid(~ class)

rsd_hn_box <- ggplot(result_hn, aes(x = class, y = RSD, fill = class)) +
  geom_boxplot(alpha = 0.5) + theme_bw(20) +
  labs(y = "RSD (%)", title = "HILIC ESI-") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  ) +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black"))

rsd_results_hn <- dcast(result_hn, class ~ Variable, value.var = "RSD")
rsd_results_hn_final <- dcast(result_hn, class ~ Variable, value.var = "RSD")
rsd_results_hn <- as.data.frame(t(rsd_results_hn))
rsd_results_hn <- rsd_results_hn %>% janitor::row_to_names(row_number = 1)
rsd_results_hn <- sapply(rsd_results_hn, as.numeric)
rsd_results_hn <- as_tibble(rsd_results_hn, .name_repair = "unique")

##------------------------------------------------------------------------------
# C18 pos
rpp_dat <- C18_Pos$data
rpp_dat$class <- C18_Pos$sample_meta$class
rpp_dat_melt <- melt(rpp_dat)

# Using dplyr and tidyr to calculate %RSD across groups for multiple columns
result_rpp <- rpp_dat %>%
  gather(key = "Variable", value = "Value", -class) %>%
  group_by(class, Variable) %>%
  summarise(RSD = rsd(Value)) %>%
  ungroup()

rsd_rpp_density <- ggplot(result_rpp, aes(x = RSD, fill = class, color =
                                            class), na.rm = TRUE) +
  geom_density(alpha = 0.1) + theme_bw(20) +
  labs(x = "RSD (%)", title = "C18 ESI+") +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  scale_color_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  facet_grid(~ class)

rsd_rpp_box <- ggplot(result_rpp, aes(x = class, y = RSD, fill = class)) +
  geom_boxplot(alpha = 0.5) + theme_bw(20) +
  labs(y = "RSD (%)", title = "C18 ESI+") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  ) +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black"))

rsd_results_rpp <- dcast(result_rpp, class ~ Variable, value.var = "RSD")
rsd_results_rpp_final <- dcast(result_rpp, class ~ Variable, value.var = "RSD")
rsd_results_rpp <- as.data.frame(t(rsd_results_rpp))
rsd_results_rpp <- rsd_results_rpp %>% janitor::row_to_names(row_number = 1)
rsd_results_rpp <- sapply(rsd_results_rpp, as.numeric)
rsd_results_rpp <- as_tibble(rsd_results_rpp, .name_repair = "unique")

##------------------------------------------------------------------------------
# C18 neg
rpn_dat <- C18_Neg$data
rpn_dat$class <- C18_Neg$sample_meta$class
rpn_dat_melt <- melt(rpn_dat)

# Using dplyr and tidyr to calculate %RSD across groups for multiple columns
result_rpn <- rpn_dat %>%
  gather(key = "Variable", value = "Value", -class) %>%
  group_by(class, Variable) %>%
  summarise(RSD = rsd(Value)) %>%
  ungroup()

rsd_rpn_density <- ggplot(result_rpn, aes(x = RSD, fill = class, color =
                                            class), na.rm = TRUE) +
  geom_density(alpha = 0.1) + theme_bw(20) +
  labs(x = "RSD (%)", title = "C18 ESI-") +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  scale_color_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black")) +
  facet_grid(~ class)

rsd_rpn_box <- ggplot(result_rpn, aes(x = class, y = RSD, fill = class)) +
  geom_boxplot(alpha = 0.5) + theme_bw(20) +
  labs(y = "RSD (%)", title = "C18 ESI-") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(
      size = rel(0.8),
      angle = 45,
      hjust = 1,
      vjust = 1.1
    )
  ) +
  scale_fill_manual(values = c("#386cb0", "#ef3b2c", "#7fc97f", "black"))

rsd_results_rpn <- dcast(result_rpn, class ~ Variable, value.var = "RSD")
rsd_results_rpn_final <- dcast(result_rpn, class ~ Variable, value.var = "RSD")
rsd_results_rpn <- as.data.frame(t(rsd_results_rpn))
rsd_results_rpn <- rsd_results_rpn %>% janitor::row_to_names(row_number = 1)
rsd_results_rpn <- sapply(rsd_results_rpn, as.numeric)
rsd_results_rpn <- as_tibble(rsd_results_rpn, .name_repair = "unique")

##------------------------------------------------------------------------------
rsd_results_hp %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")
rsd_results_hn %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")
rsd_results_rpp %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")
rsd_results_rpn %>%
  get_summary_stats(Chlor_Meth, Meth_Water, Meth_ACN, QC, type = "robust")

(rsd_hp_density / rsd_hn_density / rsd_rpp_density / rsd_rpn_density) +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  '/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/rsd_density.pdf',
  width = 10,
  height = 15
)

((rsd_hp_box | rsd_hn_box | rsd_rpp_box | rsd_rpn_box)  /
    (
      tic_Hilic_Pos_pre  |
        tic_Hilic_Neg_pre | tic_C18_Pos_pre | tic_C18_Neg_pre
    ) /
    (
      post_tic_Hilic_Pos  |
        post_tic_Hilic_Neg | post_tic_C18_Pos | post_tic_C18_Neg
    )
) +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 10))

ggsave(
  '/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/tic_processed.pdf',
  width = 12,
  height = 10.5
)


# mapping the data frames onto the list
data_frames <- list(
  "RSD HILIC ESI+" = rsd_results_hp_final,
  "RSD HILIC ESI-" = rsd_results_hn_final,
  "RSD C18 ESI+" = rsd_results_rpp_final,
  "RSD C18 ESI-" = rsd_results_rpn_final
)

# writing the list of data frames onto the xlsx file
write.xlsx(
  data_frames,
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/Supplementary Table S2.xlsx"
)




(rsd_hp_box | rsd_hn_box) / (rsd_rpp_box | rsd_rpn_box) +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.text = element_text(size = 10))
