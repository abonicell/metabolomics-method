suppressPackageStartupMessages({
  # Bioconductor packages
  library(structToolbox)
  library(pmp)
  library(ropls)
  library(BiocFileCache)
  library(SDAMS)
  
  # CRAN libraries
  library(ComplexHeatmap)
  library(ComplexUpset)
  library(ggpubr)
  library(gridExtra)
  library(patchwork)
  library(openxlsx)
  library(struct)
  library(pheatmap)
  library(rstatix)
  library(reshape2)
  library(UpSetR)
  library(rtemis)
})

## ---------------------------------------------------------------
# Load data
data_hp <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/data_hp.csv",
  row.names = 1,
  check.names = FALSE
)
variable_meta_hp <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/variable_meta_hp.csv",
  row.names = 1,
  check.names = FALSE
)

data_hn <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/data_hn.csv",
  row.names = 1,
  check.names = FALSE
)
variable_meta_hn <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/variable_meta_hn.csv",
  row.names = 1,
  check.names = FALSE
)

data_rpp <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/data_rpp.csv",
  row.names = 1,
  check.names = FALSE
)
variable_meta_rpp <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/variable_meta_rpp.csv",
  row.names = 1,
  check.names = FALSE
)

data_rpn <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/data_rpn.csv",
  row.names = 1,
  check.names = FALSE
)
variable_meta_rpn <- read.csv(
  "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/data/variable_meta_rpn.csv",
  row.names = 1,
  check.names = FALSE
)

# metadata are the same for all datasets
sample_meta <- read.csv("sample_meta.csv")

## ---------------------------------------------------------------
# Create DatasetExperiment objects for each dataset and mode
DE_Hilic_Pos <- DatasetExperiment(
  data = t(data_hp),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_hp, row.names = colnames(t(data_hp))),
  description = 'Metabolomics Testing',
  name = "Hilic ESI+"
)

DE_Hilic_Neg <- DatasetExperiment(
  data = t(data_hn),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_hn, row.names = colnames(t(data_hn))),
  description = 'Metabolomics Testing',
  name = "Hilic ESI-"
)

DE_C18_Pos <- DatasetExperiment(
  data = t(data_rpp),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_rpp, row.names = colnames(t(data_rpp))),
  description = 'Metabolomics Testing',
  name = "C18 ESI+"
)

DE_C18_Neg <- DatasetExperiment(
  data = t(data_rpn),
  sample_meta = sample_meta,
  variable_meta = as.data.frame(variable_meta_rpn, row.names = colnames(t(data_rpn))),
  description = 'Metabolomics Testing',
  name = "C18 ESI-"
)

## ---------------------------------------------------------------
# Convert sample metadata columns to factors for all DatasetExperiments
for (DE in list(DE_Hilic_Pos, DE_Hilic_Neg, DE_C18_Pos, DE_C18_Neg)) {
  DE$sample_meta$batch <- factor(DE$sample_meta$batch)
  DE$sample_meta$type <- factor(DE$sample_meta$type)
  DE$sample_meta$class <- factor(DE$sample_meta$class)
  DE$sample_meta$extraction_type <- factor(DE$sample_meta$extraction_type)
}

## ---------------------------------------------------------------
# Matrix processing - set up cleaning steps

# Blank filter
blk_filter <- blank_filter(
  fold_change = 15,
  blank_label = "Blank",
  qc_label = "QC",
  factor_name = 'type'
)

# Missing value feature filter
perc_features <- mv_feature_filter(
  threshold = 60,
  method = 'across',
  factor_name = 'type'
)

# RSD QC by feature filter
qc_features <- rsd_filter(
  rsd_threshold = 30,
  qc_label = 'QC',
  factor_name = 'type'
)

## ---------------------------------------------------------------
# Apply blank filter
blk_Hilic_Pos <- model_apply(blk_filter, DE_Hilic_Pos) |> predicted()
blk_Hilic_Neg <- model_apply(blk_filter, DE_Hilic_Neg) |> predicted()
blk_C18_Pos <- model_apply(blk_filter, DE_C18_Pos) |> predicted()
blk_C18_Neg <- model_apply(blk_filter, DE_C18_Neg) |> predicted()

# Apply missing data filter
blk_perc_Hilic_Pos <- model_apply(perc_features, blk_Hilic_Pos) |> predicted()
blk_perc_Hilic_Neg <- model_apply(perc_features, blk_Hilic_Neg) |> predicted()
blk_perc_C18_Pos <- model_apply(perc_features, blk_C18_Pos) |> predicted()
blk_perc_C18_Neg <- model_apply(perc_features, blk_C18_Neg) |> predicted()

# Apply RSD QC filter
blk_perc_qc_Hilic_Pos <- model_apply(qc_features, blk_perc_Hilic_Pos)
Hilic_Pos <- predicted(blk_perc_qc_Hilic_Pos)

blk_perc_qc_Hilic_Neg <- model_apply(qc_features, blk_perc_Hilic_Neg)
Hilic_Neg <- predicted(blk_perc_qc_Hilic_Neg)

blk_perc_qc_C18_Pos <- model_apply(qc_features, blk_perc_C18_Pos)
C18_Pos <- predicted(blk_perc_qc_C18_Pos)

blk_perc_qc_C18_Neg <- model_apply(qc_features, blk_perc_C18_Neg)
C18_Neg <- predicted(blk_perc_qc_C18_Neg)

## ---------------------------------------------------------------
# Total ion count plot setup
C <- tic_chart(
  factor_name = 'class',
  run_order = 'run_order',
  connected = FALSE
)

tic_Hilic_Pos_pre <- chart_plot(C, Hilic_Pos) + theme_bw(18) + ggtitle("HILIC ESI+")
tic_Hilic_Neg_pre <- chart_plot(C, Hilic_Neg) + theme_bw(18) + ggtitle("HILIC ESI-")
tic_C18_Pos_pre <- chart_plot(C, C18_Pos) + theme_bw(18) + ggtitle("C18 ESI+")
tic_C18_Neg_pre <- chart_plot(C, C18_Neg) + theme_bw(18) + ggtitle("C18 ESI-")

## ---------------------------------------------------------------
# PCA and data imputation

# KNN imputation
process <- knn_impute(neighbours = 5, sample_max = 100)

Hilic_Pos <- model_apply(process, Hilic_Pos)
Hilic_Neg <- model_apply(process, Hilic_Neg)
C18_Pos <- model_apply(process, C18_Pos)
C18_Neg <- model_apply(process, C18_Neg)

# Get imputed data matrices
Hilic_Pos <- predicted(Hilic_Pos)
Hilic_Neg <- predicted(Hilic_Neg)
C18_Pos <- predicted(C18_Pos)
C18_Neg <- predicted(C18_Neg)

# PCA model
PCA <- mean_centre() + structToolbox::PCA(number_components = 5)

PCA_Hilic_Pos <- model_apply(PCA, Hilic_Pos)
PCA_Hilic_Neg <- model_apply(PCA, Hilic_Neg)
PCA_C18_Pos <- model_apply(PCA, C18_Pos)
PCA_C18_Neg <- model_apply(PCA, C18_Neg)

# PCA scores plot setup
score <- pca_scores_plot(
  factor_name = 'class',
  points_to_label = "outliers",
  ellipse = 'sample',
  label_size = 5
)

pca_Hilic_Pos <- chart_plot(score, PCA_Hilic_Pos[2]) + theme_bw(16) + ggtitle("PCA Hilic ESI+")
pca_Hilic_Neg <- chart_plot(score, PCA_Hilic_Neg[2]) + theme_bw(16) + ggtitle("PCA Hilic ESI-")
pca_C18_Pos <- chart_plot(score, PCA_C18_Pos[2]) + theme_bw(16) + ggtitle("PCA C18 ESI+")
pca_C18_Neg <- chart_plot(score, PCA_C18_Neg[2]) + theme_bw(16) + ggtitle("PCA C18 ESI-")

(pca_Hilic_Pos | pca_Hilic_Neg | pca_C18_Pos | pca_C18_Neg) +
  plot_annotation(tag_levels = c('A')) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

## ---------------------------------------------------------------
# Boxplots to check data distributions before and after normalization

box_sample <- DatasetExperiment_boxplot(
  factor_name = 'class',
  number = 25,
  by_sample = TRUE,
  per_class = TRUE
)

pre_box_sample_Hilic_Pos <- chart_plot(box_sample, Hilic_Pos) + theme_bw(18) + ggtitle("HILIC ESI+")
pre_box_sample_Hilic_Neg <- chart_plot(box_sample, Hilic_Neg) + theme_bw(18) + ggtitle("HILIC ESI-")
pre_box_sample_C18_Pos <- chart_plot(box_sample, C18_Pos) + theme_bw(18) + ggtitle("C18 ESI+")
pre_box_sample_C18_Neg <- chart_plot(box_sample, C18_Neg) + theme_bw(18) + ggtitle("C18 ESI-")

## ---------------------------------------------------------------
# Data normalization and transformation

process <- pqn_norm(qc_label = 'QC', factor_name = 'type') +
  glog_transform(qc_label = 'QC', factor_name = 'class')

Hilic_Pos <- model_apply(process, Hilic_Pos)
Hilic_Neg <- model_apply(process, Hilic_Neg)
C18_Pos <- model_apply(process, C18_Pos)
C18_Neg <- model_apply(process, C18_Neg)

# Get transformed data matrices
Hilic_Pos <- predicted(Hilic_Pos)
Hilic_Neg <- predicted(Hilic_Neg)
C18_Pos <- predicted(C18_Pos)
C18_Neg <- predicted(C18_Neg)

box_sample_Hilic_Pos <- chart_plot(box_sample, Hilic_Pos) + theme_bw(18) + ggtitle("HILIC ESI+")
box_sample_Hilic_Neg <- chart_plot(box_sample, Hilic_Neg) + theme_bw(18) + ggtitle("HILIC ESI-")
box_sample_C18_Pos <- chart_plot(box_sample, C18_Pos) + theme_bw(18) + ggtitle("C18 ESI+")
box_sample_C18_Neg <- chart_plot(box_sample, C18_Neg) + theme_bw(18) + ggtitle("C18 ESI-")

# Arrange boxplots in a grid layout
(
  (pre_box_sample_Hilic_Pos / pre_box_sample_Hilic_Neg / pre_box_sample_C18_Pos / pre_box_sample_C18_Neg) |
    (box_sample_Hilic_Pos / box_sample_Hilic_Neg / box_sample_C18_Pos / box_sample_C18_Neg)
) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom")

ggsave(
  '/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/box_distribution.pdf',
  width = 13,
  height = 19
)

## ---------------------------------------------------------------
# Total ion count plots after processing
C <- tic_chart(
  factor_name = 'class',
  run_order = 'run_order',
  connected = FALSE
)

tic_Hilic_Pos <- chart_plot(C, Hilic_Pos) + theme_bw(18) + ggtitle("HILIC ESI+")
tic_Hilic_Neg <- chart_plot(C, Hilic_Neg) + theme_bw(18) + ggtitle("HILIC ESI-")
tic_C18_Pos <- chart_plot(C, C18_Pos) + theme_bw(18) + ggtitle("C18 ESI+")
tic_C18_Neg <- chart_plot(C, C18_Neg) + theme_bw(18) + ggtitle("C18 ESI-")

# Arrange and save TIC plots
(
  (tic_Hilic_Pos / tic_Hilic_Neg) |
    (tic_C18_Pos / tic_C18_Neg)
) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom")

ggsave(
  '/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/tic_plots.pdf',
  width = 13,
  height = 10
)
