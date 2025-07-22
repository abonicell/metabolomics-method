library(dplyr)
library(ComplexUpset)
library(ComplexHeatmap)
library(grid)
library(gridExtra)

# Helper function to process one dataset based on mode and models
process_dataset <- function(DE_model, mode_levels, blk_filter, perc_features, qc_features) {
  TT <- filter_smeta(mode = 'include', factor_name = 'group', levels = mode_levels)
  TT <- model_apply(TT, DE_model)
  predicted_data <- predicted(TT)
  
  predicted_data <- model_apply(blk_filter, predicted_data) %>% predicted()
  predicted_data <- model_apply(perc_features, predicted_data) %>% predicted()
  predicted_data <- model_apply(qc_features, predicted_data) %>% predicted()
  
  sample_meta_cols <- c('batch', 'type', 'class', 'extraction_type')
  for (col in sample_meta_cols) {
    predicted_data$sample_meta[[col]] <- factor(predicted_data$sample_meta[[col]])
  }
  
  return(predicted_data)
}

# Drift correction for Sample type only
drift_correction <- function(predicted_data_list) {
  TT <- filter_smeta(mode = 'include', factor_name = 'type', levels = 'Sample')
  corrected_list <- lapply(predicted_data_list, function(data) {
    TT <- model_apply(TT, data)
    predicted(TT)
  })
  return(corrected_list)
}

# Create UpSet plot
make_upset_plot <- function(predicted_data, title) {
  lt_rpn <- list(
    Chlor_Meth = c(predicted_data[[1]]),
    Meth_Water = c(predicted_data[[2]]),
    Meth_ACN = c(predicted_data[[3]])
  )
  
  m <- make_comb_mat(lt_rpn)
  ss <- set_size(m)
  cs <- comb_size(m)
  
  p <- UpSet(m,
             set_order = order(ss),
             comb_order = order(comb_degree(m), -cs),
             top_annotation = HeatmapAnnotation(
               "Intersection Size" = anno_barplot(cs,
                                                  ylim = c(0, max(cs)*1.1),
                                                  border = FALSE,
                                                  gp = gpar(fill = "black"),
                                                  height = unit(4, "cm")),
               annotation_name_side = "left",
               annotation_name_rot = 90),
             left_annotation = rowAnnotation(
               "Set Size" = anno_barplot(ss,
                                         border = FALSE,
                                         gp = gpar(fill = c(Chlor_Meth = "#386cb0",
                                                            Meth_Water = "#7fc97f",
                                                            Meth_ACN = "#ef3b2c")),
                                         width = unit(4, "cm")),
               set_name = anno_text(set_name(m),
                                    location = 0.5,
                                    just = "center",
                                    width = max_text_width(set_name(m)) + unit(4, "mm"))
             ),
             right_annotation = NULL,
             show_row_names = FALSE
  )
  
  p <- draw(p) +
    grid.text(title, x = 0.13, y = 0.95,
              gp = gpar(fontsize = 20))
  
  od <- column_order(p)
  decorate_annotation("Intersection Size", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
              default.units = "native", just = c("left", "bottom"),
              gp = gpar(fontsize = 12, col = "#404040"), rot = 45)
  })
  
  return(p)
}

# === PROCESSING ALL DATASETS ===

# HILIC ESI+
chlor_meth_pos <- process_dataset(DE_Hilic_Pos, c('A','QC'), blk_filter, perc_features, qc_features)
meth_water_pos <- process_dataset(DE_Hilic_Pos, c('B','QC'), blk_filter, perc_features, qc_features)
meth_ACN_pos <- process_dataset(DE_Hilic_Pos, c('C','QC'), blk_filter, perc_features, qc_features)

corrected_pos <- drift_correction(list(chlor_meth_pos, meth_water_pos, meth_ACN_pos))
hp <- make_upset_plot(corrected_pos, "HILIC ESI+")

# HILIC ESI-
chlor_meth_neg <- process_dataset(DE_Hilic_Neg, c('A','QC'), blk_filter, perc_features, qc_features)
meth_water_neg <- process_dataset(DE_Hilic_Neg, c('B','QC'), blk_filter, perc_features, qc_features)
meth_ACN_neg <- process_dataset(DE_Hilic_Neg, c('C','QC'), blk_filter, perc_features, qc_features)

corrected_neg <- drift_correction(list(chlor_meth_neg, meth_water_neg, meth_ACN_neg))
hn <- make_upset_plot(corrected_neg, "HILIC ESI-")

# C18 ESI+
chlor_meth_C18_pos <- process_dataset(DE_C18_Pos, c('A','QC'), blk_filter, perc_features, qc_features)
meth_water_C18_pos <- process_dataset(DE_C18_Pos, c('B','QC'), blk_filter, perc_features, qc_features)
meth_ACN_C18_pos <- process_dataset(DE_C18_Pos, c('C','QC'), blk_filter, perc_features, qc_features)

corrected_C18_pos <- drift_correction(list(chlor_meth_C18_pos, meth_water_C18_pos, meth_ACN_C18_pos))
rpp <- make_upset_plot(corrected_C18_pos, "C18 ESI+")

# C18 ESI-
chlor_meth_C18_neg <- process_dataset(DE_C18_Neg, c('A','QC'), blk_filter, perc_features, qc_features)
meth_water_C18_neg <- process_dataset(DE_C18_Neg, c('B','QC'), blk_filter, perc_features, qc_features)
meth_ACN_C18_neg <- process_dataset(DE_C18_Neg, c('C','QC'), blk_filter, perc_features, qc_features)

corrected_C18_neg <- drift_correction(list(chlor_meth_C18_neg, meth_water_C18_neg, meth_ACN_C18_neg))
rpn <- make_upset_plot(corrected_C18_neg, "C18 ESI-")

# === COMBINE ALL PLOTS INTO ONE FIGURE ===

# Use grid graphics viewport layout for 2x2
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

# Helper to print a plot in a grid cell
print_in_viewport <- function(plot, row, col) {
  vp <- viewport(layout.pos.row = row, layout.pos.col = col)
  pushViewport(vp)
  grid.draw(plot)
  popViewport()
}

print_in_viewport(hp, 1, 1)
print_in_viewport(hn, 1, 2)
print_in_viewport(rpp, 2, 1)
print_in_viewport(rpn, 2, 2)
