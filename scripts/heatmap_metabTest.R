# Load and prepare Hilic_Pos heatmap
set.seed(10)
Hilic_Pos_tab <- Hilic_Pos$data[1:15, ]
Hilic_Pos_hm <- as.data.frame(t(Hilic_Pos_tab))
rownames(Hilic_Pos_hm) <- Hilic_Pos$variable_meta$Compound

annotation_col <- data.frame(Extraction = factor(rep(c("Chlor_Meth", "Meth_Water", "Meth_ACN"), each = 5)))
rownames(annotation_col) <- rownames(Hilic_Pos_tab)

ann_colors <- list(Extraction = c(Chlor_Meth = "#386cb0", Meth_Water = "#7fc97f", Meth_ACN = "#ef3b2c"))

Hilic_Pos_pheat <- ggplotify::as.ggplot(pheatmap(
  Hilic_Pos_hm,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  border_color = NA,
  show_rownames = FALSE,
  cutree_cols = 2, cutree_rows = 2,
  legend = FALSE, annotation_legend = FALSE,
  main = "Hilic ESI+"
))

# Load and prepare Hilic_Neg heatmap
Hilic_Neg_tab <- Hilic_Neg$data[1:15, ]
Hilic_Neg_hm <- as.data.frame(t(Hilic_Neg_tab))
rownames(Hilic_Neg_hm) <- Hilic_Neg$variable_meta$Compound
rownames(annotation_col) <- rownames(Hilic_Neg_tab)

Hilic_Neg_pheat <- ggplotify::as.ggplot(pheatmap(
  Hilic_Neg_hm,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  border_color = NA,
  show_rownames = FALSE,
  cutree_cols = 2, cutree_rows = 2,
  legend = FALSE, annotation_legend = FALSE
))

# Load and prepare C18_Pos heatmap
C18_Pos_tab <- C18_Pos$data[1:15, ]
C18_Pos_hm <- as.data.frame(t(C18_Pos_tab))
rownames(C18_Pos_hm) <- C18_Pos$variable_meta$Compound
rownames(annotation_col) <- rownames(C18_Pos_tab)

C18_Pos_pheat <- ggplotify::as.ggplot(pheatmap(
  C18_Pos_hm,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  border_color = NA,
  show_rownames = FALSE,
  cutree_cols = 2, cutree_rows = 2,
  legend = FALSE, annotation_legend = FALSE,
  main = "C18 ESI+"
))

# Load and prepare C18_Neg heatmap
C18_Neg_tab <- C18_Neg$data[1:15, ]
C18_Neg_hm <- as.data.frame(t(C18_Neg_tab))
rownames(C18_Neg_hm) <- C18_Neg$variable_meta$Compound
rownames(annotation_col) <- rownames(C18_Neg_tab)

C18_Neg_pheat <- ggplotify::as.ggplot(pheatmap(
  C18_Neg_hm,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "row",
  border_color = NA,
  show_rownames = FALSE,
  cutree_cols = 2, cutree_rows = 2,
  annotation_legend = TRUE,
  main = "C18 ESI-"
))

# Combine all heatmaps into one figure
ggarrange(
  Hilic_Pos_pheat, C18_Pos_pheat, C18_Neg_pheat,
  common.legend = TRUE, legend = "bottom",
  nrow = 1, ncol = 3, widths = c(1, 1, 1.6)
)

# Save final figure
ggsave('/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/figures/heat.pdf', width = 10, height = 7)
