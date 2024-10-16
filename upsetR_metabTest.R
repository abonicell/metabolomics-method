## ---------------------------------------------------------------------------
# load Hilic ESI+
TT = filter_smeta(mode = 'include',factor_name = 'group',levels = c('A','QC'))  

# apply model
TT = model_apply(TT, DE_Hilic_Pos)
chlor_meth = predicted(TT)
chlor_meth

blk_chlor_meth <- model_apply(blk_filter, chlor_meth) %>%
  predicted()

blk_perc_chlor_meth <- model_apply(perc_features, blk_chlor_meth) %>%
  predicted()

chlor_meth <- model_apply(qc_features, blk_perc_chlor_meth) %>%
  predicted()

chlor_meth$sample_meta$batch = factor(chlor_meth$sample_meta$batch)
chlor_meth$sample_meta$type = factor(chlor_meth$sample_meta$type)
chlor_meth$sample_meta$class = factor(chlor_meth$sample_meta$class)
chlor_meth$sample_meta$extraction_type = factor(chlor_meth$sample_meta$extraction_type)

# meth_water
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('B','QC'))  

# apply model
TT = model_apply(TT, DE_Hilic_Pos)
meth_water = predicted(TT)
meth_water

blk_meth_water <- model_apply(blk_filter, meth_water) %>%
  predicted()

blk_perc_meth_water <- model_apply(perc_features, blk_meth_water) %>%
  predicted()

meth_water <- model_apply(qc_features, blk_perc_meth_water) %>%
  predicted()

meth_water$sample_meta$batch = factor(meth_water$sample_meta$batch)
meth_water$sample_meta$type = factor(meth_water$sample_meta$type)
meth_water$sample_meta$class = factor(meth_water$sample_meta$class)
meth_water$sample_meta$extraction_type = factor(meth_water$sample_meta$extraction_type)

# meth_ACN
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('C','QC'))  

# apply model
TT = model_apply(TT, DE_Hilic_Pos)
meth_ACN = predicted(TT)
meth_ACN

blk_meth_ACN <- model_apply(blk_filter, meth_ACN) %>%
  predicted()

blk_perc_meth_ACN <- model_apply(perc_features, blk_meth_ACN) %>%
  predicted()

meth_ACN <- model_apply(qc_features, blk_perc_meth_ACN) %>%
  predicted()

# drift correction
# convert to factors
meth_ACN$sample_meta$batch = factor(meth_ACN$sample_meta$batch)
meth_ACN$sample_meta$type = factor(meth_ACN$sample_meta$type)
meth_ACN$sample_meta$class = factor(meth_ACN$sample_meta$class)
meth_ACN$sample_meta$extraction_type = factor(meth_ACN$sample_meta$extraction_type)

# isolate samples only
TT = filter_smeta(mode = 'include', factor_name = 'type', levels = 'Sample')  

# apply model
TT = model_apply(TT, chlor_meth)
chlor_meth = predicted(TT)
chlor_meth

# apply model
TT = model_apply(TT, meth_water)
meth_water = predicted(TT)
meth_water

# apply model
TT = model_apply(TT, meth_ACN)
meth_ACN = predicted(TT)
meth_ACN

chlor_meth_upsetr <- chlor_meth[[2]]
meth_water_upsetr <- meth_water[[2]]
meth_ACN_upsetr <- meth_ACN[[2]]

lt_rpn = list(Chlor_Meth= c(chlor_meth_upsetr),
              Meth_Water = c(meth_water_upsetr),
              Meth_ACN = c(meth_ACN_upsetr))

m = make_comb_mat(lt_rpn)
m

ss = set_size(m)
cs = comb_size(m)

hp = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Set Size" = anno_barplot(ss, 
                                       border = FALSE, 
                                       gp = gpar(fill = c(  Chlor_Meth = "#386cb0" ,
                                                            Meth_Water = "#7fc97f",
                                                            Meth_ACN = "#ef3b2c")), 
                                       width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE) +
  grid.text("HILIC ESI+", x = 0.15, y = 0.95,
            gp = gpar(fontsize = 20))


hp = draw(hp) +
  grid.text("HILIC ESI+", x = 0.13, y = 0.95,
            gp = gpar(fontsize = 20))
od = column_order(hp)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 12, col = "#404040"), rot = 45)
})

vp1 = grid.grab()

## ---------------------------------------------------------------------------
# load Hilic ESI-
TT = filter_smeta(mode = 'include',factor_name = 'group',levels = c('A','QC'))  

# apply model
TT = model_apply(TT, DE_Hilic_Neg)
chlor_meth = predicted(TT)
chlor_meth

blk_chlor_meth <- model_apply(blk_filter, chlor_meth) %>%
  predicted()

blk_perc_chlor_meth <- model_apply(perc_features, blk_chlor_meth) %>%
  predicted()

chlor_meth <- model_apply(qc_features, blk_perc_chlor_meth) %>%
  predicted()

chlor_meth$sample_meta$batch = factor(chlor_meth$sample_meta$batch)
chlor_meth$sample_meta$type = factor(chlor_meth$sample_meta$type)
chlor_meth$sample_meta$class = factor(chlor_meth$sample_meta$class)
chlor_meth$sample_meta$extraction_type = factor(chlor_meth$sample_meta$extraction_type)

# meth_water
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('B','QC'))  

# apply model
TT = model_apply(TT, DE_Hilic_Neg)
meth_water = predicted(TT)
meth_water

blk_meth_water <- model_apply(blk_filter, meth_water) %>%
  predicted()

blk_perc_meth_water <- model_apply(perc_features, blk_meth_water) %>%
  predicted()

meth_water <- model_apply(qc_features, blk_perc_meth_water) %>%
  predicted()

meth_water$sample_meta$batch = factor(meth_water$sample_meta$batch)
meth_water$sample_meta$type = factor(meth_water$sample_meta$type)
meth_water$sample_meta$class = factor(meth_water$sample_meta$class)
meth_water$sample_meta$extraction_type = factor(meth_water$sample_meta$extraction_type)

# meth_ACN
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('C','QC'))  

# apply model
TT = model_apply(TT, DE_Hilic_Neg)
meth_ACN = predicted(TT)
meth_ACN

blk_meth_ACN <- model_apply(blk_filter, meth_ACN) %>%
  predicted()

blk_perc_meth_ACN <- model_apply(perc_features, blk_meth_ACN) %>%
  predicted()

meth_ACN <- model_apply(qc_features, blk_perc_meth_ACN) %>%
  predicted()

# drift correction
# convert to factors
meth_ACN$sample_meta$batch = factor(meth_ACN$sample_meta$batch)
meth_ACN$sample_meta$type = factor(meth_ACN$sample_meta$type)
meth_ACN$sample_meta$class = factor(meth_ACN$sample_meta$class)
meth_ACN$sample_meta$extraction_type = factor(meth_ACN$sample_meta$extraction_type)

# isolate samples only
TT = filter_smeta(mode = 'include', factor_name = 'type', levels = 'Sample')  

# apply model
TT = model_apply(TT, chlor_meth)
chlor_meth = predicted(TT)
chlor_meth

# apply model
TT = model_apply(TT, meth_water)
meth_water = predicted(TT)
meth_water

# apply model
TT = model_apply(TT, meth_ACN)
meth_ACN = predicted(TT)
meth_ACN

chlor_meth_upsetr <- chlor_meth[[1]]
meth_water_upsetr <- meth_water[[1]]
meth_ACN_upsetr <- meth_ACN[[1]]

lt_rpn = list(Chlor_Meth= c(chlor_meth_upsetr),
              Meth_Water = c(meth_water_upsetr),
              Meth_ACN = c(meth_ACN_upsetr))

m = make_comb_mat(lt_rpn)
m

ss = set_size(m)
cs = comb_size(m)

hn = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Set Size" = anno_barplot(ss, 
                                       border = FALSE, 
                                       gp = gpar(fill = c(  Chlor_Meth = "#386cb0" ,
                                                            Meth_Water = "#7fc97f",
                                                            Meth_ACN = "#ef3b2c")), 
                                       width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)


hn = draw(hn)+
  grid.text("HILIC ESI-", x = 0.13, y = 0.95,
            gp = gpar(fontsize = 20))
od = column_order(hn)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 12, col = "#404040"), rot = 45)
})

vp2 = grid.grab()

## ---------------------------------------------------------------------------
# load C18 ESI+
TT = filter_smeta(mode = 'include',factor_name = 'group',levels = c('A','QC'))  

# apply model
TT = model_apply(TT, DE_C18_Pos)
chlor_meth = predicted(TT)
chlor_meth

blk_chlor_meth <- model_apply(blk_filter, chlor_meth) %>%
  predicted()

blk_perc_chlor_meth <- model_apply(perc_features, blk_chlor_meth) %>%
  predicted()

chlor_meth <- model_apply(qc_features, blk_perc_chlor_meth) %>%
  predicted()

chlor_meth$sample_meta$batch = factor(chlor_meth$sample_meta$batch)
chlor_meth$sample_meta$type = factor(chlor_meth$sample_meta$type)
chlor_meth$sample_meta$class = factor(chlor_meth$sample_meta$class)
chlor_meth$sample_meta$extraction_type = factor(chlor_meth$sample_meta$extraction_type)

# meth_water
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('B','QC'))  

# apply model
TT = model_apply(TT, DE_C18_Pos)
meth_water = predicted(TT)
meth_water

blk_meth_water <- model_apply(blk_filter, meth_water) %>%
  predicted()

blk_perc_meth_water <- model_apply(perc_features, blk_meth_water) %>%
  predicted()

meth_water <- model_apply(qc_features, blk_perc_meth_water) %>%
  predicted()

meth_water$sample_meta$batch = factor(meth_water$sample_meta$batch)
meth_water$sample_meta$type = factor(meth_water$sample_meta$type)
meth_water$sample_meta$class = factor(meth_water$sample_meta$class)
meth_water$sample_meta$extraction_type = factor(meth_water$sample_meta$extraction_type)

# meth_ACN
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('C','QC'))  

# apply model
TT = model_apply(TT, DE_C18_Pos)
meth_ACN = predicted(TT)
meth_ACN

blk_meth_ACN <- model_apply(blk_filter, meth_ACN) %>%
  predicted()

blk_perc_meth_ACN <- model_apply(perc_features, blk_meth_ACN) %>%
  predicted()

meth_ACN <- model_apply(qc_features, blk_perc_meth_ACN) %>%
  predicted()

# drift correction
# convert to factors
meth_ACN$sample_meta$batch = factor(meth_ACN$sample_meta$batch)
meth_ACN$sample_meta$type = factor(meth_ACN$sample_meta$type)
meth_ACN$sample_meta$class = factor(meth_ACN$sample_meta$class)
meth_ACN$sample_meta$extraction_type = factor(meth_ACN$sample_meta$extraction_type)

# isolate samples only
TT = filter_smeta(mode = 'include', factor_name = 'type', levels = 'Sample')  

# apply model
TT = model_apply(TT, chlor_meth)
chlor_meth = predicted(TT)
chlor_meth

# apply model
TT = model_apply(TT, meth_water)
meth_water = predicted(TT)
meth_water

# apply model
TT = model_apply(TT, meth_ACN)
meth_ACN = predicted(TT)
meth_ACN

chlor_meth_upsetr <- chlor_meth[[1]]
meth_water_upsetr <- meth_water[[1]]
meth_ACN_upsetr <- meth_ACN[[1]]

lt_rpn = list(Chlor_Meth= c(chlor_meth_upsetr),
              Meth_Water = c(meth_water_upsetr),
              Meth_ACN = c(meth_ACN_upsetr))

m = make_comb_mat(lt_rpn)
m

ss = set_size(m)
cs = comb_size(m)

rpp = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Set Size" = anno_barplot(ss, 
                                       border = FALSE, 
                                       gp = gpar(fill = c(  Chlor_Meth = "#386cb0" ,
                                                            Meth_Water = "#7fc97f",
                                                            Meth_ACN = "#ef3b2c")), 
                                       width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)


rpp = draw(rpp) +
  grid.text("C18 ESI+", x = 0.13, y = 0.95,
            gp = gpar(fontsize = 20))
od = column_order(rpp)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 12, col = "#404040"), rot = 45)
})

vp3 = grid.grab()
## ---------------------------------------------------------------------------
# load C18 ESI-
TT = filter_smeta(mode = 'include',factor_name = 'group',levels = c('A','QC'))  

# apply model
TT = model_apply(TT, DE_C18_Neg)
chlor_meth = predicted(TT)
chlor_meth

blk_chlor_meth <- model_apply(blk_filter, chlor_meth) %>%
  predicted()

blk_perc_chlor_meth <- model_apply(perc_features, blk_chlor_meth) %>%
  predicted()

chlor_meth <- model_apply(qc_features, blk_perc_chlor_meth) %>%
  predicted()

chlor_meth$sample_meta$batch = factor(chlor_meth$sample_meta$batch)
chlor_meth$sample_meta$type = factor(chlor_meth$sample_meta$type)
chlor_meth$sample_meta$class = factor(chlor_meth$sample_meta$class)
chlor_meth$sample_meta$extraction_type = factor(chlor_meth$sample_meta$extraction_type)

# meth_water
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('B','QC'))  

# apply model
TT = model_apply(TT, DE_C18_Neg)
meth_water = predicted(TT)
meth_water

blk_meth_water <- model_apply(blk_filter, meth_water) %>%
  predicted()

blk_perc_meth_water <- model_apply(perc_features, blk_meth_water) %>%
  predicted()

meth_water <- model_apply(qc_features, blk_perc_meth_water) %>%
  predicted()

meth_water$sample_meta$batch = factor(meth_water$sample_meta$batch)
meth_water$sample_meta$type = factor(meth_water$sample_meta$type)
meth_water$sample_meta$class = factor(meth_water$sample_meta$class)
meth_water$sample_meta$extraction_type = factor(meth_water$sample_meta$extraction_type)

# meth_ACN
TT = filter_smeta(mode = 'include', factor_name = 'group', levels = c('C','QC'))  

# apply model
TT = model_apply(TT, DE_C18_Neg)
meth_ACN = predicted(TT)
meth_ACN

blk_meth_ACN <- model_apply(blk_filter, meth_ACN) %>%
  predicted()

blk_perc_meth_ACN <- model_apply(perc_features, blk_meth_ACN) %>%
  predicted()

meth_ACN <- model_apply(qc_features, blk_perc_meth_ACN) %>%
  predicted()

# drift correction
# convert to factors
meth_ACN$sample_meta$batch = factor(meth_ACN$sample_meta$batch)
meth_ACN$sample_meta$type = factor(meth_ACN$sample_meta$type)
meth_ACN$sample_meta$class = factor(meth_ACN$sample_meta$class)
meth_ACN$sample_meta$extraction_type = factor(meth_ACN$sample_meta$extraction_type)

# isolate samples only
TT = filter_smeta(mode = 'include', factor_name = 'type', levels = 'Sample')  

# apply model
TT = model_apply(TT, chlor_meth)
chlor_meth = predicted(TT)
chlor_meth

# apply model
TT = model_apply(TT, meth_water)
meth_water = predicted(TT)
meth_water

# apply model
TT = model_apply(TT, meth_ACN)
meth_ACN = predicted(TT)
meth_ACN

chlor_meth_upsetr <- chlor_meth[[1]]
meth_water_upsetr <- meth_water[[1]]
meth_ACN_upsetr <- meth_ACN[[1]]

lt_rpn = list(Chlor_Meth= c(chlor_meth_upsetr),
              Meth_Water = c(meth_water_upsetr),
              Meth_ACN = c(meth_ACN_upsetr))

m = make_comb_mat(lt_rpn)
m

ss = set_size(m)
cs = comb_size(m)

rpn = UpSet(m, 
           set_order = order(ss),
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Set Size" = anno_barplot(ss, 
                                       border = FALSE, 
                                       gp = gpar(fill = c(  Chlor_Meth = "#386cb0" ,
                                                            Meth_Water = "#7fc97f",
                                                            Meth_ACN = "#ef3b2c")), 
                                       width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)


rpn = draw(rpn) +
  grid.text("C18 ESI-", x = 0.13, y = 0.95,
            gp = gpar(fontsize = 20))
od = column_order(rpn)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("left", "bottom"), 
            gp = gpar(fontsize = 12, col = "#404040"), rot = 45)
})


## ---------------------------------------------------------------
hp_dat <- Hilic_Pos$data
hn_dat <- Hilic_Neg$data
C18_Pos_dat <- C18_Pos$data
C18_Neg_dat <- C18_Neg$data

## full sample
y <- list(
  "HILIC ESI+" = colnames(hp_dat), 
  "HILIC ESI-" = colnames(hn_dat), 
  "C18 ESI+" = colnames(C18_Pos_dat),
  "C18 ESI-" = colnames(C18_Neg_dat))

library(ggvenn)
venn <- ggvenn(
  y, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 5, fill_alpha = 0.2
)

ggsave('venn.pdf', width = 7, height = 7)

v.table <- gplots::venn(y, show.plot=FALSE)
intersections <- attr(v.table,"intersections")
print(intersections$`Hilic ESI-:C18 ESI-`)

vp4 = grid.grab()
grid.arrange(vp1,vp2,vp3,vp4)

