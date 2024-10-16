## ---------------------------------------------------------------
# ANOVA
# summary of filtered data
Hilic_Pos_uni <- Hilic_Pos_filtered$data
Hilic_Pos_uni$class <- Hilic_Pos_filtered$sample_meta$class
Hilic_Pos_uni_melted <- melt(Hilic_Pos_uni)

Hilic_Pos_uni_melted %>%
  group_by(variable) %>%
  anova_test(value ~ class ) %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_Hilic_Pos.xlsx")

Hilic_Pos_uni_melted %>%
  group_by(variable) %>%
  tukey_hsd(value ~ class,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/tukey_Hilic_Pos.xlsx")

# summary of filtered data
Hilic_Neg_uni <- Hilic_Neg_filtered$data
Hilic_Neg_uni$class <- Hilic_Neg_filtered$sample_meta$class
Hilic_Neg_uni_melted <- melt(Hilic_Neg_uni)

Hilic_Neg_uni_melted %>%
  group_by(variable) %>%
  anova_test(value ~ class ) %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_Hilic_Neg.xlsx")

Hilic_Neg_uni_melted %>%
  group_by(variable) %>%
  tukey_hsd(value ~ class,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/tukey_Hilic_Neg.xlsx")


# summary of filtered data
C18_Pos_uni <- C18_Pos_filtered$data
C18_Pos_uni$class <- C18_Pos_filtered$sample_meta$class
C18_Pos_uni_melted <- melt(C18_Pos_uni)

C18_Pos_uni_melted %>%
  group_by(variable) %>%
  anova_test(value ~ class ) %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_C18_Pos.xlsx")

C18_Pos_uni_melted %>%
  group_by(variable) %>%
  tukey_hsd(value ~ class,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/tukey_C18_Pos.xlsx")

# summary of filtered data
C18_Neg_uni <- C18_Neg_filtered$data
C18_Neg_uni$class <- C18_Neg_filtered$sample_meta$class
C18_Neg_uni_melted <- melt(C18_Neg_uni)

C18_Neg_uni_melted %>%
  group_by(variable) %>%
  anova_test(value ~ class ) %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_C18_Neg.xlsx")

C18_Neg_uni_melted %>%
  group_by(variable) %>%
  tukey_hsd(value ~ class,
            p.adjust.method = "fdr") %>%
  add_significance() %>%
  write.xlsx("/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/tukey_C18_Neg.xlsx")








