library(reshape2)
library(dplyr)
library(rstatix)    # For anova_test, tukey_hsd, add_significance
library(openxlsx)   # For write.xlsx

# Helper function to perform ANOVA and Tukey HSD tests and export results
run_univariate_stats <- function(filtered_data, output_prefix) {
  data <- filtered_data$data
  data$class <- filtered_data$sample_meta$class
  data_melted <- melt(data, id.vars = "class", variable.name = "variable", value.name = "value")
  
  # ANOVA
  anova_res <- data_melted %>%
    group_by(variable) %>%
    anova_test(value ~ class) %>%
    add_significance()
  
  write.xlsx(anova_res, paste0(output_prefix, "_anova.xlsx"))
  
  # Tukey HSD post-hoc test with FDR correction
  tukey_res <- data_melted %>%
    group_by(variable) %>%
    tukey_hsd(value ~ class, p.adjust.method = "fdr") %>%
    add_significance()
  
  write.xlsx(tukey_res, paste0(output_prefix, "_tukey.xlsx"))
}

# Apply the function to each dataset
run_univariate_stats(Hilic_Pos_filtered, "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_Hilic_Pos")
run_univariate_stats(Hilic_Neg_filtered, "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_Hilic_Neg")
run_univariate_stats(C18_Pos_filtered, "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_C18_Pos")
run_univariate_stats(C18_Neg_filtered, "/Users/andreabonicelli/Documents/GitHub/metabolomics-method/scripts/tables/anova_C18_Neg")
