#Get ready
library(jsonlite)
library(plotly)
library(tidyr)
library(dplyr)
library(data.table)
library(htmlwidgets)
library(RColorBrewer)
library(lme4)
library(emmeans)
library(ggtext)
setwd("/data/gpfs/projects/punim1869/users/yunhongh1/workspace/intergrain/NARO_3D/CT_cb/rinfo")
# json xyz - polylines
extract_polylines <- function(data, polylines_list) {
  if (!is.null(data$polyline)) {
    polylines_list[[length(polylines_list) + 1]] <- data$polyline
  }
  
  if (is.list(data)) {
    for (subdata in data) {
      if (is.list(subdata)) {
        polylines_list <- extract_polylines(subdata, polylines_list)
      }
    }
  }
  
  return(polylines_list)
}
# read in data files
files_and_conditions <- list(
  # Surplus
  list("U17P2_35DAS.rinfo", "Surplus_Vlamingh_U17P2_35DAS.txt"),
  list("U19P2_35DAS.rinfo", "Surplus_Vlamingh_U19P2_35DAS.txt"),
  list("U20P6_35DAS.rinfo", "Surplus_Vlamingh_U20P6_35DAS.txt"),
  
  list("U19P5_35DAS.rinfo", "Surplus_Buff_U19P5_35DAS.txt"),
  list("U21P3_35DAS.rinfo", "Surplus_Buff_U21P3_35DAS.txt"),
  list("U21P1_35DAS.rinfo", "Surplus_Buff_U21P1_35DAS.txt"),
  
  list("U17P1_35DAS.rinfo", "Surplus_Roe_U17P1_35DAS.txt"),
  list("U17P4_35DAS.rinfo", "Surplus_Roe_U17P4_35DAS.txt"),
  list("U20P4_35DAS.rinfo", "Surplus_Roe_U20P4_35DAS.txt"),
  
  list("U17P5_35DAS.rinfo", "Surplus_LaTrobe_U17P5_35DAS.txt"),
  list("U19P1_35DAS.rinfo", "Surplus_LaTrobe_U19P1_35DAS.txt"),
  list("U24P1_35DAS.rinfo", "Surplus_LaTrobe_U24P1_35DAS.txt"),
  # Deficient
  list("U24P6_35DAS.rinfo", "Deficient_Vlamingh_U24P6_35DAS.txt"),
  list("U21P2_35DAS.rinfo", "Deficient_Vlamingh_U21P2_35DAS.txt"),
  list("U24P4_35DAS.rinfo", "Deficient_Vlamingh_U24P4_35DAS.txt"),
  
  list("U17P3_35DAS.rinfo", "Deficient_Buff_U17P3_35DAS.txt"),
  list("U21P4_35DAS.rinfo", "Deficient_Buff_U21P4_35DAS.txt"),
  list("U24P5_35DAS.rinfo", "Deficient_Buff_U24P5_35DAS.txt"),
  
  list("U19P6_35DAS.rinfo", "Deficient_Roe_U19P6_35DAS.txt"),
  list("U24P2_35DAS.rinfo", "Deficient_Roe_U24P2_35DAS.txt"),
  list("U24P3_35DAS.rinfo", "Deficient_Roe_U24P3_35DAS.txt"),
  
  list("U17P6_35DAS.rinfo", "Deficient_LaTrobe_U17P6_35DAS.txt"),
  list("U20P3_35DAS.rinfo", "Deficient_LaTrobe_U20P3_35DAS.txt"),
  list("U20P5_35DAS.rinfo", "Deficient_LaTrobe_U20P5_35DAS.txt")
)
# calculate the length of each root with 22-unit intervals --> so in total 30 layers (deepest 644, made it to 660 for easy calculation (660/30=22...)
calculate_root_length <- function(polyline, voxel_resolution = 0.3) {
  # Convert polyline to a data frame
  root_data <- as.data.frame(polyline)
  colnames(root_data) <- c("x", "y", "z")
  
  # Calculate the euclidean distance between points
  root_data_shifted <- root_data[-1, ]  # root i+1
  root_data_original <- root_data[-nrow(root_data), ]  # root i
  segment_lengths <- sqrt(
    (root_data_shifted$x - root_data_original$x)^2 + 
      (root_data_shifted$y - root_data_original$y)^2 + 
      (root_data_shifted$z - root_data_original$z)^2
  ) * voxel_resolution / 10  # Convert voxel to cm
  
  root_data$segment_length <- c(NA, segment_lengths)
  
  max_x <- max(root_data$x, na.rm = TRUE)
  breaks <- seq(0, 660, by = 22)  # every 22
  
  # Add depth layer to the root data
  root_data$depth_layer <- cut(root_data$x, 
                               breaks = breaks, 
                               include.lowest = TRUE, 
                               right = FALSE)  # Ensure correct intervals for depth
  
  return(list(root_data = root_data, total_length = sum(segment_lengths, na.rm = TRUE)))
}

# output
final_data_table <- data.frame()

# Loop all & calculate root lengths, and group into by 22-unit depth layer
for (file_info in files_and_conditions) {
  
  json_file <- file_info[[1]]
  data_table_file <- file_info[[2]]
  
  json_data <- fromJSON(json_file, simplifyDataFrame = FALSE)
  data_table <- fread(data_table_file)
  
  polylines_list <- list()
  polylines_list <- extract_polylines(json_data, polylines_list)
  
  total_root_length <- 0  #start from length=0
  depth_layer_lengths <- data.frame()  
  
  polylines_df <- lapply(1:length(polylines_list), function(i) {
    polyline <- polylines_list[[i]]
    root_type <- data_table$Root_Type[i]  # Get the Root_Type as a scalar value
    color <- get_root_color(root_type)  # color based on root type
    
    if (!is.null(polyline) && ncol(polyline) == 3) {
      # Calculate root length for each polyline
      root_length_data <- calculate_root_length(polyline)
      
      total_root_length <<- total_root_length + root_length_data$total_length
      
      depth_layer_summary <- root_length_data$root_data %>%
        group_by(depth_layer) %>%
        summarise(total_length = sum(segment_length, na.rm = TRUE))
      
      # Combine the lengths by depth layers across all polylines
      depth_layer_lengths <<- bind_rows(depth_layer_lengths, depth_layer_summary)
      
      return(list(
        df = data.frame(x = polyline[, 1], y = polyline[, 2], z = polyline[, 3]),
        color = color
      ))
    } else {
      return(NULL)  
    }
  })
  
  polylines_df <- polylines_df[!sapply(polylines_df, is.null)]
  
  # Group the total length by depth layer
  final_depth_layer_summary <- depth_layer_lengths %>%
    group_by(depth_layer) %>%
    summarise(total_length = sum(total_length, na.rm = TRUE)) %>%
    mutate(cumulative_length = cumsum(total_length))  
  
  # Fill in missing depth layers (no root grow anymore) with 0 for consistency across all replicates (NA is difficult for plotting?i guess)
  all_depth_layers <- data.frame(depth_layer = levels(cut(seq(0, 660, by = 22), breaks = seq(0, 660, by = 22), include.lowest = TRUE, right = FALSE)))
  final_depth_layer_summary <- full_join(all_depth_layers, final_depth_layer_summary, by = "depth_layer") %>%
    replace_na(list(total_length = 0, cumulative_length = 0))
  
  # Done - output
  final_data_table <- bind_rows(final_data_table, 
                                final_depth_layer_summary %>%
                                  mutate(genotype = strsplit(data_table_file, "_")[[1]][2],
                                         condition = strsplit(data_table_file, "_")[[1]][1],
                                         replicate = strsplit(data_table_file, "_")[[1]][3]))
}

head(final_data_table)

#Sort the order... top to bottom depth; G1-G4
depth_levels <- c("[0,22)", "[22,44)", "[44,66)", "[66,88)", "[88,110)", "[110,132)", "[132,154)", "[154,176)", 
                  "[176,198)", "[198,220)", "[220,242)", "[242,264)", "[264,286)", "[286,308)", "[308,330)", 
                  "[330,352)", "[352,374)", "[374,396)", "[396,418)", "[418,440)", "[440,462)", "[462,484)", 
                  "[484,506)", "[506,528)", "[528,550)", "[550,572)", "[572,594)", "[594,616)", "[616,638)", "[638,660]")
final_data_table$depth_layer <- factor(final_data_table$depth_layer, levels = depth_levels, ordered = TRUE)
final_data_table <- final_data_table[order(final_data_table$depth_layer), ]
final_data_table$genotype <- factor(final_data_table$genotype, levels = c("Vlamingh", "Buff", "Roe", "LaTrobe"))

# extract the numeric depth_layer
final_data_table <- final_data_table %>%
  mutate(
    depth_numeric = (as.numeric(factor(depth_layer)) - 1) * 0.66 + 0.66,  # Ensure depth_numeric is calculated; otherwise 
    cumulative_length = ifelse(cumulative_length == 0, NA, cumulative_length)  # OK, NA is good...f, remeber to go back fix it in the future
  )

# Plot
ggplot(final_data_table, aes(x = depth_numeric, y = cumulative_length, color = genotype, group = interaction(genotype, replicate))) +
  geom_line(alpha = 0.7) +  # Keep line transparency
  labs(
    title = "Cumulative Root Length by Depth for Each Replicate",
    x = "Soil Depth (cm)",
    y = "Cumulative Root Length (mm)"
  ) +
  theme_minimal() +
  facet_wrap(~ condition) +
  scale_color_manual(values = c("Vlamingh" = "#F9766D", "Buff" = "#7DAF00", "Roe" = "#00C6CC", "LaTrobe" = "#BD77F2")) +
  scale_x_continuous(breaks = seq(0.66, max(final_data_table$depth_numeric, na.rm = TRUE), by = 0.66))

#Stats
##lmm
lmm_interaction <- lmer(total_length ~ genotype * condition * depth_layer + (1 | replicate), data = final_data_table)
final_data_table$predict_total_length <- predict(lmm_interaction)
wide_data <- final_data_table %>%
  pivot_wider(names_from = depth_layer, values_from = c(total_length, cumulative_length, predict_total_length))
write.table(wide_data, "length_per_layer.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
##Multiple comparison
est_by_factor <- emmeans(lmm_interaction, pairwise ~ genotype * condition * depth_layer)
contrasts_diff_genotype_same_condition_same_layer <- contrast(est_by_factor, method = "pairwise", by = c("condition", "depth_layer"))
contrasts_diff_genotype_same_condition_same_layer_results <- summary(contrasts_diff_genotype_same_condition_same_layer)
contrasts_same_genotype_diff_condition_same_layer <- contrast(est_by_factor, method = "pairwise", by = c("genotype", "depth_layer"))
contrasts_same_genotype_diff_condition_same_layer_results <- summary(contrasts_same_genotype_diff_condition_same_layer)
##LASSO
data <- fread("/data/gpfs/projects/punim1869/users/yunhongh1/workspace/intergrain/NARO_3D/D2.txt") #deficient
data <- fread("/data/gpfs/projects/punim1869/users/yunhongh1/workspace/intergrain/NARO_3D/S2.txt") #surplus
data <- as.data.frame(data)
above_ground_traits <- c("Shootdryweight", "NofT", "CN", "N_content") #maybe add C% or N%? 
X <- as.matrix(data[, above_ground_traits])
root_traits <- grep("total_length_", names(data), value = TRUE)
for (root_trait in root_traits) {
  Y <- data[[root_trait]]  # Now Y is a root trait
  # LOOCV 
  lasso_model <- cv.glmnet(X, Y, alpha = 1, nfolds = nrow(X))
  best_lambda <- lasso_model$lambda.min
  
  coefficients <- coef(lasso_model, s = best_lambda)
  important_vars <- data.frame(Variable = rownames(coefficients), Coefficient = coefficients[, 1])
  
  important_vars <- important_vars[important_vars$Coefficient != 0, ]
  important_vars <- important_vars[important_vars$Variable != "(Intercept)", ]
  
  selected_vars <- colnames(X)[colnames(X) %in% important_vars$Variable]
  
  if (length(selected_vars) > 1) {  # If we have more than 1 important variable
    
    # Pearson
    correlation_results <- rcorr(as.matrix(data[, selected_vars]), Y)
    corr_matrix <- correlation_results$r
    p_matrix <- correlation_results$P
    
    diag(p_matrix) <- NA
    
    p_matrix_significance <- p_matrix
    p_matrix_significance[p_matrix < 0.01] <- '**'
    p_matrix_significance[p_matrix >= 0.01 & p_matrix < 0.05] <- '*'
    p_matrix_significance[p_matrix >= 0.05] <- ''
    
    print(p_matrix_significance)
    
    for (var in selected_vars) {
      df_plot <- data.frame(AboveGroundTrait = data[[var]], RootTrait = Y)
      
      significance_label <- as.character(p_matrix_significance[var, "y"]) 
      
      lm_fit <- lm(RootTrait ~ AboveGroundTrait, data = df_plot)
      intercept <- coef(lm_fit)[1]
      slope <- coef(lm_fit)[2]
      equation <- paste0("y = ", round(intercept, 2), " + ", round(slope, 2), "*x")
      
      equation_with_significance <- paste0(
        ifelse(significance_label != "", 
               paste0("<span style='color:red;'>", significance_label, "</span> "), ""),
        equation)
      
      p <- ggplot(df_plot, aes(x = AboveGroundTrait, y = RootTrait)) +
        geom_point(size = 3, color = "blue") +  
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +  
        labs(title = paste("Regression Line for", root_trait, "and", var, "with significance"),
             x = var, y = root_trait) +
        geom_richtext(aes(x = Inf, y = min(df_plot$RootTrait), 
                          label = equation_with_significance), hjust = 1.2, vjust = -1.2, size = 6, fill = NA, label.color = NA) +  
        theme_minimal()
      
      print(p)
      
      # Save
      filename <- paste0("regression_plot_", var, "_vs_", root_trait, ".png")
      ggsave(filename = filename, plot = p, width = 8, height = 6)
    }
    
  } else {
    message(paste("No important variable for root trait:", root_trait))
  }
}
##plot all
for (root_trait in root_traits) {
  Y <- data[[root_trait]]  # Now Y is a root trait
  
  # LOOCV 
  lasso_model <- cv.glmnet(X, Y, alpha = 1, nfolds = nrow(X))
  best_lambda <- lasso_model$lambda.min
  
  coefficients <- coef(lasso_model, s = best_lambda)
  important_vars <- data.frame(Variable = rownames(coefficients), Coefficient = coefficients[, 1])
  
  important_vars <- important_vars[important_vars$Coefficient != 0, ]
  important_vars <- important_vars[important_vars$Variable != "(Intercept)", ]
  
  selected_vars <- colnames(X)[colnames(X) %in% important_vars$Variable]
  
  # Check if important variables were selected
  if (length(selected_vars) > 1) {  
    # Code for plotting with significance, as before
    
    correlation_results <- rcorr(as.matrix(data[, selected_vars]), Y)
    corr_matrix <- correlation_results$r
    p_matrix <- correlation_results$P
    
    diag(p_matrix) <- NA
    
    p_matrix_significance <- p_matrix
    p_matrix_significance[p_matrix < 0.01] <- '**'
    p_matrix_significance[p_matrix >= 0.01 & p_matrix < 0.05] <- '*'
    p_matrix_significance[p_matrix >= 0.05] <- ''
    
    print(p_matrix_significance)
    
    for (var in selected_vars) {
      df_plot <- data.frame(AboveGroundTrait = data[[var]], RootTrait = Y)
      
      significance_label <- as.character(p_matrix_significance[var, "y"]) 
      
      lm_fit <- lm(RootTrait ~ AboveGroundTrait, data = df_plot)
      intercept <- coef(lm_fit)[1]
      slope <- coef(lm_fit)[2]
      equation <- paste0("y = ", round(intercept, 2), " + ", round(slope, 2), "*x")
      
      equation_with_significance <- paste0(
        ifelse(significance_label != "", 
               paste0("<span style='color:red;'>", significance_label, "</span> "), ""),
        equation)
      
      p <- ggplot(df_plot, aes(x = AboveGroundTrait, y = RootTrait)) +
        geom_point(size = 3, color = "blue") +  
        geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +  
        labs(title = paste("Regression Line for", root_trait, "and", var),
             x = var, y = root_trait) +
        geom_richtext(aes(x = Inf, y = min(df_plot$RootTrait), 
                          label = equation_with_significance), hjust = 1.2, vjust = -1.2, size = 6, fill = NA, label.color = NA) +  
        theme_minimal()
      
      print(p)
      
      # Save
      filename <- paste0("regression_plot_", var, "_vs_", root_trait, ".png")
      ggsave(filename = filename, plot = p, width = 8, height = 6)
    }
    
  } else {  
    message(paste("Not important variable for root trait:", root_trait))

##Alternative way - plot as forest plot showing coefficients
root_depth_layers <- grep("total_length_", names(data), value = TRUE)
depth_levels <- paste0("total_length_", c("[638,660]", "[616,638)", "[594,616)", "[572,594)", "[550,572)", 
                                          "[528,550)", "[506,528)", "[484,506)", "[462,484)", "[440,462)", 
                                          "[418,440)", "[396,418)", "[374,396)", "[352,374)", "[330,352)", 
                                          "[308,330)", "[286,308)", "[264,286)", "[242,264)", "[220,242)", 
                                          "[198,220)", "[176,198)", "[154,176)", "[132,154)", "[110,132)", 
                                          "[88,110)", "[66,88)", "[44,66)", "[22,44)", "[0,22)"))

results_df <- data.frame(
  DepthLayer = depth_levels,
  Variable = NA,
  Coefficient = 0,
  Importance = "Not Important",
  Significance = "Not Significant"
)
for (root_trait in root_depth_layers) {
  Y <- data[[root_trait]]  # Current soil depth layer
  
  # LASSO with LOOCV
  lasso_model <- cv.glmnet(X, Y, alpha = 1, nfolds = nrow(X))
  best_lambda <- lasso_model$lambda.min
  coefficients <- coef(lasso_model, s = best_lambda)
  important_vars <- data.frame(Variable = rownames(coefficients), Coefficient = as.vector(coefficients[, 1]))
  important_vars <- important_vars[important_vars$Variable != "(Intercept)", ]
  
  if (nrow(important_vars[important_vars$Coefficient != 0, ]) > 0) {
    for (var in important_vars$Variable[important_vars$Coefficient != 0]) {
      corr_test <- rcorr(data[[var]], Y)
      p_value <- corr_test$P[1, 2]
      
      results_df <- rbind(results_df, data.frame(
        DepthLayer = root_trait,
        Variable = var,
        Coefficient = important_vars$Coefficient[important_vars$Variable == var],
        Importance = "Important",
        Significance = ifelse(p_value < 0.05, "Significant", "Not Significant")
      ))
    }
  }
}

results_df <- results_df[!(duplicated(results_df$DepthLayer) & results_df$Coefficient == 0), ]

# order depth
results_df$DepthLayer <- factor(results_df$DepthLayer, levels = depth_levels)
depth_labels <- paste0(seq(19.8, by = -0.66, length.out = length(depth_levels)), "-", 
                       seq(20.46, by = -0.66, length.out = length(depth_levels)))

# plot
forest_plot <- ggplot(results_df, aes(x = DepthLayer, y = Coefficient)) +
  geom_point(aes(color = interaction(Importance, Significance)), size = 3) +
  geom_text(aes(label = Variable), vjust = -0.5, size = 3) +  # Add trait labels above coefficient
  geom_errorbar(aes(ymin = Coefficient - 0.05, ymax = Coefficient + 0.05), width = 0.2) +
  scale_color_manual(values = c(
    "Important.Significant" = "red",
    "Important.Not Significant" = "black",
    "Not Important.Significant" = "blue",
    "Not Important.Not Significant" = "grey"
  )) +
  scale_x_discrete(labels = depth_labels) +
  labs(title = "Forest Plot of LASSO Coefficients Across Soil Depth Layers",
       x = "Soil Depth Layer",
       y = "Coefficient") +
  theme_minimal() +
  coord_flip()

print(forest_plot)









