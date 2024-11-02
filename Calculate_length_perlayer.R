#Get ready
library(jsonlite)
library(plotly)
library(tidyr)
library(dplyr)
library(data.table)
library(htmlwidgets)
library(RColorBrewer)
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
