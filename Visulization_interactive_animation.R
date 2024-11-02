#Load toolkits & Set up environment
library(jsonlite)
library(plotly)
library(dplyr)
library(data.table)
library(htmlwidgets)
library(RColorBrewer)
setwd("/data/gpfs/projects/punim1869/users/yunhongh1/workspace/intergrain/NARO_3D/CT_cb/rinfo")
# Define the root type color mapping
root_type_colors <- list(
  "SNR" = "#FF4500",
  "LNR" = "#1E90FF",  
  "SWR" = "#FFA500",  
  "LWR" = "#228B22"   
)
get_root_color <- function(root_type) {
  return(root_type_colors[[root_type]])
}

# JSON -> Polylines
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

# List of data files and corresponding conditions
files_and_conditions <- list(
  # Surplus with 3 reps for each genotype
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
  
  # Deficient with 3 reps for each genotype
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

for (file_info in files_and_conditions) {
  
  # json + root type data file
  json_file <- file_info[[1]]
  data_table_file <- file_info[[2]]
  json_data <- fromJSON(json_file, simplifyDataFrame = FALSE)
  data_table <- fread(data_table_file)
  
  # Extract polylines
  polylines_list <- list()
  polylines_list <- extract_polylines(json_data, polylines_list)
  
  # Data frame for each polyline then assign the root type
  polylines_df <- lapply(1:length(polylines_list), function(i) {
    polyline <- polylines_list[[i]]
    root_type <- data_table$Root_Type[i]  # Root_Type as scalar value
    color <- get_root_color(root_type)  # root type -> colors
    
    if (!is.null(polyline) && ncol(polyline) == 3) {
      return(list(
        df = data.frame(x = polyline[, 1], y = polyline[, 2], z = polyline[, 3]),
        color = color
      ))
    } else {
      return(NULL)  
    }
  })
  
  polylines_df <- polylines_df[!sapply(polylines_df, is.null)]
  
  # plot - interactive plotting 
  p <- plot_ly(type = 'scatter3d', mode = 'lines')
  
  # get colours for each ployline (root)
  for (polyline_obj in polylines_df) {
    poly_df <- polyline_obj$df
    color <- polyline_obj$color
    
    if (!is.null(poly_df)) {
      p <- p %>%
        add_trace(
          x = poly_df$x, 
          y = poly_df$y, 
          z = poly_df$z, 
          type = 'scatter3d', 
          mode = 'lines',
          line = list(color = color, width = 2) 
        )
    }
  }
  
  # save in html (interactive)
  file_name_html <- paste0(tools::file_path_sans_ext(data_table_file), ".html")
  
  p <- p %>%
    layout(scene = list(
      xaxis = list(title = 'X', range = c(0, 700)),  # Set fixed range for X axis
      yaxis = list(title = 'Y', range = c(0, 700)),  # Set fixed range for Y axis
      zaxis = list(title = 'Z', range = c(0, 700))   # Set fixed range for Z axis
    ))
  
  saveWidget(p, file_name_html)
  
  cat("Plot saved as:", file_name_html, "\n")  # check if finished
}

#plot - Animations - choose good ones for PPT

  # update p to animation
  p <- plot_ly()
  num_lines <- length(polylines_df)
  colors <- brewer.pal(min(9, num_lines), "Set1") 
  
  # Get the maximum X value for the animation step range
  max_x <- max(sapply(polylines_df, function(df) max(df$df$x, na.rm = TRUE)))
  
  for (i in seq_along(polylines_df)) {
    root_data <- polylines_df[[i]]$df
    color <- polylines_df[[i]]$color
    
    # sort by x (depth - so strange... change coordination x to z in the future...)
    root_data <- root_data[order(root_data$x), ]
    
    for (step in seq_len(max_x)) {
      if (step <= max(root_data$x)) {
        frame_x <- root_data$x[root_data$x <= step]
        frame_y <- root_data$y[seq_along(frame_x)]
        frame_z <- root_data$z[seq_along(frame_x)]
        
        p <- add_trace(
          p,
          type = 'scatter3d',
          mode = 'lines',
          x = frame_x,
          y = frame_y,
          z = frame_z,
          line = list(color = color, width = 5),
          frame = step,  
          showlegend = FALSE
        )
      } else {
        # Add full line once the frame exceeds the max x for that root
        p <- add_trace(
          p,
          type = 'scatter3d',
          mode = 'lines',
          x = root_data$x,
          y = root_data$y,
          z = root_data$z,
          line = list(color = color, width = 5),
          frame = step,  
          showlegend = FALSE
        )
      }
    }
  }
  
  # max layer - 644 so 700 for consistancy and look more real as a pot (small barley root :( )
  p <- layout(
    p,
    scene = list(
      xaxis = list(title = "X", range = c(0, 700)),
      yaxis = list(title = "Y", range = c(0, 700)),
      zaxis = list(title = "Z", range = c(0, 700))
    ),
    title = paste("3D Root Growth Animation for", tools::file_path_sans_ext(basename(data_table_file)))
  ) %>%
    animation_opts(
      frame = 100, transition = 0, redraw = TRUE
    ) %>%
    animation_slider(
      currentvalue = list(prefix = "Growth Step ")
    )
  
  # Save each animated plot as an HTML file
  file_name_html <- paste0(tools::file_path_sans_ext(data_table_file), "_animation.html")
  saveWidget(p, file_name_html)
  
  cat("Animated plot saved as:", file_name_html, "\n")  # check if any not saved
}
