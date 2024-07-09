##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
## 1.1 Loading the libraries and defining ROI ----
# ---------------------------------------------------------------- - - -
library(terra)
library(stringr)
library(ggplot2)
library(dplyr) # Added for data manipulation
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 ROI ----
# ---------------------------------------------------------------- - - -
SITES = c("Paragominas") # If you want to loop through multiple sites, add them to this vector

##%######################################################%##
#                                                          #
#                    2. Loading data                      ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  # Site
  Site = SITES[i]
  cat("\n Starting processing for", Site, "\n")
  
  ## 2.1 Loading data and output directories ----
  # ---------------------------------------------------------------- - - -
  ROI <- list.files(str_glue("data/{Site}"), pattern = "shp", full.names = TRUE)
  EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")
  DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_AObs_deg/DEG_pct")
  OUT_FOLDER <- str_glue("results/{Site}/Productivity/Plots/{Site}_Line_area_km2")
  
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  ## 2.2 Charging the ROI ----
  # ---------------------------------------------------------
  roi <- vect(ROI)
  
  ## 2.3 Loading the rasters ----
  # ---------------------------------------------------------
  # Loading EVI rasters
  evi_rasters <- list(
    evi_class_2000_2005 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2000_2005.tif"))),
    evi_class_2003_2008 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2003_2008.tif"))),
    evi_class_2006_2011 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2006_2011.tif"))),
    evi_class_2009_2014 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2009_2014.tif"))),
    evi_class_2012_2017 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2012_2017.tif"))),
    evi_class_2015_2020 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2015_2020.tif")))
  )
  
  # Loading degradation rasters (accumulated and time window)
  deg_acc_rasters <- list(
    deg_1990_2005 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2005.tif"))),
    deg_1990_2008 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2008.tif"))),
    deg_1990_2011 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2011.tif"))),
    deg_1990_2014 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2014.tif"))),
    deg_1990_2017 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2017.tif"))),
    deg_1990_2020 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2020.tif")))
  )
  
  deg_tw_rasters <- list(
    deg_2000_2005 = rast(file.path(DEG_FOLDER, str_glue("{Site}_TW_deg_2000_2005.tif"))),
    deg_2003_2008 = rast(file.path(DEG_FOLDER, str_glue("{Site}_TW_deg_2003_2008.tif"))),
    deg_2006_2011 = rast(file.path(DEG_FOLDER, str_glue("{Site}_TW_deg_2006_2011.tif"))),
    deg_2009_2014 = rast(file.path(DEG_FOLDER, str_glue("{Site}_TW_deg_2009_2014.tif"))),
    deg_2012_2017 = rast(file.path(DEG_FOLDER, str_glue("{Site}_TW_deg_2012_2017.tif"))),
    deg_2015_2020 = rast(file.path(DEG_FOLDER, str_glue("{Site}_TW_deg_2015_2020.tif")))
  )
  
  ## 2.4 Cropping rasters ----
  # ---------------------------------------------------------
  crop_with_roi <- function(raster, roi) {
    roi_proj <- project(roi, crs(raster))
    cropped_raster <- crop(raster, roi_proj)
    return(mask(cropped_raster, roi_proj))
  }
  
  evi_rasters <- lapply(evi_rasters, crop_with_roi, roi)
  deg_acc_rasters <- lapply(deg_acc_rasters, crop_with_roi, roi)
  deg_tw_rasters <- lapply(deg_tw_rasters, crop_with_roi, roi)
  
  # Renaming the raster values
  for (name in names(evi_rasters)) {
    names(evi_rasters[[name]]) <- "values"
  }
  
  for (name in names(deg_acc_rasters)) {
    names(deg_acc_rasters[[name]]) <- "values"
  }
  
  for (name in names(deg_tw_rasters)) {
    names(deg_tw_rasters[[name]]) <- "values"
  }
  
  ##%######################################################%##
  #                                                          #
  #                3. Masking rasters                        ----
  #                                                          #
  ##%######################################################%##
  ## 3.1 Function to create EVI masks ----
  # ---------------------------------------------------------------- - - -
  # Masking positive values
  mask_positive <- function(raster) {
    subst(raster, c(0, -1), NA, others = TRUE)
  }
  
  # Masking negative values
  mask_negative <- function(raster) {
    subst(raster, c(0, 1), NA, others = TRUE)
  }
  
  ## 3.2 Applying masks to rasters ----
  # ---------------------------------------------------------------- - - -
  apply_masks <- function(evi_raster) {
    masked_pos <- mask(evi_raster, mask_positive(evi_raster), maskvalues = c(0, NA), updatevalue = NA)
    masked_neg <- mask(evi_raster, mask_negative(evi_raster), maskvalues = c(0, NA), updatevalue = NA)
    
    list(pos = masked_pos, neg = masked_neg)
  }
  
  ##%######################################################%##
  #                                                          #
  #                4. Calculating surface                    ----
  #                                                          #
  ##%######################################################%##
  ## 4.1 Function to calculate the area in km² for each category ----
  # ---------------------------------------------------------------- - - -
  # Function to calculate the area in square meters of a single pixel at a given latitude
  area_of_pixel <- function(lat, res) {
    earth_radius <- 6378137 # Earth's radius in meters
    lat_rad <- lat * pi / 180 # Convert latitude to radians
    meter_per_deg_long <- pi * earth_radius * cos(lat_rad) / 180 # Meters per degree longitude
    meter_per_deg_lat <- 111320 # Approximate meters per degree latitude
    area_m2 <- (meter_per_deg_long * res[1]) * (meter_per_deg_lat * res[2])
    return(area_m2)
  }
  
  # Calculate area in km² for EVI categories
  calculate_area_evi <- function(raster) {
    center_lat <- (ymax(raster) + ymin(raster)) / 2
    resolution <- res(raster)
    cell_area_m2 <- area_of_pixel(center_lat, resolution)
    
    print(paste("Cell Area m2:", cell_area_m2))
    
    freq <- freq(raster)
    print(freq)
    
    freq_df <- as.data.frame(freq)
    freq_df$area_m2 <- freq_df$count * cell_area_m2
    freq_df$area_km2 <- freq_df$area_m2 / 1e6
    return(freq_df)
  }
  
  # Function to calculate the area in km² for each category of degradation (full pixel)
  calculate_area_full_pixel <- function(raster) {
    res_degrees <- res(raster) # Get the resolution in degrees
    res_meters <- res_degrees * 111320 # Convert resolution from degrees to meters
    cell_area_m2 <- res_meters[1] * res_meters[2]  # Calculate the area of each cell in square meters
    count_cells <- sum(values(raster > 0), na.rm = TRUE)  # Count the number of cells with value greater than 0
    total_area_km2 <- (count_cells * cell_area_m2) / 1e6  # Convert m² to km²
    
    return(total_area_km2)
  }
  
  ##%######################################################%##
  #                                                          #
  #                5. Preparing data for plotting            ----
  #                                                          #
  ##%######################################################%##
  ## 5.1 Function to prepare EVI data ----
  # ---------------------------------------------------------------- - - -
  # Prepare EVI data for plotting
  prepare_evi_data <- function(masked_list, period) {
    pos_area <- calculate_area_evi(masked_list$pos)
    pos_area$category <- "Positive"
    pos_area$period <- period
    
    neg_area <- calculate_area_evi(masked_list$neg)
    neg_area$category <- "Negative"
    neg_area$period <- period
    
    rbind(pos_area, neg_area)
  }
  
  ## 5.2 Function to prepare degradation data ----
  # ---------------------------------------------------------------- - - -
  # Prepare degradation data for plotting
  prepare_degradation_data <- function(rasters_acc, rasters_tw, periods) {
    data_list_acc <- list()
    data_list_tw <- list()
    
    for (i in 1:length(rasters_acc)) {
      area_acc <- calculate_area_full_pixel(rasters_acc[[i]])
      data_list_acc[[i]] <- data.frame(period = periods[i], area_km2 = area_acc, type = "Accumulated")
    }
    
    for (i in 1:length(rasters_tw)) {
      area_tw <- calculate_area_full_pixel(rasters_tw[[i]])
      data_list_tw[[i]] <- data.frame(period = periods[i], area_km2 = area_tw, type = "Window")
    }
    
    data_acc <- do.call(rbind, data_list_acc)
    data_tw <- do.call(rbind, data_list_tw)
    
    data <- rbind(data_acc, data_tw)
    return(data)
  }
  
  # Define periods for EVI and degradation
  periods_evi <- c("2000-2005", "2003-2008", "2006-2011", "2009-2014", "2012-2017", "2015-2020")
  periods_degradation <- c("2000-2005", "2003-2008", "2006-2011", "2009-2014", "2012-2017", "2015-2020")
  
  # Apply masks to EVI rasters and prepare data for plotting
  masked_evi_rasters <- lapply(evi_rasters, apply_masks)
  line_plot_data_evi <- mapply(prepare_evi_data, masked_evi_rasters, periods_evi, SIMPLIFY = FALSE)
  line_plot_data_evi <- do.call(rbind, line_plot_data_evi)
  line_plot_data_evi <- line_plot_data_evi[line_plot_data_evi$category %in% c("Positive", "Negative"), ]
  line_plot_data_evi$category <- factor(line_plot_data_evi$category, levels = c("Positive", "Negative"))
  
  # Prepare degradation data for plotting
  line_plot_data_degradation <- prepare_degradation_data(deg_acc_rasters, deg_tw_rasters, periods_degradation)
  
  ##%######################################################%##
  #                                                          #
  #                6. Creating the combined plot             ----
  #                                                          #
  ##%######################################################%##
  combined_line_plot <- function(data_evi, data_degradation, output_file) {
    data_evi$period <- as.factor(data_evi$period)
    data_degradation$period <- as.factor(data_degradation$period)
    
    ggplot() +
      geom_line(data = data_evi, aes(x = period, y = area_km2, color = category, group = category), linewidth = 3, alpha = 0.6, lineend = "round") +
      geom_line(data = data_degradation, aes(x = period, y = area_km2, color = type, group = type), linewidth = 3, alpha = 0.6, lineend = "round") +
      scale_color_manual(name = " ",
                         values = c("Positive" = "#1a9850", "Negative" = "#d73027", "Accumulated" = "#1f78b4", "Window" = "#6a3d9a"),
                         labels = c("Positive" = "Tendance EVI positive", "Negative" = "Tendance EVI Negative", "Accumulated" = "Deg. accumulée (1990 - )", "Window" = "Deg. par période")) +
      labs(title = "Surface des tendances de l'EVI et de la dégradation",
           x = "Période",
           y = "Surface (km²)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 20)),
        axis.title.x = element_text(size = 20, margin = margin(t = 10)),
        axis.title.y = element_text(size = 20, margin = margin(r = 10)),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 18, hjust = 0.5), 
        legend.position = "bottom",
        legend.box.spacing = unit(0.5, 'cm'),
        legend.spacing.x = unit(1, 'cm'), 
        legend.key.width = unit(1.5, 'cm'),  
        legend.key.height = unit(0.5, 'cm'), 
        legend.margin = margin(t = 10, b = 10, l = 10, r = 10)
      ) +
      scale_x_discrete(expand = c(0, 0.2)) +  # Adds small extra space at the ends of the X axis
      scale_y_continuous(expand = c(0.05, 0)) +  # Adds small extra space at the ends of the Y axis
      guides(color = guide_legend(override.aes = list(linetype = c(rep(1, 4)))))
    
    ggsave(output_file, width = 16, height = 6, dpi = 300)  # Increase height to prevent line cutting
  }
  
  # Create the combined plot
  combined_line_plot(line_plot_data_evi, line_plot_data_degradation, file.path(OUT_FOLDER, "01_EVI_Degradation_Area.png"))
  
  ##%######################################################%##
  #                                                          #
  #                7. Exporting the data to CSV              ----
  #                                                          #
  ##%######################################################%##
  write.csv(line_plot_data_evi, file.path(OUT_FOLDER, "01_EVI_area_km2_data.csv"), row.names = FALSE)
  write.csv(line_plot_data_degradation, file.path(OUT_FOLDER, "01_Degradation_area_km2_data.csv"), row.names = FALSE)
  
  cat("\n End of processing for", Site, "\n")
}
