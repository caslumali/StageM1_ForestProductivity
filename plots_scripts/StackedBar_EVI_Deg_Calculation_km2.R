##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
## 1.1 Loading the libraries and defining ROI ----
# ---------------------------------------------------------
# Install packages if not already installed
required_packages <- c("terra", "stringr", "ggplot2", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(terra)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Define terra options
terraOptions(tempdir = "D:/temp")

## 1.2 ROI ----
# ---------------------------------------------------------
SITES <- c("Paragominas") # Add more sites as needed

##%######################################################%##
#                                                          #
#                    2. Loading data                      ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  Site <- SITES[i]
  cat("\n Starting processing for", Site, "\n")
  
  ## 2.1 Loading data and output directories ----
  # ---------------------------------------------------------
  ROI <- list.files(str_glue("data/{Site}"), pattern = "shp", full.names = TRUE)
  EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")
  DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_AObs_deg/DEG_pct")
  OUT_FOLDER <- str_glue("results/{Site}/Productivity/Plots/{Site}_StackedBar_area_km2")
  
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  ## 2.2 Charging the ROI ----
  # ---------------------------------------------------------
  roi <- vect(ROI)
  
  ## 2.3 Loading the rasters ----
  # ---------------------------------------------------------
  evi_rasters <- list(
    evi_class_2000_2005 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2000_2005.tif"))),
    evi_class_2003_2008 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2003_2008.tif"))),
    evi_class_2006_2011 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2006_2011.tif"))),
    evi_class_2009_2014 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2009_2014.tif"))),
    evi_class_2012_2017 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2012_2017.tif"))),
    evi_class_2015_2020 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2015_2020.tif")))
  )
  
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
  
  ##%######################################################%##
  #                                                          #
  #                3. Calculate Statistics                  ----
  #                                                          #
  ##%######################################################%##
  cat("\n Calculating EVI stats for", Site, "\n")
  
  area_of_pixel <- function(lat, res) {
    earth_radius <- 6378137
    lat_rad <- lat * pi / 180
    meter_per_deg <- pi * earth_radius * cos(lat_rad) / 180
    area_m2 <- (meter_per_deg * res[1]) * (111320 * res[2])
    return(area_m2)
  }
  
  calculate_evi_stats <- function(raster, raster_name) {
    freq_values <- freq(raster)
    total_pixels <- sum(freq_values$count, na.rm = TRUE)
    pos_pixels <- sum(freq_values$count[freq_values$value == 1], na.rm = TRUE)
    neg_pixels <- sum(freq_values$count[freq_values$value == -1], na.rm = TRUE)
    nonsig_pixels <- sum(freq_values$count[freq_values$value == 0], na.rm = TRUE)
    
    center_lat <- (ymax(raster) + ymin(raster)) / 2
    resolution <- res(raster)
    pixel_area_m2 <- area_of_pixel(center_lat, resolution)
    
    total_area_km2 <- (total_pixels * pixel_area_m2) / 1e6
    pos_area_km2 <- (pos_pixels * pixel_area_m2) / 1e6
    neg_area_km2 <- (neg_pixels * pixel_area_m2) / 1e6
    nonsig_area_km2 <- (nonsig_pixels * pixel_area_m2) / 1e6
    
    results <- data.frame(
      Raster = raster_name,
      Total_Pixels = total_pixels,
      Total_Area_Km2 = total_area_km2,
      Positive_Pixels = pos_pixels,
      Positive_Area_Km2 = pos_area_km2,
      Negative_Pixels = neg_pixels,
      Negative_Area_Km2 = neg_area_km2,
      NonSig_Pixels = nonsig_pixels,
      NonSig_Area_Km2 = nonsig_area_km2
    )
    
    print(results)
    evi_stats_all <<- rbind(evi_stats_all, results)
  }
  
  calculate_deg_stats <- function(raster, raster_name, type) {
    freq_values <- freq(raster)
    total_pixels <- sum(freq_values$count, na.rm = TRUE)
    nodeg_pixels <- sum(freq_values$count[freq_values$value == 0], na.rm = TRUE)
    deg_pixels <- total_pixels - nodeg_pixels
    
    center_lat <- (ymax(raster) + ymin(raster)) / 2
    resolution <- res(raster)
    pixel_area_m2 <- area_of_pixel(center_lat, resolution)
    
    total_area_km2 <- (total_pixels * pixel_area_m2) / 1e6
    nodeg_area_km2 <- (nodeg_pixels * pixel_area_m2) / 1e6
    deg_area_km2 <- (deg_pixels * pixel_area_m2) / 1e6
    
    results <- data.frame(
      Raster = raster_name,
      Total_Pixels = total_pixels,
      Total_Area_Km2 = total_area_km2,
      NoDeg_Pixels = nodeg_pixels,
      NoDeg_Area_Km2 = nodeg_area_km2,
      Deg_Pixels = deg_pixels,
      Deg_Area_Km2 = deg_area_km2,
      Type = type
    )
    
    print(results)
    if (type == "acc") {
      acc_deg_stats_all <<- rbind(acc_deg_stats_all, results)
    } else {
      tw_deg_stats_all <<- rbind(tw_deg_stats_all, results)
    }
  }
  
  evi_stats_all <- data.frame()
  acc_deg_stats_all <- data.frame()
  tw_deg_stats_all <- data.frame()
  
  for (name in names(evi_rasters)) {
    calculate_evi_stats(evi_rasters[[name]], name)
  }
  
  for (name in names(deg_acc_rasters)) {
    calculate_deg_stats(deg_acc_rasters[[name]], name, "acc")
  }
  
  for (name in names(deg_tw_rasters)) {
    calculate_deg_stats(deg_tw_rasters[[name]], name, "tw")
  }
  
  write.csv(evi_stats_all, file.path(OUT_FOLDER, "EVI_stats_all.csv"), row.names = FALSE)
  write.csv(acc_deg_stats_all, file.path(OUT_FOLDER, "Acc_deg_stats_all.csv"), row.names = FALSE)
  write.csv(tw_deg_stats_all, file.path(OUT_FOLDER, "Tw_deg_stats_all.csv"), row.names = FALSE)
  
  ##%######################################################%##
  #                                                          #
  #                 4. Stacked Bar plots                    ----
  #                                                          #
  ##%######################################################%##
  stacked_bar_plot <- function(data, output_file, title, x_labels_order, x_labels) {
    ggplot(data, aes(x = factor(Raster, levels = x_labels_order), y = Area_Km2, fill = category)) +
      geom_bar(stat = "identity", position = "stack", alpha = 0.6, width = 0.6) +  # Diminuir a largura das barras
      scale_fill_manual(
        values = c("Positive_Area_Km2" = "#1a9850", "Negative_Area_Km2" = "#d73027", "NonSig_Area_Km2" = "grey", "NoDeg_Area_Km2" = "grey", "Deg_Area_Km2" = "#1f78b4"), # blue - "#1f78b4" | purple = "#6a3d9a"
        labels = c("Positive_Area_Km2" = "EVI positive", "Negative_Area_Km2" = "EVI Negative", "NonSig_Area_Km2" = "Non significatif", "NoDeg_Area_Km2" = "Non dégradée", "Deg_Area_Km2" = "Degradée")
      ) +
      labs(title = title, x = "Période", y = "Area (km²)") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 20)),
        axis.title.x = element_text(size = 18, margin = margin(t = 10)),
        axis.title.y = element_text(size = 18, margin = margin(r = 10)),
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
      scale_x_discrete(labels = x_labels, expand = c(0.09, 0.01))  +
      scale_y_continuous(expand = c(0.05, 0)) +
      guides(fill = guide_legend(override.aes = list(alpha = 0.6)))
    
    ggsave(output_file, width = 14, height = 6, dpi = 300)
  }
  
  evi_plot_data <- evi_stats_all %>%
    gather(key = "category", value = "Area_Km2", Positive_Area_Km2, Negative_Area_Km2, NonSig_Area_Km2) %>%
    mutate(category = factor(category, levels = c("NonSig_Area_Km2", "Negative_Area_Km2", "Positive_Area_Km2")))
  
  evi_periods <- c("evi_class_2000_2005", "evi_class_2003_2008", "evi_class_2006_2011", "evi_class_2009_2014", "evi_class_2012_2017", "evi_class_2015_2020")
  evi_labels <- c("2000 à 2005", "2003 à 2008", "2006 à 2011", "2009 à 2014", "2012 à 2017", "2015 à 2020")
  stacked_bar_plot(evi_plot_data, file.path(OUT_FOLDER, "EVI_stacked_bar_plot.png"), "Surface de la tendance de l'EVI (km²)", evi_periods, evi_labels)
  
  acc_deg_plot_data <- acc_deg_stats_all %>%
    gather(key = "category", value = "Area_Km2", NoDeg_Area_Km2, Deg_Area_Km2) %>%
    mutate(category = factor(category, levels = c("NoDeg_Area_Km2", "Deg_Area_Km2")))
  
  
  acc_deg_periods <- c("deg_1990_2005", "deg_1990_2008", "deg_1990_2011", "deg_1990_2014", "deg_1990_2017", "deg_1990_2020")
  acc_deg_labels <- c("1990 à 2005", "1990 à 2008", "1990 à 2011", "1990 à 2014", "1990 à 2017", "1990 à 2020")
  stacked_bar_plot(acc_deg_plot_data, file.path(OUT_FOLDER, "Acc_deg_stacked_bar_plot.png"), "Surface de pixels présentant une certaine dégradation (accumulée)", acc_deg_periods, acc_deg_labels)
  
  tw_deg_plot_data <- tw_deg_stats_all %>%
    gather(key = "category", value = "Area_Km2", NoDeg_Area_Km2, Deg_Area_Km2) %>%
    mutate(category = factor(category, levels = c("NoDeg_Area_Km2", "Deg_Area_Km2")))
  
  tw_deg_periods <- c("deg_2000_2005", "deg_2003_2008", "deg_2006_2011", "deg_2009_2014", "deg_2012_2017", "deg_2015_2020")
  tw_deg_labels <- c("2000 à 2005", "2003 à 2008", "2006 à 2011", "2009 à 2014", "2012 à 2017", "2015 à 2020")
  stacked_bar_plot(tw_deg_plot_data, file.path(OUT_FOLDER, "Tw_deg_stacked_bar_plot.png"), "Surface de pixels présentant une certaine dégradation (par période)", tw_deg_periods, tw_deg_labels)
  
  cat("\n End of processing for", Site, "\n")
}
