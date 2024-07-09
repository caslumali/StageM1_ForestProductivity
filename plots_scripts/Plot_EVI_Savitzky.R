##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
# Defining ROI variable
SITES = c("Paragominas") # If you want to loop # 'Cotriguacu', 'Guaviare', 'MDD', 'Paragominas'
# Site <- "Paragominas"
WRITE_INTERMEDIATE <- TRUE

## 1.1 Charging the libraries ----
# ---------------------------------------------------------------- - - -
library(terra)
library(data.table)
library(stringr)
library(ggplot2)
terraOptions(tempdir="D:/temp")
terraOptions()

##%######################################################%##
#                                                          #
#                    2. Prepare data                       ----
#                                                          #
##%######################################################%##
for ( i in 1:length(SITES)){
  # Site
  Site = SITES[i]
  cat("\n Starting processing for", Site, "\n")
  
  ## 2.1 Defining data and output directories ----
  # ---------------------------------------------------------------- - - -
  DATA_FOLDER <- stringr::str_glue("results/{Site}/Productivity/EVI")
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/Plots/Savitzky-Golay_cleaned")
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  ## 2.2 Load the names and renaming the rasters if metrics were not generated ----
  # ---------------------------------------------------------------- - - -
  file_names <- list.files(stringr::str_glue("data/{Site}/MODIS_NASA/MODIS_EVI"))
  doy_values <- str_extract(file_names, "doy[0-9]{7}")
  date_files <- as.Date(doy_values, format = "doy%Y%j")
  
  # Manually add the missing date for Paragominas
  if (Site == 'Paragominas'){
    missing_date <- as.Date("2005-08-29") # 241º dia de 2005
    if (!missing_date %in% date_files) {
      date_files <- sort(c(date_files, missing_date))
    }
  }
  
  evi_fc_masked <- rast(file.path(DATA_FOLDER, "evi_fc_masked.tif"))
  names(evi_fc_masked) <- date_files
  
  evi_smoothed <- rast(file.path(DATA_FOLDER, "evi_smoothed.tif"))
  names(evi_smoothed) <- date_files

  
  ## 2.4 Load the metrics if they have already been generated ----
  # ---------------------------------------------------------------- - - -
  mean_evi_original <- readRDS(file.path(OUT_FOLDER, str_glue("{Site}_mean_evi_original.rds")))
  mean_evi_smoothed <- readRDS(file.path(OUT_FOLDER, str_glue("{Site}_mean_evi_smoothed.rds")))
  
  ##%######################################################%##
  #                                                          #
  #                     3. EVI smoothed                     ----
  #                                                          #
  ##%######################################################%##
  cat("\n Calculating means and plotting for", Site, "\n")
  
  # 3.1 Formatting smoothed EVI data for plotting ----
  # ---------------------------------------------------------------- - - -
  # Extract the mean values for each original layer
  mean_evi_original <- sapply(1:nlyr(evi_fc_masked), function(i) {
    mean(values(evi_fc_masked[[i]]), na.rm = TRUE)
  })
  
  # Extract the mean EVI value for each smoothed layer
  mean_evi_smoothed <- sapply(1:nlyr(evi_smoothed), function(i) {
    mean(values(evi_smoothed[[i]]), na.rm = TRUE)
  })
  
  if(WRITE_INTERMEDIATE) {
    saveRDS(mean_evi_original, file.path(OUT_FOLDER, str_glue("{Site}_mean_evi_original.rds")))
    saveRDS(mean_evi_smoothed, file.path(OUT_FOLDER, str_glue("{Site}_mean_evi_smoothed.rds")))
  }
  
  # Create a data frame for plotting with dates and EVI values
  df_evi_original <- data.table(
    Date = as.Date(names(evi_fc_masked), format = "%Y-%m-%d"),
    Original_EVI = mean_evi_original
  )
  
  # Create a data frame for the smoothed line
  df_evi_smoothed <- data.table(
    Date = as.Date(names(evi_smoothed), format = "%Y-%m-%d"),
    Smoothed_EVI = mean_evi_smoothed
  )
  
  ## 3.2 Plotting and saving in French ----
  # ---------------------------------------------------------------- - - -
  evi_smoothed_plot_fr <- ggplot() +
    geom_point(data = df_evi_original, aes(x = Date, y = Original_EVI, color = "EVI moyen original"), 
               shape = 1, size = 2, stroke = 1) +
    geom_line(data = df_evi_smoothed, aes(x = Date, y = Smoothed_EVI, color = "EVI moyen lissé"), 
              linewidth = 1) +
    labs(title = "MODIS EVI 250m",
         subtitle = paste({Site} , " - période 2000 - 2020"),
         y = "EVI", x = "", color = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(color = "#666666", size = 22, hjust = 0.5, face = "bold",
                                    margin = margin(t = 5, b = 10), vjust = 2),
          plot.subtitle = element_text(color = "#666666", size = 16, hjust = 0.5,
                                       margin = margin(-10, 0, 10, 0)),
          axis.title.y = element_text(face = "bold", size = 20, margin = margin(l = 10, r = 10)),
          axis.text.x = element_text(vjust = 0.5, size = 12, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.position = "bottom", 
          legend.text = element_text(size = 14),
          legend.margin = margin(t = 0, b = 0, l = 0, r = 0), 
          legend.box.margin = margin(t = -10, b = 5, l = 0, r = 0),
          legend.spacing = unit(0.2, 'cm'),
          legend.key.width = unit(1.5, 'lines')) + 
    scale_color_manual(values = c("EVI moyen original" = "#51A3DB", "EVI moyen lissé" = "#176396"))
  print(evi_smoothed_plot_fr) # Display the plot
  
  ggsave(file.path(OUT_FOLDER, "evi_smoothed_fr.png"), evi_smoothed_plot_fr, width = 30,
         height = 12, units = "cm", dpi = 300)
  
  
  ## 3.3 Plotting and saving in English ----
  # ---------------------------------------------------------------- - - -
  evi_smoothed_plot_eng <- ggplot() +
    geom_point(data = df_evi_original, aes(x = Date, y = Original_EVI, color = "Original avg. EVI"), 
               shape = 1, size = 2, stroke = 1) +
    geom_line(data = df_evi_smoothed, aes(x = Date, y = Smoothed_EVI, color = "Smoothed avg. EVI"), 
              linewidth = 1) +
    labs(title = "MODIS EVI 250m",
         subtitle = paste({Site} , " - 2000 - 2020 period"),
         y = "EVI", x = "", color = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(color = "#666666", size = 22, hjust = 0.5, face = "bold",
                                    margin = margin(t = 5, b = 10), vjust = 2),
          plot.subtitle = element_text(color = "#666666", size = 16, hjust = 0.5,
                                       margin = margin(-10, 0, 10, 0)),
          axis.title.y = element_text(face = "bold", size = 20, margin = margin(l = 10, r = 10)),
          axis.text.x = element_text(vjust = 0.5, size = 12, hjust = 1),
          axis.text.y = element_text(size = 12),
          legend.position = "bottom", 
          legend.text = element_text(size = 14),
          legend.margin = margin(t = 0, b = 0, l = 0, r = 0), 
          legend.box.margin = margin(t = -10, b = 5, l = 0, r = 0),
          legend.spacing = unit(0.2, 'cm'),
          legend.key.width = unit(1.5, 'lines')) + 
    scale_color_manual(values = c("Original avg. EVI" = "#51A3DB", "Smoothed avg. EVI" = "#176396"))
  print(evi_smoothed_plot_eng)
  
  # Save the plot
  ggsave(file.path(OUT_FOLDER, "evi_smoothed_eng.png"), evi_smoothed_plot_eng, width = 30,
         height = 12, units = "cm", dpi = 300)
}
