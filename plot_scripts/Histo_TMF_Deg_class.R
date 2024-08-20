##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
## 1.1 Charging the libraries and defining ROI ----
# ---------------------------------------------------------------- - - -
library(terra)
library(stringr)
library(ggplot2)
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 ROI ----
# ---------------------------------------------------------------- - - -
SITES = c("Paragominas")

##%######################################################%##
#                                                          #
#                    2. Charging data                      ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  Site = SITES[i]
  
  cat("\n Starting processing for", Site, "\n")
  
  ## 2.1 Charging data and output directories ----
  # ---------------------------------------------------------------- - - -   
  EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")
  
  # DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_deg/DEG_class")
  DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_AnnualObs/DEG_class")
  
  # Defining and creating output folders
  # OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/Plots/{Site}_Histo_DegClass")
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/Plots/{Site}_Histo_AnnualObsClass")
  
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUT_FOLDER, "Histo_AccDeg"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUT_FOLDER, "Histo_TwDeg"), showWarnings = FALSE, recursive = TRUE)
  
  
  ## 2.2 Charging the rasters ----
  # ---------------------------------------------------------------- - - -
  evi_rasters <- list(
    evi_class_2000_2005 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2000_2005.tif"))),
    evi_class_2000_2005 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2000_2005.tif"))),
    evi_class_2003_2008 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2003_2008.tif"))),
    evi_class_2006_2011 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2006_2011.tif"))),
    evi_class_2009_2014 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2009_2014.tif"))),
    evi_class_2012_2017 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2012_2017.tif"))),
    evi_class_2015_2020 = rast(file.path(EVI_FOLDER, str_glue("{Site}_evi_class_2015_2020.tif")))
  )
  
  deg_acc_rasters <- list(
    deg_1990_1999 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_1999.tif"))),
    deg_1990_2005 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_2005.tif"))),
    deg_1990_2008 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_2008.tif"))),
    deg_1990_2011 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_2011.tif"))),
    deg_1990_2014 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_2014.tif"))),
    deg_1990_2017 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_2017.tif"))),
    deg_1990_2020 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_class_1990_2020.tif")))
  )
  
  deg_tw_rasters <- list(
    deg_1990_1999 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_1990_1999.tif"))),
    deg_2000_2005 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_2000_2005.tif"))),
    deg_2003_2008 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_2003_2008.tif"))),
    deg_2006_2011 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_2006_2011.tif"))),
    deg_2009_2014 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_2009_2014.tif"))),
    deg_2012_2017 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_2012_2017.tif"))),
    deg_2015_2020 = rast(file.path(DEG_FOLDER, str_glue("{Site}_Tw_deg_class_2015_2020.tif")))
  )
  
  # Renaming the raster values
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
  ## 3.1 Function to creating evi mask ----
  # ---------------------------------------------------------------- - - -
  mask_positive <- function(raster) {
    subst(raster, c(0, -1), NA, others = TRUE)
  }
  
  mask_negative <- function(raster) {
    subst(raster, c(0, 1), NA, others = TRUE)
  }
  
  mask_ns <- function(raster) {
    subst(raster, c(-1, 1), NA, others = TRUE)
  }
  
  ## 3.2 Applying mask to rasters ----
  # ---------------------------------------------------------------- - - - 
  apply_masks <- function(evi_raster, change_raster) {
    masked_pos <- mask(change_raster, mask_positive(evi_raster), maskvalues = c(0, NA), updatevalue = NA)
    masked_neg <- mask(change_raster, mask_negative(evi_raster), maskvalues = c(0, NA), updatevalue = NA)
    masked_ns <- mask(change_raster, mask_ns(evi_raster), maskvalues = c(0, NA), updatevalue = NA)
    
    list(pos = masked_pos, neg = masked_neg, ns = masked_ns)
  }
  
  # Scrolling through the rasters and applying masks
  masked_deg_acc_rasters <- mapply(apply_masks, evi_rasters, deg_acc_rasters, SIMPLIFY = FALSE)
  masked_deg_tw_rasters <- mapply(apply_masks, evi_rasters, deg_tw_rasters, SIMPLIFY = FALSE)
  
  
  ##%######################################################%##
  #                                                          #
  #               4. Creating histograms                     ----
  #                                                          #
  ##%######################################################%##
  ## 4.1 Define colors for positive and negative values ----
  # ---------------------------------------------------------------- - - -
  positive_colors <- c("0" = "lightgrey",
                       "25" = "#a6d96a",
                       "50" = "#66bd63",
                       "75" = "#1a9850",
                       "100" = "#006837")
  
  negative_colors <- c("0" = "lightgrey",
                       "25" = "#fee08b",
                       "50" = "#fdae61",
                       "75" = "#f46d43",
                       "100" = "#d73027")
  
  ## 4.2 List to names the histograms plots and files ----
  # ---------------------------------------------------------------- - - -
  acc_titles <- c("1990-1999", "1990-2005", "1990-2008", "1990-2011", "1990-2014", "1990-2017", "1990-2020")
  tw_titles <- c("1990-1999", "2000-2005", "2003-2008", "2006-2011", "2009-2014", "2012-2017", "2015-2020")
  
  histo_deg_acc_pos_names <- c(str_glue("01-{Site}_acc_deg_1990_1999_pos"),
                               str_glue("03-{Site}_acc_deg_1990_2005_pos"),
                               str_glue("05-{Site}_acc_deg_1990_2008_pos"),
                               str_glue("07-{Site}_acc_deg_1990_2011_pos"),
                               str_glue("09-{Site}_acc_deg_1990_2014_pos"),
                               str_glue("11-{Site}_acc_deg_1990_2017_pos"),
                               str_glue("13-{Site}_acc_deg_1990_2020_pos"))
  
  histo_deg_acc_neg_names <- c(str_glue("02-{Site}_acc_deg_1990_1999_neg"),
                               str_glue("04-{Site}_acc_deg_1990_2005_neg"),
                               str_glue("06-{Site}_acc_deg_1990_2008_neg"),
                               str_glue("08-{Site}_acc_deg_1990_2011_neg"),
                               str_glue("10-{Site}_acc_deg_1990_2014_neg"), 
                               str_glue("12-{Site}_acc_deg_1990_2017_neg"),
                               str_glue("14-{Site}_acc_deg_1990_2020_neg"))
  
  histo_deg_tw_pos_names <- c(str_glue("01-{Site}_tw_deg_1990_1999_pos"),
                              str_glue("03-{Site}_tw_deg_2000_2005_pos"),
                              str_glue("05-{Site}_tw_deg_2003_2008_pos"),
                              str_glue("07-{Site}_tw_deg_2006_2011_pos"),
                              str_glue("09-{Site}_tw_deg_2009_2014_pos"),
                              str_glue("11-{Site}_tw_deg_2012_2017_pos"),
                              str_glue("13-{Site}_tw_deg_2015_2020_pos"))
  
  histo_deg_tw_neg_names <- c(str_glue("02-{Site}_tw_deg_1990_1999_neg"),
                              str_glue("04-{Site}_tw_deg_2000_2005_neg"),
                              str_glue("06-{Site}_tw_deg_2003_2008_neg"),
                              str_glue("08-{Site}_tw_deg_2006_2011_neg"),
                              str_glue("10-{Site}_tw_deg_2009_2014_neg"), 
                              str_glue("12-{Site}_tw_deg_2012_2017_neg"),
                              str_glue("14-{Site}_tw_deg_2015_2020_neg"))
  
  
  ## 4.3 Function to create the histograms ----
  # ---------------------------------------------------------------- - - -
  create_histogram <- function(df, title, output_folder, filename_base, colors) {
    if (sum(df$values, na.rm = TRUE) > 0) {
      # Create a factor with the specific change values
      df$values <- factor(df$values, levels = c(0, 25, 50, 75, 100))
      
      histo <- ggplot(df, aes(x = values)) +
        geom_bar(aes(fill = values), color = "darkgrey", alpha = 0.7) +
        scale_fill_manual(values = colors, guide = "none") +
        labs(title = title, x = "Percentage of Change", y = "Frequency of Pixels") +
        theme_minimal(base_family = "Helvetica") +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(0, 0, 20, 0)),
          plot.subtitle = element_text(margin = margin(0, 0, 10, 0)),
          axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
          axis.title.y = element_text(margin = margin(0, 15, 0, 0)),
          axis.text = element_text(color = "black"),
          axis.line = element_line(color = "black")
        )
      
      filename = paste0(output_folder, "/", filename_base, ".png")
      ggsave(filename, plot = histo, width = 8, height = 5, units = "in", dpi = 300)
    } else {
      histo <- ggplot(data.frame(x = c(0, 1)), aes(x)) + 
        geom_blank() +
        labs(title = title, x = "Percentage of Change", y = "Frequency of Pixels") +
        theme_minimal(base_family = "Helvetica") +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(0, 0, 20, 0))) +
        annotate("text", x = 0.5, y = 0.5, label = "No data for the specified conditions", size = 6, hjust = 0.5)
      
      filename = paste0(output_folder, "/", filename_base, ".png")
      ggsave(filename, plot = histo, width = 8, height = 5, units = "in", dpi = 300)
      cat("No valid data to plot for ", title, ", empty plot created.\n")
    }
  }
  
  ## 4.4 Looping for creating the histograms ----
  # ---------------------------------------------------------------- - - -
  for (i in seq_along(masked_deg_acc_rasters)) {
    df_pos <- as.data.frame(masked_deg_acc_rasters[[i]]$pos)
    df_neg <- as.data.frame(masked_deg_acc_rasters[[i]]$neg)
    
    create_histogram(df_pos, paste("Positive EVI -", acc_titles[i]), file.path(OUT_FOLDER, "Histo_AccDeg"), histo_deg_acc_pos_names[i], positive_colors)
    create_histogram(df_neg, paste("Negative EVI -", acc_titles[i]), file.path(OUT_FOLDER, "Histo_AccDeg"), histo_deg_acc_neg_names[i], negative_colors)
  }
  
  for (i in seq_along(masked_deg_tw_rasters)) {
    df_pos <- as.data.frame(masked_deg_tw_rasters[[i]]$pos)
    df_neg <- as.data.frame(masked_deg_tw_rasters[[i]]$neg)
    
    create_histogram(df_pos, paste("Positive EVI -", tw_titles[i]), file.path(OUT_FOLDER, "Histo_TwDeg"), histo_deg_tw_pos_names[i], positive_colors)
    create_histogram(df_neg, paste("Negative EVI -", tw_titles[i]), file.path(OUT_FOLDER, "Histo_TwDeg"), histo_deg_tw_neg_names[i], negative_colors)
  }
  
  cat("\n End of processing for", Site, "\n")
}
