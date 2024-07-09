##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
## 1.1 Loading the libraries and defining ROI ----
# ---------------------------------------------------------
library(terra)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr) # Added for data manipulation
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 ROI ----
# ---------------------------------------------------------
SITES = c("Paragominas") # If you want to loop # 'Cotriguacu', 'Guaviare', 'MDD', 'Paragominas'

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
  # ---------------------------------------------------------
  DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_AObs_deg/DEG_pct")
  
  # Defining and creating output folder
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/Plots/{Site}_Histo_DegPct")
  
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  ## 2.2 Loading the rasters ----
  # ---------------------------------------------------------
  deg_rasters <- list(
    deg_1990_2005 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2005.tif"))),
    deg_1990_2008 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2008.tif"))),
    deg_1990_2011 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2011.tif"))),
    deg_1990_2014 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2014.tif"))),
    deg_1990_2017 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2017.tif"))),
    deg_1990_2020 = rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2020.tif")))
  )
  
  # Renaming the raster values
  for (name in names(deg_rasters)) {
    names(deg_rasters[[name]]) <- "values"
  }
  
  ##%######################################################%##
  #                                                          #
  #                3. Creating histograms                   ----
  #                                                          #
  ##%######################################################%##
  create_histogram <- function(raster, title) {
    values <- values(raster, na.rm = TRUE)
    values_df <- data.frame(values = values)
    
    hist_plot <- ggplot(values_df, aes(x = values)) +
      geom_histogram(binwidth = 1, fill = "#1f78b4", alpha = 0.6, color = "black") +
      labs(title = title,
           x = "% de dégradation",
           y = "Fréquence") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 20)),
        axis.title.x = element_text(size = 14, margin = margin(t = 10)),
        axis.title.y = element_text(size = 14, margin = margin(r = 10)),
        axis.text = element_text(size = 12)
      )
    
    return(hist_plot)
  }
  
  hist_plots <- list()
  hist_plots[[1]] <- create_histogram(deg_rasters$deg_1990_2005, "1990 à 2005")
  hist_plots[[2]] <- create_histogram(deg_rasters$deg_1990_2008, "1990 à 2008")
  hist_plots[[3]] <- create_histogram(deg_rasters$deg_1990_2011, "1990 à 2011")
  hist_plots[[4]] <- create_histogram(deg_rasters$deg_1990_2014, "1990 à 2014")
  hist_plots[[5]] <- create_histogram(deg_rasters$deg_1990_2017, "1990 à 2017")
  hist_plots[[6]] <- create_histogram(deg_rasters$deg_1990_2020, "1990 à 2020")
  
  combined_plot <- arrangeGrob(grobs = hist_plots, ncol = 2, 
                               top = textGrob("Histogrammes du pourcentage de dégradation", gp = gpar(fontsize = 18, fontface = "bold"), vjust = 0),
                               padding = unit(3, "lines"))
  
  ggsave(file.path(OUT_FOLDER, "combined_histograms.png"), combined_plot, width = 13, height = 15, dpi = 300)
  
  cat("\n End of processing for", Site, "\n")
}
