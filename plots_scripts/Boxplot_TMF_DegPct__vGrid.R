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
library(gridExtra)
library(grid)
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 ROI ----
# ---------------------------------------------------------------- - - -
SITES = c("Paragominas") # If you want to loop # 'Cotriguacu', 'Guaviare', 'MDD', 'Paragominas'

##%######################################################%##
#                                                          #
#                    2. Charging data                      ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  # Site
  Site = SITES[i]
  cat("\n Starting processing for", Site, "\n")
  
  ## 2.1 Charging data and output directories ----
  # ---------------------------------------------------------------- - - -   
  EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")
  DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_AObs_deg/DEG_pct")
  
  # Defining and creating output folder
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/Plots/{Site}_Bplots_Memoire")
  
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUT_FOLDER, "Bplot_AccDeg"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUT_FOLDER, "Bplot_TwDeg"), showWarnings = FALSE, recursive = TRUE)
  
  ## 2.2 Charging the rasters ----
  # ---------------------------------------------------------------- - - -
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
  ## 3.1 Function to creating evi mask 
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
  
  # Iterating through rasters and applying masks
  masked_deg_acc_rasters <- mapply(apply_masks, evi_rasters, deg_acc_rasters, SIMPLIFY = FALSE)
  masked_deg_tw_rasters <- mapply(apply_masks, evi_rasters, deg_tw_rasters, SIMPLIFY = FALSE)
  
  ##%######################################################%##
  #                                                          #
  #                4. Creating boxplots                      ----
  #                                                          #
  ##%######################################################%##
  ## 4.1 Preparing data for boxplots
  # ---------------------------------------------------------------- - - -
  prepare_boxplot_data <- function(masked_list, period) {
    pos_df <- as.data.frame(masked_list$pos)
    pos_df$category <- "Positive"
    pos_df$period <- period
    
    neg_df <- as.data.frame(masked_list$neg)
    neg_df$category <- "Negative"
    neg_df$period <- period
    
    ns_df <- as.data.frame(masked_list$ns)
    ns_df$category <- "Not Significant"
    ns_df$period <- period
    
    rbind(pos_df, neg_df, ns_df)
  }
  
  periods_accum <- c("1990-2005", "1990-2008", "1990-2011", "1990-2014", "1990-2017", "1990-2020")
  periods_tw <- c("2000-2005", "2003-2008", "2006-2011", "2009-2014", "2012-2017", "2015-2020")
  
  boxplot_deg_acc_data <- mapply(prepare_boxplot_data, masked_deg_acc_rasters, periods_accum, SIMPLIFY = FALSE)
  boxplot_deg_tw_data <- mapply(prepare_boxplot_data, masked_deg_tw_rasters, periods_tw, SIMPLIFY = FALSE)
  
  ## 4.2 Creating and saving individual boxplots ----
  # ---------------------------------------------------------------- - - -
  create_boxplot <- function(data, period) {
    data$category <- factor(data$category, levels = c("Positive", "Negative", "Not Significant"))
    
    ggplot(data, aes(x = category, y = values, fill = category)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.6) +
      scale_fill_manual(values = c("Positive" = "#1a9850", "Negative" = "#d73027", "Not Significant" = "#fee08b")) +
      labs(title = period,
           x = "Tendance de l'EVI",
           y = "% de dégradation") +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold", size = 16, margin = margin(0, 0, 5, 0)),
        axis.title.x = element_text(size = 14, margin = margin(5, 0, 0, 0)),
        axis.title.y = element_text(size = 14, margin = margin(0, 5, 0, 0)),
        axis.text = element_text(color = "black", size = 12),
        axis.line = element_line(color = "black"),
        legend.position = "none",
        axis.text.x = element_blank()  # Removendo texto das categorias
      )
  }
  
  ## 4.3 Creating and combining plots ----
  # ---------------------------------------------------------------- - - -
  create_combined_plot <- function(data_list, periods, output_file) {
    plots <- lapply(seq_along(data_list), function(i) {
      create_boxplot(data_list[[i]], periods[i])
    })
    
    combined_plot <- arrangeGrob(grobs = plots, ncol = 2,
                                 top = textGrob("Pourcentage de dégradation selon la tendance de l'EVI", gp = gpar(fontsize = 18, fontface = "bold"), vjust = 0),
                                 padding = unit(1, "lines"))
    
    # combined_plot <- arrangeGrob(grobs = plots, ncol = 2)
                                 
    
    
    # Custom legend
    legend_data <- data.frame(category = c("Positive", "Negative", "Non significatif"), x = 1:3, y = 1)
    
    legend_plot <- ggplot(legend_data, aes(x = factor(x), y = y, fill = category)) +
      geom_boxplot(aes(ymin = y - 0.1, lower = y - 0.05, middle = y, upper = y + 0.05, ymax = y + 0.1, alpha = 0.6), stat = "identity", width = 0.2) +
      scale_fill_manual(values = c("Positive" = "#1a9850", "Negative" = "#d73027", "Non significatif" = "#fee08b")) +
      theme_void() +
      theme(legend.position = "none") +
      geom_text(aes(label = ifelse(category == "Positive", category, "")), hjust = -0.75, vjust = 0.5, size = 5) +
      geom_text(aes(label = ifelse(category == "Negative", category, "")), hjust = -0.7, vjust = 0.5, size = 5) +
      geom_text(aes(label = ifelse(category == "Non significatif", category, "")), hjust = -0.4, vjust = 0.5, size = 5) +
      scale_x_discrete(limits = c("1", "2", "3")) +
      scale_y_continuous(expand = c(0, 0))
    
    legend_grob <- ggplotGrob(legend_plot)
    
    # Add blank space between the plots and the legend
    blank_grob <- grid::rectGrob(gp = grid::gpar(col = NA, fill = NA))
    
    # Adjust the final plot to include more space after the legend
    final_plot <- arrangeGrob(combined_plot, blank_grob, legend_grob, ncol = 1, 
                              heights = unit.c(unit(1, "npc") - unit(4, "lines"), unit(1, "lines"), unit(2, "lines")))
    
    ggsave(output_file, final_plot, width = 10, height = 12, dpi = 300)
  }
  
  
  create_combined_plot(boxplot_deg_acc_data, periods_accum, file.path(OUT_FOLDER, "Bplot_AccDeg_combined.png"))
  create_combined_plot(boxplot_deg_tw_data, periods_tw, file.path(OUT_FOLDER, "Bplot_TwDeg_combined.png"))
  
  cat("\n End of processing for", Site, "\n")
}
