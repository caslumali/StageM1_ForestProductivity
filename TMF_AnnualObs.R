##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
## 1.1 Charging the libraries ----
# ---------------------------------------------------------------- - - -
library(terra)
library(data.table)
library(stringr)
library(zoo)
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 Setting environment ----
# ---------------------------------------------------------------- - - -
SITES = c("Paragominas") # If you want to loop # 'Cotriguacu', 'Guaviare', 'MDD', 'Paragominas'
# Site <- "Paragominas"
ROI_BUFFER <- 2000 # meters if ROI CRS is lat/lon, projection units otherwise
MASK_ROI <- TRUE # if TRUE, will mask using ROI after cropping
WRITE_INTERMEDIATE <- TRUE # for write intermediate data
WRITE_ESSENTIALS <- TRUE # for write essential data
YEARS <- 1990:2020 # naming the layers stacks
SAVE_CLASSIFIED <- TRUE # If TRUE, save the classified raster (0, 25, 50, 75, 100);

##%######################################################%##
#                                                          #
#                    2. Prepare data                       ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  Site <- SITES[i]
  cat("\n Starting processing for", Site, "\n")
  
  ## 2.1 Defining data and output directories ----
  # ---------------------------------------------------------------- - - -
  ROI <- list.files(stringr::str_glue("data/{Site}"), pattern = "shp", full.names = TRUE) # Load ROI shapefile
  
  # Open Annual Change TMF Annual Disruptions
  ANNUAL_OBS <- list.files(stringr::str_glue("data/{Site}/TMF_AD"), pattern = "AnnualDisruptionObs.*\\.tif$", full.names = TRUE) # Load annual disruption files
  
  # Charging EVI data
  # EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")
  EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI_HQ")
  evi_smoothed <- rast(file.path(EVI_FOLDER, "evi_smoothed.tif")) # Load EVI data
  
  # Defining and creating output folders
  # OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/TMF_AnnualObs")
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/TMF_AnnualObs_max_HQ")
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE) # Create output directory
  
  DEG_PCT_FOLDER <- file.path(OUT_FOLDER, "DEG_pct")
  dir.create(DEG_PCT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  DEG_CLASS_FOLDER <- file.path(OUT_FOLDER, "DEG_class")
  dir.create(DEG_CLASS_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  ## 2.2 Charging the ROI ----
  # ---------------------------------------------------------------- - - -
  roi <- if (is.character(ROI)) vect(ROI) else ROI
  
  ##%######################################################%##
  #                                                          #
  #               3. TMF Annual Disruption                   ----
  #                                                          #
  ##%######################################################%##
  ## 3.1 Creates .vrt based on one or more tif files ----
  # ---------------------------------------------------------------- - - -
  if (is.character(ANNUAL_OBS) && length(ANNUAL_OBS) > 0) {
    # Create a single virtual raster tile (VRT) from all .tif files
    VRT_PATH <- paste0(OUT_FOLDER, "/Annual_Disruption.vrt")
    terra::vrt(ANNUAL_OBS, VRT_PATH, overwrite = TRUE)
    message("Single VRT created successfully at: ", VRT_PATH)
  } else {
    message("No .tif files found in the specified directory.")
  }
  
  # Load the created VRT for further use
  if (file.exists(VRT_PATH)) {
    annual_obs <- terra::rast(VRT_PATH)
    message("VRT loaded successfully.")
  } else {
    stop("Failed to find the VRT file.")
  }
  
  # Get only rasters from 1990 to 2020
  annual_obs <- annual_obs[[9:39]]
  
  # Crop TMF Annual Disruption
  roi <- project(roi, annual_obs) # project ROI to match EVI MODIS CRS
  roi_buf <- buffer(roi, ROI_BUFFER) # create buffer
  set.crs(roi_buf, crs(roi)) # ensure buffered ROI has correct CRS
  annual_obs <- crop(annual_obs, roi_buf, mask = MASK_ROI)
  
  # Transform everything that is greater than 0 into 1, keeping 0 and NA as they are
  subs_fun <- function(x) ifelse(is.na(x), NA, ifelse(x == 0, 0, 1))
  deg <- app(annual_obs, fun = subs_fun)
  
  # Saving the masked raster
  if (WRITE_INTERMEDIATE) {
    writeRaster(annual_obs, filename = file.path(OUT_FOLDER, "annual_obs.tif"), overwrite = TRUE, datatype = "FLT4S")
    writeRaster(deg, filename = file.path(OUT_FOLDER, "deg.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  ##%######################################################%##
  #                                                          #
  #           4. Calculate annual degradation            ----
  #                                                          #
  ##%######################################################%##
  cat("\n Calculation degradation for", Site, "\n")
  
  ## 4.1 Creating EVI mask and retrieving the resolution for aggregate and resample ----
  # ---------------------------------------------------------------- - - -
  evi_mask <- evi_smoothed[[1]]
  evi_mask[] <- ifelse(!is.na(evi_mask[]), 1, NA)
  res_modis <- res(evi_mask)
  res_tmf <- res(deg)
  
  # Calculating the aggregation factors (x and y)
  fator_agg_x <- ceiling(res_modis[1] / res_tmf[1])
  fator_agg_y <- ceiling(res_modis[2] / res_tmf[2])
  
  ## 4.2 Function to aggregate and resample the change rasters ----
  # ---------------------------------------------------------------- - - -
  pixel_count <- fator_agg_x * fator_agg_y
  empty_modis <- rast(evi_mask, nlyrs = 1, vals = NA)
  
  aggregate_and_resample <- function(raster, factor_x, factor_y, empty_modis) {
    agg_raster <- aggregate(raster, fact = c(factor_x, factor_y), fun = sum, na.rm = TRUE)
    terra::resample(agg_raster, empty_modis, method = "near")
  }
  
  ## 4.3 Calculate degradation in each layer ----
  # ---------------------------------------------------------------- - - -
  annual_deg_agg <- aggregate_and_resample(deg, fator_agg_x, fator_agg_y, empty_modis)
  annual_deg_pct <- annual_deg_agg / pixel_count * 100
  annual_deg_pct <- mask(annual_deg_pct, evi_mask, maskvalue = NA)
  names(annual_deg_pct) <- YEARS
  
  if (WRITE_INTERMEDIATE) {
    writeRaster(annual_deg_pct, filename = file.path(OUT_FOLDER, "annual_deg_pct.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  
  
  ##%######################################################%##
  #                                                          #
  #           5. Calculate degradation by periods            ----
  #                                                          #
  ##%######################################################%##
  ## 5.1 Defining periods in a list ----
  # ---------------------------------------------------------------- - - -
  periods <- list(
    "1990_1999" = 1:10,  
    "2000_2005" = 11:16, 
    "2003_2008" = 14:19, 
    "2006_2011" = 17:22, 
    "2009_2014" = 20:25, 
    "2012_2017" = 23:28, 
    "2015_2020" = 26:31
  )
  

  ## 5.2 Applying function to take max degradation ----
  # ---------------------------------------------------------------- - - -
  tw_deg_list <- lapply(periods, function(range) {
    app(annual_deg_pct[[range]], fun = max, na.rm = TRUE)
  })
  
  tw_deg <- rast(tw_deg_list)
  
  # Save raster with max degradation by period
  if (WRITE_INTERMEDIATE) {
    writeRaster(tw_deg, filename = file.path(OUT_FOLDER, "tw_deg.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  ## 5.3 Defining matrix to classify the degradation ----
  # ---------------------------------------------------------------- - - -
  rcl <- matrix(c(0, 0, 0,
                  0.1, 25, 25,
                  25, 50, 50,
                  50, 75, 75,
                  75, 100, 100), ncol = 3, byrow = TRUE)
  
  tw_deg_class <- classify(tw_deg, rcl, include.lowest = TRUE)
  
  if (WRITE_ESSENTIALS) {
    if (SAVE_CLASSIFIED) {
      writeRaster(tw_deg_class, filename = file.path(DEG_CLASS_FOLDER, str_glue("{Site}_TW_deg_class.tif")),
                  overwrite = TRUE, datatype = "FLT4S")
    } else {
      writeRaster(tw_deg, filename = file.path(DEG_PCT_FOLDER, str_glue("{Site}_TW_deg.tif")),
                  overwrite = TRUE, datatype = "FLT4S")
    }
  }
  
  if (WRITE_INTERMEDIATE) {
    for (period in names(periods)) {
      if (SAVE_CLASSIFIED) {
        writeRaster(classify(tw_deg_list[[period]], rcl, include.lowest = TRUE), 
                    filename = file.path(DEG_CLASS_FOLDER, str_glue("{Site}_TW_deg_class_{period}.tif")),
                    overwrite = TRUE, datatype = "FLT4S")
      } else {
        writeRaster(tw_deg_list[[period]], 
                    filename = file.path(DEG_PCT_FOLDER, str_glue("{Site}_TW_deg_{period}.tif")),
                    overwrite = TRUE, datatype = "FLT4S")
      }
    }
  }
  
  ## 5.4 Defining accumulated periods in a list ----
  # ---------------------------------------------------------------- - - -
  periods_acc <- list(
    "1990_1999" = 1:10,  
    "1990_2005" = 1:16, 
    "1990_2008" = 1:19, 
    "1990_2011" = 1:22, 
    "1990_2014" = 1:25, 
    "1990_2017" = 1:28, 
    "1990_2020" = 1:31
  )

  ## 5.5 Applying function to take max degradation ----
  # ---------------------------------------------------------------- - - -
  acc_deg_list <- lapply(periods_acc, function(range) {
    app(annual_deg_pct[[range]], fun = max, na.rm = TRUE)
  })
  acc_deg <- rast(acc_deg_list)
  
  # Save raster with max degradation by period
  if (WRITE_INTERMEDIATE) {
    writeRaster(acc_deg, filename = file.path(OUT_FOLDER, "acc_deg.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  # Classify the raster
  acc_deg_class <- classify(acc_deg, rcl, include.lowest = TRUE)
  
  # Saving the rasters
  if (WRITE_ESSENTIALS) {
    if (SAVE_CLASSIFIED) {
      writeRaster(acc_deg_class, filename = file.path(DEG_CLASS_FOLDER, str_glue("{Site}_ACC_deg_class.tif")),
                  overwrite = TRUE, datatype = "FLT4S")
    } else {
      writeRaster(acc_deg, filename = file.path(DEG_PCT_FOLDER, str_glue("{Site}_ACC_deg.tif")),
                  overwrite = TRUE, datatype = "FLT4S")
    }
  }
  
  if (WRITE_INTERMEDIATE) {
    for (period in names(periods_acc)) {
      if (SAVE_CLASSIFIED) {
        writeRaster(classify(acc_deg_list[[period]], rcl, include.lowest = TRUE), 
                    filename = file.path(DEG_CLASS_FOLDER, str_glue("{Site}_ACC_deg_class_{period}.tif")),
                    overwrite = TRUE, datatype = "FLT4S")
      } else {
        writeRaster(acc_deg_list[[period]], 
                    filename = file.path(DEG_PCT_FOLDER, str_glue("{Site}_ACC_deg_{period}.tif")),
                    overwrite = TRUE, datatype = "FLT4S")
      }
    }
  }
  
  cat("\n End of processing for", Site, "\n")
}

