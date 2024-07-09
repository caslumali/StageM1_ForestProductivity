##%######################################################%##
#                                                          #
#                     1. Settings                          ----
#                                                          #
##%######################################################%##
## 1.1 Charging the libraries ----
# ---------------------------------------------------------------- - - -
library(terra)
library(stringr)
library(bitops)
library(signal)
library(zoo)
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 Setting environment ----
# ---------------------------------------------------------------- - - -
# SITES = c("Paragominas") # If you want to loop # 'Cotriguacu', 'Guaviare', 'MDD', 'Paragominas'
Site <- "Paragominas"

ROI_BUFFER <- 2000 # meters if ROI CRS is lat/lon, projection units otherwise
MASK_ROI <- TRUE # if TRUE, will mask using ROI after cropping
WRITE_INTERMEDIATE <- TRUE # for write intermediate data
WRITE_ESSENTIALS <- TRUE #for write essential data


## 1.3 Setting TW parameters ----
# ---------------------------------------------------------------- - - -
# Uncomment/comment the line below to analyze EVI trends in different periods
# TW <- list("2000_2010" = 1:230, "2011_2020" = 231:460) # to analyze two periods
# TW <- list("2000_2013" = 1:299, "2014_2020" = 300:460) # to analyze two periods

# Define temporal windows for EVI
# For looping to create 5-year time windows rolling every 3
TW <- list()
tw_years <- seq(2000, 2015, by=3)
for (i in 1:length(tw_years)) {
  start_year <- tw_years[i]
  start_index <- 1 + 23 * (start_year - 2000)  # Calculate the initial index
  end_index <- start_index + 115 - 1           # Calculate the final index
  TW_name <- paste(start_year, start_year + 5, sep="_")  # Name of the time window
  TW[[TW_name]] <- start_index:end_index
}

##%######################################################%##
#                                                          #
#                    2. Prepare data                       ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  Site = SITES[i]
  cat("\n Starting processing MODIS-EVI for", Site, "\n")
  
  ## 2.1 Defining folder data and output directories ----
  # ---------------------------------------------------------------- - - -   
  ROI <- list.files(stringr::str_glue("data/{Site}"), pattern = "shp", full.names = TRUE) # filename or SpatVector
  ANNUAL_CHANGE <- list.dirs(stringr::str_glue("data/{Site}/TMF_AC"), recursive = FALSE, full.names = TRUE) # filename(s) or SpatRaster
  EVI_BIMONTHLY <- stringr::str_glue("data/{Site}/MODIS_NASA/MODIS_EVI")
  MODIS_QA <- stringr::str_glue("data/{Site}/MODIS_NASA/MODIS_QA") # Pixel Reliability
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/EVI_HQ")
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  OLS_FOLDER <- file.path(OUT_FOLDER, "OLS_outputs")
  dir.create(OLS_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
  
  ## 2.2 Charging the ROI ----
  # ---------------------------------------------------------------- - - -
  if (is.character(ROI)) {
    roi <- vect(ROI)
  } else {
    roi <- ROI
  }
  
  ##%######################################################%##
  #                                                          #
  #               3. TMF data preparation                    ----
  #                                                          #
  ##%######################################################%##
  ## 3.1 Creates .vrt based on one (T1.vrt) or more tif files (Annual_change.vrt) ----
  # ---------------------------------------------------------------- - - -
  if (is.character(ANNUAL_CHANGE) && length(ANNUAL_CHANGE) > 0) {
    vrt_paths <- character(length(ANNUAL_CHANGE))
    for (i in seq_along(ANNUAL_CHANGE)) {
      tif_files <- list.files(path = ANNUAL_CHANGE[i], pattern = "JRC_TMF_AnnualChange_v1_2020_.*\\.tif$", full.names = TRUE, recursive = TRUE)
      if (length(tif_files) > 0) {
        vrt_path <- paste0(ANNUAL_CHANGE[i], "/T", i, ".vrt")
        terra::vrt(tif_files, vrt_path, options = "-separate", overwrite = TRUE)
        vrt_paths[i] <- vrt_path
      }
    }
    if (any(nzchar(vrt_paths))) {
      terra::vrt(vrt_paths, paste0(OUT_FOLDER, "/Annual_change_2020.vrt"), overwrite = TRUE)
    } else {
      message("No .tif files found in the specified directories.")
    }
  } else {
    message("No annual change data specified or the path is not character.")
  }
  
  ANNUAL_CHANGE <- paste0(OUT_FOLDER, "/Annual_change_2020.vrt")
  if (file.exists(ANNUAL_CHANGE)) {
    annual_change <- terra::vrt(ANNUAL_CHANGE)
  } else {
    message("Final VRT was not created or not found.")
  }
  
  ANNUAL_CHANGE <- list.files(OUT_FOLDER, pattern = "Annual_change_2020.vrt", full.names = TRUE) # filename(s) or SpatRaster
  annual_change <- vrt(ANNUAL_CHANGE)
  
  # Crop annual change data
  roi <- project(roi, annual_change) # project ROI to match annual change CRS
  roi_buf <- buffer(roi, ROI_BUFFER) # create buffer
  set.crs(roi_buf, crs(roi)) # ensure buffered ROI has correct CRS
  annual_change <- crop(annual_change, roi_buf, mask = MASK_ROI) 
  
  # Reclassify Annual Change to create forest cover mask
  fc_2020 <- subst(annual_change, c(NA, 1, 2), c(NA, 1, 1), others = 0)
  
  # Saving the forest cover raster
  if (WRITE_INTERMEDIATE) {
    writeRaster(fc_2020, filename = file.path(OUT_FOLDER, "fc_2020.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  ##%######################################################%##
  #                                                          #
  #               4. Charging MODIS data                     ----
  #                                                          #
  ##%######################################################%##
  ## 4.1 Create EVI .vrt based tif files ----
  # ---------------------------------------------------------------- - - -
  if (is.character(EVI_BIMONTHLY)) {
    if (length(EVI_BIMONTHLY) == 1) {
      terra::vrt(list.files(path = EVI_BIMONTHLY, pattern = ".tif$", full.names = TRUE, recursive = TRUE), paste0(OUT_FOLDER, "/EVI_bimonthly.vrt"), options = "-separate", overwrite = TRUE)
    }  
  } else {
    evi_bimonthly <- EVI_BIMONTHLY
  }
  
  EVI_BIMONTHLY <- list.files(OUT_FOLDER, pattern = "EVI_bimonthly.vrt", full.names = TRUE) # filename(s) or SpatRaster
  evi_bimonthly <- vrt(EVI_BIMONTHLY)
  
  
  ## 4.2 Create Pixel Realiability .vrt based tif files ----
  # ---------------------------------------------------------------- - - -
  if (is.character(MODIS_QA)) {
    if (length(MODIS_QA) == 1) {
      terra::vrt(list.files(path = MODIS_QA, pattern = ".tif$", full.names = TRUE, recursive = TRUE), paste0(OUT_FOLDER, "/MODIS_QualityAssessment.vrt"), options = "-separate", overwrite = TRUE)
    }  
  } else {
    modis_qa <- MODIS_QA
  }
  
  MODIS_QA <- list.files(OUT_FOLDER, pattern = "MODIS_QualityAssessment.vrt", full.names = TRUE) # filename(s) or SpatRaster
  modis_qa <- vrt(MODIS_QA)
  
  ##%######################################################%##
  #                                                          #
  #                 5. Cleaning MODIS data                   ----
  #                                                          #
  ##%######################################################%##
  cat("\n Cleaning MODIS-EVI data for", Site, "\n")
  
  
  ## 5.1 Filtering MODIS with Quality Assessment data ----
  # ---------------------------------------------------------------- - - -
  QA_mask_fun <- function(x) {
    vi_quality <- bitwAnd(x, 0x03) %in% c(0, 1, 2) # Apply mask to VI quality bits (bits 0-1)
    vi_usefulness <- bitwAnd(x %/% 4, 0x0F) %in% c(0, 1, 2, 4, 8, 9, 10, 12) # Apply mask to VI usefulness bits (bits 2-5)
    return(vi_quality & vi_usefulness) # Returns the logical combination of the two masks
  }
  
  qa_mask <- app(modis_qa, QA_mask_fun)
  
  ## 5.2 Masking MODIS with Quality Assessment data ----
  # ---------------------------------------------------------------- - - -
  # Process paragominas before to solve missing raster problem
  if (Site == "Paragominas") {
    evi_useful <- mask(evi_bimonthly, qa_mask, maskvalue = 0) # Apply mask to EVI raster
    
    index1 <- 116 # Image from 2005_225
    index2 <- 117 # Image from 2005_257
    
    mean_fun <- function(x) { # Function to calculate mean between the two rasters handling NAs
      if (is.na(x[1]) && !is.na(x[2])) return(x[2])
      else if (!is.na(x[1]) && is.na(x[2])) return(x[1])
      else if (!is.na(x[1]) && !is.na(x[2])) return((x[1] + x[2]) / 2)
      else return(NA)
    }
    
    raster_mean <- app(c(evi_useful[[index1]], evi_useful[[index2]]), fun = mean_fun) # Calculate the mean
    paragominas_pt1 <- evi_useful[[1:(index2-1)]] # Take the fist stack
    paragominas_pt2 <- evi_useful[[index2:nlyr(evi_useful)]] # Take the second stack
    evi_useful <- c(paragominas_pt1, raster_mean, paragominas_pt2) # Put the mean raster in the middle
    
  } else {
    evi_useful <- mask(evi_bimonthly, qa_mask, maskvalue = 0) # Apply mask to EVI raster
  }
  
  
  ## 5.3 EVI data cleaning and interpolation ----
  # ---------------------------------------------------------------- - - -
  max_consecutive_na <- function(x) { # Function to calculate the maximum consecutive NAs
    rle_na <- rle(is.na(x))
    max(rle_na$lengths[rle_na$values])
  }
  
  # Create mask with a pixel has more than 12 consecutive NA
  mask_consecutive_na <- app(evi_useful, fun = max_consecutive_na)
  discard_mask <- mask_consecutive_na < 12
  evi_cleaned <- mask(evi_useful, discard_mask, maskvalue = FALSE) # Apply mask and discard inadequate pixels
  
  # Crop and interpolation of EVI useful data
  roi <- project(roi, evi_cleaned) # project ROI to match EVI MODIS CRS
  roi_buf <- buffer(roi, ROI_BUFFER) # create buffer
  set.crs(roi_buf, crs(roi)) # ensure buffered ROI has correct CRS
  evi_cleaned <- crop(evi_cleaned, roi_buf, mask = FALSE)
  evi_cleaned <- crop(evi_cleaned, roi_buf, mask = MASK_ROI)
  evi_filled <- approximate(evi_cleaned, method = "linear") # Linear interpolation of remaining NA values
  
  if (WRITE_INTERMEDIATE) { # Save the masked raster
    writeRaster(evi_filled, filename = file.path(OUT_FOLDER, "evi_filled.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  ## 5.4 Masking EVI data based on 2020 forest cover ----
  # ---------------------------------------------------------------- - - -
  res_evi <- res(evi_filled) # Retrieving the resolution of both rasters in degrees
  res_fc <- res(fc_2020)
  fator_agg_x <- ceiling(res_evi[1] / res_fc[1]) # Calculating the aggregation factors (x and y)
  fator_agg_y <- ceiling(res_evi[2] / res_fc[2])
  
  fc_2020_agg <- aggregate(fc_2020, fact = c(fator_agg_x, fator_agg_y), fun = sum, na.rm = TRUE) # Aggregating using the computed factors
  num_cells <- fator_agg_x * fator_agg_y
  fc_perc <- fc_2020_agg / num_cells # Calculating the correct percentage of forest cover
  fc_mask <- ifel(fc_perc >= 0.95, 1, 0) # creating a mask using a 95% threshold
  
  fc_mask_rsp <- terra::resample(fc_mask, evi_filled, method = "near") # Resampling and applying the mask
  evi_fc_masked <- mask(evi_filled, fc_mask_rsp, maskvalue = 0)
  
  if (WRITE_INTERMEDIATE) { # Saving the masked raster
    writeRaster(evi_fc_masked, filename = file.path(OUT_FOLDER, "evi_fc_masked.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  
  ## 5.5 Smoothing NDVI data ----
  # ---------------------------------------------------------------- - - -
  EVI_smoothed <- function(raster, window_size, poly_order) { # Function for smoothing the EVI with Savitzky-Golay, maintaining the NA
    smoothing_fun <- function(v) {
      if (all(is.na(v))) {
        return(rep(NA, length(v)))
      } # Check if all values are NA
      v <- na.spline(v, na.rm = FALSE) # Spline interpolation for NA
      padded <- c(rev(v[1:window_size]), v, rev(v[(length(v) - window_size + 1):length(v)])) # Apply mirroring to data to fix the border problem 
      smoothed_padded <- sgolayfilt(padded, p = poly_order, n = window_size)
      smoothed <- smoothed_padded[(window_size + 1):(length(smoothed_padded) - window_size)] # Remove padding
      return(smoothed)
    }
    evi_smoothed <- app(raster, smoothing_fun) # Applies the smoothing function to the raster
    return(evi_smoothed)
  }
  
  evi_smoothed <- EVI_smoothed(evi_fc_masked, window_size = 9, poly_order = 3) # Applies the smoothing function
  
  if (WRITE_ESSENTIALS) { # Saving the raster
    writeRaster(evi_smoothed, filename = file.path(OUT_FOLDER, "evi_smoothed.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  
  evi_smoothed <- rast(file.path(OUT_FOLDER, "evi_smoothed.tif"))
  
  ##%######################################################%##
  #                                                          #
  #               6. Calculating EVI trends                  ----
  #                                                          #
  ##%######################################################%##
  ## 6.1 Analyze Trend function ----
  # ---------------------------------------------------------------- - - -
  EVI_trend <- function(evi_data, time_window = NULL) {
    if (is.null(time_window)) { # If 'time_window' is not provided, use the entire dataset
      time_window <- 1:nlyr(evi_data)
    }
    time <- 1:length(time_window) # Setting the time vector for the specific window
    evi_subset <- subset(evi_data, time_window)
    
    slope_fun <- function(x) {
      if (all(is.na(x))) return(NA_real_)
      lm(x ~ time)$coefficients[2] # Extracting the slope coefficient
    }
    
    pvalue_fun <- function(x) {
      if (all(is.na(x))) return(NA_real_)
      summary(lm(x ~ time))$coefficients["time", "Pr(>|t|)"] # P-value for 'time'
    }
    
    slope <- app(evi_subset, slope_fun)
    pvalue <- app(evi_subset, pvalue_fun)
    
    return(list(slope = slope, pvalue = pvalue))
  }
  
  
  ## 6.2 Classify Trend function ----
  # ---------------------------------------------------------------- - - -
  EVI_classify <- function(pvalue_raster, slope_raster, strict_pvalue = 0.001) {
    class_raster <- pvalue_raster # Create a new raster for classifications based on pvalue_raster
    values(class_raster) <- 0 # Initialize with 0 (not significant)
    
    pvalues <- values(pvalue_raster) # Get raster values for processing
    slopes <- values(slope_raster)
    
    values(class_raster)[slopes > 0 & pvalues <= strict_pvalue] <- 1  # Positive and significant
    values(class_raster)[slopes < 0 & pvalues <= strict_pvalue] <- -1 # Negative and significant
    is_na <- is.na(pvalues)
    values(class_raster)[is_na] <- NA # Keep NA where pvalue is NA
    
    return(class_raster)
  }
  
  
  ## 6.3 Calculating Transition Matrix ----
  # ---------------------------------------------------------------- - - -
  calculate_transition <- function(raster1, raster2, output_folder, transition_name) {
    values1 <- values(raster1) # Get the values of the rasters
    values2 <- values(raster2)
    
    transition_class <- rep(NA, length(values1)) # Create a vector to store transitions
    
    # Classify transitions
    transition_class[values1 == -1 & values2 == -1] <- -4 # from negative to negative
    transition_class[values1 == 1 & values2 == -1] <- -3  # from positive to negative
    transition_class[values1 == -1 & values2 == 0] <- -2  # from negative to not significant
    transition_class[values1 == 0 & values2 == -1] <- -1  # from not significant to negative
    transition_class[values1 == 0 & values2 == 0] <- 0    # from not significant to not significant
    transition_class[values1 == 0 & values2 == 1] <- 1    # from not significant to positive
    transition_class[values1 == 1 & values2 == 0] <- 2    # from positive to not significant
    transition_class[values1 == -1 & values2 == 1] <- 3   # from negative to positive
    transition_class[values1 == 1 & values2 == 1] <- 4    # from positive to positive
    
    transition_raster <- raster1 # Create a transition raster with the calculated values
    values(transition_raster) <- transition_class
    
    transition_table <- table(transition_class, useNA = "ifany") # Create the transition table
    
    if (WRITE_INTERMEDIATE) { # Save the transition matrix as CSV
      write.csv(as.data.frame(transition_table), file.path(output_folder, paste0("transition_", transition_name, ".csv")))
    }
    
    if (WRITE_ESSENTIALS) { # Save the transition raster
      writeRaster(transition_raster, filename = file.path(output_folder, str_glue("{Site}_{transition_name}.tif")), overwrite = TRUE, datatype = "FLT4S")
    }
    
    return(list(transition_matrix = transition_table, transition_raster = transition_raster))
  }
  
  ## 6.4 Applying trend, classify, and transition functions within existing loop ----
  # ---------------------------------------------------------------- - - -
  classified_rasters <- list()
  transition_rasters <- list()
  
  if (exists("TW") && length(TW) > 0) {
    for (i in 1:length(TW)) {
      tw_name <- names(TW)[i]
      cat("\n Calculating EVI trends for", tw_name, "in", Site, "\n")
      time_window <- TW[[tw_name]]
      results_evi_tw <- EVI_trend(evi_smoothed, time_window)
      
      # Saving the trend results
      if (WRITE_INTERMEDIATE) {
        writeRaster(results_evi_tw$slope, file.path(OLS_FOLDER, str_glue("evi_slope_{tw_name}.tif")),
                    overwrite = TRUE, datatype="FLT4S")
        writeRaster(results_evi_tw$pvalue, file.path(OLS_FOLDER, str_glue("evi_pvalue_{tw_name}.tif")),
                    overwrite = TRUE, datatype="FLT4S")
      }
      
      classified_raster_tw <- EVI_classify(results_evi_tw$pvalue, results_evi_tw$slope) # Classify the rasters
      classified_rasters[[tw_name]] <- classified_raster_tw
      if (WRITE_ESSENTIALS) {
        writeRaster(classified_raster_tw, file.path(OUT_FOLDER, str_glue("{Site}_evi_class_{tw_name}.tif")), overwrite = TRUE, datatype = "FLT4S")
      }
    }
    cat("\n Creating transition matrix and raster for", Site, "\n")
    for (j in 1:(length(classified_rasters) - 1)) { # Generate transition rasters and count significant transitions
      transition_name <- paste(names(TW)[j], names(TW)[j + 1], sep = "_to_")
      transition_results <- calculate_transition(classified_rasters[[j]], classified_rasters[[j + 1]], OUT_FOLDER, transition_name)
      transition_rasters[[transition_name]] <- transition_results$transition_raster
    }
  } else {
    cat("\n Calculating EVI trends for the entire dataset in", Site, "\n")
    results_evi <- EVI_trend(evi_smoothed)
    
    # Saving the trend results
    if (WRITE_INTERMEDIATE) {
      writeRaster(results_evi$slope, file.path(OLS_FOLDER, str_glue("evi_slope.tif")),
                  overwrite = TRUE, datatype="FLT4S")
      writeRaster(results_evi$pvalue, file.path(OLS_FOLDER, str_glue("evi_pvalue.tif")),
                  overwrite = TRUE, datatype="FLT4S")
    }
    
    classified_raster <- EVI_classify(results_evi$pvalue, results_evi$slope)
    if (WRITE_ESSENTIALS) {
      writeRaster(classified_raster, file.path(OUT_FOLDER, str_glue("{Site}_evi_class_all.tif")), overwrite = TRUE, datatype = "FLT4S")
    }
  }
  cat("\n End of processing for", Site, "\n")
}
