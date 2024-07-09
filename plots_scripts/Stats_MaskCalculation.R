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
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 Setting environment ----
# ---------------------------------------------------------------- - - -
# SITES = c("Paragominas") # If you want to loop # 'Cotriguacu', 'Guaviare', 'MDD', 'Paragominas'
Site <- "Paragominas"


ROI_BUFFER <- 2000 # meters if ROI CRS is lat/lon, projection units otherwise
MASK_ROI <- TRUE # if TRUE, will mask using ROI after cropping
WRITE_INTERMEDIATE <- TRUE # for write intermediate data
WRITE_ESSENTIALS <- TRUE # for write essential data

##%######################################################%##
#                                                          #
#                    2. Prepare data                       ----
#                                                          #
##%######################################################%##
for (i in 1:length(SITES)) {
  Site
  # Site = SITES[i]
  cat("\nStarting processing MODIS-EVI for", Site, "\n")
  
  ## 2.1 Defining folder data and output directories ----
  # ---------------------------------------------------------------- - - -   
  ROI <- list.files(stringr::str_glue("data/{Site}"), pattern = "shp", full.names = TRUE) # filename or SpatVector
  ANNUAL_CHANGE <- list.dirs(stringr::str_glue("data/{Site}/TMF_AC"), recursive = FALSE, full.names = TRUE) # filename(s) or SpatRaster
  EVI_BIMONTHLY <- stringr::str_glue("data/{Site}/MODIS_NASA/MODIS_EVI")
  MODIS_QA <- stringr::str_glue("data/{Site}/MODIS_NASA/MODIS_QA") # Pixel Reliability
  
  OUT_FOLDER <- stringr::str_glue("results/{Site}/Productivity/Plots/Mask_stats")
  dir.create(OUT_FOLDER, showWarnings = FALSE, recursive = TRUE)
  
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
  annual_change <- crop(annual_change, roi, mask = MASK_ROI) 
  
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
  
  ## 5.4 Crop filtre for stats calculation ----
  # ---------------------------------------------------------------- - - -
  roi <- project(roi, discard_mask) # project ROI to match EVI MODIS CRS
  useful_px_mask <- crop(discard_mask, roi, mask = MASK_ROI)
  
  # Saving the masks for verification
  if (WRITE_INTERMEDIATE) {
    writeRaster(useful_px_mask, filename = file.path(OUT_FOLDER, "useful_px_mask.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
 
  ##%######################################################%##
  #                                                          #
  #                 6. Calculate Statistics                  ----
  #                                                          #
  ##%######################################################%##
  
  # Function to calculate the area in square meters of a single pixel at a given latitude
  area_of_pixel <- function(lat, res) {
    earth_radius <- 6378137 # Earth's radius in meters
    lat_rad <- lat * pi / 180 # Convert latitude to radians
    meter_per_deg <- pi * earth_radius * cos(lat_rad) / 180 # Meters per degree longitude
    area_m2 <- (meter_per_deg * res[1]) * (111320 * res[2]) # 111320 meters per degree latitude
    return(area_m2)
  }
  
  ## Function to calculate pixel statistics and area ----
  # ---------------------------------------------------------------- - - -
  calculate_stats <- function(raster, output_file) {
    freq_values <- freq(raster) # Frequency of values
    total_pixels <- sum(freq_values$count) # Total number of pixels
    false_pixels <- freq_values$count[freq_values$value == 0] # Total number of FALSE (0) pixels
    percent_false <- (false_pixels / total_pixels) * 100 # Percentage of FALSE pixels
    
    # Calculate the latitude at the center of the raster for area calculation
    center_lat <- (ymax(raster) + ymin(raster)) / 2
    resolution <- res(raster)
    pixel_area_m2 <- area_of_pixel(center_lat, resolution) # Area of a single pixel in square meters
    
    false_area_ha <- (false_pixels * pixel_area_m2) / 10000 # Convert to hectares
    total_area_ha <- (total_pixels * pixel_area_m2) / 10000 # Convert to hectares
    
    # Create a data frame to store the results
    results <- data.frame(
      Total_Pixels = total_pixels,
      False_Pixels = false_pixels,
      Percent_False = percent_false,
      False_Area_Ha = false_area_ha,
      Total_Area_Ha = total_area_ha
    )
    
    # Print results to console
    print(results)
    
    # Write results to CSV
    write.csv(results, file = output_file, row.names = FALSE)
  }
  
  # Calculate stats for useful_px_mask
  calculate_stats(useful_px_mask, file.path(OUT_FOLDER, "useful_px_mask_stats.csv"))
  
  # Calculate stats for fc_2020
  calculate_stats(fc_2020, file.path(OUT_FOLDER, "fc_2020_stats.csv"))
  
  
  cat("\n End of processing for", Site, "\n")

}

  