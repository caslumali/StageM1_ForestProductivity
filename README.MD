### README for Forest Productivity and Degradation Monitoring Project

## Dissertation

The full dissertation titled "Productivité végétale et dégradation forestière : Étude à Paragominas dans le cadre du projet TerrAmaz" is available in this repository. You can view or download it [here](./UE801_Rapport_LucasLima.pdf).


---

#### Overview

This repository contains the complete set of scripts and tools used in the master's thesis project titled **"Forest Productivity and Degradation Monitoring: A Case Study in Paragominas within the TerrAmaz Project"** by **Lucas Pereira das Neves Souza Lima**. The project was conducted as part of the "Master 1 Géomatique" program at the Université de Toulouse, under the supervision of CIRAD (Centre de coopération Internationale en Recherche Agronomique pour le Développement).

The study aimed to evaluate forest productivity in Amazonian territories using Enhanced Vegetation Index (EVI) data derived from MODIS satellite imagery. The project also compared these productivity trends with forest degradation data provided by the Joint Research Centre (JRC) of the European Commission.

This repository is structured into several directories and scripts, each responsible for different aspects of data processing, analysis, and visualization.

---

### Repository Structure

#### 1. **`data_analyses` Directory**

The data_analyses directory contains the primary scripts used throughout the research project. These scripts handle data processing, analysis, and the main statistical evaluations.

- **`EVI_Trend.R`**: Script to compute EVI trends over time using a linear regression model on MODIS data.
- **`TMF_AnnualObs.R`**: Handles the processing of Tropical Moist Forest (TMF) data to assess annual forest disruptions.
- **`ShinyApp_EVI_AccDeg.R`**: A Shiny application designed for interactive exploration of EVI trends and forest degradation data.

#### 2. **`plot_scripts` Directory**

This directory contains scripts specifically designed for generating plots and calculating key statistics related to the research.

- **`Boxplot_TMF_DegPct__vGrid.R`**: Produces boxplots comparing the percentage of forest degradation against EVI trends.
- **`StackedBar_EVI_Deg_Calculation_km2.R`**: Creates stacked bar plots showing the area of different EVI trends and degradation classes over time.
- **`PlotLine_EVI_Deg_area_km2.R`**: Script to plot the area of EVI trends and degraded areas over time.
- **`Histo_TMF_Deg_class.R`**: Generates histograms to visualize the distribution of forest degradation classes.
- **`Plot_EVI_Savitzky.R`**: Applies the Savitzky-Golay filter to smooth the EVI data and plots the results.
- **`Stats_MaskCalculation.R`**: Computes various statistics related to forest cover and degradation based on masking techniques.

#### 3. **`GEE_exportation` Directory**

This directory contains Jupyter Notebooks used for interacting with Google Earth Engine (GEE) to download and visualize relevant satellite data.

- **`GEE_visualisation_EVI.ipynb`**: Notebook for visualizing EVI data in GEE.
- **`GEE_export_Mapbiomas_Amazonia.ipynb`**: Exports MapBiomas data specific to the Amazon region from GEE.
- **`GEE_export_NDVI-EVI.ipynb`**: Notebook to download NDVI and EVI time series data from GEE.
- **`GEE_TMF_AnnualDisruptionsObs.ipynb`**: Exports TMF data related to annual forest disruptions observed in the Amazon.

---

### Installation and Requirements

To run the scripts provided in this repository, you will need the following software and packages:

- **R** (version 4.3.3 or higher)
  - Packages: `bitops`, `data.table`, `dplyr`, `ggplot2`, `leaflet`, `plotly`, `shiny`, `signal`, `stringr`, `terra`, `tidyr`, `zoo`
- **QGIS** (version 3.34.4 or higher) for map visualization.
- **Python** (version 3.12 or higher) with `ee` and `geemap` packages for Google Earth Engine integration.

**Jupyter Notebooks**: Ensure that you have `jupyter` installed to run the `.ipynb` files in the `GEE_exportation` directory.

---

### Usage

1. **Data Processing**: Start with the scripts in the root directory to process the MODIS EVI and TMF data.
2. **Plotting**: Use the scripts in the `plot_scripts` directory to generate visualizations and key statistics.
3. **Google Earth Engine**: Leverage the notebooks in the `GEE_exportation` directory to download and visualize additional satellite data.

---

### License

This project is licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license. This means you are free to share, copy, and redistribute the material in any medium or format and adapt, remix, transform, and build upon the material for any purpose, even commercially, as long as proper credit is given to the author.

---

### Contact

For any questions or further information, please contact **Lucas Lima** at [caslumali@gmail.com].
