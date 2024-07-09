##%######################################################%##
#                                                          #
#                 1. EVI Trend Analysis                    ----
#                                                          #
##%######################################################%##
## 1.1 Charging the libraries
# ---------------------------------------------------------------- - - -
library(shiny)
library(terra)
library(plotly)
library(stringr)
library(leaflet)
terraOptions(tempdir="D:/temp")
terraOptions()

## 1.2 Settings ----
# ---------------------------------------------------------------- - - -
Site <- "Paragominas"
# EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")
EVI_FOLDER <- str_glue("results/{Site}/Productivity/EVI")

# DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_deg/DEG_pct")
DEG_FOLDER <- str_glue("results/{Site}/Productivity/TMF_AObs_deg/DEG_pct")
TW <- list()
tw_years <- seq(2000, 2015, by=3) 
for (i in 1:length(tw_years)) {
  start_year <- tw_years[i]
  start_index <- 1 + 23 * (start_year - 2000)  # Calculate the initial index
  end_index <- start_index + 115 - 1           # Calculate the final index
  TW_name <- paste(start_year, start_year + 5, sep="_")  # Name of the time window
  TW[[TW_name]] <- start_index:end_index
}

## 1.3 Charging data ----
# ---------------------------------------------------------------- - - -
roi <- vect(list.files(stringr::str_glue("data/{Site}"), pattern = "shp", full.names = TRUE))
evi_smoothed <- rast(file.path(EVI_FOLDER, "evi_smoothed.tif"))
evi_class_2000_2005 <- rast(file.path(EVI_FOLDER, stringr::str_glue("{Site}_evi_class_2000_2005.tif")))
evi_class_2003_2008 <- rast(file.path(EVI_FOLDER, stringr::str_glue("{Site}_evi_class_2003_2008.tif")))
evi_class_2006_2011 <- rast(file.path(EVI_FOLDER, stringr::str_glue("{Site}_evi_class_2006_2011.tif"))) 
evi_class_2009_2014 <- rast(file.path(EVI_FOLDER, stringr::str_glue("{Site}_evi_class_2009_2014.tif")))
evi_class_2012_2017 <- rast(file.path(EVI_FOLDER, stringr::str_glue("{Site}_evi_class_2012_2017.tif")))
evi_class_2015_2020 <- rast(file.path(EVI_FOLDER, stringr::str_glue("{Site}_evi_class_2015_2020.tif")))

deg_1990_2005 <- rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2005.tif")))
deg_1990_2008 <- rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2008.tif")))
deg_1990_2011 <- rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2011.tif")))
deg_1990_2014 <- rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2014.tif")))
deg_1990_2017 <- rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2017.tif")))
deg_1990_2020 <- rast(file.path(DEG_FOLDER, str_glue("{Site}_ACC_deg_1990_2020.tif")))

##%######################################################%##
#                                                          #
#                2. Setting user interface                 ----
#                                                          #
##%######################################################%##
# CSS and HTML for user interface ---
# ---------------------------------------------------------------- - - -
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        height: 90vh;
        max-width: 100%;
        max-height: 90vh;
        overflow: auto;
      }
      .flex-container {
        display: flex;
        flex-direction: column;
        align-items: center;
        max-height: 80%;
        overflow: auto; 
      }
      #timeWindowContainer {
        width: 100%; 
        text-align: center; 
      }
      #timeWindow {
        width: auto; 
        display: inline-block;
        margin-bottom: 5px;
      }
      #rasterMap {
        max-width: 90%;
        max-height: 300px;
      }
      .main-title {
        font-size: 1.0em;
        margin-bottom: 5px;
        text-align: center;
        width: 100%; 
        margin-top: 5px; 
      }
    "))
  ),
  div(class='main-title', titlePanel("Analyse des tendances de l'EVI")),
  div(class = 'flex-container', # Div to centralize the items
      selectInput("timeWindow", "Choisissez la période :", choices = names(TW)),
      leafletOutput("map", height = "350px"), 
      plotlyOutput("regressionPlot"),
      tableOutput("degradationValues")
  )
)

##%######################################################%##
#                                                          #
#                 3. Setting Shiny server                  ----
#                                                          #
##%######################################################%##
## 3.1 Server function ----
# ---------------------------------------------------------------- - - -
server <- function(input, output, session) {
  
  # Initialize the reactive variable to store the pixel data
  pixelData <- reactiveVal(NULL)
  
  # Initialize the reactive variable to store the degradation data
  degradationData <- reactiveVal(NULL)
  
  ### 3.1.1 Click definition ----
  # ---------------------------------------------------------------- - - -
  # Observe clicks on the leaflet map
  observe({
    click <- input$map_click
    if (is.null(click)) return()
    
    # Transform lat/lng coordinates into raster coordinates
    xy <- cbind(click$lng, click$lat)
    pixel_values <- extract(evi_smoothed, xy)
    pixelData(pixel_values) # Update reactive data
    
    # Extract degradation values for all periods
    deg_values <- sapply(list(deg_1990_2005, deg_1990_2008, deg_1990_2011, deg_1990_2014, deg_1990_2017, deg_1990_2020), function(rast) {
      extract(rast, xy)
    })
    degradationData(deg_values) # Update reactive data
  })
  
  ### 3.1.2 Setting color palettes ----
  # ---------------------------------------------------------------- - - -
  colors_class <- c("#FF0000", # Class -1 red
                    "#FFFFFF", # Class 0 white
                    "#00FF00") # Class 1 green
  
  color_pal_class <- colorFactor(palette = colors_class, domain = c(-1, 0, 1), na.color = "transparent")
  
  color_pal_deg <- colorBin(palette = c("#FFFFFF", "#d0e6f7ff", "#a7cdee", "#7fb4e5", "#559adb", 
                                        "#2979b9", "#23659a", "#1c507a", "#163b5b", "#11263c"), 
                            bins = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
                            na.color = "transparent")
  
  color_pal_evi <- colorNumeric(palette = c("transparent", "green", "yellow", "red"), 
                                domain = c(NA, range(values(evi_smoothed), na.rm = TRUE)),
                                na.color = "transparent")
  
  ### 3.1.3 Setting Leaflet mapping ----
  # ---------------------------------------------------------------- - - -
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery", options = providerTileOptions(opacity = 0.8)) %>% 
      addProviderTiles("OpenStreetMap.Mapnik", group = "OSM", options = providerTileOptions(opacity = 0.8)) %>%  
      addPolygons(data = roi, fillColor = "transparent", color = "blue", weight = 2, opacity = 1, group = "ROI") %>%
      addRasterImage(evi_smoothed, colors = color_pal_evi, opacity = 0.7,
                     project = FALSE, group = "EVI") %>%
      addRasterImage(evi_class_2000_2005, colors = color_pal_class, opacity = 0.8, project = FALSE, group = "EVI Class 2000-2005") %>%
      addRasterImage(evi_class_2003_2008, colors = color_pal_class, opacity = 0.8, project = FALSE, group = "EVI Class 2003-2008") %>%
      addRasterImage(evi_class_2006_2011, colors = color_pal_class, opacity = 0.8, project = FALSE, group = "EVI Class 2006-2011") %>%
      addRasterImage(evi_class_2009_2014, colors = color_pal_class, opacity = 0.8, project = FALSE, group = "EVI Class 2009-2014") %>%
      addRasterImage(evi_class_2012_2017, colors = color_pal_class, opacity = 0.8, project = FALSE, group = "EVI Class 2012-2017") %>%
      addRasterImage(evi_class_2015_2020, colors = color_pal_class, opacity = 0.8, project = FALSE, group = "EVI Class 2015-2020") %>%
      addRasterImage(deg_1990_2005, colors = color_pal_deg, opacity = 0.8, project = FALSE, group = "Degradation 1990-2005") %>%
      addRasterImage(deg_1990_2008, colors = color_pal_deg, opacity = 0.8, project = FALSE, group = "Degradation 1990-2008") %>%
      addRasterImage(deg_1990_2011, colors = color_pal_deg, opacity = 0.8, project = FALSE, group = "Degradation 1990-2011") %>%
      addRasterImage(deg_1990_2014, colors = color_pal_deg, opacity = 0.8, project = FALSE, group = "Degradation 1990-2014") %>%
      addRasterImage(deg_1990_2017, colors = color_pal_deg, opacity = 0.8, project = FALSE, group = "Degradation 1990-2017") %>%
      addRasterImage(deg_1990_2020, colors = color_pal_deg, opacity = 0.8, project = FALSE, group = "Degradation 1990-2020") %>%
      addLayersControl(
        overlayGroups = c("EVI",
                          "EVI Class 2000-2005",
                          "EVI Class 2003-2008",
                          "EVI Class 2006-2011",
                          "EVI Class 2009-2014",
                          "EVI Class 2012-2017",
                          "EVI Class 2015-2020",
                          "Degradation 1990-2005",
                          "Degradation 1990-2008",
                          "Degradation 1990-2011",
                          "Degradation 1990-2014",
                          "Degradation 1990-2017",
                          "Degradation 1990-2020"),
        baseGroups = c("Esri World Imagery", "OSM"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })
  
  # Hide desired layers after initial map loading
  observe({
    leafletProxy("map") %>%
      hideGroup("EVI Class 2000-2005") %>%
      hideGroup("EVI Class 2003-2008") %>%
      hideGroup("EVI Class 2006-2011") %>%
      hideGroup("EVI Class 2009-2014") %>%
      hideGroup("EVI Class 2012-2017") %>%
      hideGroup("EVI Class 2015-2020") %>%
      hideGroup("Degradation 1990-2005") %>%
      hideGroup("Degradation 1990-2008") %>%
      hideGroup("Degradation 1990-2011") %>%
      hideGroup("Degradation 1990-2014") %>%
      hideGroup("Degradation 1990-2017") %>%
      hideGroup("Degradation 1990-2020")
  })
  
  ### 3.3.4 Calculating the regression ----
  # ---------------------------------------------------------------- - - -
  output$regressionPlot <- renderPlotly({
    req(pixelData()) # Ensures that pixel values are available
    
    selected_TW <- TW[[input$timeWindow]]
    years <- seq(from = as.numeric(substr(input$timeWindow, 1, 4)), length.out = length(selected_TW), by = 1/23)
    
    values <- as.numeric(pixelData()[selected_TW])
    
    if (!all(is.na(values))) {
      data_for_plot <- data.frame(year = years, EVI_value = values)
      model <- lm(EVI_value ~ year, data = data_for_plot)
      data_for_plot$predicted_EVI_value <- predict(model, newdata = data_for_plot)
      
      # Calculating confidence interval
      conf_int <- predict(model, newdata = data_for_plot, interval = "confidence")
      data_for_plot$lower_ci <- conf_int[, "lwr"]
      data_for_plot$upper_ci <- conf_int[, "upr"]
      
      plot_ly() %>%
        add_ribbons(data = data_for_plot, x = ~year, ymin = ~lower_ci, ymax = ~upper_ci, name = "Intervalle de confiance",
                    fillcolor = 'rgba(204,204,204,0.5)', line = list(color = 'transparent')) %>%
        add_lines(data = data_for_plot, x = ~year, y = ~predicted_EVI_value, name = "Ligne de tendance",
                  line = list(color = 'blue')) %>%
        add_markers(data = data_for_plot, x = ~year, y = ~EVI_value, name = "EVI observé",
                    mode = "markers", marker = list(color = 'darkorange', size = 10)) %>%
        layout(yaxis = list(title = "Valeur de l'EVI",
                            titlefont = list(size = 16, color = 'black', family = 'Arial', style = 'bold'),
                            standoff = 15),  # Ajuste a distância como desejar
               xaxis = list(title = 'Year',
                            titlefont = list(size = 16, color = 'black', family = 'Arial', style = 'bold'),
                            standoff = 10))  # Ajuste a distância como desejar
    } else {
      return(NULL)
    }
  })
  
  ### 3.3.5 Displaying degradation values ----
  # ---------------------------------------------------------------- - - -
  output$degradationValues <- renderTable({
    req(degradationData())
    deg_values <- degradationData()
    periods <- c("1990-2005", "1990-2008", "1990-2011", "1990-2014", "1990-2017", "1990-2020")
    data.frame("Période" = periods, "Dégradation " = as.numeric(deg_values))
  })
}

##%######################################################%##
#                                                          #
#                 4. Calling UI/Server                     #
#                                                          #
##%######################################################%##
# 4. Calling UI/Server ----
# ---------------------------------------------------------------- - - -
shinyApp(ui, server)
