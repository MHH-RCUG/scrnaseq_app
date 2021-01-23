library(shinydashboard) #Basic
library(shinyjs) #Javascript
library(shinycssloaders) #Loaders for plots
library(Seurat) #Plots
library(ggplot2) #Plots
library(readxl) #Read Excel file
library(markdown) #Read Markdown
library(shinyalert) #Alerts
library(zip) #zip

script.dir <- dirname(sys.frame(1)$ofile)
setwd(script.dir) #set wd to file path

runApp(
  appDir = "shiny/",
  launch.browser = TRUE
)

