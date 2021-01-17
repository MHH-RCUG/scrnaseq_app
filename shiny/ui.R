library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(Seurat)
library(ggplot2)
library(readxl)
library(markdown)
library(shinyalert)
library(zip)

# Define UI for application that draws a histogram
ui <- tagList(
  dashboardPage(
    dashboardHeader(title = "scrnaseq_app"),
    dashboardSidebar(
      sidebarMenu(
        id = "tabs",
        menuItem("Input Data", tabName = "datainput", icon = icon("upload")),
        menuItemOutput("settings"),
        menuItemOutput("plots"),
        menuItemOutput("download"),
        menuItem("Help", tabName = "help", icon = icon("question")),
        menuItem("About", tabName = "about", icon = icon("info"))
      )#sidebarMenu
    ),#dashboardSidebar

    dashboardBody(
      tags$head(tags$style(HTML('.skin-blue .wrapper {background-color: #ecf0f5;}'))),

      useShinyjs(),
      useShinyalert(),
      tabItems(
        source("ui-tab-datainput.R", local = TRUE)$value,
        source("ui-tab-settings.R", local = TRUE)$value,
        tabItem(tabName = "featureplots", uiOutput("ui_feature")),
        tabItem(tabName = "ridgeplots_raw", uiOutput("ui_ridge_raw")),
        tabItem(tabName = "ridgeplots_norm", uiOutput("ui_ridge_norm")),
        tabItem(tabName = "violinplots_raw", uiOutput("ui_vln_raw")),
        tabItem(tabName = "violinplots_norm", uiOutput("ui_vln_norm")),
        tabItem(
          tabName = "dotplot",
          shinycssloaders::withSpinner(
            plotOutput("plot_dotplot")
          )
        ),
        tabItem(
          tabName = "heatmap",
          shinycssloaders::withSpinner(
            plotOutput("plot_heatmap")
          )
        ),
        source("ui-tab-download.R", local = TRUE)$value,
        source("ui-tab-help.R", local = TRUE)$value
        )#tabItems
      )#dashboardbody
    ),#dashboardPage

    tags$footer(
      wellPanel(
        HTML(
        '
        <p align="center" width="4">scrnaseq_app</p>
        <p align="center" width="4">Github: <a href="https://github.com/MHH-RCUG/scrnaseq_app">https://github.com/MHH-RCUG/scrnaseq_app</a></p>
        <p align="center" width="4">Created by: <a href="https://www.mhh.de/genomics/kontakt">MHH RCUG</a> </p>
        '
        )#HTML
        ),#wellPanel
      )#footer
  )#UI
