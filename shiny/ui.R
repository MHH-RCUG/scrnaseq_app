
source("functions_plotting.R")

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
        # Sourcing files in order to associate them with a "tabName"
        source("ui-tab-datainput.R", local = TRUE)$value,
        source("ui-tab-settings-general.R", local = TRUE)$value,
        source("ui-tab-settings-colors.R", local = TRUE)$value,
        source("ui-tab-settings-heatmap.R", local = TRUE)$value,

        # "tabItem"s for plots
        tabItem(
          tabName = "featureplots",
          uiOutput("ui_feature")
          ),
        tabItem(
          tabName = "ridgeplots_raw",
          uiOutput("ui_ridge_raw")
          ),
        tabItem(
          tabName = "ridgeplots_norm",
          uiOutput("ui_ridge_norm")
          ),
        tabItem(
          tabName = "violinplots_raw",
          uiOutput("ui_vln_raw")
          ),
        tabItem(
          tabName = "violinplots_norm",
          uiOutput("ui_vln_norm")
          ),
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

  # Footer shown at the bottom of every page
    tags$footer(
      wellPanel(
        HTML(
        '
        <p align="center" width="4">scrnaseq_app</p>
        <p align="center" width="4">Github: <a href="https://github.com/MHH-RCUG/scrnaseq_app">https://github.com/MHH-RCUG/scrnaseq_app</a></p>
        <p align="center" width="4">Created by: <a href="https://www.mhh.de/genomics/kontakt">MHH RCUG</a> </p>
        <p align="center" width="4"><a href="https://www.mhh.de/genomics"> <img src="RCUG_Logo.png" height="15%" width="15%" alt="RCUG"></a><a href="https://www.mhh.de/"> <img src="Logo_engl_schwarz.png" height="15%" width="15%" alt="RCUG"></a></p>
        '
        )#HTML
      )#wellPanel
    )#footer
  )#UI
