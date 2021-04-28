# Libraries
library(shiny)
library(shinyjs)
library(shinycssloaders)
library(shinythemes)
library(shinyalert)
library(Seurat)
library(ggplot2)
library(readxl)

# Increase max upload size
options(shiny.maxRequestSize = 10000*1024^2)

source("www/functions_plotting.R")

# Define global variables
features_names_ids = NULL
stored_FeaturePlots = NULL
stored_RidgePlotRaws = NULL
stored_RidgePlotNorms = NULL
stored_ViolinPlotRaws = NULL
stored_ViolinPlotNorms = NULL
stored_DotPlot = NULL
stored_Heatmap = NULL

# Define parameters
param = list()
param$col = "palevioletred"
param$col_cluster = "empty"

# UI------------------------------------------------------------------------------------------------
ui = tagList(
    # Initialze the use of shinyjs and shinyalert
    shinyjs::useShinyjs(),
    shinyalert::useShinyalert(),

    navbarPage(
        # Loading themes with shinythemes
        theme = shinythemes::shinytheme("flatly"),

        # The title of the app
        title = "scrnaseq_app",
        id = "tabs",

        # UI-Data Input---------------------------------------------------------------------------
        # Defining the tab "Input Data"
        tabPanel(
            title = "Data Input",
            icon = icon("upload"),
            fluidPage(
                fluidRow(
                    column(
                        width = 12,
                        align = "center",
                        tags$h3("1. Upload"),
                        fileInput(
                            width = "50%",
                            inputId = "file_rds",
                            label = "Upload Seurat file: (.rds)",
                            accept = ".rds",
                            buttonLabel = "Browse",
                            placeholder = "Please upload .rds file!"
                        ),
                        uiOutput(
                            outputId = "ui_datainput_excel"
                        ),
                        uiOutput(
                            outputId = "ui_datainput_buttons"
                        )
                    )
                )
            ),
        ),

        # UI-Settings-----------------------------------------------------------------------------
        tabPanel(
            title = "Settings",
            icon = icon("cog"),

            sidebarLayout(
                sidebarPanel(
                    tags$h1("Test")
                ),
                mainPanel(
                    fluidRow(
                        column(
                            width = 12,
                            splitLayout(
                                cellWidths = c("25%","25%"),
                                numericInput(
                                    inputId = "x_axis",
                                    label = "X-Axis (px):",
                                    value = 1024,
                                    min = 1,
                                    max = 3000,
                                    #width = "50%"
                                ),
                                numericInput(
                                    inputId = "y_axis",
                                    label = "Y-Axis (px):",
                                    value = 576,
                                    min = 1,
                                    max = 3000,
                                    #width = "50%"
                                )
                            ),
                            numericInput(
                              inputId = "res",
                              label = "Resolution of plot, in pixels per inch:",
                              value = 144,
                              min = 1,
                              max = 3000
                            ),
                        ),
                    ),
                    fluidRow(
                        column(
                            width = 12,
                            ## Action buttons to set the aspect ratio
                            actionButton(inputId = "ratio_512",
                                         label = "512x288 px",
                                         class = "btn btn-outline-secondary"),
                            actionButton(inputId = "ratio_768",
                                         label = "768x432 px",
                                         class = "btn btn-outline-secondary"),
                            actionButton(inputId = "ratio_1024",
                                         label = "1024x576 px",
                                         class = "btn btn-outline-secondary"),
                            actionButton(inputId = "ratio_1152",
                                         label = "1152x648 px",
                                         class = "btn btn-outline-secondary"),
                            actionButton(inputId = "ratio_1280",
                                         label = "1280x720 px",
                                         class = "btn btn-outline-secondary"),
                            actionButton(inputId = "ratio_1920",
                                         label = "1920x1080 px",
                                         class = "btn btn-outline-secondary")
                        )
                    ),
                    tags$div()
                )
            )
        ),
        # UI-Plots--------------------------------------------------------------------------------
        tabPanel(
            title = "Plots",
            icon = icon("bar-chart"),
            
            # Panel with tabs, to switch between plot types
            tabsetPanel(
                tabPanel(
                    title = "Feature Plots",
                    uiOutput(
                      outputId = "ui_feature"
                    )
                ),
                tabPanel(
                    title = "Ridge Plots (Raw)",
                    uiOutput(
                        outputId = "ui_ridge_raw"
                    )
                ),
                tabPanel(
                    title = "Ridge Plots (Norm)",
                    uiOutput(
                        outputId = "ui_ridge_norm"
                    )
                ),
                tabPanel(
                    title = "Violin Plots (Raw)",
                    uiOutput(
                        outputId = "ui_vln_raw"
                    )
                ),
                tabPanel(
                    title = "Violin Plots (Norm)",
                    uiOutput(
                        outputId = "ui_vln_norm"
                    )
                ),
                tabPanel(
                    title = "Dotplot",
                    shinycssloaders::withSpinner(
                        plotOutput("plot_dotplot")
                    )
                ),
                tabPanel(
                    title = "Heatmap",
                    shinycssloaders::withSpinner(
                        plotOutput("plot_heatmap")
                    )
                )
            )
        ),
        # UI-Download-----------------------------------------------------------------------------
        tabPanel(
            title = "Download",
            icon = icon("download")
        ),
        # UI-Help---------------------------------------------------------------------------------
        tabPanel(
            title = "Help",
            icon = icon("question"),
            tabsetPanel(
                tabPanel("1"),
                tabPanel("2")
            ),
        )
    )
)

# Server--------------------------------------------------------------------------------------------
server = function(input, output, session) {

  source("www/global.R", local = TRUE)

    # Server-Hide Tabs----------------------------------------------------------------------------
    # Hide unused tabs in the beginning
    hideTab(
        inputId = "tabs",
        target = "Settings"
    )
    hideTab(
        inputId = "tabs",
        target = "Plots"
    )
    hideTab(
        inputId = "tabs",
        target = "Download"
    )

    # Server-Data Input---------------------------------------------------------------------------
    # Upload .rds file
    sc = reactive({
        inFile = input$file_rds
        if (is.null(inFile)) {
            tmp = NULL
        } else {
            shinyalert(
                title = "Please wait!",
                text = "Upload complete! Please wait while the data is being processed!",
                size = "s",
                type = "info",
                showConfirmButton = FALSE
            )
            tmp = readRDS(inFile$datapath)
        }
        return(tmp)
    })

    # Process uploaded .rds file
    observeEvent(sc(), {
        if (!is.null(sc())) {
            features_names_ids <<-
                paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][, 1], sep = "")
            param$col_clusters <<- as.vector(sc()@misc[["colors"]][["seurat_clusters"]])
        }
        output$ui_datainput_excel = renderUI(
            tagList(
                tags$hr(),
                tags$h3("2. Selection"),
                selectizeInput(
                    width = "50%",
                    inputId = "select_genes",
                    label = "Select genes:",
                    choices = NULL,
                    multiple = TRUE
                ),
                fileInput(
                    width = "50%",
                    inputId = "xlsx_file",
                    label = "Select through excel file: (.xlsx)",
                    accept = ".xlsx",
                    buttonLabel = "Browse",
                    placeholder = "No file selected"
                )
            )
        )
        updateSelectizeInput(
            inputId = "select_genes",
            label = "Select Genes:",
            choices = features_names_ids,
            server = TRUE
        )
        shinyjs::delay(500, shinyjs::runjs("swal.close();"))
    })

    # Excel file upload
    excel_genes = reactive({
        inFile = input$xlsx_file
        if (is.null(inFile)) {
            tmp = NULL
        } else {
            tmp = readxl::read_excel(
                path = inFile$datapath,
                sheet = 1,
                range = cell_cols("A:B"),
                col_names = FALSE
            )
            Sys.sleep(1)
            shinyalert(
                title = "Please wait!",
                text = "Upload complete! Please wait while the data is being processed!",
                size = "s",
                type = "info",
                showConfirmButton = FALSE
            )
        }
        return(tmp)
    })

    # Excel file processing
    observeEvent(excel_genes(), {
        tryCatch(
            {
                excel_list = unlist(excel_genes()[, 1])
                x = features_names_ids[unlist(lapply(excel_list, function(one_gene)
                    grep(one_gene, features_names_ids)))]
            },
            error = function(e){
                message("Excel file upload: There was an error")
                message(e)
                shinyalert(
                    title = "Error!",
                    text = "The Excel file was empty!",
                    size = "s",
                    type = "error",
                    showConfirmButton = TRUE
                )
            },
            warning = function(w){
                message("Excel file upload: There was a warning message.")
                message(w)
            },
            finally = {
                suppressWarnings(
                    updateSelectizeInput(
                        inputId = "select_genes",
                        label = "Select Genes:",
                        choices = features_names_ids,
                        selected = x,
                        server = TRUE
                    )
                )
                shinyjs::delay(500, shinyjs::runjs("swal.close();"))
            }
        )
    })

    # Render UI for buttons when genes are selected
    observeEvent(input$select_genes,{
        output$ui_datainput_buttons = renderMenu(
            tagList(
                tags$hr(),
                tags$h3("3. Render"),
                actionButton(
                    width = "50%",
                    inputId = "render_with_defaults",
                    label = "Render with default settings!",
                    class = "btn btn-success btn-lg btn-block"
                ),
                actionButton(
                    width = "50%",
                    inputId = "goto_settings",
                    label = "Go to settings!",
                    class = "btn btn-success btn-lg btn-block"
                )
            )
        )
    })
    
    # Excecute when button "render with default settings" is pressed
    observeEvent(input$render_with_defaults,{
        render_plots()
        showTab(
            inputId = "tabs",
            target = "Settings"
        )
        showTab(
            inputId = "tabs",
            target = "Plots",
            select = TRUE
        )
        showTab(
            inputId = "tabs",
            target = "Download"
        )
    })

    # Action button to go to settings
    observeEvent(input$goto_settings,{
        showTab(
            inputId = "tabs",
            target = "Settings",
            select = TRUE,
        )
    })

    # Server-Settings-----------------------------------------------------------------------------


    # Server-Plots--------------------------------------------------------------------------------

}

# Run the application
shinyApp(ui = ui, server = server)
