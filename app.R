# Author: Marius Rueve

# ! This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.

# ! For convenience click "Run External" to run in you default browser

# Global --------------------
# Load packages
library(shiny)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(markdown)
library(readxl)
library(shinyalert)
library(shinyjs)
library(zip)

# Set max upload size to 300 MB
options(shiny.maxRequestSize = 5000 * 1024 ^ 2)

# Parameter list for convenience
param = list()
param$col = "palevioletred"

# Load data if upload is broken
#sc = readRDS(file = "pbmc_2020-06-09.rds")
#features_names_ids = paste(rownames(sc[["RNA"]][[]]), "_", sc[["RNA"]][[]][,1], sep = "")

# Empty object so Input UI gets rendered
features_names_ids = NULL

# Define global variable to store plots
stored_FeaturePlots = NULL
stored_RidgePlotRaws = NULL
stored_RidgePlotNorms = NULL
stored_ViolinPlotRaws = NULL
stored_ViolinPlotNorms = NULL
stored_DotPlot = NULL
stored_Heatmap = NULL

# UI --------------------
ui = fluidPage(

  #theme = shinytheme("paper"),

  # Place RCUG logo at top right
  img(
    src = "RCUG_Logo.png",
    height = "10%",
    width = "10%",
    align = "right"
  ),

  # Place MHH logo at top right
  img(
    src = "Logo_engl_schwarz.png",
    height = "15%",
    width = "15%",
    align = "right"
  ),

  # Application title
  titlePanel("scrnaseq_app"),

  # Sidebar
  sidebarPanel(
    width = 3,
    h4("1. Upload:"),

    # Upload of .rds file
    fileInput(
      inputId = "rds_file",
      label = "Upload Seurat file: (.rds)",
      accept = ".rds",
      buttonLabel = "Browse",
      placeholder = "Please upload .rds file!"
    ),

    # UI will be rendered once the .rds file has been uploaded
    uiOutput("insert_ui"),
  ),

  # Main panel
  mainPanel(

    useShinyjs(),
    useShinyalert(),

    # Set tabs
    tabsetPanel(
      type = "tabs",
      tabPanel("FeaturePlot", uiOutput("ui_feature")),
      tabPanel("RidgePlot (Raw)", uiOutput("ui_ridge_raw")),
      tabPanel("RidgePlot (Norm)", uiOutput("ui_ridge_norm")),
      tabPanel("ViolinPlot (Raw)", uiOutput("ui_vln_raw")),
      tabPanel("ViolinPlot (Norm)", uiOutput("ui_vln_norm")),
      tabPanel("DotPlot", plotOutput("plot_dotplot")),
      tabPanel("Heatmap", plotOutput("plot_heatmap")),
      tabPanel("Help", includeMarkdown("HELP.md"))
    )
  )
)


# Server --------------------
server = function(input, output, session) {

  output$insert_ui = renderUI({
    req(input$rds_file)
    tagList(
      tags$hr(),
      h4("2. Select genes:"),
      selectInput(
        "genes",
        "Select through list:",
        features_names_ids,
        multiple = TRUE
      ),

      # Select genes
      actionButton("clear_selection", "Clear selection"),
      fileInput(
        inputId = "xlsx_file",
        label = "Select through excel file: (.xlsx)",
        accept = ".xlsx",
        buttonLabel = "Browse",
        placeholder = "No file selected"
      ),

      tags$hr(),

      fluidRow(
        column(
          width = 6,
          numericInput(
            inputId = "x_axis",
            label = "X-Axis (px):",
            value = 1280,
            min = 1,
            max = 3000
          )
        ),

        column(
          width = 6,
          ofset = 3,
          numericInput(
            inputId = "y_axis",
            label = "Y-Axis (px):",
            value = 720,
            min = 1,
            max = 3000
          )
        )
      ),
      
      actionButton(inputId = "preset_512",
                   label = "512x288 px"),
      
      actionButton(inputId = "preset_1024",
                   label = "1024x576 px"),

      actionButton(inputId = "restore_axes",
                   label = "Default settings"),

      tags$hr(),

      h4("3. Download:"),

      textInput(
        inputId = "archive_download",
        label = "Enter name of archive (.zip):",
        value = paste0("Download", "_", Sys.Date())
      ),

      fluidRow(
        column(
          5,
          checkboxInput("check_featureplot", "FeaturePlot", value = TRUE)
        ),
        column(
          5,
          checkboxInput("check_ridgeplot_raw", "RidgePlot (Raw)", value = TRUE)
        ),
        column(
          5,
          checkboxInput("check_ridgeplot_norm", "RidgePlot (Norm)", value = TRUE)
        ),
        column(
          5,
          checkboxInput("check_vlnplot_raw", "ViolinPlot (Raw)", value = TRUE)
        ),
        column(
          5,
          checkboxInput("check_vlnplot_norm", "ViolinPlot (Norm)", value = TRUE)
        ),
        column(5, checkboxInput("check_dotplot", "DotPlot", value = TRUE))
      ),

      downloadButton("download_plots")
    )
  })

  # Upload ====================
  # .rds file upload
  sc = reactive({
    inFile = input$rds_file
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

  observeEvent(sc(), {
    if (!is.null(sc())) {
      features_names_ids <<-
        paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][, 1], sep = "")
      #print(features_names_ids)
      updateSelectInput(session, "genes", "Select Genes:", features_names_ids)
    }
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

  # Observe event: when excel file has been uploaded and read
  observeEvent(excel_genes(), {
    if (length(excel_genes()) == 0) {
      shinyjs::runjs("swal.close();")
      shinyjs::delay(1000)
      shinyalert(
        title = "Error!",
        text = "The Excel file was empty!",
        size = "s",
        type = "error",
        showConfirmButton = TRUE
      )
      return(NULL)
    } else{
      excel_list = unlist(excel_genes()[, 1])
      x = features_names_ids[unlist(lapply(excel_list, function(one_gene)
        grep(one_gene, features_names_ids)))]
      updateSelectInput(
        session = session,
        inputId = "genes",
        label = "Select Genes:",
        choices = features_names_ids,
        selected = x
      )
      shinyjs::delay(500, shinyjs::runjs("swal.close();"))
    }
  })

  # If Button "clear selection" is pressed the selection is null
  observeEvent(input$clear_selection, {
    updateSelectInput(session, "genes", "Select Genes:", features_names_ids)
  })
  
  observeEvent(input$preset_1024, {
    updateNumericInput(
      session,
      "x_axis",
      "X-Axis (px):",
      value = 1024,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session,
      "y_axis",
      "Y-Axis (px):",
      value = 576,
      min = 1,
      max = 3000
    )
  })
  
  observeEvent(input$preset_512, {
    updateNumericInput(
      session,
      "x_axis",
      "X-Axis (px):",
      value = 512,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session,
      "y_axis",
      "Y-Axis (px):",
      value = 288,
      min = 1,
      max = 3000
    )
  })

  # Button to restore default settings for axes
  observeEvent(input$restore_axes, {
    updateNumericInput(
      session,
      "x_axis",
      "X-Axis (px):",
      value = 1280,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session,
      "y_axis",
      "Y-Axis (px):",
      value = 720,
      min = 1,
      max = 3000
    )
  })


  # Plots =================================

  # Rendering Plots, location of Seurat functions
  observeEvent(input$genes, {
    for (i in 1:length(input$genes)) {
      # Feature Plots
      stored_FeaturePlots[[i]] <<- # <<- for global assignments
        Seurat::FeaturePlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          cols = c("lightgrey", param$col),
          label = TRUE
        ) +
        theme_light() + theme(panel.border = element_blank())

      # Ridge Plots Raw
      stored_RidgePlotRaws[[i]] <<-
        Seurat::RidgePlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          assay = "RNA",
          slot = "counts"
        ) +
        theme_light() + theme(panel.border = element_blank()) +
        labs(color = "Cell identity", fill = "Cell identity") +
        ylab("Cluster")

      # Ridge Plots Norm
      stored_RidgePlotNorms[[i]] <<-
        Seurat::RidgePlot(sc(),
                          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
                          slot = "data") +
        theme_light() + theme(panel.border = element_blank()) +
        labs(color = "Cell identity", fill = "Cell identity") +
        ylab("Cluster")

      # Violin Plots Raw
      stored_ViolinPlotRaws[[i]] <<-
        Seurat::VlnPlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          assay = "RNA",
          slot = "counts",
          pt.size = 0.2
        ) +
        theme_light() + theme(panel.border = element_blank()) +
        labs(color = "Cell identity", fill = "Cell identity") +
        xlab("Cluster")

      # Violin Plots Norm
      stored_ViolinPlotNorms[[i]] <<-
        Seurat::VlnPlot(sc(),
                        features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
                        pt.size = 0.2) +
        theme_light() + theme(panel.border = element_blank()) +
        labs(color = "Cell identity", fill = "Cell identity") +
        xlab("Cluster")

      # DotPlot
      stored_DotPlot <<-
        Seurat::DotPlot(
          sc(),
          features = unlist(strsplit(input$genes, "_"))[c(T, F)],
          cols = c("lightgrey", param$col)
        ) +
        theme_light() + theme(panel.border = element_blank()) +
        ylab("Cluster") +
        theme(axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = .5
        ))

      tryCatch({
        # Heatmap
        stored_Heatmap <<-
          Seurat::DoHeatmap(
            sc(),
            features = unlist(strsplit(input$genes, "_"))[c(T, F)],
            group.colors = param$col_clusters,
            slot = "counts"
          ) +
          NoLegend()
      },
      warning = function(w){
          showNotification(
            ui = "One or more of the selected genes could not be found in the 
                  default search locations. These gene(s) will not show up in the Heatmap!",
            duration = NULL,
            id = "heatmap_warning",
            type = "warning"
          )
      }
      )
    }
  })


  # UI FeaturePlot #############################
  output$ui_feature = renderUI({
    tmp = list()

    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      tmp[[i]] =  plotOutput(
        outputId = paste0("plot_feature", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(tmp)
  })


  # UI RidgePlot Raw #############################
  output$ui_ridge_raw = renderUI({
    tmp = list()

    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      tmp[[i]] =  plotOutput(
        outputId = paste0("plot_ridge_raw", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(tmp)
  })


  # UI RidgePlot Normalised #############################
  output$ui_ridge_norm = renderUI({
    tmp = list()

    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      tmp[[i]] =  plotOutput(
        outputId = paste0("plot_ridge_norm", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(tmp)
  })


  # UI ViolinPlot Raw #############################
  output$ui_vln_raw = renderUI({
    tmp = list()

    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      tmp[[i]] =  plotOutput(
        outputId = paste0("plot_vln_raw", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(tmp)
  })


  # UI ViolinPlot Normalised #############################
  output$ui_vln_norm = renderUI({
    tmp = list()

    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      tmp[[i]] =  plotOutput(
        outputId = paste0("plot_vln_norm", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(tmp)
  })


  # renderPlots() --> plotOutput() --> renderUI()
  observe({
    for (i in 1:length(input$genes)) {
      local({
        # because expressions are evaluated at app init
        ii = i

        # This chunks takes the stored plots from the global variable and renders a reactive plots that is
        # suitable for assigning to an `output` slot; passing vector to UI
        output[[paste0("plot_feature", ii)]] = renderPlot(stored_FeaturePlots[[ii]])

        output[[paste0("plot_ridge_raw", ii)]] = renderPlot(stored_RidgePlotRaws[[ii]])

        output[[paste0("plot_ridge_norm", ii)]] = renderPlot(stored_RidgePlotNorms[[ii]])

        output[[paste0("plot_vln_raw", ii)]] = renderPlot(stored_ViolinPlotRaws[[ii]])

        output[[paste0("plot_vln_norm", ii)]] = renderPlot(stored_ViolinPlotNorms[[ii]])

        output$plot_dotplot = renderPlot(stored_DotPlot,
                                         width = input$x_axis,
                                         height = input$y_axis)

        output$plot_heatmap = renderPlot(stored_Heatmap,
                                         width = input$x_axis,
                                         height = input$y_axis)

      })
    }
  })


  # Downloads =================================
  output$download_plots = downloadHandler(
    filename = function() {
      # Takes the filename from the input of the user
      paste0(input$archive_download, ".zip")
    },
    content = function(file) {
      #temporarily switch to the temp dir, in case you do not have write permission to the current working directory
      owd = setwd(tempdir())
      on.exit(setwd(owd))

      # Idea to check if all devices are off on exit
      # on.exit({
      #  setwd(owd)
      #  while length(dev.list())>1; do
      #    dev.off()
      # })

      files = NULL
      on.exit(unlink(files))

      if (length(input$genes) == 0) {
        # Message if download fails due to no genese selected
        shinyalert(
          title = "Download failed!",
          text = "Failed to download plots! Please make sure to upload .rds file and select genes before downloading!",
          size = "s",
          type = "error"
        )
        return(NULL)
      } else{
        # Message that the archive is being created and the download will start soon
        shinyalert(
          title = "Please wait!",
          text = "The creation of files (.png and .pdf) can take a while. The download will start, once everything is ready.",
          size = "s",
          type = "info",
          showConfirmButton = FALSE
        )
      }


      # Download FeaturPlot #############################
      if (input$check_featureplot == TRUE) {
        # PNG
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("FeaturePlot_", input$genes[[i]], ".png")

          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(stored_FeaturePlots[[i]])
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF
        fileName_pdf = "FeaturePlot.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / 96),
          height = (input$y_axis / 96)
        )
        for (i in 1:length(input$genes)) {
          print(stored_FeaturePlots[[i]])
        }
        dev.off()
        files = c(fileName_pdf, files)
      }


      # Download RidgePlot Raw #############################
      if (input$check_ridgeplot_raw == TRUE) {
        # PNG
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("RidgePlot_Raw_", input$genes[[i]], ".png")

          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(stored_RidgePlotRaws[[i]])
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF
        fileName_pdf = "RidgePlot_Raw.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / 96),
          height = (input$y_axis / 96)
        )
        for (i in 1:length(input$genes)) {
          print(stored_RidgePlotRaws[[i]])
        }
        dev.off()
        files = c(fileName_pdf, files)
      }


      # Download RidgePlot Normalised #############################
      if (input$check_ridgeplot_norm == TRUE) {
        # PNG
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("RidgePlot_Norm_", input$genes[[i]], ".png")

          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(stored_RidgePlotNorms[[i]])
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF
        fileName_pdf = "RidgePlot_Norm.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / 96),
          height = (input$y_axis / 96)
        )
        for (i in 1:length(input$genes)) {
          print(stored_RidgePlotNorms[[i]])
        }
        dev.off()
        files = c(fileName_pdf, files)
      }


      # Download ViolinPlot Raw #############################
      if (input$check_vlnplot_raw == TRUE) {
        # PNG ViolinPlot Raw
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("ViolinPlot_Raw_", input$genes[[i]], ".png")

          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(stored_ViolinPlotRaws[[i]])
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF ViolinPlot Raw
        fileName_pdf = "ViolinPlot_Raw.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / 96),
          height = (input$y_axis / 96)
        )
        for (i in 1:length(input$genes)) {
          print(stored_ViolinPlotRaws[[i]])
        }
        dev.off()
        files = c(fileName_pdf, files)
      }


      # Download ViolinPlot Normalised #############################
      if (input$check_vlnplot_norm == TRUE) {
        # PNG ViolinPlot Normalised
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("ViolinPlot_Norm_", input$genes[[i]], ".png")

          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(stored_ViolinPlotNorms[[i]])
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF ViolinPlot Normalised
        fileName_pdf = "ViolinPlot_Norm.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / 96),
          height = (input$y_axis / 96)
        )
        for (i in 1:length(input$genes)) {
          print(stored_ViolinPlotNorms[[i]])
        }
        dev.off()
        files = c(fileName_pdf, files)
      }


      # Download DotPlot #############################
      if (input$check_dotplot == TRUE) {
        # PNG
        fileName_png = "DotPlot.png"

        png(fileName_png,
            width = input$x_axis,
            height = input$y_axis)
        print(stored_DotPlot)
        dev.off()
        files = c(fileName_png, files)
        # PDF
        fileName_pdf = "DotPlot.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / 96),
          height = (input$y_axis / 96)
        )
        print(stored_DotPlot)
        dev.off()
        files = c(fileName_pdf, files)
      }


      # Create zip file for Download, uses array of files
      zip(zipfile = file, files =  files)
      shinyjs::runjs("swal.close();")
    },
    contentType = "application/zip"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
