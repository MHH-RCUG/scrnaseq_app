# Author: Marius Rueve
# Global ---------------------------------
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
options(shiny.maxRequestSize = 300 * 1024 ^ 2)

# Parameter list for convenience
param = list()
param$col = "palevioletred"

# Load data if upload is broken
#sc = readRDS(file = "pbmc_2020-06-09.rds")
#features_names_ids = paste(rownames(sc[["RNA"]][[]]), "_", sc[["RNA"]][[]][,1], sep = "")

# Empty object so Input UI gets rendered
features_names_ids = NULL

# Define global variable for plots
plots_FeaturePlot = NULL
plots_RidgePlotRaw = NULL
plots_RidgePlotNorm = NULL
plots_ViolinPlotRaw = NULL
plots_ViolinPlotNorm = NULL
plots_DotPlot = NULL
plots_Heatmap = NULL

# UI ---------------------------------
ui = fluidPage(
  
  #theme = shinytheme("paper"),
  img(src = "Logo_engl_schwarz.png", height="15%", width="15%", align = "right"),
  titlePanel("scrnaseq_app"),
  
  sidebarPanel(
    width = 3,
    h4("1. Upload:"),
    
    fileInput(
      inputId = "rds_file",
      label = "Upload Seurat file: (.rds)",
      accept = ".rds",
      buttonLabel = "Browse",
      placeholder = "Please upload .rds file!"
    ), # fileInput: Upload of .rds file
    
    uiOutput("insert_ui"),
  ), #sidebarPanel
  
  # Tabs
  mainPanel(
    useShinyjs(),
    useShinyalert(),
    tabsetPanel(
      type = "tabs",
      tabPanel("FeaturePlot", uiOutput("ui_feature")),
      tabPanel("RidgePlot (Raw)", uiOutput("ui_ridge_raw")),
      tabPanel("RidgePlot (Norm)", uiOutput("ui_ridge_norm")),
      tabPanel("ViolinPlot (Raw)", uiOutput("ui_vln_raw")),
      tabPanel("ViolinPlot (Norm)", uiOutput("ui_vln_norm")),
      tabPanel("DotPlot", plotOutput("plot_dotplot")),
      tabPanel("Heatmap", plotOutput("plot_heatmap")),
      tabPanel("Help", includeMarkdown("README.md"))
    ) #tabPanel: UI for the creation of tabs and which content to display
  ) #mainPanel
) #fluidPage


# Server ---------------------------------
server = function(input, output, session) {
  
  output$insert_ui = renderUI({
    req(input$rds_file)
    tagList(
      tags$hr(),
      h4("2. Select genes:"),
      selectInput("genes", "Select through list:", features_names_ids, multiple = TRUE), # Select genes 
      actionButton("clear_selection", "Clear selection"),
      fileInput(
        inputId = "xlsx_file",
        label = "Select through excel file: (.xlsx)",
        accept = ".xlsx",
        buttonLabel = "Browse",
        placeholder = "No file selected"
      ), #fileInput: uplaod of an excel file (.xlsx) to select genes by ensemblID
      
      tags$hr(),
      
      fluidRow(column(
        6,
        numericInput(
          "x_axis",
          "X-Axis (px):",
          value = 1280,
          min = 1,
          max = 3000
        )
      ), #column: value of x-axis
      column(
        6,
        ofset = 3,
        numericInput(
          "y_axis",
          "Y-Axis (px):",
          value = 720,
          min = 1,
          max = 3000
        )
      )), #column: value of y-axis
      actionButton("restore_axes", "Default settings"), # Restore default values of axes (px)
      tags$hr(),
      h4("3. Download:"),
      textInput(
        "archive_download",
        "Enter name of archive (.zip):",
        value = paste0("Download", "_", Sys.Date())
      ), #textInput: give name of .zip file
      fluidRow(
        column(5,checkboxInput("check_featureplot", "FeaturePlot", value = TRUE)),
        column(5,checkboxInput("check_ridgeplot_raw", "RidgePlot (Raw)", value = TRUE)),
        column(5,checkboxInput("check_ridgeplot_norm", "RidgePlot (Norm)", value = TRUE)),
        column(5,checkboxInput("check_vlnplot_raw", "ViolinPlot (Raw)", value = TRUE)),
        column(5,checkboxInput("check_vlnplot_norm", "ViolinPlot (Norm)", value = TRUE)),
        column(5, checkboxInput("check_dotplot", "DotPlot", value = TRUE))
        ), #fluidRow: Cointains checkboxes wether or not to download certain plot types
      downloadButton("download_plots")
    ) #tagList
  }) #renderUI: dynamically expands the UI when the user uploads an .rds file
  
  # Upload =================================
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
  
  observeEvent(sc(),{
    if (!is.null(sc())) {
      features_names_ids <<- paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][, 1], sep = "")
      #print(features_names_ids)
      updateSelectInput(session, "genes", "Select Genes:", features_names_ids)
    }
    shinyjs::delay(500,shinyjs::runjs("swal.close();"))
  })
  
  # Excel file upload
  excel_genes = reactive({
    inFile = input$xlsx_file
    if (is.null(inFile)) {
      tmp = NULL
    } else {
      tmp = readxl::read_excel(path = inFile$datapath,
                               sheet = 1,
                               range = cell_cols("A:B"),
                               col_names = FALSE)
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
  
  observeEvent(excel_genes(),{
    if(length(excel_genes()) == 0){
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
    }else{
      excel_list = unlist(excel_genes()[,1])
      x = features_names_ids[unlist(lapply(excel_list, function(one_gene)grep(one_gene, features_names_ids)))]
      updateSelectInput(session = session,
                        inputId = "genes",
                        label = "Select Genes:",
                        choices = features_names_ids, 
                        selected = x)
      # updateTextInput(session = session,
      #                 inputId = "archive_download",
      #                 label = "Enter name of archive (.zip):",
      #                 value = input$xlsx_file$name
      #                 )
      shinyjs::delay(500,shinyjs::runjs("swal.close();"))
    }
  }) #observerEvent
  
  observeEvent(input$clear_selection, {
    updateSelectInput(session, "genes", "Select Genes:", features_names_ids)
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
  ## Rendering Plots
  observeEvent(input$genes, {
    for (i in 1:length(input$genes)) {
      
      # Feature Plots
      plots_FeaturePlot[[i]] <<- # <<- for global assignments
        Seurat::FeaturePlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          cols = c("lightgrey", param$col),
          label = TRUE
        ) +
        theme_light() + theme(panel.border = element_blank())
      
      # Ridge Plots Raw
      plots_RidgePlotRaw[[i]] <<- 
        Seurat::RidgePlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          assay = "RNA",
          slot = "counts"
        ) + 
        theme_light() + theme(panel.border = element_blank()) +
        labs(color="Cell identity", fill="Cell identity") + 
        ylab("Cluster")
      
      # Ridge Plots Norm
      plots_RidgePlotNorm[[i]] <<- 
        Seurat::RidgePlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          slot = "data"
        ) + 
        theme_light() + theme(panel.border = element_blank()) +
        labs(color="Cell identity", fill="Cell identity") + 
        ylab("Cluster")
      
      # Violin Plots Raw
      plots_ViolinPlotRaw[[i]] <<-
        Seurat::VlnPlot( 
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          assay = "RNA",
          slot = "counts",
          pt.size = 0.2
        ) + 
        theme_light() + theme(panel.border = element_blank()) +
        labs(color="Cell identity", fill="Cell identity") + 
        xlab("Cluster")
      
      # Violin Plots Norm
      plots_ViolinPlotNorm[[i]] <<-
        Seurat::VlnPlot(
          sc(),
          features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
          pt.size = 0.2
        ) + 
        theme_light() + theme(panel.border = element_blank()) +
        labs(color="Cell identity", fill="Cell identity") + 
        xlab("Cluster")
      
      # Dot Plot
      plots_DotPlot <<-
        Seurat::DotPlot(
          sc(),
          features = unlist(strsplit(input$genes, "_"))[c(T, F)],
          cols = c("lightgrey", param$col)
        ) + 
        theme_light() + theme(panel.border = element_blank()) + 
        ylab("Cluster") + 
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))
      
      # Heatmap
      plots_Heatmap <<-
        Seurat::DoHeatmap(
          sc(),
          features=unlist(strsplit(input$genes, "_"))[c(T, F)],
          group.colors=param$col_clusters,
          slot = "counts") +
        NoLegend()
    }
  })
  
  ### UI FeaturePlot ############################# 
  output$ui_feature = renderUI({
    
    out_feature = list()

    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      out_feature[[i]] =  plotOutput(
        outputId = paste0("plot_feature", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(out_feature)
  })
  
  ### UI RidgePlot Raw ############################# 
  output$ui_ridge_raw = renderUI({
    out_ridge_raw = list()
    
    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      out_ridge_raw[[i]] =  plotOutput(
        outputId = paste0("plot_ridge_raw", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(out_ridge_raw)
  })
  
  ### UI RidgePlot Normalised ############################# 
  output$ui_ridge_norm = renderUI({
    out_ridge_norm = list()
    
    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      out_ridge_norm[[i]] =  plotOutput(
        outputId = paste0("plot_ridge_norm", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(out_ridge_norm)
  })
  
  ### UI ViolinPlot Raw ############################# 
  output$ui_vln_raw = renderUI({
    out_vln_raw = list()
    
    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      out_vln_raw[[i]] =  plotOutput(
        outputId = paste0("plot_vln_raw", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(out_vln_raw)
  })
  
  ### UI ViolinPlot Normalised ############################# 
  output$ui_vln_norm = renderUI({
    out_vln_norm = list()
    
    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      out_vln_norm[[i]] =  plotOutput(
        outputId = paste0("plot_vln_norm", i),
        width = paste0(input$x_axis, "px"),
        height = paste0(input$y_axis, "px")
      )
    }
    return(out_vln_norm)
  })
  
  ### renderPlots() --> plotOutput() --> renderUI()
  observe({
    for (i in 1:length(input$genes)) {
      local({
        #because expressions are evaluated at app init
        ii = i
        output[[paste0("plot_feature", ii)]] = renderPlot(
            plots_FeaturePlot[[ii]]
        )
        
        output[[paste0("plot_ridge_raw", ii)]] = renderPlot(
          plots_RidgePlotRaw[[ii]]
        )
        
        output[[paste0("plot_ridge_norm", ii)]] = renderPlot(
          plots_RidgePlotNorm[[ii]]
        )
        
        output[[paste0("plot_vln_raw", ii)]] = renderPlot(
          plots_ViolinPlotRaw[[ii]]
        )
        
        output[[paste0("plot_vln_norm", ii)]] = renderPlot(
          plots_ViolinPlotNorm[[ii]]
        )
        
        output$plot_dotplot = renderPlot(
          plots_DotPlot,
          width = input$x_axis,
          height = input$y_axis
        )
        
        output$plot_heatmap = renderPlot(
          plots_Heatmap,
          width = input$x_axis,
          height = input$y_axis
        )
        
      })
    }
  })

  # Downloads =================================
  
  output$download_plots = downloadHandler(
    filename = function() {
      paste0(input$archive_download, ".zip")
    },
    content = function(file) {
      # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
      # owd = setwd(tempdir())
      # on.exit(setwd(owd))
      
      # Idea to check if all devices are off on exit
      # on.exit({
      #  setwd(owd)
      #  while length(dev.list())>1; do
      #    dev.off()
      # })

      files = NULL
      on.exit(unlink(files))
      
      if(length(input$genes) == 0){
        shinyalert(
          title = "Download failed!",
          text = "Failed to download plots! Please make sure to upload .rds file and select genes before downloading!",
          size = "s", 
          type = "error"
        )
        return(NULL)
      }else{
        shinyalert(
          title = "Please wait!",
          text = "The creation of files (.png and .pdf) can take a while. The download will start, once everything is ready.",
          size = "s",
          type = "info",
          showConfirmButton = FALSE
        )
      }
      
      ### Download FeaturPlot ############################# 
      if (input$check_featureplot == TRUE) {
        # PNG
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("FeaturePlot_", input$genes[[i]], ".png")
          
          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(
            plots_FeaturePlot[[i]]
          )
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF
        fileName_pdf = "FeaturePlot.pdf"
        pdf(file = fileName_pdf,
            width = (input$x_axis/96),
            height = (input$y_axis/96))
        for (i in 1:length(input$genes)) {
          print(
            plots_FeaturePlot[[i]]
          )
        }
        dev.off()
        files = c(fileName_pdf, files)
      }
      
      ### Download RidgePlot Raw ############################# 
      if (input$check_ridgeplot_raw == TRUE) {
        # PNG
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("RidgePlot_Raw_", input$genes[[i]], ".png")
          
          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(
            plots_RidgePlotRaw[[i]]
          )
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF
        fileName_pdf = "RidgePlot_Raw.pdf"
        pdf(file = fileName_pdf,
            width = (input$x_axis/96),
            height = (input$y_axis/96))
        for (i in 1:length(input$genes)) {
          print(
            plots_RidgePlotRaw[[i]]
          )
        }
        dev.off()
        files = c(fileName_pdf, files)
      }

      ### Download RidgePlot Normalised #############################
      if (input$check_ridgeplot_norm == TRUE) {
        # PNG
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("RidgePlot_Norm_", input$genes[[i]], ".png")
          
          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(
            plots_RidgePlotNorm[[i]]
          )
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF
        fileName_pdf = "RidgePlot_Norm.pdf"
        pdf(file = fileName_pdf,
            width = (input$x_axis/96),
            height = (input$y_axis/96))
        for (i in 1:length(input$genes)) {
          print(
            plots_RidgePlotNorm[[i]]
          )
        }
        dev.off()
        files = c(fileName_pdf, files)
      }

      ### Download ViolinPlot Raw ############################# 
      if (input$check_vlnplot_raw == TRUE) {
        # PNG ViolinPlot Raw
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("ViolinPlot_Raw_", input$genes[[i]], ".png")
          
          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(
            plots_ViolinPlotRaw[[i]]
          )
          dev.off()
          files = c(fileName_png, files)
        }
        
        # PDF ViolinPlot Raw
        fileName_pdf = "ViolinPlot_Raw.pdf"
        pdf(file = fileName_pdf,
            width = (input$x_axis/96),
            height = (input$y_axis/96))
        for (i in 1:length(input$genes)) {
          print(
            plots_ViolinPlotRaw[[i]]
          )
        }
        dev.off()
        files = c(fileName_pdf, files)
      }

      ### Download ViolinPlot Normalised ############################# 
      if (input$check_vlnplot_norm == TRUE) {
        # PNG ViolinPlot Normalised
        for (i in 1:length(input$genes)) {
          fileName_png = paste0("ViolinPlot_Norm_", input$genes[[i]], ".png")
          
          png(fileName_png,
              width = input$x_axis,
              height = input$y_axis)
          print(
            plots_ViolinPlotNorm[[i]]
          )
          dev.off()
          files = c(fileName_png, files)
        }
        # PDF ViolinPlot Normalised
        fileName_pdf = "ViolinPlot_Norm.pdf"
        pdf(file = fileName_pdf,
            width = (input$x_axis/96),
            height = (input$y_axis/96))
        for (i in 1:length(input$genes)) {
          print(
            plots_ViolinPlotNorm[[i]]
          )
        }
        dev.off()
        files = c(fileName_pdf, files)
      }
      
      ### Download DotPlot ############################# 
      if (input$check_dotplot == TRUE) {
        # PNG
        fileName_png = "DotPlot.png"
        
        png(fileName_png,
            width = input$x_axis,
            height = input$y_axis)
        print(
          plots_DotPlot
        )
        dev.off()
        files = c(fileName_png, files)
        
        # PDF
        fileName_pdf = "DotPlot.pdf"
        pdf(file = fileName_pdf,
            width = (input$x_axis/96),
            height = (input$y_axis/96))
        print(
          plots_DotPlot
        )
        dev.off()
        files = c(fileName_pdf, files)
      }
      #print(files)
      # Create zip file for Download, uses array of files
      zip(zipfile = file, files =  files)
      shinyjs::runjs("swal.close();")
    },
    contentType = "application/zip"
  )
}

# Run the application
shinyApp(ui = ui, server = server)
