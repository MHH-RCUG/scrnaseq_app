# Author: Marius Rüve

#Load packages
library(shiny)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(markdown)

# Set max upload size to 300 MB
options(shiny.maxRequestSize=300*1024^2)

# Parameter list for convenience
param = list()
param$col = "palevioletred"

# Load data if upload is broken
#sc = readRDS(file = "pbmc_2020-06-09.rds")
#features_names_ids = paste(rownames(sc[["RNA"]][[]]), "_", sc[["RNA"]][[]][,1], sep = "")

# Empty object so Input UI gets rendered
features_names_ids = NULL



# Define UI for application that draws a histogram
ui = fluidPage(
  theme = shinytheme("paper"),
  img(src = "MHH.png", align = "right"),
  titlePanel("scrnaseq_app"),
  
  sidebarPanel(width = 3,
    h5("Plots:"),
    fileInput("rds_file", "Choose Seurat file:", accept = ".rds", buttonLabel = "Browse..."), # File upload
    selectInput("genes", "Genes:", features_names_ids, multiple = TRUE), # Select genes 
    
    fluidRow(column(6, numericInput("x_axis", "X-Axis (px):", value = 1024, min = 1, max = 3000)),
             column(6, ofset = 3, numericInput("y_axis", "Y-Axis (px):", value = 576, min = 1, max = 3000))
    ),
    
    actionButton("restore_axes", "Default settings for axes"), # Restore default values of axes (px)
    
    # br() extend spacing between elements
    br(),
    br(),
    
    h5("Download:"),
    textInput('archive_download',"Enter name of archive (.zip):", value = paste0("Download", "_", Sys.Date()), placeholder = paste0("Download", "_", Sys.Date())),
    fluidRow(
      column(4, checkboxInput("check_featureplot", "FeaturePlot", value = TRUE)),
      column(4, checkboxInput("check_ridgeplot", "RidgePlot", value = TRUE)),
      column(4, checkboxInput("check_dotplot", "DotPlot", value = TRUE)),
      column(4, checkboxInput("check_vlnplot", "ViolinPlot", value = TRUE))
    ),

    downloadButton('download_plots')
    # conditionalPanel( condition = "input.rds_file$datapath",
    #                   )
  ),
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("FeaturePlot", uiOutput('ui_feature')),
                tabPanel("RidgePlot", uiOutput('ui_ridge')),
                tabPanel("DotPlot", plotOutput('plot_dotplot')),
                tabPanel("ViolinPlot", plotOutput('ui_vlnplot')),
                tabPanel("Help", includeMarkdown("include.md"))
    )
  )
)

server = function(input, output, session) {
  
  # File upload
  sc = reactive({
    inFile = input$rds_file
    if (is.null(inFile)) {
      tmp_sc = NULL
    } else {
      tmp_sc = readRDS(inFile$datapath)
    }
    tmp_sc
  })
  
  observe({
    if(!is.null(sc())){
      updateSelectInput(session, "genes","Genes:",paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][,1], sep = ""))
    }
  })
  
  # Button to restore default settings for axes
  observeEvent(input$restore_axes,{
    updateNumericInput(session,"x_axis", "X-Axis (px):", value = 1024, min = 1, max = 3000)
    updateNumericInput(session,"y_axis", "X-Axis (px):", value = 576, min = 1, max = 3000)
  })

  # Feature Plot
  output$ui_feature = renderUI({
      out_feature = list()

        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out_feature[[i]] =  plotOutput(outputId = paste0("plot_feature",i),
                                           width = paste0(input$x_axis, "px"),
                                           height = paste0(input$y_axis, "px"))
        }
        return(out_feature)
    })
    observe({
        for (i in 1:length(input$genes)){
            local({  #because expressions are evaluated at app init
                ii = i
                output[[paste0('plot_feature',ii)]] = renderPlot({
                  return(Seurat::FeaturePlot(sc(), features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], cols=c("lightgrey", param$col), label=TRUE, combine=FALSE))
                })
            })
        }
    })

    #Ridge Plot
    output$ui_ridge = renderUI({
        out_ridge = list()
        
        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out_ridge[[i]] =  plotOutput(outputId = paste0("plot_ridge",i),
                                         width = paste0(input$x_axis, "px"),
                                         height = paste0(input$y_axis, "px"))
        }  
        return(out_ridge) 
    })
    observe({  
        for (i in 1:length(input$genes)){  
            local({  #because expressions are evaluated at app init
                ii = i 
                output[[paste0('plot_ridge',ii)]] = renderPlot({ 
                    return(Seurat::RidgePlot(sc(), features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], combine=FALSE))
                })
            })
        }                                  
    })
    
    # DotPlot
    output$plot_dotplot = renderPlot({
      Seurat::DotPlot(sc(), features=unlist(strsplit(input$genes,"_"))[c(T,F)], cols=c("lightgrey", param$col))
    })
    
    # VlnPlot
    # output$ui_vlnplot = renderUI({
    #   out_vln = list()
    # 
    #   if (length(input$genes)==0){return(NULL)}
    #   for (i in 1:length(input$genes)){
    #     out_vln[[i]] =  plotOutput(outputId = paste0("plot_vln",i),
    #                                  width = paste0(input$x_axis, "px"),
    #                                  height = paste0(input$y_axis, "px"))
    #   }
    #   return(out_vln)
    # })
    # observe({
    #   for (i in 1:length(input$genes)){
    #     local({  #because expressions are evaluated at app init
    #       ii = i
    #       output[[paste0('plot_vln',ii)]] = renderPlot({
    #         return(Seurat::VlnPlot(sc, features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], combine=FALSE))
    #       })
    #     })
    #   }
    # })
    
    # Download
    output$download_plots = downloadHandler(
        filename = function() {
          paste0(input$archive_download, ".zip")
        },
        content = function(file) {
          
          # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
          owd = setwd(tempdir())
          on.exit(setwd(owd))
          
          # Idea to check if all devices are off on exit
          # on.exit({
          #  setwd(owd)
          #  while length(dev.list())>1; do
          #    dev.off()
          # })
        
          files = NULL;
          
          # Download FeaturPlot
          if (input$check_featureplot == TRUE) {
            # PNG
            for (i in 1:length(input$genes)){
              fileName_png = paste0("FeaturePlot_",input$genes[[i]],".png")
              
              png(fileName_png,width = input$x_axis, height = input$y_axis)
              print(Seurat::FeaturePlot(sc(), features=unlist(strsplit(input$genes[[i]],"_"))[c(T,F)], cols=c("lightgrey", param$col), label=TRUE, combine=FALSE))
              dev.off()
              
              files = c(fileName_png,files)
            }
            # PDF
            fileName_pdf = "FeaturePlots.pdf"
            pdf(file = fileName_pdf, width = 16, height = 9)
            for(i in 1:length(input$genes)){
              print(Seurat::FeaturePlot(sc(), features=unlist(strsplit(input$genes[[i]],"_"))[c(T,F)], cols=c("lightgrey", param$col), label=TRUE, combine=FALSE))
            }
            dev.off()
            files = c(fileName_pdf,files)
          }
          
          # Download RidgePlot
          if (input$check_ridgeplot == TRUE) {
            # PNG
            for (i in 1:length(input$genes)){
              fileName_png = paste0("RidgePlot_",input$genes[[i]],".png")

              png(fileName_png,width = input$x_axis, height = input$y_axis)
              print(Seurat::RidgePlot(sc(), features=unlist(strsplit(input$genes[[i]],"_"))[c(T,F)], combine=FALSE))
              dev.off()

              files = c(fileName_png,files)
            }
            # PDF
            fileName_pdf = "RidgePlot.pdf"
            pdf(file = fileName_pdf, width = 16, height = 9)
            for(i in 1:length(input$genes)){
              print(Seurat::RidgePlot(sc(), features=unlist(strsplit(input$genes[[i]],"_"))[c(T,F)], combine=FALSE))
            }
            dev.off()
            files = c(fileName_pdf,files)
          }
          
          # Download DotPlot
          if (input$check_dotplot == TRUE) {
            # PNG
            fileName_png = "DotPlot.png"

            png(fileName_png,width = input$x_axis, height = input$y_axis)
            print(Seurat::DotPlot(sc(), features=unlist(strsplit(input$genes,"_"))[c(T,F)], cols=c("lightgrey", param$col)))
            dev.off()

            files = c(fileName_png,files)
            
            # PDF
            fileName_pdf = "DotPlot.pdf"
            pdf(file = fileName_pdf, width = 16, height = 9)
            print(Seurat::DotPlot(sc(), features=unlist(strsplit(input$genes,"_"))[c(T,F)], cols=c("lightgrey", param$col)))
            dev.off()
            files = c(fileName_pdf,files)
          }
          # Create zip file for Download, uses array of files
          zip(file,files)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
