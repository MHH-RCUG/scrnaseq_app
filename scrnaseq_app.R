# Author: Marius Rueve

#Load packages
library(shiny)
library(Seurat)
library(ggplot2)
library(shinythemes)
library(markdown)
library(readxl)

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



########################################
# UI
########################################
ui = fluidPage(
  theme = shinytheme("paper"),
  img(src = "MHH.png", align = "right"),
  titlePanel("scrnaseq_app"),
  
  sidebarPanel(width = 3,
    h5("Plots:"),
    fileInput("rds_file", "Select Seurat file:", accept = ".rds", buttonLabel = "Browse..."), # File upload
    
    fileInput("xlsx_file", "Select Excel file:", accept = ".xlsx", buttonLabel = "Browse..."),
    checkboxInput("header", "Check if column names are given: (Header)", TRUE),
    selectInput("genes", "Genes:", features_names_ids, multiple = TRUE), # Select genes

    fluidRow(column(6, numericInput("x_axis", "X-Axis (px):", value = 1024, min = 1, max = 3000)),
             column(6, ofset = 3, numericInput("y_axis", "Y-Axis (px):", value = 576, min = 1, max = 3000))
    ),

    actionButton("restore_axes", "Default settings for axes"), # Restore default values of axes (px)

    # br() extend spacing between elements
    br(),
    br(),

    h5("Download:"),
    textInput("archive_download","Enter name of archive (.zip):", value = paste0("Download", "_", Sys.Date()), placeholder = paste0("Download", "_", Sys.Date())),
    fluidRow(
      column(5, checkboxInput("check_featureplot", "FeaturePlot", value = TRUE)),
      column(5, checkboxInput("check_ridgeplot_raw", "RidgePlot (Raw)", value = TRUE)),
      column(5, checkboxInput("check_ridgeplot_norm", "RidgePlot (Norm)", value = TRUE)),
      column(5, checkboxInput("check_vlnplot_raw", "ViolinPlot (Raw)", value = TRUE)),
      column(5, checkboxInput("check_vlnplot_norm", "ViolinPlot (Norm)", value = TRUE)),
      column(5, checkboxInput("check_dotplot", "DotPlot", value = TRUE))
    ),

    downloadButton("download_plots")
    # conditionalPanel( condition = "input.rds_file$datapath",
    #                   )
  ),
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("FeaturePlot", uiOutput("ui_feature")),
                tabPanel("RidgePlot (Raw)", uiOutput("ui_ridge_raw")),
                tabPanel("RidgePlot (Norm)", uiOutput("ui_ridge_norm")),
                tabPanel("ViolinPlot (Raw)", uiOutput("ui_vln_raw")),
                tabPanel("ViolinPlot (Norm)", uiOutput("ui_vln_norm")),
                tabPanel("DotPlot", plotOutput("plot_dotplot")),
                tabPanel("Help", includeMarkdown("include.md"))
    )
  )
)



########################################
# SERVER
########################################
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
      features_names_ids = paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][,1], sep = "")
      #print(features_names_ids)
      updateSelectInput(session, "genes","Genes:",features_names_ids)
    }
  })
  
  # Excel file input
  excel_genes = reactive({
    inFile = input$xlsx_file
    if (is.null(inFile)) {
      tmp = NULL
    } else {
      tmp = read_excel(path = inFile$datapath,
                       sheet = 1,
                       col_names = input$header)
    }
    tmp
  })
  
  observe({
    if(!is.null(sc()) & !is.null(excel_genes())){
      print(excel_genes())
      
      excel_ensemblids = excel_genes()[,1]
      
      print(excel_ensemblids)
      
      print(grepl(features_names_ids[1,], excel_ensemblids()[1,]))
    }
  })
  
  # observe({
  # 
  #   req(input$xlsx_file)
  # 
  #   inFile <- input$xlsx_file
  # 
  #   df_excel_file = read_excel(path = inFile$datapath,
  #                              sheet = 1,
  #                              col_names = input$header)
  # 
  #   excel_ensemblids = as.vector(df_excel_file[,1])
  #   
  #   #print(excel_ensemblids)
  #   
  #   x = vector()
  #   
  #   print(excel_ensemblids[[1]])
  #   
  #   #print(grepl(features_names_ids[1,], excel_ensemblids[1,]))
  # 
  #   for(i in excel_ensemblids){
  #     for(ii in features_names_ids){
  #       if(grep(i, ii, ignore.case = TRUE)){
  #         print("grepl")
  #       }
  #       print(i)
  #       print(ii)
  #     }
  #   }
  # 
  #   updateSelectInput(session = session,
  #                     inputId = "genes",
  #                     label = "Genes:",
  #                     choices = features_names_ids,
  #                     selected = x)
  # 
  # })
  
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
                output[[paste0("plot_feature",ii)]] = renderPlot({
                  return(Seurat::FeaturePlot(sc(), features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], cols=c("lightgrey", param$col), label=TRUE, combine=FALSE))
                })
            })
        }
    })

    #Ridge Plot Raw
    output$ui_ridge_raw = renderUI({
        out_ridge_raw = list()
        
        if (length(input$genes)==0){return(NULL)}
        for (i in 1:length(input$genes)){
            out_ridge_raw[[i]] =  plotOutput(outputId = paste0("plot_ridge_raw",i),
                                         width = paste0(input$x_axis, "px"),
                                         height = paste0(input$y_axis, "px"))
        }  
        return(out_ridge_raw) 
    })
    observe({  
        for (i in 1:length(input$genes)){  
            local({  #because expressions are evaluated at app init
                ii = i 
                output[[paste0("plot_ridge_raw",ii)]] = renderPlot({ 
                    return(Seurat::RidgePlot(sc(), 
                                             features = unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], 
                                             assay = "RNA", 
                                             slot = "counts", 
                                             combine=FALSE))
                })
            })
        }                                  
    })
    
    #Ridge Plot Normalised
    output$ui_ridge_norm = renderUI({
      out_ridge_norm = list()
      
      if (length(input$genes)==0){return(NULL)}
      for (i in 1:length(input$genes)){
        out_ridge_norm[[i]] =  plotOutput(outputId = paste0("plot_ridge_norm",i),
                                     width = paste0(input$x_axis, "px"),
                                     height = paste0(input$y_axis, "px"))
      }  
      return(out_ridge_norm) 
    })
    observe({  
      for (i in 1:length(input$genes)){  
        local({  #because expressions are evaluated at app init
          ii = i 
          output[[paste0("plot_ridge_norm",ii)]] = renderPlot({ 
            return(Seurat::RidgePlot(sc(), 
                                     features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)],
                                     slot = "data",
                                     combine=FALSE))
          })
        })
      }                                  
    })
    
    #Violin plot of raw gene expression counts
    output$ui_vln_raw = renderUI({
      out_vln_raw = list()
      
      if (length(input$genes)==0){return(NULL)}
      for (i in 1:length(input$genes)){
        out_vln_raw[[i]] =  plotOutput(outputId = paste0("plot_vln_raw",i),
                                          width = paste0(input$x_axis, "px"),
                                          height = paste0(input$y_axis, "px"))
      }  
      return(out_vln_raw) 
    })
    observe({  
      for (i in 1:length(input$genes)){  
        local({  #because expressions are evaluated at app init
          ii = i 
          output[[paste0("plot_vln_raw",ii)]] = renderPlot({ 
            return(Seurat::VlnPlot(sc(),
                                   features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)],
                                   slot = "data",
                                   combine=FALSE,
                                   pt.size = 0.2
                                   )
                   )
          })
        })
      }                                  
    })
    
    #Violin plot of normalised gene expression data
    output$ui_vln_norm = renderUI({
      out_vln_norm = list()
      
      if (length(input$genes)==0){return(NULL)}
      for (i in 1:length(input$genes)){
        out_vln_norm[[i]] =  plotOutput(outputId = paste0("plot_vln_norm",i),
                                       width = paste0(input$x_axis, "px"),
                                       height = paste0(input$y_axis, "px"))
        }
      return(out_vln_norm)
      })
    observe({  
      for (i in 1:length(input$genes)){  
        local({  #because expressions are evaluated at app init
          ii = i 
          output[[paste0("plot_vln_norm",ii)]] = renderPlot({ 
            return(Seurat::VlnPlot(sc(),
                                   features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)],
                                   combine=FALSE,
                                   pt.size = 0.2
                                   )
                   )
            })
          })
        }
      })
    
    # DotPlot
    output$plot_dotplot = renderPlot({
      Seurat::DotPlot(sc(), features=unlist(strsplit(input$genes,"_"))[c(T,F)], cols=c("lightgrey", param$col))
    })
    
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
          
          # Download RidgePlot Raw
          if (input$check_ridgeplot_raw == TRUE) {
            # PNG
            for (i in 1:length(input$genes)){
              fileName_png = paste0("RidgePlot_Raw_",input$genes[[i]],".png")

              png(fileName_png,width = input$x_axis, height = input$y_axis)
              print(Seurat::RidgePlot(sc(), 
                                      features = unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], 
                                      assay = "RNA", 
                                      slot = "counts", 
                                      combine=FALSE))
              dev.off()

              files = c(fileName_png,files)
            }
            # PDF
            fileName_pdf = "RidgePlot_Raw.pdf"
            pdf(file = fileName_pdf, width = 16, height = 9)
            for(i in 1:length(input$genes)){
              print(Seurat::RidgePlot(sc(), 
                                      features = unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], 
                                      assay = "RNA", 
                                      slot = "counts", 
                                      combine=FALSE))
            }
            dev.off()
            files = c(fileName_pdf,files)
          }
          
          # Download RidgePlot Normalised
          if (input$check_ridgeplot_norm == TRUE) {
            # PNG
            for (i in 1:length(input$genes)){
              fileName_png = paste0("RidgePlot_Norm_",input$genes[[i]],".png")
              
              png(fileName_png,width = input$x_axis, height = input$y_axis)
              print(Seurat::RidgePlot(sc(), 
                                      features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)],
                                      slot = "data",
                                      combine=FALSE))
              dev.off()
              
              files = c(fileName_png,files)
            }
            # PDF
            fileName_pdf = "RidgePlot_Norm.pdf"
            pdf(file = fileName_pdf, width = 16, height = 9)
            for(i in 1:length(input$genes)){
              print(Seurat::RidgePlot(sc(), 
                                      features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)],
                                      slot = "data",
                                      combine=FALSE))
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
