#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(Seurat)
library(ggplot2)

options(shiny.maxRequestSize=300*1024^2)
#source("FeaturePlot.R")

param = list()
param$col = "palevioletred"

#sc = readRDS(file = "pbmc_2020-06-09.rds")
#sc = c("test1","test2","test3")
features_names_ids = NULL
#features_names_ids = paste(rownames(sc[["RNA"]][[]]), "_", sc[["RNA"]][[]][,1], sep = "")
test = FALSE

# Define UI for application that draws a histogram
ui = fluidPage(
        titlePanel("scrnaseq"),

        sidebarPanel(
            fileInput("rds_file","Choose Seurat file:", accept = ".rds", buttonLabel = "Browse..."),
            selectInput("genes", "Genes:", features_names_ids, multiple = TRUE),
            numericInput("x_axis", "X-Axis (px):", value = 1024, min = 1, max = 3000),
            numericInput("y_axis", "Y-Axis (px):", value = 576, min = 1, max = 3000),
            textInput('archive_download',"Enter name of archive (.zip):", value = paste0("Download", "_", Sys.Date()), placeholder = paste0("Download", "_", Sys.Date())),
            checkboxInput("check_featureplot", "FeaturePlot", value = TRUE),
            checkboxInput("check_ridgeplot", "RidgePlot", value = TRUE),
            checkboxInput("check_dotplot", "DotPlot", value = TRUE),
            checkboxInput("check_vlnplot", "ViolinPlot", value = TRUE),
            downloadButton('download_plots')
            # conditionalPanel( condition = "input.rds_file$datapath",
            #                   )
        ),
        
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("FeaturePlot", uiOutput('ui_feature')),
                        #tabPanel("Feature", plotOutput('plot_feature_test')),
                        tabPanel("RidgePlot", uiOutput('ui_ridge')),
                        tabPanel("DotPlot", plotOutput('plot_dotplot')),
                        tabPanel("ViolinPlot", plotOutput('ui_vlnplot'))
            )
        )
)


server = function(input, output, session) {
  
  sc <- reactive({
    inFile <- input$rds_file
    if (is.null(inFile)) {
      d <- NULL
    } else {
      d <- readRDS(inFile$datapath)
    }
    d
  })
  
  # output$contents <- renderTable({
  #   myData()
  # })
  
  observe({
    if(!is.null(sc())){
      updateSelectInput(session, "genes",
                        label = "Gene(s):",
                        choices = paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][,1], sep = ""))
    }
  })
  
  renderUI({
    
  })
  #sc = c("test1","test2","test3")

  # fileInput = eventReactive(input$rds_file, {
  #   sc = readRDS(input$rds_file$datapath)
  #   print("done")
  # 
  #   features_names_ids = paste(rownames(sc[["RNA"]][[]]), "_", sc[["RNA"]][[]][,1], sep = "")
  #   updateSelectInput(session, "genes", label = features_names_ids, multiple = TRUE)
  #   return(features_names_ids)
  # })
  # 
  # reactive({
  #   features_names_ids = fileInput()
  # })

  ########################################
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
                  return(Seurat::FeaturePlot(sc(), features=unlist(strsplit(input$genes[[ii]],"_"))[c(T,F)], cols=c("lightgrey", param$col), combine=FALSE))
                })
            })
        }
    })

    ########################################
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
    
    ########################################
    # DotPlot
    output$plot_dotplot = renderPlot({
      Seurat::DotPlot(sc(), features=unlist(strsplit(input$genes,"_"))[c(T,F)], cols=c("lightgrey", param$col))
    })
    
    ########################################
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
    
    ########################################
    # Download
    output$download_plots = downloadHandler(
        filename = function() {
          paste0(input$archive_download, ".zip")
        },
        content = function(file) {
          #ggsave(file, plot = out_feature[[1]], device = "png")

          owd <- setwd(tempdir())
          on.exit(setwd(owd))
          
          # on.exit({
          #  setwd(owd)
          #  while length(dev.list())>1; do
          #    dev.off()
          # })
        
          files <- NULL;
          
          if (input$check_featureplot == TRUE) {
            for (i in 1:length(input$genes)){
              #print(paste0(input$archive_download,"_0",i,".png"))
              fileName <- paste0("FeaturePlot_",input$genes[[i]],".png")
              
              png(fileName,width = input$x_axis, height = input$y_axis)
              print(Seurat::FeaturePlot(sc(), features=unlist(strsplit(input$genes[[i]],"_"))[c(T,F)], cols=c("lightgrey", param$col), combine=FALSE))
              dev.off()
              
              files <- c(fileName,files)
            }
          }
          
          if (input$check_ridgeplot == TRUE) {
            for (i in 1:length(input$genes)){
              #print(paste0(input$archive_download,"_0",i,".png"))
              fileName <- paste0("RidgePlot_",input$genes[[i]],".png")

              png(fileName,width = input$x_axis, height = input$y_axis)
              print(Seurat::RidgePlot(sc(), features=unlist(strsplit(input$genes[[i]],"_"))[c(T,F)], combine=FALSE))
              dev.off()

              files <- c(fileName,files)
            }
          }
          
          if (input$check_dotplot == TRUE) {
            #print(paste0(input$archive_download,"_0",i,".png"))
            fileName <- paste0("DotPlot",".png")

            png(fileName,width = input$x_axis, height = input$y_axis)
            print(Seurat::DotPlot(sc(), features=unlist(strsplit(input$genes,"_"))[c(T,F)], cols=c("lightgrey", param$col)))
            dev.off()

            files <- c(fileName,files)
          }
          
          zip(file,files)
        }
    )
    
    ########################################
    # Test
    # output$plot_feature_test = renderPlot({
    #   p = Seurat::FeaturePlot(sc, features=input$genes, cols=c("lightgrey", param$col), combine=FALSE)
    #   names(p) = input$genes
    #   for (i in input$genes) p[[i]] = PlotMystyle(p[[i]], title=i)
    #   patchwork::wrap_plots(p, ncol=1)
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
