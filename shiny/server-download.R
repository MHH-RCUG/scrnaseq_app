output$download_plots = downloadHandler(
  filename = function() {
    # Takes the filename from the input of the user
    paste0(input$archive_download, ".zip")
  },
  content = function(file) {
    #temporarily switch to the temp dir, in case you do not have write permission to the current working directory
    owd = setwd(tempdir())
    on.exit(setwd(owd))
    #on.exit(unlink(files))

    files = NULL

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

    # Download Heatmap #############################
    if (input$check_heatmap == TRUE) {
      # PNG
      fileName_png = "Heatmap.png"

      png(fileName_png,
          width = input$x_axis,
          height = input$y_axis)
      print(stored_Heatmap)
      dev.off()
      files = c(fileName_png, files)
      # PDF
      fileName_pdf = "Heatmap.pdf"
      pdf(
        file = fileName_pdf,
        width = (input$x_axis / 96),
        height = (input$y_axis / 96)
      )
      print(stored_Heatmap)
      dev.off()
      files = c(fileName_pdf, files)
    }

    # Create zip file for Download, uses array of files
    shinyjs::runjs("swal.close();")
    zip(zipfile = file, files =  files)
  },
  contentType = "application/zip"
)
