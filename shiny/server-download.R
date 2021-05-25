# Wirt functions for png tiff and pdf to reduce repi
output$download_plots <- downloadHandler(
  filename = function() {
    # Takes the filename from the input of the user
    paste0(input$archive_download, ".zip")
  },
  content = function(file) {
    # temporarily switch to the temp dir, in case you do not have write permission to the current working directory
    owd <- setwd(tempdir())
    on.exit(setwd(owd))
    # on.exit(unlink(files))

    files <- NULL

    if (length(input$genes) == 0) {
      # Message if download fails due to no genese selected
      shinyalert(
        title = "Download failed!",
        text = "Failed to download plots! Please make sure to upload .rds file and select genes before downloading!",
        size = "s",
        type = "error"
      )
      return(NULL)
    } else {
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
    if (input$check_featureplot) {
      for (i in 1:length(input$genes)) {
        # PNG FeaturPlot
        if (input$check_png) {
          fileName_png <- paste0("FeaturePlot_", input$genes[[i]], ".png")
          png(
            filename = fileName_png,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_FeaturePlots[[i]])
          dev.off()
          files <- c(fileName_png, files)
        }

        # TIFF FeaturPlot
        if (input$check_tiff) {
          fileName_tiff <- paste0("FeaturePlot_", input$genes[[i]], ".tiff")
          tiff(
            filename = fileName_tiff,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_FeaturePlots[[i]])
          dev.off()
          files <- c(fileName_tiff, files)
        }
      }

      # PDF FeaturPlot
      if (input$check_pdf) {
        fileName_pdf <- "FeaturePlot.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        for (i in 1:length(input$genes)) {
          print(stored_FeaturePlots[[i]])
        }
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }


    # Download RidgePlot Raw #############################
    if (input$check_ridgeplot_raw) {
      for (i in 1:length(input$genes)) {
        # PNG RidgePlot Raw
        if (input$check_png) {
          fileName_png <- paste0("RidgePlot_Raw_", input$genes[[i]], ".png")
          png(fileName_png,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_RidgePlotRaws[[i]])
          dev.off()
          files <- c(fileName_png, files)
        }

        # TIFF RidgePlot Raw
        if (input$check_tiff) {
          fileName_tiff <- paste0("RidgePlot_Raw_", input$genes[[i]], ".tiff")
          tiff(
            filename = fileName_tiff,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_RidgePlotRaws[[i]])
          dev.off()
          files <- c(fileName_tiff, files)
        }
      }
      # PDF
      if (input$check_pdf) {
        fileName_pdf <- "RidgePlot_Raw.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        for (i in 1:length(input$genes)) {
          print(stored_RidgePlotRaws[[i]])
        }
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }


    # Download RidgePlot Normalised #############################
    if (input$check_ridgeplot_norm) {
      for (i in 1:length(input$genes)) {
        # PNG RidgePlot Norm
        if (input$check_png) {
          fileName_png <- paste0("RidgePlot_Norm_", input$genes[[i]], ".png")
          png(
            filename = fileName_png,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_RidgePlotNorms[[i]])
          dev.off()
          files <- c(fileName_png, files)
        }

        # TIFF RidgePlot Norm
        if (input$check_tiff) {
          fileName_tiff <- paste0("RidgePlot_Norm_", input$genes[[i]], ".tiff")
          tiff(
            filename = fileName_tiff,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_RidgePlotNorms[[i]])
          dev.off()
          files <- c(fileName_tiff, files)
        }
      }
      # PDF RidgePlot Norm
      if (input$check_pdf) {
        fileName_pdf <- "RidgePlot_Norm.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        for (i in 1:length(input$genes)) {
          print(stored_RidgePlotNorms[[i]])
        }
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }


    # Download ViolinPlot Raw #############################
    if (input$check_vlnplot_raw) {
      for (i in 1:length(input$genes)) {
        # PNG ViolinPlot Raw
        if (input$check_png) {
          fileName_png <- paste0("ViolinPlot_Raw_", input$genes[[i]], ".png")
          png(
            filename = fileName_png,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_ViolinPlotRaws[[i]])
          dev.off()
          files <- c(fileName_png, files)
        }

        # TIFF ViolinPlot Raw
        if (input$check_tiff) {
          fileName_tiff <- paste0("ViolinPlot_Raw_", input$genes[[i]], ".tiff")
          tiff(
            filename = fileName_tiff,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_ViolinPlotRaws[[i]])
          dev.off()
          files <- c(fileName_tiff, files)
        }
      }

      # PDF ViolinPlot Raw
      if (input$check_pdf) {
        fileName_pdf <- "ViolinPlot_Raw.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        for (i in 1:length(input$genes)) {
          print(stored_ViolinPlotRaws[[i]])
        }
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }

    # Download ViolinPlot Normalised #############################
    if (input$check_vlnplot_norm) {
      for (i in 1:length(input$genes)) {
        # PNG ViolinPlot Normalised
        if (input$check_png) {
          fileName_png <- paste0("ViolinPlot_Norm_", input$genes[[i]], ".png")
          png(
            filename = fileName_png,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_ViolinPlotNorms[[i]])
          dev.off()
          files <- c(fileName_png, files)
        }

        # TIFF ViolinPlot Normalised
        if (input$check_tiff) {
          fileName_tiff <- paste0("ViolinPlot_Norm_", input$genes[[i]], ".tiff")
          tiff(
            filename = fileName_tiff,
            width = input$x_axis,
            height = input$y_axis,
            res = input$res
          )
          print(stored_ViolinPlotNorms[[i]])
          dev.off()
          files <- c(fileName_tiff, files)
        }
      }

      # PDF ViolinPlot Normalised
      if (input$check_pdf) {
        fileName_pdf <- "ViolinPlot_Norm.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        for (i in 1:length(input$genes)) {
          print(stored_ViolinPlotNorms[[i]])
        }
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }

    # Download DotPlot #############################
    if (input$check_dotplot) {
      # PNG DotPlot
      if (input$check_png) {
        fileName_png <- "DotPlot.png"
        png(
          filename = fileName_png,
          width = input$x_axis,
          height = input$y_axis,
          res = input$res
        )
        print(stored_DotPlot)
        dev.off()
        files <- c(fileName_png, files)
      }

      # TIFF DotPlot
      if (input$check_tiff) {
        fileName_tiff <- "DotPlot.tiff"
        tiff(
          filename = fileName_tiff,
          width = input$x_axis,
          height = input$y_axis,
          res = input$res
        )
        print(stored_DotPlot)
        dev.off()
        files <- c(fileName_tiff, files)
      }

      # PDF
      if (input$check_pdf) {
        fileName_pdf <- "DotPlot.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        print(stored_DotPlot)
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }


    # Download Heatmap #############################
    if (input$check_heatmap == TRUE) {
      # PNG Heatmap
      if (input$check_png) {
        fileName_png <- "Heatmap.png"
        png(fileName_png,
          width = input$x_axis,
          height = input$y_axis,
          res = input$res
        )
        print(stored_Heatmap)
        dev.off()
        files <- c(fileName_png, files)
      }

      # TIFF Heatmap
      if (input$check_tiff) {
        fileName_tiff <- "Heatmap.tiff"
        tiff(
          filename = fileName_tiff,
          width = input$x_axis,
          height = input$y_axis,
          res = input$res
        )
        print(stored_Heatmap)
        dev.off()
        files <- c(fileName_tiff, files)
      }

      # PDF Heatmap
      if (input$check_pdf) {
        fileName_pdf <- "Heatmap.pdf"
        pdf(
          file = fileName_pdf,
          width = (input$x_axis / input$res),
          height = (input$y_axis / input$res)
        )
        print(stored_Heatmap)
        dev.off()
        files <- c(fileName_pdf, files)
      }
    }

    # Create zip file for Download, uses array of files
    shinyjs::runjs("swal.close();")
    zip(zipfile = file, files = files)
  },
  contentType = "application/zip"
)
