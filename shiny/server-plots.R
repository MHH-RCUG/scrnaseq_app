# renderPlots() --> plotOutput() --> renderUI()

render_plots = function(){
  shinyalert(
    title = "Please wait!",
    text = "The creation of plots can take a while.",
    size = "s",
    type = "info",
    showConfirmButton = FALSE
  )

  output$plots = renderMenu({
    menuItem("Plots", tabName = "plots", icon = icon("bar-chart"), startExpanded = TRUE,
             menuSubItem("Feature Plots", tabName = "featureplots"),
             menuSubItem("Ridge Plots (Raw)", tabName = "ridgeplots_raw"),
             menuSubItem("Ridge Plots (Norm)", tabName = "ridgeplots_norm"),
             menuSubItem("Violin Plots (Raw)", tabName = "violinplots_raw"),
             menuSubItem("Violin Plots (Norm)", tabName = "violinplots_norm"),
             menuSubItem("Dot Plot", tabName = "dotplot"),
             menuSubItem("Heatmap", tabName = "heatmap")
    )#menuItem
  })#renderMenu

  output$download = renderMenu({
    menuItem("Download", tabName = "download", icon = icon("download"))
  })#renderMenu

  for (i in 1:length(input$genes)) {
    # Feature Plots
    stored_FeaturePlots[[i]] <<- # <<- for global assignments
      Seurat::FeaturePlot(
        sc(),
        features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
        col=c("lightgrey", param$col),
        label = TRUE
      ) +
      AddStyle(title=input$genes[[i]])

    # Ridge Plots Raw
    stored_RidgePlotRaws[[i]] <<-
      Seurat::RidgePlot(
        sc(),
        features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
        assay = "RNA",
        slot = "counts"
      ) +
      AddStyle(title=input$genes[[i]],
               legend_title="Cluster",
               fill=param$col_clusters,
               ylab="Cluster")

    # Ridge Plots Norm
    stored_RidgePlotNorms[[i]] <<-
      Seurat::RidgePlot(
        sc(),
        features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
        slot = "data"
      ) +
      AddStyle(title=input$genes[[i]],
               legend_title="Cluster",
               fill=param$col_clusters,
               ylab="Cluster")

    # Violin Plots Raw
    stored_ViolinPlotRaws[[i]] <<-
      Seurat::VlnPlot(
        sc(),
        features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
        assay = "RNA",
        slot = "counts",
        pt.size = 0.2
      ) +
      AddStyle(title=input$genes[[i]],
               legend_title="Cluster",
               fill=param$col_clusters,
               xlab="Cluster")

    # Violin Plots Norm
    stored_ViolinPlotNorms[[i]] <<-
      Seurat::VlnPlot(
        object = sc(),
        features = unlist(strsplit(input$genes[[i]], "_"))[c(T, F)],
        pt.size = 0.2
      ) +
      AddStyle(title = input$genes[[i]],
               legend_title = "Cluster",
               fill = param$col_clusters,
               xlab = "Cluster")
  }#for

  # DotPlot
  stored_DotPlot <<-
    Seurat::DotPlot(
      sc(),
      features = unlist(strsplit(input$genes, "_"))[c(T, F)],
      cols = c("lightgrey", param$col)
    ) +
    AddStyle(title = paste0("Test"), ylab="Cluster", legend_position = "bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
    guides(size = guide_legend(order = 1))

  tryCatch({
    # Heatmap
    stored_Heatmap <<-
      Seurat::DoHeatmap(
        sc(),
        features = unlist(strsplit(input$genes, "_"))[c(T, F)],
        group.colors = param$col_clusters,
        slot = input$heatmap_slot,
        assay = input$heatmap_assay
      ) +
      NoLegend() +
      theme(axis.text.y = element_blank())
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

  for (i in 1:length(input$genes)) {
    local({
      # Without it, the value of i in the renderPlot() will be the same across all instances,
      # because of when the expression is evaluated.
      ii = i

      # This chunks takes the stored plots from the global variable and renders a reactive plots that is
      # suitable for assigning to an `output` slot; passing vector to UI
      output[[paste0("plot_feature", ii)]] = renderPlot(stored_FeaturePlots[[ii]],
                                                        res = input$res,
                                                        width = input$x_axis,
                                                        height = input$y_axis)

      output[[paste0("plot_ridge_raw", ii)]] = renderPlot(stored_RidgePlotRaws[[ii]],
                                                          res = input$res,
                                                          width = input$x_axis,
                                                          height = input$y_axis)

      output[[paste0("plot_ridge_norm", ii)]] = renderPlot(stored_RidgePlotNorms[[ii]],
                                                           res = input$res,
                                                           width = input$x_axis,
                                                           height = input$y_axis)

      output[[paste0("plot_vln_raw", ii)]] = renderPlot(stored_ViolinPlotRaws[[ii]],
                                                        res = input$res,
                                                        width = input$x_axis,
                                                        height = input$y_axis)

      output[[paste0("plot_vln_norm", ii)]] = renderPlot(stored_ViolinPlotNorms[[ii]],
                                                         res = input$res,
                                                         width = input$x_axis,
                                                         height = input$y_axis)
    })#local
  }#for

  # UI FeaturePlot #############################
  output$ui_feature = renderUI({
    tmp = list()
    if (length(input$genes) == 0) {
      return(NULL)
    }
    for (i in 1:length(input$genes)) {
      tmp[[i]] =  shinycssloaders::withSpinner(
        plotOutput(
          outputId = paste0("plot_feature", i),
          width = paste0(input$x_axis, "px"),
          height = paste0(input$y_axis, "px")
        )
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
      tmp[[i]] = shinycssloaders::withSpinner(
        plotOutput(
          outputId = paste0("plot_ridge_raw", i),
          width = paste0(input$x_axis, "px"),
          height = paste0(input$y_axis, "px")
        )
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
      tmp[[i]] = shinycssloaders::withSpinner(
        plotOutput(
          outputId = paste0("plot_ridge_norm", i),
          width = paste0(input$x_axis, "px"),
          height = paste0(input$y_axis, "px")
        )
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
      tmp[[i]] = shinycssloaders::withSpinner(
        plotOutput(
          outputId = paste0("plot_vln_raw", i),
          width = paste0(input$x_axis, "px"),
          height = paste0(input$y_axis, "px")
        )
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
      tmp[[i]] = shinycssloaders::withSpinner(
        plotOutput(
          outputId = paste0("plot_vln_norm", i),
          width = paste0(input$x_axis, "px"),
          height = paste0(input$y_axis, "px")
        )
      )
    }
    return(tmp)
  })

  output$plot_dotplot = renderPlot(stored_DotPlot,
                                   width = input$x_axis,
                                   height = input$y_axis,
                                   res = input$res)

  output$plot_heatmap = renderPlot(stored_Heatmap,
                                   width = input$x_axis,
                                   height = input$y_axis,
                                   res = input$res)

  shinyjs::runjs("swal.close();")
  Sys.sleep(0.5)
  updateTabItems(session = session,
                 inputId = "tabs",
                 selected = "featureplots")
}#funcion: render_plots

observeEvent(input$renderPlots, {
  render_plots()
})#observeEvent
