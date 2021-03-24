tabItem(
  tabName = "download",
  fluidRow(
    column(
      width = 6,

      box(
        width = 12,
        title = "Select Plots to download",
        solidHeader = TRUE,
        status = "primary",
        id = "download_Plotselection",

        fluidRow(
          column(
            width = 6,
            checkboxInput("check_featureplot", "FeaturePlot", value = TRUE),
            checkboxInput("check_ridgeplot_raw", "RidgePlot (Raw)", value = TRUE),
            checkboxInput("check_ridgeplot_norm", "RidgePlot (Norm)", value = TRUE),
            checkboxInput("check_heatmap", "Heatmap", value = TRUE)
          ),

          column(
            width = 6,
            ofset = 3,
            checkboxInput("check_vlnplot_raw", "ViolinPlot (Raw)", value = TRUE),
            checkboxInput("check_vlnplot_norm", "ViolinPlot (Norm)", value = TRUE),
            checkboxInput("check_dotplot", "DotPlot", value = TRUE)
          ),
        ),
      ),#box

      box(
        width = 12,
        title = "Enter name of archive (.zip):",
        solidHeader = TRUE,
        status = "primary",

        textInput(
          inputId = "archive_download",
          label = "",
          value = paste0("scrnaseq_download", "_", Sys.Date())
        ),
      ),#box

      box(
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        downloadButton("download_plots", class = "btn btn-success btn-lg btn-block")
      )#box
    )#column
  ),#fluidRow
)#tabItem
