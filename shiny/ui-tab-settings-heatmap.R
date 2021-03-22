tabItem(
  tabName = "settings_heatmap",
  
  fluidRow(
    column(
      width = 6,
      box(
        width = 12,
        title = "Heatmap settings",
        solidHeader = TRUE,
        status = "primary",
        id = "heatmap_settings",
        
        selectInput(
          inputId = "heatmap_slot",
          label = "Select data slot to use",
          choices = c("counts", "data"),
          selected = "counts"
        ),
        
        selectInput(
          inputId = "heatmap_assay",
          label = "Assay to pull from",
          choices = c("RNA", "SCT"),
          selected = "RNA"
        )
        
      ),#box_heatmap
    )#column
  )#fluidRow
)#tabItem