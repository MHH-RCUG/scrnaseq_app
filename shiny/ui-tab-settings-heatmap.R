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
      
      box(
        width = 12,
        title = "Info",
        solidHeader = TRUE,
        status = "primary",

        p(strong("slot: "),"Data slot to use, choose from 'counts', 'data', or 'scale.data'"),
        p(strong("assay: "),"Assay to pull from")
      ),
      
      box(
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        id = "render_plots",
        actionButton(inputId = "renderPlots",
                     label = "Render Plots!",
                     class = "btn btn-success btn-lg btn-block")
      )#box
    )#column
  )#fluidRow
)#tabItem