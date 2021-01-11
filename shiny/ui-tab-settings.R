tabItem(
  tabName = "settings",
  fluidRow(
    column(
      width = 8,

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
      
      box(
        width = 12,
        title = "Change size of plots (in pixel)",
        solidHeader = TRUE,
        status = "primary",
        id = "changepixel",

        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = "x_axis",
              label = "X-Axis (px):",
              value = 1024,
              min = 1,
              max = 3000
            )
          ),

          column(
            width = 6,
            ofset = 3,
            numericInput(
              inputId = "y_axis",
              label = "Y-Axis (px):",
              value = 576,
              min = 1,
              max = 3000
            )
          )
        ),

        actionButton(inputId = "ratio_512",
                     label = "512x288 px",
                     class = "btn btn-outline-secondary"),
        actionButton(inputId = "ratio_1024",
                     label = "1024x576 px",
                     class = "btn btn-outline-secondary")
      ),#box,

      box(
        width = 12,
        title = "Change resolution of plots",
        solidHeader = TRUE,
        status = "primary",
        id = "uploadbox",

        numericInput(
          inputId = "res",
          label = "Resolution of plot, in pixels per inch:",
          value = 144,
          min = 1,
          max = 3000
        ),

        actionButton(inputId = "resolution_70",
                     label = "70",
                     class = "btn btn-outline-primary"),

        actionButton(inputId = "resolution_144",
                     label = "144",
                     class = "btn btn-outline-primary"),
      ),#box

      box(
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        id = "render_plots",
        actionButton(inputId = "renderPlots",
                     label = "Render Plots!",
                     class = "btn btn-success btn-lg btn-block")
      ),#box
    )#column
  )#fluidRow
)#tabItem
