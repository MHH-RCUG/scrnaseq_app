tabItem(
  tabName = "settings_general",
  fluidRow(
    column(
      width = 6,
      box(
        width = 12,
        title = "Change Size (in pixel) and Resolution (in dpi) of Plots",
        solidHeader = TRUE,
        status = "primary",
        id = "changepixel",

        fluidRow(
          column(
            width = 6,
            numericInput(
              inputId = "x_axis",
              label = "X-Axis (px):",
              value = 512,
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
              value = 288,
              min = 1,
              max = 3000
            )
          )
        ),
        
        numericInput(
          inputId = "res",
          label = "Resolution of plot, in pixels per inch:",
          value = 70,
          min = 1,
          max = 3000
        )
      ),#box,

      box(
        width = 12,
        title = "Presets",
        solidHeader = TRUE,
        status = "primary",
        
        actionButton(inputId = "preset_512x288_70",
                     label = "512 x 288; 70 dpi",
                     class = "btn btn-outline-secondary"),
        actionButton(inputId = "preset_768x432_110",
                     label = "768 x 432; 110 dpi",
                     class = "btn btn-outline-secondary"),
        actionButton(inputId = "preset_1280x720_144",
                     label = "1280 x 720; 144 dpi",
                     class = "btn btn-outline-secondary"),
        actionButton(inputId = "preset_1920x1080_300",
                     label = "1920 x 1080; 300 dpi",
                     class = "btn btn-outline-secondary"),
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
