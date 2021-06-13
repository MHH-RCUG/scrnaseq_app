tabItem(
  tabName = "settings_colors",
  
  fluidRow(
    column(
      width = 6,
      box(
        width = 12,
        title = "Select colors",
        solidHeader = T,
        status = "primary",
        id = "selectcolors",
        
        selectizeInput(
          inputId = "colors",
          label = "Select colors for plots:",
          choices = c(1,2,3),
          multiple = TRUE
        ),#selectInput
      ),#box
      
      box(
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        id = "render_plots",
        actionButton(inputId = "renderPlots",
                     label = "Render Plots!",
                     class = "btn btn-success btn-lg btn-block")
      )#box
    )
  )
)