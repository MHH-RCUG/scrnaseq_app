tabItem(
  tabName = "settings_colors",
  
  fluidRow(
    column(
      width = 8,
      
      box(
        width = 12,
        title = "Select colors",
        solidHeader = T,
        status = "primary",
        id = "selectcolors",
        
        selectInput(
          inputId = "colors",
          label = "Select colors for plots:",
          choices = c(1,2,3),
          multiple = TRUE
        ),#selectInput
      )#box
    )
  )
)