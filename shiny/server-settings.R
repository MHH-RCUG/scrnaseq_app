tabItem(
  tabName = "settings",

  observeEvent(input$ratio_512, {
    updateNumericInput(
      session,
      "x_axis",
      "X-Axis (px):",
      value = 512,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session,
      "y_axis",
      "Y-Axis (px):",
      value = 288,
      min = 1,
      max = 3000
    )
  }),

  observeEvent(input$ratio_1024, {
    updateNumericInput(
      session = session,
      inputId = "x_axis",
      label = "X-Axis (px):",
      value = 1024,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session = session,
      inputId = "y_axis",
      label = "Y-Axis (px):",
      value = 576,
      min = 1,
      max = 3000
    )
  }),
  
  observeEvent(input$ratio_1152, {
    updateNumericInput(
      session = session,
      inputId = "x_axis",
      label = "X-Axis (px):",
      value = 1152,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session = session,
      inputId = "y_axis",
      label = "Y-Axis (px):",
      value = 648,
      min = 1,
      max = 3000
    )
  }),
  
  observeEvent(input$ratio_1280, {
    updateNumericInput(
      session = session,
      inputId = "x_axis",
      label = "X-Axis (px):",
      value = 1280,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session = session,
      inputId = "y_axis",
      label = "Y-Axis (px):",
      value = 720,
      min = 1,
      max = 3000
    )
  }),
  
  observeEvent(input$ratio_1920, {
    updateNumericInput(
      session = session,
      inputId = "x_axis",
      label = "X-Axis (px):",
      value = 1920,
      min = 1,
      max = 3000
    )
    updateNumericInput(
      session = session,
      inputId = "y_axis",
      label = "Y-Axis (px):",
      value = 1080,
      min = 1,
      max = 3000
    )
  }),

  observeEvent(input$resolution_70, {
    updateNumericInput(
      session = session,
      inputId = "res",
      label = "Resolution of plot, in pixels per inch:",
      value = 70,
      min = 1,
      max = 300
    )
  }),
  
  observeEvent(input$resolution_110, {
    updateNumericInput(
      session = session,
      inputId = "res",
      label = "Resolution of plot, in pixels per inch:",
      value = 110,
      min = 1,
      max = 300
    )
  }),

  observeEvent(input$resolution_144, {
    updateNumericInput(
      session = session,
      inputId = "res",
      label = "Resolution of plot, in pixels per inch:",
      value = 144,
      min = 1,
      max = 300
    )
  })
)
