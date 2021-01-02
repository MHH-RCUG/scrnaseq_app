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

  observeEvent(input$resolution_70, {
    updateNumericInput(
      session = session,
      inputId = "res",
      label = "Resolution: (higher values = bigger text)",
      value = 70,
      min = 1,
      max = 300
    )
  }),

  observeEvent(input$resolution_144, {
    updateNumericInput(
      session = session,
      inputId = "res",
      label = "Resolution: (higher values = bigger text)",
      value = 144,
      min = 1,
      max = 300
    )
  })
)
