update_x_axis = function(input_x){
  updateNumericInput(
    session,
    "x_axis",
    "X-Axis (px):",
    value = input_x,
    min = 1,
    max = 3000
  )
}

update_y_axis = function(input_y){
  updateNumericInput(
    session,
    "y_axis",
    "Y-Axis (px):",
    value = input_y,
    min = 1,
    max = 3000
  )
}

update_resolution = function(input_res){
  updateNumericInput(
    session = session,
    inputId = "res",
    label = "Resolution of plot, in pixels per inch:",
    value = input_res,
    min = 1,
    max = 300
  )
}

tabItem(
  tabName = "settings",

  #Buttons for aspect ratio
  observeEvent(input$ratio_512, {
    update_x_axis(512)
    update_y_axis(288)
  }),
  
  observeEvent(input$ratio_768, {
    update_x_axis(768)
    update_y_axis(432)
  }),

  observeEvent(input$ratio_1024, {
    update_x_axis(1024)
    update_y_axis(576)
  }),

  observeEvent(input$ratio_1152, {
    update_x_axis(1152)
    update_y_axis(648)
  }),

  observeEvent(input$ratio_1280, {
    update_x_axis(1280)
    update_y_axis(720)
  }),

  observeEvent(input$ratio_1920, {
    update_x_axis(1920)
    update_y_axis(1080)
  }),
  
  
  #Buttons for resolution/ ppi
  observeEvent(input$resolution_70, {
    update_resolution(70)
  }),

  observeEvent(input$resolution_110, {
    update_resolution(110)
  }),

  observeEvent(input$resolution_144, {
    update_resolution(144)
  }),
  
  observeEvent(input$resolution_300, {
    update_resolution(300)
  })
)
