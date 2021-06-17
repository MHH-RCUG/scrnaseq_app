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

  #Buttons for presets
  observeEvent(input$preset_512x288_70, {
    update_x_axis(512)
    update_y_axis(288)
    update_resolution(70)
  }),
  
  observeEvent(input$preset_768x432_110, {
    update_x_axis(768)
    update_y_axis(432)
    update_resolution(110)
  }),

  observeEvent(input$preset_1280x720_144, {
    update_x_axis(1280)
    update_y_axis(720)
    update_resolution(144)
  }),

  observeEvent(input$preset_1920x1080_300, {
    update_x_axis(1920)
    update_y_axis(1080)
    update_resolution(300)
  }),
)
