tabItem(
  tabName = "datainput",
  fluidRow(
    column(
      width = 6,
      box(
        width = 12,
        title = "Upload Data (.rds)",
        solidHeader = T,
        status = "primary",
        id = "uploadbox",

        fileInput(
          inputId = "file_rds",
          label = "Upload Seurat file: (.rds)",
          accept = ".rds",
          buttonLabel = "Browse",
          placeholder = "Please upload .rds file!"
        )#fileInput
      ),#box

      conditionalPanel(
        "output.fileUploaded == true",
        box(
          width = 12,
          title = "Select genes",
          solidHeader = T,
          status = "primary",
          id = "selectgenes",

          selectInput(
            inputId = "genes",
            label = "Select through list:",
            choices = c(1, 2, 3),
            multiple = TRUE
          ),#selectInput

          fileInput(
            inputId = "xlsx_file",
            label = "Select through excel file: (.xlsx)",
            accept = ".xlsx",
            buttonLabel = "Browse",
            placeholder = "No file selected"
          ),#fileInput
          
          box(
            width = 12,
            title = "Example",
            solidHeader = FALSE,
            #status = "info",
            collapsible = TRUE,
            collapsed = TRUE,

            img(
              src = "excel_file_template.png",
              align = "right"
            )
          )#box
        )#box
      ),#conditionalPanel

      conditionalPanel(
        "output.gene_selected == true",
        box(
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          id = "datainput_buttons",
          actionButton(inputId = "goto_settings",
                       label = "Next!",
                       class = "btn btn-success btn-lg btn-block"),
          actionButton(inputId = "render_with_defaults",
                       label = "Render with default settings!",
                       class = "btn btn-success btn-lg btn-block")
        )#box
      )#conditionalPanel
    )#column
    )#fluidRow
  )#tabItem
