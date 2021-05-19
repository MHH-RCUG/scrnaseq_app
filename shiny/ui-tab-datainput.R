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
          label = "Upload Seurat File: (.rds)",
          accept = ".rds",
          buttonLabel = "Browse",
          placeholder = "Please upload .rds file!"
        )#fileInput
      ),#box

      conditionalPanel(
        "output.fileUploaded == true",
        box(
          width = 12,
          title = "Select Genes",
          solidHeader = T,
          status = "primary",
          id = "selectgenes",

          selectizeInput(
            inputId = "genes",
            label = "Select Genes Manually:",
            choices = NULL,
            multiple = TRUE
          ),#selectInput
          
          tags$hr(),
          
          selectizeInput(
            inputId = "select_gene_list",
            label = "Select through Gene List:",
            choices = NULL,
            selected = NULL,
            multiple = FALSE
          ),
          actionButton(
            inputId = "button_gene_list",
            label = "Select List"
          ),
          
          tags$hr(),

          fileInput(
            inputId = "xlsx_file",
            label = "Select through Excel File: (.xlsx)",
            accept = ".xlsx",
            buttonLabel = "Browse",
            placeholder = "No file selected"
          ),#fileInput
          
          box(
            width = 12,
            title = "Example Excel File",
            solidHeader = FALSE,
            #status = "info",
            collapsible = TRUE,
            collapsed = TRUE,

            img(
              src = "excel_file_template.png",
              align = "left"
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
          actionButton(inputId = "render_with_defaults",
                       label = "Render with Default Settings",
                       class = "btn btn-success btn-lg btn-block"),
          actionButton(inputId = "goto_settings",
                       label = "Advanced Settings",
                       class = "btn btn-secondary btn-lg btn-block")
        )#box
      )#conditionalPanel
    )#column
    )#fluidRow
  )#tabItem
