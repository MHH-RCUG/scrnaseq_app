sc = reactive({
  inFile = input$file_rds
  if (is.null(inFile)) {
    tmp = NULL
  } else {
    shinyalert(
      title = "Please wait!",
      text = "Upload complete! Please wait while the data is being processed!",
      size = "s",
      type = "info",
      showConfirmButton = FALSE
    )
    tmp = readRDS(inFile$datapath)
  }
  return(tmp)
})

observeEvent(sc(), {
  if (!is.null(sc())) {
    features_names_ids <<-
      paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][, 1], sep = "")
    updateSelectInput(session, "genes", "Select Genes:", features_names_ids)
  }
  shinyjs::delay(500, shinyjs::runjs("swal.close();"))
})#observeEvent

# Excel file upload
excel_genes = reactive({
  inFile = input$xlsx_file
  if (is.null(inFile)) {
    tmp = NULL
  } else {
    tmp = readxl::read_excel(
      path = inFile$datapath,
      sheet = 1,
      range = cell_cols("A:B"),
      col_names = FALSE
    )
    Sys.sleep(1)
    shinyalert(
      title = "Please wait!",
      text = "Upload complete! Please wait while the data is being processed!",
      size = "s",
      type = "info",
      showConfirmButton = FALSE
    )
  }
  return(tmp)
})

observeEvent(excel_genes(), {
  tryCatch(
    expr = {
      excel_list = unlist(excel_genes()[, 1])
      x = features_names_ids[unlist(lapply(excel_list, function(one_gene)
        grep(one_gene, features_names_ids)))]
      updateSelectInput(
        session = session,
        inputId = "genes",
        label = "Select Genes:",
        choices = features_names_ids,
        selected = x
      )
    },
    error = function(e){
      shinyalert(
        title = "Error!",
        text = "The Excel file was empty!",
        size = "s",
        type = "error",
        showConfirmButton = TRUE
      )
    },
    warning = function(w){
      message("Excel file: There was a warning message.")
    },
    finally = {
      message("Excel file: tryCatch is finished.")
      shinyjs::delay(500, shinyjs::runjs("swal.close();"))
    }
  )
})#observeEvent

output$fileUploaded = reactive({
  return(!is.null(sc()))
})#reactive
outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)

observeEvent(input$goto_settings,{
  output$settings = renderMenu({
    menuItem("Settings", tabName = "settings", icon = icon("cog"), selected = TRUE)
  })
  updateTabItems(session = session,
                 inputId = "tabs",
                 selected = "settings")
})