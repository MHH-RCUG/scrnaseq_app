# rds file upload
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

# rds file processing
observeEvent(sc(), {
  if (!is.null(sc())) {
    features_names_ids <<-
      paste(rownames(sc()[["RNA"]][[]]), "_", sc()[["RNA"]][[]][, 1], sep = "")
    param$col <<- as.vector(sc()@misc[["colors"]][["seurat_clusters"]])
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

# Excel file processing
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

# Returns true if .rds file has been processed and thus making gene slection available
output$fileUploaded = reactive({
  return(!is.null(sc()))
})#reactive
outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)

output$gene_selected = reactive({
  return(!is.null(input$genes))
})#reactive
outputOptions(output, 'gene_selected', suspendWhenHidden = FALSE)

# actionButton renders settings tab and selects it
observeEvent(input$goto_settings,{
  output$settings = renderMenu({
    menuItem("Settings", tabName = "settings", icon = icon("cog"))
  })
  updateTabItems(session = session,
                 inputId = "tabs",
                 selected = "settings")
})
