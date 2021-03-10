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
    param$col_clusters <<- as.vector(sc()@misc[["colors"]][["seurat_clusters"]])
    #param$col_cluster <<- sc@misc$colors$seurat_clusters
    
    suppressWarnings(
      updateSelectInput(
        session = session, 
        inputId = "genes",
        label = "Select Genes:",
        choices = features_names_ids
      ) 
    )
    
    updateSelectInput(
      session = session,
      inputId = "colors",
      label = "Select colors for plots:",
      choices = param$col_clusters,
      selected = param$col_clusters
      )
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
    {
      excel_list = unlist(excel_genes()[, 1])
      x = features_names_ids[unlist(lapply(excel_list, function(one_gene)
        grep(one_gene, features_names_ids)))]
    },
    error = function(e){
      message("Excel file upload: There was an error")
      message(e)
      shinyalert(
        title = "Error!",
        text = "The Excel file was empty!",
        size = "s",
        type = "error",
        showConfirmButton = TRUE
      )
    },
    warning = function(w){
      message("Excel file upload: There was a warning message.")
      message(w)
    },
    finally = {
      suppressWarnings(
        updateSelectInput(
          session = session,
          inputId = "genes",
          label = "Select Genes:",
          choices = features_names_ids,
          selected = x
        )
      )
      
      message("Excel file upload: tryCatch is finished.")
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
