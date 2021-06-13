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
      paste(rownames(sc()[["RNA"]][[]]), "_", sc()@assays[["RNA"]]@meta.features[["feature_id"]], sep = "")
    param$col_clusters <<- as.vector(sc()@misc[["colour_lists"]][["seurat_clusters"]])
    gene_lists <<- names(sc()@misc[["gene_lists"]])
    heatmap_assays <<- names(sc()@assays)
    #print(gene_lists)

    suppressWarnings(
      updateSelectizeInput(
        session = session,
        inputId = "genes",
        choices = features_names_ids,
        server = TRUE
      )
    )
    updateSelectizeInput(
      session = session,
      inputId = "colors",
      choices = param$col_clusters,
      selected = param$col_clusters
      )
    updateSelectizeInput(
      session = session,
      inputId = "select_gene_list",
      choices = c("", gene_lists),
      selected = NULL,
      server = TRUE
    )
    updateSelectInput(
      session = session,
      inputId = "heatmap_assay",
      choices = heatmap_assays,
      selected = heatmap_assays[1]
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
      message = toString(e)
      shinyalert(
        title = "Error!",
        text = e,
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
        updateSelectizeInput(
          session = session,
          inputId = "genes",
          choices = features_names_ids,
          selected = x,
          server = TRUE
        )
      )
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

observeEvent(input$button_gene_list,{
  shinyalert(
    title = "Please wait!",
    text = "Upload complete! Please wait while the data is being processed!",
    size = "s",
    type = "info",
    showConfirmButton = FALSE
  )
  
  tmp = NULL
  tryCatch(
    {
      tmp = features_names_ids[unlist(lapply(sc()@misc[["gene_lists"]][[input$select_gene_list]], function(one_gene)
        grep(paste0("^",one_gene,"_"), features_names_ids)))]
    },
    error = function(e){
      message("Error in selecting marker genes")
      message(e)
    },
    warning = function(w){
      message("Warning in selecting marker genes")
      message(w)
    },
    finally = {
      if(length(tmp)==0){
        shinyjs::delay(500, shinyjs::runjs("swal.close();"))
        shinyalert(
          title = "Error!",
          text = "Error: Marker genes could not be selected",
          size = "s",
          type = "error",
          showConfirmButton = TRUE
        )
      } else {
        updateSelectizeInput(
          session = session,
          inputId = "genes",
          choices = features_names_ids,
          selected = tmp,
          server = TRUE
        )
        shinyjs::delay(500, shinyjs::runjs("swal.close();"))
      }
    }
  )
})

observeEvent(input$select_marker_genes, {
  marker_genes_list = NULL
  shinyalert(
    title = "Please wait!",
    text = "Upload complete! Please wait while the data is being processed!",
    size = "s",
    type = "info",
    showConfirmButton = FALSE
  )
  tryCatch(
    {
      marker_genes_list = features_names_ids[unlist(lapply(sc()@misc[["gene_lists"]][["CC_S_phase"]], function(one_gene)
        grep(paste0("^",one_gene,"_"), features_names_ids)))]
    },
    error = function(e){
      message("Error in selecting marker genes")
      message(e)

    },
    warning = function(w){
      message("Warning in selecting marker genes")
      message(w)
    },
    finally = {
      if(length(marker_genes_list)==0){
        shinyjs::delay(500, shinyjs::runjs("swal.close();"))
        shinyalert(
          title = "Error!",
          text = "Error: Marker genes could not be selected",
          size = "s",
          type = "error",
          showConfirmButton = TRUE
        )
      } else {
        updateSelectizeInput(
          session = session,
          inputId = "genes",
          choices = features_names_ids,
          selected = marker_genes_list,
          server = TRUE
        )
        shinyjs::delay(500, shinyjs::runjs("swal.close();"))
      }
    }
  )
})

render_tab_settings = function(){
  output$settings = renderMenu({
    menuItem("Settings", tabName = "settings", icon = icon("cog"), startExpanded = TRUE,
             menuSubItem("General", tabName = "settings_general"),
             menuSubItem("Colors", tabName = "settings_colors"),
             menuSubItem("Heatmap", tabName = "settings_heatmap")
    )
  })
  updateTabItems(
    session = session,
    inputId = "tabs",
    selected = "settings_general"
  )
}

# actionButton renders settings tab and selects it
observeEvent(input$goto_settings,{
  render_tab_settings()
})

observeEvent(input$render_with_defaults,{
  render_tab_settings()
  render_plots()
})
