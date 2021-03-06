options(shiny.maxRequestSize = 10000*1024^2) #10GB
param = list()
param$col = "palevioletred"
param$col_cluster = "empty"

# Empty object so Input UI gets rendered
features_names_ids = NULL
gene_lists = NULL
heatmap_slots = NULL
heatmap_assays = NULL

# Define global variable to store plots
stored_FeaturePlots = NULL
stored_RidgePlotRaws = NULL
stored_RidgePlotNorms = NULL
stored_ViolinPlotRaws = NULL
stored_ViolinPlotNorms = NULL
stored_DotPlot = NULL
stored_Heatmap = NULL

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  source("server-datainput.R", local = TRUE)
  source("server-settings.R", local = TRUE)
  source("server-plots.R", local = TRUE)
  source("server-download.R", local = TRUE)
  }
