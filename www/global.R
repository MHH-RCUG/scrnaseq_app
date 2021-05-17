seurat_functions = function(gene_list, seurat_object, heatmap_slot, heatmap_assay){
    for (i in 1:length(gene_list)){
        stored_FeaturePlots[[i]] <<- # <<- for global assignments
            Seurat::FeaturePlot(
                seurat_object,
                features = unlist(strsplit(gene_list[[i]], "_"))[c(T, F)],
                col = c("lightgrey", param$col),
                label = TRUE
            ) +
            AddStyle(title = gene_list[[i]])
        stored_RidgePlotRaws[[i]] <<-
            Seurat::RidgePlot(
                seurat_object,
                features = unlist(strsplit(gene_list[[i]], "_"))[c(T, F)],
                assay = "RNA",
                slot = "counts"
            ) +
            AddStyle(title = gene_list[[i]],
                     legend_title = "Cluster",
                     fill = param$col_clusters,
                     ylab = "Cluster")
        stored_RidgePlotNorms[[i]] <<-
            Seurat::RidgePlot(
                seurat_object,
                features = unlist(strsplit(gene_list[[i]], "_"))[c(T, F)],
                slot = "data"
            ) +
            AddStyle(title = gene_list[[i]],
                     legend_title = "Cluster",
                     fill = param$col_clusters,
                     ylab = "Cluster")
        stored_ViolinPlotRaws[[i]] <<-
            Seurat::VlnPlot(
                seurat_object,
                features = unlist(strsplit(gene_list[[i]], "_"))[c(T, F)],
                assay = "RNA",
                slot = "counts",
                pt.size = 0.2
            ) +
            AddStyle(title = gene_list[[i]],
                     legend_title = "Cluster",
                     fill = param$col_clusters,
                     xlab = "Cluster")
        stored_ViolinPlotNorms[[i]] <<-
            Seurat::VlnPlot(
                object = seurat_object,
                features = unlist(strsplit(gene_list[[i]], "_"))[c(T, F)],
                pt.size = 0.2
            ) +
            AddStyle(title = gene_list[[i]],
                     legend_title = "Cluster",
                     fill = param$col_clusters,
                     xlab = "Cluster")
    }
    stored_DotPlot <<-
        Seurat::DotPlot(
            seurat_object,
            features = unlist(strsplit(gene_list, "_"))[c(T, F)],
            cols = c("lightgrey", param$col)
        ) +
        AddStyle(title = paste0("Test"), ylab = "Cluster", legend_position = "bottom") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
        guides(size = guide_legend(order = 1))
    tryCatch({
        stored_Heatmap <<-
            Seurat::DoHeatmap(
                seurat_object,
                features = unlist(strsplit(gene_list, "_"))[c(T, F)],
                group.colors = param$col_clusters,
                slot = heatmap_slot,
                assay = heatmap_assay
            ) +
            NoLegend() +
            theme(axis.text.y = element_blank())
    },
    warning = function(w){
        showNotification(
            ui = "One or more of the selected genes could not be found in the
                  default search locations. These gene(s) will not show up in the Heatmap!",
            duration = NULL,
            id = "heatmap_warning",
            type = "warning"
        )
    })
}