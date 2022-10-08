# simulate scRNA-seq data 
simulate_clusters <- function(ref.data = NULL, 
                              clust.n = NULL, 
                              genes.p = NULL, 
                              genes.mean.FC = NULL, 
                              genes.sd.FC = NULL) {
  # simulate dataset 
  scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref_dataset,
                                                          sceUMI = TRUE,
                                                          useUMI = TRUE,
                                                          protocol = "droplet",
                                                          numCells = clust.n,
                                                          popHet = c(1, 1), 
                                                          usePops = list(propGenes = genes.p,
                                                                         fc_mean = genes.mean.FC,
                                                                         fc_sd = genes.sd.FC))
  sim_data <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.data)
  # prepare Seurat object
  sim_data <- Seurat::CreateSeuratObject(counts = counts(sim_data), 
                                         meta.data = (as.data.frame(SingleCellExperiment::colData(sim_data))), 
                                         min.cells = 3, 
                                         min.features = 200) %>% 
              Seurat::NormalizeData(verbose = FALSE) %>%
              Seurat::FindVariableFeatures(verbose = FALSE) %>% 
              Seurat::ScaleData(verbose = FALSE) %>% 
              Seurat::RunPCA(features = Seurat::VariableFeatures(.), verbose = FALSE) %>% 
              Seurat::RunUMAP(reduction = "pca", dims = 1:30, dim.embed = 2, verbose = FALSE)
  return(sim_data)
}
