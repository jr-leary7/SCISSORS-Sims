# fetch lung reference dataset
get_lung_ref <- function() {
  ref <- scRNAseq::ZilionisLungData(which = "human", filter = TRUE)
  ref <- ref[rowSums(BiocGenerics::counts(ref) > 0) >= 3, ]  # genes found in >= 3 cells
  return(ref)
}

# fetch pancreas reference dataset 
get_panc_ref <- function() {
  ref <- scRNAseq::BaronPancreasData(which = "human")
  ref <- ref[rowSums(BiocGenerics::counts(ref) > 0) >= 3, ]  # genes found in >= 3 cells
  return(ref)
}

# simulate scRNA-seq data 
simulate_clusters <- function(ref.data = NULL, 
                              clust.n = NULL, 
                              genes.n = NULL, 
                              genes.p = NULL, 
                              genes.mean.FC = NULL, 
                              genes.sd.FC = NULL) {
  # simulate dataset 
  scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.data,
                                                          sceUMI = TRUE,
                                                          useUMI = TRUE,
                                                          protocol = "droplet",
                                                          numCells = clust.n, 
                                                          numGenes = genes.n, 
                                                          popHet = c(1, 1), 
                                                          usePops = list(propGenes = genes.p,
                                                                         fc_mean = genes.mean.FC,
                                                                         fc_sd = genes.sd.FC))
  sim_data <- scaffold::simulateScaffold(scaffoldParams = scaffold_params, originalSCE = ref.data)
  # prepare Seurat object
  sim_data <- Seurat::CreateSeuratObject(counts = BiocGenerics::counts(sim_data), 
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
