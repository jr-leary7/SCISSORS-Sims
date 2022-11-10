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
                              genes.p = NULL, 
                              genes.mean.FC = NULL, 
                              genes.sd.FC = NULL) {
  # simulate dataset 
  scaffold_params <- scaffold::estimateScaffoldParameters(sce = ref.data,
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
  sim_data <- Seurat::CreateSeuratObject(counts = BiocGenerics::counts(sim_data), 
                                         meta.data = (as.data.frame(SingleCellExperiment::colData(sim_data))), 
                                         min.cells = 3, 
                                         min.features = 200) %>% 
              Seurat::NormalizeData(verbose = FALSE) %>%
              Seurat::FindVariableFeatures(verbose = FALSE) %>% 
              Seurat::ScaleData(verbose = FALSE) %>% 
              Seurat::RunPCA(features = Seurat::VariableFeatures(.), verbose = FALSE) %>% 
              Seurat::RunUMAP(reduction = "pca", 
                              dims = 1:30,
                              n.components = 2, 
                              metric = "cosine", 
                              verbose = FALSE) %>% 
              Seurat::FindNeighbors(nn.method = "annoy", 
                                    annoy.metric = "cosine", 
                                    verbose = FALSE) %>% 
              Seurat::FindClusters(resolution = 0.3, 
                                   algorithm = 1, 
                                   verbose = FALSE, 
                                   random.seed = 312)
  return(sim_data)
}

# run SCISSORS clustering
evaluate_SCISSORS <- function(seu.obj = NULL, 
                              orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # evaluate method
  res <- SCISSORS::ReclusterCells(seurat.object = seu.obj, 
                                  which.clust = unique(seu.obj$seurat_clusters), 
                                  resolution.vals = seq(0.1, 0.7, by = c(0.1)), 
                                  k.vals = c(10, 25, 50, 75), 
                                  merge.clusters = TRUE, 
                                  use.sct = FALSE, 
                                  n.HVG = 2000, 
                                  n.PC = 30, 
                                  cutoff.score = 0.25, 
                                  random.seed = 312)
  # collate results
  clustering_res <- data.frame(parameter = NA_real_, 
                               parameter_type = NA_character_, 
                               ari = mclust::adjustedRandIndex(orig.clusters, res$seurat_clusters), 
                               method = "SCISSORS")
  return(clustering_res)
}

# run Seurat Louvain clustering 
evaluate_Seurat_Louvain <- function(seu.obj = NULL, 
                                    res.vals = NULL, 
                                    orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(res.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # evaluate method
  reclust_results <- purrr::map_dbl(res.vals, function(x) {
    seu_obj_reclust <- Seurat::FindClusters(seu.obj, 
                                            resolution = x, 
                                            algorithm = 1, 
                                            random.seed = 312, 
                                            verbose = FALSE)
    ARI <- mclust::adjustedRandIndex(orig.clusters, seu_obj_reclust$seurat_clusters)
    return(ARI)
  })
  # collate results
  clustering_res <- data.frame(parameter = res.vals, 
                               parameter_type = "Resolution", 
                               ari = reclust_results, 
                               method = "Louvain (Seurat)")
  return(clustering_res)
}

# run Seurat Leiden clustering 
evaluate_Seurat_Leiden <- function(seu.obj = NULL, 
                                   res.vals = NULL, 
                                   orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(res.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # set Python virtual environment
  reticulate::use_virtualenv("/nas/longleaf/home/jrleary/Python", required = TRUE)
  # evaluate method
  reclust_results <- purrr::map_dbl(res.vals, function(x) {
    seu_obj_reclust <- Seurat::FindClusters(seu.obj, 
                                            resolution = x, 
                                            algorithm = 4, 
                                            random.seed = 312, 
                                            verbose = FALSE)
    ARI <- mclust::adjustedRandIndex(orig.clusters, seu_obj_reclust$seurat_clusters)
    return(ARI)
  })
  # collate results
  clustering_res <- data.frame(parameter = res.vals, 
                               parameter_type = "Resolution", 
                               ari = reclust_results, 
                               method = "Louvain (Seurat)")
  return(clustering_res)
}

# run hierarchical clustering 
evaluate_hclust <- function(seu.obj = NULL, 
                            k.vals = NULL, 
                            orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(k.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # evaluate method
  hclust_tree <- hclust(dist(seu.obj@reductions$pca@cell.embeddings[, 1:30]), method = "ward.D2")
  hclust_ari <- purrr::map_dbl(k.vals, function(x) {
    hclust_res <- cutree(hclust_tree, k = x)
    ARI <- mclust::adjustedRandIndex(orig.clusters, hclust_res)
    return(ARI)
  })
  # collate results
  clustering_res <- data.frame(parameter = k.vals, 
                               parameter_type = "Cut K", 
                               ari = hclust_ari, 
                               method = "Hierarchical (Ward)")
  return(clustering_res)
}

# run k-means clustering 
evaluate_kmeans <- function(seu.obj = NULL, 
                            k.vals = NULL, 
                            orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(k.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # evaluate method
  kmeans_ari <- purrr::map_dbl(k.vals, function(x) {
    clust_res <- kmeans(seu.obj@reductions$pca@cell.embeddings[, 1:30],
                        centers = x,
                        nstart = 5,
                        algorithm = "Hartigan-Wong")
    ARI <- mclust::adjustedRandIndex(orig.clusters, clust_res$cluster)
    return(ARI)
  })
  # collate results
  clustering_res <- data.frame(parameter = k.vals, 
                               parameter_type = "K Centroids", 
                               ari = kmeans_ari, 
                               method = "K-means (Hartigan-Wong)")
  return(clustering_res)
}

evaluate_dbscan <- function(seu.obj = NULL, 
                            orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # choose possible values of epsilon via segmented regression
  Y <- sort(dbscan::kNNdist(seu.obj@reductions$pca@cell.embeddings[, 1:30], k = 10))
  X <- 1:ncol(seu.obj)
  base_lm <- lm(Y ~ X)
  breakpoint_vals <- purrr::map(c(1:4), function(x) {
    seg_reg <- segmented::segmented.lm(base_lm, npsi = x)
    breakpoints <- as.numeric(seg_reg$psi[, 2])
    return(breakpoints)
  })
  breakpoint_vals <- unique(unlist(breakpoint_vals))
  eps_vals <- purrr::map_dbl(breakpoint_vals, function(x) {
    eps <- Y[which.min(abs(X - x))]
    return(eps)
  })
  eps_vals <- unique(round(eps_vals))
  # evaluate method 
  dbscan_ari <- purrr::map_dbl(eps_vals, function(x) {
    dens_clust <- dbscan(scale(seu.obj@reductions$pca@cell.embeddings[, 1:30], scale = FALSE), 
                         eps = x, 
                         minPts = 5, 
                         borderPoints = TRUE)
    ARI <- mclust::adjustedRandIndex(orig.clusters, dens_clust$cluster)
    return(ARI)
  })
  # collate results
  clustering_res <- data.frame(parameter = eps_vals, 
                               parameter_type = "Epsilon", 
                               ari = dbscan_ari, 
                               method = "DBSCAN")
  return(clustering_res)
}
# run all clustering methods for a given dataset 
evaluate_clustering_all <- function(sim.data = NULL) {
  # check inputs
  if (is.null(sim.data)) { stop("sim.data must be non-NULL.") }
  # evaluate all methods 
  SCISSORS_res <- evaluate_SCISSORS(seu.obj = sim.data, orig.clusters = sim.data$cellPopulation)
  Louvain_res <- evaluate_Seurat_Louvain(seu.obj = sim.data, 
                                         res.vals = seq(0.1, 1.5, by = 0.1), 
                                         orig.clusters = sim.data$cellPopulation)
  Leiden_res <- evaluate_Seurat_Leiden(seu.obj = sim.data, 
                                       res.vals = seq(0.1, 1.5, by = 0.1), 
                                       orig.clusters = sim.data$cellPopulation)
  kmeans_res <- evaluate_kmeans(seu.obj = sim.data, 
                                k.vals = c(2:10), 
                                orig.clusters = sim.data$cellPopulation)
  hclust_res <- evaluate_hclust(seu.obj = sim.data, 
                                k.vals = c(2:10), 
                                orig.clusters = sim.data$cellPopulation)
  dbscan_res <- evaluate_dbscan(seu.obj = sim.data, 
                                orig.clusters = sim.data$cellPopulation)
  # collate results
  clustering_res_all <- purrr::reduce(list(SCISSORS_res, Louvain_res, Leiden_res, kmeans_res, hclust_res, dbscan_res), rbind)
  return(clustering_res_all)
}
