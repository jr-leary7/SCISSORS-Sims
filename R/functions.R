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
              Seurat::RunTSNE(reduction = "pca", 
                              dims = 1:30, 
                              dim.embed = 2, 
                              seed.use = 312, 
                              verbose = FALSE) %>% 
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
                              cutoff.vals = NULL, 
                              orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(cutoff.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # evaluate method
  seu.obj <- Seurat::FindNeighbors(seu.obj, 
                                   k.param = sqrt(ncol(seu.obj)), 
                                   nn.method = "annoy", 
                                   annoy.metric = "cosine", 
                                   verbose = FALSE) %>% 
             Seurat::FindClusters(resolution = 0.25, 
                                  algorithm = 1, 
                                  verbose = FALSE, 
                                  random.seed = 312)
  SCISSORS_res <- purrr::map(cutoff.vals, function(x) {
    start_time <- Sys.time()
    reclust_res <- SCISSORS::ReclusterCells(seurat.object = seu.obj, 
                                            auto = TRUE, 
                                            use.parallel = TRUE,
                                            n.cores = 2, 
                                            resolution.vals = seq(0.1, 0.7, by = c(0.1)), 
                                            k.vals = c(10, 25, 40, 60), 
                                            merge.clusters = FALSE, 
                                            redo.embedding = FALSE, 
                                            nn.metric = "cosine", 
                                            use.sct = FALSE, 
                                            n.HVG = 2000, 
                                            n.PC = 30, 
                                            cutoff.score = x, 
                                            random.seed = 312)
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    # integrate with previous clustering
    if (length(reclust_res) == 1) {
      if (length(unique(reclust_res$seurat_clusters)) == 1L) {
        res <- seu.obj 
      } else {
        res <- SCISSORS::IntegrateSubclusters(seu.obj, reclust.results = reclust_res) 
      }
    } else {
      if (all(purrr::map_int(reclust_res, \(x) length(unique(x$seurat_clusters))) == 1L)) {
        res <- seu.obj
      } else {
        res <- SCISSORS::IntegrateSubclusters(seu.obj, reclust.results = reclust_res)
      }
    }
    # metrics
    ARI <- mclust::adjustedRandIndex(orig.clusters, res$seurat_clusters)
    NMI <- aricode::NMI(orig.clusters, 
                        res$seurat_clusters, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = as.integer(res$seurat_clusters))[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = cutoff.vals, 
                               parameter_type = "Silhouette Cutoff", 
                               ari = purrr::map_dbl(SCISSORS_res, \(x) x$ARI),
                               nmi = purrr::map_dbl(SCISSORS_res, \(x) x$NMI), 
                               sil = purrr::map_dbl(SCISSORS_res, \(x) x$sil), 
                               runtime = purrr::map_dbl(SCISSORS_res, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(SCISSORS_res, \(x) x$runtime_units),
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
  reclust_results <- purrr::map(res.vals, function(x) {
    start_time <- Sys.time()
    seu_obj_reclust <- Seurat::FindClusters(seu.obj, 
                                            resolution = x, 
                                            algorithm = 1, 
                                            random.seed = 312, 
                                            verbose = FALSE)
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    ARI <- mclust::adjustedRandIndex(orig.clusters, seu_obj_reclust$seurat_clusters)
    NMI <- aricode::NMI(orig.clusters, 
                        seu_obj_reclust$seurat_clusters, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = as.integer(seu_obj_reclust$seurat_clusters))[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = res.vals, 
                               parameter_type = "Resolution", 
                               ari = purrr::map_dbl(reclust_results, \(x) x$ARI),
                               nmi = purrr::map_dbl(reclust_results, \(x) x$NMI), 
                               sil = purrr::map_dbl(reclust_results, \(x) x$sil), 
                               runtime = purrr::map_dbl(reclust_results, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(reclust_results, \(x) x$runtime_units),
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
  reclust_results <- purrr::map(res.vals, function(x) {
    start_time <- Sys.time()
    seu_obj_reclust <- Seurat::FindClusters(seu.obj, 
                                            resolution = x, 
                                            algorithm = 4, 
                                            random.seed = 312, 
                                            verbose = FALSE)
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    ARI <- mclust::adjustedRandIndex(orig.clusters, seu_obj_reclust$seurat_clusters)
    NMI <- aricode::NMI(orig.clusters, 
                        seu_obj_reclust$seurat_clusters, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = as.integer(seu_obj_reclust$seurat_clusters))[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = res.vals, 
                               parameter_type = "Resolution", 
                               ari = purrr::map_dbl(reclust_results, \(x) x$ARI),
                               nmi = purrr::map_dbl(reclust_results, \(x) x$NMI), 
                               sil = purrr::map_dbl(reclust_results, \(x) x$sil), 
                               runtime = purrr::map_dbl(reclust_results, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(reclust_results, \(x) x$runtime_units),
                               method = "Leiden (Seurat)")
  return(clustering_res)
}

# run hierarchical clustering 
evaluate_hclust <- function(seu.obj = NULL, 
                            k.vals = NULL, 
                            orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(k.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # evaluate method
  hclust_ari <- purrr::map(k.vals, function(x) {
    start_time <- Sys.time()
    hclust_tree <- hclust(SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), method = "ward.D2")
    hclust_res <- cutree(hclust_tree, k = x)
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    ARI <- mclust::adjustedRandIndex(orig.clusters, hclust_res)
    NMI <- aricode::NMI(orig.clusters, 
                        hclust_res, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = hclust_res)[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = k.vals, 
                               parameter_type = "Cut K", 
                               ari = purrr::map_dbl(hclust_ari, \(x) x$ARI),
                               nmi = purrr::map_dbl(hclust_ari, \(x) x$NMI), 
                               sil = purrr::map_dbl(hclust_ari, \(x) x$sil), 
                               runtime = purrr::map_dbl(hclust_ari, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(hclust_ari, \(x) x$runtime_units),
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
  kmeans_ari <- purrr::map(k.vals, function(x) {
    start_time <- Sys.time()
    clust_res <- kmeans(Seurat::Embeddings(seu.obj, "pca")[, 1:30],
                        centers = x,
                        nstart = 5,
                        algorithm = "Hartigan-Wong")
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    ARI <- mclust::adjustedRandIndex(orig.clusters, clust_res$cluster)
    NMI <- aricode::NMI(orig.clusters, 
                        clust_res$cluster, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = clust_res$cluster)[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = k.vals, 
                               parameter_type = "K Centroids", 
                               ari = purrr::map_dbl(kmeans_ari, \(x) x$ARI),
                               nmi = purrr::map_dbl(kmeans_ari, \(x) x$NMI), 
                               sil = purrr::map_dbl(kmeans_ari, \(x) x$sil), 
                               runtime = purrr::map_dbl(kmeans_ari, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(kmeans_ari, \(x) x$runtime_units),
                               method = "K-means (Hartigan-Wong)")
  return(clustering_res)
}

# density-based clustering 
evaluate_dbscan <- function(seu.obj = NULL, 
                            orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  # estimate possible values of epsilon via repeated segmented regression
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
  dbscan_ari <- purrr::map(eps_vals, function(x) {
    start_time <- Sys.time()
    dens_clust <- dbscan::dbscan(scale(seu.obj@reductions$pca@cell.embeddings[, 1:30], scale = FALSE), 
                                 eps = x, 
                                 minPts = 5, 
                                 borderPoints = TRUE)
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    ARI <- mclust::adjustedRandIndex(orig.clusters, dens_clust$cluster)
    NMI <- aricode::NMI(orig.clusters, 
                        dens_clust$cluster, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = dens_clust$cluster)[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = eps_vals, 
                               parameter_type = "Epsilon", 
                               ari = purrr::map_dbl(dbscan_ari, \(x) x$ARI),
                               nmi = purrr::map_dbl(dbscan_ari, \(x) x$NMI), 
                               sil = purrr::map_dbl(dbscan_ari, \(x) x$sil), 
                               runtime = purrr::map_dbl(dbscan_ari, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(dbscan_ari, \(x) x$runtime_units),
                               method = "DBSCAN")
  return(clustering_res)
}

# GiniClust3
evaluate_giniclust <- function(seu.obj = NULL, 
                               res.vals = NULL, 
                               orig.clusters = NULL) {
  # check inputs
  if (is.null(seu.obj) | is.null(res.vals) | is.null(orig.clusters)) { stop("arguments must be non-NULL.") }
  reticulate::use_virtualenv("/nas/longleaf/home/jrleary/Python", required = TRUE)
  reticulate::source_python("./Python/run_GiniClust3.py")
  # evaluate method
  gc3_ari <- purrr::map(res.vals, function(x) {
    start_time <- Sys.time()
    gc3_clust <- run_gc3(count_mat=as.matrix(seu.obj@assays$RNA@counts), res_val=x)
    end_time <- Sys.time()
    time_diff <- end_time - start_time
    runtime <- as.numeric(time_diff[[1]])
    runtime_units <- attributes(time_diff)$units
    ARI <- mclust::adjustedRandIndex(orig.clusters, gc3_clust$finalCluster)
    NMI <- aricode::NMI(orig.clusters, 
                        gc3_clust$finalCluster, 
                        variant = "sqrt")
    sil <- mean(cluster::silhouette(dist = SCISSORS::CosineDist(Seurat::Embeddings(seu.obj, "pca")[, 1:30]), x = as.integer(gc3_clust$finalCluster))[, 3])
    return(list(ARI = ARI, 
                NMI = NMI, 
                sil = sil, 
                runtime = runtime, 
                runtime_units = runtime_units))
  })
  # collate results
  clustering_res <- data.frame(parameter = res.vals, 
                               parameter_type = "Resolution", 
                               ari = purrr::map_dbl(gc3_ari, \(x) x$ARI),
                               nmi = purrr::map_dbl(gc3_ari, \(x) x$NMI), 
                               sil = purrr::map_dbl(gc3_ari, \(x) x$sil), 
                               runtime = purrr::map_dbl(gc3_ari, \(x) x$runtime), 
                               runtime_units = purrr::map_chr(gc3_ari, \(x) x$runtime_units),
                               method = "GiniClust3")
  return(clustering_res)
}

# run all clustering methods for a given dataset 
evaluate_clustering_all <- function(sim.data = NULL) {
  # check inputs
  if (is.null(sim.data)) { stop("sim.data must be non-NULL.") }
  # evaluate all methods 
  SCISSORS_res <- evaluate_SCISSORS(seu.obj = sim.data, 
                                    cutoff.vals = c(0.15, 0.25, 0.35), 
                                    orig.clusters = sim.data$cellPopulation)
  Louvain_res <- evaluate_Seurat_Louvain(seu.obj = sim.data, 
                                         res.vals = seq(0.1, 1.3, by = 0.1), 
                                         orig.clusters = sim.data$cellPopulation)
  Leiden_res <- evaluate_Seurat_Leiden(seu.obj = sim.data, 
                                       res.vals = seq(0.1, 1.3, by = 0.1), 
                                       orig.clusters = sim.data$cellPopulation)
  kmeans_res <- evaluate_kmeans(seu.obj = sim.data, 
                                k.vals = c(3:10), 
                                orig.clusters = sim.data$cellPopulation)
  hclust_res <- evaluate_hclust(seu.obj = sim.data, 
                                k.vals = c(3:10), 
                                orig.clusters = sim.data$cellPopulation)
  dbscan_res <- evaluate_dbscan(seu.obj = sim.data, 
                                orig.clusters = sim.data$cellPopulation)
  giniclust3_res <- evaluate_giniclust(seu.obj = sim.data,
                                       res.vals = seq(0.1, 1.3, by = 0.1), 
                                       orig.clusters = sim.data$cellPopulation)
  # collate results
  clustering_res_all <- purrr::reduce(list(SCISSORS_res, Louvain_res, Leiden_res, kmeans_res, hclust_res, dbscan_res, giniclust3_res), 
                                      rbind)
  return(clustering_res_all)
}
