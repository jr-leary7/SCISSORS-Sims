library(future)
library(targets)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)

tar_option_set(packages = c("scaffold", "SingleCellExperiment", "Seurat", "BiocGenerics", "magrittr"),
               format = "rds", 
               garbage_collection = TRUE, 
               error = "continue")

# options(clustermq.scheduler = "slurm")
# options(clustermq.template = "clustermq.tmpl")

tar_source("R/functions.R")

# Replace the target list below with your own:
list(
  ##### reference datasets #####
  tar_target(panc_ref, get_panc_ref()),
  tar_target(lung_ref, get_lung_ref()), 
  
  ##### pancreas sim datasets #####
  
  ### 1,000 cells
  tar_target(sim_panc_ncell1000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(500, 400, 100), 
                                                 genes.p = c(0, 0.25, 0.1), 
                                                 genes.mean.FC = c(0, 0.7, 1.4), 
                                                 genes.sd.FC = c(0, 0.4, 0.15))), 
  tar_target(sim_panc_ncell1000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(200, 200, 200, 200, 200), 
                                                 genes.p = c(0, 0.1, 0.3, 0.4, 0.05), 
                                                 genes.mean.FC = c(0, 1.2, 0.9, 2.5, 1.4), 
                                                 genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8))), 
  tar_target(sim_panc_ncell1000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(250, 250, 150, 150, 100, 50, 50), 
                                                 genes.p = c(0, 0.1, 0.3, 0.4, 0.05, 0.2, 0.3), 
                                                 genes.mean.FC = c(0, 1.2, 0.8, 2.3, 1.4, 1.5, 0.6), 
                                                 genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8, 0.4, 0.4))), 
  ### 3,000 cells
  tar_target(sim_panc_ncell3000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(1500, 1200, 300), 
                                                 genes.p = c(0, 0.25, 0.1), 
                                                 genes.mean.FC = c(0, 0.7, 1.4), 
                                                 genes.sd.FC = c(0, 0.4, 0.15))), 
  tar_target(sim_panc_ncell3000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(600, 950, 250, 800, 400), 
                                                 genes.p = c(0, 0.15, 0.8, 0.4, 0.3), 
                                                 genes.mean.FC = c(0, 0.75, 3.1, 2.5, 1.4), 
                                                 genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8))), 
  tar_target(sim_panc_ncell3000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(600, 700, 50, 150, 500, 450, 550), 
                                                 genes.p = c(0, 0.1, 0.3, 0.4, 0.05, 0.2, 0.15), 
                                                 genes.mean.FC = c(0, 1.2, 3.2, 2.5, 0.4, 0.4, 1.5), 
                                                 genes.sd.FC = c(0, 0.2, 0.5, 0.2, 0.75, 0.4, 0.2))), 
  ### 5,000 cells
  tar_target(sim_panc_ncell5000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(2000, 2800, 200), 
                                                 genes.p = c(0, 0.25, 0.4), 
                                                 genes.mean.FC = c(0, 1.1, 0.2), 
                                                 genes.sd.FC = c(0, 0.4, 0.1))), 
  tar_target(sim_panc_ncell5000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(3000, 800, 200, 700, 300), 
                                                 genes.p = c(0, 0.1, 0.5, 0.2, 0.07), 
                                                 genes.mean.FC = c(0, 0.2, 1.1, 4.1, 2.4), 
                                                 genes.sd.FC = c(0, 0.25, 0.25, 0.4, 0.1))), 
  tar_target(sim_panc_ncell5000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                 clust.n = c(100, 900, 50, 950, 1000, 600, 400), 
                                                 genes.p = c(0, 0.25, 0.3, 0.4, 0.03, 0.2, 0.15), 
                                                 genes.mean.FC = c(0, 1.2, 3.2, 2.5, 3.4, 0.4, 1.3), 
                                                 genes.sd.FC = c(0, 0.25, 0.5, 0.2, 0.5, 0.4, 0.1))), 
  
  ##### lung sim datasets #####
  
  ### 1,000 cells
  tar_target(sim_lung_ncell1000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(500, 400, 100), 
                                                 genes.p = c(0, 0.25, 0.1), 
                                                 genes.mean.FC = c(0, 0.7, 1.4), 
                                                 genes.sd.FC = c(0, 0.4, 0.15))), 
  tar_target(sim_lung_ncell1000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(200, 200, 200, 200, 200), 
                                                 genes.p = c(0, 0.1, 0.3, 0.4, 0.05), 
                                                 genes.mean.FC = c(0, 1.2, 0.9, 2.5, 1.4), 
                                                 genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8))), 
  tar_target(sim_lung_ncell1000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(250, 250, 150, 150, 100, 50, 50), 
                                                 genes.p = c(0, 0.1, 0.3, 0.4, 0.05, 0.2, 0.3), 
                                                 genes.mean.FC = c(0, 1.2, 0.8, 2.3, 1.4, 1.5, 0.6), 
                                                 genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8, 0.4, 0.4))), 
  ### 3,000 cells
  tar_target(sim_lung_ncell3000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(1500, 1200, 300), 
                                                 genes.p = c(0, 0.25, 0.1), 
                                                 genes.mean.FC = c(0, 0.7, 1.4), 
                                                 genes.sd.FC = c(0, 0.4, 0.15))), 
  tar_target(sim_lung_ncell3000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(600, 950, 250, 800, 400), 
                                                 genes.p = c(0, 0.15, 0.8, 0.4, 0.3), 
                                                 genes.mean.FC = c(0, 0.75, 3.1, 2.5, 1.4), 
                                                 genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8))), 
  tar_target(sim_lung_ncell3000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(600, 700, 50, 150, 500, 450, 550), 
                                                 genes.p = c(0, 0.1, 0.3, 0.4, 0.05, 0.2, 0.15), 
                                                 genes.mean.FC = c(0, 1.2, 3.2, 2.5, 0.4, 0.4, 1.5), 
                                                 genes.sd.FC = c(0, 0.2, 0.5, 0.2, 0.75, 0.4, 0.2))), 
  ### 5,000 cells
  tar_target(sim_lung_ncell5000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(2000, 2800, 200), 
                                                 genes.p = c(0, 0.25, 0.4), 
                                                 genes.mean.FC = c(0, 1.1, 0.2), 
                                                 genes.sd.FC = c(0, 0.4, 0.1))), 
  tar_target(sim_lung_ncell5000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(3000, 800, 200, 700, 300), 
                                                 genes.p = c(0, 0.1, 0.5, 0.2, 0.07), 
                                                 genes.mean.FC = c(0, 0.2, 1.1, 4.1, 2.4), 
                                                 genes.sd.FC = c(0, 0.25, 0.25, 0.4, 0.1))), 
  tar_target(sim_lung_ncell5000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                 clust.n = c(100, 900, 50, 950, 1000, 600, 400), 
                                                 genes.p = c(0, 0.25, 0.3, 0.4, 0.03, 0.2, 0.15), 
                                                 genes.mean.FC = c(0, 1.2, 3.2, 2.5, 3.4, 0.4, 1.3), 
                                                 genes.sd.FC = c(0, 0.25, 0.5, 0.2, 0.5, 0.4, 0.2))), 
  
  ##### clustering evaluation #####
  
  ### pancreas sim datasets 
  # 1,000 cells
  tar_target(clustres_panc_ncell1000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell1000_nclust3)), 
  tar_target(clustres_panc_ncell1000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell1000_nclust5)), 
  tar_target(clustres_panc_ncell1000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell1000_nclust7)), 
  # 3,000 cells
  tar_target(clustres_panc_ncell3000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell3000_nclust3)), 
  tar_target(clustres_panc_ncell3000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell3000_nclust5)), 
  tar_target(clustres_panc_ncell3000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell3000_nclust7)), 
  # 5,000 cells
  tar_target(clustres_panc_ncell5000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell5000_nclust3)), 
  tar_target(clustres_panc_ncell5000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell5000_nclust5)), 
  tar_target(clustres_panc_ncell5000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell5000_nclust7)), 
  ### lung datasets 
  # 1,000 cells
  tar_target(clustres_lung_ncell1000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell1000_nclust3)), 
  tar_target(clustres_lung_ncell1000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell1000_nclust5)), 
  tar_target(clustres_lung_ncell1000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell1000_nclust7)), 
  # 3,000 cells
  tar_target(clustres_lung_ncell3000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell3000_nclust3)), 
  tar_target(clustres_lung_ncell3000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell3000_nclust5)), 
  tar_target(clustres_lung_ncell3000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell3000_nclust7)), 
  # 5,000 cells
  tar_target(clustres_lung_ncell5000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell5000_nclust3)), 
  tar_target(clustres_lung_ncell5000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell5000_nclust5)), 
  tar_target(clustres_lung_ncell5000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell5000_nclust7)), 
  
  ##### simulation analysis 
  tar_render(simulation_summary, "SCISSORS_Simulation_Summary.Rmd")
)
