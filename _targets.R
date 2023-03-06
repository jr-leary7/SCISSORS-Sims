library(future)
library(targets)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 100000 * 1024^2)  # 100GB max size 

tar_option_set(packages = c("stats", 
                            "mclust", 
                            "Seurat", 
                            "dbscan", 
                            "cluster", 
                            "aricode", 
                            "scaffold", 
                            "magrittr", 
                            "segmented",
                            "tidyverse",
                            "clusterSim", 
                            "reticulate", 
                            "BiocGenerics", 
                            "SingleCellExperiment"),
               format = "rds", 
               garbage_collection = TRUE, 
               error = "continue")

tar_source("R/functions.R")

list(
  ##### reference datasets #####
  tar_target(panc_ref, get_panc_ref()),
  tar_target(lung_ref, get_lung_ref()), 
  
  ##### pancreas simulations #####
  
  ### 1,000 cells
  tar_target(sim_panc_ncell1000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(500, 400, 100), 
                                                           genes.p = c(0.05, 0.25, 0.01), 
                                                           genes.mean.FC = c(0.5, 0.7, 0.4), 
                                                           genes.sd.FC = c(0.2, 0.4, 0.15))), 
  tar_target(sim_panc_ncell1000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(200, 200, 200, 200, 200), 
                                                           genes.p = c(0.08, 0.1, 0.15, 0.1, 0.2), 
                                                           genes.mean.FC = c(0.4, 0.9, 1.5, 1.1, 0.5), 
                                                           genes.sd.FC = c(0.5, 0.3, 0.2, 0.6, 0.4))), 
  tar_target(sim_panc_ncell1000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(250, 250, 150, 150, 100, 50, 50), 
                                                           genes.p = c(0.05, 0.1, 0.3, 0.4, 0.15, 0.05, 0.3), 
                                                           genes.mean.FC = c(0.5, 0.8, 0.8, 2.3, 1.4, 1.5, 1.9), 
                                                           genes.sd.FC = c(0.25, 0.5, 0.5, 0.2, 0.8, 0.25, 0.5))), 
  ### 3,000 cells
  tar_target(sim_panc_ncell3000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(1500, 1200, 300), 
                                                           genes.p = c(0.1, 0.05, 0.05), 
                                                           genes.mean.FC = c(0.2, 0.1, 0.5), 
                                                           genes.sd.FC = c(0.7, 0.4, 0.42))), 
  tar_target(sim_panc_ncell3000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(600, 950, 250, 800, 400), 
                                                           genes.p = c(0, 0.15, 0.8, 0.4, 0.3), 
                                                           genes.mean.FC = c(0, 0.75, 3.1, 2.5, 1.4), 
                                                           genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8))), 
  tar_target(sim_panc_ncell3000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(600, 700, 50, 150, 500, 450, 550), 
                                                           genes.p = c(0.05, 0.1, 0.05, 0.4, 0.05, 0.2, 0.9), 
                                                           genes.mean.FC = c(0.75, 0.4, 3.2, 2.5, 0.2, 0.4, 1.5), 
                                                           genes.sd.FC = c(0.2, 0.05, 0.2, 0.2, 0.75, 0.4, 0.55))), 
  ### 5,000 cells
  tar_target(sim_panc_ncell5000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(2000, 2800, 200), 
                                                           genes.p = c(0, 0.25, 0.4), 
                                                           genes.mean.FC = c(0, 1.1, 0.2), 
                                                           genes.sd.FC = c(0, 0.4, 0.1))), 
  tar_target(sim_panc_ncell5000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(3000, 800, 200, 700, 300), 
                                                           genes.p = c(0.03, 0.1, 0.15, 0.2, 0.03), 
                                                           genes.mean.FC = c(0.15, 0.2, 3.1, 4.1, 0.18), 
                                                           genes.sd.FC = c(0.08, 0.25, 1.1, 0.4, 0.1))), 
  tar_target(sim_panc_ncell5000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                           clust.n = c(100, 900, 50, 950, 1000, 600, 400), 
                                                           genes.p = c(0.08, 0.25, 0.4, 0.4, 0.03, 0.2, 0.15), 
                                                           genes.mean.FC = c(1.3, 1.2, 2.4, 0.3, 1.5, 0.4, 1.3), 
                                                           genes.sd.FC = c(0.07, 0.5, 0.2, 0.2, 0.5, 0.4, 0.6))), 
  ### 10,000 cells
  tar_target(sim_panc_ncell10000_nclust3, simulate_clusters(ref.data = panc_ref, 
                                                            clust.n = c(5000, 3500, 1500), 
                                                            genes.p = c(0.05, 0.1, 0.25), 
                                                            genes.mean.FC = c(0.7, 0.3, 0.1), 
                                                            genes.sd.FC = c(0.7, 0.3, 0.4))), 
  tar_target(sim_panc_ncell10000_nclust5, simulate_clusters(ref.data = panc_ref, 
                                                            clust.n = c(200, 800, 2000, 6000, 1000), 
                                                            genes.p = c(0.04, 0.15, 0.1, 0.04, 0.07), 
                                                            genes.mean.FC = c(2.3, 0.4, 0.4, 0.08, 1.3), 
                                                            genes.sd.FC = c(0.08, 0.5, 0.1, 0.3, 0.13))), 
  tar_target(sim_panc_ncell10000_nclust7, simulate_clusters(ref.data = panc_ref, 
                                                            clust.n = c(300, 700, 50, 2950, 4500, 500, 1000), 
                                                            genes.p = c(0, 0.5, 0.12, 0.12, 0.03, 0.14, 0.19), 
                                                            genes.mean.FC = c(0, 2.1, 1.8, 0.07, 3.8, 0.8, 1.7), 
                                                            genes.sd.FC = c(0, 0.5, 1.2, 1.2, 0.4, 0.9, 1.4))), 
  
  ##### lung simulations #####
  
  ### 1,000 cells
  tar_target(sim_lung_ncell1000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(500, 400, 100), 
                                                           genes.p = c(0.05, 0.25, 0.05), 
                                                           genes.mean.FC = c(0.2, 0.7, 0.2), 
                                                           genes.sd.FC = c(0.03, 0.4, 0.03))), 
  tar_target(sim_lung_ncell1000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(200, 200, 200, 200, 200), 
                                                           genes.p = c(0.02, 0.02, 0.3, 0.4, 0.05), 
                                                           genes.mean.FC = c(0.3, 0.3, 0.9, 2.5, 1.4), 
                                                           genes.sd.FC = c(0.05, 0.05, 0.5, 0.2, 0.8))), 
  tar_target(sim_lung_ncell1000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(250, 250, 150, 150, 100, 50, 50), 
                                                           genes.p = c(0.03, 0.1, 0.3, 0.4, 0.05, 0.03, 0.3), 
                                                           genes.mean.FC = c(0.25, 1.2, 0.8, 2.3, 1.4, 0.25, 0.4),
                                                           genes.sd.FC = c(0.05, 0.3, 0.5, 0.2, 0.8, 0.01, 0.4))), 
  ### 3,000 cells
  tar_target(sim_lung_ncell3000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(1500, 1200, 300), 
                                                           genes.p = c(0.02, 0.25, 0.02), 
                                                           genes.mean.FC = c(0.75, 0.7, 1.8), 
                                                           genes.sd.FC = c(0.35, 0.4, 0.35))), 
  tar_target(sim_lung_ncell3000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(600, 950, 250, 800, 400), 
                                                           genes.p = c(0, 0.15, 0.8, 0.4, 0.3), 
                                                           genes.mean.FC = c(0, 0.75, 1.1, 0.5, 1.4), 
                                                           genes.sd.FC = c(0, 0.3, 0.5, 0.2, 0.8))), 
  tar_target(sim_lung_ncell3000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(600, 700, 50, 150, 500, 450, 550), 
                                                           genes.p = c(0.04, 0.05, 0.07, 0.4, 0.05, 0.2, 0.15), 
                                                           genes.mean.FC = c(0.25, 0.3, 3.4, 2.5, 0.4, 0.4, 1.5), 
                                                           genes.sd.FC = c(0.15, 0.17, 0.05, 0.2, 0.75, 0.4, 0.2))), 
  ### 5,000 cells
  tar_target(sim_lung_ncell5000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(2000, 2800, 200), 
                                                           genes.p = c(0, 0.25, 0.4), 
                                                           genes.mean.FC = c(0, 1.1, 0.2), 
                                                           genes.sd.FC = c(0, 0.4, 0.1))), 
  tar_target(sim_lung_ncell5000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(3000, 800, 200, 700, 300), 
                                                           genes.p = c(0.05, 0.1, 0.1, 0.2, 0.05), 
                                                           genes.mean.FC = c(0.3, 0.2, 1.8, 4.1, 1.8), 
                                                           genes.sd.FC = c(0.05, 0.25, 0.15, 0.4, 0.1))), 
  tar_target(sim_lung_ncell5000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                           clust.n = c(100, 1350, 100, 950, 1500, 600, 400), 
                                                           genes.p = c(0.1, 0.03, 0.05, 0.02, 0.03, 0.2, 0.15), 
                                                           genes.mean.FC = c(0.3, 0.1, 2.6, 0.8, 0.4, 0.07, 3.2), 
                                                           genes.sd.FC = c(0.07, 0.13, 0.1, 0.02, 0.5, 0.4, 0.4))), 
  ### 10,000 cells
  tar_target(sim_lung_ncell10000_nclust3, simulate_clusters(ref.data = lung_ref, 
                                                            clust.n = c(1000, 7000, 2000), 
                                                            genes.p = c(0, 0.05, 0.1), 
                                                            genes.mean.FC = c(0, 0.17, 1.05), 
                                                            genes.sd.FC = c(0, 0.2, 0.6))), 
  tar_target(sim_lung_ncell10000_nclust5, simulate_clusters(ref.data = lung_ref, 
                                                            clust.n = c(300, 700, 5500, 2500, 1000), 
                                                            genes.p = c(0, 0.1, 0.07, 0.13, 0.3), 
                                                            genes.mean.FC = c(0, 0.95, 0.4, 2.6, 1.5), 
                                                            genes.sd.FC = c(0, 0.5, 1.4, 0.4, 0.05))), 
  tar_target(sim_lung_ncell10000_nclust7, simulate_clusters(ref.data = lung_ref, 
                                                            clust.n = c(75, 525, 400, 2700, 1300, 3500, 1500), 
                                                            genes.p = c(0.03, 0.5, 0.34, 0.11, 0.03, 0.15, 0.19), 
                                                            genes.mean.FC = c(1.1, 1.57, 0.17, 2.0, 1.1, 0.75, 0.75), 
                                                            genes.sd.FC = c(0.4, 0.9, 1.1, 1.2, 0.5, 0.1, 0.4))), 
  
  ##### clustering evaluation #####
  
  ### pancreas sim datasets
  # 1,000 cells
  tar_target(clustres_panc_ncell1000_nclust3, evaluate_clustering_all(sim.data = sim_panc_ncell1000_nclust3)), 
  tar_target(clustres_panc_ncell1000_nclust5, evaluate_clustering_all(sim.data = sim_panc_ncell1000_nclust5)), 
  tar_target(clustres_panc_ncell1000_nclust7, evaluate_clustering_all(sim.data = sim_panc_ncell1000_nclust7)), 
  # 3,000 cells
  tar_target(clustres_panc_ncell3000_nclust3, evaluate_clustering_all(sim.data = sim_panc_ncell3000_nclust3)), 
  tar_target(clustres_panc_ncell3000_nclust5, evaluate_clustering_all(sim.data = sim_panc_ncell3000_nclust5)), 
  tar_target(clustres_panc_ncell3000_nclust7, evaluate_clustering_all(sim.data = sim_panc_ncell3000_nclust7)), 
  # 5,000 cells
  tar_target(clustres_panc_ncell5000_nclust3, evaluate_clustering_all(sim.data = sim_panc_ncell5000_nclust3)), 
  tar_target(clustres_panc_ncell5000_nclust5, evaluate_clustering_all(sim.data = sim_panc_ncell5000_nclust5)), 
  tar_target(clustres_panc_ncell5000_nclust7, evaluate_clustering_all(sim.data = sim_panc_ncell5000_nclust7)), 
  # 10,000 cells 
  tar_target(clustres_panc_ncell10000_nclust3, evaluate_clustering_all(sim.data = sim_panc_ncell10000_nclust3)), 
  tar_target(clustres_panc_ncell10000_nclust5, evaluate_clustering_all(sim.data = sim_panc_ncell10000_nclust5)), 
  tar_target(clustres_panc_ncell10000_nclust7, evaluate_clustering_all(sim.data = sim_panc_ncell10000_nclust7)), 
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
  # 10,000 cells 
  tar_target(clustres_lung_ncell10000_nclust3, evaluate_clustering_all(sim.data = sim_lung_ncell10000_nclust3)), 
  tar_target(clustres_lung_ncell10000_nclust5, evaluate_clustering_all(sim.data = sim_lung_ncell10000_nclust5)), 
  tar_target(clustres_lung_ncell10000_nclust7, evaluate_clustering_all(sim.data = sim_lung_ncell10000_nclust7)), 
  
  ##### summary & analysis documents
  tar_render(simulation_qc, "./Reports/SCISSORS_Simulated_Data_QC.Rmd"), 
  tar_render(simulation_summary, "./Reports/SCISSORS_Simulation_Summary.Rmd")
)
