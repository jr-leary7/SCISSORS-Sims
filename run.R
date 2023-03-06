setwd("/work/users/j/r/jrleary/scRNAseq/SCISSORS_Project/SCISSORS-Sims/")

library(targets)
library(tarchetypes)

targets::tar_make_future(workers = 4)

# command 
# sbatch -t 48:00:00 -c 2 --mem=200G -J SCISSORS_Sim  --mail-type=ALL --mail-user=jrleary@live.unc.edu --wrap="module load r; Rscript run.R"
# sbatch -t 48:00:00 -c 1 --mem=100G -J SCISSORS_Sim  --mail-type=ALL --mail-user=jrleary@live.unc.edu --wrap="module load r; Rscript run.R"
# monitoring 
# targets::tar_watch(level_separation = 1200, seconds = 120, seconds_max = 360)
