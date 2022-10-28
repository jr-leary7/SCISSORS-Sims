# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

setwd("/pine/scr/j/r/jrleary/R_Projects/SCISSORS-Sims")
library(targets)
library(tarchetypes)

# targets::tar_make()
# tar_make_clustermq(workers = 3) # nolint
targets::tar_make_future(workers = 3) # nolint
