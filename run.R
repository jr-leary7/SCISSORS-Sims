setwd("/pine/scr/j/r/jrleary/R_Projects/SCISSORS-Sims")

library(targets)
library(tarchetypes)

# targets::tar_make()
targets::tar_make_future(workers = 4)
