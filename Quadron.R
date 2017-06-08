################################################################################
# Requires the libraries "doMC", "foreach" and "itertools".                    #
# If not already installed in R, you can install those by typing:              #
# install.packages(c("doMC", "foreach" and "itertools"))                       #
################################################################################
#setwd("./lib")
#source("bitcompile.R")
#rm(list=ls())
#setwd("..")

print("NOTE: Loading Quadron core...", quote=FALSE)
load("Quadron.lib")

Quadron(FastaFile= "test.fasta", 
        OutFile  = "out.txt",
        nCPU     = 4,
        NonCanonical = TRUE,
        SeqPartitionBy = 1000000,
        ReturnOnlyNC = FALSE)
        
#file.remove("Quadron.lib")
#rm(list=ls())
################################################################################


