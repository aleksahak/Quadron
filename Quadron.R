################################################################################
# Requires the libraries "doMC", "foreach" and "itertools".                    #
# If not already installed in R, you can install those by typing:              #
# install.packages(c("doMC", "foreach" and "itertools"))                       #
################################################################################
print("NOTE: Loading Quadron core...", quote=FALSE)
load("Quadron.lib")

Quadron(FastaFile= "test.fasta", 
        OutFile  = "out.txt",
        nCPU     = 1)
################################################################################
