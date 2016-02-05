################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
#  This function takes a character line and decomposes it cutting the spaces.  #
################################################################################
linesplit <- function(line) {

  Line <- strsplit(line," ")[[1]]     # only compatibel with a single-line input
  #Line <- unlist(strsplit(line," ")) # could also parse/merge many lines, 
  #                                      depreciated since 21 Aug 2015.
  return( Line[which(Line!="")] )

}
################################################################################