################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
## This function takes a vector of characters and/or strings and trims those  ##
## from both sides, excluding all the terminal excl.char repeats (if aplcbl.).##
## The function returns a vector of the trimmed results, where some strings   ##
## can be fully untouched because of not containing the specified excl.char.  ##
## For the input entrie(s) that are fully a repeat of excl.char, the function ##
## returns "". The max.length.repeats argument specifies the maximum length   ##
## of the repeats of excl.char that shoud be trimmed from either sides.       ##
################################################################################
exclude.terminal.char <- function(excl.char="G", max.length.repeats=100,
                                  seqvec=c("GGGATGCATTTGG","GGAAAA","TTT")){
    
    G.start.ind <- NULL
    
    for(iter in 1:max.length.repeats){
      starts <- substr(start=1, stop=1, x=seqvec)
      ends   <- substr(start=nchar(seqvec), stop=nchar(seqvec), x=seqvec)
      G.start.ind <- which(starts==excl.char)
      if(length(G.start.ind)!=0){
        seqvec[G.start.ind] <- substr(start=2, stop=100, x=seqvec[G.start.ind])
      } else {
        break
      }
    }
    
    for(iter in 1:max.length.repeats){
      ends   <- substr(start=nchar(seqvec), stop=nchar(seqvec), x=seqvec)
      G.end.ind <- which(ends==excl.char)
      if(length(G.end.ind)!=0){
        seqvec[G.end.ind] <- substr(start=1, stop=(nchar(seqvec[G.end.ind])-1),
        x=seqvec[G.end.ind])
      } else {
        break
      }
    }
    
    return(seqvec)
    
}
################################################################################