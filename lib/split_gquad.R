################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
## This function takes a G-quadruplex with minimum G-3 (GGG) intervals, then  ##
## trimming the remaining loops, so that they don't start or end with Gs.     ##
## The function returns a vector of the trimmed loops. If one or many of the  ##
## loops are just another repeat of Gs, "" will be returned in place of that. ##
## REQUIRES: exclude.terminal.char() in exclude_terminal_char.R file!!!       ##
################################################################################
split.gquad <- function(seq="GGGGATGCATTTGGGAAAAAGGGGTGTTT", spl="GGG"){

  if(!exists("exclude.terminal.char")){
    stop("split.gquad: failed dependency: load exclude.terminal.char().")
  }
  splitseq <- unlist(strsplit(seq, spl))
  splitseq <- splitseq[-which(splitseq=="")]
  if(length(splitseq)==0){
      return(c("","",""))
      # stop("split.gquad: splitseq returned character(0).")
  } else {
    return(exclude.terminal.char(excl.char="G", seqvec=splitseq,
                                 max.length.repeats=12))
  }
  
}
################################################################################
