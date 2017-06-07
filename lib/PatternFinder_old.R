################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
# This function takes a DNA sequence and looks for the requested motifs.       #
#                                                                              #
# seq - a string of characters. The merged sequence in capital letters should  #
#       be provided, such as "ATGG...". The allowed characters are N, A, T, G, #
#       C, U.                                                                  #
#                                                                              #
################################################################################
PatternFinder <- function(seq=seq,
                          plus.strand.pattern="([G]{3}[NATGCU]{1,12}){3,}[G]{3}",
                          minus.strand.pattern="([C]{3}[NATGCU]{1,12}){3,}[C]{3}"){

  plus.strand  <- gregexpr(text = seq, pattern = plus.strand.pattern)
  minus.strand <- gregexpr(text = seq, pattern = minus.strand.pattern)

  start.pos  <- c( as.vector(plus.strand[[1]]), as.vector(minus.strand[[1]]) )
  seq.length <- c( attr(plus.strand[[1]], "match.length"), attr(minus.strand[[1]], "match.length") )
  strand     <- c( rep("+",length(plus.strand[[1]])), rep("-",length(minus.strand[[1]])) )

  # Checking whether there are 0 returns (-1 by gregexpr)
  rm.ind <- which(start.pos==-1)
  if(length(rm.ind)!=0){
    start.pos  <- start.pos[-rm.ind]
    seq.length <- seq.length[-rm.ind]
    strand     <- strand[-rm.ind]
  }

  if(length(start.pos)!=0){ # there ARE detected occurrences

    # sorting the results in the order of their start.pos:
    new.order  <- order(start.pos)

    start.pos  <- start.pos[new.order]
    seq.length <- seq.length[new.order]
    strand     <- strand[new.order]
    #num.occ    <- length(start.pos)
    #num.occ.plus  <- length(which(strand=="+"))
    #num.occ.minus <- length(which(strand=="-"))
    sequence  <- sapply(1:length(start.pos), FUN=function(i){
                        substr(seq, start=start.pos[i], stop=(start.pos[i]+seq.length[i]-1) )
                       }, simplify=TRUE, USE.NAMES=FALSE)

  } else { # there are NO detected occurrences
    start.pos     <- 0
    seq.length    <- 0
    strand        <- 0
    #num.occ       <- 0
    #num.occ.plus  <- 0
    #num.occ.minus <- 0
    sequence      <- 0
  }

  QP.RESULTS <- NULL
  QP.RESULTS$start.pos     <- start.pos
  QP.RESULTS$seq.length    <- seq.length
  QP.RESULTS$strand        <- strand
  QP.RESULTS$sequence      <- sequence
  #QP.RESULTS$num.occ       <- num.occ
  #QP.RESULTS$num.occ.plus  <- num.occ.plus
  #QP.RESULTS$num.occ.minus <- num.occ.minus

  #return( QP.RESULTS )
  return( as.data.frame(QP.RESULTS, stringsAsFactors=FALSE) )

}
## FUNCTION ####################################################################
