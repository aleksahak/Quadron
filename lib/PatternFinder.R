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
PatternFinder <- function(seq=seq, str.pattern="([C]{3}[NATGCU]{1,12}){3,}[C]{3}"){

  all.strands  <- gregexpr(text = seq, pattern = str.pattern)

  QP.RESULTS <- NULL
  QP.RESULTS$start.pos  <- as.vector(all.strands[[1]])
  QP.RESULTS$seq.length <- attr(all.strands[[1]], "match.length")

  # Checking whether there are 0 returns (-1 by gregexpr)
  rm.ind <- which(QP.RESULTS$start.pos==-1)
  if(length(rm.ind)!=0){
    QP.RESULTS$start.pos  <- QP.RESULTS$start.pos[-rm.ind]
    QP.RESULTS$seq.length <- QP.RESULTS$seq.length[-rm.ind]
  }

  if(length(QP.RESULTS$start.pos)!=0){ # there ARE detected occurrences

    # sorting the results in the order of their start.pos:
    new.order  <- order(QP.RESULTS$start.pos)
    QP.RESULTS$start.pos  <- QP.RESULTS$start.pos[new.order]
    QP.RESULTS$seq.length <- QP.RESULTS$seq.length[new.order]

    QP.RESULTS$sequence  <- sapply(1:length(QP.RESULTS$start.pos),
                                   FUN=function(i){
                                     substr(
                                       seq,
                                       start=QP.RESULTS$start.pos[i],
                       stop=(QP.RESULTS$start.pos[i]+QP.RESULTS$seq.length[i]-1)
                                     )
                                   }, simplify=TRUE, USE.NAMES=FALSE)

    QP.RESULTS$strand <- sapply(QP.RESULTS$sequence,
                                FUN=function(i){
                                  if( substr(i, start=1, stop=1)=="G"){
                                    return("+")
                                  } else {
                                    return("-")
                                  }
                                }, simplify=TRUE, USE.NAMES=FALSE)

  } else { # there are NO detected occurrences
    QP.RESULTS$start.pos     <- NULL
    QP.RESULTS$seq.length    <- NULL
    QP.RESULTS$strand        <- NULL
    QP.RESULTS$sequence      <- NULL
  }

  #return( QP.RESULTS )
  return( as.data.frame(QP.RESULTS, stringsAsFactors=FALSE) )

}
## FUNCTION ####################################################################
