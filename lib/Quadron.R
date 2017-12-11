################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
# Requires the libraries "doMC", "foreach", "itertools", "xgboost" (0.4-4),    #
# "caret" and "plyr".                                                          #
# If not already installed in R, you can install those by typing:              #
# install.packages(c("doMC", "foreach", "itertools", "plyr", "caret"))         #
# Specific steps are needed to install the xgboost version 0.4-4, as detailed  #
# in the Readme file.                                                          #
# The default fastread==TRUE option in readfasta requires "data.table" library.#
################################################################################
Quadron <- function(FastaFile    = "test.fasta",
                    OutFile      = "out.txt",
                    nCPU         = 4,
                    parsed.seq   = "",
                    SeqPartitionBy = 1000000){
################################################################################

# These lines can become part of the accepted arguments, after the NC develop-
# ment is done.
ReturnOnlyNC = FALSE
NonCanonical = FALSE
# if NonCanonical==TRUE, will also require bioconductor package IRanges:
#> source("http://bioconductor.org/biocLite.R")
#> biocLite("IRanges")

  if(ReturnOnlyNC==TRUE & NonCanonical==FALSE){
    stop("Quadron: ReturnOnlyNC can be activated if NonCanonical is TRUE.")
  }

  suppressWarnings(suppressPackageStartupMessages(library(doMC)))
  suppressWarnings(suppressPackageStartupMessages(library(foreach)))
  suppressWarnings(suppressPackageStartupMessages(library(itertools)))
  registerDoMC(cores = nCPU)
  
  if(packageDescription("xgboost")$Version!="0.4-4"){
    stop("Quadron: our model is robust and reproducible with the\n       xgboost version 0.4-4. Install the specific version of xgboost\n       as described in the Quadron documentation.")
    # install.packages(
    #  "http://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_0.4-4.tar.gz",
    #  repos=NULL,
    #  type="source"
    # )
    #
    # https://support.rstudio.com/hc/en-us/articles/219949047-Installing-older-versions-of-packages
  }
  
  info <- INFOline(OUT=info, msg=
  "NOTE: *:)* Sequence-Based Prediction of DNA Quadruplex Structures *(:*",
  initial=TRUE)

  info <- INFOline(OUT=info, msg=
  paste("NOTE: Date - ", date(), sep=""))

  info <- INFOline(OUT=info, msg=
  "NOTE: Parsing the sequence...")
  if(parsed.seq==""){
    seq <- readfasta(filename=FastaFile, case="UPPER", split=FALSE, fastread=TRUE)
  } else {
    seq <- NULL
    seq$seq <- paste(parsed.seq, collapse="")
    seq$length <- nchar(seq$seq)
  }

  info <- INFOline(OUT=info, msg=
  paste("NOTE: The digested sequence is of ", seq$length, "-nt length.", sep=""))

  info <- INFOline(OUT=info, msg=
  "NOTE: Scanning the sequence for G4 motifs...")

  if(NonCanonical==TRUE){
    info <- INFOline(OUT=info, msg=
    "NOTE: Relaxing the criteria out of the canonical G4 scope...")

    suppressPackageStartupMessages(library(IRanges))
    #### FOREACH EXECUTION #########
    QP <- foreach(j=isplitVector(1:seq$length, chunks=ceiling(seq$length/SeqPartitionBy)),
                  .combine="rbind", .inorder=TRUE) %dopar% {
            rng <- range(j)
            # Adding 15nt tolerance range for retrieving complete G4s from split parts:
            tol.rng <- c( max(1,(rng[1]-15)), min((rng[2]+15),seq$length) )

            qp <- PatternFinder(seq=substr(seq$seq, start=tol.rng[1], stop=tol.rng[2]),
               str.pattern="(([G]{3}[NATGCU]{1,12}){3,}[G]{3})|(([C]{3}[NATGCU]{1,12}){3,}[C]{3})")

            qp.nc <- PatternFinder(seq=substr(seq$seq, start=tol.rng[1], stop=tol.rng[2]),
               str.pattern="(([G]{2,}[NATGCU]{1,12}){3,}[G]{2,})|(([C]{2,}[NATGCU]{1,12}){3,}[C]{2,})")

            # Identifying overlaps between QP.nc and QP ************************
            if(length(qp$start.pos)!=0){
             if(length(qp.nc$start.pos)!=0){

              query=RangedData( IRanges(start = qp.nc$start.pos,
                                        end = qp.nc$start.pos+qp.nc$seq.length-1),
                                space = qp.nc$strand,
                                query.ind = 1:length(qp.nc$start.pos) )

              subject=RangedData( IRanges(start = qp$start.pos,
                                          end = qp$start.pos+qp$seq.length-1),
                                  space = qp$strand,
                                  subject.ind = 1:length(qp$start.pos) )

              ol <- as.matrix(findOverlaps(query=query, subject=subject, type="any",
                                           maxgap=0L, minoverlap=1L, select="all"))

              ol <- cbind(query=query[ol[,"queryHits"], ]$query.ind,
                          subject=subject[ol[,"subjectHits"], ]$subject.ind)

              if(!is.na(ol[,"query"][1])){ # there is at least one overlap
                # eliminating rows in non-canonical qp.nc that overlap with qp
                qp.nc <- qp.nc[-ol[,"query"],]
              }

              #----------------------
              if(ReturnOnlyNC==TRUE){
                qp <- qp.nc
              } else {
                # merging qp with qp.nc with reordering via genomic coordinates
                qp <- rbind(qp, qp.nc)
                reordering  <- order(qp$start.pos)
                qp <- qp[reordering,]
              }
              #----------------------

             }
             # shifting the coordinate system as soon as at least qp is not NULL df
             qp$start.pos <- qp$start.pos+tol.rng[1]-1
            }
            #*******************************************************************

            return(qp)
          }
    #### FOREACH EXECUTION DONE ####

  } else {

    #### FOREACH EXECUTION #########
    QP <- foreach(j=isplitVector(1:seq$length, chunks=ceiling(seq$length/SeqPartitionBy)),
                  .combine="rbind", .inorder=TRUE) %dopar% {
            rng <- range(j)
            # Adding 15nt tolerance range for retrieving complete G4s from split parts:
            tol.rng <- c( max(1,(rng[1]-15)), min((rng[2]+15),seq$length) )
            qp <- PatternFinder(seq=substr(seq$seq, start=tol.rng[1], stop=tol.rng[2]),
            str.pattern="(([G]{3}[NATGCU]{1,12}){3,}[G]{3})|(([C]{3}[NATGCU]{1,12}){3,}[C]{3})")
            qp$start.pos <- qp$start.pos+tol.rng[1]-1
            return(qp)
          }
    #### FOREACH EXECUTION DONE ####

  }

  # will use only:
  # QP$sequence
  # QP$start.pos
  # QP$seq.length
  # QP$strand

  ##############################################################################
  QP.num.occ <- length(QP$start.pos)

  if(QP.num.occ!=0){ # hence there are PQSs found.

  info <- INFOline(OUT=info, msg=
  paste("NOTE: Extracting features using ", nCPU, " processing core(s).", sep=""))

  #### FOREACH EXECUTION #########
  # Creating the feature holding data frame (RES), which will correspond to each
  # found PQS. In case the flanks are incomplete because of terminal location of
  # G4, the corresponding row will be populated by NAs only.

  #RES=foreach(i = 1:pqs.len, .combine="rbind") %dopar% {
  #pqs.len=100  ### *** ### TEST LINE
  RES <- foreach(ch=isplitVector(1:QP.num.occ, chunks=ceiling(QP.num.occ/1000)),
                 .combine="rbind", .inorder=TRUE) %dopar% {

    count <- 1
    CHRES <- as.data.frame(matrix(NA, ncol=length(feature.names.for.df), nrow=length(ch)))
    names(CHRES) <- feature.names.for.df

    for(i in ch){
      #*# print(i, quote=F)

      # Getting the GENOMIC sequences for PQS and flanks. Note, that
      # here 5' and 3' assignment naming will be correct for only + strand.
      pqs.seq    <-  QP$sequence[i]
      flank5.seq <- substr( seq$seq,
                            start = (QP$start.pos[i]-50),
                            stop  = (QP$start.pos[i]-1) )
      flank3.seq <- substr( seq$seq,
                            start = (QP$start.pos[i]+QP$seq.length[i]),
                            stop  = (QP$start.pos[i]+QP$seq.length[i]+49) )

      #-- Continuing only if both flanks are fully present.
      if(nchar(flank5.seq)==50 & nchar(flank3.seq)==50){

        # Correcting the sequences for - strand G4s.
        if(QP$strand[i]=="-"){
          pqs.seq    <- get.revcomp(pqs.seq)
          plus.flank5.seq  <- flank5.seq
          flank5.seq <- get.revcomp(flank3.seq)
          flank3.seq <- get.revcomp(plus.flank5.seq)
          rm(plus.flank5.seq)
        }

        # Extracting the required features.
        if(NonCanonical==TRUE){
          dfr <- FeatureExtractorQr(seq=pqs.seq, type="PQS", O=TRUE, OT=TRUE, G4=TRUE, NC=TRUE)
        } else {
          dfr <- FeatureExtractorQr(seq=pqs.seq, type="PQS", O=TRUE, OT=TRUE, G4=TRUE, NC=FALSE)
        }
        dfr <- data.frame(dfr)

        df.f5 <- FeatureExtractorQr(seq=flank5.seq, type="FLANK5", O=TRUE, OT=TRUE, G4=FALSE)
        df.f5 <- data.frame(df.f5)
        #*# names(df.f5) <- paste("FLANK5.",names(df.f5),sep="")

        df.f3 <- FeatureExtractorQr(seq=flank3.seq, type="FLANK3", O=TRUE, OT=TRUE, G4=FALSE)
        df.f3 <- data.frame(df.f3)
        #*# names(df.f3) <- paste("FLANK3.",names(df.f3),sep="")

        dfr <- cbind(dfr, df.f5, df.f3)
        CHRES[count,] <- dfr
      }
      #-- End of "continuing only if both flanks are fully present".

      count <- count+1
    }

    return(CHRES)

  }
  #### FOREACH EXECUTION DONE ####

    if(dim(RES)[1]==1 & all(is.na(RES[1,]))==TRUE){
      info <- INFOline(OUT=info, msg=
      "NOTE: There is only a single motif detected, with insufficient flanks.")
      info <- INFOline(OUT=info, msg=
      "NOTE: Quadron core will not be executed.")
      PRED <- NA

    } else {

      info <- INFOline(OUT=info, msg=
      "NOTE: Pre-processing the extracted features...")
      # load("./ModelProc.env")
      RES[,feature.names.for.df] <- t( (t(RES[,feature.names.for.df]) - medians[feature.names.for.df]) / sdevs[feature.names.for.df] )

      #info <- INFOline(OUT=info, msg=
      #"NOTE: Loading Quadron core...")
      #load("./QuadronML")
      suppressPackageStartupMessages(library(xgboost))
      suppressPackageStartupMessages(library(plyr))

      info <- INFOline(OUT=info, msg=
      "NOTE: Executing the Quadron core...")
      PRED <- rep(NA, QP.num.occ)
      non.na.rows <- which(!is.na(RES[,1]))
      PRED[non.na.rows] <- predict(QuadronML, newdata=RES[non.na.rows,])

    }

  } else { # no PQS is found
    PRED <- NA
    QP   <- NULL
    QP$start.pos     <- NA
    QP$seq.length    <- NA
    QP$strand        <- NA
    QP$sequence      <- NA
  }
  ##############################################################################

  info <- INFOline(OUT=info, msg=
  "NOTE: Formatting and saving the results...")

  write(info, file=OutFile)

  out.start <- getSTART()
  write(out.start, file=OutFile, append=TRUE)

  out.data <- paste("DATA:",
                    format(QP$start.pos),
                    format(QP$strand),
                    format(QP$seq.length),
                    format(round(PRED,2), nsmal=2),
                    QP$sequence,
                    sep=" ")

  write(out.data, file=OutFile, append=TRUE)

  out.end <- getEND()
  write(out.end, file=OutFile, append=TRUE)

  print("NOTE: Quadron is done!", quote=FALSE)
  return(c(info, out.end))

}
################################################################################
