################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
# Requires split.gquad, exclude.terminal.char, linesplit                       #
# and the read in dna.loopfold object, in case G4 is TRUE.                     #
#                                                                              #
# Returns 3 types of features:                                                 #
# O  -- features that describe the general characteristics of the whole sequence
# OT -- triad count of the whole sequence                                      #
# G4 -- features that describe the engulfed quadruplex,                        #
#       assuming that the whole sequence is a quadruplex                       #
################################################################################

FeatureExtractorQr <- function(seq="GGGGTGGTGGGCCCGCGGGGCGCGGGGAAGCGGGGGAGGAGGGCAGGG",
                               O=TRUE, OT=TRUE, G4=TRUE,
                               type="PQS", # "FLANK5", "FLANK3"
                               NC = FALSE
                              ){

  # O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O
  if(O==TRUE){
    seq.split <- strsplit(seq, split="")[[1]]
    O.Length  <- length(seq.split)            #--$
    O.Gcont   <- sum(seq.split=="G")/O.Length #--$
    O.Ccont   <- sum(seq.split=="C")/O.Length #--$
    O.Acont   <- sum(seq.split=="A")/O.Length #--$
  }
  # O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O O

  # OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT
  if(OT==TRUE){
    kmers <- table(substring( seq, first=1:O.Length, last=((1:O.Length)+3-1) ))

    if(type=="PQS"){
      kmer.names <- c("AAG", "AGA", "AGC", "AGG", "ATA", "CAG", "CCC", "CGG", "CTG",
                      "GAA", "GAG", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG",
                      "GGT", "GTG", "GTT", "TGC", "TGG", "TGT", "TTG")
    }
    if(type=="FLANK5"){
      kmer.names <- c("AAA", "AAG", "AAT", "ACA", "ACC", "AGA", "AGC", "AGG", "ATG",
                      "CAA", "CAC", "CAG", "CCA", "CCC", "CCG", "CCT", "CGC", "CGG",
                      "CTC", "CTG", "CTT", "GAA", "GAG", "GCA", "GCC", "GCG", "GCT",
                      "GGA", "GGC", "GGG", "GGT", "GTG", "TCA", "TCC", "TCT", "TGA",
                      "TGC", "TGG", "TGT", "TTT")
      rm(O.Length)
    }
    if(type=="FLANK3"){
      kmer.names <- c("AAA", "AAG", "ACA", "AGA", "AGC", "AGG", "AGT", "ATA", "ATG",
                      "CAG", "CCA", "CCC", "CCG", "CCT", "CGC", "CTG", "GAA", "GAG",
                      "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT",
                      "GTG", "TAC", "TAG", "TAT", "TCT", "TGA", "TGC", "TGG", "TGT",
                      "TTT")
      rm(O.Length)
    }

    # Reordering kmers in kmer.names order, producing NAs for triads absent in kmers.
    kmers <- as.vector(kmers[match(kmer.names, names(kmers))])
    # Replacing NAs with 0s.
    kmers[is.na(kmers)] <- 0

    for(k in 1:length(kmer.names)){
       eval(parse(text=paste("OT.triad.",kmer.names[k]," <- kmers[",k,"]",sep="")))
    }
  }
  # OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT OT

  # G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4
  if(G4==TRUE){
    # a vector where each entry is a string of a loop sequence
    if(NC==TRUE){
      nonG3plus.loops     <- split.gquad(seq, spl="GG") # that would be non G2 loops
    } else {
      nonG3plus.loops     <- split.gquad(seq, spl="GGG")
    }

    #*#- # the loop merged and split vector of individual characters in loops
    #*#- nonG3plus.seq.split <- unlist(strsplit(nonG3plus.loops, split=""))
    #*#- G4.nonG3plus.Length <- length(nonG3plus.seq.split) # total length of the loops

    G4.nonG3plus.Nloops <- length(nonG3plus.loops)     # number of loops

    #F7.nonG3plus.Nloops.with.G2 <- length(grep("GG", nonG3plus.loops))
    #F8.nonG3plus.Nloops.with.C3 <- length(grep("CCC", nonG3plus.loops))

    Loop123.lengths <- nchar(nonG3plus.loops)[1:3]
    Loop123.lengths[is.na(Loop123.lengths)] <- 0

    G4.nonG3plus.Loop1.Length <- Loop123.lengths[1]
    G4.nonG3plus.Loop2.Length <- Loop123.lengths[2]
    G4.nonG3plus.Loop3.Length <- Loop123.lengths[3]

    G4.Loop1.efe <- 0
    if(G4.nonG3plus.Loop1.Length>=5 & G4.nonG3plus.Loop1.Length<=12){
      pl <- which(dna.loopfold == nonG3plus.loops[1])
      if(length(pl)!=0){G4.Loop1.efe <- as.numeric(linesplit(dna.loopfold[pl+2])[6])}
    }

    G4.Loop2.efe <- 0
    if(G4.nonG3plus.Loop2.Length>=5 & G4.nonG3plus.Loop2.Length<=12){
      pl <- which(dna.loopfold == nonG3plus.loops[2])
      if(length(pl)!=0){G4.Loop2.efe <- as.numeric(linesplit(dna.loopfold[pl+2])[6])}
    }

    G4.Loop3.efe <- 0
    if(G4.nonG3plus.Loop3.Length>=5 & G4.nonG3plus.Loop3.Length<=12){
      pl <- which(dna.loopfold == nonG3plus.loops[3])
      if(length(pl)!=0){G4.Loop3.efe <- as.numeric(linesplit(dna.loopfold[pl+2])[6])}
    }

  }
  # G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4 G4

  RESULT <- NULL
  #Embeded UTIL_CollectObjects code:
  #-------------------
  for( pf in c("O.","OT.","G4.")[c(O,OT,G4)] ){

    compon   <- grep(pf, ls(), value=TRUE, fixed=TRUE) # vector of character(0)
    if(length(compon) != 0){
      for(h in compon){
        eval(parse(text=paste("RESULT$",h,"<-",h,sep="")))
      }
    }

  }
  #-------------------

  return( RESULT )
  # THE ORDER IN THE RESULT SHOULD ALWAYS BE KEPT THE SAME!

}
################################################################################
# # 119 FEATURES
# 1- "O.Length"
# 3- c("O.Acont", "O.Ccont", "O.Gcont")
# 25- c("OT.triad.AAG", "OT.triad.AGA", "OT.triad.AGC", "OT.triad.AGG", "OT.triad.ATA", "OT.triad.CAG", "OT.triad.CCC", "OT.triad.CGG", "OT.triad.CTG", "OT.triad.GAA", "OT.triad.GAG", "OT.triad.GCA", "OT.triad.GCC", "OT.triad.GCG", "OT.triad.GCT", "OT.triad.GGA", "OT.triad.GGC", "OT.triad.GGG", "OT.triad.GGT", "OT.triad.GTG", "OT.triad.GTT", "OT.triad.TGC", "OT.triad.TGG", "OT.triad.TGT", "OT.triad.TTG")
# 3- c("G4.Loop1.efe", "G4.Loop2.efe", "G4.Loop3.efe")
# 3- c("G4.nonG3plus.Loop1.Length", "G4.nonG3plus.Loop2.Length", "G4.nonG3plus.Loop3.Length")
# 1- "G4.nonG3plus.Nloops"
# 3- c("FLANK5.O.Acont", "FLANK5.O.Ccont", "FLANK5.O.Gcont")
# 40- c("FLANK5.OT.triad.AAA", "FLANK5.OT.triad.AAG", "FLANK5.OT.triad.AAT", "FLANK5.OT.triad.ACA", "FLANK5.OT.triad.ACC", "FLANK5.OT.triad.AGA", "FLANK5.OT.triad.AGC", "FLANK5.OT.triad.AGG", "FLANK5.OT.triad.ATG", "FLANK5.OT.triad.CAA", "FLANK5.OT.triad.CAC", "FLANK5.OT.triad.CAG", "FLANK5.OT.triad.CCA", "FLANK5.OT.triad.CCC", "FLANK5.OT.triad.CCG", "FLANK5.OT.triad.CCT", "FLANK5.OT.triad.CGC", "FLANK5.OT.triad.CGG", "FLANK5.OT.triad.CTC", "FLANK5.OT.triad.CTG", "FLANK5.OT.triad.CTT", "FLANK5.OT.triad.GAA", "FLANK5.OT.triad.GAG", "FLANK5.OT.triad.GCA", "FLANK5.OT.triad.GCC", "FLANK5.OT.triad.GCG", "FLANK5.OT.triad.GCT", "FLANK5.OT.triad.GGA", "FLANK5.OT.triad.GGC", "FLANK5.OT.triad.GGG", "FLANK5.OT.triad.GGT", "FLANK5.OT.triad.GTG", "FLANK5.OT.triad.TCA", "FLANK5.OT.triad.TCC", "FLANK5.OT.triad.TCT", "FLANK5.OT.triad.TGA", "FLANK5.OT.triad.TGC", "FLANK5.OT.triad.TGG", "FLANK5.OT.triad.TGT", "FLANK5.OT.triad.TTT")
# 3- c("FLANK3.O.Acont", "FLANK3.O.Ccont", "FLANK3.O.Gcont")
# 37- c("FLANK3.OT.triad.AAA", "FLANK3.OT.triad.AAG", "FLANK3.OT.triad.ACA", "FLANK3.OT.triad.AGA", "FLANK3.OT.triad.AGC", "FLANK3.OT.triad.AGG", "FLANK3.OT.triad.AGT", "FLANK3.OT.triad.ATA", "FLANK3.OT.triad.ATG", "FLANK3.OT.triad.CAG", "FLANK3.OT.triad.CCA", "FLANK3.OT.triad.CCC", "FLANK3.OT.triad.CCG", "FLANK3.OT.triad.CCT", "FLANK3.OT.triad.CGC", "FLANK3.OT.triad.CTG", "FLANK3.OT.triad.GAA", "FLANK3.OT.triad.GAG", "FLANK3.OT.triad.GAT", "FLANK3.OT.triad.GCA", "FLANK3.OT.triad.GCC", "FLANK3.OT.triad.GCG", "FLANK3.OT.triad.GCT", "FLANK3.OT.triad.GGA", "FLANK3.OT.triad.GGC", "FLANK3.OT.triad.GGG", "FLANK3.OT.triad.GGT", "FLANK3.OT.triad.GTG", "FLANK3.OT.triad.TAC", "FLANK3.OT.triad.TAG", "FLANK3.OT.triad.TAT", "FLANK3.OT.triad.TCT", "FLANK3.OT.triad.TGA", "FLANK3.OT.triad.TGC", "FLANK3.OT.triad.TGG", "FLANK3.OT.triad.TGT", "FLANK3.OT.triad.TTT")
#
# triad.namesG4 <- c("OT.triad.AAG", "OT.triad.AGA", "OT.triad.AGC", "OT.triad.AGG",
#                    "OT.triad.ATA", "OT.triad.CAG", "OT.triad.CCC", "OT.triad.CGG",
#                    "OT.triad.CTG", "OT.triad.GAA", "OT.triad.GAG", "OT.triad.GCA",
#                    "OT.triad.GCC", "OT.triad.GCG", "OT.triad.GCT", "OT.triad.GGA",
#                    "OT.triad.GGC", "OT.triad.GGG", "OT.triad.GGT", "OT.triad.GTG",
#                    "OT.triad.GTT", "OT.triad.TGC", "OT.triad.TGG", "OT.triad.TGT",
#                    "OT.triad.TTG")
# triad.namesF5 <- c("FLANK5.OT.triad.AAA", "FLANK5.OT.triad.AAG", "FLANK5.OT.triad.AAT",
#                    "FLANK5.OT.triad.ACA", "FLANK5.OT.triad.ACC", "FLANK5.OT.triad.AGA",
#                    "FLANK5.OT.triad.AGC", "FLANK5.OT.triad.AGG", "FLANK5.OT.triad.ATG",
#                    "FLANK5.OT.triad.CAA", "FLANK5.OT.triad.CAC", "FLANK5.OT.triad.CAG",
#                    "FLANK5.OT.triad.CCA", "FLANK5.OT.triad.CCC", "FLANK5.OT.triad.CCG",
#                    "FLANK5.OT.triad.CCT", "FLANK5.OT.triad.CGC", "FLANK5.OT.triad.CGG",
#                    "FLANK5.OT.triad.CTC", "FLANK5.OT.triad.CTG", "FLANK5.OT.triad.CTT",
#                    "FLANK5.OT.triad.GAA", "FLANK5.OT.triad.GAG", "FLANK5.OT.triad.GCA",
#                    "FLANK5.OT.triad.GCC", "FLANK5.OT.triad.GCG", "FLANK5.OT.triad.GCT",
#                    "FLANK5.OT.triad.GGA", "FLANK5.OT.triad.GGC", "FLANK5.OT.triad.GGG",
#                    "FLANK5.OT.triad.GGT", "FLANK5.OT.triad.GTG", "FLANK5.OT.triad.TCA",
#                    "FLANK5.OT.triad.TCC", "FLANK5.OT.triad.TCT", "FLANK5.OT.triad.TGA",
#                    "FLANK5.OT.triad.TGC", "FLANK5.OT.triad.TGG", "FLANK5.OT.triad.TGT",
#                    "FLANK5.OT.triad.TTT")
# triad.namesF3 <- c("FLANK3.OT.triad.AAA", "FLANK3.OT.triad.AAG", "FLANK3.OT.triad.ACA",
#                    "FLANK3.OT.triad.AGA", "FLANK3.OT.triad.AGC", "FLANK3.OT.triad.AGG",
#                    "FLANK3.OT.triad.AGT", "FLANK3.OT.triad.ATA", "FLANK3.OT.triad.ATG",
#                    "FLANK3.OT.triad.CAG", "FLANK3.OT.triad.CCA", "FLANK3.OT.triad.CCC",
#                    "FLANK3.OT.triad.CCG", "FLANK3.OT.triad.CCT", "FLANK3.OT.triad.CGC",
#                    "FLANK3.OT.triad.CTG", "FLANK3.OT.triad.GAA", "FLANK3.OT.triad.GAG",
#                    "FLANK3.OT.triad.GAT", "FLANK3.OT.triad.GCA", "FLANK3.OT.triad.GCC",
#                    "FLANK3.OT.triad.GCG", "FLANK3.OT.triad.GCT", "FLANK3.OT.triad.GGA",
#                    "FLANK3.OT.triad.GGC", "FLANK3.OT.triad.GGG", "FLANK3.OT.triad.GGT",
#                    "FLANK3.OT.triad.GTG", "FLANK3.OT.triad.TAC", "FLANK3.OT.triad.TAG",
#                    "FLANK3.OT.triad.TAT", "FLANK3.OT.triad.TCT", "FLANK3.OT.triad.TGA",
#                    "FLANK3.OT.triad.TGC", "FLANK3.OT.triad.TGG", "FLANK3.OT.triad.TGT",
#                    "FLANK3.OT.triad.TTT")
# paste(sapply(triad.namesF3,
#              FUN=function(i){strsplit(i,".", fixed=TRUE)[[1]][4]},
#              USE.NAMES=FALSE, simplify=TRUE), collapse=", ")
# # Below is the complete set of 64 triads:
# kmer.names <- c( "AAA","AAC","AAG","AAT","ACA","ACC","ACG",
#                  "ACT","AGA","AGC","AGG","AGT","ATA","ATC",
#                  "ATG","ATT","CAA","CAC","CAG","CAT","CCA",
#                  "CCC","CCG","CCT","CGA","CGC","CGG","CGT",
#                  "CTA","CTC","CTG","CTT","GAA","GAC","GAG",
#                  "GAT","GCA","GCC","GCG","GCT","GGA","GGC",
#                  "GGG","GGT","GTA","GTC","GTG","GTT","TAA",
#                  "TAC","TAG","TAT","TCA","TCC","TCG","TCT",
#                  "TGA","TGC","TGG","TGT","TTA","TTC","TTG",
#                  "TTT" )
