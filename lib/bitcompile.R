################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
LibPath      = "."
################################################################################

suppressWarnings(suppressPackageStartupMessages(library(compiler)))
source(paste(LibPath,"/INFOline.R",sep=""))
source(paste(LibPath,"/readfasta.R",sep=""))
source(paste(LibPath,"/PatternFinder.R",sep=""))
source(paste(LibPath,"/get_revcomp.R",sep=""))
source(paste(LibPath,"/linesplit.R",sep=""))
source(paste(LibPath,"/exclude_terminal_char.R",sep=""))
source(paste(LibPath,"/split_gquad.R",sep=""))
source(paste(LibPath,"/FeatureExtractorQr.R",sep=""))
source(paste(LibPath,"/getSTART.R",sep=""))
source(paste(LibPath,"/getEND.R",sep=""))
source(paste(LibPath,"/Quadron.R",sep=""))

INFOline           <- cmpfun(INFOline,     options=list(suppressUndefined=TRUE))
readfasta          <- cmpfun(readfasta,    options=list(suppressUndefined=TRUE))
PatternFinder      <- cmpfun(PatternFinder, options=list(suppressUndefined=TRUE))
get.revcomp        <- cmpfun(get.revcomp,  options=list(suppressUndefined=TRUE))
linesplit          <- cmpfun(linesplit,    options=list(suppressUndefined=TRUE))
exclude.terminal.char <- cmpfun(exclude.terminal.char, options=list(suppressUndefined=TRUE))
split.gquad        <- cmpfun(split.gquad, options=list(suppressUndefined=TRUE))
FeatureExtractorQr <- cmpfun(FeatureExtractorQr, options=list(suppressUndefined=TRUE))
getSTART           <- cmpfun(getSTART,     options=list(suppressUndefined=TRUE))
getEND             <- cmpfun(getEND,       options=list(suppressUndefined=TRUE))
Quadron            <- cmpfun(Quadron,      options=list(suppressUndefined=TRUE))

feature.names.for.df <- c('O.Acont', 'O.Ccont', 'O.Gcont', 'O.Length',
                          'OT.triad.AAG', 'OT.triad.AGA', 'OT.triad.AGC',
                          'OT.triad.AGG', 'OT.triad.ATA', 'OT.triad.CAG',
                          'OT.triad.CCC', 'OT.triad.CGG', 'OT.triad.CTG',
                          'OT.triad.GAA', 'OT.triad.GAG', 'OT.triad.GCA',
                          'OT.triad.GCC', 'OT.triad.GCG', 'OT.triad.GCT',
                          'OT.triad.GGA', 'OT.triad.GGC', 'OT.triad.GGG',
                          'OT.triad.GGT', 'OT.triad.GTG', 'OT.triad.GTT',
                          'OT.triad.TGC', 'OT.triad.TGG', 'OT.triad.TGT',
                          'OT.triad.TTG',
                          'G4.Loop1.efe', 'G4.Loop2.efe', 'G4.Loop3.efe',
                          'G4.nonG3plus.Loop1.Length',
                          'G4.nonG3plus.Loop2.Length',
                          'G4.nonG3plus.Loop3.Length',
                          'G4.nonG3plus.Nloops',
                          'FLANK5.O.Acont', 'FLANK5.O.Ccont', 'FLANK5.O.Gcont',
                          'FLANK5.OT.triad.AAA', 'FLANK5.OT.triad.AAG', 'FLANK5.OT.triad.AAT',
                          'FLANK5.OT.triad.ACA', 'FLANK5.OT.triad.ACC', 'FLANK5.OT.triad.AGA',
                          'FLANK5.OT.triad.AGC', 'FLANK5.OT.triad.AGG', 'FLANK5.OT.triad.ATG',
                          'FLANK5.OT.triad.CAA', 'FLANK5.OT.triad.CAC', 'FLANK5.OT.triad.CAG',
                          'FLANK5.OT.triad.CCA', 'FLANK5.OT.triad.CCC', 'FLANK5.OT.triad.CCG',
                          'FLANK5.OT.triad.CCT', 'FLANK5.OT.triad.CGC', 'FLANK5.OT.triad.CGG',
                          'FLANK5.OT.triad.CTC', 'FLANK5.OT.triad.CTG', 'FLANK5.OT.triad.CTT',
                          'FLANK5.OT.triad.GAA', 'FLANK5.OT.triad.GAG', 'FLANK5.OT.triad.GCA',
                          'FLANK5.OT.triad.GCC', 'FLANK5.OT.triad.GCG', 'FLANK5.OT.triad.GCT',
                          'FLANK5.OT.triad.GGA', 'FLANK5.OT.triad.GGC', 'FLANK5.OT.triad.GGG',
                          'FLANK5.OT.triad.GGT', 'FLANK5.OT.triad.GTG', 'FLANK5.OT.triad.TCA',
                          'FLANK5.OT.triad.TCC', 'FLANK5.OT.triad.TCT', 'FLANK5.OT.triad.TGA',
                          'FLANK5.OT.triad.TGC', 'FLANK5.OT.triad.TGG', 'FLANK5.OT.triad.TGT',
                          'FLANK5.OT.triad.TTT',
                          'FLANK3.O.Acont', 'FLANK3.O.Ccont', 'FLANK3.O.Gcont',
                          'FLANK3.OT.triad.AAA', 'FLANK3.OT.triad.AAG', 'FLANK3.OT.triad.ACA',
                          'FLANK3.OT.triad.AGA', 'FLANK3.OT.triad.AGC', 'FLANK3.OT.triad.AGG',
                          'FLANK3.OT.triad.AGT', 'FLANK3.OT.triad.ATA', 'FLANK3.OT.triad.ATG',
                          'FLANK3.OT.triad.CAG', 'FLANK3.OT.triad.CCA', 'FLANK3.OT.triad.CCC',
                          'FLANK3.OT.triad.CCG', 'FLANK3.OT.triad.CCT', 'FLANK3.OT.triad.CGC',
                          'FLANK3.OT.triad.CTG', 'FLANK3.OT.triad.GAA', 'FLANK3.OT.triad.GAG',
                          'FLANK3.OT.triad.GAT', 'FLANK3.OT.triad.GCA', 'FLANK3.OT.triad.GCC',
                          'FLANK3.OT.triad.GCG', 'FLANK3.OT.triad.GCT', 'FLANK3.OT.triad.GGA',
                          'FLANK3.OT.triad.GGC', 'FLANK3.OT.triad.GGG', 'FLANK3.OT.triad.GGT',
                          'FLANK3.OT.triad.GTG', 'FLANK3.OT.triad.TAC', 'FLANK3.OT.triad.TAG',
                          'FLANK3.OT.triad.TAT', 'FLANK3.OT.triad.TCT', 'FLANK3.OT.triad.TGA',
                          'FLANK3.OT.triad.TGC', 'FLANK3.OT.triad.TGG', 'FLANK3.OT.triad.TGT',
                          'FLANK3.OT.triad.TTT')

dna.loopfold <- readLines(paste(LibPath,"/dna_loopfold.txt",sep=""))
load(paste(LibPath,"/ModelProc.env",sep=""))

system(paste("cat ",LibPath,"/QuadronML.* > ",LibPath,"/QuadronML",sep=""))
load(paste(LibPath,"/QuadronML",sep=""))
file.remove(paste(LibPath,"/QuadronML",sep=""))

save(list=c("INFOline",
            "readfasta",
            "PatternFinder",
            "get.revcomp",
            "linesplit",
            "exclude.terminal.char",
            "split.gquad",
            "FeatureExtractorQr",
            "getSTART",
            "getEND",
            "Quadron",
            "feature.names.for.df",
            "dna.loopfold",
            "FeatureNames",
            "medians",
            "sdevs",
            "QuadronML"), file="../Quadron.lib")

#file.copy(from=paste(LibPath,"/QuadronML",sep=""),
#          to="../QuadronML")
################################################################################
