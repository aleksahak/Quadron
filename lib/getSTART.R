################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
getSTART <- function(){

  return( c(
    "HEADER: **********************************************************************",
    "HEADER: POS - position of the genomic start for the PQS_L12 potential",
    "HEADER:       quadruplex sequence with canonical G3+ tracks and maximum",
    "HEADER:       loop size of 12 nt. The position corresponds to 5'-start for",
    "HEADER:       the motifs in + strand and to 3'-end for the ones in - strand.",
    "HEADER: STR - strand where the PQS_L12 is located. The supplied sequence is",
    "HEADER:       treated as +.",
    "HEADER: L   - length of the retrieved PQS_L12 motif.",
    "HEADER: Q   - Quadron prediction of the corresponding G4-seq mismatch level",
    "HEADER:       for a polymerase stalling at quadruplex sites. NA indicates that",
    "HEADER:       the PQS_L12 is to close to the sequence termini for the 50-nt",
    "HEADER:       flanks to be analysed, as required for Quadron predictions.",
    "HEADER: SEQUENCE - PQS_L12 sequence in the supplied + strand. If, for a ",
    "HEADER:       particular PQS_L12, str is +, then this is the PQS_L12 sequence",
    "HEADER:       in 5'-3' direction. Otherwise, for - strand motifs, SEQUENCE",
    "HEADER:       would indicate the reverse complementar of the actual quadruplex",
    "HEADER:       sequence (i.e. complementar to its 3'-5'span).",
    "HEADER: **********************************************************************",
    "HEADER: POS STR L    Q SEQUENCE",
    "HEADER: **********************************************************************")
  )

}
################################################################################
