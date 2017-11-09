################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
getSTART <- function(){

  return( c(
    "HEADER: **********************************************************************",
    "HEADER: POS - position of the genomic start for the PQS_L12, i.e. potential",
    "HEADER:       quadruplex sequence with canonical G3+ tracts and maximum",
    "HEADER:       loop size of 12 nt. The genomic start position thus",
    "HEADER        corresponds to 5'-start for PQS_L12 sequence, if that",
    "HEADER:       is in the '+' strand, and to 3'-end for the one in '-' strand.",
    "HEADER: STR - strand where the PQS_L12 is located. The supplied sequence is",
    "HEADER:       treated as '+'.",
    "HEADER: L   - length of the retrieved PQS_L12 motif.",
    "HEADER: Q   - Quadron prediction of the corresponding G4-seq mismatch level",
    "HEADER:       for a polymerase stalling at quadruplex sites. NA indicates that",
    "HEADER:       the PQS_L12 is too close to the sequence termini for the 50-nt",
    "HEADER:       flanks to be analysed, as required for Quadron predictions.",
    "HEADER:       Q values above 19 indicate that the corresponding PQS_L12 is",
    "HEADER:       a highly stable G-quadruplex.",
    "HEADER: SEQUENCE - PQS_L12 sequence in the supplied + strand. If, for a ",
    "HEADER:       particular PQS_L12, str is '+', then this is the PQS_L12 sequence",
    "HEADER:       in 5'-3' direction. Otherwise, for '-' strand motifs, SEQUENCE",
    "HEADER:       would be reverse complementary to the actual PQS_L12 (i.e.",
    "HEADER:       complementary to its 3'-5'span).",
    "HEADER: **********************************************************************",
    "HEADER: POS STR L    Q SEQUENCE",
    "HEADER: **********************************************************************")
  )

}
################################################################################
