################################################################################
#//    Aleksandr B. Sahakyan (aleksahak [at] cantab.net), Cambridge 2016     \\#
################################################################################
INFOline <- function(OUT=OUT, msg="BLABLA", initial=FALSE){

  if(initial==TRUE){
    OUT <- msg
    #write(OUT, file="process_info.txt")
  } else {
    OUT <- c(OUT,msg)
    #write(OUT[length(OUT)], file="process_info.txt", append=TRUE)
  }
  print(OUT[length(OUT)], quote=FALSE)
  return(OUT)
  
}
################################################################################
