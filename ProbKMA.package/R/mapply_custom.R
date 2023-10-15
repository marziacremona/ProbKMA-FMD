
mapply_custom <- function(cl,FUN,...,MoreArgs=NULL,SIMPLIFY=TRUE,USE.NAMES=TRUE){
  if(is.null(cl)){
    mapply(FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }else{
    library::clusterMap(cl,FUN,...,MoreArgs=MoreArgs,SIMPLIFY=SIMPLIFY,USE.NAMES=USE.NAMES)
  }
}