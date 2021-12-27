#' @title process your specturm with different kinds of prep-processing methods
#'
#' @description We provide five kinds of pre-processing methods to help compare different results among prep.msc() prep.norm() prep.ref2km() prep.savgol() prep.snv()
#'
#'
#' @param x
#' specturm data, can be produced by read_nir()
#'
#' @param y
#' reference value.
#' @param cv
#'  cross valiation number
#' @param ncomp
#' 	maximum number of components to calculate.You need to choose how many compents for each pre-processing methods.
#'
#' @usage
#' preprocess_compare(x,y,cv = NULL,ncomp = NULL)
#'
#'
#'
#'
#' @export
preprocess_compare <- function(x = NULL,y =NULL ,cv = NULL,ncomp=NULL){
  a <- list(prep.msc(x),
    prep.norm(x),
    prep.savgol(x),
    prep.snv(x ),
    prep.ref2km(x)
    )
  b <-ncomp
  names(a) <- c("prep.msc","prep.norm","prep.savgol","prep.snv"," prep.ref2km")
  ml_pingguosuan_cat <- list(prep.msc= NULL,prep.norm = NULL,prep.savgol =NULL,prep.snv = NULL, prep.ref2km = NULL)


  for (i in 1:5){
    cat(paste(names(a)[i]," 's summary",sep = ""))
    ml_pingguosuan_cat[[i]]<-pls(a[[i]],
      as.vector(y),
      ncomp = b[i] ,
      cv = cv)
    cat(summary(ml_pingguosuan_cat[[i]]))
  }
  return(lapply(X = ml_pingguosuan_cat ,FUN = function(X){
    rbind(X$calres$r2[X$ncomp.selected],
      X$calres$rmse[X$ncomp.selected],
      X$cvres$r2[X$ncomp.selected],
      X$cvres$rmse[X$ncomp.selected],
      X$ncomp.selected)
  })  %>% data.frame(row.names = c("Rcal","RSMEC","Rcv","RMSECV","selected component")))
}
