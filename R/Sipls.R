#' @title apply Synergy interval PLSPLS for specturm data
#'
#' @description we seperated specturm into differents parts, and just randonly selected differents numbers of region to PLS modeling. And most improtant, we make sure every kinds of combination is used to modeling and record the results. NOTE：this function can be time consuming.
#'
#' @importFrom stringi stri_c
#' @param x
#' specturm data
#' @param y
#' reference value
#' @param nint
#'  numbers of intervals that has been devided. It should be divisible by numbers of ncol(x)
#' @param input
#' How many intervals would be apply to the models. Round and it's a string of numbers
#'
#' @usage
#' sipls(spectrum,
#' y,
#' cv = NULL,
#' nint = NULL,
#' input = NULL
#' )
#' @importFrom magrittr %>%
#'@export
Sipls <- function(spectrum = NULL,y = NULL,nint = NULL,input = NULL,cv){
  x = spectrum
  n = nint
  m = input
  a <- data.frame()
  for (i in n) {

    for (j in m) {
      cat(paste("划分区间数量：",i,"\n选取输入的数据区间数量:",j,"\n",sep = ""))
      ma_subset <- combn(i,j)
      ## 对数据集进行分类标签

      x = spectrum[,1:8000]
      x <-  t(x) %>% data.frame()
      x$tag = NA
      x$tag <- rep(1:i,each = 8000/i)

      for( k in 1:ncol(ma_subset)){
        cat(paste("目前输入的区间数",stringi::stri_c(ma_subset[,k] ,collapse = ","),"\n",sep = ""))

        spe_subset<- x[x[["tag"]] %in% ma_subset[,k],]
        ml<- pls(x = prep.snv(t(spe_subset[,-ncol(x)])),y = y,cv=cv,ncomp = 20)

        b <- data.frame(ml$res$cal$r2[[ml$res$cal$ncomp.selected]],
          ml$res$cal$rmse[[ml$res$cal$ncomp.selected]],
          ml$res$cv$r2[[ml$res$cv$ncomp.selected]],
          ml$res$cv$rmse[[ml$res$cv$ncomp.selected]],
          ml$res$cv$ncomp.selected,
          i,
          stringi::stri_c(ma_subset[,k] ,collapse = ","))


        a<- rbind(a,b)

      }
    }
  }

  colnames(a) <- c("Rcal","RMSE","Rcv","RMSECV","ncomp","interval","selected interval")
  return(a[order(a[,3],decreasing = T),])
}
