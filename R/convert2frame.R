#'@title convert HPLC data into a dataframe
#'
#' @description There are superfluous information among the HPLC data from machine. This function can extract the main information we need to next analyis.
#'
#'
#'
#' @name convert2frame
#' @usage convert2frame(data,
#' dis,
#' npeak)
#' @format
#'
#' @export
convert2frame<- function(x ,dis = NULL,npeak = NULL)
{

  ## 如何检查:用nrow/dis = nsamples
  nsamples <-  (nrow(x)/6) %>% round()
  idx <-seq(1,by = dis ,length.out = nsamples)
  peaktime <- matrix(data = x[,3],ncol =dis,byrow = TRUE)
  peak<- matrix(data = x[,4],ncol = dis,byrow = TRUE)
  samplenames<- x[,1][idx]
  dat2<- data.frame(samplenames,peak,peaktime)
  dat2<- dat2[,-which(dat2[1,] == ""|dat2[1,] == "面积"|dat2[1,] == "保留时间")]
  colnames(dat2) <- c("samplenames",paste("peak",1:npeak,sep = ""),paste("time",1:npeak,sep = ""))
  dat2
}

