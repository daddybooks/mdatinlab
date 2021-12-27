#' @title a function used to read data from NIR machine data
#'
#' @description This function can automacially read all the csvfiles
#'
#'
#' @param path
#' the csvfiles' file path
#'
#' @usage
#' raed_nir(path)
#' @return
#'  returns a data frame with all spectra information in the ".csv" files.
#' @export
read_nir <- function(path){
  ## 首先进入csv文件所在的文件夹，PATH，其次对读取文件进行处理
  csvfiles<- list.files(path = path)[list.files(path = path) %>% grep(pattern = ".*.csv")]
  z <- list()
  for (i in csvfiles) {
    z[[i]] <- read.csv(paste(path,i,sep = ""),fill = T,stringsAsFactors = T,header = F)[,2]
  }
  dat<- data.frame(z)
  ## 转置
  dat1 <- t(dat)
  ## 转置之后第一列第二列全部是单位A,以及其他信息 可以去除
  dat1<- dat1[,-(1:2)]
  ## 完成一个frame 之后，要数字格式
  dat1 <- matrix(data =(as.vector(dat1) %>% as.numeric()),nrow =nrow(dat1),ncol = ncol(dat1) )
  ## 重新名命样本名
  rownames(dat1) <- ((csvfiles) %>% gsub(pattern = ".csv",replacement = ""))
  return(dat1)
}
