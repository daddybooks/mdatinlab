#' @title
#' check before the orginal data is converted because uncorectly copy may introdcing wrong names with
#'
#' @param x  original data from Waters software Empower.
#' @param dis distance between each sample in original data frame normally 6 or 12.
#' @param last The last sample's name check if they were identity
#'
#'
#'
#' @export

check_nsample<- function(x,dis= NULL,last = NULL)
{
  samplename<- (x[,1] %>% matrix(ncol = dis,byrow = T))[,1]
  ## 下机数据中的样本名是否和当前的一致？如果不一致，用户需要检查哪里出现问题，检查的依据是什么呢？

  ## 获得最后一个样本名
  last_samplename = samplename[(nrow(x)/dis )%>% round]
  if (last_samplename == last) {
    print("last sample name is indentity, you can go on")
  }  else {
    cat("处理得到的最后一个样本名是：",last_samplename,"请确认重哪一个样本名出现问题\t----------------------")
    samplename
  }}










