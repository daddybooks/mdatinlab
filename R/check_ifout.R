#' @title check if there are outliers exist in our dataframe
#'
#' @description by using the funtion outlier.IQR, we mainly extract colume of time and point out which value is outlier.
#'
#'
#' @param x
#' A dataframe after converted.
#'
#' @usage
#' read_ifout(x)
#' @export
check_ifout<- function(x){
  outlier.IQR <- function(x, multiple = 1.5, replace = FALSE, revalue = NA) {
    q <- quantile(x, na.rm = TRUE) #四分位间距3倍间距以外的认为是离群值
    IQR <- q[4] - q[2]
    x1 <- which(x < q[2] - multiple * IQR | x > q[4] + multiple * IQR)
    x2 <- x[x1]
    if (length(x2) > 0) outlier <- data.frame(location = x1, value = x2)
    else outlier <- data.frame(location = 0, value = 0)
    if (replace == TRUE) {
      x[x1] <- revalue
    }
    return(list(new.value = x, outlier = outlier))
  }
  ### 通过四分位数来检测是否存在离群值
  for( i in 1:(grep("time",colnames(x)) %>% length())) {
    (x[,grep("time",colnames(x))][,i] %>%
        as.numeric() %>%
        outlier.IQR())$outlier %>% print}
}

