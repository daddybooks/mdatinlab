#' @title Check if there are factitious error causing duplication.
#'
#' @description
#' apply by using \code{duplicated()} funtion and it will check on every row
#'
#' @param x
#' a data frame after converted
#'
#'
#' @usage
#' check_ifdup(x)
#' @importFrom magrittr %>%
#' @export
check_ifdup<- function(x){
  if(any(x %>% duplicated())){
         print("there is a duplicated copy pls check")
     }
  else {


             print("OK!")
}}
