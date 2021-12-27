#' @title read the hplc data from Waters machine by copying with hands.
#'
#' @param file  file path.
#' @examples
#' read_hplc("./file/path/to/your/data")
#'
#'@export

read_hplc<- function (file, header = F, sep = ",", quote = "\"", dec = ".",
  fill = TRUE, comment.char = "", ...)
  read.table(file = file, header = header, sep = sep, quote = quote,
    dec = dec, fill = fill, comment.char = comment.char, ...)
