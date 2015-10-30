#' @useDynLib MACode
#' @export
doStuff <- function() {
    .C('dostuff', arg=as.integer(0), PACKAGE='MACode')$arg
}
