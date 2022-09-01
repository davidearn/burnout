##' Transparent colour
##' @description Given a colour, change it to a transparent version
##'     with a given \code{alpha} level
##' @param col as in \code{\link{col2rgb}}: vector of any of the three
##'     kinds of R color specifications, i.e., either a colour name
##'     (as listed by \code{colors()}), a hexadecimal string of the
##'     form \code{"##rrggbb"} or \code{"##rrggbbaa"} (see
##'     \code{\link{rgb}}), or a positive integer \code{i} meaning
##'     \code{\link{palette}}\code{()[i]}.
##' @param alpha opacity level on scale from \code{0} to \code{1}
##'     (i.e., \code{1} means completely opaque)
##' @seealso \code{\link{col2rgb}}, \code{\link{rgb}}
##' @importFrom grDevices col2rgb
##' @importFrom grDevices rgb
##' @export
transparent_colour <- function(col,alpha=150/255) {
    alpha <- round(alpha * 255)
    v <- col2rgb(col)[,1] ## color as rgb vector
    tcol <- rgb(v["red"],v["green"],v["blue"],alpha=alpha,maxColorValue = 255)
    return(tcol)
}