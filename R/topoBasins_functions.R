# Topographical and basin associated operations



#' Delineate a watershed based on ldd and an outlet point
#'
#' @description
#' This function create basin mask map based on user inputs of ldd and XY location of the outlet.
#'
#' @param ldd Input flow direction map (Only PCRaster ldd is accepted)
#' @param lat,lon latitude and longitude of the outlet in the ldd reference system
#' @param loud `logical` if TRUE, the search list (including cell ids) will be printed every iteration
#' @return a named list
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
maskmapFromXY <- function(ldd, lat, lon, loud = TRUE) {

  getUpstream <- function(x, cell) {
    ngbrs <- raster::adjacent(x, cell, directions = 8, pairs = FALSE)
    # order 7 4 1 9 6 3 8 2
    tmpdf <- data.frame(nvals = x[ngbrs],
                        rvals = c(3, 6, 9, 1, 4, 7, 2, 8))
    return(ngbrs[tmpdf$nvals == tmpdf$rvals])
  }

  result <- srchList <- raster::cellFromXY(ldd, xy = c(lon, lat))


  while(length(srchList) > 0) {
    if(loud) print(srchList)
    newCells <- getUpstream(x = ldd, cell = srchList[1])
    srchList <- c(srchList, newCells)
    srchList <- srchList[-1]

    result <- c(result, newCells)

  }
  rout <- ldd
  rout[] <- 0
  rout[result] <- 1

  return(rout)
}

