# Topographical and basin associated operations

#' Reclassifies  flow-direction grid to a ldd
#'
#' @description
#' This function converts the flow-direction raster (by ArcGIS) to a ldd (by PCRaster).
#'
#' @param flowdir Input flow-direction map
#' @param noflow a value defined for sinks (default = 0)
#' @return a ldd `RasterLayer` object
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
flowdir2ldd <- function(flowdir, noflow = 0) {
  ctbl <- setNames(c(1, 2, 4, 8, 16, 32, 64, 128, noflow),
                   nm = c(6, 3, 2, 1, 4, 7, 8, 9, 5))

  lddout <- flowdir
  lddout[] <- 0

  for(i in seq_len(length(ctbl))) {
    lddout[ctbl[i] == flowdir] <- as.numeric(names(ctbl)[i])
  }
  return(lddout)
}

#' Point & click to zoom in on a RasterLayer object
#'
#' @description
#' This function zoom in to a RasterLayer using R `base` graphics and the `raster` libraries.
#' Calling the function plots the raster layer. Zoom in is triggered by a user point and click.
#' The scope of the output raster is defined by the `size` argument, resulting in a window of area `(2*size + 1)^2`.
#'
#' @param r Input flow direction map (Only PCRaster ldd is accepted)
#' @param size pre-defined distance (grid-cell) for defining a zoom-in window (default = 15)
#' @return a `RasterLayer` object
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
zoomClip <- function(r, size = 15) {
  plot(r)
  cid <- raster::click(r, n = 1, cell = TRUE, show = FALSE)
  mx <- matrix(1, nrow = 1 + 2 * size, ncol = 1 + 2 * size)
  mx[size, size] <- 0
  ngbrs <- raster::adjacent(r, cell = cid$cell, directions = mx,
                            include = TRUE, pairs = FALSE)

  return(raster::crop(r, raster::extentFromCells(r, cells = ngbrs)))
}

#' Delineate a watershed based on ldd and an outlet point
#'
#' @description
#' This function create basin mask map based on user inputs of ldd and XY location of the outlet.
#'
#' @param ldd Input ldd map
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

