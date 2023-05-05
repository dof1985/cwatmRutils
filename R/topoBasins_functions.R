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

#' Delineate a watershed based on ldd and an outlet point
#'
#' @description
#' This function create basin mask map based on user inputs of ldd and XY location of the outlet.
#'
#' @param ldd Input ldd map
#' @param lat,lon latitude and longitude of the outlet in the ldd reference system
#' @param crop if set to 'TRUE' the resulting `RasterLayer` is cropped to the extent of the mask
#' @param loud `logical` if TRUE, the search list (including cell ids) will be printed every iteration
#' @return a `RasterLayer` of the basin mask
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
maskmapFromXY <- function(ldd, lat, lon, crop = FALSE, loud = TRUE) {

  getUpstream <- function(x, cell) {
    ngbrs <- raster::adjacent(x, cell, directions = 8, pairs = FALSE)
    # order 7 4 1 9 6 3 8 2

    if(length(ngbrs) < 8) {
      ngbrs <- unlist(lapply(c(1, 2, 3, 7, 8, 9, 4, 6), function(i) {
        v <- rep(NA, 9)
        v[i] <- 1

        mx <- matrix(v, nrow = 3, ncol = 3)
        mx[2, 2] <- 0
        tmp <- raster::adjacent(x, cell, directions = mx, pairs = FALSE)
        if(length(tmp) == 0) {
          return(NA)
        } else {
          return(tmp)
        }

      }))
    }
    tmpdf <- data.frame(nvals = x[ngbrs],
                        rvals = c(3, 6, 9, 1, 4, 7, 2, 8))
    return(ngbrs[tmpdf$nvals == tmpdf$rvals])
  }

  result <- srchList <- raster::cellFromXY(ldd, xy = c(lon, lat))


  while(length(srchList) > 0) {
    if(loud) print(srchList)
    newCells <- getUpstream(x = ldd, cell = srchList[1])
    newCells <- newCells[!is.na(newCells)]
    srchList <- c(srchList, newCells)
    srchList <- srchList[-1]

    result <- c(result, newCells)

  }
  rout <- ldd
  rout[] <- 0
  rout[result] <- 1

  if(crop) {
    maskExtent <- raster::extentFromCells(ldd, result)
    rout <- raster::crop(rout, maskExtent)
  }

  return(rout)
}

