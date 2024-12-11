# Topographical and basin associated operations


#' Identify inlets into desalinated spatial units
#'
#' @description
#' Utilize the ldd and the stream-network to identify the inlets into a deliniated spatial units, e.g., river basins or administrative units.
#'
#' @param r a spatial unit `RasterLayer` object
#' @param network a stream-network `RasterLayer` object
#' @param ldd a ldd `RasterLayer` object
#' @return a named `list` containing a inlet `RasterLayer` and `data.frame` objects
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
getInlets <- function(r, network,ldd) {
  tmp <- raster::as.matrix(r)
  ldd_ <- raster::crop(ldd, r)
  network <- raster::crop(network, r)

  ## get dimensions
  dim_ <- dim(tmp)

  ## shifted matrices
  r_NA <- matrix(NA, nrow = dim_[1], ncol = 1)
  c_NA <- matrix(NA, nrow = 1, ncol = dim_[2])

  tmp_r <- cbind(r_NA, tmp[ , -dim_[2]])
  tmp_l <- cbind(tmp[ , -1], r_NA)

  tmp_d <- rbind(c_NA, tmp[-dim_[1], ])
  tmp_u <- rbind(tmp[-1, ], c_NA)

  con1 <- (tmp - tmp_r) !=0
  con2 <- (tmp - tmp_l) !=0
  con3 <- (tmp - tmp_d) !=0
  con4 <- (tmp - tmp_u) !=0

  con1[is.na(con1)] <- 0
  con2[is.na(con2)] <- 0
  con3[is.na(con3)] <- 0
  con4[is.na(con4)] <- 0

  ## ldd overlay - from cells

  r1 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 1) *(con3 * con2) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  r2 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 2) *(con3) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  r3 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 3) *(con3 * con1) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  r4 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 4) *(con2) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)

  r6 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 6) *(con1) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  r7 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 7) *(con4 * con2) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  r8 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 8) *(con4) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  r9 <- raster::raster(raster::as.matrix(network) * (raster::as.matrix(ldd_ == 9) *(con4 * con1) > 0) * 1,
                       xmn = r@extent@xmin, xmx = r@extent@xmax,
                       ymn = r@extent@ymin, ymx = r@extent@ymax,
                       crs = r@crs)
  inlets <- reg
  inlets[] <- 0
  ## get downstream
  inlts_tbl <- do.call("rbind", lapply(c(1:4, 6:9), function(dir) {
    tmpr <- get(paste0("r", dir))
    tmp <- raster::as.data.frame(tmpr, xy = TRUE)
    tmp <- tmp[!is.na(tmp$layer) & tmp$layer == 1, ]
    tmp$cell <- raster::cellFromXY(tmpr, tmp[c("x", "y")])
    tmp$toreg <- unlist(lapply(tmp$cell, function(ic) {
      tocell <- raster::adjacent(tmpr, cell = ic, dir = 8)[ which(c(7, 4, 1, 9, 6, 3, 8, 2) %in% dir), 2]
      reg[tocell]
    }))
    return(tmp)
  }))

  inlets[inlts_tbl$cell] <- inlts_tbl$toreg


  return(list("inletGrid" = inlets,
              "inletTbl" = inlts_tbl))
  # transform to right left

}



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

