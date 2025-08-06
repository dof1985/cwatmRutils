# Topographical and basin associated operations

#' Get upstream cell indices
#'
#' @description
#' Identify the upstream cell flowing into each grid cell, according to a ldd.
#'
#' @param x a ldd `RasterLayer` object
#' @return a `data.frame` associating the id of each grid cell with its one or more direct upstream cells
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
getUpstream <- function(x) {
  names(x) <- "layer"
  downdir <- data.frame("val" = c(1:4, 6:9),
                        "row" = c(1, 1, 1, 0, 0, -1, -1, -1),
                        "col" = c(-1, 0, 1, - 1, 1, -1, 0, 1))

  allcells <- 1:raster::ncell(x)
  upc <- do.call("rbind", (lapply(downdir$val, function(dir) {
    upldd <- upcells <- x
    upldd[] <- upcells[] <- NA
    cells_ <- allcells
    upcells[cells_] <- raster::cellFromRowCol(x,
                                              row = raster::rowFromCell(x, cells_) - downdir[downdir$val == dir, "row"],
                                              col = raster::colFromCell(x, cells_) - downdir[downdir$val == dir, "col"])
    df <- data.frame("cells" = cells_,
                     "upcells" = raster::as.data.frame(upcells)$layer)
    df <- df[x[df$upcells] == dir, ]
    # upldd[df$cells] <- x[df$upcells]
    # (upldd == 1)[859]

    return(df)
  })))
  return(upc)
  #upc maps upstream table ids for each cell. If you also add a map of discharge - this can be easily sum, so each cell would have all upstream inflows
}


#' Identify inlets into desalinated spatial units
#'
#' @description
#' Utilize the ldd and the stream-network to identify the inlets into a delineated spatial units, e.g., river basins or administrative units.
#'
#' @param r a spatial unit `RasterLayer` object
#' @param network a stream-network `RasterLayer` object
#' @param ldd a ldd `RasterLayer` object
#' @param plot a logical variable. If `TRUE`, a plot of the inlet points, spatial units, and river network is plotted `RasterLayer` object
#' @return a `data.frame` object with the inflow point location, cell ids, and associated from- and to-zones
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
getInlets <- function (r, network, ldd, plot_ = FALSE) {
  cell_streams <- raster::Which(network == 1, cells = TRUE)
  ups_ <- na.omit(cwatmRutils::getUpstream(ldd))
  ups_ <- ups_[ups_$cells %in% cell_streams, ]
  ups_ <- ups_[r[ups_$cells] != r[ups_$upcells], ]
  ups_$fromreg <- r[ups_$upcells]
  ups_$toreg <- r[ups_$cells]

  ups_$to_x <- raster::xFromCell(r, ups_$cells)
  ups_$to_y <- raster::yFromCell(r, ups_$cells)

  if(plot_) {
    raster::plot(r)
    raster::plot(cwatmRutils::val_to_na(network > 0), col = "black", add =  TRUE, legend = FALSE)
    raster::plot(sp::SpatialPoints(ups_[c("to_x", "to_y")]), col = "red", add = TRUE, legend = FALSE)
  }
  return(ups_)
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

  # v2 - optimized - run vectors of cells simulatnously.
  # GET SHIFT MATRICES
  tmp <- raster::as.matrix(ldd)
  dim_ <- dim(tmp)
  r_NA <- matrix(NA, nrow = dim_[1], ncol = 1)
  c_NA <- matrix(NA, nrow = 1, ncol = dim_[2])
  tmp_r <- cbind(r_NA, tmp[, -dim_[2]])
  tmp_l <- cbind(tmp[, -1], r_NA)
  tmp_d <- rbind(c_NA, tmp[-dim_[1], ])
  tmp_u <- rbind(tmp[-1, ], c_NA)
  tmp_dr <- cbind(r_NA, tmp_d[, -dim(tmp_d)[2]])
  tmp_dl <- cbind(tmp_d[, -1], r_NA)
  tmp_ur <- cbind(r_NA, tmp_u[, -dim(tmp_u)[2]])
  tmp_ul <- cbind(tmp_u[, -1], r_NA)

  # create shift rasters

  sr1 <- raster::raster(tmp_dl, xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr2 <- raster::raster(tmp_d,  xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr3 <- raster::raster(tmp_dr,  xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr4 <- raster::raster(tmp_l, xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr6 <- raster::raster(tmp_r,  xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr7 <- raster::raster(tmp_ul,  xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr8 <- raster::raster(tmp_u, xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)
  sr9 <- raster::raster(tmp_ur, xmn = ldd@extent@xmin, xmx = ldd@extent@xmax,
                        ymn = ldd@extent@ymin, ymx = ldd@extent@ymax, crs = ldd@crs)

  conds_ <- setNames(lapply(c(1:4, 6:9), function(i) {
    get(paste0("sr", i)) == i
  }), nm = paste0("cond", c(1:4, 6:9)))


  #cell = srchList[1]
  getUpstream <- function(srchCells,  x, conds) {

    newcells <- integer(0)

    newcells <- c(newcells, conds$cond1[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells) - 1, col = raster::colFromCell(ldd, srchCells) + 1))
    newcells <- c(newcells, conds$cond2[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells) - 1, col = raster::colFromCell(ldd, srchCells)))
    newcells <- c(newcells, conds$cond3[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells) - 1, col = raster::colFromCell(ldd, srchCells) - 1))
    newcells <- c(newcells, conds$cond4[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells), col = raster::colFromCell(ldd, srchCells) + 1))
    newcells <- c(newcells, conds$cond6[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells), col = raster::colFromCell(ldd, srchCells) - 1))
    newcells <- c(newcells, conds$cond7[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells) + 1, col = raster::colFromCell(ldd, srchCells) + 1))
    newcells <- c(newcells, conds$cond8[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells) + 1, col = raster::colFromCell(ldd, srchCells)))
    newcells <- c(newcells, conds$cond9[srchCells] * raster::cellFromRowCol(ldd, row = raster::rowFromCell(ldd, srchCells) + 1, col = raster::colFromCell(ldd, srchCells) - 1))

    newcells <- newcells[newcells != 0]

    return(newcells)
  }


  result <- srchList <- raster::cellFromXY(ldd, xy = c(lon, lat))


  while (length(srchList) > 0) {
    if (loud)
      print(srchList)
    newCells <- getUpstream(srchCells = srchList, x = ldd, conds = conds_)
    newCells <- newCells[!is.na(newCells)]
    newCells <- newCells[!newCells %in% result]
    srchList <- newCells
    #srchList <- srchList[-1]
    result <- c(result, newCells)
  }
  rout <- ldd
  rout[] <- 0
  rout[result] <- 1
  if (crop) {
    maskExtent <- raster::extentFromCells(ldd, result)
    rout <- raster::crop(rout, maskExtent)
  }
  return(rout)

}
