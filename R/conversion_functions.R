# conversions



#' Gets Information and metadata from netcdf
#'
#' @description
#' This function extracts dimensions and attributes from 'NetCDF' file.
#'
#' @param pth a path to a NetCDF file with the extension '.nc' or '.nc4'
#' @param dim `logical` if TRUE, dimensions the files are returned
#' @param attrs `logical` if TRUE, variable and global attributes are returned
#' @return a named list
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
ncdfInfo <- function(pth, dim = TRUE, attrs = FALSE) {
  tmp <- ncdf4::nc_open(pth)
  info <- list()
  vars <- names(tmp$var)
  if(dim) {
    dims <- names(tmp$dim)
    info <- c(info, "vars"= list(vars), "dims" = list(dims))
    #if(getVals) {
    #  #getVals = TRUE
    #  for(dim1 in dims) {
    #    tmpList <- tmp$dim[[dim1]]$vals
    #  }
    }

  if(attrs) {
    varattr <- do.call("rbind", lapply(vars, function(v) {
      outdf <- as.data.frame(t(as.data.frame(unlist(ncdf4::ncatt_get(tmp, varid = v)))))
      row.names(outdf) <- v
      return(outdf)
    }))

    globattr <- as.data.frame(t(as.data.frame(unlist(ncdf4::ncatt_get(tmp, varid = 0)))))
    row.names(globattr) <- NULL

    info <- c(info, "varAttributes" = varattr, "globalAttributes" = globattr)
  }
  ncdf4::nc_close(tmp)
  return(info)
}


#' Imports and converts NetCDF files to 'RasterLayer','RasterStack', and 'data.frame'
#'
#' @description
#' This function imports netCDF files to R 'RasterLayer','RasterStack', and 'data.frame'.
#'
#' @details
#' The function can imports 2 or 3-dimensional single/multi-variable NetCDF files into R. It provides users
#' with temporal and spatial subset and summary options.
#'
#' The `spatial` argument  accepts a point coordinate data.frame (x, y) or a `RasterLayer` mask as input. If points are provided, the function
#' returns a data.frame with the value for each point, time, and variable. A `RasterLayer` input return the masked map by default. The user
#' can use a weight `RasterLayer` as a spatial mask, so the result is the multiplication of the values by the weights.
#'
#' The `time` argument accepts one or two (from, to) points in time as a `numeric`/`integer` index or as `Date`. In the former case, the i to j th
#' time points are being extracted from the file. When `Date` is provided it is first being converted to days since the defined `origin`, which is
#' later searched within the values of the file's time dimension. `origin` should fit units defined in the file's time dimension.
#'
#' Summarizing the result can be applied spatially with the `fun` and `...` argument, or temporally with the `temporal_fun` argument.
#' Temporal summary applies pre-defined statistical transformation (e.g., sum, mean, sd, and cv, coefficient of variance) to every grid-cell on the x-y plane across time points,
#' thus it converts 3 dimensional array to 2 dimensional.
#' Spatial summary applies user defined statistical transformation (e.g., sum, mean, sd) to every time point, resulting in a `data.frame` with the
#' var, time, and summarized value.
#' The user can apply both summaries at the same time.
#'
#' @param pth a path to a NetCDF file with the extension `.nc` or `.nc4`
#' @param flip `character` "x", "y" or NULL. If not set to `NULL`, the resulting array is being flipped to the defined direction.
#' @param transpose `logical` If `TRUE`, the resulting array is being transposed (defaults to `TRUE`).
#' @param time If not set to `NULL`, defines temporal subset (see Details).
#' @param origin temporal origin of the input NetCDF, defaults to `1901-01-01` (Optional; see Details).
#' @param spatial If not set to `NULL`, defines spatial subset (see Details).
#' @param varName If not set to `NULL`, specifies the variable to load.
#' @param fun function for spatial summarize.
#' @param temporal_fun `character` One of the following: c("sum", "mean", "sd", "cv")
#' @param crs proj4 string input to `raster::crs()` used to construct the output `RasterLayer`
#' @param ... additional arguments to function provided to the `fun` argument.


#' @param attrs `logical` if TRUE, variable and global attributes are returned
#' @return a `data.frame`, a `RasterLayer`, a `RasterStack`, or a named `list`
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export

ncdf2raster <- function(pth, flip = NULL, transpose = FALSE, time = NULL, origin = "1901-01-01", spatial = NULL,
                        varName = NULL, fun = NULL, temporal_fun = NULL, crs = "+init=EPSG:4326", ...) {

  ## input validation
  if(!is.null(time)) {
    stopifnot("'time' argument should be of class 'Date', 'integer' or 'numeric'" = any(c("integer", "Date", "numeric") %in% class(time)))

    stopifnot("'time' argument should be of length 1 or 2" = length(time) == 1 || length(time) ==2)

    if(length(time) == 2) {
      stopifnot("The first member of the 'time' argument should be smaller than its second member" = time[2] > time[1])
    }



  }

  if(!is.null(spatial)) {
    stopifnot("'spatial' argument should be of class 'RasterLayer' or 'data.frame'" = any(c("RasterLayer", "data.frame") %in% class(spatial)))

    if(class(spatial) == "data.frame") {
      stopifnot("Points coordinates columun should recieve one of the following names: c('x', 'y'), c('X', 'Y'), c('lon', 'lat')" = all(c(any(c("lat", "Y", "y") %in% names(spatial)), any(c("lon", "X", "x") %in% names(spatial)))))
    }


  }

  if(!is.null(varName)) {
    stopifnot("'varName' should be of class 'character'" = class(varName) %in% "character")
  }


  if(!is.null(temporal_fun)) {
    stopifnot("'temporal_fun' should be of class 'character'" = class(temporal_fun) %in% "character")
    stopifnot("'temporal_fun' can recieve one of the following: 'sum', 'mean', 'sd', 'cv'" = temporal_fun %in% c("sum", "mean", "sd", "cv"))
  }

  ## functions
  getAxis <- function(array, idx, axis) {
    ndim <- length(dim(array))
    idx_list <- lapply(seq_len(ndim), function(d) {
      if(d == axis) {
        return(idx)
      } else {
        return(seq_len(dim(array)[d]))
      }
    })
    return(do.call("[", c(list(array), idx_list)))
  }


  # open file
  tmp <- ncdf4::nc_open(pth)

  dim_nmes <- c(grep("[Ll]atitude|[Ll]at", names(tmp$dim), value = TRUE),
                grep("[Ll]ongitue|[Ll]on", names(tmp$dim), value = TRUE))
  if(length(dim_nmes) == 0) {
    dim_nmes <- c(grep("[Yy]", names(tmp$dim), value = TRUE),
                  grep("[Xx]", names(tmp$dim), value = TRUE))
  }

  # get dim x, dim y
  y <- tmp$dim[[dim_nmes[1]]]$vals
  x <- tmp$dim[[dim_nmes[2]]]$vals

  resx <- x[2] - x[1]
  resy <- abs(y[2] - y[1])

  # set temporal dim
  timeExists <- "time" %in% names(tmp$dim)
  tempnm <- NULL

  if(timeExists && is.null(time)) {
    temp <- tmp$dim$time$vals
    s_time <- 1
    e_time <- length(temp)
    tempnm <- temp
  }

  if(!is.null(time)) {
    errmsg = "'time' argument included, but input data has no time dimension"
    stopifnot(errmsg = timeExists)

    temp <- tmp$dim$time

    # time inputs asDate
    if(class(time) %in% "Date") {
      ncheck <- length(time)
      time <- which(temp$vals %in% as.numeric(time - as.Date(origin)))
      stopifnot('The temporal coverage could not be fully retrived' = (length(time) == ncheck))
    }

    s_time <- time[1]
    e_time <- s_time
    if(length(time) == 2) e_time <- time[2] - s_time + 1

    tempnm <- temp$vals[s_time:(s_time + e_time - 1)]

  }

  s_x <- 1
  c_x <- -1
  s_y <- 1
  c_y <- -1

  # set spatial mask
  spatExists <- !is.null(spatial)
  isPts <- FALSE
  isMask <- FALSE

  if(spatExists) {
    isPts <- class(spatial) %in% "data.frame"
    isMask <- class(spatial) %in% "RasterLayer"
  }

  # pts
  if(isPts) {
    x_idx <- na.omit(match(c("lon", "X", "x"), names(spatial)))
    y_idx <- na.omit(match(c("lat", "Y", "y"), names(spatial)))

    x_loc <- unlist(lapply(spatial[ , x_idx], function(x1) {
      which.min(abs(x1 - x))
    }))

    y_loc <- unlist(lapply(spatial[ , y_idx], function(y1) {
      which.min(abs(y1 - y))
    }))

    #mask2array <- as.matrix(spatial)
    mask2Extent <- c(min(x[x_loc]), max(x[x_loc]), min(y[y_loc]), max(y[y_loc]))

    s_x <- which.min(abs(mask2Extent[1] - x))
    e_x <- which.min(abs(mask2Extent[2] - x))

    s_y <- which.min(abs(mask2Extent[4] - y))
    e_y <- which.min(abs(mask2Extent[3]- y))

    c_x <- e_x - s_x + 1
    c_y <- e_y - s_y + 1
  }

  # msk
  if(isMask) {


    mask2Extent <- raster::extentFromCells(spatial, raster::Which(!is.na(spatial), cell = TRUE))

    mask_count_x <- (mask2Extent@xmax - mask2Extent@xmin) / resx
    mask_count_y <- (mask2Extent@ymax - mask2Extent@ymin) / resy

    s_x <- which(min(abs(mask2Extent@xmin - x)) == abs(mask2Extent@xmin - x))
    s_x <- s_x[length(s_x)]
    # correct if mask is the same size as input ncdf
    e_x <- s_x + min(mask_count_x, length(tmp$dim$lon$vals) - 1)

    if(e_x > length(x)) {
      e_x <- e_x - 1
      s_x <- s_x - 1
    }
    #(x[e_x] - x[s_x]) / (e_x - s_x)

    s_y <- which(min(abs(mask2Extent@ymax - y)) == abs(mask2Extent@ymax - y))
    s_y <- s_y[length(s_y)]
    # correct if mask is the same size as input ncdf
    e_y <- s_y + min(mask_count_y, length(tmp$dim$lat$vals) - 1)

    if(e_y > length(y)) {
      e_y <- e_y - 1
      s_y <- s_y - 1
    }

    #(y[e_y] - y[s_y]) / (e_y - s_y)

    c_x <- e_x - s_x + 1
    c_y <- e_y - s_y + 1
  }


  varid <- names(tmp$var)
  if(is.null(varName)) {
    if(length(varid) > 1) warning(sprintf("varName is NULL, attempting to extract %s variables: %s;",
                                          length(varid),
                                          paste0(varid, collapse = ", ")))
  } else {
    varid <- varName
  }

  from <- c(s_x, s_y)
  counts <- c(c_x, c_y)


  if(timeExists) {
    # Search position of time in dim names

    if(tmp$var[[varid[1]]]$dim[[1]]$name == "time") {
      from <- c(s_time, from)
      counts <- c(e_time, counts)
    } else {
      from <- c(from, s_time)
      counts <- c(counts, e_time)
    }

  }
  out_ds <- setNames(lapply(varid, function(varid) {

    arr <-  tryCatch({
      ncdf4::ncvar_get(tmp, varid = varid, start = from, count = counts)
      }, error =  function(e) {
        from <- c(s_x, s_y)
        counts <- c(c_x, c_y)
        timeExists <- FALSE
        ncdf4::ncvar_get(tmp, varid = varid, start = from, count = counts)
             })
    #arr <- ncdf4::ncvar_get(tmp, varid = varid, start = from, count = counts)
    arrDims <- dim(arr)
    time_arrDim <- NULL
    if(timeExists) {
      time_arrDim <- 3
      if(tmp$var[[varid[1]]]$dim[[1]]$name == "time") time_arrDim <- 1
      #time_arrDim <- which(names(tmp$dim) %in% "time")
      # if(is.null(time)) {
      #   time_arrDim <- length(arrDims)
      # } else if (length(time) > 1) {
      #   time_arrDim <- length(arrDims)
      # }
    }

        temporal_sum <- FALSE
    if(!is.null(temporal_fun) && !is.null(time_arrDim) && !isPts) { # ignore points
      n <- dim(arr)[time_arrDim]
      rast_tmp <- raster::stack(lapply(seq_len(n), function(i) {
        raster::raster(getAxis(array = arr, idx = i, axis = time_arrDim))
      }))
      # 'sum', 'mean', 'sd', 'cv'

      naMask <- is.na(rast_tmp[[1]])
      if(temporal_fun == "sum") rast_tmp <- sum(rast_tmp, na.rm = TRUE)
      if(temporal_fun == "mean") rast_tmp <- sum(rast_tmp, na.rm = TRUE) / n
      if(temporal_fun == "sd") {
        m <-  sum(rast_tmp, na.rm = TRUE) / n
        rast_tmp <- sum((rast_tmp - m) ^ 2, na.rm = TRUE) / n
      }
      if(temporal_fun == "cv") {
        m <-  sum(rast_tmp, na.rm = TRUE) / n
        rast_tmp <- sum((rast_tmp - m) ^ 2, na.rm = TRUE) / n
        rast_tmp <- rast_tmp / m
      }

      rast_tmp[naMask] <- NA

      arr <- as.array(matrix(getValues(rast_tmp), nrow = rast_tmp@nrows, ncol = rast_tmp@ncols, byrow = TRUE))



      tempnm <- NULL
      # arrDims <- dim(arr)
      # if(is.null(time)) {
      #   time_arrDim <- length(arrDims)
      # } else if (length(time) > 1) {
      #   time_arrDim <- length(arrDims)
      # }
      temporal_sum <- TRUE
    }

    if(isMask) {
      xmn = x[s_x] - 0.5 * resx
      xmx = x[e_x] + 0.5 * resx
      ymn = y[e_y] - 0.5 * resy
      ymx = y[s_y]  + 0.5 * resy

      # check dims

      tmprast <- raster::crop(raster::extend(spatial, raster::extent(xmn, xmx, ymn, ymx)),
                              raster::extent(xmn, xmx, ymn, ymx))
      mask2array <- matrix(raster::getValues(tmprast), byrow = TRUE, nrow = tmprast@nrows, ncol = tmprast@ncols)
      #mask2array <- as.matrix(raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx)))
      if(transpose) mask2array <- raster::t(mask2array)
      if(!is.null(time_arrDim) && !temporal_sum) mask2array <- replicate(n = dim(arr)[time_arrDim], expr = mask2array, simplify = "array")

      if(all(dim(arr) == dim(mask2array))) {
        arr <- mask2array * arr
      } else {
        xmn = x[s_x] - 0.5 * resx
        xmx = x[e_x]# + 0.5 * resx
        ymn = y[e_y]# - 0.5 * resy
        ymx = y[s_y]#  + 0.5 * resy

        tmprast <- raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx))
        mask2array <- matrix(raster::getValues(tmprast), byrow = TRUE, nrow = tmprast@nrows, ncol = tmprast@ncols)
        #mask2array <- as.matrix(raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx)))
        if(transpose) mask2array <- raster::t(mask2array)
        if(!is.null(time_arrDim) && !temporal_sum) mask2array <- replicate(n = dim(arr)[time_arrDim], expr = mask2array, simplify = "array")
        if(all(dim(arr) == dim(mask2array))) {
          arr <- mask2array * arr
        } else {
          xmn = x[s_x]# - 0.5 * resx
          xmx = x[e_x] #+ 0.5 * resx
          ymn = y[e_y] - 0.5 * resy
          ymx = y[s_y] # + 0.5 * resy

          tmprast <- raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx))
          mask2array <- matrix(raster::getValues(tmprast), byrow = TRUE, nrow = tmprast@nrows, ncol = tmprast@ncols)
          #mask2array <- as.matrix(raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx)))
          if(transpose) mask2array <- raster::t(mask2array)
          if(!is.null(time_arrDim) && !temporal_sum) mask2array <- replicate(n = dim(arr)[time_arrDim], expr = mask2array, simplify = "array")
          if(all(dim(arr) == dim(mask2array))) {
          arr <- mask2array * arr
          } else {
          xmn = x[s_x] #- 0.5 * resx
          xmx = x[e_x] #+ 0.5 * resx
          ymn = y[e_y]# - 0.5 * resy
          ymx = y[s_y]#  + 0.5 * resy

          tmprast <- raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx))
          mask2array <- matrix(raster::getValues(tmprast), byrow = TRUE, nrow = tmprast@nrows, ncol = tmprast@ncols)
          #mask2array <- as.matrix(raster::crop(spatial, raster::extent(xmn, xmx, ymn, ymx)))
          if(transpose) mask2array <- raster::t(mask2array)
          if(!is.null(time_arrDim) && !temporal_sum) mask2array <- replicate(n = dim(arr)[time_arrDim], expr = mask2array, simplify = "array")
          arr <- mask2array * arr
        }}}}


    iter <- 1

    if(!(isPts && nrow(spatial) == 1)) {
      if(!is.null(time_arrDim) && !temporal_sum) iter <- seq_len(dim(arr)[time_arrDim])


      outr <- setNames(lapply(iter, function(l) {
      if(!is.null(time_arrDim)) {
        arr2rast <- as.matrix(getAxis(array = arr, idx = l, axis = time_arrDim))
        #if(time_arrDim == 2) arr2rast <- raster::t(arr2rast)
      } else {
        arr2rast <- as.matrix(arr)
        #if(isPts) arr2rast <- raster::t(arr2rast)
      }

      if(transpose) arr2rast <- raster::t(arr2rast)

      if(isPts) {
        x_ext <- x_loc - min(x_loc) + 1
        y_ext <- y_loc - min(y_loc) + 1
        if(transpose) {
          y_ext <- x_loc - min(x_loc) + 1
          x_ext <- y_loc - min(y_loc) + 1
        }

        do.call("rbind", lapply(seq_len(length(x_ext)), function(i) {
          dfpts <- data.frame("x" = x[x_loc[i]],
                              "y" = y[y_loc[i]],
                              "var"= varid, stringsAsFactors = FALSE)
          if(!is.null(time_arrDim)) dfpts$time <- tempnm[l]
          dfpts$value <-  arr2rast[x_ext[i], y_ext[i]]
          return(dfpts)
        }))


      } else {

        # coords
        xmn = min(x) - 0.5 * resx
        xmx = max(x) + 0.5 * resx
        ymn = min(y) - 0.5 * resy
        ymx = max(y) + 0.5 * resy

        if(isMask) {
          xmn = x[s_x]# - 0.5 * resx
          xmx = x[e_x]# + 0.5 * resx
          ymn = y[e_y]# - 0.5 * resy
          ymx = y[s_y]#  + 0.5 * resy

        }


        if(!is.null(fun)) {
          dfout <- data.frame("var"= varid, stringsAsFactors = FALSE)
          if(!is.null(time_arrDim) && is.null(temporal_fun)) dfout$time <- tempnm[l]
          dfout$value <- do.call(fun, list(as.numeric(arr2rast), ...))

          return(dfout)
        }
        #print(dim(arr2rast))
        #print((xmx - xmn) / resx)
        #print((ymx - ymn) / resy)

        rast <- raster::raster(arr2rast, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = raster::crs(crs))
        if(isMask) rast@extent <- mask2Extent
        if(!is.null(flip)) rast <- raster::flip(rast, direction = flip)

        return(rast)
      }

    }), nm = tempnm)
    } else {
      outr <- data.frame("x" = x[x_loc],
                          "y" = y[y_loc],
                          "var"= varid,
                          "value" = as.numeric(arr),
                          "time" = tempnm,
                          stringsAsFactors = FALSE)
    }
    if(length(outr) == 1) outr <- outr[[1]]
    if(isPts && !is.null(time_arrDim) && nrow(spatial) > 1) outr <- do.call("rbind", outr)
    if(!is.null(fun) && is.null(temporal_fun) && !is.null(time_arrDim)) outr <- do.call("rbind", outr)
    return(outr)
  }), nm = varid)

  if((isPts | !is.null(fun)) && length(varid) > 1) out_ds <- do.call("rbind", out_ds)
  if(class(out_ds) %in% "data.frame") row.names(out_ds) <- NULL
  if(class(out_ds) %in% "list") {
    if(length(varid) > 1 && ((length(tempnm) == 1) | is.null(tempnm))) {
      out_ds <- raster::stack(out_ds)
    } else if(length(varid) > 1) {
      out_ds <- setNames(lapply(tempnm, function(timename) {
        tmp <- raster::stack(lapply(varid, function(varname) {
          out_ds[[varname]][[as.character(timename)]]
        }))
        names(tmp) <- varid
        return(tmp)
      }), nm = tempnm)
    } else if(length(varid) == 1) {
      out_ds <- out_ds[[1]]
    }
  }

  ncdf4::nc_close(tmp)
  return(out_ds)
}






#' Converts raster to netCDF
#'
#' @description
#' This function converts a `RasterLayer`, `RasterStack` or a `list` of these classes to a `netCDF` file, and exports it to a pre-defined output path.
#'
#'
#' @details
#' Inputs can represent a single map (`RasterLayer`), a stacked map (e.g., with multiple variables; `RasterStack`),
#' a time-series (`list` of maps) or a time-series of stacked maps.
#' For each variable, a `name` must be provided, so the length of the `name` argument must be equal to the number of stacked layers in the input.
#' In case a time-series is used (e.g., input raster is a `list`), one shall provide a `time` argument, either as a `numeric` vector or as a `character` vector.
#' In the latter option, the `character` input will be converted to `numeric` representing the number of days since a user defined `origin`. Both the `time`
#' and `origin` arguments should be of format `%Y-%m-%d`.
#' For spatio-temporal data, a list of `RasterLayer`/`RasterStack` objects is expected. Each list item stands for a point in time, and
#' each `RasterLayer` stands for the spatial dimensions and data of a variable. If multiple variable are required, a `RasterStack` can be used, so
#' each layer stands for a variable.
#' The user must provide the following attributes: `name`, `unit`, `prec`, and `missing_value`. Other optional attributes include `longname` (for each variable),
#' or the global variables `title`, `author`, `institute`, `source`, and `description`.
#'
#' @param rast_in an input `RasterLayer`, `RasterStack` or a `list` (see Details)
#' @param path_out `character` path for the output file with the extension `.nc` or `.nc4` (if `is_ncdf4 = TRUE`)
#' @param name `character` variable name(s)
#' @param unit `character` variable attributes describing the unit
#' @param is_ncdf4 `logical`, if `TRUE` the output is a netCDF 4 file
#' @param prec `character` variable type, either `integer` or `float`
#' @param missing_value missing value
#' @param time `character` representation of a `Date` with format `%Y-%m-%d` or  `numeric` value
#' @param origin `character` representation of an origin `Date`, used to convert `character` time argument to `numeric`.
#' @param axis_names `character` vector naming the axis names. Provides flexibility in cases the third axis is not time.
#' @param longname `character` variable long name(s)
#' @param title `character` a global attribute
#' @param author `character` a global attribute
#' @param institute `character` a global attribute
#' @param source `character` a global attribute
#' @param description `character` a global attribute
#' @return NULL
#'
#' @examples
#' \dontrun{
#' PUT EXAMPLE HERE
#' }
#' @export
raster2ncdf <- function(rast_in, path_out, name, unit, is_ncdf4 = FALSE, prec = "float", missing_value = 32000,
                        time = NULL, origin = "1901-01-01", time_unit = "days",
                        longname = NULL, title =NULL, author = NULL, institute = NULL,
                        source = NULL, description = NULL, axis_names = c("lon", "lat", "time")) {

  # flip and transpose are not used - currently  flip = NULL, transpose = FALSE,

  if(class(rast_in) %in% "RasterBrick") stop("RasterBrick intputs are not supported, use RasterStack")
  # check input rast_in
  if(class(rast_in) %in% "RasterLayer") {
    #print("Simple one map to nc")
    r_template <- rast_in[[1]]
    stopifnot("'name' argument shall be of the same length as the number of variable" = length(name) == raster::nlayers(rast_in))
    rast_list <- list(rast_in)
  }

  if(class(rast_in) %in% "RasterStack") {
    #print("Multi-variable map to nc")
    r_template <- rast_in[[1]]
    stopifnot("'name' argument shall be of the same length as the number of variable" = length(name) == raster::nlayers(rast_in))
    stopifnot("'name' should hold unique values" = length(unique(name)) == raster::nlayers(rast_in))
    rast_list <- rast_in
  }

  if(class(rast_in) %in% "list") {

    r_template <- rast_in[[1]][[1]]

    stopifnot("time-series input requires a numeric/character 'time' vector" = !is.null(time))
    stopifnot("time-series input requires a numeric/character 'time' vector" = class(time) %in% c("character", "numeric", "integer"))

    if(is.character(time)) {
      time <- as.numeric(as.Date(time) - as.Date(origin))
      #print(time)
    }

    if(class(rast_in[[1]]) %in% "RasterLayer") {
      #print("Time-series of simple map to nc")
      rast_list <- list(rast_in)
    }

    if(class(rast_in[[1]]) %in% "RasterStack") {
      #print("Time-series of multi-variable map to nc")
      stopifnot("'name' argument shall be of the same length as the number of variable" = length(name) == raster::nlayers(rast_in[[1]]))

      rast_list <- lapply(seq_along(name), function(i) {
        lapply(seq_len(length(rast_in)), function(t) {
          rast_in[[t]][[i]]
        })
      })
    }
  }

  # build x, y dimensions
  xvals <- raster::xFromCol(r_template, seq_len(r_template@ncols))
  yvals <- raster::yFromRow(r_template, seq_len(r_template@nrows))

  nx <- r_template@ncols
  ny <- r_template@nrows

  x_def <- c(axis_names[1], "degrees_east")
  y_def <- c(axis_names[2], "degrees_north")

  r_crs <- sf::st_as_text(sf::st_crs(r_template@crs))
  axis_str <- substr(r_crs,  regexec("AXIS\\[", r_crs)[[1]], nchar(r_crs))
  if(regexec("Latitude", axis_str)[[1]] == -1) {
    # set x, y_def to metric
    x_def <- c(axis_names[1], "meters")
    y_def <- c(axis_names[2], "meters")
  }

  x_lon <- ncdf4::ncdim_def(x_def[1], x_def[2], xvals)
  y_lat <- ncdf4::ncdim_def(y_def[1], y_def[2], yvals)

  dims <- list(x_lon, y_lat)

  nt <- 1
  # if time included add time dimension
  if(!is.null(time)) {
    nt <- length(time)
    t_time <- ncdf4::ncdim_def(axis_names[3], units = sprintf("%s since %s", time_unit, origin), calendar =  "standard", vals = 0, unlim = TRUE)
    dims <- list(x_lon, y_lat, t_time)
  }

  if(is.null(longname)) longname <- name

  var_tmp <- lapply(seq_along(name), function(ivar) {
    ncdf4::ncvar_def(name = name[ivar], units = unit, dim = dims, prec = prec, missval = missing_value)

  })

  ncnew <- ncdf4::nc_create(path_out, var_tmp, force_v4=is_ncdf4)
  for(ivar in seq_along(name)) {
    rast_list_var <- rast_list[[ivar]]

    # put longname
    ncdf4::ncatt_put(ncnew, var_tmp[[ivar]], "long_name", longname[ivar])

    for(itime in seq_len(nt)) {

      r_towrite <- raster::t(rast_list_var[[itime]])
      arr_towrite <- matrix(getValues(r_towrite), nrow = r_towrite@nrows, ncol = r_towrite@ncols, byrow = TRUE)

      if(!is.null(time)) {
        ncdf4::ncvar_put(nc = ncnew, varid = var_tmp[[ivar]], vals = arr_towrite, start = c(1, 1, itime), count = c(nx, ny, 1))
        ncdf4::ncvar_put(nc = ncnew, varid = t_time, vals = time[itime], start = itime, count = 1)
      } else {
        ncdf4::ncvar_put(nc = ncnew, varid = var_tmp[[ivar]], vals = arr_towrite, start = c(1, 1), count = c(nx, ny))
      }

      ncdf4::nc_sync(ncnew)

    }
  }

  ## write attributes

  if(!is.null(title)) ncdf4::ncatt_put(ncnew, 0, "title", title)
  if(!is.null(author)) ncdf4::ncatt_put(ncnew, 0, "author", author)
  if(!is.null(description)) ncdf4::ncatt_put(ncnew, 0, "description", description)
  if(!is.null(institute)) ncdf4::ncatt_put(ncnew, 0, "institute", institute)
  if(!is.null(source)) ncdf4::ncatt_put(ncnew, 0, "source", source)

  ncdf4::nc_close(ncnew)
}
