#' Make Ecologically-based Low Fractal-Dimensional polygons
#'
#' Fires have key features that are different than "normal" ecological zonations.
#' These include: 1) the property that they spread horizontally across a landscape,
#' regardless of human-mapped ecological areas. 2) Furthermore, even when
#' homogeneous fire regions are mapped, the homogeneousness may be because of
#' fuel changes that occupy the average condition for an ecological area, not
#' necessarily the boundary conditions. 3) Small areas that encompass an area
#' that is only "several times larger than the largest fires" will be too small
#' to capture the potential for multiple very large fires.
#' 4) Furthermore, very large ecological
#' regions may encompass too much variation in fire-fuel-climate relationships
#' that any estimates of fire behaviour will not be accurate across the whole
#' region. 5) No matter what area is chosen to Thus, to estimate fire behaviour across a large region like Canada,
#' we need polygons that are not too small (>> largest fires), not too big (<<), are not constrained
#' by ecological zones defined by vegetation, have buffers.
#'
#' With these five points in mind, we define *Ecologically-based
#' Low Fractal-dimensional polygons, "ELFs", that are constrained by large
#' ecological zones (Canadian Ecozones), are "larg-ish", more "blobby" (low fractal dimension)
#' than other ecological zonation in mountains (which tend to follow elevation
#' contours), and are each buffered (default 20km).
#' @param x A polygon layer (`SpatVector` or `sf`) or gridded layer (`SpatRaster`)
#' that maps the entire area of interest for the ELFs. If missing, it will try to use
#' `scfmutils::prepInputsFireRegimePolys(type = "FRU")` if `scfmutils` is installed.
#' Currently, only options with a column called `FRU` will work.
#' @param desiredBuffer The distance in `m` that the buffers should be. Default is `20000`.
#' @param maxArea The area in `m^2` that will be split if the ecozones are
#' larger than this. Default is `2.4e+11` which is about 500km x 500km plus the buffer.
#' This means that the ELFs will be smaller than this, with lower limit being the
#' individual ecozones plus buffer.
#' @param destinationPath A path where any downloads should be put.
#' @param singleSpatVector Logical. If `TRUE`, the function returns a single
#'   `SpatVector` of all ELFs (rather than the per-ELF list described in
#'   `Value`). Default `FALSE`.
#' @param useCache Logical. If `TRUE`, then some internal steps will be cached.
#'   If `FALSE`, then some (or most) will not be cached.
#'
#' @return A list of length 2, with `rasWhole` and `rasCentered`. Each of these
#'   will have the complete list of ELFs. `rasWhole` is the whole map of Canada
#'   with the individual ELF highlighted, and `rasCentered` is *just* the individual ELF
#'   centred and reprojected so its projected parameters are: long0 is the x mid,
#'   and lat0 = ymin and lat1 = ymax of the individual ELF. This will create ELFs
#'   with the least amount of pixel deformation.
#' @export
#' @importFrom stats setNames
makeELFs <- function(x, desiredBuffer = 20000,
                     maxArea = 2.4e+11, destinationPath = ".", singleSpatVector = FALSE,
                     useCache = TRUE) {
  ## keep many small terra intermediates in memory rather than spawning scratch tifs in
  ## terra's tmpdir -- those scratch tifs are the dominant FD_SETSIZE-exhausting holders
  ## (see climateData::makeClimateCluster diagnostics)
  tO <- list(todisk = FALSE, memmax = -1)
  forExit <- terra::terraOptions(print = FALSE)[names(tO)]
  do.call(terra::terraOptions, tO)
  on.exit(do.call(terra::terraOptions, forExit), add = TRUE)


  if (missing(x)) {
    if (!requireNamespace("scfmutils")) stop("Please install PredictiveEcology/scfmutils ",
                                        "or supply a x")
    x <-
      scfmutils::prepInputsFireRegimePolys(type = "FRU", destinationPath = destinationPath)
  }
  digNFP <- reproducible::.robustDigest(x)
  if (is(x, "sf")) {
    dv <- terra::vect(x) |>
      reproducible::Cache(omitArgs = "x", .cacheExtra = digNFP,
                          useCache = useCache)
  } else if (is(x, "SpatRaster")) {
    dv <- {
      `>=`(x, 0) |>
      terra::as.polygons() |>
        terra::buffer(4*res(x)[1]) |>
        terra::buffer(-4*res(x)[1]) } |>
      Cache(.cacheExtra = digNFP, omitArgs = "x",
            .functionName = "Create polygon from input SpatRaster")
  }
  
  # remove tiny NA holes --> this should not exist because this is coming from SCANFI files;
  #    but there were tiny NA slivers that were being magnified by the aggregate for ELFs.
  # Raster-only: when `x` is (or defaults to) an sf of fire regime polygons there are no NA holes
  # to fill, and focal() has no method for it -- `no method found for signature sf, data.frame`.
  if (is(x, "SpatRaster")) {
    x <- {
      focal(x, w = 5, fun = "modal", na.policy = "only", na.rm = TRUE) |>
        terra::mask(dv)
    } |> Cache()
  }
  
  ecoNames <- c(
    "zone", "region",
    "province"
    , "district"
  )
  makeEcoURLs <- function(ecoNames) {
    vapply(ecoNames, FUN.VALUE = character(1), function(en)
      paste0(paste0("https://sis.agr.gc.ca/cansis/nsdb/ecostrat/", en, "/eco", en, "_shp.zip")))
  }
  urls <- makeEcoURLs(ecoNames)
  ecosLCC <- Map(nam = ecoNames, url = urls, function(nam, url) {
    reproducible::prepInputs(url = url,
                             destinationPath = destinationPath,
                             projectTo = x,
                             ## always a SpatVector: an sf gives a different 5 km grid below
                             fun = terra::vect,
                             useCache = useCache) |>
      reproducible::Cache(omitArgs = "projectTo",
                          .cacheExtra = list(x = digNFP),
                          useCache = useCache, .functionName = paste0("prepInputs_", nam))
  })
  tmpl <- terra::rast(ecosLCC[[1]], resolution = 5000)
  if (is(x, "SpatRaster")) {
    hf <- terra::project(x, tmpl)
  } else {
    hf <- terra::rasterize(x, tmpl, field = "FRU") |>
      reproducible::Cache(omitArgs = "x", .cacheExtra = digNFP,
                          useCache = useCache)
  }

  hf <- terra::trim(hf)

  ecosR <- Map(ecoLCC = ecosLCC, nam = names(ecosLCC),
               function(ecoLCC, nam) {
                 terra::rasterize(ecoLCC, tmpl, field = toupper(paste0("eco", substr(nam, 1, 7)))) |>
                   reproducible::postProcess(maskTo = hf)
               }) |>
    reproducible::Cache(omitArgs = "maskTo", .functionName = "rasterize_eco*_files",
                        .cacheExtra = list(hf = attr(hf, "tags")),
                        useCache = useCache)

  ecopR <- ecosR$province
  # Some are on the edge, e.g., in the tundra --> remove if less than 100 pixels
  fre <- terra::freq(ecopR)
  fre <- fre[fre$count > 100, ]
  ecopRseg <- terra::segregate(ecopR)
  if (terra::is.factor(ecopR)) {
    categories <- terra::cats(ecopR)[[1]]
    ## the first column holds the cell values, whatever its name: "ID" when made here, "value" when the
    ## raster comes back from a .tif (a cache restore)
    names(ecopRseg) <- categories[match(names(ecopRseg), categories[[1]]), "ECOPROVINC"]
  } else {
    stop("ecoprovince raster has no categories: segment names must be the '4.3'-style ecoprovince names, got: ",
         paste(head(names(ecopRseg), 5), collapse = ", "))
    categories <- data.frame(names(ecopRseg)) |> setNames(names(ecopR))
  }
  
  ord <- order(as.numeric(intersect(categories$ECOPROVINC, fre$value)))
  ecopRseg <- ecopRseg[[names(ecopRseg) %in% fre$value]][[ord]]

  # plotStack(ecopRseg)
  message("Starting to buffer ELFs")
  # out <- bufferOut(dv, desiredBuffer = desiredBuffer, # spatRas = x,
  #                  useCache = useCache) |>
  #   reproducible::Cache(omitArgs = "v", .cacheExtra = digNFP,
  #                       useCache = useCache)
  out2 <- mergeAndSplitRas(ecopRseg, ecosLCC$province, maxArea = maxArea,
                           useCache = useCache, destinationPath = destinationPath) |>
    reproducible::Cache(omitArgs = "ecopRseg",
                        .cacheExtra = list(attr(ecosR, "tags"), deparse(mergeAndSplitRas), deparse(split_poly)))

  out3 <- lapply(out2, function(x) as.list(segregateKeepNames(x, omitClasses = 0))) |>
    unlist(recursive = FALSE)
  names(out3) <- sapply(out3, names)

  ELFs <- bufferOut(spatRasSeg = out3, mask = hf,
                    desiredBuffer = desiredBuffer,
                    useCache = useCache)
  # ELFs <- bufferOut(spatRasSeg = out3, mask = out$rasWhole[[1]],
  #                   desiredBuffer = desiredBuffer,
  #                   useCache = useCache) |>
  #   reproducible::Cache(omitArgs = c("spatRasSeg", "mask"),
  #                       .cacheExtra = list(x = digNFP,
  #                                          attr(out2, "tags"),
  #                                          attr(out, "tags")
  #                       ),
  #                       useCache = useCache)
  ELFsPoly <- NULL
  # if (isTRUE(singleSpatVector)) {
  
  # ELFs <- makeELFs(destinationPath = "inputs") # recursive if singleSpatVector is TRUE
  whole <- ELFs$rasWhole$`4.1`
  whole[whole[] >= 0] <- 1
  outerbound <- terra::as.polygons(whole)
  
  allVec <- Map(id = ELFs$rasWhole, nam = names(ELFs$rasWhole), function(id, nam) {
    a <- id
    a[a[] == 0] <- NA
    vec <- terra::as.polygons(a)
    bb <- try(vec[,"ID"] <- nam)
    if (is(bb, "try-error"))
      stop("could not assign ID column '", nam, "' to the polygon layer: ", conditionMessage(attr(bb, "condition")))
    vec[, "buffer"] <- vec[, nam]
    vec[,nam] <- NULL
    vec
    # terra::plot(vec, col = c("turquoise", "yellow")[vec$`4.1`], alpha = 0.5)
  })
  ELFsPoly <- Reduce(rbind, allVec)
  
  # }
  
  # ELFs2 <- Map(e = ELFs, function(e) {
  #   ff <- Filenames(e)
  #   nzff <- nzchar(ff)
  #   if (any(nzff)) {
  #     whnzff <- which(nzff)
  #     e[whnzff] <- Map(r = e[whnzff], function(r) {r[] <- terra::values(r, mat = FALSE); r})
  #   }
  #   e
  # })
  
  append(ELFs, list(poly = ELFsPoly))
}


bufferOut <- function(v, spatRasSeg, spatRas, mask, field = "FRU", desiredBuffer = 20000,
                      useCache = TRUE) {

  ## see makeELFs(); applied here too so direct callers of bufferOut benefit
  tO <- list(todisk = FALSE, memmax = -1)
  forExit <- terra::terraOptions(print = FALSE)[names(tO)]
  do.call(terra::terraOptions, tO)
  on.exit(do.call(terra::terraOptions, forExit), add = TRUE)

  if (missing(spatRasSeg) && missing(spatRas)) {
    tmpl <- terra::rast(v, resolution = 5000)
    spatRas <- {tmpl |>
        terra::rasterize(v, y = _, field = "FRU", touches = TRUE)} |>
      reproducible::Cache(omitArgs = "x",
                          useCache = useCache)
  }
  if (missing(spatRasSeg)) {
    spatRasSeg <- terra::segregate(spatRas)
  }

  crsThis <- if (is(spatRasSeg, "list")) terra::crs(spatRasSeg[[1]]) else terra::crs(spatRasSeg)

  vCentred <- list()
  r <- list()
  ca <- list()
  if (!missing(v)) {
    vdf <- as.data.frame(v)
    layers <- sort(vdf[[field]])
  } else {
    layers <- names(spatRasSeg)
  }

  lostPixels <- list()
  for (i in layers) {
    if (missing(v)) {
      tmpsrs <- spatRasSeg[[i]]
      tmpsrs[!tmpsrs[] > 0] <- NA
      nf <- terra::as.polygons(tmpsrs)
      nf[[field]] <- names(tmpsrs)
    } else {
      nf <- v[vdf[[field]] %in% i, 1]
    }

    # center the individual polygon
    nfll <- terra::project(nf, "epsg:4943")
    exts <- round(terra::ext(nfll))[]
    middle <- mean(c(exts[["xmin"]],exts[["xmax"]]))
    prj <- paste0("+proj=lcc +lat_0=", exts[["ymin"]], " +lon_0=",middle," +lat_1=",exts[["ymax"]],
                  " +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
    vCentred[[i]] <- terra::project(nfll, prj)

    b41 <- vCentred[[i]]
    
    tmpl <- terra::rast(b41, resolution = 5000)
    largestBuffer <- 8e4
    tmpl1 <- terra::extend(tmpl, largestBuffer/terra::res(tmpl)[1])
    b41Orig <- terra::rasterize(b41, y = tmpl, field = field)
    b41Orig <- as.numeric(b41Orig) # remove factor
    b41Extended <- terra::rasterize(b41, y = tmpl1, field = field)
    b41Extended <- as.numeric(b41Extended) # remove factor
    b41ToUse <- b41Extended
    b41 <- b41ToUse
    if (is.factor(b41)) {
      cat <- cats(b41)[[1]][[1]]
      b41[!b41[] %in% cat] <- NA
    } else {
      # b41[!b41[] > as.integer(i)] <- NA
    }
    # do buffer on the centred polygon
    buffer_filled <- function(x, width, original = x) {
      # 1. data footprint as polygon, holes filled
      mask <- terra::ifel(is.na(x[[1]]), NA, 1L)
      p    <- terra::as.polygons(mask, dissolve = TRUE)
      p    <- terra::fillHoles(p)                  # <-- closes ALL internal holes
      
      # 2. buffer the filled footprint
      pb   <- terra::buffer(p, width = width)
      
      # 3. rasterize buffer as 1, then overwrite original footprint as 2
      template <- terra::rast(x)
      r    <- terra::rasterize(pb, template, field = 1L, background = NA)
      r[!is.na(original[[1]])] <- 2L               # original region = 2
      
      terra::trim(r)
    }
    b41ddd <- buffer_filled(b41, desiredBuffer)
    
    # terra::buffer is failing for me: it just makes whole map "TRUE"; can't diagnose; claude neither
    #   this is a workaround
    #terra::toMemory(b41)
    # b41ddd <- terra::buffer(b41, largestBuffer) # goes well past the edge of the map
    # b41ddd[b41ddd[] %in% TRUE] <- NA
    # b41ddd <- terra::buffer(b41ddd, largestBuffer - desiredBuffer)
    # b41ddd[] <- b41ddd[] %in% FALSE
    # b41ddd[b41ddd[] %in% FALSE] <- NA
    # b41ddd[] <- b41ddd + 0L # convert to numeric/integer
    # b41ddd[!is.na(b41ToUse[])] <- 2L
    # b41ddd <- terra::trim(b41ddd)

    names(b41ddd) <- i
    if (is.factor(b41Orig)) {
      cat <- cats(b41Orig)[[1]][[1]]
      b41OrigNew <- terra::resample(b41Orig, b41ddd, method = "near")
      b41ddd[b41OrigNew[] %in% cat & !is.na(b41OrigNew[])] <- 2
    } else {
      # b41ddd[!b41Orig[] > as.integer(i)] <- NA
    }

    r[[i]] <- b41ddd

    # Add it to the national map
    b41b <- terra::project(b41ddd, crsThis, method = "near")
    b41b <- terra::resample(b41b, spatRasSeg[[i]], method = "near")
    ## the core is exactly the ELF's own cells, so the cores tile the map without overlaps or gaps;
    ## the centred map projected back gives only the buffer ring ("near" shifts its edges a cell)
    ca[[i]] <- spatRasSeg[[i]]
    ca[[i]][b41b > 0] <- 1L
    ca[[i]][spatRasSeg[[i]] > 0] <- 2L
    if (missing(mask)) {
      mask <- spatRasSeg[[i]]
    }


    # need to mask 2x maybe
    r[[i]] <- reproducible::maskTo(r[[i]], mask, verbose = FALSE, touches = FALSE)
    if (NROW(terra::freq(is.na(r[[i]]))) > 1) {
      ca[[i]] <- reproducible::maskTo(ca[[i]], mask, verbose = FALSE, touches = FALSE)

      patchs <- terra::patches(r[[i]], allowGaps = FALSE, values = FALSE, directions = 8)
      tab <- terra::freq(patchs)

      if (NROW(tab) > 1) {
        tooSmall <- tab$count < 500
        if (!all(tooSmall)) {
          if (any(tooSmall)) {

            r[[i]][patchs[] %in% tab$value[tooSmall]] <- NA
            r[[i]] <- terra::trim(r[[i]])
            a <- ca[[i]]
            a[ca[[i]] == 0] <- NA
            patchsCA <- terra::patches(a, allowGaps = FALSE, values = FALSE, directions = 8)
            # patchsCA <- terra::patches(ca[[i]] > 0, allowGaps = FALSE, values = FALSE, directions = 8)
            tab <- terra::freq(patchsCA)
            tooSmall <- tab$count < 500
            wh <- which(patchsCA[] %in% tab$value[tooSmall])
            dtLost <- data.table(pixelID = wh, value = terra::values(ca[[i]])[wh])

            lostPixels[[i]] <- dtLost
            ca[[i]][wh] <- 0
          }
        }
      }
      # r[[i]] <- maskTo(r[[i]], mask, verbose = FALSE, touches = FALSE)
      # ca[[i]] <- maskTo(ca[[i]], mask, verbose = FALSE, touches = FALSE)
    }
    message("Done bufferOut for ", i)
  }

  ll <- moveSliversToOtherELFs(lostPixels, ca, i, r)
  # ll3 <- ll
  destinationPath <- unique(dirname(Filenames(spatRasSeg)))
  # There may be tempfiles along the way; remove them
  destinationPath <- grep(basename(tempdir()), destinationPath, invert = TRUE, value = TRUE)
  destinationPath <- grep(terra::terraOptions()$tempdir, destinationPath, invert = TRUE, value = TRUE)
  ELFpath <- file.path(destinationPath, "ELFs_final")
  unlink(ELFpath, recursive = TRUE)
  dir.create(ELFpath, recursive = TRUE, showWarnings = FALSE)
  message("Writing ELF rasters to disk")
  rr <- Map(out = ll, namOut = names(ll), function(out, namOut) {
    Map(inner = out, nam = names(out), function(inner, nam) {
      fn <- file.path(ELFpath, paste0(namOut, "_", nam, ".tif"))
      unlink(fn)
      terra::writeRaster(inner, filename = fn, overwrite = TRUE)
    })
  })
  list(rasCentered = rr$r, rasWhole = rr$ca)
}

segregateKeepNames <- function(ecopR, omitClasses, classes = NULL) {
  uniquVals <- unique(ecopR[]) |> na.omit()
  if (!missing(omitClasses)) {
    classes <- setdiff(uniquVals, omitClasses)
  }
  if (length(setdiff(uniquVals, 0)) > 1) {
    origFilenames <- Filenames(ecopR)
    ecopRseg <- terra::segregate(ecopR, classes = classes)
    if (is.factor(ecopR)) {
      names(ecopRseg) <- cats(ecopR)[[1]][match(as.numeric(names(ecopRseg)), cats(ecopR)[[1]][[1]]), 2]
      dt <- data.table(nam = names(ecopRseg), num = as.numeric(names(ecopRseg)))
    } else {
      nam <- names(ecopR)
      namOrig <- names(ecopRseg)
      names(ecopRseg) <- paste0(nam, ".", names(ecopRseg))
      dt <- data.table(nam = names(ecopRseg), num = as.numeric(namOrig))
    }
    data.table::setorder(dt, "num")
    ecopR <- ecopRseg[[match(dt$nam, names(ecopRseg))]]
  }
  if (terra::inMemory(ecopR))
    ecopR <- writeRaster(ecopR, overwrite = TRUE,
                         filename = file.path(dirname(origFilenames),
                                              paste0(paste(names(ecopR), collapse = "_"), ".tif")))
  ecopR

}

#' @importFrom sf st_as_sf st_sample st_geometry st_intersection st_area
split_poly <- function(sf_poly, n_areas) {
  # Create random points
  wasTerra <- FALSE
  if (is(sf_poly, "SpatVector")) {
    sf_poly <- sf::st_as_sf(sf_poly)
    wasTerra <- TRUE
  }
  points_rnd <- sf::st_sample(sf_poly, size = 10000)
  # k-means clustering
  points <- do.call(rbind, sf::st_geometry(points_rnd)) %>%
    as.data.frame() %>% setNames(c("lon","lat"))
  k_means <- kmeans(points, centers = n_areas)
  # Create voronoi polygons
  if (!requireNamespace("dismo"))
    stop("Please install `dismo` package to proceed")
  voronoi_polys <- dismo::voronoi(k_means$centers, ext = sf_poly)
  # Clip to sf_poly
  terra::crs(voronoi_polys) <- terra::crs(sf_poly)
  voronoi_sf <- sf::st_as_sf(voronoi_polys)
  equal_areas <- sf::st_intersection(voronoi_sf, sf_poly)
  equal_areas$area <- sf::st_area(equal_areas)
  if (wasTerra)
    equal_areas <- terra::vect(equal_areas)
  
  # Put them in xmin to xmax, ymin to ymax order
  mins <- Map(ind = seq(NROW(equal_areas)), function(ind) {
    cbind(xmin = terra::xmin(equal_areas[ind, ]), ymin = terra::ymin(equal_areas[ind, ]))
  }) |> do.call(args = _, rbind)
  ord <- order(mins[, "xmin"], mins[, "ymin"])
  equal_areas <- equal_areas[ord, ]
  equal_areas[, "id"] <- seq(NROW(equal_areas))
  
  return(equal_areas[ord, ])
}

mergeAndSplitRas <- function(ecopRseg, ecopLCC, maxArea = 2.4e+11,
                             field = "ECOPROVINC", useCache = TRUE,
                             destinationPath) { # FRU 27 ... or FRU 26 is 4.02075e+11
  # Require::Require(c(sf, dismo))
  ## split_poly() is random (sf::st_sample(), kmeans()). Set the RNG here so the split does not depend
  ## on what ran before. Seed 123 then 31 uniform draws is the state the ELFs used since 2026-10-05
  ## were split with, so those maps do not change.
  withr::local_seed(123, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                    .rng_sample_kind = "Rejection")
  stats::runif(31)
  # Each ecoprovince is written to its own .tif below, named by the province code ("4.1.tif").
  # Those used to go straight into destinationPath, scattering ~36 anonymous files through a shared
  # inputs directory and letting successive calls overwrite one another. Give each call its own
  # subdirectory, keyed on the inputs that determine the output.
  # NOTE: this protects the files in destinationPath only. reproducible's cache-RESTORE path collapses
  # these to `cachePath/<basename>.tif`, dropping the cacheId, so distinct calls still collide there;
  # that is a reproducible bug, not something this subdirectory can fix.
  # .robustDigest() on a list returns one digest PER ELEMENT, so collapse to a single hash.
  dpKey <- reproducible::.robustDigest(
    paste(unlist(reproducible::.robustDigest(list(ecopRseg, ecopLCC, maxArea, field))),
          collapse = "")
  )
  destinationPath <- reproducible::checkPath(
    file.path(destinationPath, paste0("ELFs_", substr(dpKey, 1, 10))), create = TRUE
  )
  prov1 <- as.character(sort(as.numeric(unique(names(ecopRseg)))))
  out <- Map(prov = prov1, function(prov) {
    provChar <- as.character(prov)
    tmp <- ecopLCC[ecopLCC$ECOPROVINC == provChar, ]
    if (!is(tmp, "SpatVector")) {
      tmp <- terra::vect(tmp)
    }
    tmp[, "AREA"] <- terra::expanse(tmp)
    needUpdate <- FALSE
    # merge multiples into 1
    if (NROW(tmp) > 1) {
      tmp <- terra::hull(tmp)
      terra::values(tmp) <- data.frame(AREA = terra::expanse(tmp))
      # tmp[, "AREA"] <- terra::expanse(tmp)
      needUpdate <- TRUE
      tmp[, field] <- paste0(provChar, seq_len(NROW(tmp)))
    }
    tooBig <- tmp$AREA > maxArea
    # split too large into 2 or more
    if (any(tooBig)) {
      numAreas <- tmp$AREA/ maxArea
      message("Splitting ", provChar, " in ", ceiling(numAreas))
      tmp <- split_poly(tmp, n_areas = ceiling(numAreas))
      tmp[, "AREA"] <- NULL
      needUpdate <- TRUE
      tmp[, field] <- tmp[, "id"]
    }
    if (needUpdate) {
      tmpRas <- try(terra::rasterize(tmp, ecopRseg[[provChar]], field = field, touches = TRUE))
      ## only this province's cells: the hull of a multi-polygon province (and a split piece's
      ## touching edge) also covers its neighbours' cells
      whChange <- tmpRas[] > 0 & ecopRseg[[provChar]][] > 0
      a <- ecopRseg[[provChar]]
      # Next line doesn't work --> it pulls it into memory
      # set.values(a, cells = which(whChange %in% 1), values(tmpRas, mat = FALSE)[whChange %in% 1])
      a[whChange] <- tmpRas[whChange]
    } else {
      a <- ecopRseg[[provChar]]
    }
    # Need to write it explicitly because it keeps the original "fullRaster" as its .tif.
    # destinationPath is the per-call subdirectory established above, so `provChar` alone is a
    # sufficient (and readable) filename within it.
    tmp <- writeRaster(a, filename = file.path(destinationPath, paste0(provChar, ".tif")),
                       overwrite = TRUE)
    message("Done ", prov)
    tmp
  })
  out
}

moveSliversToOtherELFs <- function(lostPixels, ca, i, r) {

  message("moving slivers to neighbouring ELF")
  if (NROW(unlist(lostPixels))) {
    if (is.null(names(lostPixels))) {
      hasNames <- FALSE
      haveLostPixels <- which(lengths(lostPixels) > 0)
    } else {
      hasNames <- TRUE
      haveLostPixels <- names(lostPixels)
    }
    for (lp in haveLostPixels) {
      numOverlap <- list()
      for (i in seq(ca)) {
        whPixelIDsHere <- which(ca[[i]][] > 0)
        lostPixelsElsewhere2 <- lostPixels[[lp]][value == 2]$pixelID %in% whPixelIDsHere
        lostPixelsElsewhere1 <- lostPixels[[lp]][value == 1]$pixelID %in% whPixelIDsHere
        if (any(lostPixelsElsewhere2)) {
          numOverlap[[i]] <- sum(lostPixelsElsewhere2)
          print(paste(i, ": ", sum(lostPixelsElsewhere2), "pixels"))
        }
      }

      # If the pixels that were lost do not show up in any other ELF blob, they are lost forever
      #   In the examples I saw, it was tiny islands off the coast of Newfoundland
      if (NROW(numOverlap)) {
        addTo <- which(sapply(numOverlap, function(x) !is.null(x)))
        if (hasNames) {
          addTo <- names(ca)[addTo]
        }

        a <- list()
        sums <- list()
        for (addToInd in addTo) {
          curPixelVal <- ca[[addToInd]][] != 2
          arr <- try(a[[addToInd]] <- ca[[addToInd]])
          if (is(arr, "try-error"))
            stop("could not build the pixel array: ", conditionMessage(attr(arr, "condition")))
          a[[addToInd]][lostPixels[[lp]]$pixelID] <- pmax(terra::values(a[[addToInd]])[lostPixels[[lp]]$pixelID], lostPixels[[lp]]$value, na.rm = TRUE)
          theA <- terra::freq(a[[addToInd]])
          # theA <- lapply(a, function(x) if (!is.null(x)) terra::freq(x))
          ar1 <- theA[-1, ]
          br1 <- terra::freq(ca[[addToInd]])[-1,]
          if (NROW(ar1) != NROW(br1)) {
            missings <- setdiff(br1$value, ar1$value)
            a2 <- data.frame(layer = ar1$layer[1], value = missings, count = br1$count[br1$value %in% missings])
            ar1 <- rbind(ar1, a2)
            ar1 <- ar1[order(ar1$value),]
          }
          sums[[addToInd]] <- sum(ar1 - br1)
        }
        addTo <- which.min(unlist(lapply(sums, function(ss) if (is.null(ss)) 1e7 else ss )))
        if (hasNames)
          addTo <- names(addTo)
        newVals <- pmax(terra::values(ca[[addTo]])[lostPixels[[lp]]$pixelID], lostPixels[[lp]]$value, na.rm = TRUE)
        ca[[addTo]][lostPixels[[lp]]$pixelID] <- newVals

        # dealt with the lp-th element in lostPixels; remove it: next line doesn't work. It shrinks the list.
        # lostPixels[[lp]] <- NULL
        a <- ca[[addTo]]
        a[] <- NA
        a[lostPixels[[lp]]$pixelID] <- newVals
        # a[a[] == 0] <- NA
        a <- terra::trim(a)
        a <- terra::project(a, terra::crs(r[[addTo]]), method = "near")
        bb <- terra::resample(a, r[[addTo]], method = "near")
        whVals <- which(terra::values(bb) > 0)
        r[[addTo]][whVals] <- terra::values(bb)[whVals]
      } else {
        numLostPixelsForever <- try(sum(lostPixels[[lp]]$value))
        if (is(numLostPixelsForever, "try-error"))
          stop("could not count lost pixels: ", conditionMessage(attr(numLostPixelsForever, "condition")))
        if (is.null(names(ca))) {
          message("From ELF ", lp, ", lost ", numLostPixelsForever," isolated pixels that do not exist in another ELF")
        }
      }
      nam <- if (is.null(names(lostPixels))) seq(lostPixels) else names(lostPixels)
      lostPixels <- Map(na = nam, function(na) if (na == lp) NULL else lostPixels[[na]])
    }
  }
  list(ca = ca, r = r)
}



#' Find a module in a project's module list, seeing through parent modules
#'
#' A project's `modules` may list a parent module (e.g. `PredictiveEcology/fireSense@development`)
#' instead of its children. This searches the listed modules and, recursively, the child modules
#' of any listed parent, read from `modulePath` with [SpaDES.core::moduleMetadata()].
#'
#' @param modules Character vector of module names or specs (`"owner/repo@branch"`).
#' @param modulePath Directory holding the module folders (children included).
#' @param pattern Regular expression matched against module names (no owner or branch).
#'
#' @return The plain name of the first matching module; stops if none matches.
#' @importFrom SpaDES.core moduleMetadata
#' @importFrom Require extractPkgName
#' @keywords internal
.findModuleInProject <- function(modules, modulePath, pattern) {
  toName <- function(x) unname(Require::extractPkgName(unname(as.character(x))))
  searched <- character(0)
  search <- function(mods) {
    for (m in setdiff(toName(mods), searched)) {
      searched <<- c(searched, m)
      if (grepl(pattern, m)) return(m)
      if (file.exists(file.path(modulePath, m, paste0(m, ".R")))) {
        kids <- SpaDES.core::moduleMetadata(module = m, path = modulePath,
                                            defineModuleListItems = "childModules")$childModules
        kids <- unlist(kids)
        if (length(kids)) {
          out <- search(kids)
          if (!is.null(out)) return(out)
        }
      }
    }
    NULL
  }
  out <- search(modules)
  if (is.null(out))
    stop("No module matching '", pattern, "' among the project's modules or their child modules; ",
         "searched: ", paste(searched, collapse = ", "), call. = FALSE)
  out
}

#' Run the ELFs module and return the ELF names or maps
#'
#' Extracting only the "ELF" module, runs it via \code{SpaDES.core::simInitAndSpades2()} and
#' returns the vector of polygon IDs from \code{spreadFitPreRun}, or the ELF maps. The run is not
#' cached here: the module caches the national ELF list itself, keyed on its \code{ELFversion}
#' parameter alone and shared through the cloud cache, so a run only downloads it.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Restricts \code{preRunSetupProject$modules} to entries matching \code{"ELFs"}.
#'   \item Initializes and runs the simulation with \code{SpaDES.core::simInitAndSpades2()}.
#'   \item Returns \code{sim$spreadFitPreRun[[fireSenseUtils::polygonIDTxt]]}, or the names or maps
#'         chosen by \code{whatOut}.
#' }
#'
#' The function expects \code{preRunSetupProject} to be a SpaDES project list with at least:
#' \itemize{
#'   \item \code{$modules}: character vector of module names.
#'   \item \code{$paths$modulePath}: base directory containing module folders.
#' }
#'
#' @param preRunSetupProject list. A SpaDES project object prepared for
#'   \code{SpaDES.core::simInitAndSpades2()}, containing \code{$modules} and
#'   \code{$paths$modulePath}.
#' @param urlELFresults Unused and deprecated; kept so existing calls still work. `runELFs()` no longer
#'   reads this Google Drive file or uploads to it.
#' @param whatOut string. One of `c("fittedNamesOnly", "allNames", "maps" )`. Partial
#'   matching is active (i.e., `fit` will do `fittedNamesOnly`). Will
#'   default to the first one if more than one supplied. If not `fittedNamesOnly` nor `allNames`,
#'   then will do `maps`. `fittedNamesOnly` will return the vector of the `.ELFind` names
#'   that have `SpreadFit`` values in the Cloud Geo Cache. `allNames` will return all .ELFind names,
#'   even if they have no fits. `maps` or not one of the first two options,
#'   will return a list of length 2, each with length `NROW(names(.ELFinds))`. These
#'   are individual rasters, one per .ELFind, either with the entire map of 
#'   Canada (`rasWhole`) or with a projection centered on the centroid of the ELF (`rasCentred`).
#' 
#'
#' @return
#' A vector corresponding to \code{sim$spreadFitPreRun[[fireSenseUtils::polygonIDTxt]]}.
#' Currently, no high arctic
#' ELFs are returned (e.g., Ecoprovinces starting with either 1. or 2.), nor the Prairie
#' grassland ELFs (ecoprovince 10.2, i.e. 10.2 and 10.2.*); see [ELFsNotFitted()].
#'
#' @examples
#' \dontrun{
#' # Assuming 'prj' is a valid SpaDES project list:
#' # prj <- somePreparationFunction(...)
#' ids <- runELFs(prj)
#' }
#'
#' @importFrom SpaDES.core simInitAndSpades2
#' @export
#'
#' @seealso
#' \code{\link[SpaDES.core]{simInitAndSpades2}}
#'
#' @author
#' Based on the provided code snippet.
#'
runELFs <- function(
    preRunSetupProject,
    whatOut = c("fittedNamesOnly", "allNames", "maps"),
    urlELFresults = paste0(
      "https://drive.google.com/file/d/",
      "1tkC944mPzR9-y-qCMDB5o2cAR_1MoDz4/view?usp=drive_link"
    )) {
  # The project may list a parent module (e.g. PredictiveEcology/fireSense@development) rather
  # than fireSense_ELFs itself; only the ELFs module is run here.
  preRunSetupProject$modules <- .findModuleInProject(
    preRunSetupProject$modules, preRunSetupProject$paths$modulePath, pattern = "^fireSense_ELFs$")
  sim <- SpaDES.core::simInitAndSpades2(preRunSetupProject)

  # Arctic and Prairie grassland ELFs are out permanently (ELFsNotFitted())
  arcticELFs <- .ELFsNotFittedPattern()
  # An ELF with no SCANFI tree species is NOT excluded: fireSense fits nonForest fuel
  # classes, and the dataPrep chain now carries speciesLayers = NULL and a zero-row
  # sppEquiv through (LandR loadSCANFISpeciesLayers/assertSpeciesLayers/sppHarmonize,
  # Biomass_speciesData, Biomass_borealDataPrep no-species branch, fireSenseUtils
  # cohortsToFuelClasses/assessFuelClasses, the Biomass_* init guards).
  # Excluding them here was a workaround for that chain stopping; the chain no longer stops.
  # ELFs with no fire in the fitted years, counted from the fire records by
  # fireSense_ELFs (ELFfitStatus() -> ELFsExcluded()). This is the same judgement
  # dataPrepFit makes only at the end of a cold run ("no ignitions present",
  # fireSense_dataPrepFit.R:1595-1597); making it here keeps them out of the queue.
  # Data-driven, and complementary to the Arctic pattern above: an Arctic
  # exclusion is not a statement about fire.
  zeroFireELFs <- if (is.null(sim$ELFsExcluded)) character(0) else sim$ELFsExcluded
  excludeELFs <- paste(c(arcticELFs,
                         if (length(zeroFireELFs))
                           paste0("^", gsub(".", "\\.", zeroFireELFs, fixed = TRUE), "$")),
                       collapse = "|")
  if (grepl("^fit|^all", whatOut[1])) {
    if (grepl("^all", whatOut[1])){
      .ELFinds <- names(sim$ELFs$rasCentered)
    } else {
      .ELFinds <- sim$spreadFitPreRun[[polygonIDTxt]]
    }
    .ELFinds <- grep(excludeELFs, invert = TRUE, value = TRUE, .ELFinds)
  } else {
    for (i in names(sim$ELFs)) {
      .ELFinds <- grep(excludeELFs, value = TRUE, names(sim$ELFs[[i]]))
      sim$ELFs[[i]][.ELFinds] <- NULL
    }
    .ELFinds <- sim$ELFs
  }

  #remove Arctic that is far from treeline
  .ELFinds
}



#' Identify and rasterize Ecologically-based Low Fractal-dimensional objects within a study area
#'
#' Intersects a study area with a set of Ecologically-based Low
#' Fractal-dimensional (ELF) polygons, identifies which ELFs overlap the study
#' area, and constructs a single categorical raster combining each ELF's
#' per-cell classification. A cell in an ELF's core (`2`) is labelled with that
#' ELF. A cell in the buffer (`1`) of one or more ELFs is labelled with the ELF,
#' among those, whose core is nearest. All other in-domain cells are `0`.
#'
#' The function assumes that `ELFs` (a list-like object with `rasWhole` rasters
#' and an `ID` field) and `rastTemplate` are available in the calling
#' environment, along with `inputPath`. ELF rasters are expected to contain the
#' integer values `0`, `1`, and `2`, where `2` indicates the core area of the ELF
#' and `1` indicates the buffer for that ELF.
#'
#' @param studyArea A `SpatVector` (or object coercible to one) defining the
#'   area of interest. Will be reprojected to the CRS of `ELFs` internally.
#' @param inputPath Character. Directory in which to store any downloads
#'   produced by `ELFtemplateRaster()` when `ELFsPolygon` is not supplied.
#' @param ELFsRaster Optional named list (one entry per ELF) of `SpatRaster`
#'   objects containing the integer values `0`, `1`, and `2`, accessed via
#'   `ELFsRaster$rasWhole`. If `NULL`, the per-ELF rasters must be loaded by
#'   the caller before invoking this function.
#' @param ELFsPolygon Optional `SpatVector` of ELF polygons. If `NULL`
#'   (the default), ELFs are generated by calling `makeELFs()` with a 20 km
#'   buffer around `rastTemplate`.
#'
#' @return A categorical `SpatRaster` with a `value`/`ELFind` levels table,
#'   where `ELFind` gives the ELF identifier (or `"none"` for background).
#'
#' @details
#' The three layers are combined with `terra::cover()` so that:
#' \itemize{
#'   \item cores override everything else,
#'   \item a buffer cell takes the ELF with the nearest core (`terra::distance()`), and
#'   \item the `0` background fills remaining in-domain cells.
#' }
#'
#' @seealso [makeELFs()], [terra::distance()], [terra::mosaic()],
#'   [terra::cover()]
#'
#' @export
ELFsInStudyArea <- function(studyArea, inputPath, ELFsRaster = NULL, ELFsPolygon = NULL) {
  if (is.null(ELFsPolygon)) {
    rastTemplate <- ELFtemplateRaster(inputPath)
    ELFsPolygon <- fireSenseUtils::makeELFs(rastTemplate, desiredBuffer = 20000, 
                                           destinationPath = inputPath, singleSpatVector = TRUE) |>
      Cache(omitArgs = "nationalForestPolygon",
            .cacheExtra = list(rt = attr(rastTemplate, "tags"),
                               makeELFs = fireSenseUtils::makeELFs,
                               bufferOutFn = bufferOut,
                               splitPoly = split_poly,
                               merge = mergeAndSplitRas))
    ELFsPolygon <- ELFsPolygon$poly
    
  }
  ELFsPolygon <- terra::intersect(ELFsPolygon, terra::project(studyArea, terra::crs(ELFsPolygon)))
  ELFsHere <- unique(ELFsPolygon$ID)
  
  
  # Pick the ones you want
  keep <- ELFsHere
  sub  <- ELFsRaster$rasWhole[keep]

  codes <- seq_along(sub)
  names(codes) <- keep

  # the ELF rasters on one grid: their union; 0 wherever any ELF has a value
  zero_bg <- terra::mosaic(terra::sprc(lapply(sub, function(r) r * 0)), fun = "min")
  onGrid <- terra::rast(lapply(sub, function(r) terra::crop(terra::extend(r, zero_bg), zero_bg)))

  # a core (2) is labelled with its ELF; cores do not overlap
  isCore <- onGrid == 2
  cores <- terra::mask(terra::which.max(isCore), any(isCore), maskvalues = FALSE)

  # a buffer (1) cell is labelled with the ELF, among those it is buffer of, whose core is nearest
  distToCore <- terra::rast(lapply(seq_along(sub), function(i) {
    terra::ifel(onGrid[[i]] == 1, terra::distance(terra::ifel(isCore[[i]], 1, NA)), NA)
  }))
  buffers <- terra::which.min(distToCore)

  # --- Stack with priority: cores > buffers > 0 background ---
  out <- terra::cover(cores, terra::cover(buffers, zero_bg))

  # Attach ELFind labels
  lvls <- data.frame(
    value  = c(0, codes),
    ELFind = c("none", names(codes)),
    stringsAsFactors = FALSE
  )
  levels(out) <- lvls
  ELFSout <- list(rast = out, poly = ELFsPolygon)
  ELFSout
}



#' Build a Canada-wide raster template for ELF construction
#'
#' Downloads the SCANFI Douglas-fir raster (by default the v2 2020 file on the
#' Arbutus object store), aggregates it by a factor of 8, turns it into a 1/`NA`
#' mask, loads it into memory, and writes it to disk as a reusable template for
#' building Ecologically-based Low Fractal-dimensional (ELF) objects. Both the
#' remote download and the aggregation step are cached via `reproducible::Cache()`,
#' keyed on `url` and the remote file's hash, so the template is only rebuilt when
#' the source changes.
#'
#' @param inputPath Character. Directory in which to store the downloaded
#'   source raster and the resulting `rastTemplate_Canada.tif`.
#' @param url Character. The source raster. It is used only for its footprint
#'   and grid, so any SCANFI layer with the national extent will do. The default
#'   is the v2 Douglas-fir 2020 file, which has the same extent, 30 m
#'   resolution and CRS as the v1.2 file on `ftp.maps.canada.ca` it replaced,
#'   and is smaller (163 MB vs 226 MB).
#'
#' @return A `SpatRaster` pointing to `rastTemplate_Canada.tif` inside
#'   `inputPath`. The raster covers Canada at roughly 8x the native SCANFI
#'   resolution; it is `1` where the aggregated source has data (including zero
#'   cover) and `NA` elsewhere, and is suitable as a template for downstream ELF
#'   construction.
#'
#' @details
#' The Douglas-fir layer is used purely for its spatial footprint and grid, so
#' the returned values are a mask, not cover. At `fact = 8`, the v1.2 file has
#' 156,026,794 non-`NA` cells and the v2 file 156,482,565; 460,121 cells are `NA`
#' only in v1.2 and 4,350 only in v2, so ELF maps and `rasterToMatch` built from
#' the template differ by about 0.3% of cells between the two.
#' Remote metadata is checked at most once per week
#' (`notOlderThan = Sys.time() - 60*60*24*7`) to avoid unnecessary network
#' calls.
#'
#' @seealso [reproducible::prepInputs()], [reproducible::Cache()],
#'   [terra::aggregate()], [terra::writeRaster()]
#'
#' @export
ELFtemplateRaster <- function(inputPath,
                              url = "https://object-arbutus.alliancecan.ca/529feda02d404836b0e91b1d360b4df0:predictiveecology/SCANFI_v2/2020/SCANFI_spsCC_PSEU_MEN_2020_v2_20260119.tif") {
  rastTemplate <- {
    hash <- utils::getFromNamespace("getRemoteMetadata", "reproducible")(isGDurl = FALSE, url = url) |>
      Cache(notOlderThan = Sys.time() - 60*60*24*7)
    {
      prepInputs(url = url, destinationPath = inputPath) |>
        terra::aggregate(fact = 8, na.rm = TRUE) |>
        is.na() |>
        terra::ifel(NA, 1) |>
        terra::toMemory()
    } |> Cache(omitArgs = TRUE, .cacheExtra = list(url, hash$remoteHash))
  }

  rastTemplate <- terra::writeRaster(
    rastTemplate, filename = file.path(inputPath, "rastTemplate_Canada.tif"), 
    overwrite = TRUE) 
  rastTemplate
  
}
