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
makeELFs <- function(x, desiredBuffer = 20000,
                     maxArea = 2.4e+11, destinationPath = ".", singleSpatVector = FALSE,
                     useCache = TRUE) {
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
                             useCache = useCache) |>
      reproducible::Cache(omitArgs = "projectTo",
                          .cacheExtra = list(x = digNFP),
                          useCache = useCache, .functionName = paste0("prepInputs_", nam))
  })
  tmpl <- terra::rast(ecosLCC[[1]], res = 5000)
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
  categories <- terra::cats(ecopR)[[1]]
  names(ecopRseg) <- categories[match(names(ecopRseg), categories$ID), "ECOPROVINC"]
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
    reproducible::Cache(omitArgs = "ecopRseg", .cacheExtra = attr(ecosR, "tags"))

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
  if (isTRUE(singleSpatVector)) {

    ELFs <- makeELFs(destinationPath = "inputs")
    whole <- ELFs$rasWhole$`4.1`
    whole[whole[] >= 0] <- 1
    outerbound <- terra::as.polygons(whole)

    allVec <- Map(id = ELFs$rasWhole, nam = names(ELFs$rasWhole), function(id, nam) {
      a <- id
      a[a[] == 0] <- NA
      vec <- terra::as.polygons(a)
      vec[,"ID"] <- nam
      vec[, "buffer"] <- vec[, nam]
      vec[,nam] <- NULL
      vec
      # terra::plot(vec, col = c("turquoise", "yellow")[vec$`4.1`], alpha = 0.5)
    })
    ELFs <- Reduce(rbind, allVec)

  }
  
  # ELFs2 <- Map(e = ELFs, function(e) {
  #   ff <- Filenames(e)
  #   nzff <- nzchar(ff)
  #   if (any(nzff)) {
  #     whnzff <- which(nzff)
  #     e[whnzff] <- Map(r = e[whnzff], function(r) {r[] <- terra::values(r, mat = FALSE); r})
  #   }
  #   e
  # })
  
  ELFs
}


bufferOut <- function(v, spatRasSeg, spatRas, mask, field = "FRU", desiredBuffer = 20000,
                      useCache = TRUE) {

  if (missing(spatRasSeg) && missing(spatRas)) {
    tmpl <- terra::rast(v, res = 5000)
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
    tmpl <- terra::rast(b41, res = 5000)
    largestBuffer <- 8e4
    tmpl1 <- terra::extend(tmpl, largestBuffer/terra::res(tmpl)[1])
    b41Orig <- terra::rasterize(b41, y = tmpl, field = field)
    b41Extended <- terra::rasterize(b41, y = tmpl1, field = field)
    b41ToUse <- b41Extended


    b41 <- b41ToUse
    if (is.factor(b41)) {
      cat <- cats(b41)[[1]]$ID
      b41[!b41[] %in% cat] <- NA
    } else {
      # b41[!b41[] > as.integer(i)] <- NA
    }
    # do buffer on the centred polygon

    b41 <- terra::buffer(b41, largestBuffer) # goes well past the edge of the map
    b41[b41[] %in% TRUE] <- NA
    b41 <- terra::buffer(b41, largestBuffer - desiredBuffer)
    b41[] <- b41[] %in% FALSE
    b41[b41[] %in% FALSE] <- NA
    b41[] <- b41 + 0L # convert to numeric/integer
    # b41[b41ToUse[] %in% as.integer(i)] <- 2L
    b41[!is.na(b41ToUse[])] <- 2L
    b41 <- terra::trim(b41)

    names(b41) <- i
    if (is.factor(b41Orig)) {
      cat <- cats(b41Orig)[[1]]$ID
      b41OrigNew <- terra::resample(b41Orig, b41, method = "near")
      b41[b41OrigNew[] %in% cat & !is.na(b41OrigNew[])] <- 2
    } else {
      # b41[!b41Orig[] > as.integer(i)] <- NA
    }

    r[[i]] <- b41

    # Add it to the national map
    b41b <- terra::project(b41, crsThis, method = "near")
    b41b <- terra::resample(b41b, spatRasSeg[[i]], method = "near")
    ca[[i]] <- spatRasSeg[[i]]
    ca[[i]][b41b > 0] <- b41b
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
    print(paste0("Done ", i))
  }

  ll <- moveSliversToOtherELFs(lostPixels, lp, ca, i, r)
  # ll3 <- ll
  destinationPath <- unique(dirname(Filenames(spatRasSeg)))
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
      names(ecopRseg) <- cats(ecopR)[[1]][match(as.numeric(names(ecopRseg)), cats(ecopR)[[1]]$ID), 2]
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
      whChange <- tmpRas[] > 0
      a <- ecopRseg[[provChar]]
      # Next line doesn't work --> it pulls it into memory
      # set.values(a, cells = which(whChange %in% 1), values(tmpRas, mat = FALSE)[whChange %in% 1])
      a[whChange] <- tmpRas[whChange]
    } else {
      a <- ecopRseg[[provChar]]
    }
    # Need to write it explicitly because it keeps the original "fullRaster" as its .tif
    tmp <- writeRaster(a, filename = file.path(destinationPath, paste0(provChar, ".tif")),
                       overwrite = TRUE)
    message("Done ", prov)
    tmp
  })
  out
}

moveSliversToOtherELFs <- function(lostPixels, lp, ca, i, r) {

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
          if (is(arr, "try-error")) browser()
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
        if (is(numLostPixelsForever, "try-error")) browser()
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



#' Run ELFs modules and optionally update googledrive png
#'
#' Extracting only the "ELF" module, runs the projectmodule
#' via \code{SpaDES.core::simInitAndSpades2()} with caching (optionally
#' updates a Google Drive file if the current user is \code{"emcintir"}), and
#' returns the vector of polygon IDs from \code{spreadFitPreRun}.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Restricts \code{preRunSetupProject$modules} to entries matching \code{"ELFs"}.
#'   \item Collects all \code{.R} source files under each module directory (recursively),
#'         excluding any under \code{tests/}, and uses these as part of the cache key.
#'   \item Retrieves remote metadata (hash) from a Google Drive file URL and assigns
#'         the remote hash to \code{preRunSetupProject$params$fireSense_ELFs$hashSpreadFitRemoteFile}.
#'   \item Initializes and runs the simulation with \code{SpaDES.core::simInitAndSpades2()},
#'         wrapping the call in \code{Cache()} so that the run is skipped unless inputs
#'         (source files and remote hash) have changed.
#'   \item If the current \code{SpaDES.project::user()} is \code{"emcintir"}, updates
#'         the Google Drive file at the provided URL with any output files whose names
#'         contain \code{"ELF"}, caching the update by \code{cacheId(sim)}.
#'   \item Returns \code{sim$spreadFitPreRun$polygonID}.
#' }
#'
#' Internal `Cache`ing relies on two extra keys:
#' \itemize{
#'   \item \code{src}: paths to module \code{.R} files (excluding \code{tests/}).
#'   \item \code{remoteHash}: the remote metadata hash obtained from the URL.
#' }
#'
#' The function expects \code{preRunSetupProject} to be a SpaDES project list with at least:
#' \itemize{
#'   \item \code{$modules}: character vector of module names.
#'   \item \code{$paths$modulePath}: base directory containing module folders.
#'   \item \code{$params$fireSense_ELFs}: a list where \code{hashSpreadFitRemoteFile} can be set.
#' }
#'
#' @param preRunSetupProject list. A SpaDES project object prepared for
#'   \code{SpaDES.core::simInitAndSpades2()}, containing \code{$modules},
#'   \code{$paths$modulePath}, and \code{$params$fireSense_ELFs}.
#' @param urlELFresults A googledrive url where the ELF data.frame that contains
#'   geometries and fitted SpreadFit parameters and several other 
#'   attributes.
#' @param onlyFittedELFs logical. If `TRUE`, the default, then this will only
#'   return the ELF indices that have been successfully fitted via `fireSense_SpreadFit`.
#'   Whether individual ELFs have been fitted or not will be determined based on
#'   downloading the `urlELFresults`.
#'
#' @return
#' A vector corresponding to \code{sim$spreadFitPreRun$polygonID}. Currently, no high arctic
#' ELFs are returned (e.g., Ecoprovinces starting with either 1. or 2.)
#'
#' @section Side Effects:
#' Updates a Google Drive file when the current user is \code{"emcintir"} and when
#' the run produces outputs with names containing \code{"ELF"}.
#'
#' @section Notes:
#' \itemize{
#'   \item The remote file URL is hard-coded as:
#'         \code{"https://drive.google.com/file/d/1tkC944mPzR9-y-qCMDB5o2cAR_1MoDz4/view?usp=drive_link"}.
#'   \item \code{Cache()} and \code{cacheId()} are assumed to be available (typically from
#'         \pkg{reproducible} in SpaDES workflows).
#'   \item \code{asPath()} is expected to be available in your environment (commonly from
#'         \pkg{reproducible}); adjust if your project defines it elsewhere.
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'prj' is a valid SpaDES project list:
#' # prj <- somePreparationFunction(...)
#' ids <- runELFs(prj)
#' }
#'
#' @importFrom SpaDES.core simInitAndSpades2
#' @importFrom reproducible Cache
#' @importFrom googledrive drive_update
#' @importFrom SpaDES.project user
#' @export
#'
#' @seealso
#' \code{\link[SpaDES.core]{simInitAndSpades2}},
#' \code{\link[reproducible]{Cache}},
#' \code{\link[googledrive]{drive_update}}
#'
#' @author
#' Based on the provided code snippet.
#'
runELFs <- function(preRunSetupProject, onlyFittedELFs = TRUE,
                    urlELFresults = "https://drive.google.com/file/d/1tkC944mPzR9-y-qCMDB5o2cAR_1MoDz4/view?usp=drive_link") {
  preRunSetupProject$modules <- grep("ELFs", preRunSetupProject$modules, value = TRUE)
  # For digesting; i.e., whether it needs to be re-run
  srcFiles <- asPath(dir(file.path(preRunSetupProject$paths$modulePath, preRunSetupProject$modules), 
                         pattern = "\\.R$", recursive = TRUE, full.names = TRUE) |> 
                       grep(pattern = "tests\\/", invert = TRUE, value = TRUE))
  md <- reproducible:::getRemoteMetadata(url = urlELFresults)
  preRunSetupProject$params$fireSense_ELFs$hashSpreadFitRemoteFile <- md$remoteHash
  
  # Run the module
  sim <- SpaDES.core::simInitAndSpades2(preRunSetupProject) |> 
    Cache(.cacheExtra = list(src = srcFiles, remoteHash = md$remoteHash),
          omitArgs = "l")
  
  # Upload if it is emcintir; and if it has changed
  if (SpaDES.project::user() %in% "emcintir") {
    cid <- cacheId(sim)
    ll <- googledrive::drive_update(file = urlELFresults,
                                    media = grep("ELF", sim@outputs$file, value = TRUE)) |> 
      Cache(omitArgs = c("file", "media"), .cacheExtra = list(cacheId = cid))
  }
  # })
  if (isTRUE(onlyFittedELFs)) {
    .ELFinds <- sim$spreadFitPreRun$polygonID
  } else {
    .ELFinds <- names(sim$ELFs$rasCentered)
  }
  #remove Arctic that is far from treeline
  .ELFinds <- grep("^1\\.|^2\\.", invert = TRUE, value = TRUE, .ELFinds)

  .ELFinds
}
