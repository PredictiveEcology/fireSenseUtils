## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  eval = FALSE ## TODO: re-enable, but note that vignette takes a while to build
)

## ----packages-----------------------------------------------------------------
#  library(data.table)
#  library(fasterize)
#  library(ggplot2)
#  library(raster)
#  library(reproducible)
#  library(sf)

## ----compare_LCC_and_fire_data------------------------------------------------
#  dPath <- tempdir()
#  lccUrl <- paste0("ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/",
#                   "LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip")
#  lcc <- prepInputs(url = lccUrl, destinationPath = dPath,
#                    targetFile = "LCC2005_V1_4a.tif", alsoExtract = NA)
#  
#  fireUrl <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"
#  firePolygons <- prepInputs(url = fireUrl, destinationPath = dPath, fun = "sf::read_sf")
#  firePolygons <- sf::st_transform(firePolygons, crs = crs(lcc))
#  
#  #rasterize firePolygons
#  firePolygons$dummyVar <- 1
#  firePolygons <- firePolygons[firePolygons$YEAR > 2004, ]
#  fireRaster <- fasterize(sf = firePolygons, raster = lcc, field = "dummyVar")
#  fireLoc <- 1:ncell(lcc)
#  fireLoc <- fireLoc[!is.na(getValues(fireRaster))]
#  
#  burned <- data.table(pixelID = 1:ncell(lcc), lcc = getValues(lcc))
#  burned[pixelID %in% fireLoc, burn := 1]
#  burnCalc <- burned[, .(available = .N, burned = sum(burn, na.rm = TRUE)), .(lcc)]
#  
#  burnCalc[, "percentBurned" := round(burned/available * 100, digits = 3)]
#  setkey(burnCalc, lcc)
#  burnCalc

## ----fire, echo=FALSE---------------------------------------------------------
#  ggplot(data = burnCalc, aes(x = lcc, y = percentBurned)) +
#    geom_bar(stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
#    ylab("percent burned (%)") +
#    scale_x_continuous("LCC", labels = as.character(burnCalc$lcc), breaks = burnCalc$lcc) +
#    theme_bw()

## ----eval = FALSE-------------------------------------------------------------
#  ## original fireSense_dataPrepFit defaults
#  LCC2005_nonFlam <- c(0, 25, 30, 33, 36, 37, 38, 39)

