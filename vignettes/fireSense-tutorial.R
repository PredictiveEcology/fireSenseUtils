## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  eval = FALSE ## TODO: temporary until I can rework this vignette to actually run
)

library(reproducible)
Require(c("magrittr", "data.table", "dplyr", "pemisc", "PtProcess", "quickPlot",
          "raster", "rgeos", "sp", "SpaDES.core", "SpaDES.tools", "spatstat"))

options(reproducible.useCache = FALSE) # Turn off caching

moduleDir <- checkPath(file.path(tempdir(), "modules"), create = TRUE)
outputDir <- checkPath(file.path(tempdir(), "outputs"), create = TRUE)

setPaths(
  modulePath = moduleDir,
  outputPath = outputDir
)

## TODO: download modules to moduleDir -- use git or downloadModule() ???

normalizeRaster <- function(x) { ## TODO: use pemisc version
  # Normalize raster values
  (x - minValue(x)) / (maxValue(x) - minValue(x))
}

## ----create_data, message=FALSE-----------------------------------------------
#  # Create a raster template describing the study area
#  raster_width <- raster_height <- 100
#  r_template <- raster(matrix(nrow = raster_height, ncol = raster_width))
#  
#  # Create rasters describing environmental drivers controlling fire
#  # Here we use for example the land cover and the weather
#  
#  set.seed(123)
#  
#  ## Landcover -- these should be "proportion of that single cover type", e.g., Proportion Deciduous
#  landTypeOne <- gaussMap(r_template, scale = 300, var = 300)
#  landTypeOne <- normalizeRaster(landTypeOne)
#  landTypeTwo <- 1 - landTypeOne
#  
#  ## Weather
#  weather <- gaussMap(r_template, scale = 300, var = 3000)
#  
#  # Create low resolution rasters of the environmental drivers
#  # These are used for fitting the ignition model. The main advantages
#  # of fire modeling using low resolution are described in Marchal et al. 2017.
#  # In short, unlike high-resolution modeling, fire locations need not be known
#  # precisely and do not require data at scales where local homogeneity of
#  # landcover can be assumed.
#  
#  landTypeOneLowRes <- aggregate(landTypeOne, fact = 25, fun = mean)
#  landTypeTwoLowRes <- 1 - landTypeOneLowRes
#  
#  weatherLowRes <- aggregate(weather, fact = 25, fun = mean)
#  
#  # Let's use a known pattern / equation
#  nbFires <- 500
#  fireLocations <- rpoint(nbFires, f = function(x, y, ...) 300 * y) # Introduce some spatial variability: more fires in the North
#  fireLocations <- as(fireLocations, "SpatialPoints")
#  
#  #landcoverDataLowRes <- extract(landTypeOneLowRes, fireLocations, cellnumbers = TRUE) %>%
#  #  as_tibble %>%
#  #  rename(cell_id = 1, landtype_1_pp = 2) %>% # Rename columns using column index
#  #  distinct(cell_id, landtype_1_pp) %>%
#  #  mutate(landtype_2_pp = 1 - landtype_1_pp)  # Add a column with the 2nd landcover type
#  
#  #weatherDataLowRes <- extract(weatherLowRes, fireLocations, cellnumbers = TRUE) %>%
#  #  as_tibble %>%
#  #  rename(cell_id = 1, weather = 2) %>% # Rename columns using column index
#  #  distinct(cell_id, weather)
#  
#  # Create the dataset necessary to fit the statistical model of fire ignitions
#  #dataFireSense_IgnitionFit <- landcoverDataLowRes %>%
#  #  group_by(cell_id) %>%
#  #  summarise(n_fires = n()) %>%
#  #  right_join(landcoverDataLowRes) %>%
#  #  left_join(weatherDataLowRes) %>%
#  #  mutate(n_fires = ifelse(is.na(n_fires), 0, n_fires))
#  
#  # Converting the above from tibble to data.table
#  dataFireSense_IgnitionFit <- data.table(
#    extract(landTypeOneLowRes, fireLocations, cellnumbers = TRUE),
#    landtype_2_pp = extract(landTypeTwoLowRes, fireLocations, cellnumbers = FALSE),
#    weather = extract(weatherLowRes, fireLocations, cellnumbers = FALSE)
#  )
#  setnames(dataFireSense_IgnitionFit, old = "layer", new = "landtype_1_pp")
#  dataFireSense_IgnitionFit <- dataFireSense_IgnitionFit[, list(
#    n_fires = .N, landtype_1_pp = landtype_1_pp[1],
#    landtype_2_pp = landtype_2_pp[1], weather = weather[1]
#  ), by = cells]
#  
#  dataFireSense_EscapeFit <- data.table(
#    landtype_1_pp = extract(landTypeOne, fireLocations, cellnumbers = FALSE),
#    landtype_2_pp = extract(landTypeTwo, fireLocations, cellnumbers = FALSE),
#    weather = extract(weather, fireLocations, cellnumbers = FALSE)
#  )
#  
#  # Create the dataset necessary to fit the statistical model fire escapes
#  #dataFireSense_EscapeFit <- tibble(
#  #  landtype_1_pp = extract(landTypeOne, fireLocations),
#  #  landtype_2_pp = extract(landTypeTwo, fireLocations),
#  #  weather = extract(weather, fireLocations)
#  #)
#  
#  b1 <- .01
#  b2 <- -1
#  b3 <- 1
#  
#  escapeProb <- with(dataFireSense_EscapeFit, binomial()$linkinv(weather * b1 + landtype_1_pp * b2 + landtype_2_pp * b3))
#  escaped <- rbinom(n = nbFires, size = 1, prob = escapeProb)
#  
#  dataFireSense_EscapeFit[, `:=`(escaped = escaped, n_fires = 1)]
#  # dataFireSense_EscapeFit <- mutate(dataFireSense_EscapeFit, escaped = escaped, n_fires = 1)
#  
#  # Create the dataset necessary to fit the statistical model of the fire size distribution
#  # Calculate mean landcover and weather conditions around the locations of escaped fires
#  buf_around_loc_escaped <- gBuffer(fireLocations[as.logical(escaped)], byid = TRUE, width = .03)
#  
#  landtypeTypeOneBufMn <- extract(landTypeOne, buf_around_loc_escaped) %>%
#    lapply(mean) %>%
#    unlist
#  
#  landtypeTypeTwoBufMn <- extract(landTypeTwo, buf_around_loc_escaped) %>%
#    lapply(mean) %>%
#    unlist
#  
#  weatherBufMn <- extract(weather, buf_around_loc_escaped) %>%
#    lapply(mean) %>%
#    unlist
#  
#  dataFireSense_SizeFit <- data.table(
#    landtype_1_pp = landtypeTypeOneBufMn,
#    landtype_2_pp = landtypeTypeTwoBufMn,
#    weather = weatherBufMn
#  )
#  
#  b1_l <- -.03
#  b2_l <- 4
#  b3_l <- -2
#  
#  lambda <- with(dataFireSense_SizeFit, exp(weather * b1_l + landtype_1_pp * b2_l + landtype_2_pp * b3_l))
#  
#  b1_t <- .02
#  b2_t <- 4
#  b3_t <- 2
#  
#  theta <- with(dataFireSense_SizeFit, exp(weather * b1_t + landtype_1_pp * b2_t + landtype_2_pp * b3_t))
#  
#  fireSize <- round(rtappareto(n = length(lambda), lambda = lambda, theta = theta, a = 1))
#  
#  dataFireSense_SizeFit[, fire_size := fireSize]
#  # dataFireSense_SizeFit <- mutate(dataFireSense_SizeFit, fire_size = fireSize)
#  
#  # Create the fire attribute dataset that describes the starting locations
#  # and the size of the fires to be spread. This is needed to fit the statistical model of spread probabilities
#  fireAttributesFireSense_SpreadFit <- SpatialPointsDataFrame(fireLocations[as.logical(escaped)],
#                                                              data = data.frame(size = fireSize))

## ----fit_ignition_model, message=FALSE----------------------------------------
#  modules <- list("fireSense_IgnitionFit")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_IgnitionFit module inputs
#  objects <- list(dataFireSense_IgnitionFit = dataFireSense_IgnitionFit)
#  
#  # Define fireSense_IgnitionFit module parameters
#  formula <- formula(n_fires ~ landtype_1_pp:weather + landtype_2_pp:weather - 1)
#  family <- poisson(link = "identity")
#  ub <- list(coef = 1)
#  dataObjName <- "dataFireSense_IgnitionFit"
#  trace <- 1
#  iterDEoptim <- 60
#  iterNlminb <- 100
#  cores <- 50
#  
#  parameters <- list(
#    fireSense_IgnitionFit = list(
#      formula = formula,           # Formula of the statistical model
#      family = family,             # Distribution family, here the negative binomial distribution
#      ub = ub,                     # Upper bounds for coefficients to be estimated
#      data = dataObjName,          # Name of the data.frame containing the variables in the statistical model. By default "dataFireSense_IgnitionFit"
#      trace = trace,               # Print progress every 10 iterations
#      iterDEoptim = iterDEoptim,   # Integer defining the maximum number of iterations allowed for the DEoptim optimizer
#      iterNlminb = iterNlminb,     # Number of trials, or searches, to be performed by the nlminb optimizer in order to find the best solution
#      cores = cores              # Number of logical cores to use for parallel computation during optimization
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    params = parameters,
#    objects = objects,
#    times = times
#  )
#  
#  fireSense_IgnitionFitted <- sim$fireSense_IgnitionFitted # Extract the fitted model from the sim object

## ----fit_escape_model, message=FALSE------------------------------------------
#  modules <- list("fireSense_EscapeFit")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_EscapeFit module inputs
#  objects <- list(dataFireSense_EscapeFit = dataFireSense_EscapeFit)
#  
#  # Define fireSense_EscapeFit module parameters
#  formula <- formula(cbind(escaped, n_fires - escaped) ~ landtype_1_pp + landtype_2_pp + weather - 1)
#  dataObjName <- "dataFireSense_EscapeFit"
#  
#  parameters <- list(
#    fireSense_EscapeFit = list(
#      formula = formula, # Formula of the statistical model
#      data = dataObjName # Name of the data.frame containing the variables in the statistical model. By default "dataFireSense_EscapeFit"
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    params = parameters,
#    objects = objects,
#    times = times
#  )
#  
#  fireSense_EscapeFitted <- sim$fireSense_EscapeFitted # Extract the fitted model from the sim object

## ----fit_size_distribution_model, message=FALSE-------------------------------
#  modules <- list("fireSense_SizeFit")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_SizeFit module inputs
#  objects <- list(dataFireSense_SizeFit = dataFireSense_SizeFit)
#  
#  # Define fireSense_SizeFit module parameters
#  formula <- formula(fire_size ~ landtype_1_pp + landtype_2_pp + weather - 1)
#  dataObjName <- "dataFireSense_SizeFit"
#  
#  parameters <- list(
#    fireSense_SizeFit = list(
#      formula = list(     # Formulas of the statistical model
#        beta = formula,   # The formulas for the beta and theta parameters of the tapered Pareto.
#        theta = formula   # They do not have to be identical, even if it's the case here
#      ),
#      data = dataObjName, # Name of the data.frame containing the variables in the statistical model. By default "dataFireSense_SizeFit"
#      link = list(
#        beta = "log",     # Link function for the beta parameter of the tapered Pareto
#        theta = "log"     # Link function for the theta parameter of the tapered Pareto
#      ),
#      a = 1,              # Lower truncation point a of the tapered Pareto
#      ub = list(
#        beta = 10,
#        theta = 10
#      ),
#      itermax = 2000,
#      trace = 100
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    params = parameters,
#    objects = objects,
#    times = times
#  )
#  
#  fireSense_SizeFitted <- sim$fireSense_SizeFitted # Extract the fitted model from the sim object

## ----produce_maps_tp_params, message=FALSE------------------------------------
#  modules <- list("fireSense_SizePredict")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_SizePredict module inputs
#  objects <- list(
#    fireSense_SizeFitted = fireSense_SizeFitted,
#    landtype_1_pp = landTypeOne,
#    landtype_2_pp = landTypeTwo,
#    weather = weather
#  )
#  
#  # Define fireSense_SizePredict module outputs
#  outputs <- rbind(
#    data.frame(
#      file = "fireSense_SizePredicted_Beta.tif",
#      fun = "writeRaster",
#      objectName = "fireSense_SizePredicted_Beta",
#      package = "raster",
#      saveTime = 1
#    ),
#    data.frame(
#      file = "fireSense_SizePredicted_Theta.tif",
#      fun = "writeRaster",
#      objectName = "fireSense_SizePredicted_Theta",
#      package = "raster",
#      saveTime = 1
#    )
#  )
#  
#  # Define fireSense_SizePredict module parameters
#  parameters <- list(
#    fireSense_SizePredict = list(
#      data = c("landtype_1_pp", "landtype_2_pp", "weather"),
#      modelName = "fireSense_SizeFitted" # This is the default
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    objects = objects,
#    outputs = outputs,
#    params = parameters,
#    times = times
#  )

## ----fit_spread_model, message=FALSE------------------------------------------
#  modules <- list("fireSense_SpreadFit")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_SpreadFit module inputs
#  inputs <- rbind(
#    # tapered Pareto's beta
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_SizePredicted_Beta",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "beta",
#      loadTime = 1
#    ),
#    # tapered Pareto's theta
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_SizePredicted_Theta",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "theta",
#      loadTime = 1
#    )
#  )
#  
#  objects <- list(fireAttributesFireSense_SpreadFit = fireAttributesFireSense_SpreadFit)
#  
#  # Define fireSense_SpreadFit module parameters
#  formula <- formula(~ I(1/beta) + log(theta) - 1)
#  
#  parameters <- list(
#    fireSense_SpreadFit = list(
#      formula = formula, # Formula of the statistical model
#      data = c("beta", "theta"),
#      lower = c(.01, 0, .1, .3, .001, .001),
#      upper = c(.20, .1, 10, 4., .300, .300),
#      cores = 8
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    inputs = inputs,
#    modules = modules,
#    objects = objects,
#    params = parameters,
#    times = times
#  )
#  
#  fireSense_SpreadFitted <- sim$fireSense_SpreadFitted # Extract the fitted model from the sim object

## ----predict_ignition_rates, message=FALSE------------------------------------
#  modules <- list("fireSense_IgnitionPredict")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_IgnitionPredict module inputs
#  objects <- list(
#    fireSense_IgnitionFitted = fireSense_IgnitionFitted,
#    landtype_1_pp = landTypeOne,
#    landtype_2_pp = landTypeTwo,
#    weather = weather
#  )
#  
#  # Define fireSense_IgnitionPredict module outputs
#  outputs <- rbind(
#    data.frame(
#      file = paste0("fireSense_IgnitionPredicted.tif"),
#      fun = "writeRaster",
#      objectName = "fireSense_IgnitionPredicted",
#      package = "raster",
#      saveTime = 1
#    )
#  )
#  
#  # Define fireSense_IgnitionPredict module parameters
#  parameters <- list(
#    fireSense_IgnitionPredict = list(
#      data = c("landtype_1_pp", "landtype_2_pp", "weather"),
#      modelName = "fireSense_IgnitionFitted" # This is the default
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    objects = objects,
#    outputs = outputs,
#    params = parameters,
#    times = times
#  )

## ----predict_escape_probabilities, message=FALSE------------------------------
#  modules <- list("fireSense_EscapePredict")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_EscapePredict module inputs
#  objects <- list(
#    fireSense_EscapeFitted = fireSense_EscapeFitted,
#    landtype_1_pp = landTypeOne,
#    landtype_2_pp = landTypeTwo,
#    weather = weather
#  )
#  
#  # Define fireSense_EscapePredict module outputs
#  outputs <- rbind(
#    data.frame(
#      file = paste0("fireSense_EscapePredicted.tif"),
#      fun = "writeRaster",
#      objectName = "fireSense_EscapePredicted",
#      package = "raster",
#      saveTime = 1
#    )
#  )
#  
#  # Define fireSense_EscapePredict module parameters
#  parameters <- list(
#    fireSense_EscapePredict = list(
#      data = c("landtype_1_pp", "landtype_2_pp", "weather"),
#      modelName = "fireSense_EscapeFitted" # This is the default
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    objects = objects,
#    outputs = outputs,
#    params = parameters,
#    times = times
#  )

## ----predict_spread_probabilities, message=FALSE------------------------------
#  modules <- list("fireSense_SpreadPredict")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense_SpreadPredict module inputs
#  # tapered Pareto's beta
#  inputs <- rbind(
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_SizePredicted_Beta",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "beta",
#      loadTime = 1
#    ),
#    # tapered Pareto's theta
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_SizePredicted_Theta",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "theta",
#      loadTime = 1
#    )
#  )
#  
#  objects <- list(fireSense_SpreadFitted = fireSense_SpreadFitted)
#  
#  # Define fireSense_SpreadPredict module outputs
#  outputs <- rbind(
#    data.frame(
#      file = paste0("fireSense_SpreadPredicted.tif"),
#      fun = "writeRaster",
#      objectName = "fireSense_SpreadPredicted",
#      package = "raster",
#      saveTime = 1
#    )
#  )
#  
#  # Define fireSense_SpreadPredict module parameters
#  parameters <- list(
#    fireSense_SpreadPredict = list(
#      data = c("landtype_1_pp", "landtype_2_pp", "weather"),
#      modelName = "fireSense_SpreadFitted" # This is the default
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    modules = modules,
#    objects = objects,
#    outputs = outputs,
#    params = parameters,
#    times = times
#  )

## ----run_fireSense, message=FALSE---------------------------------------------
#  modules <- list("fireSense")
#  
#  times <- list(start = 1, end = 1)
#  
#  # Define fireSense module inputs
#  inputs <- rbind(
#    # Fire ignition rates
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_IgnitionPredicted",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "ignitionProbRaster",
#      loadTime = 1
#    ),
#    # Fire escape probabilities
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_EscapePredicted",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "escapeProbRaster",
#      loadTime = 1
#    ),
#    # Fire spread probabilities
#    data.frame(
#      files = dir(
#        getPaths()$outputPath,
#        pattern = "fireSense_SpreadPredicted",
#        full.names = TRUE
#      ),
#      functions = "raster::raster",
#      objectName = "spreadProbRaster",
#      loadTime = 1
#    )
#  )
#  
#  # Run the simulation
#  sim <- simInitAndSpades(
#    inputs = inputs,
#    modules = modules,
#    objects = objects,
#    times = times
#  )
#  
#  burnMap <- sim$burnMap
#  
#  Plot(burnMap)

