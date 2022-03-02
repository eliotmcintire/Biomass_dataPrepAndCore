## install/load required packages
if (!require("Require")) {install.packages("Require"); require("Require")}
Require("PredictiveEcology/SpaDES.install")
installSpaDES()
Require("SpaDES.core")

## environment setup -- all the functions below rely on knowing where modules
## are located via this command
setPaths(cachePath = "cache",
         inputPath = "inputs",
         modulePath = "modules",
         outputPath = "outputs")

## modules
moduleGitRepos <- c('Biomass_borealDataPrep',
                    'Biomass_core')
getModule(moduleGitRepos, overwrite = FALSE)
modules <- extractPkgName(moduleGitRepos)

## packages that are required by modules
makeSureAllPackagesInstalled()
Require::Require(c("PredictiveEcology/LandR@development (>= 1.0.7.9007)"))

## Module documentation -- please go to these pages to read about each module
##  In some cases, there will be important defaults that a user should be aware of
##  or important objects (like studyArea) that may be essential
## browseURL('https://github.com/PredictiveEcology/Biomass_borealDataPrep/blob/main/Biomass_borealDataPrep.rmd')
## browseURL('https://github.com/PredictiveEcology/Biomass_core/blob/main/Biomass_core.rmd')


## simulation initialization. These may not be appropriate start and end times for one
## or more of the modules, e.g., they may only be defined with calendar dates
simTimes <- list(start = 2011, end = 2021)



fixRTM <- function(x) {
  x <- raster::raster(x)
  x[!is.na(x[])] <- 1
  RIArtm3 <- terra::rast(x)
  aaa <- terra::focal(RIArtm3, fun = "sum", na.rm = TRUE, w = 5)
  RIArtm2 <- raster::raster(x)
  RIArtm2[aaa[] > 0] <- 1
  RIArtm4 <- terra::rast(RIArtm2)
  bbb <- terra::focal(RIArtm4, fun = "sum", na.rm = TRUE, w = 5)
  ccc <- raster::raster(bbb)[] > 0 & !is.na(x[])
  RIArtm2[ccc] <- 1
  RIArtm2[!ccc & !is.na(x[])] <- 0
  sa <- sf::st_as_sf(stars::st_as_stars(RIArtm2), as_points = FALSE, merge = TRUE)
  sa <- sf::st_buffer(sa, 0)
  sa <- sf::as_Spatial(sa)
  return(sa)
}
SA_ERIntersect <- function(x, studyArea) {
  x <- sf::st_read(x)
  sa_sf <- sf::st_as_sf(studyArea)
  ecoregions <- sf::st_transform(x, sf::st_crs(sa_sf))
  studyAreaER <- sf::st_intersects(ecoregions, sa_sf, sparse = FALSE)
  sf::as_Spatial(ecoregions[studyAreaER,])
}
# googledrive::drive_auth(cache = ".secret")
studyArea <- Cache(prepInputs, url = "https://drive.google.com/file/d/1h7gK44g64dwcoqhij24F2K54hs5e35Ci/view?usp=sharing",
                   destinationPath = Paths$inputPath,
                   fun = fixRTM, overwrite = TRUE)
studyAreaER <- Cache(prepInputs, url =  "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip",
                     destinationPath = Paths$inputPath, fun = quote(SA_ERIntersect(x = targetFilePath, studyArea)),
                     overwrite = TRUE)

speciesToUse <- c("Abie_las", "Betu_pap", "Pice_gla", "Pice_mar", "Pinu_con",
                  "Popu_tre", "Pice_eng")
speciesToUse
speciesNameConvention <- LandR::equivalentNameColumn(speciesToUse, LandR::sppEquivalencies_CA)
sppEquiv <- LandR::sppEquivalencies_CA[LandR::sppEquivalencies_CA[[speciesNameConvention]] %in% speciesToUse,]

# Assign a colour convention for graphics for each species
sppColorVect <- LandR::sppColors(sppEquiv, speciesNameConvention, palette = "Set1")
sppColorVect <- c(sppColorVect, "Mixed" = "#F781BF")


## Set module parameters -- to see options, look at module documentation
parameters = list(
  Biomass_borealDataPrep  = list(.useCache = "init",
                                 subsetDataBiomassModel = 50)
)

## Objects to provide to simInit from e.g., .GlobalEnv,
objects <- list(studyArea = studyArea, studyAreaLarge = studyArea,
                studyAreaANPP = studyAreaER,
                sppEquiv = sppEquiv,
                sppColorVect = sppColorVect)

## Set module parameters -- to see options, look at module documentation
parameters = list(
  Biomass_borealDataPrep  = list(subsetDataBiomassModel = 50),
  Biomass_core  = list(
    vegLeadingProportion = 0
  )
)


## Module inputs to load on a schedule. There may be specific inputs required by modules that are not
## supplied by default
inputs = data.frame()

## Module outputs to save on a schedule. There may be specific outputs required by modules that are not
## supplied by default
outputs = data.frame()

## See ?simInit for all options
mySim <- simInit(
  times = simTimes,
  modules = modules,
  params = parameters,
  inputs = inputs,
  objects = objects,
  outputs = outputs
)

## run the simulation
mySimOut <- spades(copy(mySim), events = "init", debug = 1) # optionally, use `spades(Copy(mySim))`
mySimOut1 <- spades(copy(mySimOut), debug = 1) # optionally, use `spades(Copy(mySim))`


