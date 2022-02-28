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

## Module documentation -- please go to these pages to read about each module
##  In some cases, there will be important defaults that a user should be aware of
##  or important objects (like studyArea) that may be essential
## browseURL('https://github.com/PredictiveEcology/Biomass_borealDataPrep/blob/main/Biomass_borealDataPrep.rmd')
## browseURL('https://github.com/PredictiveEcology/Biomass_core/blob/main/Biomass_core.rmd')


## simulation initialization. These may not be appropriate start and end times for one
## or more of the modules, e.g., they may only be defined with calendar dates
simTimes <- list(start = 0, end = 10)

## Set module parameters -- to see options, look at module documentation
parameters = list(
  Biomass_borealDataPrep  = list(),
  Biomass_core  = list()
)

## Objects to provide to simInit from e.g., .GlobalEnv,
objects = list()

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
mySimOut <- spades(mySim, debug = 1) # optionally, use `spades(Copy(mySim))`


