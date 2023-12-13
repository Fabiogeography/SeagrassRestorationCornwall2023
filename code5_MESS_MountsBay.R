##############################################################
##### CALCULATE MESS EXTRAPOLATION MASK FOR PLYMOUTH #########
##### Written by: Regan Early ################################
##### Written on: 8th June 2020 ##############################
##### Modified on: December 2023   ###########################
##############################################################

.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(terra)
library(raster) ## Note MESS doesn't run with terra, so must install an old version or raster instead
library(MASS)
library(dismo) ## MESS with rasters

##### Data #####
vars <- c("bathymetry", "slope", "avg_sst", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_coar_perc", "eunis_sand_perc")

wd.dat <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/combined"
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS/SPATIAL"

### Spatial model
load(paste0(wd.out,"/multivar4i_bestModel_MountsBay_rac_m2")) ## multivar4i.rac

### Environmental rasters
bathymetry <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/bathymetry/bath_cornwall_30m.tif")
slope <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/slope/slope_cornwall.tif")
avg_sst <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/sst/avg_sst_cornwall.tif")
mlw_dist <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/mlw/distance_to_mlw_cornwall.tif")
avg_spm <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/spm/spm_avg_cornwall.tif") ## is this correct?
covar_spm <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/spm/spm_covar_cornwall.tif")
expo <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/exposure/exposure_cornwall.tif")
eunis_other_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_1_pres_abs.tif")
eunis_coar_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_3_pres_abs.tif")
eunis_sand_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_4_pres_abs.tif")

## Identify the coarse variables, which need resampling to the fine variables
for (v in vars) {
  print(paste(v, res(get(v))))
}

avg_sst <- resample(avg_sst, bathymetry)
avg_spm <- resample(avg_spm, bathymetry)
covar_spm <- resample(covar_spm, bathymetry)
expo <- resample(expo, bathymetry)

env.ras <- c(bathymetry, slope, avg_sst, mlw_dist, avg_spm, covar_spm, expo,
         eunis_other_perc, eunis_coar_perc, eunis_sand_perc)
names(env.ras) <- vars

## Make all nodata values of EUNIS categories 0 (so that can remove NAs later)
env.ras$eunis_other_perc[is.na(env.ras$eunis_other_perc)] <- 0
env.ras$eunis_coar_perc[is.na(env.ras$eunis_coar_perc)] <- 0
env.ras$eunis_sand_perc[is.na(env.ras$eunis_sand_perc)] <- 0

## Create matrix of env conditions for the cells where presence/absence recorded
dat.rac <- read.csv(paste0(wd.out,"/dat_rac_MountsBay_m2.csv"))
dat.rac <- vect(dat.rac) ## Make spatial using raster
crs(dat.rac) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ## Asign the WGS84 projection
dat.rac <- terra::project(dat.rac, env.ras) ## project to match the raster

env.pa.mat <- terra::extract(env.ras, dat.rac)
env.pa.mat <- subset(env.pa.mat, select=-ID)

env.ras.r <- raster::stack(env.ras)

###### Calculate MESS #####
## MESS identifies areas where the SDM predictions are extrapolated beyond the environmental conditions used to construct the SDM, i.e. no-analogue environments
## Negative values are sites where at least one variable has a value that is outside the range of environments used in construction.
## The values in the MESS are influenced by the full distribution of the data used to construct SDMs, so that sites within the environmental range of the  construction data, but in relatively unusual environments, will have a smaller value than those in very common environments.
mess.out <- mess(x=env.ras.r, v=env.pa.mat, full=TRUE)
names(mess.out) <- c(vars, "rmess")
plot(mess.out[["rmess"]])
writeRaster(mess.out$rmess, paste0(wd.out,"/MESS_MountsBay_m2.tif"))

for(i in vars) {
  writeRaster(mess.out[[i]], paste0(wd.out,"/MESS_MountsBay_m2_",i,".tif"))
}
