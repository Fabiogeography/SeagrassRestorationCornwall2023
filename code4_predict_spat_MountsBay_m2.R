##############################################################
##### PREDICT SEAGRASS SUITABILITY FOR RESTORATION ###########
##### Written by: Regan Early ################################
##### Written on: 21st November 2023 ##############################
##### Modified on: December 2023  ###########################
##############################################################

.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(merTools) ## predictInterval - bootstrapped confidence intervals for mixed models in a reasonable time. averageObs.
library(terra)

##### Data #####
vars <- c("bathymetry", "slope", "avg_sst", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_coar_perc", "eunis_sand_perc")

wd.dat <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/combined"
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS/SPATIAL"

### Read in the standardised data with spatial autocovariate
dat.rac <- read.csv(paste0(wd.out,"/dat_rac_MountsBay_m2.csv"))

destdize <- read.csv(paste0(wd.dat, "/allAreas_means_sds.csv"))

### Spatial model
load(paste0(wd.out,"/multivar4i_bestModel_MountsBay_rac_m2"))
final <- multivar4i.rac

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

##### Predict (can't do the inverse link process or the GLMM process direct to raster) #####
### Put environmental data through same standardisation process as the original environmental data
for (n in names(env.ras)) {
  std <- destdize[,n] ## values of model input environment data to standardise raster values
  env.ras[[n]] <- (env.ras[[n]] - std[1]) / std[2] ## subtract original mean and divide by original standard deviation
}

## Add autocovariate term
rac <- rast(paste0(wd.out,"/rac_MountsBay_m2.tif"))
names(rac) <- "rac"
env.ras <- c(env.ras, rac)

## Save standardised and rac rasters together
writeRaster(env.ras, paste0(wd.dat,"/env_standardised_MountsBay_m2.tif")) ## Saves all layers

## Make all nodata values of EUNIS categories 0 (so that can remove NAs later)
env.ras$eunis_other_perc[is.na(env.ras$eunis_other_perc)] <- 0
env.ras$eunis_coar_perc[is.na(env.ras$eunis_coar_perc)] <- 0
env.ras$eunis_sand_perc[is.na(env.ras$eunis_sand_perc)] <- 0

## Extract data from raster to data.frame
env <- as.data.frame(terra::extract(env.ras, 1:ncell(env.ras)))

## Add the random effect value
env$id_600m <- 3334 # Derived from mean(final@frame$id_600m) ## @frame is the data associated with a merMod object

### Make the prediction for the final model, including 95% confidence intervals
p <- predictInterval(final, env, type="probability", which="fixed", level=0.95, n.sims=100, include.resid.var=F, stat="median")

### Put the predictions in a raster
## Mean prediction
p.ras <- env.ras[[1]]
names(p.ras) <- "prediction"
values(p.ras) <- p$fit
writeRaster(p.ras, file=paste0(wd.out,"/pred_rac_MountsBay_m2.tiff"))

## Lower 95% prediction
p.l95 <- env.ras[[1]]
names(p.l95) <- "l95"
values(p.l95) <- p$lwr
writeRaster(p.l95, file=paste0(wd.out,"/pred_rac_MountsBay_m2_l95.tiff"))

## Upper 95% prediction
p.u95 <- env.ras[[1]]
names(p.u95) <- "u95"
values(p.u95) <- p$upr
writeRaster(p.u95, file=paste0(wd.out,"/pred_rac_MountsBay_m2_u95.tiff"))
##############