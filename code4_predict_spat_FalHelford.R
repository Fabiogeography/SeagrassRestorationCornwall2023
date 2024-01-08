##############################################################
##### PREDICT SEAGRASS SUITABILITY FOR RESTORATION ###########
##### Written by: Regan Early ################################
##### Written on: 21st November 2023 ##############################
##### Modified on: 21st December 2023 by Shari Mang ###############
##############################################################


install.packages("merTools", dependencies=TRUE, repos='http://cran.rstudio.com/')
library(merTools) ## predictInterval - bootstrapped confidence intervals for mixed models in a reasonable time. averageObs.
library(terra)
library(here)
conflicted::conflict_prefer("here", "here")

##### Data #####
vars <- c("bathymetry", "slope", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_litt_perc", "eunis_coar_perc", "eunis_sand_perc", "eunis_mix_perc") 
wd.dat <- here("variables/falhel_only/")
wd.out <- here("results/multivar_mod/spatial/")

### Read in the standardised data with spatial autocovariate
dat.rac <- read.csv(paste0(wd.out,"dat_rac.csv"))
destdize <- read.csv(paste0(wd.dat, "/allAreas_means_sds.csv"))

### Spatial model
load(paste0(wd.out,"multivar4i_bestModel_FalHelford_rac")) # multivar4i.rac


### Environmental rasters
bathymetry <- terra::rast(here("variables/bathymetry/bath_cornwall_30m.tif"))
slope <- terra::project(terra::rast(here("variables/slope/slope_cornwall.tif")), crs(bathymetry)) ## only one not in BNG, not sure why
mlw_dist <- terra::rast(here("variables/mlw/distance_to_mlw_cornwall.tif"))
avg_spm <- terra::rast(here("variables/spm/spm_avg_cornwall.tif"))
covar_spm <- terra::rast(here("variables/spm/spm_covar_cornwall.tif"))
expo <- terra::rast(here("variables/exposure/exposure_cornwall.tif"))
eunis_other_perc <- terra::rast(here("variables/eunis/study_area/cornwall_eunis_cat_1_pres_abs.tif"))
eunis_litt_perc <- terra::rast(here("variables/eunis/study_area/cornwall_eunis_cat_2_pres_abs.tif"))
eunis_coar_perc <- terra::rast(here("variables/eunis/study_area/cornwall_eunis_cat_3_pres_abs.tif"))
eunis_sand_perc <- terra::rast(here("variables/eunis/study_area/cornwall_eunis_cat_4_pres_abs.tif"))
eunis_mix_perc <- terra::rast(here("variables/eunis/study_area/cornwall_eunis_cat_6_pres_abs.tif"))


## Identify the coarse variables, which need resampling to the fine variables
for (v in vars) {
  print(paste(v, res(get(v))))
}
avg_spm <- terra::resample(avg_spm, bathymetry)
covar_spm <- terra::resample(covar_spm, bathymetry)
expo <- terra::resample(expo, bathymetry)

env.ras <- c(bathymetry, slope, mlw_dist, avg_spm, covar_spm, expo,
         eunis_other_perc, eunis_litt_perc, eunis_coar_perc, eunis_sand_perc, eunis_mix_perc) 
names(env.ras) <- vars


##### Predict (can't do the inverse link process or the GLMM process direct to raster) #####
### Put environmental data through same standardisation process as the original environmental data
for (n in names(env.ras)) {
  std <- destdize[,n] ## values of model input environment data to standardise raster values
  env.ras[[n]] <- (env.ras[[n]] - std[1]) / std[2] ## subtract original mean and divide by original standard deviation
}

## Add autocovariate term
rac <- terra::rast(paste0(wd.out,"rac_FalHelford.tif"))
names(rac) <- "rac"
env.ras <- c(env.ras, rac)

## Save standardised and rac rasters together
writeRaster(env.ras, paste0(wd.dat,"env_rac_standardised.tif")) ## Saves all layers

## Make all nodata values of EUNIS categories 0 (so that can remove NAs later)
env.ras$eunis_other_perc[is.na(env.ras$eunis_other_perc)] <- 0
env.ras$eunis_litt_perc[is.na(env.ras$eunis_litt_perc)] <- 0
env.ras$eunis_coar_perc[is.na(env.ras$eunis_coar_perc)] <- 0
env.ras$eunis_sand_perc[is.na(env.ras$eunis_sand_perc)] <- 0
env.ras$eunis_mix_perc[is.na(env.ras$eunis_mix_perc)] <- 0

## Extract data from raster to data.frame
env <- as.data.frame(terra::extract(env.ras, 1:ncell(env.ras)))

## Add the random effect value
mean(multivar4i.rac@frame$id_600m) # 3417.318  # @frame is the data associated with a merMod object
env$id_600m <- 3417  



### Make the prediction for the final model, including 95% confidence intervals
p <- predictInterval(multivar4i.rac, env, type="probability", which="fixed", level=0.95, n.sims=100, include.resid.var=F, stat="median")

### Put the predictions in a raster
## Mean prediction
p.ras <- env.ras[[1]]
names(p.ras) <- "prediction"
values(p.ras) <- p$fit
writeRaster(p.ras, file=paste0(wd.out,"/pred_rac_FalHelford.tiff"))

## Lower 95% prediction
p.l95 <- env.ras[[1]]
names(p.l95) <- "l95"
values(p.l95) <- p$lwr
writeRaster(p.l95, file=paste0(wd.out,"/pred_rac_FalHelford_l95.tiff"))

## Upper 95% prediction
p.u95 <- env.ras[[1]]
names(p.u95) <- "u95"
values(p.u95) <- p$upr
writeRaster(p.u95, file=paste0(wd.out,"/pred_rac_FalHelford_u95.tiff"))
##############
