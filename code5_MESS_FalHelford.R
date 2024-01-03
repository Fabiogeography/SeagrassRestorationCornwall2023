##############################################################
##### CALCULATE MESS EXTRAPOLATION MASK FOR PLYMOUTH #########
##### Written by: Regan Early ################################
##### Written on: 8th June 2020 ##############################
##### Modified on: December 2023 by Shari Mang ###############
##############################################################

#.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(terra)
library(raster) ## Note MESS doesn't run with terra, so must install an old version or raster instead
library(MASS)
library(dismo) ## MESS with rasters
library(here)
conflicted::conflict_prefer("here", "here")

##### Data #####
vars <- c("bathymetry", "slope", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_litt_perc", "eunis_coar_perc", "eunis_sand_perc", "eunis_mix_perc") 
wd.dat <- here("variables/falhel_only/")
wd.out <- here("results/multivar_mod/spatial/")
bng_epsg <- 27700

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


## Make all nodata values of EUNIS categories 0 (so that can remove NAs later)
env.ras$eunis_other_perc[is.na(env.ras$eunis_other_perc)] <- 0
env.ras$eunis_litt_perc[is.na(env.ras$eunis_litt_perc)] <- 0
env.ras$eunis_coar_perc[is.na(env.ras$eunis_coar_perc)] <- 0
env.ras$eunis_sand_perc[is.na(env.ras$eunis_sand_perc)] <- 0
env.ras$eunis_mix_perc[is.na(env.ras$eunis_mix_perc)] <- 0


## Create matrix of env conditions for the cells where presence/absence recorded
dat.rac <- read.csv(paste0(wd.out,"dat_rac.csv"))
projcrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" # Assign the WGS84 projection
dat.rac <- vect(dat.rac, crs = projcrs) 
dat.rac <- terra::project(dat.rac, env.ras) ## project to match the raster

env.pa.mat <- terra::extract(env.ras, dat.rac)
env.pa.mat <- subset(env.pa.mat, select=-ID)

# convert to raster object 
env.ras.r <- raster::stack(env.ras)


###### Calculate MESS #####
## MESS identifies areas where the SDM predictions are extrapolated beyond the environmental conditions used to construct the SDM, i.e. no-analogue environments
## Negative values are sites where at least one variable has a value that is outside the range of environments used in construction.
## The values in the MESS are influenced by the full distribution of the data used to construct SDMs, so that sites within the environmental range of the  construction data, but in relatively unusual environments, will have a smaller value than those in very common environments.
mess.out <- mess(x=env.ras.r, v=env.pa.mat, full=TRUE)
names(mess.out) <- c(vars, "rmess")
plot(mess.out[["rmess"]])
writeRaster(mess.out$rmess, paste0(wd.out,"MESS_FalHelford.tif"), overwrite = TRUE)

for(i in vars) {
  writeRaster(mess.out[[i]], paste0(wd.out,"MESS_FalHelford_",i,".tif"))
}

# Mask of negative values --> indicates dissimilar environments 
mess_mask <- mess.out$rmess < 0
plot(mess_mask)
writeRaster(mess_mask, paste0(wd.out,"MESS_mask_negval.tif"), overwrite = TRUE)


# Load alt mess rmess that doesn't have values on land
mess.alt <- terra::rast(paste0(wd.out,"MESS_FalHelford_alt.tif"))
mess_positive <- mess.alt
mess_positive[mess_positive < 0] <- NA
mess_positive <- subst(mess_positive, 0, NA) # don't care about zeros either
plot(mess_positive)
writeRaster(mess_positive, paste0(wd.out, "MESS_FalHelford_positive.tif"), overwrite = TRUE)


