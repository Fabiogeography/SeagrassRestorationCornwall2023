##############################################################
##### INCLUDE SPATIAL AUTOCOVARIATE IN MULTIVARATE SEAGRASS MODEL ################################
##### Written by: Regan Early ################################
##### Written on: 21st Nov 2023 ##############################
##### Modified on: December 2023   ###########################
##############################################################

.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(Rcpp); library(terra) ## raster operations, focal
library(car) ## Variance Inflation Factor, vif
library(lme4)
library(tidyverse)
library(AICcmodavg) ## predictSE.mer
library(DHARMa) ## residual diagnostics for hierarchical (multi-level/mixed) regression models. https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#interpreting-residuals-and-recognizing-misspecification-problems
library(partR2)
library(pROC) ## rocr
library(interactions)

##### Data #####
vars <- c("bathymetry", "slope", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_litt_perc", "eunis_coar_perc", "eunis_sand_perc", "eunis_mix_perc") 

wd.dat <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/combined"
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS"

### Read in the standardised data
dat <- read.csv(paste0(wd.dat,"/cornwall_sg_env_unique_30m_stnd.csv")) ## Data already standardized
dat <- dat[dat$site=="FalHelford",]

dat <- dat %>% ## Models will not run with NA values in data
  drop_na()

### Non-spatial model
load(paste0(wd.out,"/multivar4i_bestModel_FalHelford"))
final <- multivar4i.best

### Environmental rasters
bathymetry <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/bathymetry/bath_cornwall_30m.tif")
slope <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/slope/slope_cornwall.tif")
max_sst <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/sst/max_sst_cornwall.tif")
avg_sst <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/sst/avg_sst_cornwall.tif")
mlw_dist <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/mlw/distance_to_mlw_cornwall.tif")
avg_spm <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/spm/spm_avg_cornwall.tif") ## is this correct?
covar_spm <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/spm/spm_covar_cornwall.tif")
expo <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/exposure/exposure_cornwall.tif")
eunis_other_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_1_pres_abs.tif")
eunis_litt_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_2_pres_abs.tif")
eunis_coar_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_3_pres_abs.tif")
eunis_sand_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_4_pres_abs.tif")
eunis_mix_perc <- rast("F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/eunis/study_area/cornwall_eunis_cat_6_pres_abs.tif")

## Identify the coarse variables, which need resampling to the fine variables
for (v in vars) {
  print(paste(v, res(get(v))))
}

max_sst <- resample(max_sst, bathymetry)
avg_sst <- resample(avg_sst, bathymetry)
avg_spm <- resample(avg_spm, bathymetry)
covar_spm <- resample(covar_spm, bathymetry)
expo <- resample(expo, bathymetry)

env <- c(bathymetry, slope, mlw_dist, avg_spm, covar_spm, expo,
         eunis_other_perc, eunis_litt_perc, eunis_coar_perc, eunis_sand_perc, eunis_mix_perc) 
names(env) <- vars

##### Extract and map spatial autocorrelation in the residuals - GLM #####
## Set up a blank rasterfile
rast <- (env$slope / env$slope) - 1

## Extract residuals from the GLM
resids <- residuals(final)

##### Alternative way to calculate residuals #####
resids.sim <- simulateResiduals(final) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.
summary(resids.sim$scaledResiduals)
plot(resids, resids.sim[["scaledResiduals"]])

xy_resids <- as.data.frame(cbind(dat$lon, dat$lat, resids))
colnames(xy_resids) <- c("lon","lat", "resids")

xy_resids <- vect(xy_resids) ## Make spatial
crs(xy_resids) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" ## Asign the WGS84 projection
xy_resids <- project(xy_resids, rast) ## project to match the raster

plot(rast); points(xy_resids) ## Seems right

## Map the residuals
# rast[xy_resids] <- xy_resids$resids ## Fails in R 4.3.2, possibly due to updated terra package
cellnum <- terra::cellFromXY(rast, geom(xy_resids)[,c("x","y")])
rast[cellnum] <- xy_resids$resids
names(rast) <- "resids"
summary(rast, size=ncell(rast)) ## Ensure this matches the range of values in xy_resids

writeRaster(rast, paste0(wd.out,"/SPATIAL/resids_FalHelford.tif"))

#####  Use residual autocovariateas a predictor #####
## Focal operations: ngb is neighbourhood size, set to 3 by 3 cells; fun is function, 
## here the mean value within the defined neighbourhood

## Calculate the moving window values for the neighbourhood of focal cells. 
rac <- focal(rast, w=matrix(1/9,nrow=3,ncol=3), fun="sum") ## 3x3 mean filter. Something weird happening with this on 1st Dec 2023 (agter R upgrade)- use the saved raster instead
plot(rac)
plot(xy_resids, add=T)
writeRaster(rac, paste0(wd.out,"/SPATIAL/rac_FalHelford.tif"))

## Extract the autocovariate to a vector and add to the presence/absence locations
rac <- terra::extract(rac, xy_resids)$resids
rac[is.na(rac)] <- 0 ## remove NAs
dat.rac <- cbind(dat, rac)
write.csv(dat.rac, paste0(wd.out,"/SPATIAL/dat_rac.csv"), row.names=F)

## fit the RAC model using the original environmental variables and the residuals autocovariate
multivar4i.rac <- update(final, .~. + rac , data=dat.rac) ## failed to converge
AIC(final, multivar4i.rac) ## Model with spatial covariate is better, delta AIC -2426.
vif(multivar4i.rac) ## same as previously - sand and mix are a bit dubious. Bathymetry VIF is high because of colinearity with the quadratic term and interaction - not of concern. 

save(multivar4i.rac, file=paste0(wd.out,"/SPATIAL/multivar4i_bestModel_FalHelford_rac")) 

## Remove sea surface temperature variables as interested in extrapolating model beyond SST raster.
multivar4i.rac.noSST <- update(final, .~. + rac -max_sst -avg_sst, data=dat.rac) ## failed to converge
AIC(final, multivar4i.rac.noSST) ## Model with spatial covariate is better, delta AIC -2426.
vif(multivar4i.rac.noSST) ## same as previously - sand and mix are a bit dubious

save(multivar4i.rac.noSST, file=paste0(wd.out,"/SPATIAL/multivar4i_bestModel_FalHelford_rac_noSST")) 

### multivar4i.rac did not converge. However, if all optimizers converge to values that are practically equivalent, then the model fit is good enough 
multivar4i.rac.all <- allFit(multivar4i.rac)
ss <- summary(multivar4i.rac.all)
ss$ fixef               ## fixed effects
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors of the random effect covariance matrix
ss$ which.OK            ## which fits worked
### Values are all very similar so ignore the convergence warnings for now at least

##### How does the model with teh spatial autocovariate perform? #####
## Compare parameter estimates with and without spatial autocovariate ***and SST***
params <- cbind(c(fixef(final),NA), fixef(multivar4i.rac)) ## avg_spm, eunis_other, eunis_sand, have steeper effects. This is driven by removal of SST, rather than the autocovariate.

write.csv(as.data.frame(fixef(multivar4i.rac)), paste0(wd.out, "/SPATIAL/params_FalHelford.csv"), row.names=T)

## Inspect (scaled) residuals
resids.rac <- simulateResiduals(multivar4i.rac) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.
plot(resids.rac) ## Looks pretty good

## Residual metrics
##  If you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesn’t necessarily mean that you need to change your model. 
## Look at the magnitude of the residual pattern to decide if the diagnostic is of concern.

## Dispersion. Over/under dispersion is common in binomial models
## If overdispersion is present, the main effect is that confidence intervals tend to be too narrow, and p-values to small, leading to inflated type I error.
## The opposite is true for underdispersion, i.e. the main issue of underdispersion is that you loose power.
par(mfrow=c(1,2))
testDispersion(multivar4i.rac) ## Resids are not significantly over or under--dispersed. 
testZeroInflation(multivar4i.rac) ## Only very slight zero-inflation (1.4, 1.3 without SST).

## Plot residuals against all predictors to see what could be causing the trend in variance
## If the predictor is a factor, or if there is just a small number of observations on the x axis, plotResiduals will plot a box plot with additional tests instead of a scatter plot.
resids.orig <- simulateResiduals(final) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.

par(mfrow=c(1,2))
for(v in vars) {
  x <- dat[,v]
  par(mfrow=c(1,2))
  plotResiduals(resids.orig, form=x, xlab=v)
  plotResiduals(resids.rac, form=x, xlab=v)
}
## When autocovariate is in the model, there is no more trend in residuals for explanatory variables than in the original model

## Spatial autocorrelation
par(mfrow=c(1,2))
testSpatialAutocorrelation(resids.rac, x=dat$lon, y=dat$lat) ## p-value significant so spatial autocorrelation remains

##### Plot explanatory relationships #####
### Calculate the means and standard deviations of the unstandardised data. Don't need to run if run once already, just load csv
# dat.unstd <- read.csv(paste0(wd.dat,"/cornwall_sg_env_prop_presabs_unique_30m.csv")) ## Shari's data_prep code tells me this is the unstandardised data used to make the standardised data used in modelling
# 
# means <- apply(dat.unstd[,vars], 2, function(x) mean(x, na.rm=T))
# sds <- apply(dat.unstd[,vars], 2, function(x) sd(x, na.rm=T))
# 
# destdize <- rbind(means, sds)
# write.csv(destdize, paste0(wd.dat, "/allAreas_means_sds.csv"), row.names=F)
destdize <- read.csv(paste0(wd.dat, "/allAreas_means_sds.csv"))

##### Plot graphs #####

### Make data frame of explanatory variables on which to make predictions
vars <- names(fixef(multivar4i.rac))
vars <- vars[vars!="(Intercept)"]
vars <- vars[!vars %in% grep("I", vars, value = T)] ## linear terms only
vars <- vars[!vars %in% grep(":", vars, value = T)] ## main effects terms only

mns <- apply(dat.rac[,vars], 2, mean)

## y.limits for graphs
ylims <- rep(1, length(vars))
names(ylims) <- vars

##### Final model #####
jpeg(paste0(wd.out,"/SPATIAL/MULTIVAR4i_results_FalHelford_1.jpg"), width=600, height=400)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1)) # c(bottom, left, top, right)
for (i in vars[1:4]) { ## Run through vars in groups of fours. Rename jpeg with appropriate number
  newdat <- data.frame(matrix(data=mns, byrow=T, nrow=10000, ncol=length(vars), dimnames=(list(NULL,vars))))
  newdat[,i] <- seq(min(dat.rac[,i]), max(dat.rac[,i]), length.out=10000)
  newdat$id_600m <- rep_len(dat.rac$id_600m, nrow(newdat))
  
  p <- as.data.frame(predictSE.mer(mod=multivar4i.rac, newdata=newdat, se.fit = TRUE,
                type="response", level = 0, print.matrix = TRUE)) ## Based on fixed effects only. The current version of the function only supports predictions for the populations excluding random effects (i.e., level = 0).
  
  
  ### Calculate the raw (destandardised) values of the explanatory variable
  ## Multiply by standard deviation (row 2 of destdize) and add mean (row 1 of destdize)
  x.std <- seq(min(newdat[,i]), max(newdat[,i]), length.out=6)
  x.destd <- signif((x.std * destdize[2,i]) + destdize[1,i], digits=4) ## Use signif rather than round as it uses significant figures rather than decimal places
  x.destd[abs(x.destd)<0.001] <- 0 ## The significant figures function gives some tiny values for the percentage cover variables - convert to 0
  
  plot(newdat[,i], p$fit, type="n", xlab=i, ylab="Suitability", xaxt="n", yaxt="n", ylim=c(0,ylims[i])) ## Note that the input data is a proportion, not a presence/absence, so 'suitability' is a better name for the y axis than Prob. Occupied.
  axis(side=1, at=x.std, labels=x.destd) ## Add x axis on original scale
  axis(side=2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1))
  l.ci <- p$fit - (1.96*p$se.fit)
  u.ci <- p$fit + (1.96*p$se.fit)
  polygon(c(newdat[,i], rev(newdat[,i])), c(c(l.ci, rev(u.ci))), col = 'grey80', border = NA)
  lines(newdat[,i], p$fit)

}

dev.off()

##### Interaction plots ##### 
### bathymetry:covar_spm interaction
jpeg(paste0(wd.out,"/SPATIAL/MULTIVAR4i_results_FalHelford_interact_bathyCovarSPM_noSST.jpg"), width=600, height=400)
interactions::interact_plot(multivar4i.rac, pred=bathymetry, modx=covar_spm, interval=T,
                            y.label="Prob. Occupied", legend.main="Coef. Var. SPM")#, rug=T)
dev.off()

##### Calculate AUC by cross-validation #####
AUC <- vector()

a <- names(fixef(multivar4i.rac))
a <- a[a!="(Intercept)"]
a <- paste(a, collapse=" + ")
f <- as.formula(paste0("cbind(total_pres, total_abs) ~ ", a, " + (1|id_600m)"))

for(i in 1:10) { ## 10-fold cross-validation
  n.calib <- trunc(nrow(dat.rac)*0.7) ## Number of data rows to use for calibration
  dat.calib <- dat.rac[sample(nrow(dat.rac), n.calib), ] ## Randomly selected calibration data
  dat.valid <- dat.rac[!(rownames(dat.rac) %in% rownames(dat.calib)),] ## Validation data, independent of calibration data
  
  m.calib <- glmer(f, family=binomial, dat.calib, na.action=na.fail) ## Model failed to converge
  p.valid <- predict(m.calib, dat.valid, re.form=NA, type="response")
  
  roc_obj <- roc(dat.valid$present, p.valid)
  AUC <- c(AUC, auc(roc_obj))
}
mean(AUC) ## 
sd(AUC) ## 

##### Calculate partial R-squared of each variable in the final model #####
### https://cran.r-project.org/web/packages/partR2/vignettes/Using_partR2.html
Rsquare <- partR2(multivar4i.rac, partvars=vars, max_level=1, # the argument max_level=1 calculates part R2 for each predictor, but not for all their combinations. We can specify this with .
                 R2_type = "marginal", nboot = 10) ## Marginal R2 refers to the variance explained by fixed effect predictors relative to the total variance in the respons

Rsquare$R2   # Partial R2s Very low when random intercept of 600m grid-cell not included
write.csv(Rsquare$R2, paste0(wd.out,"/SPATIAL/Rsq_FalHelford.csv"), row.names=F)

library(patchwork)
p1 <- forestplot(Rsquare, type = "R2", text_size = 10)
p2 <- forestplot(Rsquare, type = "IR2", text_size = 10) ## Inclusive R2. This is SC^2 * R2_full.
p3 <- forestplot(Rsquare, type = "SC", text_size = 10) ## Structure coefficients are the correlation between a predictor and the predicted response. Squared structure coefficients indicate how much of the regression effect can be attributed to a given predictor
p4 <- forestplot(Rsquare, type = "BW", text_size = 10) ## Standardised model estimates (beta weights) for fixed effects. Beta weights for Gaussian models are calculated as beta * sd(x)/sd(y), with beta being the estimated slope of a fixed effect for predictor x and response y. Beta weight for Non-Gaussian models are calculated as beta * sd(x). Beta weights for interactions or polynomial terms are not informative at the moment and we recommend users to standardise variables themselves before fitting the model and to look at the model estimates (Ests) instead of beta weights (BW) in the partR2 output. See vignette for details.
(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")
## SCs are better than beta weights when there is multicollinearity in the explanatory variables. A low β weight and a high structure coefficient indicates multicollinearity, and beta weight is masked by multicollinearity
## Confidence Intervals are 95% by default
