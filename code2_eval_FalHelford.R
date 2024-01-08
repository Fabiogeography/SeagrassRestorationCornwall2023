##############################################################
##### EVALUATE SEAGRASS MODEL ################################
##### Written by: Regan Early ################################
##### Written on: 28th May 2020 ##############################
##### Modified on: December 2023 by Shari Mang  ##############
##############################################################

pacman::p_install(interactions, partR2) # both unavailable for this version of R

#.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(ggplot2) 
library(MuMIn)
library(lme4) ## Ensure Matrix package is up to date to support lme4 after R version 4.3.2
library(boot) ## Calculate AUC using cross-validation
library(pROC)
# library(rsq) ## perhaps not available on cluster? --> R can't find
library(cli); library(scales); library(vctrs); library(lifecycle); library(tidyverse) ## Issues witn  namespace mean I need to load some packages before tidyverse. Shouldn't be necessary normally
library(interactions) ## interact_plot
library(AICcmodavg) ## predictSE.mer
library(MLmetrics) ## Gini coefficient
library(partR2)
library(DHARMa) ## residual diagnostics for hierarchical (multi-level/mixed) regression models. https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#interpreting-residuals-and-recognizing-misspecification-problems
library(glmmTMB)
library(here)
conflicted::conflict_prefer("here", "here")


##### Data #####
vars <- c("bathymetry", "slope", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_litt_perc", "eunis_coar_perc", "eunis_sand_perc", "eunis_mix_perc") 

wd.dat <- here("variables/falhel_only/")
wd.out <- here("results/multivar_mod/")
wd.plot <- here("results/multivar_mod/plots/")

dat <- read.csv(paste0(wd.dat, "falhel_sg_env_prop_occ_unique_30m_stnd.csv")) ## Data already standardized
# only data for fal and helford
colnames(dat) # sea surface temp variables not included in this data

dat <- dat %>% ## Models will not run with NA values in data
  drop_na()
View(dat)


### Model
load(paste0(wd.out,"multivar4i_bestModel_FalHelford"))
final <- multivar4i.best

## Inspect (scaled) residuals
resids <- DHARMa::simulateResiduals(final) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.
png(paste0(wd.plot, "final_residuals.png"), width = 900, height = 700)
plot(resids)
dev.off()
## Residual metrics
## If you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. 
## That, however, doesn’t necessarily mean that you need to change your model. 
## Look at the magnitude of the residual pattern to decide if there is a problem

## Dispersion. Over/under dispersion is common in binomial models
## If overdispersion is present, the main effect is that confidence intervals tend to be too narrow, and p-values to small, leading to inflated type I error.
## The opposite is true for underdispersion, i.e. the main issue of underdispersion is that you loose power.
png(paste0(wd.plot, "final_disp.png"), width = 800, height = 800)
testDispersion(final) 
dev.off()
# Significant under dispersion (dispersion = 1.5532)
png(paste0(wd.plot, "final_zeroInf.png"), width = 800, height = 800)
testZeroInflation(final) 
dev.off()
## Zero-inflation is present.
png(paste0(wd.plot, "final_quantiles.png"), width = 800, height = 800)
testQuantiles(final) 
dev.off()

## Try refitting model with a zero-inflation term

f <- as.formula("cbind(total_pres, total_abs) ~ bathymetry + slope + mlw_dist +  
  avg_spm + covar_spm + expo + eunis_other_perc + eunis_litt_perc +  
  eunis_coar_perc + eunis_sand_perc + eunis_mix_perc + I(mlw_dist^2) +  
  bathymetry * mlw_dist + (1 | id_600m)") ## Formula used in final
final.zi <- glmmTMB(f, family=binomial, data=dat, na.action=na.fail, ziformula = ~1 )
summary(final.zi)
testZeroInflation(final.zi) ## Zero-inflation is removed

params <- cbind(fixef(final), fixef(final.zi)[[1]]) 
## Fairly similar  
  # eunis_coarse sign differs (- in model with zero-inf term)

resids.zi <- simulateResiduals(final.zi) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.
png(paste0(wd.plot, "final_zero-inf-mod_residuals.png"), width = 900, height = 700)
plot(resids.zi)
dev.off()
# less dispersion on qq plot; resid vs pred still curvey but now pulls towards 1 more
png(paste0(wd.plot, "final_zero-inf-mod_disp.png"), width = 800, height = 800)
testDispersion(final.zi) 
dev.off()
## Significant but small under-dispersion: 1.5043 (1.5532 in other model). This reduces model power, rather than suggest spurious relationships.
png(paste0(wd.plot, "final_zero-inf-mod_quantiles.png"), width = 800, height = 800)
testQuantiles(final.zi) 
dev.off()
## Sig. p-value means there is heteroscedasticity and there is an upward trend in the quantile splines that wasn't present in the model without dispersion




## Plot residuals against all predictors to see what could be causing the trend in variance
## If the predictor is a factor, or if there is just a small number of observations on the x axis, plotResiduals will plot a box plot with additional tests instead of a scatter plot.
for(v in vars) {
  x <- dat[,v]
  par(mfrow=c(1,2))
  plotResiduals(resids, form=x, xlab=v)
  plotResiduals(resids.zi, form=x, xlab=v)
}
## When zero-inflation is in the model, there is slightly more trend in residuals for some explanatory variables, but also less significant quantile splines for others. 

#### Later investigations found that incorporating a spatial autocovariate in the model removed most of the zero-inflation.
### Given the expected importance of zero inflation in the dataset, it was demeed preferable to use the spatial autocovariate than pursue the zero-inflated model.

## Spatial autocorrelation
testSpatialAutocorrelation(resids, x=dat$lon, y=dat$lat) ## p-value significant so spatial autocorrelation exists
testSpatialAutocorrelation(resids.zi, x=dat$lon, y=dat$lat) ## p-value significant so spatial autocorrelation exists

##### Plot explanatory relationships #####
### Calculate the means and standard deviations of the unstandardised data. Don't need to run if run once already, just load csv
dat.unstd <- read.csv(paste0(wd.dat, "falhel_sg_env_prop_occ_unique_30m.csv")) ## Shari's data_prep code tells me this is the unstandardised data used to make the standardised data used in modelling

means <- apply(dat.unstd[,vars], 2, function(x) mean(x, na.rm=T))
sds <- apply(dat.unstd[,vars], 2, function(x) sd(x, na.rm=T))

destdize <- rbind(means, sds)
write.csv(destdize, paste0(wd.dat, "/allAreas_means_sds.csv"), row.names=F)
destdize <- read.csv(paste0(wd.dat, "/allAreas_means_sds.csv"))


##### Plot graphs #####

### Make data frame of explanatory variables on which to make predictions
vars <- names(fixef(final))
vars <- vars[vars!="(Intercept)"]
vars <- vars[!vars %in% grep("I", vars, value = T)] ## linear terms only
vars <- vars[!vars %in% grep(":", vars, value = T)] ## main effects terms only

mns <- apply(dat[,vars], 2, mean)

## y.limits for graphs
ylims <- rep(1, length(vars))
names(ylims) <- vars

##### Final model #####
jpeg(paste0(wd.out,"/MULTIVAR4i_results_FalHelford_1.jpg"), width=600, height=400)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1)) # c(bottom, left, top, right)
for (i in vars[1:4]) { ## Run through vars in groups of fours. Rename jpeg with appropriate number
  newdat <- data.frame(matrix(data=mns, byrow=T, nrow=10000, ncol=length(vars), dimnames=(list(NULL,vars))))
  newdat[,i] <- seq(min(dat[,i]), max(dat[,i]), length.out=10000)
  newdat$id_600m <- rep_len(dat$id_600m, nrow(newdat))
  
 p <- as.data.frame(predictSE.mer(mod=final, newdata=newdat, se.fit = TRUE,
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
par(mfrow=c(1,1))
##### Interaction plots ##### 
### bathymetry:covar_spm interaction
jpeg(paste0(wd.out,"/MULTIVAR4i_results_FalHelford_interact_bathyCovarSPM.jpg"), width=600, height=400)
i <- interactions::interact_plot(final, pred=bathymetry, modx=covar_spm, interval=T,
                            y.label="Prob. Occupied", legend.main="Coef. Var. SPM")#, rug=T)
dev.off()