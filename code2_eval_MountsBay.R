##############################################################
##### EVALUATE SEAGRASS MODEL ################################
##### Written by: Regan Early ################################
##### Written on: 28th May 2020 ##############################
##### Modified on: December 2023   ###########################
##############################################################

.libPaths("C:/SOFTWARE/R/R-4.3.2/library")
library(cli); library(scales); library(vctrs); library(lifecycle); library(tidyverse) ## Issues witn  namespace mean I need to load some packages before tidyverse. Shouldn't be necessary normally
library(ggplot2) 
library(MuMIn)
library(lme4) ## Ensure Matrix package is up to date to support lme4 after R version 4.3.2
library(boot) ## Calculate AUC using cross-validation
library(pROC)
library(rsq) ## perhaps not available on cluster?
library(interactions) ## interact_plot
library(AICcmodavg) ## predictSE.mer
library(MLmetrics) ## Gini coefficient
library(partR2)
library(DHARMa) ## residual diagnostics for hierarchical (multi-level/mixed) regression models. https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html#interpreting-residuals-and-recognizing-misspecification-problems

##### Data #####
vars <- c("bathymetry", "slope", "avg_sst", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_coar_perc", "eunis_sand_perc")

wd.dat <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/combined"
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS"

dat <- read.csv(paste0(wd.dat,"/cornwall_sg_env_unique_30m_stnd.csv")) ## Data already standardized
dat <- dat[dat$site=="MountsBay",]

dat <- dat %>% ## Models will not run with NA values in data
  drop_na()

### Model
load(paste0(wd.out,"/multivar4i_bestModel_MountsBay"))
load(paste0(wd.out,"/multivar4i_2ndbestModel_MountsBay"))

## Inspect (scaled) residuals
resids.best <- simulateResiduals(multivar4i.best) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.
plot(resids.best) ## Look really good

resids.2ndbest <- simulateResiduals(multivar4i.2ndbest) ## calculate scaled residuals. a scaled residual value of 0.5 means that half of the simulated data are higher than the observed value, and half of them lower. Min / max values are 0 / 1.
plot(resids.2ndbest) ## Look amazing

## Residual metrics
##  If you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesn’t necessarily mean that you need to change your model. 
## Look at the magnitude of the residual pattern

## Dispersion. Over/under dispersion is common in binomial models
## If overdispersion is present, the main effect is that confidence intervals tend to be too narrow, and p-values to small, leading to inflated type I error.
## The opposite is true for underdispersion, i.e. the main issue of underdispersion is that you loose power.
testDispersion(resids.best) ## Resids are not significantly over or under--dispersed. 
testZeroInflation(resids.best) ## No Zero-inflation is present. Don't worry!

testDispersion(resids.2ndbest) ## Resids are not significantly over or under--dispersed. 
testZeroInflation(resids.2ndbest) ## No Zero-inflation is present. Don't worry!

## Spatial autocorrelation
testSpatialAutocorrelation(resids.best, x=dat$lon, y=dat$lat) ## p-value significant so spatial autocorrelation exists
testSpatialAutocorrelation(resids.2ndbest, x=dat$lon, y=dat$lat) ## p-value significant so spatial autocorrelation exists

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
mns <- apply(dat[,vars], 2, mean)

## y.limits for graphs
ylims <- rep(1, length(vars))
names(ylims) <- vars

##### Final model #####
jpeg(paste0(wd.out,"/MULTIVAR4i_results_MountsBay_model2_1.jpg"), width=600, height=400)
par(mfrow=c(2,2))
par(mar=c(4,4,2,1)) # c(bottom, left, top, right)
for (i in vars[1:4]) { ## Run through vars in groups of fours. Rename jpeg with appropriate number
  newdat <- data.frame(matrix(data=mns, byrow=T, nrow=10000, ncol=length(vars), dimnames=(list(NULL,vars))))
  newdat[,i] <- seq(min(dat[,i]), max(dat[,i]), length.out=10000)
  newdat$id_600m <- rep_len(dat$id_600m, nrow(newdat))
  
  p <- as.data.frame(predictSE.mer(mod=multivar4i.2ndbest, newdata=newdat, se.fit = TRUE,
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
### expo:eunis_other_perc interaction in best model
jpeg(paste0(wd.out,"/MULTIVAR4i_results_MountsBay_interact_expoEunisOther_m1.jpg"), width=600, height=400)
interactions::interact_plot(multivar4i.best, pred=expo, modx=eunis_other_perc, interval=T,
                                 y.label="Prob. Occupied", legend.main="Eunis Other")#, rug=T)
dev.off()

### bathymetry:avg_spm interaction
jpeg(paste0(wd.out,"/MULTIVAR4i_results_MountsBay_interact_bathyAvgSPM_m2.jpg"), width=600, height=400)
interactions::interact_plot(multivar4i.2ndbest, pred=bathymetry, modx=avg_spm, interval=T,
                            y.label="Prob. Occupied", legend.main="Avg. SPM")#, rug=T)
dev.off()


##### Calculate Gini's coefficient by cross-validation #####
## Gini coefficient can compare two continuous predictions: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.stonybrook.edu/commcms/irpe/reports/_presentations/DataMiningDemystified_Galambos_2016_06_02.pdf
## Also called Somer's D

### Best model
gini <- vector()

a <- names(fixef(multivar4i.best))
a <- a[a!="(Intercept)"]
a <- paste(a, collapse=" + ")
f <- as.formula(paste0("present ~ ", a, " + (1|id_600m)"))

for(i in 1:10) { ## 10-fold cross-validation
  n.calib <- trunc(nrow(dat)*0.7) ## Number of data rows to use for calibration
  dat.calib <- dat[sample(nrow(dat), n.calib), ] ## Randomly selected calibration data
  dat.valid <- dat[!(rownames(dat) %in% rownames(dat.calib)),] ## Validation data, independent of calibration data
  
  m.calib <- lme4::glmer(f, family=binomial, dat.calib, na.action=na.fail) ## Model failed to converge
  p.valid <- predict(m.calib, dat.valid, re.form=NA, type="response")
  
  g <- Gini(p.valid, dat.valid$prop_pres)
  gini <- c(gini, g)
}
mean(gini) ## 0.7242677
sd(gini) ## 0.03190105

### Second best model
gini <- vector()

a <- names(fixef(multivar4i.2ndbest))
a <- a[a!="(Intercept)"]
a <- paste(a, collapse=" + ")
f <- as.formula(paste0("present ~ ", a, " + (1|id_600m)"))

for(i in 1:10) { ## 10-fold cross-validation
  n.calib <- trunc(nrow(dat)*0.7) ## Number of data rows to use for calibration
  dat.calib <- dat[sample(nrow(dat), n.calib), ] ## Randomly selected calibration data
  dat.valid <- dat[!(rownames(dat) %in% rownames(dat.calib)),] ## Validation data, independent of calibration data
  
  m.calib <- glmer(f, family=binomial, dat.calib, na.action=na.fail) ## Model failed to converge
  p.valid <- predict(m.calib, dat.valid, re.form=NA, type="response")
  
  g <- Gini(p.valid, dat.valid$prop_pres)
  gini <- c(gini, g)
}
mean(gini) ## 0.7323672
sd(gini) ## 0.02155791

##### Calculate partial R-squared of each variable in the final model #####
### https://cran.r-project.org/web/packages/partR2/vignettes/Using_partR2.html

### Best model
Rsquare_m1 <- partR2(multivar4i.best, partvars=vars, max_level=1, # the argument max_level=1 calculates part R2 for each predictor, but not for all their combinations. We can specify this with .
                 R2_type = "marginal", nboot = 10) ## Marginal R2 refers to the variance explained by fixed effect predictors relative to the total variance in the respons

Rsquare_m1$R2   # Partial R2s Very low when random intercept of 600m grid-cell not included

### 2nd best model
Rsquare_m2 <- partR2(multivar4i.2ndbest, partvars=vars, max_level=1, # the argument max_level=1 calculates part R2 for each predictor, but not for all their combinations. We can specify this with .
                  R2_type = "marginal", nboot = 10) ## Marginal R2 refers to the variance explained by fixed effect predictors relative to the total variance in the respons

Rsquare_m2$R2   # Partial R2s Very low when random intercept of 600m grid-cell not included

## Combine and ouptut
Rsq_both <- cbind(Rsquare_m1$R2, Rsquare_m2$R2)
write.csv(Rsq_both, paste0(wd.out,"/Rsq_MountsBay.csv"), row.names=F)

library(patchwork)
p1 <- forestplot(Rsquare_m2, type = "R2", text_size = 10)
p2 <- forestplot(Rsquare_m2, type = "IR2", text_size = 10) ## Inclusive R2. This is SC^2 * R2_full.
p3 <- forestplot(Rsquare_m2, type = "SC", text_size = 10) ## Structure coefficients are the correlation between a predictor and the predicted response. Squared structure coefficients indicate how much of the regression effect can be attributed to a given predictor
p4 <- forestplot(Rsquare_m2, type = "BW", text_size = 10) ## Standardised model estimates (beta weights) for fixed effects. Beta weights for Gaussian models are calculated as beta * sd(x)/sd(y), with beta being the estimated slope of a fixed effect for predictor x and response y. Beta weight for Non-Gaussian models are calculated as beta * sd(x). Beta weights for interactions or polynomial terms are not informative at the moment and we recommend users to standardise variables themselves before fitting the model and to look at the model estimates (Ests) instead of beta weights (BW) in the partR2 output. See vignette for details.
(p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = "A")
## SCs are better than beta weights when there is multicollinearity in the explanatory variables. A low β weight and a high structure coefficient indicates multicollinearity, and beta weight is masked by multicollinearity

