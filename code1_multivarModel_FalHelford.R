########################################################################################
##### Combine selected variables to make multivariate models for seagrass #####
##### Model Fal & Helford only 
##### Written by: Regan Early
##### Written on: 13th November 2023
##### Modified: 
########################################################################################

.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(tidyverse)
library(lme4) # Mixed models
library(MuMIn) # dredge() stdize() model.sel()
library(car) ## Variance Inflation Factor, vif
library(snow) ## cluster

##### Set Up #####
wd.dat <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/combined"
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS"

dat <- read.csv(paste0(wd.dat,"/cornwall_sg_env_unique_30m_stnd.csv")) ## Data already standardized
dat <- dat[dat$site=="FalHelford",]

dat <- dat %>% ## Models will not run with NA values in data
  drop_na()

### Response variable is number of presence and number of absence in each 30m grid cell
# Response is the number of presence and absence points within a 30m grid cell
# Bind response var into single matrix
# y <- cbind(dat$total_pres, dat$total_abs)

######################################################
##################### multivar1 ######################
######################################################

### Just use all variables - not that many
vars <- c("bathymetry", "slope", "max_sst", "avg_sst", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_litt_perc", "eunis_coar_perc", "eunis_sand_perc", "eunis_mix_perc") ## Do not include "eunis_mud_perc" (no points within this habitat)
vars <- paste(vars, collapse=" + ")

multivar1 <- glmer(as.formula(paste0("cbind(total_pres, total_abs) ~ ", vars, " + (1|id_600m)")), family=binomial, dat, na.action=na.fail) ## Model failed to converge
summary(multivar1) ## SD of coarse grid-cell is 0.895 This is in the ball park than the estimates of the fixed effects, so retain coarse grid-cell
summary(multivar1)$AICtab ## 6917

### multivar1 did not converge. However, if all optimizers converge to values that are practically equivalent, then the model fit is good enough 
multivar1.all <- allFit(multivar1)
ss <- summary(multivar1.all)
ss$ fixef               ## fixed effects
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors of the random effect covariance matrix
ss$ which.OK            ## which fits worked
### Values are all very similar (except with Nelder Mead) so ignore the convergence warnings for now.

## Compare models with and without random effect of coarse grid-cell.
multivar1.noRFX <- glm(as.formula(paste0("cbind(total_pres, total_abs) ~ ", vars)), family=binomial, dat, na.action=na.fail) 
AIC(multivar1, multivar1.noRFX) ## The model with random effects has notably lower AIC than the model without. 
cbind(fixef(multivar1), coefficients(multivar1.noRFX))
## The coefficients of the fixed effects are fairly similar - eunis_other changes coefficient, but effect size is very small. Random effect might not affect parameter estimates much, though could affect model fit.

### Check collinearity.
vif(multivar1) 
## Try to keep both eunis_sand_perc and eunis_mix_perc - they both have important effects in univariate models.
## covar_spm and avg_spm seem to covary a bit. Keep both for now.

#####################################################
#################### multivar2 ######################
#####################################################

### Dredge the multivariate base model. Repects marginality constraints. Ranks by AICc by default.
## Set up the cluster. Number of cores limited by computing resources - best on server
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 30), type = clusterType))
clusterExport(clust, "dat") ## Send the data to each node of the cluster
clusterCall(clust, function() library(MuMIn)) ## Call the function to load the library on each node of the cluster
clusterCall(clust, function() library(lme4)) ## Call the function to load the library on each node of the cluster

## Dredge the mixed model
multivar2 <- dredge(multivar1, beta="sd", evaluate=T, trace=2, cluster=clust) ## Dredge the models without quadratic terms as otherwise the process is too lengthy. 
## Must be done on cluster in order to complete in reasonable time frame

colnames(multivar2)[1] <- "Intercept" ## write.csv doesn't like the column name having parentheses
write.csv(multivar2, file=paste0(wd.out,"/multivar2_FalHelford.csv"), row.names=F)

multivar2[multivar2$delta<=2,] ## 3 models fall into the best model subset. The best model has a delta AIC of -0.14. EUNIS other and max_sst are the variables not in all models.

multivar2.best <- get.models(multivar2, subset=1)[[1]]
save(multivar2.best, file=paste0(wd.out,"/multivar2_bestModel_FalHelford"))
# load(paste0(wd.out,"/multivar2_bestModel_FalHelford")) ## Dredging had to be done on the cluster. Load here.

### Investigate the best model. Check VIF of single best model.
vif(multivar2.best, singular.ok = TRUE) ## Try to keep both eunis_sand_perc and eunis_mix_perc - they both have important effects in univariate models.

##################################
##########  multivar2q ###########
##################################

### Add quadratic terms of the variables added in multivar2
## Use multivar1 as a model containing all variables was included in the best model subset

### Check quadratic terms ('new.vars') of existing main effects ('orig.vars')
multivar2.vars <- names(summary(multivar1)$coefficients[-1,1])
multivar2.vars <- paste(multivar2.vars, collapse=" + ")

## Make the models to be compared
multivar2q.vars <- paste0("I(", names(summary(multivar1)$coefficients[-1,1]), "^2)")
for (i in 1:length(multivar2q.vars)) {
  assign(paste0("f",i), as.formula(paste0("cbind(total_pres, total_abs) ~ ",multivar2.vars, " + ", multivar2q.vars[i], " + (1|id_600m)")))
  assign(paste0("m",i), glmer(get(paste0("f",i)), family=binomial, dat, na.action=na.fail))
} 

## Compare models using AICc
multivar2q.vars.sel <- model.sel(mget(c(paste0("m",1:length(multivar2q.vars)), "multivar1"))) 

## Get the explanatory variables in the models below the delta AICc threshold.
multivar2q <- get.models(multivar2q.vars.sel, subset=delta<2) ## There was one best model, m1, that included all linear terms and the quadratic term of bathymetry (delta AIC -9.42)

names(multivar2q) <- "multivar2q"

## Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(multivar2q)) {print(summary(multivar2q[[i]]))} 

## Check VIF
for(i in 1:length(multivar2q)) {print(vif(multivar2q[[i]], singular.ok = TRUE))} ## Try to keep both eunis_sand_perc and eunis_mix_perc - they both have important effects in univariate models.

## Check that a non-mixed model wouldn't perform better
multivar2q.vars.all <- unique(unlist(lapply(multivar2q , function(x) {rownames(summary(x)$coefficients)})))
multivar2q.vars.all <- multivar2q.vars.all[multivar2q.vars.all!="(Intercept)"]
multivar2q.vars.all <- paste(multivar2q.vars.all, collapse=" + ")

f <- as.formula(paste0("cbind(total_pres, total_abs) ~ ", multivar2q.vars.all))
multivar2q.noRFX <- glm(f, family=binomial, dat, na.action=na.fail)
model.sel(multivar2q[[1]], multivar2q.noRFX) ## Model with random effects has a substantially lower AIC than teh model without (delta AIC -757). Continue with mixed model. 

## Save model
multivar2q.best <- multivar2q[[1]]
save(multivar2q.best, file=paste0(wd.out,"/multivar2q_bestModel_FalHelford"))

### Mixed model did not converge. Check results are similar using all optimisers.
multivar2q.best.all <- allFit(multivar2q.best)
ss <- summary(multivar2q.best.all)
ss$ fixef               ## fixed effects - estimates all very similar. OK to proceed with this model.
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors of the random effect covariance matrix
ss$ which.OK            ## which fits worked
### Values are all very similar (except with Nelder Mead) so ignore the convergence warnings.


###########################################################
######################### multivar4i #######################
###########################################################

### Check interactions between linear terms of all variables retained in multivar2q
multivar4 <- multivar2q.best

### Identify all the variables in multivar4 best model
a <- names(fixef(multivar4)) ## Get all the variables in teh best model 
a <- a[a!="(Intercept)"]
al <- a[!a %in% grep("I", a, value = T)] ## linear terms only

## Create a vector of all the interaction combos
ints <- as.data.frame(lapply(1:2, function(y) combn(al, y)))
ints <- as.vector(apply(ints, 2, paste, collapse=" * "))
ints <- ints[-c(1:length(a))] ## Remove combos of the same variable

## Fit all models with one interaction
multivar4.vars <- paste(a, collapse=" + ")
fs <- lapply(ints, function(x) {as.formula(paste0("cbind(total_pres, total_abs) ~ ", multivar4.vars, " + ", x, " + (1|id_600m)"))}) ## Make one formula for each interaction added to the existing fixed effects
multivar4i <- lapply(fs, function(x) {glmer(x, family=binomial, dat, na.action=na.fail)}) ## Make models of the above formulas. Many failed convergences. "fixed-effect model matrix is rank deficient so dropping 1 column / coefficient"
names(multivar4i) <- ints

multivar4i.sel <- model.sel(c(multivar4i, multivar4)) ## 8 warnings: 1: In vcov.merMod(object, use.hessian = use.hessian) : variance-covariance matrix computed from finite-difference Hessian is not positive definite or contains NA values: falling back to var-cov estimated from RX
## Single best model with delta AIC of -177
multivar4i  <- get.models(multivar4i.sel , subset=delta<2) ## Single best model containing interaction between batymetry amd covariance in spm. SD of random effect is ~1, and fixed effect estimates are between 0.4 and -8. So random effect still seems reasonable to include (best model without random effect contained the same linear and quadratic effects byt no interactions)
multivar4i.best <- multivar4i[[1]]

### Compare to model without random effect of coarse grid-cell
a <- names(fixef(multivar4i.best)) ## Get all the variables in teh best model 
a <- a[a!="(Intercept)"]
a <- paste0(a, collapse=" + ")

f <- as.formula(paste0("cbind(total_pres, total_abs) ~ ", a))
m <- glm(f, family=binomial, dat, na.action=na.fail) 

## Compare models with fixed and random effects.
AIC(multivar4i.best) ## 6279.291
AIC(m) ## 7295.989
## Model with random effect has much lower AIC (delta AIC ~100). Models fairly similar (estimates same order of magnitude and sign), but differences in parameter estimates may well affect prediction. 

## Check VIF
vif(multivar4i.best, singular.ok = TRUE) ## Sand and mix still a bit high.

## Save the best model 
save(multivar4i.best, file=paste0(wd.out,"/multivar4i_bestModel_FalHelford"))
# load(paste0(wd.out,"/multivar4i_bestModel_FalHelford"))

#################################################
###### Produce & save final model ###############
#################################################
### 4i produced a single best model so use that

final <- multivar4i.best

final.coefs <- fixef(final)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

## Write outputs
write.csv(final.coefs, file=paste0(wd.out,"/multivar4i_bestModel_FalHelford_coefs.csv"), row.names=F)

