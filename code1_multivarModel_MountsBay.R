########################################################################################
##### Combine environmental variables to make multivariate models for seagrass #####
##### Model Mounts Bay only 
##### Written by: Regan Early
##### Written on: 13th November 2023
########################################################################################



.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(tidyverse)
library(lme4) # Mixed models
library(MuMIn) # dredge() stdize() model.sel()
library(car) ## Variance Inflation Factor, vif
library(snow) ## for running on cluster

##### Set Up #####
wd.dat <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/DATA/variables/combined" ## Directory containing environmental variables
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS" ## Directory where results will be stored

dat <- read.csv(paste0(wd.dat,"/cornwall_sg_env_unique_30m_stnd.csv")) ## Data already standardized
dat <- dat[dat$site=="MountsBay",]

dat <- dat %>% ## Models will not run with NA values in data
  drop_na()
# Response variable is whether the 30m grid cell is occupied or not ('Presence/Absence'), i.e. column "present"

######################################################
##################### multivar1 ######################
######################################################

### Combine all variables in a single model. Linear terms only.
vars <- c("bathymetry", "slope", "max_sst", "avg_sst", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_coar_perc", "eunis_sand_perc") ## Do not include: "eunis_litt_perc" (very few values), "eunis_mix_perc" (only one value), "eunis_mud_perc" (no values within this region)
vars <- paste(vars, collapse=" + ")

multivar1 <- glmer(as.formula(paste0("present ~ ", vars, " + (1|id_600m)")), family=binomial, dat, na.action=na.fail) ## Model failed to converge. See line 42.

summary(multivar1) ## SD of coarse grid-cell is 1.44 This is in the ball park of the estimates of some fixed effects, so retain coarse grid-cell as random effect
summary(multivar1)$AICtab ## 2056

### multivar1 did not converge. However, if all optimizers converge to values that are practically equivalent, then the model fit is good enough 
multivar1.all <- allFit(multivar1)
ss <- summary(multivar1.all)
ss$ fixef               ## fixed effects
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors of the random effect covariance matrix
ss$ which.OK            ## which fits worked
### Differences in fixed effects are noticable

## Compare models with and without random effect of coarse grid-cell.
multivar1.noRFX <- glm(as.formula(paste0("present ~ ", vars)), family=binomial, dat, na.action=na.fail) 
AIC(multivar1, multivar1.noRFX) ## The model with random effects has notably lower AIC than the model without, so retain random effect.
cbind(fixef(multivar1), coefficients(multivar1.noRFX))
## The coefficients of the fixed effects are fairly similar - expo changes coefficient, but effect size is very small. Random effect might not affect parameter estimates much.

### Check collinearity
vif(multivar1) 
## max_sst has vif of 8.4. avg_sst 6.7. Try removing max_sst
multivar1b <- update(multivar1, ~. - max_sst) ## Model failed to converge
summary(multivar1b) ## SD of coarse grid-cell is 1.6 This is the ball park of the estimates of some fixed effects, so retain coarse grid-cell
summary(multivar1b)$AICtab ## 2058. Slightly higher than with max_sst but within delta AIC of 2.
vif(multivar1b) ## Fine

### multivar1 did not converge. However, if all optimizers converge to values that are practically equivalent, then the model fit is good enough 
multivar1b.all <- allFit(multivar1b)
ss <- summary(multivar1b.all)
ss$ fixef               ## fixed effects
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors of the random effect covariance matrix
ss$ which.OK            ## which fits worked
### Differences in fixed effects noticable, avg_sst and avg_spm especially. However, given AIC values retain random effect.

## Update vars
vars <- c("bathymetry", "slope", "avg_sst", "mlw_dist", "avg_spm", "covar_spm", "expo",
          "eunis_other_perc", "eunis_coar_perc", "eunis_sand_perc")
vars <- paste(vars, collapse=" + ")

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
multivar2 <- dredge(multivar1b, beta="sd", evaluate=T, trace=2, cluster=clust) ## Dredge the models without quadratic terms as otherwise the process is too lengthy. 
save(multivar2, file=paste0(wd.out,"/multivar2_MountsBay"))
## Must be done on cluster in order to complete in reasonable time frame
# load(paste0(wd.out,"/multivar2_MountsBay")) ## Dredging had to be done on the cluster. Load here.

colnames(multivar2)[1] <- "Intercept" ## write.csv doesn't like the column name having parentheses
write.csv(multivar2, file=paste0(wd.out,"/multivar2_MountsBay.csv"), row.names=F)

multivar2_subset <- multivar2[multivar2$delta<=2,] ## 14 models fall into the best model subset. The best model has a delta AIC of -0.23. EUNIS_coar, avg_sst, mlw_dist are the variables not in all models.
save(multivar2_subset, file=paste0(wd.out,"/multivar2_bestSubset_MountsBay"))
# load(paste0(wd.out,"/multivar2_bestSubset_MountsBay")) ## Dredging had to be done on the cluster. Load here.

multivar2.best <- get.models(multivar2, subset=1)[[1]]
save(multivar2.best, file=paste0(wd.out,"/multivar2_bestModel_MountsBay"))
# load(paste0(wd.out,"/multivar2_bestModel_MountsBay")) ## Dredging had to be done on the cluster. Load here.

### Investigate the best model. Check VIF of single best model.
vif(multivar2.best, singular.ok = TRUE) ## OK

##################################
##########  multivar2q ###########
##################################

### Test addition of quadratic terms of the variables individually
## Use multivar1b as a model containing all of variables in that model was included in the best model subset

### Check quadratic terms ('new.vars') of existing main effects ('orig.vars')
multivar2.vars <- names(summary(multivar1b)$coefficients[-1,1])
multivar2.vars <- paste(multivar2.vars, collapse=" + ")

## Make the models to be compared
multivar2q.vars <- paste0("I(", names(summary(multivar1)$coefficients[-1,1]), "^2)")
for (i in 1:length(multivar2q.vars)) {
  assign(paste0("f",i), as.formula(paste0("present ~ ",multivar2.vars, " + ", multivar2q.vars[i], " + (1|id_600m)")))
  assign(paste0("m",i), glmer(get(paste0("f",i)), family=binomial, dat, na.action=na.fail))
} 

## Compare models using AICc
multivar2q.vars.sel <- model.sel(mget(c(paste0("m",1:length(multivar2q.vars)), "multivar1b"))) 

## Get the explanatory variables in the models below the delta AICc threshold.
multivar2q <- get.models(multivar2q.vars.sel, subset=delta<2) ## There was one best model, m1, that included the quadratic term of bathymetry (delta AIC -226).

names(multivar2q) <- "multivar2q"

## Inspect models (ask whether fixed effect relationships meaningful)
for(i in 1:length(multivar2q)) {print(summary(multivar2q[[i]]))} 

## Check VIF
for(i in 1:length(multivar2q)) {print(vif(multivar2q[[i]], singular.ok = TRUE))} ## ok

## Check that a non-mixed model wouldn't perform better
multivar2q.vars.all <- unique(unlist(lapply(multivar2q , function(x) {rownames(summary(x)$coefficients)})))
multivar2q.vars.all <- multivar2q.vars.all[multivar2q.vars.all!="(Intercept)"]
multivar2q.vars.all <- paste(multivar2q.vars.all, collapse=" + ")

f <- as.formula(paste0("present ~ ", multivar2q.vars.all))
multivar2q.noRFX <- glm(f, family=binomial, dat, na.action=na.fail)
model.sel(multivar2q[[1]], multivar2q.noRFX) ## Model with random effects has a substantially lower AIC than the model without (delta AIC -155). Continue with mixed model. 

## Save model
multivar2q.best <- multivar2q[[1]]
save(multivar2q.best, file=paste0(wd.out,"/multivar2q_bestModel_MountsBay"))

### Mixed model did not converge. Check results are similar using all optimisers.
multivar2q.best.all <- allFit(multivar2q.best)
ss <- summary(multivar2q.best.all)
ss$ fixef               ## fixed effects - estimates all very similar. OK to proceed with this model.
ss$ llik                ## log-likelihoods
ss$ sdcor               ## SDs and correlations
ss$ theta               ## Cholesky factors of the random effect covariance matrix
ss$ which.OK            ## which fits worked
### Values are now fairly very similar so ignore the convergence warnings

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
fs <- lapply(ints, function(x) {as.formula(paste0("present ~ ", multivar4.vars, " + ", x, " + (1|id_600m)"))}) ## Make one formula for each interaction added to the existing fixed effects
multivar4i <- lapply(fs, function(x) {glmer(x, family=binomial, dat, na.action=na.fail)}) ## Make models of the above formulas. Many failed convergences. 
names(multivar4i) <- ints

multivar4i.sel <- model.sel(c(multivar4i, multivar4)) ## 2 warnings: 1: In vcov.merMod(object, use.hessian = use.hessian) : variance-covariance matrix computed from finite-difference Hessian is not positive definite or contains NA values: falling back to var-cov estimated from RX
## two best models both with delta AIC of -~ -6
multivar4i  <- get.models(multivar4i.sel , subset=delta<2) ## expo:eunis_other_perc and bathymetry:avg_spm. Both highly significant. SD of random effect is ~1.3, and fixed effect estimates are between 0.1 and -8. So random effect still seems reasonable to include.)
multivar4i.best <- multivar4i[[1]]

### Compare to model without random effect of coarse grid-cell
a <- names(fixef(multivar4i.best)) ## Get all the variables in teh best model 
a <- a[a!="(Intercept)"]
a <- paste0(a, collapse=" + ")

f <- as.formula(paste0("present ~ ", a))
m <- glm(f, family=binomial, dat, na.action=na.fail) 

## Compare models with and without random effects.
AIC(multivar4i.best) ## 1761
AIC(m) ## 1920
## Model with random effect has much lower AIC (delta AIC ~100) so retain. Models fairly similar (estimates same order of magnitude and sign).

## Check VIF
vif(multivar4i.best, singular.ok = TRUE) ## ok

## Save the two best models
save(multivar4i.best, file=paste0(wd.out,"/multivar4i_bestModel_MountsBay"))
multivar4i.2ndbest <- multivar4i[[2]]; save(multivar4i.2ndbest, file=paste0(wd.out,"/multivar4i_2ndbestModel_MountsBay"))
save(multivar4i, file=paste0(wd.out,"/multivar4i_MountsBay"))

#################################################
###### Produce & save final model ###############
#################################################
load(paste0(wd.out,"/multivar4i_MountsBay"))
final <- model.avg(multivar4i, subset = delta<2) ## Average the models with delta AIC < value.

final.coefs <- final$coefficients["full",] ## Average parameter estimates over all models (value will be 0 in models wehre the parameter wasn't included)
vars <- names(final.coefs)
coefs <- as.vector(final.coefs)
final.coefs <- data.frame(Variables=vars, MeanCoefficients=coefs)
rownames(final.coefs) <- vars

### Output the coefficients in a table
## The variable names are the same in both models
vars <- names(unlist(fixef(multivar4i[[1]])))
coefs <- as.vector(unlist(fixef(multivar4i[[1]])))
d1 <- data.frame(Variables=vars, Coefficients1=coefs)

vars <- names(unlist(fixef(multivar4i[[2]])))
coefs <- as.vector(unlist(fixef(multivar4i[[2]])))
d2 <- data.frame(Variables=vars, Coefficients2=coefs)

all.coefs <- Reduce(function(x,y) merge(x=d1, y=d2, by = "Variables", all=T), 
                    list(d1, d2)) ## Can presumably be modified to deal with all elements in a list.
colnames(all.coefs) <- c("Variables", "CoefsModel1", "CoefsModel2")

## Write outputs
save(final, file=paste0(wd.out,"/multivar4i_final_MountsBay"))
write.csv(all.coefs, file=paste0(wd.out,"/multivar4i_final_MountsBay_coefs.csv"), row.names=F)
