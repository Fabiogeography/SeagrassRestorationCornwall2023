##############################################################
##### CALCULATE AVERAGE OF TWO BEST SPATIAL MODELS ################################
##### Use the 2nd best model #####
##### Written by: Regan Early ################################
##### Written on: 21st Nov 2023 ##############################
##### Modified on: 12th Dec 2023   ###########################
##############################################################

.libPaths("C:/SOFTWARE/R-4.3.2/library")
library(MuMIn) ## model.sel
library(lme4)

### Models
wd.out <- "F:/NON_PROJECT/SEAGRASS/CORNWALL_COUNCIL/RESULTS/MULTIVARIATE_MODELS/SPATIAL"

load(paste0(wd.out,"/multivar4i_bestModel_MountsBay_rac_m1")) 
m1 <- multivar4i.rac

load(paste0(wd.out,"/multivar4i_bestModel_MountsBay_rac_m2")) 
m2 <- multivar4i.rac

rm(multivar4i.rac)

##### Find best model #####
model.sel(list(m1, m2)) ## m2 is the single best model, with a delta AIC of 17.3

write.csv(as.data.frame(fixef(m2)), paste0(wd.out, "/params_m2.csv"), row.names=T)

