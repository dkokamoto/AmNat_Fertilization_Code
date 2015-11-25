
##############################################
###  Script to estimate GLMMs using lme4   ###
###  Author:  D.K. Okamoto                 ###
##############################################

### load packages ###
library(lme4)

### read in data
fert <- read.csv("Okamoto_Fert.csv",header= T)  

### ensure replicate is a factor 
fert$Replicate <- factor(fert$Replicate)

### add individual ID to include random normal variance on the logit scale
fert$INDID <- factor(1:180)

### fit the generalized linear mixed effects models for total fertilization ###
fit1 <- glmer(TFert/100~poly(log10(S0),3)+log10(E0)+  ### fixed effects
              (1|Replicate)+                           ### random effects 
              (1|INDID),                               ### error variance
              family=binomial(link= "logit"),
              data=fert,weights= rep(100,180))
fit2 <-  update(fit1,.~.-log10(E0))
fit3 <-  update(fit1,.~.-poly(log10(S0),3))

### perform likelihood ratio tests ###
anova(fit1,fit2)
anova(fit1,fit3)

### fit the generalized linear mixed effects models for polyspermy ###
fit1p <- glmer(PFert/100~log10(S0)+log10(E0)+ ### fixed effects
               (1|Replicate)+                   ### random effects 
               (1|INDID),                       ### error variance 
               data=fert,
               family= "binomial",
               weights= rep(100,nrow(fert)),
               verbose= F)

fit2p <-  update(fit1p,.~.-log10(E0))
fit3p <-  update(fit1p,.~.-log10(S0))

### perform likelihood ratio tests ###
anova(fit1p,fit2p)
anova(fit1p,fit3p)

### extract coefficients for model of total fertilization
coef <- data.frame(fixef(fit1))
varcomp <- data.frame(summary(fit1)$varcor)[,c(1,5)]

### extract coefficients for model of polyspermy
coefp <- data.frame(fixef(fit1p))
varcompp <- data.frame(summary(fit1p)$varcor)[,c(1,5)]

### compute confidence intervals
profile <- confint(fit1)  ### total fertilization
profilep <- confint(fit1p) ### polyspermy

