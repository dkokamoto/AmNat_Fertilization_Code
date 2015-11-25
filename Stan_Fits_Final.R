######################################################
###  Script to sample model posteriors in Stan,    ###
###  calculate summaries and generate WAIC         ###
###  Author:  D.K. Okamoto                         ###
######################################################

require(rstan);library(gaussquad)
set_cppo(mode = "fast")

### read in data - must be saved as 'Okamoto_Fert.csv' ###
fert <-read.csv("Okamoto_Fert.csv")
attach(fert)

### load models provided in the file "Stan_Models_Final.R"
source("Stan_Models_Final.R")

### set number of replciates (R) and total number of samples (N) ###
R <- 8
N <- 180
  
### design matrices for collision rates 
Dmat1 <- cbind(D1,D2,D3,D4)  ### different collision rates by egg dilution
Dmat2 <- cbind(rep(1,180),rep(0,180),rep(0,180),rep(0,180))  ### one common collision rate 

### MCMC details ###
n.chains =3
n.iter =1000
n.burnin=500
set.seed <- 1234

### numerical integration times and weights
n.points <- 10  ### total number of points in each interval
n.intervals <- 3  ### total number of intervals
L <- n.points*n.intervals  ### total number of points in all intervals
a4 <- 7200;a3 <- 600;a2 <- 30;a1 <- 0  ### endpoints for each interval (seconds)
  
### calculate time points and weights
weights <- legendre.quadrature.rules(n.points)[[n.points]] 
  
### convert to actual time points (seconds) for numerical integration
t <- c(weights$x*(a2-a1)/2+(a2+a1)/2,
       weights$x*(a3-a2)/2+(a3+a2)/2,
       weights$x*(a4-a3)/2+(a4+a3)/2)
### calculate actual weights for numerical integration
w <- c(weights$w,weights$w,weights$w)*c(rep((a2-a1)/2,n.points),
       rep((a3-a2)/2,n.points),rep((a4-a3)/2,n.points))

### different contact times for the Millar & Anderson and Styan models###
tau1 <- rep(2310,180)  ### sperm half-life estimated in this study
tau2 <- rep(7200,180)  ### total contact time 

### Levitan's equation for half-life as a functon of sperm concentration
tau.lev <- 10^(log10(S0)*0.457+2.798) 

### data for the models with density dependent collision rates ###
   
data.ALL <-    list(Fert=MFert,S=S0,E=E0,R=R,N=N,Rep=Replicate)
data.FULL <-   append(data.ALL,list(D=Dmat1))
data.NULL <-   append(data.ALL,list(D=Dmat2))
data.Dyn <-    list(tau=rep(0,180),w=w,t=t,L=L,r_s=0.0003)       

### vector of parameter names to save ###
### for models with random collision rates
 
params.NULL <- c("b1_mu","sd_b1","b_r",
                 "delta_mu","delta_r","sd_delta",
                 "gamma_mu","gamma_r","sd_gamma",
                 "log_lik","p","p2","lambda","a_beta","b_beta")
  
### for the models with density-dependent collision rates
params <- c(params.NULL,"b2","b3","b4")

### ESTIMATE THE POSTERIORS ###

  ### dynamic model ###
  Dyn.FULL = stan(model_code = Dyn.model, 
               data=append(data.FULL,data.Dyn),
               pars=params,
               iter = n.iter, warmup= n.burnin,chains = n.chains,
               verbose = FALSE,init="random",seed= set.seed)

  Dyn.NULL = stan(model_code = Dyn.model, 
               data=append(data.NULL,data.Dyn),
               pars=params.NULL,
               iter = n.iter, warmup= n.burnin,chains = n.chains,
               verbose = FALSE,init="random",seed=set.seed)

  ### MA model with tau = estimated sperm half-life ###
  MA.FULL = stan(model_code = MA.model, 
               data=append(data.FULL,list(tau=tau1)),
               pars=params,
               iter = n.iter, warmup= n.burnin,chains = n.chains,
               verbose = FALSE,init="random",seed=set.seed)

  MA.NULL = stan(model_code = MA.model, 
               data=append(data.NULL,list(tau=tau1)),
               pars=params.NULL,
               iter = n.iter, warmup= n.burnin,chains = n.chains, 
               verbose = FALSE,init="random",seed=set.seed)

  ### MA model with tau =total contact time (7200 s) ###
  MA.FULL2 = stan(model_code = MA.model, 
               data=append(data.FULL,list(tau=tau2)),
               pars=params,
               iter = n.iter, warmup= n.burnin,chains = n.chains,
               verbose = FALSE,init="random",seed=set.seed)

  MA.NULL2 = stan(model_code = MA.model, 
                 data=append(data.NULL,list(tau=tau2)),
                 pars=params.NULL,
                 iter = n.iter, warmup= n.burnin,chains = n.chains, 
                 verbose = FALSE,init="random",seed=set.seed)


  ### MA model with tau dependent on sperm concentration  ###
  MA.FULL3 = stan(model_code = MA.model, 
                data=append(data.FULL,list(tau=tau.lev)),
                pars=params,
                iter = n.iter, warmup= n.burnin,chains = n.chains, 
                verbose = FALSE,init="random",seed=set.seed)

  MA.NULL3 = stan(model_code = MA.model, 
                data=append(data.NULL,list(tau=tau.lev)),
                pars=params.NULL,
                iter = n.iter, warmup= n.burnin,chains = n.chains, 
                verbose = FALSE,init="random",seed=set.seed)

  ### Styan model with tau = estimated sperm half-life ###
  Styan.FULL = stan(model_code = Styan.model, 
                data=append(data.FULL,list(tau=tau1)),
                pars=params,
                iter = n.iter, warmup= n.burnin,chains = n.chains,
                verbose = FALSE,init="random",seed=set.seed)

  Styan.NULL = stan(model_code = Styan.model, 
                data=append(data.NULL,list(tau=tau1)),
                pars=params.NULL,
                iter = n.iter, warmup= n.burnin,chains = n.chains, 
                verbose = FALSE,init="random",seed=set.seed)
  
  ### Styan model with tau =total contact time (7200 s) ###
  Styan.FULL2 = stan(model_code = Styan.model, 
                  data=append(data.FULL,list(tau=tau2)),
                  pars=params,
                  iter = n.iter, warmup= n.burnin,chains = n.chains, 
                  verbose = FALSE,init="random",seed=set.seed)

  Styan.NULL2 = stan(model_code = Styan.model, 
                  data=append(data.NULL,list(tau=tau2)),
                  pars=params.NULL,
                  iter = n.iter, warmup= n.burnin,chains = n.chains, 
                  verbose = FALSE,init="random",seed=set.seed)

  ### Styan model with tau dependent on sperm concentration  ###
  Styan.FULL3 = stan(model_code = Styan.model, 
                  data=append(data.FULL,list(tau=tau.lev)),
                  pars=params,
                  iter = n.iter, warmup= n.burnin,chains = n.chains, 
                  verbose = FALSE,init="random",seed=set.seed)

  Styan.NULL3 = stan(model_code = Styan.model, 
                  data=append(data.NULL,list(tau=tau.lev)),
                  pars=params.NULL,
                  iter = n.iter, warmup= n.burnin,chains = n.chains, 
                  verbose = FALSE,init="random",seed=set.seed)

### function to extract coefficient means and credible sets ###
params.summary <- function(mod){
  params <- extract(mod)
  params2 <- params[names(params)%in%c("b_r","gamma_r","delta_r")]
  params <- params[names(params)%in%c("b1_mu","b2","b3","b4","delta_mu","gamma_mu",
                                      "sd_b1","sd_delta","sd_gamma","lambda")]
    
   
  summarydf <- data.frame(list(mean=do.call(rbind,lapply(params,mean)),
    CSL= do.call(rbind,lapply(params,quantile,probs= 0.025)),
    CSU= do.call(rbind,lapply(params,quantile,probs= 0.975))))
    
  return(summarydf)
}

## apply function to extract mean and credible sets ###
param.results <- lapply(list(Millar_Anderson=MA.FULL,
                             Styan=Styan.FULL,
                             Dynamic= Dyn.FULL,
                             Styan2=Styan.FULL2,
                             MA2=MA.FULL2,
                             Styan3=Styan.FULL3,
                             MA3=MA.FULL3),
                        params.summary
                        )
param.results

### WAIC ###
  ### function to calculate WAIC 
waic <- function(mod){ 
  ### beta-binomial log-likelihood ###
  ll <- extract(mod)$log_lik  
  ### log posterior predictive loss ###
  lpd <- sum(log(colMeans(exp(ll))))  
  ### log posterior predictive loss ###
  p_waic <- sum(apply(ll,2,function(x) sum((x-mean(x))^2)/(length(x)-1))) 
  ### calculate waic and convert it to the deviance scale  ###
  return(-2*(lpd-p_waic))
}

### apply WAIC function to posterior samples
waic.results <- data.frame(lapply(list(Millar_Anderson=MA.FULL,
                                       Millar_Anderson.NULL=MA.NULL,
                                       Millar_Anderson3=MA.FULL3,
                                       Millar_Anderson3.NULL=MA.NULL3,
                                       Styan=Styan.FULL,
                                       Styan.NULL= Styan.NULL,
                                       Dynamic= Dyn.FULL,
                                       Dynamic.NULL= Dyn.NULL),
                                  waic)
                           )
waic.results

