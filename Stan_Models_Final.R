#################################################
###  script to specify models for Stan        ###
###  Author:  D.K. Okamoto                    ###
#################################################

### Dynamic model ###
Dyn.model <-  
"  
  data {
    int<lower=0> N; // number of observations points
    int<lower= 0> R; // number of replicates
    int<lower= 0> L; 	// number of interpolations for numerical integration 
    int<lower=1, upper = R> Rep[N]; // replicate ID
    int<lower= 0,upper=100> Fert[N]; // number of samples
    vector[N] E; // Egg concentration
    vector[N] S; // Sperm concentration
    matrix[N,4] D; //Design matrix for collision rates
    row_vector[L] w;  	// Gauss-Legendre quadrature weights
	  vector[L] t; 	// Gauss-Legendre quadrature evaluation points
  }

parameters {
  real<lower=0.00001,upper=0.1> b1_mu; // sperm attack rate at 1 egg per microliter
  real<lower=0.00001,upper=4> sd_b1; // sd in log sperm attack rate among pairs
  real<lower=0.00001,upper=0.1> b_r[R]; // pair level attack rate at 1 egg per microliter

  real<lower=0.01,upper=2> b2; // adjustment for 1/4 eggs per micrioliter
  real<lower=0.01,upper=2> b3; // adjustment for 1/16 eggs per micrioliter
  real<lower=0.01,upper=2> b4; // adjustment for 1/64 eggs per micrioliter

  real<lower=0.00001,upper=0.15> gamma_mu; // mean egg selectivity among pairs
  real<lower=0.00001,upper=4> sd_gamma; // sd in egg selectivity among pairs
  real<lower=0.00001,upper=0.15> gamma_r[R]; // pair level egg selectivity

  real<lower=0.00001,upper=20> delta_mu; // rate of polyspermy block
  real<lower=0.00001,upper= 4>  sd_delta; // rate of polyspermy block
  real<lower=0.00001,upper=20> delta_r[R]; // rate of polyspermy block

  real<lower=0.00001, upper= 2000> lambda; // total count (dispersion) parameters
  real<lower=0.00001,upper= .99999> p2[N]; //
}

transformed parameters {
  vector[N] b0; // estimated attack rate for each observation
  real <lower=0.00001,upper=.99999> p[N]; // monozygotic fertilization probability observation
  vector[N] log_lik; // likelihood for each observation
  vector[N] a_beta; // shape parameters 
  vector[N] b_beta;

  for (i in 1:N) {
    b0[i] <- b_r[Rep[i]]*(D[i,1]+b2*D[i,2]+b3*D[i,3]+b4*D[i,4]);     
  	
    // Numerical integration to generate predicted probability of monozygotic fertilization
	  p[i] <- fmax(fmin((w*(delta_r[Rep[i]]*(b0[i]*exp(-
      ((b0[i]*exp((-0.0003- b0[i]*E[i])*t)*gamma_r[Rep[i]]*S[i])/(-0.0003-b0[i]*E[i]))-
      (b0[i]*gamma_r[Rep[i]]*S[i])/(0.0003 + b0[i]*E[i])-
      (0.0003+b0[i]*E[i]-delta_r[Rep[i]])*t-delta_r[Rep[i]]*t).*
      (-1+exp((0.0003+b0[i]*E[i]-delta_r[Rep[i]])*t))*
      E[i]*gamma_r[Rep[i]]*S[i])/(0.0003+b0[i]*E[i]-delta_r[Rep[i]])))/E[i],
      0.9999),0.0001);

   // generate shape parameters and calculate likelihood
    a_beta[i] <-  p[i]*lambda;
    b_beta[i] <-  (1-p[i])*lambda;
    log_lik[i] <- beta_binomial_log(Fert[i],100,a_beta[i],b_beta[i]);
  }
}

model {
 lambda ~ pareto(1,1.5);
  
  // random effects
  for ( r in 1:R){
    b_r[r] ~      lognormal(log(b1_mu),sd_b1);
    gamma_r[r] ~  normal(gamma_mu,sd_gamma)T[.00001,.15];
    delta_r[r] ~  lognormal(log(delta_mu),sd_delta);
  }

  p2 ~ beta(a_beta,b_beta);
  Fert ~ binomial(100,p2); // define the likelihood
}
"

### Millar & Anderson model ###
MA.model <-  "

data {
  int<lower=0> N; // number of observations points
  int<lower= 0> R; // number of replicates
  int<lower=1, upper = R> Rep[N]; // replicate ID
  int<lower= 0,upper=100> Fert[N]; // number of samples
  vector[N] E; // Egg concentration
  vector[N] S; // Sperm concentration
  matrix[N,4] D; //Design matrix for collision rates
  real<lower=0> tau[N];
}

parameters {
real<lower=0.00001,upper=0.1> b1_mu; // sperm attack rate at 1 egg per microliter
  real<lower=0.00001,upper=4> sd_b1; // sd in log sperm attack rate among pairs
  real<lower=0> b_r[R]; // pair level attack rate at 1 egg per microliter

  real<lower=0.01,upper=2> b2; // adjustment for 1/4 eggs per micrioliter
  real<lower=0.01,upper=2> b3; // adjustment for 1/16 eggs per micrioliter
  real<lower=0.01,upper=2> b4; // adjustment for 1/64 eggs per micrioliter

  real<lower=0.00001,upper=0.15> gamma_mu; // mean egg selectivity among pairs
  real<lower=0.00001,upper=4> sd_gamma; // sd in egg selectivity among pairs
  real<lower=0.00001,upper=0.15> gamma_r[R]; // pair level egg selectivity

  real<lower=0.00001,upper=20> delta_mu; // rate of polyspermy block
  real<lower=0.00001,upper= 4>  sd_delta; // rate of polyspermy block
  real<lower=0.00001,upper=20> delta_r[R]; // rate of polyspermy block

  real<lower=0.00001, upper= 2000> lambda; // total count (dispersion) parameters
  real<lower=0.00001,upper= .99999> p2[N]; //
}

transformed parameters {
  vector[N] b0; // estimated attack rate for each observation
  real <lower=0.00001,upper=.99999> p[N]; // monozygotic fertilization probability observation
  vector[N] log_lik; // likelihood for each observation
  vector[N] a_beta; // shape parameters 
  vector[N] b_beta;
  vector[N] x1;
	vector[N] x2;
	vector[N] x3;

  for (i in 1:N) {
    b0[i] <- b_r[Rep[i]]*(D[i,1]+b2*D[i,2]+b3*D[i,3]+b4*D[i,4]);             
    x1[i] <- gamma_r[Rep[i]]*S[i]/E[i]*(1-exp(-b0[i]*E[i]*tau[i])); 
	  x2[i] <-  gamma_r[Rep[i]]*S[i]/E[i]*(1-exp(-b0[i]*E[i]*delta_r[Rep[i]]));
	  x3[i] <-  gamma_r[Rep[i]]*S[i]/E[i]*(1-exp(-b0[i]*E[i]*(tau[i]-delta_r[Rep[i]])));					
	  
    p[i] <- fmax(fmin((-((x1[i]-x3[i])*exp(-x1[i])-
      (exp(-x2[i])-exp(-x1[i]))*
      exp(b0[i]*E[i]*delta_r[Rep[i]]))),
      0.99999),0.00001);	

    // generate shape parameters and calculate likelihood
    a_beta[i] <-  p[i]*lambda;
    b_beta[i] <-  (1-p[i])*lambda;
    log_lik[i] <- beta_binomial_log(Fert[i],100,a_beta[i],b_beta[i]);
  }
}

model {
 lambda ~ pareto(1,1.5);
  
  // random effects
  for ( r in 1:R){
    b_r[r] ~      lognormal(log(b1_mu),sd_b1);
    gamma_r[r] ~  normal(gamma_mu,sd_gamma)T[.00001,.15];
    delta_r[r] ~  lognormal(log(delta_mu),sd_delta);
  }

  p2 ~ beta(a_beta,b_beta);
  Fert ~ binomial(100,p2); // define the likelihood
}
"

Styan.model <-  "

data {
  int<lower=0> N; // number of observations points
  int<lower= 0> R; // number of replicates
  int<lower=1, upper = R> Rep[N]; // replicate ID
  int<lower= 0,upper=100> Fert[N]; // number of samples
  vector[N] E; // Egg concentration
  vector[N] S; // Sperm concentration
  matrix[N,4] D; //Design matrix for collision rates
  real<lower=0> tau[N];
}

parameters {
  real<lower=0.00001,upper=0.1> b1_mu; // sperm attack rate at 1 egg per microliter
  real<lower=0.00001,upper=4> sd_b1; // sd in log sperm attack rate among pairs
  real<lower=0.00001,upper=0.1> b_r[R]; // pair level attack rate at 1 egg per microliter

  real<lower=0.01,upper=2> b2; // adjustment for 1/4 eggs per micrioliter
  real<lower=0.01,upper=2> b3; // adjustment for 1/16 eggs per micrioliter
  real<lower=0.01,upper=2> b4; // adjustment for 1/64 eggs per micrioliter

  real<lower=0.00001,upper=0.15> gamma_mu; // mean egg selectivity among pairs
  real<lower=0.00001,upper=4> sd_gamma; // sd in egg selectivity among pairs
  real<lower=0.00001,upper=0.15> gamma_r[R]; // pair level egg selectivity

  real<lower=0.00001,upper=20> delta_mu; // rate of polyspermy block
  real<lower=0.00001,upper= 4>  sd_delta; // rate of polyspermy block
  real<lower=0.00001,upper=20> delta_r[R]; // rate of polyspermy block

  real<lower=0.00001, upper= 2000> lambda; // total count (dispersion) parameters
  real<lower=0.00001,upper= .99999> p2[N]; //
}

transformed parameters {
  vector[N] b0; // estimated attack rate for each observation
  real <lower=0.00001,upper=.99999> p[N]; // monozygotic fertilization probability observation
  vector[N] log_lik; // likelihood for each observation;
  vector[N] a_beta; // shape parameters 
  vector[N] b_beta;
  vector[N] x1;

  for (i in 1:N) {
    b0[i] <- b_r[Rep[i]]*(D[i,1]+b2*D[i,2]+b3*D[i,3]+b4*D[i,4]);
    x1[i] <- gamma_r[Rep[i]]*S[i]/E[i]*(1-exp(-b0[i]*E[i]*tau[i]));
    p[i] <- fmax(fmin(1-exp(-x1[i])-(1-exp(-x1[i])-x1[i]*exp(-x1[i]))*
      (1-exp(-gamma_r[Rep[i]]*S[i]/E[i]*(1-exp(-b0[i]*E[i]*delta_r[Rep[i]])))),
      0.99999),0.00001);

// generate shape parameters and calculate likelihood
    a_beta[i] <-  p[i]*lambda;
    b_beta[i] <-  (1-p[i])*lambda;
    log_lik[i] <- beta_binomial_log(Fert[i],100,a_beta[i],b_beta[i]);
  }
}

model {
  lambda ~ pareto(1,1.5);
  
  // random effects
  for ( r in 1:R){
    b_r[r] ~      lognormal(log(b1_mu),sd_b1);
    gamma_r[r] ~  normal(gamma_mu,sd_gamma)T[.00001,.15];
    delta_r[r] ~  lognormal(log(delta_mu),sd_delta);
  }

  p2 ~ beta(a_beta,b_beta);
  Fert ~ binomial(100,p2); // define the likelihood
}
"