"
// Stan code to estimate temporal prevalence
// this model includes both true- and false-
// detection probabilities
// see Tabak et al. 2019 Ecology and Evolution 
// https://onlinelibrary.wiley.com/doi/10.1002/ece3.5558
// for details. 
data{
    int<lower=1> N;            // number of observations
    int<lower=0,upper=1> y[N]; // observations
    int<lower=1> K;            // number of demographic parameters (age and sex)
    matrix[N,K] X;             // model matrix of demographic information
    int<lower=1> N_yr;         // number of unique years
    matrix[N,N_yr] yrs;        // model matrix of years sampled
    int starts[N_yr];          // index in vector where each year starts
    int ends[N_yr];            // index in vector where each year ends
  }
parameters{
  vector[K] beta;              // coefficients for demographic effects
  real<lower=0> sigma_beta[K]; // variance for betas
  real<lower=0,upper=1> rho;   // sensitivity
  real<lower=0,upper=1> phi;   // specificity
}
transformed parameters{
  vector[N] psi;  // probability of pathogen prevalence  
  
  for(n in 1:N){
    psi[n] = inv_logit(X[n]*beta);
  }
}
model{
  // priors
  rho ~ beta(1,1); // flat priors on these parameters. Dont need to specify, but doing this for clarity
  phi ~ beta(1,1);
  sigma_beta ~ gamma(2,0.1);
  beta ~ normal(0,sigma_beta); // vague prior on betas
  
  // likelihood statement; marginalizing out the discrete parameter (z)
  for(n in 1:N){
    for(r in 1:N_yr){
      if(yrs[n,r]==1){
        target += log_mix(psi[n], bernoulli_lpmf(y[n]|rho),
                          bernoulli_lpmf(y[n]|(1-phi)));
      }
    }
  }
  
}
generated quantities{
  vector[N_yr] mean_psi_yrs; // the psi value in each year (0 for years where not sampled)
  real<lower=0> mean_psi_tot;
  real<lower=0,upper=1> p_step[N];
  vector[N] log_lik;
  
  mean_psi_tot = mean(psi);
  
  for(r in 1:N_yr){
    mean_psi_yrs[r] = mean(psi[starts[r]:ends[r]]); 
  }
  
  // calculating pi and log liklihood
  for(n in 1:N){
    for (r in 1:N_yr){
      if(yrs[n,r]==1){
        p_step[n] = step(mean_psi_yrs[r] - psi[n]);
        log_lik[n] = bernoulli_lpmf(y[n]|psi[n]*rho + (1-psi[n])*(1-phi));
      }
    }
  }
  
}"
