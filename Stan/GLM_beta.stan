// This is Stan code for a GLM with a beta-distributed response variable
// code by Mikey Tabak tabakma@gmail.com
// from: Tabak MA, Piaggio AJ, Miller RS, Sweitzer RA, Ernest HB (2017) Anthropogenic factors predict movement of an invasive species. Ecosphere 8: e01844.

data {
 int<lower=1> N;                   // sample size
 int<lower=1> K;                   // K predictors
 vector<lower=0,upper=1>[N] y;     // response 
 matrix[N,K] X;                    // predictor matrix
}
parameters {
 vector[K] theta;                  // coefficients
 real<lower=0> phi;                // dispersion parameter

}
transformed parameters{
 vector[K] beta;
 vector[N] mu;    // transformed linear predictor
 vector[N] Xbeta;                  // linear predictor

 beta = theta * 10;               // use this for logit link function
 Xbeta = X * beta;
 for (n in 1:N) {  
  mu[n] = inv_logit(Xbeta[n]); 
 }
}
model {
 vector[N] A;             // parameter for beta distn
 vector[N] B;             // parameter for beta distn
 //int y_pred[N];           // predicted y values for posterior check

 // model calculations
 A = mu * phi;
 B = (1.0 - mu) * phi;
 //for(n in 1:N){
 //y_pred[n] ~ beta(A, B); 
 //}

 // priors 
 //phi ~ cauchy(0, 5);               // different options for phi  
 phi ~ inv_gamma(.001, .001);
 //phi ~ uniform(0, 500);          // put upper on phi if using this
 //eps ~ normal(0,10000);

 for (k in 1:K){
  theta[k] ~ normal(0,100);
 }

 // likelihood
 y ~ beta(A, B);
}
generated quantities{
 vector[N] log_lik;       // for calculating waic
 real log_lik_model;      // log likelihood for the whole model
 //real y_pred;           // predicted y values for posterior check

 // predicted values based on estimated parameters
 //y_pred = beta(mu * phi, (1 - mu) * phi);

 // log likelihood of the model
 log_lik_model = beta_lpdf(y| mu*phi, (1-mu)*phi); 

 for (n in 1:N){
  log_lik[n] = beta_lpdf(y[n]| mu[n] * phi, (1.0 - mu[n]) * phi);
 }
}
