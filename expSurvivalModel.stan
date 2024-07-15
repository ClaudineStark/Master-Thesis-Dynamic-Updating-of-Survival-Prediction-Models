// In this file:
// - data block: type of data is defined that will be used in rstan
// - parameter block: parameters are defined (coefficients and log bh)
// - model block: bayesian updating model
// - generated quantities block: surivial times are sampled (needed later on)

///////////////////////////////////////////////////////////////////////////////

//define type and dimensions of data provide through rstan
data {
    int<lower=1> N; //total number of individuals
    int<lower=1> N_uncensored; //number of individuals with event                                    
    int<lower=1> N_censored; //number of censored individuals                                       
    int<lower=0> NC; //number of covariates (allowed to be 0)    
    //splitting design matrix (X) into two parts                                           
    matrix[N_censored,NC] X_censored; //with event                                
    matrix[N_uncensored,NC] X_uncensored; //censored  
    matrix[N,NC] X; //full design matrix                         
    vector<lower=0>[N_censored] times_censored; //time of censoring                          
    vector<lower=0>[N_uncensored] times_uncensored; //time of event 
    //priors
    vector[NC] prior_mean_betas;
    vector[NC] prior_sd_betas;
    real prior_mean_log_lambda;
    real prior_sd_log_lambda;                     
}

///////////////////////////////////////////////////////////////////////////////

//all parameters of use
parameters {
    vector[NC] betas; //coefficients                                    
    real log_lambda; //log baseline hazard                                 
}

///////////////////////////////////////////////////////////////////////////////

//bayesian updating model
model {
  //priors:
  betas ~ normal(prior_mean_betas,prior_sd_betas);
  log_lambda ~ normal(prior_mean_log_lambda,prior_sd_log_lambda);

  //assumption: survival times are mutually independent
  //vecotrized statements more efficient that for loop :)
  
  //likelihoods:
  target += exponential_lpdf(times_uncensored | exp(log_lambda+X_uncensored*betas));
  target += exponential_lccdf(times_censored | exp(log_lambda+X_censored*betas)); 
}

///////////////////////////////////////////////////////////////////////////////

//sample surival times (for posterior predictive checks later)
generated quantities {
  vector[N_uncensored] times_uncensored_sampled;
  {
    //code for survtimes max horizon
    real tmp;
    real max_time;
    real max_time_censored;
    max_time = max(times_uncensored);
    max_time_censored = max(times_censored);
    if(max_time_censored > max_time) max_time = max_time_censored;
        
    for(i in 1:N_uncensored) {
        tmp= max_time + 1; 
        while(tmp > max_time) {
            tmp = exponential_rng(exp(log_lambda+X_uncensored[i,]*betas));
        }
        times_uncensored_sampled[i] = tmp;
    }
  }
  
}

