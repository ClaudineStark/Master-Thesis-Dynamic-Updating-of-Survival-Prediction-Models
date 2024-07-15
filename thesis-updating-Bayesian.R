# In this file:
# Bayesian updating is performed.
# Updating model is specified in separate Stan file.

#------------------------------------------------------------------------------
#Load necessary libraries
library(rstan)

#-------------------------------------------------------------------------------
#First, the set up

#For execution on a local, multicore CPU with excess RAM
options(mc.cores = parallel::detectCores()) 
#To avoid recompilation of unchanged Stan programs
rstan_options(auto_write = TRUE) #standard needed
#The stan model is in a separate .stan file which must be loaded
sm <- stan_model(file="expSurvivalModel.stan")

CORES <- 2
CHAINS <- 2  #low chains and iterations for speed in this example
ITER <- 4000 
BURNIN <- 1000
lambdaF <- 0.9  #the forgetting factor #TODO chose!
set.seed(12345)

#original model
#model_w1

#-------------------------------------------------------------------------------
# Bayesian updating functions:

#function that performs Bayesian update based on newdata and previous model
bayes_update <- function(newdata, previousmodel){
  
  #first adapting newdata to be given to stan
  dataUpdate <- newdata
  dataUpdate <- rename(dataUpdate, sexM = sex) |>
    mutate(sexM = if_else(sexM == "M", 1, 0))
  #covariates
  covar <- c("age", "sexM", "numb_comorb", "rr_min", 
             "sp_o2", "crp_mg_l", "creatinine_ser_mumol_l")
  X <- dataUpdate[, covar] #covariate matrix
  time <- dataUpdate$survtime
  event <- dataUpdate$status
  N_censored <- length(which(dataUpdate$status==0)) #number of individuals censored
  N <- nrow(X) #total number of individuals
  
  #setting priors
  priors <- list()
  
  #when previous model is not model format but output from previous update
  if(is.data.frame(previousmodel)){ 
    coef <- c("betas[1]", "betas[2]", "betas[3]", "betas[4]", "betas[5]", 
              "betas[6]", "betas[7]")
    #extracting only coefficients (betas)
    df <- previousmodel[(rownames(previousmodel) %in% coef),]
    priors$location <- df$mean #mean
    n <- length(priors$location)
    priors$scale <- df$sd/sqrt(lambdaF) #sd
  }
  #when previous model is model from wave 1 (model format)
  else{
    priors$location <- coefficients(previousmodel) 
    n <- length(priors$location)
    priors$scale <- vector(length=n)
    for(i in 1:n){
      priors$scale[i] <- previousmodel$var[i,i]/lambdaF #variance
    }
    priors$scale <- sqrt(priors$scale) #standard deviation
  }
  
  #can set an uninformative prior on the rate parameter
  prior_int <- list(location=0, scale=2.5) #TODO (make choice)
  
  #creating Stan data
  stan_data <- list(N=N, N_uncensored=N-N_censored,
                    N_censored=N_censored, 
                    X_censored=as.matrix(X[event==0,]),
                    X_uncensored=as.matrix(X[event!=0,]), X=X,
                    times_censored=time[event==0],
                    times_uncensored=time[event!=0],
                    NC=ncol(X),
                    prior_mean_log_lambda = prior_int$location, #prior log bh
                    prior_sd_log_lambda = prior_int$scale,
                    prior_mean_betas = priors$location, #priors coefficients
                    prior_sd_betas = priors$scale)
  
  #updating
  updateBayes <- sampling(sm, 
                          data=stan_data, 
                          chains=CHAINS, 
                          iter=ITER, 
                          cores=CORES, 
                          warmup=BURNIN) 
  
  #returning updated coefficients+lambda and sampled survival times from updated model
  return(updateBayes)
  
}

#updating
bayesian <- function(data){
  previousmodel <- NULL
  update <- NULL
  
  #updating:
  for (t in 1:length(data)) {
    if(t==1) previousmodel[[t]] <- model_w1 #first update uses coeff from model wave 1
    update[[t]] <- bayes_update(newdata = data[[t]], previousmodel = previousmodel[[t]]) #new data, coeff from previous update
    update[[t]] <- summary(update[[t]])$summary
    #extracting updated mean and sd and n_eff and Rhat (all rows!)
    previousmodel[[t+1]] <- as.data.frame(update[[t]][,c(1,3,9,10)]) 
  }
  return(previousmodel) 
}



#----------------------------------------------------------------------------
#Updating coefficients and log bh: 

#updating every 100 events:
update100_fit0 <- bayesian(data = df_int100_fit0) #no sliding window
update100_fit1 <- bayesian(data = df_int100_fit1) #1 week sliding window
update100_fit2 <- bayesian(data = df_int100_fit2) #2 weeks sliding window
update100_fit4 <- bayesian(data = df_int100_fit4) #4 weeks sliding window

bayes_data_list_100 <- list(update100_fit0, update100_fit1,
                            update100_fit2, update100_fit4)
names(bayes_data_list_100) <- c("sw0", "sw1", "sw2", "sw4")
#save this output
save(bayes_data_list_100, file =  "bayes_data_list_100.rda") 

#updating every 50 events:
update50_fit0 <- bayesian(data = df_int50_fit0) 
update50_fit1 <- bayesian(data = df_int50_fit1)
update50_fit2 <- bayesian(data = df_int50_fit2)
update50_fit4 <- bayesian(data = df_int50_fit4)

bayes_data_list_50 <- list(update50_fit0, update50_fit1,
                           update50_fit2, update50_fit4)
names(bayes_data_list_50) <- c("sw0", "sw1", "sw2", "sw4")
#save this output
save(bayes_data_list_50, file = "bayes_data_list_50.rda")

#updating every 150 events:
update150_fit0 <- bayesian(data = df_int150_fit0)
update150_fit1 <- bayesian(data = df_int150_fit1)
update150_fit2 <- bayesian(data = df_int150_fit2)
update150_fit4 <- bayesian(data = df_int150_fit4)

bayes_data_list_150 <- list(update150_fit0, update150_fit1,
                           update150_fit2, update150_fit4)
names(bayes_data_list_150) <- c("sw0", "sw1", "sw2", "sw4")
save(bayes_data_list_150, file = "bayes_data_list_150.rda")


#---------------------------------------------------------------------------
#Looking at the coefficients from the updates throughout plots:

#extracting coeff from updates as dataframe 
bay_reshape <- function(data){
  df_large <- data.frame()
  for (i in 1:length(data)) {
    df_tot <- data.frame()
    data[[i]] <- data[[i]][-1] #coef original model not needed now
    for (j in 1:(length(data[[i]])-1)) {
      df_temp <- data.frame()
      df_temp <- data[[i]][[j]][1:8,]
      df_temp$update <- j
      df_tot <- rbind(df_tot, df_temp)
    }
    df_tot$slidingwindow <- i-1
    df_large <- rbind(df_large, df_tot)
  }
  df_large$update <- as.factor(df_large$update) #needed for plots
  df_large$slidingwindow <- case_when(
    df_large$slidingwindow == 3 ~ "4 weeks",
    df_large$slidingwindow == 2 ~ "2 weeks",
    df_large$slidingwindow == 1 ~ "1 week",
    TRUE ~ "0 weeks" #when slidingwindow == 0
    )
  df_large$slidingwindow <- as.factor(df_large$slidingwindow)
  return(df_large)
  #return(df_tot)
}

#extracting rows of specific coefficient (beta)
select_coeff <- function(data, coeff_name){
  df <- bay_reshape(data = data)
  coef <- c("betas[1]", "betas[2]", "betas[3]", "betas[4]", "betas[5]", 
            "betas[6]", "betas[7]", "log_lambda")
  index <- which(coef %in% coeff_name)
  df_coef <- df[seq(index, nrow(df), 8),]
}

#check convergences
test_convergence_coef <- function(data){
  result <- c()
  for (i in 1:nrow(data)) {
    if(data$Rhat[i] > 1.1) result[i] <- "bad"
    else result[i] <- "good"
  }
  ok <- which(result=="bad")
  return(ok)
}

#making plot
plot_function_bay <- function(covar, data){
  df_new <- select_coeff(data = data, coeff_name = covar)
  #change names betas to original coefficient names
  index <- which(all_coeff == covar) 
  covar <- original_names[index]
  if(is.na(covar)) covar <- "log_lambda"
  
  plt <- ggplot(data = df_new,
                aes(y= mean, x= update, col = slidingwindow))+
    geom_point(position = position_dodge(width=0.4))+
    geom_linerange(aes(ymin= mean-1.96*sd, ymax=mean+1.96*sd),
                   position = position_dodge(width = 0.4))+
    labs(y= "coefficient", x ="updates",
         title = covar)
  return(plt)
}


all_coeff <- c("betas[1]", "betas[2]", "betas[3]", "betas[4]", "betas[5]", 
                      "betas[6]", "betas[7]", "log_lambda")
original_names <- c("age", "sex", "numb_comorb", "rr_min", "sp_o2", "crp_mg_l",
                    "creatinine_ser_mumol_l")
#plots:
#50 events
plot50 <- lapply(all_coeff, function(covar) plot_function_bay(covar = covar, data = bayes_data_list_50))
plot50_allcoef <- ggarrange(plot50[[1]], plot50[[2]], plot50[[3]],
                             plot50[[4]], plot50[[5]], plot50[[6]],
                             plot50[[7]], plot50[[8]],
                             nrow = 4, ncol = 2, 
                             common.legend = TRUE, legend="bottom")
annotate_figure(plot50_allcoef, 
                top = text_grob("Coefficients throughout updates (Bayes method, periods of 50 events)", 
                                size = 14))
#100 events
plot100 <- lapply(all_coeff, function(covar) plot_function_bay(covar = covar, data = bayes_data_list_100))
plot100_allcoef <- ggarrange(plot100[[1]], plot100[[2]], plot100[[3]],
                             plot100[[4]], plot100[[5]], plot100[[6]],
                             plot100[[7]], plot100[[8]],
                             nrow = 4, ncol = 2, 
                             common.legend = TRUE, legend="bottom")
annotate_figure(plot100_allcoef, 
                top = text_grob("Coefficients throughout updates (Bayes method, periods of 100 events)", 
                                size = 14))
#150 events
plot150 <- lapply(all_coeff, function(covar) plot_function_bay(covar = covar, data = bayes_data_list_150))
plot150_allcoef <- ggarrange(plot150[[1]], plot150[[2]], plot150[[3]],
                             plot150[[4]], plot150[[5]], plot150[[6]],
                             plot150[[7]], plot150[[8]],
                             nrow = 4, ncol = 2, 
                             common.legend = TRUE, legend="bottom")
annotate_figure(plot150_allcoef, 
                top = text_grob("Coefficients throughout updates (Bayes method, periods of 150 events)", 
                                size = 14))
