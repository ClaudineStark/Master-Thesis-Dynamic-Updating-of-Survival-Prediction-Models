#This file contains the functions for the updating methods on a period
#and validation on next period

###############################################################################
#--------------- Updating + validation -----------------------------------------
################################################################################

# ----------------------------------------------------------------------------
#------ No updating ----------------------------------------------------------
# ----------------------------------------------------------------------------

#function that validates original model M0 on intervals (no updating)
no_update <- function(mod, datapred){
  summary_model_noupdate <- lapply(datapred, function(data){
    getValidation(model = mod, data = data)
  })
  return(summary_model_noupdate)
  #return(summary(mod))
}

# -----------------------------------------------------------------------------
# ----- Refitting -------------------------------------------------------------
# -----------------------------------------------------------------------------

#function that performs refitting and validation
refit <- function(datafit, datapred){
  updateRefit <- NULL #will contain updated models
  summary_model_updateRefit <- NULL
  #summary_model_noupdateprevious <- NULL

  #loop that performs refitting on intervals
  # and validates on next interval
  for (i in 1:(length(datafit)-1)) {
    #refit model (using same formula M0) on i-th interval
    updateRefit[[i]] <- coxph(fmla, 
                              data = datafit[[i]], #using i-th interval
                              method = "breslow", x=T, y=T) 
    #validation of updated model on i+1 interval (no sliding window)
    summary_model_updateRefit[[i]] <- getValidation(model = updateRefit[[i]],
                                                    data = datapred[[i+1]])
    #comparison to no update previous step
    #summary_model_noupdateprevious[[i]] <- getValidation(model = updateRefit[[i]],
    #                                                data = datapred[[i+2]])
  }
  return(summary_model_updateRefit)
  #return(updateRefit) #for plotting coefficient through updates
}

#------------------------------------------------------------------------------
# ----- Recalibration of intercept ---------------------------------------------
# -----------------------------------------------------------------------------

#function that performs recalibration of intercept and validation
recal_intercept <- function(datafit, datapred){
  updateRecal <- NULL #will conatain updated models
  summary_model_updateRecal <- NULL

  #loop that performs recalibration of intercept updates on intervals 
  # and validates on next intervals
  for (i in 1:(length(datafit)-1)) {
    #i-th update using data i-th interval
    updateRecal[[i]] <- coxph(Surv(survtime, status) ~ offset(lp_origmod),
                              data = datafit[[i]], #using i-th interval
                              method="breslow", x=T, y=T)
    
    #validation of updated model on i+1 interval (no sliding window)
    summary_model_updateRecal[[i]] <- getValidation(model = updateRecal[[i]],
                                                    data = datapred[[i+1]])
  }
  #return(updateRecal) #for retrieving individuals risk 
  return(summary_model_updateRecal)
}

#function that makes a model with only offset into a "normal" model 
remove_offset <- function(model) {
  lp <- model$linear.predictors
  model$coefficients <- c("lp_origmod" = 1)
  model$means <- c("lp_origmod" = mean(lp)) # mean is 0, looks like our old predictions were centered (need to watch out)
  terms <- formula("Surv(survtime, status) ~ lp_origmod")
  attributes(terms) <- attributes(model$terms)
  attr(terms, "variables") <- parse(text = "list(Surv(survtime, status), lp_origmod)")[[1]]
  attr(terms, "factors") <- matrix(c(0L, 1L), nrow = 2, ncol = 1, 
                                   dimnames = list(c("Surv(survtime, status)", "lp_origmod"),
                                                   c("lp_origmod")))
  attr(terms, "term.labels") <- c("lp_origmod")
  attr(terms, "order") <- c(1L)
  attr(terms, "offset") <- NULL
  attr(terms, "predvars") <- parse(text = "list(Surv(survtime, status), lp_origmod)")[[1]]
  attr(attr(terms, "dataClasses"), "names") <- c("Surv(survtime, status)", "lp_origmod")
  model$terms <- terms
  model$assign <- list("lp_origmod" = 1L)
  x <- matrix(c(lp), dimnames = list(attr(model$x, "dimnames")[[1]], c("lp_origmod")))
  attr(x, "assign") = c(1L)
  model$x <- x
  model$formula <- formula("Surv(survtime, status) ~ lp_origmod", env = globalenv())
  model$offset <- NULL
  model$call <- parse(text = 'coxph(formula = Surv(survtime, status) ~ lp_origmod, 
                      data = datafit[[i]], x = T, y = T, method = "breslow")')[[1]]
  model$loglik <- coxph(model$y ~ 1,  method = "breslow")$loglik
  return(model)
}

#------------------------------------------------------------------------------
# ---------- Bayesian updating ------------------------------------------------
#------------------------------------------------------------------------------

# Predictions of the Bayesian update

# This function gives the predicted risk on the new data at a given time
predBayes <- function (bayesModel, newdata, origModel, time) {
  newdata <- data.frame(newdata)
  if(length(time)!=1) "Function implemented for one time-point only"
  newdata$sexM <- if_else(newdata$sex == "M", 1, 0)
  coef_names <- names(origModel$coefficients)
  coeffs <- bayesModel[grep("beta",rownames(bayesModel)), "mean"]
  names(coeffs) <- coef_names
  lambda <- bayesModel[grep("lambda",rownames(bayesModel)), "mean"]
  lp <- rowSums(t(t(newdata[, coef_names])*coeffs))
  predRisk <- 1- exp(-exp(lambda)*time*exp(lp))
  
}

#This is just for the next step: we want to use predBayes on all updated models
#(not the first one, which is the original model)
ValBayes <- function(a, BayData, preddata) {
  if (any(is(BayData[[a]]) == "coxph")) pred <- BayData[[a]]
  else pred <- predBayes(BayData[[a]], preddata[[a]], model_w1, 28)
  
  getValidation(pred,preddata[[a]])
}

