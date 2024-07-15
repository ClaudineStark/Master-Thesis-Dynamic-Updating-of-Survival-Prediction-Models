################################################################################
# In this file:
# Functions for model validation are provided.
#
# - Overall performance:
#   - Brier score
# - Discrimination:
#   - Uno's c-index
#   - AUCt
# - Calibration
#   - observed-expected ratio
#   - calibration plot
#   - calibration intercept and slope

# Code based on
# - https://github.com/danielegiardiello/Prediction_performance_survival/blob/main/01_predsurv_minimal.R
# - https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC_minimal.R


#------------- Validation functions ------------------------------------------

# ------- Overall performance -------------------------------------------------
# --------- Brier score -------------------------------------------------------

#function that computes Brier score given the model and data
getBrier <- function(model, data){
  #reformulate if model has offset
  if(any(is.numeric(model))) {
    bs <- riskRegression::Score(list("numeric" = model), #model
                                formula = Surv(survtime, status) ~ 1,
                                data = data, #new data, only complete cases
                                conf.int = TRUE,
                                times = 28, #prediction horizon
                                cens.model = "km", #using KM to estimate observed
                                metrics = "brier",
                                summary = "ipa",
                                plots="cal")
  }
  if(any(grepl("coxph", is(model)))) {
  if(!is.null(model$offset)) model <- remove_offset(model)
  bs <- riskRegression::Score(list("cox" = model), #model
                              formula = Surv(survtime, status) ~ 1,
                              data = data, #new data, only complete cases
                              conf.int = TRUE,
                              times = 28, #prediction horizon
                              cens.model = "km", #using KM to estimate observed
                              metrics = "brier",
                              summary = "ipa",
                              plots="cal")
  }
  brier <- c(as.numeric(bs$Brier$score[2,3]), #estimate Brier score
             as.numeric(bs$Brier$score[2,5]), #lower
             as.numeric(bs$Brier$score[2,6])) #upper
  #estimated scale Brier score (IPA) (using formula)
  ipa <- c((1-(as.numeric(bs$Brier$score[2,3])/as.numeric(bs$Brier$score[1,3]))), #estimate
           (1-(as.numeric(bs$Brier$score[2,5])/as.numeric(bs$Brier$score[1,5]))), #lower
           (1-(as.numeric(bs$Brier$score[2,6])/as.numeric(bs$Brier$score[1,6])))) #upper
  return(rbind(brier, ipa))
}

# --------- Discrimination ----------------------------------------------------
# ------ Uno's C index --------------------------------------------------------

#function that computes Uno's c-index with 95% CI
getCindex <- function(model, data){
  #add linear predictor (bX)
  if(any(is.numeric(model))) data$lp <- model
  if(any(grepl("coxph", is(model)))) {
  data$lp <- predict(model, newdata = data, type = "lp", reference = "zero")
  if (!is.null(model$offset)) data$lp <- data$lp_origmod #use original lp
  }
  
  Uno_C <- concordance(Surv(survtime, status) ~ lp, 
                       data = data, 
                       reverse = TRUE,
                       timewt = "n/G2") #weights
  alpha <- 0.05
  res_C <- matrix(
    c(Uno_C$concordance,
      Uno_C$concordance - 
        qnorm(1 - alpha/2) * sqrt(Uno_C$var),
      Uno_C$concordance + 
        qnorm(1 - alpha/2) * sqrt(Uno_C$var)), 
    nrow = 1, ncol = 3, byrow = T,
    dimnames = list(c("Uno C"),
                    c("Estimate", "2.5 %", "97.5 %")))
  return(res_C)
}

# ------- AUCt ----------------------------------------------------------------

#function that computes time-dependent ROC curve estimation with 95% CI
getAUCt <- function(model, data){
  #add linear predictor (bX)
  if(any(is.numeric(model))) data$lp <- model
  if(any(grepl("coxph", is(model)))) {
    data$lp <- predict(model, newdata = data, type = "lp")
    if (!is.null(model$offset)) data$lp <- data$lp_origmod #use original lp
  }

  #IPCW estimation of time-dependent ROC curve
  Uno <- timeROC::timeROC(
    T = data$survtime,  #vector of event-times
    delta = data$status, #vector of event indicator
    marker = data$lp,  #marker values for which we want to compute the time-dependent ROC curves
    cause = 1, #value that indicates non-censored values
    weighting = "marginal", # > KM estimator for censoring distribution
    times = 28, #time horizon
    ROC = TRUE,
    iid = TRUE #want to compute the iid-representation of the area under 
    #time-dependent ROC curve estimator
  )
  alpha <- 0.05
  #95% CI table of estimate
  Uno_AUC_res <- c(
    "time-dependent Uno AUC" = unname(Uno$AUC[2]),
    "2.5 %" = unname(Uno$AUC["t=28"] -
                       qnorm(1 - alpha / 2) * Uno$inference$vect_sd_1["t=28"]),
    "97. 5 %" = unname(Uno$AUC["t=28"] +
                         qnorm(1 - alpha / 2) * Uno$inference$vect_sd_1["t=28"])
  )
  return(Uno_AUC_res)
}

# ---- Calibration ------------------------------------------------------------
# ----- Observed-Expected ratio -----------------------------------------------

#function that computes observed-expected ratio with 95% CI
getOEratio <- function(model, data){
  time_horizon <- 28
  #observed risk:
  observed <- summary(survfit(Surv(survtime, status)~1, data = data), 
                      times = time_horizon) #KM survival curve
  observed_t <- 1 - observed$surv #observed risk (=1-surival)
  
  #predicted risk:
  if(any(is.numeric(model))) data$predrisk <- model
  if(any(grepl("coxph", is(model)))) {
    if(!is.null(model$offset)) model <- remove_offset(model)
    data$predrisk <-  riskRegression::predictRisk(model, newdata = data,
                                                  times = time_horizon) 
  }

  expected_t <- mean(data$predrisk) #expected
  
  #ratio:
  OE_t <- observed_t / expected_t
  alpha <- .05
  OE_summary <- c(
    "OE" = OE_t,
    "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / observed$n.event)),
    "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / observed$n.event))
  )
  return(OE_summary)
}

# --- Calibration plot --------------------------------------------------------

#function that computes the calibration plot
getCalplot <- function(model, data){
  if(any(is.numeric(model))) {
    cal <-riskRegression::Score(list("numeric" = model), #model
                                formula = Surv(survtime, status) ~ 1,
                                data = data, #new data
                                conf.int = TRUE,
                                times = 28, #prediction horizon
                                cens.model = "km", #using KM to estimate observed
                                metrics = c("auc", "brier"),
                                summary = c("ipa"),
                                cause = 1,
                                plots="calibration")
  }
  if(any(grepl("coxph", is(model)))) {
    if(!is.null(model$offset)) model <- remove_offset(model)
    cal <-riskRegression::Score(list("cox" = model), #model
                              formula = Surv(survtime, status) ~ 1,
                              data = data, #new data
                              conf.int = TRUE,
                              times = 28, #prediction horizon
                              cens.model = "km", #using KM to estimate observed
                              metrics = c("auc", "brier"),
                              summary = c("ipa"),
                              cause = 1,
                              plots="calibration")
  }
  return(plotCalibration(cal, 
                         method = "nn",
                         cens.method = "pseudo", #pseudo observations approach
                         round = FALSE, #keep unique estimates
                         legend = TRUE,
                         xlab = "Estimated risks",
                         ylab = "Observed outcome proportions", #censoring estimates with KM
                         plot = TRUE)
  )
}

# --- Calibration slope and intercept -------------------------------------------------------

#function that computes the calibration slope and intercept with 95% CI
#based on: https://github.com/survival-lumc/ValidationCompRisks/blob/main/Prediction_CSC_minimal.R
getCalslopeint <- function(model, data){
  if(any(is.numeric(model))) {
    cal <-riskRegression::Score(list("numeric" = model), #model
                                formula = Surv(survtime, status) ~ 1,
                                data = data, #new data
                                conf.int = TRUE,
                                times = 28, #prediction horizon
                                cens.model = "km", #using KM to estimate observed
                                metrics = c("auc", "brier"),
                                summary = c("ipa"),
                                cause = 1,
                                plots="calibration")
  }
  if(any(grepl("coxph", is(model)))) {
  if(!is.null(model$offset)) model <- remove_offset(model)
  cal <-riskRegression::Score(list("cox" = model), #model
                              formula = Surv(survtime, status) ~ 1,
                              data = data, #new data, only complete cases
                              conf.int = TRUE,
                              times = 28, #prediction horizon
                              cens.model = "km", #using KM to estimate observed
                              #metrics = c("auc", "brier"),
                              summary = c("ipa"),
                              cause = 1,
                              plots="calibration")
  }
  pseudos <- data.frame(cal$Calibration$plotframe)
  #sometimes the pseudovalue is 1 and that creates issues:
  pseudos$risk[pseudos$risk == 1] <- 1 - 10^-6
  pseudos$cll_pred <- log(-log(1-pseudos$risk)) #change scale
  
  #Fit model for calibration intercept
  fit_cal_int <- geese(
    pseudovalue ~ offset(cll_pred),
    data = pseudos,
    id = ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog", #back to original scale
    corstr = "independence", #correlation structure
    jack = TRUE,
    na.action = na.omit)
  
  #Fit model for calibration slope
  fit_cal_slope <- geese(
    pseudovalue ~ offset(cll_pred) + cll_pred,
    data = pseudos,
    id = ID,
    scale.fix = TRUE,
    family = gaussian,
    mean.link = "cloglog",
    corstr = "independence",
    jack = TRUE,
    na.action = na.omit)
  
  #Value and CI
  sum1 <- with(
    summary(fit_cal_slope)$mean["cll_pred", ],
    c("slope" = 1 + estimate,
      "2.5 %" = 1 + (estimate - qnorm(0.975) * san.se),
      "97.5 %" = 1 + (estimate + qnorm(0.975) * san.se)))
  #Value and CI
  sum2 <- with(
    summary(fit_cal_int)$mean,
    c("intercept" = estimate,
      "2.5 %" = estimate - qnorm(0.975) * san.se,
      "97.5 %" = estimate + qnorm(0.975) * san.se))
  output <- rbind(sum2, sum1)
  colnames(output)[1] <- "estimate"
  rownames(output) <- c("intercept", "slope")
  return(output)
}

#------------------------------------------------------------------------------
# -------- Function that returns information of all validation functions -------

getValidation <- function(model, data){
  b1 <- getBrier(model = model, data = data) #Brier+IPA score
  b1_t <- t(b1)
  cindex1 <- getCindex(model = model, data = data) #C-index Uno
  auct1 <- getAUCt(model = model, data = data) #AUCt
  oe1 <- getOEratio(model = model, data = data) #OE ratio
  #calpl1 <- getCalplot(model = model, data = data) #calibration plot
  #calibration slope and intercept:
  slopeintercept1 <- getCalslopeint(model = model, data = data)
  slopeintercept1_t <- t(slopeintercept1)

  #summary of numerical estimates (with CI)
  summary_model <- data.frame("Brier score" = round(b1_t[,1], digits = 3),
                              "Scaled Brier" = round(b1_t[,2], 3),
                              "Unos c-index" = round(cindex1[1:3], 3),
                              "AUCt" = round(auct1, 3),
                              "OE ratio" = round(oe1, 3),
                              "Calibration intercept" = round(slopeintercept1_t[,1], 3),
                              "Calibration slope" = round(slopeintercept1_t[,2], 3))
  rownames(summary_model) <- c("estimate", "2.5%", "97.5%")
  summary_model_t <- t(summary_model)
  return(summary_model_t)
}

#did not end up using calibration plot for results
