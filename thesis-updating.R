#In this file:
# - model of first wave is build and inspected
# - updating periods of the second wave are created
# - refitting and recalibration is performed and validated
# - updated Bayesian models are evaluated
# - results from all methods are plot and numerically summarised
# - updated coefficients from th refitting method are plotted
# - first period of 50 events is illustared through a plot
# - visualization of individual predicted risk of refitting and recalibration are plotted

###############################################################################

#libraries needed:
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
#library(splines)
library(reshape2)
library(gt)
library(ggpubr)

#needed for validation functions:
library(riskRegression)
library(MatrixModels)
library(kableExtra)
library(dplyr)
library(tab)
library(timeROC)
library(geepack)

#############################################################################
#load data 
load("cleanscaledcut-cc-data.rda")

#load functions for validation of model
source("thesis-validation-functions.R") 

#load functions for no-updating, refitting and recalibration of the intercept
source("thesis-updating-functions.R")

#load results updated Bayesian models
load("bayes_data_list_50.rda") 
load("bayes_data_list_100.rda") 
load("bayes_data_list_150.rda") 

#############################################################################
# -------- First wave model -------------------------------------------------
#############################################################################

#covariates:
#all covariates to include in cox model
var <- c("age", "sex", "numb_comorb", "rr_min", 
         "sp_o2", "crp_mg_l", "creatinine_ser_mumol_l")
#formula
fmla <- formula(paste("Surv(survtime, status) ~ ", paste(var, collapse = "+")))
fmla

#cox model fitted on first wave data
model_w1 <- coxph(fmla, data=datasub_w1_cc, method="breslow", x=T, y=T)
summary(model_w1)
is(model_w1)

#plot of summary output model wave 1
library(tab)
tabcoxph(model_w1,
         columns = c("beta.ci", "hr.ci", "p"))

#validation model wave 1 on data wave 1
val_model_w1<- getValidation(model = model_w1, data = datasub_w1_cc)
val_model_w1

###############################################################################
# ------- Second wave --------------------------------------------------------
###############################################################################

#linear predictor for new data (wave 2) based on original model
datasub_w2_cc$lp_origmod <- predict(model_w1,
                                    newdata = datasub_w2_cc,
                                    type = "lp")
#needed for updates of recalibration method

################################################################################
# --------- Intervals/periods ----------------------------------------------------------
#creating intervals

startdate <- as.Date("2020-08-01") #start second wave
data <- datasub_w2_cc 

make_intervals <- function(data, startdate, n_events, n_interval, sliding_win, 
                           extraevents = 0) {
  
  # Arrange data by event date
  data <- arrange(data, date_event)
  enddate <- NULL
  
  # Find startdate and enddate for interval number n
  # Here I assume the new interval starts the day after the previous one ends
  for (i in 1:n_interval) { 
    data$event_count <- cumsum((data$admitted_date >= startdate) & 
                                 (data$status == 1))
    enddate <- data$date_event[data$event_count == n_events][1] 
    if (i != n_interval) startdate <- enddate + 1 
    if (extraevents !=0) enddate <- 
      data$date_event[data$event_count == n_events + extraevents][1] 
      #period of n_events+extraevents (no influence on startdate!)
  }
  
  # Select data between starting date (minus sliding window) and end date 
  data <- subset(data, admitted_date >= startdate - sliding_win 
                 & admitted_date <= enddate)
  
  #censor event after end date
  data$status <- ifelse(!is.na(data$date_event) & data$date_event > enddate,
                        0, data$status)
  #max survival time put to date of enddate
  data$survtime <- ifelse(!is.na(data$date_event) & data$date_event > enddate,
                          data$survtime - (data$date_event - enddate), data$survtime)
  return(data)
}

# -----------------------------------------------------------------------------
# 50 events:
#creating 15 (and not 16) periods, because validation is also on 15
# Intervals for model fitting: 0 week sliding window
df_int50_fit0 <- lapply(1:15, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 50, n_interval = n_int, sliding_win=0)
}
)
# Intervals for model fitting: 1 week sliding window
df_int50_fit1 <- lapply(1:15, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 50, n_interval = n_int, sliding_win=7)
}
)
# Intervals for model fitting: 2 weeks sliding window
df_int50_fit2 <- lapply(1:15, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 50, n_interval = n_int, sliding_win=14)
}
)
# Intervals for model fitting: 4 weeks sliding window
df_int50_fit4 <- lapply(1:15, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 50, n_interval = n_int, sliding_win=28)
}
)

# -----------------------------------------------------------------------------
# 100 events:
# Intervals for model fitting: 0 week sliding window 
df_int100_fit0 <- lapply(1:9, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 100, n_interval = n_int, sliding_win=0)
}
)
# Intervals for model fitting: 1 week sliding window
df_int100_fit1 <- lapply(1:9, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 100, n_interval = n_int, sliding_win=7)
}
)
# Intervals for model fitting: 2 weeks sliding window
df_int100_fit2 <- lapply(1:9, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 100, n_interval = n_int, sliding_win=14)
}
)
# Intervals for model fitting: 4 weeks sliding window
df_int100_fit4 <- lapply(1:9, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 100, n_interval = n_int, sliding_win=28)
}
)

#------------------------------------------------------------------------------
# 150 events:
# Intervals for model fitting: 0 week sliding window
df_int150_fit0 <- lapply(1:6, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 150, n_interval = n_int, sliding_win=0)
}
)
# Intervals for model fitting: 1 week sliding window
df_int150_fit1 <- lapply(1:6, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 150, n_interval = n_int, sliding_win=7)
}
)
# Intervals for model fitting: 2 weeks sliding window
df_int150_fit2 <- lapply(1:6, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 150, n_interval = n_int, sliding_win=14)
}
)
# Intervals for model fitting: 4 weeks sliding window
df_int150_fit4 <- lapply(1:6, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 150, n_interval = n_int, sliding_win=28)
}
)

# -------------------------------------------------------------------------------
# Intervals for prediction:
#no sliding window!
# 50 events -> 100 events (but startdate determined by 50 events periods)
df_int50_pred <- lapply(1:15, function(n_int) {
  #length is one period shorter
  make_intervals(datasub_w2_cc, startdate, n_events = 50, n_interval = n_int, 
                 sliding_win=0, extraevents = 50) 
  #validate on period of n_events+extraval=100 events
}
)
# 100 events
df_int100_pred <- lapply(1:9, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 100, n_interval = n_int, 
                 sliding_win=0)
}
)
# 150 events -> 100 events (but startdate determined by 150 events periods)
df_int150_pred <- lapply(1:6, function(n_int) {
  make_intervals(datasub_w2_cc, startdate, n_events = 150, n_interval = n_int, 
                 sliding_win=0, extraevents = -50) 
  #validate on period of n_events-extraevents= 100 events
}
)


################################################################################
# --------------- Applying methods --------------------------------------------
###############################################################################

#no updating
no_update_50 <- no_update(mod = model_w1, datapred = df_int50_pred)#50 events
no_update_50 <- no_update_50[-1] #validation for updates is on second interval
no_update_100 <- no_update(mod = model_w1, datapred = df_int100_pred) #100 events
no_update_100 <- no_update_100[-1] #validation for updates is on second interval
no_update_150 <- no_update(mod = model_w1, datapred = df_int150_pred) #100 events
no_update_150 <- no_update_150[-1] #validation for updates is on second interval
#will be used to compare performance

#-----------------------------------------------------------------------------

#Performing refitting method on periods of 50/100/150 events with varying 
#sliding windows, and validating on periods of 100 events (with dates depended 
#on the periods of 50/100/150 events)

#refitting: 50 events
refit_event50_sw1 <- refit(datafit = df_int50_fit1, datapred = df_int50_pred)
refit_event50_sw2 <- refit(datafit = df_int50_fit2, datapred = df_int50_pred)
refit_event50_sw4 <- refit(datafit = df_int50_fit4, datapred = df_int50_pred)
refit_event50_sw0 <- refit(datafit = df_int50_fit0, datapred = df_int50_pred)
#putting list of data together that we want to compare
refit_data_list_50 <- list(refit_event50_sw0,
                           refit_event50_sw1,
                           refit_event50_sw2,
                           refit_event50_sw4,
                           no_update_50)

#refitting: 100 events
refit_event100_sw1 <- refit(datafit = df_int100_fit1, datapred = df_int100_pred)
refit_event100_sw2 <- refit(datafit = df_int100_fit2, datapred = df_int100_pred)
refit_event100_sw4 <- refit(datafit = df_int100_fit4, datapred = df_int100_pred)
refit_event100_sw0 <- refit(datafit = df_int100_fit0, datapred = df_int100_pred)
refit_data_list_100 <- list(refit_event100_sw0,
                            refit_event100_sw1,
                            refit_event100_sw2,
                            refit_event100_sw4,
                            no_update_100)

#refitting: 150 events
refit_event150_sw1 <- refit(datafit = df_int150_fit1, datapred = df_int150_pred)
refit_event150_sw2 <- refit(datafit = df_int150_fit2, datapred = df_int150_pred)
refit_event150_sw4 <- refit(datafit = df_int150_fit4, datapred = df_int150_pred)
refit_event150_sw0 <- refit(datafit = df_int150_fit0, datapred = df_int150_pred)
refit_data_list_150 <- list(refit_event150_sw0,
                            refit_event150_sw1,
                            refit_event150_sw2,
                            refit_event150_sw4,
                            no_update_150)

#------------------------------------------------------------------------------

#Performing recalibration method on periods of 50/100/150 events with varying 
#sliding windows, and validating on periods of 100 events (with dates depended 
#on the periods of 50/100/150 events)


#recalibration of intercept with 50 events
recal_event50_sw1 <- recal_intercept(datafit = df_int50_fit1, datapred = df_int50_pred)
recal_event50_sw2 <- recal_intercept(datafit = df_int50_fit2, datapred = df_int50_pred)
recal_event50_sw4 <- recal_intercept(datafit = df_int50_fit4, datapred = df_int50_pred)
recal_event50_sw0 <- recal_intercept(datafit = df_int50_fit0, datapred = df_int50_pred)
recal_data_list_50 <- list(recal_event50_sw0,
                           recal_event50_sw1,
                           recal_event50_sw2,
                           recal_event50_sw4,
                           no_update_50)

#recalibration of intercept with 100 events
recal_event100_sw1 <- recal_intercept(datafit = df_int100_fit1, datapred = df_int100_pred)
recal_event100_sw2 <- recal_intercept(datafit = df_int100_fit2, datapred = df_int100_pred)
recal_event100_sw4 <- recal_intercept(datafit = df_int100_fit4, datapred = df_int100_pred)
recal_event100_sw0 <- recal_intercept(datafit = df_int100_fit0, datapred = df_int100_pred)
recal_data_list_100 <- list(recal_event100_sw0,
                            recal_event100_sw1,
                            recal_event100_sw2,
                            recal_event100_sw4,
                            no_update_100)

#recalibration of intercept with 150 events
recal_event150_sw1 <- recal_intercept(datafit = df_int150_fit1, datapred = df_int150_pred)
recal_event150_sw2 <- recal_intercept(datafit = df_int150_fit2, datapred = df_int150_pred)
recal_event150_sw4 <- recal_intercept(datafit = df_int150_fit4, datapred = df_int150_pred)
recal_event150_sw0 <- recal_intercept(datafit = df_int150_fit0, datapred = df_int150_pred)
recal_data_list_150 <- list(recal_event150_sw0,
                            recal_event150_sw1,
                            recal_event150_sw2,
                            recal_event150_sw4,
                            no_update_150)

#-------------------------------------------------------------------------------
# Bayesian update

# Predicting on all prediction interval except interval 1 (where we have no updates yet)
l50 <- length(df_int50_pred) - 1
Bayes_event50_sw0 <- lapply(1:l50, function(a) ValBayes(a, bayes_data_list_50[["sw0"]][-1], df_int50_pred[-1]))
Bayes_event50_sw1 <- lapply(1:l50, function(a) ValBayes(a, bayes_data_list_50[["sw1"]][-1], df_int50_pred[-1]))
Bayes_event50_sw2 <- lapply(1:l50, function(a) ValBayes(a, bayes_data_list_50[["sw2"]][-1], df_int50_pred[-1]))
Bayes_event50_sw4 <- lapply(1:l50, function(a) ValBayes(a, bayes_data_list_50[["sw4"]][-1], df_int50_pred[-1]))
Bayes_data_list_50 <- list(Bayes_event50_sw0,
                           Bayes_event50_sw1,
                           Bayes_event50_sw2,
                           Bayes_event50_sw4,
                           no_update_50)

l100 <- length(df_int100_pred) - 1
Bayes_event100_sw0 <- lapply(1:l100, function(a) ValBayes(a, bayes_data_list_100[["sw0"]][-1], df_int100_pred[-1]))
Bayes_event100_sw1 <- lapply(1:l100, function(a) ValBayes(a, bayes_data_list_100[["sw1"]][-1], df_int100_pred[-1]))
Bayes_event100_sw2 <- lapply(1:l100, function(a) ValBayes(a, bayes_data_list_100[["sw2"]][-1], df_int100_pred[-1]))
Bayes_event100_sw4 <- lapply(1:l100, function(a) ValBayes(a, bayes_data_list_100[["sw4"]][-1], df_int100_pred[-1]))
Bayes_data_list_100 <- list(Bayes_event100_sw0,
                            Bayes_event100_sw1,
                            Bayes_event100_sw2,
                            Bayes_event100_sw4,
                            no_update_100)

l150 <- length(df_int150_pred) - 1
Bayes_event150_sw0 <- lapply(1:l150, function(a) ValBayes(a, bayes_data_list_150[["sw0"]][-1], df_int150_pred[-1]))
Bayes_event150_sw1 <- lapply(1:l150, function(a) ValBayes(a, bayes_data_list_150[["sw1"]][-1], df_int150_pred[-1]))
Bayes_event150_sw2 <- lapply(1:l150, function(a) ValBayes(a, bayes_data_list_150[["sw2"]][-1], df_int150_pred[-1]))
Bayes_event150_sw4 <- lapply(1:l150, function(a) ValBayes(a, bayes_data_list_150[["sw4"]][-1], df_int150_pred[-1]))
Bayes_data_list_150 <- list(Bayes_event150_sw0,
                            Bayes_event150_sw1,
                            Bayes_event150_sw2,
                            Bayes_event150_sw4,
                            no_update_150)


#-------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#Some functions to plot the validation results:

#function that manipulates validation output such that it can be used for the 
#plots and other later use
reshape_validationdata <- function(list_data){
  df <- NULL
  #reshaping table with summary of validation measures
  for (i in 1:length(list_data)) {
    df[[i]] <- melt(data = list_data[[i]], level = 1) 
    df[[i]] <- reshape(data = df[[i]], idvar = c("Var1", "L1"), 
                       timevar = "Var2", direction = "wide")
    names(df[[i]]) <- c("measure", "update", "estimate", "CIlow", "CIhigh")
    df[[i]]$slidingwindow <- i-1 #assuming first element of list_data has no 
    #sliding window, second to 1, ...
  }
  df <- as.data.frame(do.call(rbind, df))
  df$update <- as.factor(df$update) #needed for plots
  
  df$slidingwindow <- case_when(
    df$slidingwindow == 4 ~ "no updating",
    df$slidingwindow == 3 ~ "4 weeks",
    df$slidingwindow == 2 ~ "2 weeks",
    df$slidingwindow == 1 ~ "1 week",
    TRUE ~ "0 weeks" #when slidingwindow == 0
  )
  df$slidingwindow <- as.factor(df$slidingwindow) #needed for plots
  return(df)
}

#function makes plot for each validation measure
plot_val_scale <- function(val_measure, list_data){
  df <- reshape_validationdata(list_data = list_data)
  plotrange<-NA
  xint <- NULL
  #scale of x-axis plot based on validation measure
  if(val_measure == "Scaled.Brier") plotrange <- c(0,0.5)
  else if (val_measure %in% c("Unos.c.index", "AUCt")) plotrange <- c(0.5,1)
  else if (val_measure %in% c("OE.ratio","Calibration.slope")) {
    plotrange <- c(0,2)
    xint <- 1
  }
  else{
    plotrange <- c(-1,1)
    xint <- 0
  }

  #creating plot
  plt <- ggplot(data = df[df$measure==val_measure,], 
                aes(x=estimate,y = update, 
                    col=slidingwindow)) + #comparing different sliding windows
    geom_point(position = position_dodge(width=0.4))+
    geom_linerange(aes(xmin= CIlow, xmax=CIhigh), 
                   position = position_dodge(width = 0.4))+
    geom_vline(xintercept = xint, col = "red", linetype = "dashed")+
    labs(xlab="estimate", ylab ="updates", 
         title = val_measure)+
    coord_flip(xlim = plotrange) #transpose plot (x-axis becomes y-axis) 
  return(plt)
}

#validation measure names
validation_names <- c("Scaled.Brier", "Unos.c.index", "AUCt", 
                     "OE.ratio","Calibration.intercept", "Calibration.slope")

#-----------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#numerical summary of performance measures (over all updates):

#function that computes numerical summary values for each validation measure 
#per sliding window
summary_performance <- function(val_measure, list_data){
  df <- reshape_validationdata(list_data = list_data) 
  df <- df[df$measure == val_measure,]
  sumval <- data.frame()
  if(val_measure %in%  c("Scaled.Brier", "Unos.c.index", "AUCt")) {
    sumval <- aggregate(df$estimate, list(df$slidingwindow), 
                        FUN = function(x) round(mean(x),3))
  }
  else if(val_measure == "OE.ratio"){
    sumval <- aggregate(df$estimate, list(df$slidingwindow), 
                        FUN = function(x) round(mean(abs(x-1)),3) ) 
  }
  else if(val_measure == "Calibration.slope"){
    sumval <- aggregate(df$estimate, list(df$slidingwindow), 
                        FUN = function(x) round(mean(abs(x-1)),3))
  }
  else { #calibration intercept
    sumval <- aggregate(df$estimate, list(df$slidingwindow), 
                        FUN = function(x) round(mean(abs(x-0)),3))
  }
  names(sumval) <- c("slidingwindow", "numeric_summary")
  return(sumval)
}

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#plotting results of refitting method

#plots:
#refitting:
#50 events
plots_refit_50events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = refit_data_list_50)
})
plot_refit_50 <- ggarrange(plots_refit_50events[[1]], plots_refit_50events[[2]],
                           plots_refit_50events[[3]], plots_refit_50events[[4]],
                           plots_refit_50events[[5]], plots_refit_50events[[6]],
                           nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_refit_50, 
                top = text_grob("Refitting method using periods of 50 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#100 events
plots_refit_100events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = refit_data_list_100)
})
plot_refit_100 <- ggarrange(plots_refit_100events[[1]], plots_refit_100events[[2]],
                            plots_refit_100events[[3]], plots_refit_100events[[4]],
                            plots_refit_100events[[5]], plots_refit_100events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_refit_100, 
                top = text_grob("Refitting method using periods of 100 events and varying sliding windows 
                                (validation on periods of 100 events)", 
                                size = 14))

#150 events:
plots_refit_150events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = refit_data_list_150)
})
plot_refit_150 <- ggarrange(plots_refit_150events[[1]], plots_refit_150events[[2]],
                            plots_refit_150events[[3]], plots_refit_150events[[4]],
                            plots_refit_150events[[5]], plots_refit_150events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_refit_150, 
                top = text_grob("Refitting method using periods of 150 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#numerical summaries:
#50 events 
summary_refit_50events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = refit_data_list_50)
})
names(summary_refit_50events) <- validation_names
#summary_refit_50events

#transform to table format
t <- lapply(summary_refit_50events, function(x) x[, 'numeric_summary', drop = FALSE])
t <- as.data.frame(do.call(cbind, t))
names(t) <- validation_names
t <- cbind(summary_refit_50events[[1]][1], t)
t
#transform to table latex:
gt::gt(t) |> 
  gt::as_latex()

#100 events
summary_refit_100events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = refit_data_list_100)
})
names(summary_refit_100events) <- validation_names
#summary_refit_100events

#transform to table format
t2 <- lapply(summary_refit_100events, function(x) x[, 'numeric_summary', drop = FALSE])
t2 <- as.data.frame(do.call(cbind, t2))
names(t2) <- validation_names
t2 <- cbind(summary_refit_100events[[1]][1], t2)
t2
gt::gt(t2) |> 
  gt::as_latex()

#150 events
summary_refit_150events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = refit_data_list_150)
})
names(summary_refit_150events) <- validation_names
#summary_refit_150events

#transform to table format
t3 <- lapply(summary_refit_150events, function(x) x[, 'numeric_summary', drop = FALSE])
t3 <- as.data.frame(do.call(cbind, t3))
names(t3) <- validation_names
t3 <- cbind(summary_refit_150events[[1]][1], t3)
t3
gt::gt(t3) |> 
  gt::as_latex()

#just for simple overview in console
tot_t <- list(t,t2,t3)
names(tot_t) <- c("50 events", "100 events", "150 events")
tot_t

#------------------------------------------------------------------------------
#plotting results of recalibration method

#plots:
#50 events
plots_recal_50events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = recal_data_list_50)
})
plot_recal_50 <- ggarrange(plots_recal_50events[[1]], plots_recal_50events[[2]],
                           plots_recal_50events[[3]], plots_recal_50events[[4]],
                           plots_recal_50events[[5]], plots_recal_50events[[6]],
                           nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_recal_50, 
                top = text_grob("Recalibration of intercept method using periods of 50 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#100 events
plots_recal_100events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = recal_data_list_100)
})
plot_recal_100 <- ggarrange(plots_recal_100events[[1]], plots_recal_100events[[2]],
                            plots_recal_100events[[3]], plots_recal_100events[[4]],
                            plots_recal_100events[[5]], plots_recal_100events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_recal_100, 
                top = text_grob("Recalibration of intercept method using periods of 100 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#150 events
plots_recal_150events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = recal_data_list_150)
})
plot_recal_150 <- ggarrange(plots_recal_150events[[1]], plots_recal_150events[[2]],
                            plots_recal_150events[[3]], plots_recal_150events[[4]],
                            plots_recal_150events[[5]], plots_recal_150events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_recal_150, 
                top = text_grob("Recalibration of intercept method using periods of 150 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#numerical summary:
#50 events
summary_recal_50events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = recal_data_list_50)
})
names(summary_recal_50events) <- validation_names
#transform to table format
t4 <- lapply(summary_recal_50events, function(x) x[, 'numeric_summary', drop = FALSE])
t4 <- as.data.frame(do.call(cbind, t4))
names(t4) <- validation_names
t4 <- cbind(summary_recal_50events[[1]][1], t4)
t4
gt::gt(t4) |> 
  gt::as_latex()

#100 events
summary_recal_100events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = recal_data_list_100)
})
names(summary_recal_100events) <- validation_names
#transform to table format
t5 <- lapply(summary_recal_100events, function(x) x[, 'numeric_summary', drop = FALSE])
t5 <- as.data.frame(do.call(cbind, t5))
names(t5) <- validation_names
t5 <- cbind(summary_recal_100events[[1]][1], t5)
t5
gt::gt(t5) |> 
  gt::as_latex()

#150 events
summary_recal_150events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = recal_data_list_150)
})
names(summary_recal_150events) <- validation_names
#transform to table format
t6 <- lapply(summary_recal_150events, function(x) x[, 'numeric_summary', drop = FALSE])
t6 <- as.data.frame(do.call(cbind, t6))
names(t6) <- validation_names
t6 <- cbind(summary_recal_150events[[1]][1], t6)
t6
gt::gt(t6) |> 
  gt::as_latex()

#just for simple overview in console
tot_t2 <- list(t4,t5,t6)
names(tot_t2) <- c("50 events", "100 events", "150 events")
tot_t2

#-------------------------------------------------------------------------------
#plotting results of Bayesian method

#plots:
#50 events
plots_bayes_50events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = Bayes_data_list_50)
})
plot_bayes_50 <- ggarrange(plots_bayes_50events[[1]], plots_bayes_50events[[2]],
                            plots_bayes_50events[[3]], plots_bayes_50events[[4]],
                            plots_bayes_50events[[5]], plots_bayes_50events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_bayes_50, 
                top = text_grob("Bayesian updating method using periods of 50 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#100 events
plots_bayes_100events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = Bayes_data_list_100)
})
plot_bayes_100 <- ggarrange(plots_bayes_100events[[1]], plots_bayes_100events[[2]],
                            plots_bayes_100events[[3]], plots_bayes_100events[[4]],
                            plots_bayes_100events[[5]], plots_bayes_100events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_bayes_100, 
                top = text_grob("Bayesian updating method using periods of 100 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#150 events
plots_bayes_150events <- lapply(validation_names, function(val){
  plot_val_scale(val_measure = val, list_data = Bayes_data_list_150)
})
plot_bayes_150 <- ggarrange(plots_bayes_150events[[1]], plots_bayes_150events[[2]],
                            plots_bayes_150events[[3]], plots_bayes_150events[[4]],
                            plots_bayes_150events[[5]], plots_bayes_150events[[6]],
                            nrow = 3, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plot_bayes_150, 
                top = text_grob("Bayesian updating method using periods of 150 events and varying sliding windows
                                (validation on periods of 100 events)", 
                                size = 14))

#numerical summary:
#50 events
summary_bayes_50events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = Bayes_data_list_50)
})
names(summary_bayes_50events) <- validation_names
#transform to table format
t7 <- lapply(summary_bayes_50events, function(x) x[, 'numeric_summary', drop = FALSE])
t7 <- as.data.frame(do.call(cbind, t7))
names(t7) <- validation_names
t7 <- cbind(summary_bayes_50events[[1]][1], t7)
t7
gt::gt(t7) |> 
  gt::as_latex()

#100 events
summary_bayes_100events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = Bayes_data_list_100)
})
names(summary_bayes_100events) <- validation_names
#transform to table format
t8 <- lapply(summary_bayes_100events, function(x) x[, 'numeric_summary', drop = FALSE])
t8 <- as.data.frame(do.call(cbind, t8))
names(t8) <- validation_names
t8 <- cbind(summary_bayes_100events[[1]][1], t8)
t8
gt::gt(t8) |> 
  gt::as_latex()

#150 events
summary_bayes_150events <- lapply(validation_names, function(val){
  summary_performance(val_measure = val, list_data = Bayes_data_list_150)
})
names(summary_bayes_150events) <- validation_names
#transform to table format
t9 <- lapply(summary_bayes_150events, function(x) x[, 'numeric_summary', drop = FALSE])
t9 <- as.data.frame(do.call(cbind, t9))
names(t9) <- validation_names
t9 <- cbind(summary_bayes_150events[[1]][1], t9)
t9
gt::gt(t9) |> 
  gt::as_latex()

#just for simple overview in console
tot_t3 <- list(t7,t8,t9)
names(tot_t3) <- c("50 events", "100 events", "150 events")
tot_t3


###############################################################################
# Additional information about updating periods
###############################################################################

#Closer look at first few  periods
max(df_int150_fit0[[1]]$admitted_date, na.rm = T)#-min(df_int50_pred[[1]]$admitted_date) 

max(df_int100_pred[[1]]$admitted_date, na.rm = T)

nrow(df_int50_fit0[[1]])
nrow(df_int50_pred[[2]])
nrow(df_int50_pred[[3]])
nrow(df_int50_pred[[4]])
nrow(df_int50_pred[[5]])


sum(df_int50_fit4[[15]]$status==1)
#mean number of events
#updating every 50 events:
#no sliding window
mean(c(50,54,50,52,51,51,51,51, 52, 53,53,52,55,50,51))
#1 week
mean(c(50,73,89,96, 89,75,83,84,87,78,74,75,77,79,64))
#2 weeks
mean(c(50,94,122,139,128,108,110,120,128,107,92,94,95,109,79))
#4 weeks
mean(c(50,105,162,212,220,188,165,179,198,173,140,139,153,111))
#updating every 100 events:
sum(df_int100_fit4[[9]]$status==1)
#no sliding window
mean(c(103,109,105,100,102,104,103,103,100))
#1 week
mean(c(103,154,144,122,136,130,125,125,114))
#2 weeks
mean(c(103,182,190,147,172,167,144,143,133))
#4 weeks
mean(c(103,219,278,223,231,245,190,187,184))
#updating every 150 events:
sum(df_int150_fit4[[6]]$status==1)
#no sliding window
mean(c(157,150,151,151,155,150))
#1 week
mean(c(157,196,176,190,173,180))
#2 weeks
mean(c(157,241,214,220,201,200))
#4 weeks
mean(c(157,294,291,294,249,244))

l <- lapply(df_int50_fit4, nrow)
l<- unlist(l)
l
mean(l)
summary(l)

################################################################################
# Plots of coefficients throughout time for refitting methods
################################################################################
#change updating functions to return model instead of validation values!!!!!

coeff_ci <- function(covar, data){
  list_coef <- numeric()
  list_cilow <- numeric()
  list_cihigh <- numeric()
  list_covar <- c("age", "sexM", "numb_comorb", "rr_min", 
                  "sp_o2", "crp_mg_l", "creatinine_ser_mumol_l")
  index_covar <- which(list_covar == covar)
  for (i in 1:length(data)) {
    list_coef <- cbind(list_coef, data[[i]]$coefficients[covar]) #coefficients
    list_cilow <- cbind(list_cilow, data[[i]]$coefficients[covar] 
                        - 1.96*summary(data[[i]])$coefficients[covar,"se(coef)"])
    list_cihigh <- cbind(list_cihigh, data[[i]]$coefficients[covar] 
                         + 1.96*summary(data[[i]])$coefficients[covar,"se(coef)"])
  }
  res <- rbind(1:length(data),list_coef, list_cilow, list_cihigh)
  rownames(res) <- c("update", covar, "CIlow", "CIhigh")
  res <- t(res)
  return(res)
}

reshape_validationdata2 <- function(list_data, covar){
  df <- NULL
  for (i in 1:length(list_data)) {
    df[[i]] <- melt(data = list_data[[i]], level = 1) 
    df[[i]] <- reshape(data = df[[i]], idvar = c("Var1"), 
                       timevar = "Var2", direction = "wide")
    names(df[[i]]) <- c("update", "update", covar, "CIlow", "CIhigh")
    df[[i]]$slidingwindow <- i-1 #assuming first element of list_data has no 
    #sliding window, second to 1, ...
  }
  df <- as.data.frame(do.call(rbind, df))
  df$update <- as.factor(df$update) #needed for plots
  
  df$slidingwindow <- case_when(
    df$slidingwindow == 4 ~ "no updating",
    df$slidingwindow == 3 ~ "4 weeks",
    df$slidingwindow == 2 ~ "2 weeks",
    df$slidingwindow == 1 ~ "1 week",
    TRUE ~ "0 weeks" #when slidingwindow == 0
  )
  df$slidingwindow <- as.factor(df$slidingwindow) #needed for plots
  return(df)
}

plot_function <- function(covar, res){
  res <- as.data.frame(res)
  plt <- ggplot(data = res,
                aes(y= res[,covar], x= update, col = slidingwindow))+
    geom_point(position = position_dodge(width=0.4))+
    geom_linerange(aes(ymin= CIlow, ymax=CIhigh),
                   position = position_dodge(width = 0.4))+
    labs(y= "coefficient", x ="updates",
        title = covar)#+
  return(plt)
}

make_plot_coeff <- function(data, covar){
  t1 <- lapply(data, function(dat) coeff_ci(covar, data = dat))
  t2<- reshape_validationdata2(t1, covar)
  t2<-t2[,-1]
  plot_test <- plot_function(covar, t2)
  return(plot_test)
}

var1 <- c("age", "sexM", "numb_comorb", "rr_min", 
          "sp_o2", "crp_mg_l", "creatinine_ser_mumol_l")

refit_data_list_100 <- list(refit_event100_sw0,
                            refit_event100_sw1,
                            refit_event100_sw2,
                            refit_event100_sw4,
                            list(model_w1, model_w1, model_w1, model_w1, 
                                 model_w1, model_w1, model_w1, model_w1))

refit_100_coeff <- lapply(var1, function(var2) 
  make_plot_coeff(data = refit_data_list_100, covar = var2))
refit_100_coeff[[2]]

plots_refit_100_coeff<- ggarrange(refit_100_coeff[[1]], refit_100_coeff[[2]],
                                       refit_100_coeff[[3]], refit_100_coeff[[4]],
                                       refit_100_coeff[[5]], refit_100_coeff[[6]],
                                       refit_100_coeff[[7]],
                            nrow = 4, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plots_refit_100_coeff, 
                top = text_grob("Coefficients throughout updates (refitting method, periods of 100 events)", 
                                size = 14))


refit_data_list_150 <- list(refit_event150_sw0,
                            refit_event150_sw1,
                            refit_event150_sw2,
                            refit_event150_sw4,
                            list(model_w1, model_w1, model_w1, model_w1, model_w1))

refit_data_list_50 <- list(refit_event50_sw0,
                            refit_event50_sw1,
                            refit_event50_sw2,
                            refit_event50_sw4,
                            list(model_w1, model_w1, model_w1, model_w1, 
                                 model_w1,model_w1, model_w1, model_w1, model_w1, model_w1,model_w1, model_w1, model_w1, model_w1))
refit_150_coeff <- lapply(var1, function(var2) 
  make_plot_coeff(data = refit_data_list_150, covar = var2))
plots_refit_150_coeff<- ggarrange(refit_150_coeff[[1]], refit_150_coeff[[2]],
                                  refit_150_coeff[[3]], refit_150_coeff[[4]],
                                  refit_150_coeff[[5]], refit_150_coeff[[6]],
                                  refit_150_coeff[[7]],
                                  nrow = 4, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plots_refit_150_coeff, 
                top = text_grob("Coefficients throughout updates (refitting method, periods of 150 events)", 
                                size = 14))

refit_50_coeff <- lapply(var1, function(var2) 
  make_plot_coeff(data = refit_data_list_50, covar = var2))
plots_refit_50_coeff<- ggarrange(refit_50_coeff[[1]], refit_50_coeff[[2]],
                                 refit_50_coeff[[3]], refit_50_coeff[[4]],
                                 refit_50_coeff[[5]], refit_50_coeff[[6]],
                                 refit_50_coeff[[7]],
                                  nrow = 4, ncol = 2,  common.legend = TRUE, legend="bottom")
annotate_figure(plots_refit_50_coeff, 
                top = text_grob("Coefficients throughout updates (refitting method, periods of 50 events)", 
                                size = 14))

################################################################################
# More plots
################################################################################
#visualization events throughout time:

#only individuals with event
df <- datasub_w2_cc[!is.na(datasub_w2_cc$date_event),]
#order by admission date
df <- df[order(df$admitted_date),]

#last individual admitted in interval
last_inv_admitted <- max(which(df$admitted_date==df$date_event[df$count==50]))
last_inv_admitted

#individual with event 50
event50_id <- which(df$count==50)
event50_id
#admission date of 50th event interval 1
admission_event50 <- df$admitted_date[df$count==50]
admission_event50
#date of event of 50th event in interval 1
date_event_first50 <- df$date_event[df$count==50]
date_event_first50

#first 50 events
ggplot(df[1:70,])+geom_point(aes(y=as.numeric(row.names(df))[1:70], x = admitted_date))+
  geom_point(aes(y =as.numeric(row.names(df))[1:70], x= date_event)) +
  geom_segment(aes(x = admitted_date, y = as.numeric(row.names(df))[1:70], xend = date_event, yend = as.numeric(row.names(df))[1:70])) +
  geom_hline(yintercept = last_inv_admitted, size=1) + #all individuals included before this line
  geom_hline(yintercept = event50_id, linetype = "dashed", size=1) + #50th event
  geom_vline(xintercept = date_event_first50, size=1) + #all events censored above this line
  labs(title = "Illustration first period of 50 events: admission and event dates of individuals with an event",
       y = "inidividual",
       x = "date")+
  annotate("rect", xmin = date_event_first50, xmax = date_event_first50+16.5, ymin = 0, ymax =last_inv_admitted,
           alpha = .1,fill = "blue")

#-------------------------------------------------------------------------------
#visualization of individual predicted risk of refitting and recalibration
#updating every 50 events with nos sliding window
#first update, validated on next period of 100 events

#change updating functions to return updated model only!!!

#predicted risk by model on next period of 100 events REFITTING
predrisk <-  riskRegression::predictRisk(refit_event50_sw0[[1]], newdata = df_int50_pred[[2]],
                                         times = 28)

#predicted risk by model on next period of 100 events RECALIBRATION
model <- remove_offset(recal_event50_sw0[[1]])
predrisk2 <- riskRegression::predictRisk(model, newdata = df_int50_pred[[2]],
                                         times = 28)
#plot of how predictions differ between methods
plot(predrisk-predrisk2, main="The difference in predicted risk of individuals between refitting and recalibration 
     after update first period of 50 events (no sliding window) 
     validated on individuals from next period of 100 events",
     ylab = "difference of predicted risk between refitting and recalibration",
     xlab = "individual")
abline(h=0)
#positive: predicted risk refitting higher than predicted risk recalibration
#negative: predicted risk recalibration higher than predicted risk refitting
