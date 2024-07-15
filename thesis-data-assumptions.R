##############################################################################
# In this script:
# - KM curve is plotted
# - centering and rescaling of covariates is investigated
# - extreme values are cut off 
# - reverse KM is plotted
# - proportional hazard assumption is checked
# - linearity assumption for continuous covariates is checked
# - some splines are added for assumptions to hold
# - trajectory of HR for continuous covariates with splines is plotted

#############################################################################
#first load cleaned data:
load("cleandata.rda")

#packages needed:
library(ggplot2)
library(survival)
library(survminer)
library(tidyverse)
library(splines)
library(gridExtra)

#############################################################################
# -------- Kaplan Meier curve ----------------------------------------------
# plot the survival function S(t) (Kaplan meier)
sf <- Surv(time = datasub$survtime, event = datasub$status == 1)
km <-survfit(sf ~ 1, data = datasub)
plot(km, main = "Kaplan Meier 28 day survival curve", col =2)
summary(km)

#distribution survival times
hist(datasub$survtime,prob = F, breaks = 28, xlim = c(0,28),
     main = "Distribution of survival time")
#distribution survival times for status==1
hist(datasub$survtime[datasub$status==1],prob = F, breaks = 28, xlim = c(0,28),
     main = "Distribution of survival time for patients with event")
#distribution survival times for status==0
hist(datasub$survtime[datasub$status==0],prob = F, breaks = 28, xlim = c(0,28),
     main = "Distribution of survival time for patients with no event")
#spread of discharge destinations
plot(datasub$discharge_dest, las = 1)

#############################################################################
# ------- Centering and rescaling ---------------------------------------------

#The covariates we want to use in the model are: age, sex, number of cormorbidities,
#respiratory rate, peripheral oxygen saturation, ureum, creatinine and c-reactive protein.
#We have to take a closer look at them to see whether centering or rescaling is needed.

summary(datasub[,1:8])

hist(datasub$age)
hist(scale(datasub$age))

plot(datasub$sex)

hist(datasub$numb_comorb, main = "Distribution of number of comorbidities",
     xlab = "Number of cormorbidities")

hist(datasub$RR...min.)
plot(datasub$RR...min., main = "Spread of Respiratory rate per minute (RR min)")
range(datasub$RR...min., na.rm = T)

hist(datasub$SpO2....)
plot(datasub$SpO2...., main = "Spread of values Peripheral oxygen saturation") 
range(datasub$SpO2...., na.rm = T)

hist(datasub$CRP..mg.L.)
plot(datasub$CRP..mg.L., main = "Spread of values C-reactive protein")
range(datasub$CRP..mg.L., na.rm = T)

hist(datasub$Ureum.SER..mmol.L.)
plot(datasub$Ureum.SER..mmol.L., main="Spread of values Urea")

plot(datasub$Creatinine.SER..Âµmol.L.,  main="Spread of values Creatine")
hist(datasub$Creatinine.SER..Âµmol.L.)
#some below 40:
hist(datasub$Creatinine.SER..Âµmol.L., xlim = c(0,100), breaks=1000)
range(datasub$Creatinine.SER..Âµmol.L., na.rm = T)

# ----------------------------------------------------------------------------
#extreme values correction:

#RR:
#values above 60 strange, exaggerated measurements
#cut down to 60
datasub$RR...min. <- ifelse(datasub$RR...min. >60, 60, datasub$RR...min.)

#Peripheral oxygen saturation:
#values below 60% strange
datasub$SpO2.... <- ifelse(datasub$SpO2.... < 60, 60, datasub$SpO2....)

#creatinine:
#values above 300 strange 
#cut down to 300
datasub$Creatinine.SER..µmol.L. <- ifelse(datasub$Creatinine.SER..µmol.L. > 300, 
                                          300, datasub$Creatinine.SER..µmol.L.)

#------------------------------------------------------------------------------
#closer look at amount of missing values

#missingness of RR and SpO2
sum(is.na(datasub$RR...min.))
sum(is.na(datasub$SpO2....))
sum(is.na(datasub$RR...min.)& is.na(datasub$SpO2....))
#conclusion: a lot missing in pairs

#missingness urea and creatinine
sum(is.na(datasub$Creatinine.SER..µmol.L.))
sum(is.na(datasub$Ureum.SER..mmol.L.))
#will use creatine instead of ureum in model, because much less missing data

# -----------------------------------------------------------------------------
#centering:

#centering all covariates that have large values 
#(also to make sense of baseline hazard interpretation)

reshape_vars <- c("age", "RR...min.", "SpO2....", "CRP..mg.L.", "Creatinine.SER..µmol.L.")
rescale_df <- datasub |>
  summarise_at(vars(all_of(reshape_vars)), 
               .fun = list(~ mean(.x,  na.rm = T), ~ sd(.x,  na.rm = T))) 
reshape(rescale_df, direction = "long", 
        v.names = c("mean", "sd"), times = reshape_vars, timevar = "Variable",
        varying = 1:ncol(rescale_df), sep = "_")

rescale_df <- data.frame(  
  mean = c(mean(datasub$age), 
           mean(datasub$RR...min., na.rm = T),
           mean(datasub$SpO2...., na.rm = T), 
           mean(datasub$CRP..mg.L., na.rm = T), 
           mean(datasub$Creatinine.SER..µmol.L., na.rm=T)),
  sd = c(sd(datasub$age), 
         sd(datasub$RR...min., na.rm = T), 
         sd(datasub$SpO2...., na.rm = T), 
         sd(datasub$CRP..mg.L., na.rm = T), 
         sd(datasub$Creatinine.SER..µmol.L., na.rm = T))
)
row.names(rescale_df) <- c("age", "RR...min.", "SpO2....", "CRP..mg.L.", 
                           "Creatinine.SER..µmol.L.")
rescale_df

datasub <- datasub |> mutate_at(vars(age, RR...min., SpO2...., CRP..mg.L., Creatinine.SER..µmol.L.), ~ as.numeric(scale(.x, center = T))) 
#the scale function scales for the same values as in our rescale_df

#making sure first wave and second wave data are also scaled
datasub_w1 <- datasub[datasub$admitted_year_month <= "2020-06",]
datasub_w2 <- datasub[datasub$admitted_year_month > "2020-06",]


############################################################################
# ----- reverse KM ---------------------------------------------------------
#follow up distribution
#Follow-up tells you something about the maturity of the data, 
#hence about the reliability of your results.

datasub$status <- as.factor(datasub$status)
km2 <- survfit(Surv(survtime, status==0) ~ 1, data = datasub)
ggsurvplot(km2, risk.table = T)

custom_theme <- theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsurvplot(km2, risk.table = T, risk.table.y.text = F, 
           xlim = c(0, 28), xlab = "Days since admittance", break.x.by = 7,
           ggtheme = theme_bw(), tables.theme = custom_theme)

summary(km2) 

#############################################################################
# ------- Assumptions check ------------------------------------------------
# ------- Linearity and PH assumption ----------------------------------------
#only on first wave data

#Martingale residual plots for linearity assumption
#for continuous variables
p_res <- function (data, var, coxmodel) {
  df <- data %>% select(all_of(var))
  names(df) <- "Variable"
  indx <- complete.cases(df) #making sure NAs dont give error
  df$Residual <- NA #creating new column for residuals
  #compute residuals only of rows with complete cases of variable (no Nas)
  df$Residual[indx] <- residuals(coxmodel,  type = "martingale") 
  qplot(Variable, Residual, data=df) + 
    stat_smooth(method="loess", span=0.9, se=TRUE) +
    coord_cartesian(
      ylim = c(-0.5,0.5) #zooming in (no loss of data)
    ) +
    xlab(var)
}

#log-minus-log plot for proportional hazard assumption check
#for categorical variables only
km_lin <- function (data, var) {
  df <- data %>% select(all_of(var), survtime, status)
  names(df)[1] <- "Variable"
  df$Variable <- factor(df$Variable)
  km <- survfit(Surv(survtime, status) ~ Variable, data = df)
  ggsurvplot(km, data = df, fun = function(a) log(-log(a)), 
             xlim=c(0, 18),
             ylab = paste("Log(-Log(Survival)), Variable", var),
             title = paste("Log-minus-log plot of variable",var)) #hazard
}

#use coxph instead of km, use confouders (all other covar) and strata numb_comorb
cox_lin <- function (data, var, allvar) {
  df <- data %>% select(all_of(c(allvar,var)), survtime, status)
  #df$Variable <- factor(df$Variable)
  fmla <- formula(paste("Surv(survtime, status) ~",paste(allvar, collapse = "+"),
                        "+ strata(", var, ")"))
  coxfit <- coxph(fmla, data = df)
  bh <- basehaz(coxfit, centered = F) |>
    rename(cumhaz = hazard) |>
    group_by(strata) |>
    mutate(hazard = diff(c(0,cumhaz))) |>
    ungroup()
  
  ggplot(bh, aes(x= bh$time, y = bh$cumhaz, color=strata))+ geom_line() +
    ggtitle(paste("Cumulative hazard plot of variable",var)) + 
    xlab("Time")+ylab("Cumulative hazard")#cum hazard
}

# ------------ Discrete variables ---------------------------------------------
# PH check using survival curve plots

vardisc <- c("sex")
km_lin2 <- function (var) km_lin(datasub_w1, var)
p <- lapply(vardisc, km_lin2) 
p #lines dont cross => assumption ok 

#use strata and confouders for numb_cormb covariate
test <- cox_lin(datasub_w1, "numb_comorb", c("age", "sex", "RR...min.", 
                                             "SpO2....", "CRP..mg.L.", 
                                             "Creatinine.SER..µmol.L."))
test

# ----------- Continuous variables -------------------------------------------
#linearity check using residual plots

vars <- c("age", "sex", "numb_comorb", "RR...min.", "SpO2....", "CRP..mg.L.", 
          "Creatinine.SER..µmol.L.") #all covariates to include in cox model
fmla <- formula(paste("Surv(survtime, status) ~ ", paste(vars, collapse = "+")))
c_fit <- coxph(fmla, data=datasub_w1, method="breslow")
p_res2 <- function (var) p_res(datasub_w1, var, c_fit)
varcont <- c("age", "RR...min.", "SpO2....", "CRP..mg.L.", 
             "Creatinine.SER..µmol.L.")
q <- lapply(varcont, p_res2) 

grid.arrange(q[[1]], q[[2]], q[[3]]+ xlim(c(-1.2,1)), 
             q[[4]]+ xlim(c(-1,3)), q[[5]]+ xlim(c(-1,2)), 
             ncol = 3, nrow = 2,
             top = "Continuous variables vs Martingale residuals for linearity assumption") 
#closer look at plots of RR and SpO2:
q[[2]] + xlim(c(-1.5,2))
q[[3]] + xlim(c(-1.5,1)) 

# ---------------------------------------------------------------------------
#PH check:
#all covariates
cox.zph(c_fit) #want large p-values

# -----------------------------------------------------------------------------
# ------------------------------ adding splines -------------------------------

#only spline for RR
var3 <- c("age", "sex", "numb_comorb", "ns(RR...min., df =3)", 
          "SpO2....", "CRP..mg.L.", "Creatinine.SER..Âµmol.L.")

#both splines for RR and SpO2
var4 <- c("age", "sex", "numb_comorb", "ns(RR...min., df =3)", 
          "ns(SpO2...., df=3)", "CRP..mg.L.", "Creatinine.SER..Âµmol.L.")

#both splines for Creatinine
var5 <- c("age", "sex", "numb_comorb", "RR...min.", 
          "SpO2....", "CRP..mg.L.", "ns(Creatinine.SER..Âµmol.L., df = 4)")

#var2 <- c("ns(age, df=4)", "sex", "numb_comorb", "ns(RR...min., df =4)", 
#          "ns(SpO2...., df = 3)", "CRP..mg.L.", "Creatinine.SER..?mol.L.")

#fmla2 <- formula(paste("Surv(survtime, status) ~ ", paste(var2, collapse = "+")))
fmla3 <- formula(paste("Surv(survtime, status) ~ ", 
                       paste(var3, collapse = "+")))
fmla4 <- formula(paste("Surv(survtime, status) ~ ", 
                       paste(var4, collapse = "+")))


#c_fit2 <- coxph(fmla2, data=datasub_w1, method="breslow", x = T)
#p_res2 <- function (var) p_res(as.data.frame(c_fit2$x), var, c_fit2)
#varcont <- names(c_fit2$coefficients)

c_fit3 <- coxph(fmla3, data=datasub_w1, method="breslow", x = T)
p_res3 <- function (var) p_res(as.data.frame(c_fit3$x), var, c_fit3)
varcont3 <- names(c_fit3$coefficients)

c_fit4 <- coxph(fmla4, data=datasub_w1, method="breslow", x = T)
p_res4 <- function (var) p_res(as.data.frame(c_fit4$x), var, c_fit4)
varcont4 <- names(c_fit4$coefficients)

#q2 <- lapply(varcont, p_res2) 
q3 <- lapply(varcont3, p_res3)
q4 <- lapply(varcont4, p_res4)

grid.arrange(q3[[1]], q3[[4]],q3[[5]],q3[[6]], q3[[7]], q3[[8]], q3[[9]], 
             #q3[[10]],
             ncol = 3, nrow = 3) 

#adding one more spline on SpO2 as well gives improvement
grid.arrange(q4[[1]], q4[[4]],q4[[5]],q4[[6]], q4[[7]], q4[[8]], q4[[9]], 
             q4[[10]], q4[[11]],# q4[[12]],
             ncol = 3, nrow = 3,
             top = "Continuous variables vs Martingale residuals for linearuty assumption") 
#looks better
q4[[6]]

#grid.arrange(q2[[1]], q2[[2]], q2[[3]], q2[[4]], q2[[7]], q2[[8]], q2[[9]], 
#             q2[[10]], q2[[11]], q2[[12]], q2[[13]], q2[[14]],
#             ncol = 3, nrow = 4)

#cox.zph(c_fit2)
cox.zph(c_fit3)
cox.zph(c_fit4)
#still not so good

##############################################################################
#plot hazard ratio trough range of values of each of the covariates with a spline
#to see the effect

#baseline value
df0.1 <- datasub_w1[1,] |>
  select(- c(RR...min.)) |> #covariate of interest removed
  #set other covariates to baseline value
  mutate(age = 0, sex = "F", numb_comorb = 0, CRP..mg.L. =0, 
         Creatinine.SER..?mol.L.=0, SpO2.... = 0)
#range of values selected coavariate to plot
df.1 <- data.frame(RR...min. = seq(-3.1, 4.9, 0.1)) 
#predict
sp1.1 <- predict(c_fit4, 
               newdata = expand_grid(df0.1, df.1),
               type = "risk")

plot1 <- ggplot(NULL, aes(x=df.1[,1], y= sp1.1))+
  geom_line()+
  labs(x="RR", y="HR", title = "Hazard ratio for respiratory rate")+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")

hist(datasub_w1$SpO2....)
range(datasub_w1$SpO2...., na.rm = T)

#baseline value
df0.2 <- datasub_w1[1,] |>
  select(- c(SpO2....)) |>
  #set other covariates to baseline value
  mutate(age = 0, sex = "F", numb_comorb = 0, CRP..mg.L. =0, 
         Creatinine.SER..?mol.L.=0, RR...min. = 0)
#range of values selected coavriate to plot
df.2 <- data.frame(SpO2.... = seq(-5, 1.1, 0.1)) 
#predict
sp1.2 <- predict(c_fit4, 
               newdata = expand_grid(df0.2, df.2),
               type = "risk")

plot2 <- ggplot(NULL, aes(x=df.2[,1], y= sp1.2))+
  geom_line()+
  labs(x="SpO2", y="HR", 
       title = "Hazard ratio for peripheral oxygen saturation")+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")


grid.arrange(plot1, plot2,
             ncol = 2,
             top = "Trajectory of hazard ratio for respiratory rate and peropheral oxygen saturation")

