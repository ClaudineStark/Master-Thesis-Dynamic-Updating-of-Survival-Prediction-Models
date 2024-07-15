###############################################################################
# In this file:
# - first: extreme values of variables are cut (based on clinical expertise)
# - plots of missing values
# - descriptives table of covariates
# - other plots

#load data
load("cleandata.rda")

# ----------------------------------------------------------------------------
#load packages
library(ggplot2)
library(gridExtra)
library(naniar)
library(summarytools)
library(gtsummary)

#----------------------------------------------------------------------------
# ----- cutting extreme values ------------------------------------------------
#extreme values correction
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

#----------- Missing values --------------------------------------------------
names(datasub)
#missing data in covariates (excludin ureum)
gg_miss_var(datasub[,c(1:5,7:8)], show_pct = T)+
  ggtitle("Proportion of missing data in covariates")#+
  #ylim(c(0,100))

#all covaiates, also ureum
gg_miss_var(datasub[,1:8], show_pct = T)+
  ggtitle("Proportion of missing data in covariates")

#missing data in the other variables
gg_miss_var(datasub[,9:27], show_pct = T)+
  ggtitle("Proportion of missing data in variables")


# ---------- Descriptive summaries ------------------------------------------
#number of events per admission month
table(datasub$status, by= datasub$admitted_year_month)
#number of patients per hospital
table(datasub$hospital)

#making variable that indicates from what wave people are
datasub$wave <- as.factor(ifelse(datasub$admitted_year_month <= "2020-06", "wave 1", "wave 2"))
#nice summary table of all 7 covariates!
datasub1 <- datasub[,c(1:8,28)]
datasub1$"Age" <- datasub1$age
datasub1$numb_comorb <- ifelse(datasub1$numb_comorb == 2, "2 or more", datasub1$numb_comorb)
datasub1$"Sex" <- ifelse(datasub1$sex == "F", "Female", "Male")
datasub1$"Comorbidity count" <- datasub1$numb_comorb
datasub1$"Respiratory Rate (per min)" <- datasub1$RR...min.
datasub1$"Peripheral oxygen saturations (%)" <- datasub1$SpO2....
datasub1$"C-reactive protein (mg/l)" <- datasub1$CRP..mg.L.
datasub1$"Creatinine (µmol/l)" <- datasub1$Creatinine.SER..µmol.L.
subnam <- datasub1 %>% select("Age", "Sex", "Comorbidity count", 
                              "Respiratory Rate (per min)", 
                              "Peripheral oxygen saturations (%)", 
                              "C-reactive protein (mg/l)",
                              "Creatinine (µmol/l)", wave) 
  
tbl_summary_1 <- tbl_summary(subnam, by = wave,
                             statistic =  
                               list(all_continuous() ~ "{mean} ({sd})",
                                    all_categorical() ~ "{n} ({p}%)")) %>% 
                  add_overall()
tbl_summary_1 

# ------------ Visual summaries ----------------------------------------------

#distribution admissions per month
plot1.1 <- ggplot(data=datasub, aes(x=admitted_year_month)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of hospital admissions per month")+
  xlab("year and month of admission")
plot1.1

#plot of age distribution per sex
plot2 <- ggplot(data = datasub, aes(x =age, fill = sex)) +
  geom_histogram() +ggtitle("Distribution age per sex")
plot2

#covariates plots:

#plot sex
p1 <-ggplot(data=datasub, aes(x=sex))+ 
  geom_bar()+ ggtitle("Distribution of sex")
#plot age
p2 <-ggplot(data=datasub, aes(x=age))+ 
  geom_histogram()+ ggtitle("Distribution of age")
#plot numb corm
p3 <-ggplot(data=datasub1[complete.cases(datasub1$numb_comorb),], aes(x=numb_comorb))+ 
  geom_bar( )+ ggtitle("Distribution of comorbidities")+xlab("number of comorbidities")
#plot RR...min. 
p4 <- ggplot(data = datasub, aes(x= admitted_year_month, y = RR...min. ))+
  geom_point() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of Respiratory rate (per min)") +
  xlab("admission month and year")
#plot SpO2....
p5 <- ggplot(data = datasub, aes(x= admitted_year_month, y = SpO2....))+
  geom_point() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of Peripheral oxygen sat. (%)") +
  xlab("admission month and year")
#plot CRP..mg.L.
p6 <- ggplot(data = datasub, aes(x= admitted_year_month, y = CRP..mg.L.))+
  geom_point() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of C-reactive protein (mg/l)")+
  xlab("admission month and year")
p7 <- ggplot(data = datasub, aes(x= admitted_year_month, y = Creatinine.SER..µmol.L.))+
  geom_point() +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of Creatinine (µmol/l)")+
  xlab("admission month and year")

grid.arrange(p2, p1, p3,
             nrow =2, ncol=2, 
             top="Distribution of covariates (part 1)")

grid.arrange(p4, p5, p6, p7, p2,
             nrow =3, ncol=2, 
             top="Distribution of covariates")

#---------------------------------------------------------------------------

#plot of censoring times
hist_censor <- ggplot(data = datasub[datasub$status == 0,], aes(x=survtime))+
  geom_histogram() +
  ggtitle("Distribution of censoring time") + xlab("days")

#plot of time to event
hist_time <- ggplot(data = datasub[datasub$status == 1,], aes(x=survtime))+
  geom_histogram()+
  ggtitle("Distribution of time to event") + xlab("days")

grid.arrange(hist_time, hist_censor,
             ncol = 2, top= "Distribution event and censoring time")

#distribution event rate
datasub$status <- as.factor(datasub$status)
#event = 1 = died, ICU admission or discharge to hospice (zorginstelling)
plot1 <- ggplot(data = datasub, aes(x=admitted_year_month,  after_stat(count))) +
  geom_bar(aes(fill = status), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of events through time") +
  xlab("year and month of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))
plot1

#plot of total number of events
d1 <- ggplot(data = datasub, aes(x=status, fill = status))+
  geom_bar()+
  ggtitle("Total distribution of events")+
  scale_x_discrete(name = "Event", labels=c("0" = "No", "1" = "Yes"))+
  theme(legend.position="none")
d1

#total spread of events and per month
grid.arrange(d1, plot1, ncol=2, 
             top = "Distribution of events")


#plot of total number of events
data$status<-as.factor(data$status)
plot1 <- ggplot(data = data, aes(x=admitted_year_month,  after_stat(count))) +
  geom_bar(aes(fill = status), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of events through time") +
  xlab("year and month of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))
plot1

#admission per hospital per month
data <- datasub_clean
plot1.1 <- ggplot(data = data, aes(x=admitted_year_month,  after_stat(count))) +
  geom_bar(aes(fill = hospital), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of admissions per month for each hospital hospital") +
  xlab("year and month of admission") +
  scale_fill_discrete(name = "Hospital")
plot1.1

datasub <- data[data$admitted_year_month<="2020-12",]
plot2 <- ggplot(data = datasub, 
                aes(x=admitted_year_month,  after_stat(count))) +
  geom_bar(aes(fill = status), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of events by month") +
  xlab("year and month of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))
#plot2

plot3 <- ggplot(data = datasub, 
                aes(x=admitted_date,  after_stat(count))) +
  geom_bar(aes(fill = status), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of events by day") +
  xlab("month and day of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))
#plot3

library(gridExtra)
grid.arrange(plot2, plot3,
             nrow = 2, ncol =1,
             top = "Distribution of events wave 2 (first 5 months)")


datasub$icu_adm<-as.factor(datasub$icu_adm)
datasub$died<-as.factor(datasub$died)
plot4 <- ggplot(data = datasub, 
                aes(x=admitted_date,  after_stat(count))) +
  geom_bar(aes(fill = icu_adm), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of ICU admissions by day") +
  xlab("month and day of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))+
  ylim(c(0,26))
#plot4
plot5 <- ggplot(data = datasub, 
                aes(x=admitted_date,  after_stat(count))) +
  geom_bar(aes(fill = died), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of died by day") +
  xlab("month and day of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))+
  ylim(c(0,26))
#plot5

grid.arrange(plot4, plot5,
             nrow = 2, ncol =1,
             top = "Distribution of ICU admissions and death wave 2 (first 5 months)")

plot3.1 <- ggplot(data = datasub[datasub$survtime<28 & datasub$status==0,], 
                  aes(x=admitted_date,  after_stat(count))) +
  geom_bar(aes(fill = status), position = "dodge")+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("Distribution of censored patients (with survtime <28) by day") +
  xlab("month and day of admission") +
  scale_fill_discrete(name = "Event", labels = c("no", "yes"))
#plot3.1
