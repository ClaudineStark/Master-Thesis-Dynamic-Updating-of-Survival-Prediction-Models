###############################################################################
# In this script:
# - cleandata.rda is loaded
# - extreme values are cut
# - names from dataset are "cleaned"
# - continuous variables are scaled
# - complete cases for waves 1 and 2 are selected
# - cleaned-scaled-cut data is saved (complete cases for waves 1 and 2)

#############################################################################
library(janitor)
library(dplyr)
library(ggplot2)

#load data 
load("cleandata.rda")

###############################################################################
# ----- cutting extreme values ------------------------------------------------

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

##############################################################################
# ----- cleaning column names ------------------------------------------------

names(datasub)
#making column names "nicer" - creating new dataframe with cleaned names
datasub_clean = clean_names(datasub, replace = c("µ"= "mu"))
names(datasub_clean)

##############################################################################
# ---- scaling continuous covariates -----------------------------------------
#applying to df with nicer names 

rescale_df <- data.frame(  
  # Create a df where we can find this information when needed
  mean = c(mean(datasub_clean$age), 
           mean(datasub_clean$rr_min, na.rm = T), 
           mean(datasub_clean$sp_o2, na.rm = T), 
           mean(datasub_clean$crp_mg_l, na.rm = T),
           mean(datasub_clean$creatinine_ser_mumol_l, na.rm = T)),
  sd = c(sd(datasub_clean$age), 
         sd(datasub_clean$rr_min , na.rm = T), 
         sd(datasub_clean$sp_o2, na.rm = T), 
         sd(datasub_clean$crp_mg_l, na.rm = T),
         sd(datasub_clean$creatinine_ser_mumol_l, na.rm = T))
)
row.names(rescale_df) <- c("age", "rr_min", "sp_o2", "crp_mg_l", 
                           "creatinine_ser_mumol_l")
rescale_df

datasub_clean <- datasub_clean |> 
  mutate_at(vars(age, rr_min, sp_o2, crp_mg_l, creatinine_ser_mumol_l), 
            ~ as.numeric(scale(.x, center = T))) 
#the scale function scales for the same values as in our rescale_df

#making sure first wave and second wave data are also scaled
datasub_w1 <- datasub_clean[datasub_clean$admitted_year_month <= "2020-06",]
datasub_w2 <- datasub_clean[datasub_clean$admitted_year_month > "2020-07",] 

##############################################################################
# --- complete cases wave 1 --------------------------------------------------
#covariates in model
var1 <- c("age", "sex","numb_comorb","rr_min" ,"sp_o2", "crp_mg_l", 
          "creatinine_ser_mumol_l")
#complete cases with respect to the covariates
ind1 <- complete.cases(subset(datasub_w1, select = var1))
datasub_w1_cc <- datasub_w1[ind1, ]
nrow(datasub_w1_cc) #total number

###############################################################################
# ---- complete cases wave 2 --------------------------------------------------
#using only complete cases for now:
#covariates in model
var2 <- c("age", "sex","numb_comorb","rr_min" ,"sp_o2", "crp_mg_l", 
          "creatinine_ser_mumol_l")
#complete cases with respect to the covariates
ind2 <- complete.cases(subset(datasub_w2, select = var2))
datasub_w2_cc <- datasub_w2[ind2, ]
nrow(datasub_w2_cc) #total number

###############################################################################
#some plots:

ggplot(data=datasub_w1_cc, aes(x=survtime, fill=admitted_year_month)) +
  geom_histogram() +
  labs(x = "days", fill = "admission month")

qplot(datasub_w1_cc$survtime[datasub_w1_cc$admitted_year_month == "2020-06"])


###############################################################################
# ------- Making variable that counts events in second wave ------------------
#usefull when making intervals based on number of events

#creating column that has date of event (only for patients with an event)
datasub_w2_cc$date_event <- case_when(
  #date of event 
  datasub_w2_cc$status == 1 ~ datasub_w2_cc$admitted_date + datasub_w2_cc$survtime, 
  TRUE ~ NA #for patients with no event
)
#ordering data wave 2 by event date (increasing)
datasub_w2_cc <- datasub_w2_cc[order(datasub_w2_cc$date_event,
                                     datasub_w2_cc$admitted_date),]
#when ties occur: order by admission date

datasub_w2_cc$count <- cumsum(datasub_w2_cc$status)

###############################################################################

#save cleaned, cut, scaled and complete cases (wave 1 and wave 2) data
save(datasub_clean, datasub_w1_cc, datasub_w2_cc, file = "cleanscaledcut-cc-data.rda")

