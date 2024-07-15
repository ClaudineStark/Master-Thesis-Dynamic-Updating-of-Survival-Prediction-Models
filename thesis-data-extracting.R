################################################################################
# In this script:
# - reading in the data
# - converting variables to factor and dates
# - creating new variables for:
#       - ICU admission indicator
#       - status (event) indicator
#       - survival time
# - excluding patients based on exclusion criteria
# - creating number of comorbidities count variable
# - creating subsample with only necessary covariates for future model and other
#   needed variables
# - creating subsamples with data from first and second wave

# Resulting datasets for further use: 
# - datasub: data with useful variables included
# - datasub_w1: subset of datasub with only wave 1 data
# - datasub_w2: subset of datasub with only wave 2 data

################################################################################
# ---- LOADING DATA  ---------------------------------------------------------
library(dplyr)
library(readxl)

#loading in the data:
data <- read_xlsx() #fill in path to data

# --- CONVERTING VARIABLES ---------------------------------------------------

#converting to factor
data <- data |>
  mutate_at(
    vars(sex, died, Chloroquine, Dexamethason, Nadroparine, Remdesivir, no_IC,
         no_intubate, no_resuscitate, discharge_dest),
    ~ as.factor(.))

#convert to date
data$admitted_date <- as.Date(data$admitted_date, "%Y-%m-%d %H:%M:%S") 
#extracting month and year of admission
data$admitted_year_month <- format(data$admitted_date, "%Y-%m")

#convert to date
data$died_date <- as.Date(data$died_date, "%Y-%m-%d %H:%M:%S")
#extracting month and year of death
data$died_date_month <- format(as.Date(data$died_date, "%Y-%m-%d %H:%M:%S"),
                               "%Y-%m")

#converting to date
data$date_ICU <- as.Date(data$admitted_to_IC, "%Y-%m-%d %H:%M:%S")
#covariate giving month of ICU admission
data$date_ICU_month <- format(as.Date(data$admitted_to_IC, "%Y-%m-%d %H:%M:%S"),
                              "%Y-%m")

#converting to date
data$discharge_date <- as.Date(data$discharge_date, "%Y-%m-%d %H:%M:%S")
#convert discharge date to month+year 
data$discharge_date_month <- format(as.Date(data$discharge_date, 
                                            "%Y-%m-%d %H:%M:%S"),"%Y-%m")

#converting yes/no to 1/0
data <- data |>
  mutate_at(
    vars(Chronische.Hartziekte, Hypertensie, Chronische.Longziekte, Astma, 
         Chronische.Nierziekte, Chronische.Leverziekte, 
         Chronische.neurologische.aandoening,
         HIV, Diabetes, TBC, Asplenie, Maligniteit, Hematologische.maligniteit,
         Metabool.syndroom, Niertransplantatie),
    ~ if_else(.x == "Ja", 1, 0))

# --- CREATING NEW VARIABLES --------------------------------------------------

#covariate indicates ICU admission
data$ICU_adm <- if_else(is.na(data$admitted_to_IC), 0, 1)
#data$ICU_adm <- as.factor(data$ICU_adm) #not needed? 

## STATUS
#covariate indicating icu or death or discharge to hospice 
data$status <- if_else(
  data$died == 1 | data$ICU_adm == 1 | 
    (!is.na(data$discharge_dest) & data$discharge_dest == "Zorginstelling"),
  1, # yes 
  0  # no event (right censored) 
)

## SURVIVAL TIME
#variable for date ICU admission or death or censored 
data$survtime <- case_when(
  data$status == 1 & data$ICU_adm == 1 ~ 
    as.numeric(difftime(data$date_ICU, data$admitted_date, units = "days")),
  data$status == 1 & data$died == 1 ~ 
    as.numeric(difftime(data$died_date, data$admitted_date, units = "days")),
  data$status == 1 ~  as.numeric(difftime(data$discharge_date, 
                                          data$admitted_date, units = "days")),
  !is.na(data$discharge_dest) & 
    data$discharge_dest == "Eigen woonomgeving" ~ 28.001,
  TRUE ~ as.numeric(difftime(data$discharge_date, data$admitted_date, 
                             units = "days")) 
) 


#time horizon of 28 days
#max survival time is set to 28.001 days and status is set to 0
data$status <- ifelse(data$survtime > 28.001, 0, data$status)
#data$status <- as.factor(data$status)
#also put icu indicator to 0 if icu admission was after 28 days
data$ICU_adm <- ifelse(data$survtime > 28.001, 0, data$ICU_adm) #ok?
data$survtime <- ifelse(data$survtime > 28.001, 28.001, data$survtime)


# ---- exclusion criteria ----------------------------------------------------
#removing patients transferred from other hospital
index1 <- which(data$admission_origin == "Ander ziekenhuis")
index2 <- which(data$admission_origin == "Ziekenhuis buitenland")

data <- data[-index1,]
data <- data[-index2,]

#removing non existing patient
data <- filter(data, !rowSums(is.na(data)) + rowSums(data == 0, na.rm = T) >= ncol(data)-5)

nrow(data)
#4064 individuals left

# ------ Selecting covariates to use -----------------------------------------
#names(data)

# creating count variable of comorbidities
# based on 4C comorbidity count variable
com_vars <- c("Chronische.Hartziekte", "Chronische.Longziekte",
              "Chronische.Leverziekte","Diabetes", "Maligniteit", 
              "Hematologische.maligniteit", 
              "Chronische.Nierziekte", 
              "Hypertensie", #not included in 4c because not recorded in data
              "Chronische.neurologische.aandoening", #not in 4C but better to 
              #include from clinical perspective
              "Astma") # astma also possible to included
              # obesity not in dataset
data$numb_comorb <- rowSums(data[,com_vars]) 

#making range numb_comorb 0-1-2, with 2 meaning 2 or more
data$numb_comorb <- ifelse(data$numb_comorb >2, 2, data$numb_comorb)

# selecting covariates to use in next steps of analysis
covar <- c(#covariates selected according to 4C:
           "age", 
           "sex", 
           "numb_comorb", #comorbidities (combined score)
           "RR...min.", #Respiratory  rate
           "SpO2....", #Peripheral oxygen saturations
           #Glasgow comma (not available!)
           "Ureum.SER..mmol.L.", #Urea
           "CRP..mg.L.",#c reactive protein
           "Creatinine.SER..µmol.L.", #alternative for ureum values
           #other variables also needed (dates, ...):
           "pseudo_id", 
           "admitted_date", 
           "discharge_date", 
           "admitted_loc", 
           "admission_origin", 
           "admitted_to_IC", 
           "discharged_from_IC", 
           "discharge_dest", 
           "died", 
           "died_date", 
           "hospital",
           "admitted_year_month" , 
           "died_date_month", 
           "date_ICU",                                 
           "date_ICU_month", 
           "discharge_date_month", 
           "ICU_adm", 
           "status", "survtime")

# creating dataset with only necessary covariates included
datasub <- data[,covar]


############################################################################
# ----------- Data splitting -----------------------------------------------
# splitting dataset in two the two waves

# first wave:
datasub_w1 <- datasub[datasub$admitted_year_month <= "2020-06",]
#dim(datasub_w1)

# second wave:
datasub_w2 <- datasub[datasub$admitted_year_month > "2020-06",]
#dim(datasub_w2)

###########################################################################
# ---------- Save data ----------------------------------------------------

#saving cleaned data as a rda file
save(datasub, datasub_w1, datasub_w2, file = "cleandata.rda")

