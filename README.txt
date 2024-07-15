Project: Dynamical updating

In this file an overvieuw is given of the different R scripts used in this project, with an explanation of each file.
More detailed explanation is provided in the beginning of each file.

- thesis-data-extracting.R
	In this file the original data is loaded and adapted for the analysis that will be performed later. 
	Variables are converted to factors or dates, and some new variables are created (ICU admission indicator, event indicator, 
	survival time variable and comorbidities count variable).
	Exclusion criteria is applied. 
	The final subset created (named "datasub") consist of all the necessary variables that we will need.
	Also two seperate datasets are created with data from the first and second wave only (named "datasub_w1" and "datasub_w2").
	These tree dataset are saved in cleandata.rda.

- thesis-data-descriptives.R
	In this file some descriptives plots are provided (which will be used in Descriptives chapter in thesis).

- thesis-data-assumptions.R
	In this file the assumtpions for Cox model of the first wave are checked. Data is loaded from cleandata.rda.
	The covariates in the model are also first inspected and scaled.
	Some extreme values of the covairates are cut based on expert knowledge.
	The final model, of the first wave that holds the assumtpions the best, is constructed here.  

- thesis-data-cleaning.R
	In this file the data from thesis-data-extracting is cleaned further.
	Names are cleaned, extreme values for some variables are cut and all continuous varibales are scaled. The resulting dataset is called datasub_clean.
	Two datasets are created for wave 1 and 2, both contain only the complete cases.
	This data is then saved in cleanscaledcut-cc-data.rda.

- thesis-validation-functions.R
	In this file the validation functions are provided.
	In file "thesis-updating-functions" these functions are used to evaluate the performance of the updated models.

- thesis-updating-functions.R
	In this file several functions are provided that update a model. 
	There is a function for no-updating, refitting, and recalibration of the intercept.
	Additional, functions are provided to generate predictions from Bayesian updated models.

- thesis-updating-Bayesian.R
	In this file Bayesian updating is performed (aditional Stan file required: expSurvivalModel.stan).
	The results of the updates are saved in: bayes_data_list_50.rda, bayes_data_list_100.rda, and bayes_data_list_150.rda.
	Code for plotting the coefficients and log baseline hazard are profided. 

- thesis-updating.R
	In this file the model of the first wave is again inspected.
	Updating periods of the second wave are build for different lengths and different sliding windows. 
	Also periods for validation are created.
	Then, no-updating, refitting, and recalibration are applied to the periods.
	Also the updated Bayesian models from bayes_data_list_50.rda, bayes_data_list_100.rda, and bayes_data_list_150.rda are evaluated.
	For the resuls of the performance of updated models, code for plots and numerical summaries are provided.
	Code to plot the updated coefficients from the refitting method is also provided.
	Lastly, code to plot an illustration of the first updating period of 50 events is given, as well as visualization of individual predicted risk of refitting and recalibration (of one period).

- expSurvivalModel.stan
	In this file Stan code is provided for the Bayeisan updating method.