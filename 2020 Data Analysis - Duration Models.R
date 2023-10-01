#Whatcom TAST - 2020 Analysis - Duration (surface counts) 
#Script written by Kathleen McKeegan
#Last update made 9/30/2023

## Working Directory: 

setwd("~/Desktop/TAST-  Data Cache and R Scripts/Data - csv")

## Load Packages
library('tidyverse')
library('glmmTMB')
library('sjPlot') #for random effects plots
library('DHARMa')
library('emmeans')
library('performance')

#### Upload and format 2020 TAST data ####
tast.dat <- read.csv("2020_data.csv", header=T)
tast.dat$Fish<- sub("^$", "0", tast.dat$Fish)
tast.dat$Fish[tast.dat$Fish == "Y"] <- "1"
tast.dat$Fish <- as.numeric(tast.dat$Fish)
tast.dat$Date <- as.Date(tast.dat$Date, "%m/%d/%y")
tast.dat$ID <- str_pad(tast.dat$ID, 4, pad = "0") #make all IDs a 4-digit number
tast.dat$ID[tast.dat$ID == "Sea_lion"]  <- "sea_lion" #fix typo in sea_lion impacting summary
tast.dat$ID[tast.dat$ID == "00CF"] <- "0198"
tast.dat$ID[tast.dat$ID == "0265"] <- "0072"
tast.dat$ID <- as.factor(tast.dat$ID)
tast.dat$Salmon3day <- as.character(tast.dat$Salmon3day)
tast.dat$Salmon3day <- as.numeric(tast.dat$Salmon3day)
tast.dat$Salmon5day <- as.character(tast.dat$Salmon5day)
tast.dat$Salmon5day <- as.numeric(tast.dat$Salmon5day)
tast.dat$TAST <- as.factor(tast.dat$TAST)
tast.dat$Camera <- as.factor(tast.dat$Camera)

#Remove any RS from the New IDs, keep only LS since we know they are unique individuals

tast.dat <- subset(tast.dat, ID!= '00AO'& ID!='00BA'& ID!= '00BC'& ID!='00BD'& ID!='00BF'& ID!='00BG'& 
                     ID!='00BR'& ID!='00BT'& ID!='00BW'& ID!='00BY'& ID!='00BZ'& ID!='00CA'& ID!='00CJ'&
                     ID!= '00CL'& ID!='00CT'& ID!='00CU' & ID!= '00CW' & ID!="sea_lion")

## Create dataframe with summary of occurrence and catches
# Narrow it down to the experimental window 10/25-11/25
tast.summary <- tast.dat %>%
  group_by(Date, ID, Length, Camera, Tide_Height, TAST, Salmon5day, Salmon3day) %>%
  summarise(Occurrence = n_distinct(Time),
            Catches = sum(Fish)) %>%
  filter(Date >= '2020-10-25' & Date <= '2020-11-25')

# Add presence column in case it's needed
tast.summary <- tast.summary %>%
  mutate(Presence = 1)

# Get rid of NAs in Salmon rolling average
tast.summary[is.na(tast.summary)] <- 0

# Add column for observation days across TAST status (potential offset term)
tast.summary <- tast.summary %>%
  mutate(ObsDays = if_else( TAST == 'OFF', 14, 16))

# Add column with sum total number of seals seen per observation period
tast.summary <- tast.summary %>%
  group_by(Date) %>%
  mutate(Total=length(ID))

# Add column with Observation length in hours, not minutes (in case it's needed)
tast.summary <- tast.summary %>%
  mutate(LengthHR = (Length/60))

# Add column with Occurrence as proportion (surface counts/observation length)
tast.summary <- tast.summary %>%
  mutate(propDuration = (Occurrence/Length))

# Add column with Across-years site fidelity
ID.fidelity <- read.csv("2020_fidelitytags.csv", header=T)
ID.fidelity$ID <- str_pad(ID.fidelity$ID, 4, pad = "0") #make all IDs a 4-digit number
ID.fidelity$ID <- as.factor(ID.fidelity$ID)
ID.fidelity$Tag <- as.factor(ID.fidelity$Tag)
tast.summary <- merge(tast.summary, ID.fidelity, by= "ID", all.x=TRUE)
tast.summary <- tast.summary[order(as.Date(tast.summary$Date, format="%d/%m/%Y")),]

#### GLMM Analysis for Duration ####

# Response Variable is discrete (Surface Counts) -- need poisson or negative binomial
## Surface counts = how many times a seal surfaced during observation/
## on a minute-increment (proxy for Duration in min)

# Model Selection process: 
## Step 1 - optimal combination of random effects with fully populated fixed effect term constant
## Step 2 - keep selected random effects and find optimal fixed effect combination

# Random Effects: ID, Date, number of cameras (observers present)
# Fixed Effects: TAST, tide, rolling average of Salmon, fidelity tag

#### Step 1 - fit random effects ####
hist(tast.summary2$Occurrence)


# Model1a - see is poisson distribution is ok
model1a <- glmmTMB(Occurrence ~ TAST + Tide_Height + Salmon5day + Tag + (1|ID), 
                  data=tast.summary, offset=log(LengthHR), family=poisson(link='log'))
testDispersion(model1a) #could be ok...
check_overdispersion(model1a) #nope, over dispersed. Try negative binomial

# Negative binomial will handle overdispsed data
# Response variable is zero truncated, all seals are present for at least one minute (one count)
# use zero truncated neg binom for model fitting

# Model 1b - Try ID as random with Observation offset (hr) and negative binomial distribution
model1b <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + scale(Salmon5day) + Tag + (1|ID), 
                  data=tast.summary, offset=log(LengthHR), family=truncated_nbinom2(link='log'))
testDispersion(model1b) 
check_overdispersion(model1b) #much better, use this distribution

# Model 2 - Try Date as random with Observation Offset 
model2 <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + scale(Salmon5day) + Tag + (1|Date), 
                  data=tast.summary,offset=log(LengthHR), family=truncated_nbinom2(link='log'))

# Model 3 - Try Camera as random with Observation offset 
model3 <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + scale(Salmon5day) + Tag  +(1|Camera), 
                  data=tast.summary,offset=log(LengthHR), family=truncated_nbinom2(link='log'))

# Model 4 - Try ID and Date as crossed random effects with Observation offset
model4 <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + scale(Salmon5day) + Tag + (1|ID) + (1|Date), 
                  data=tast.summary,offset=log(LengthHR), family=truncated_nbinom2(link='log'))

#Compare models with different random effects
AIC(model1b, model2, model3, model4) 
#Model4 with ID and Date as crossed random effects is the best


#### Step 2 - fit fixed effects ####

# Model 5 - just TAST
model5 <- glmmTMB(Occurrence ~ TAST + (1|ID) + (1|Date), data=tast.summary,
                  offset=log(LengthHR), family=truncated_nbinom2(link='log'))

# Model 6 - TAST and Tide 
model6 <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + (1|ID)+ (1|Date), data=tast.summary,
                  offset=log(LengthHR), family=truncated_nbinom2(link='log'))

# Model 7 - TAST and Salmon 
model7 <- glmmTMB(Occurrence ~ TAST + scale(Salmon5day) + (1|ID) + (1|Date), 
                  data=tast.summary,
                  offset=log(LengthHR), 
                  family=truncated_nbinom2(link='log'))

# Model 8 - TAST and Tag
model8 <- glmmTMB(Occurrence ~ TAST + Tag + (1|ID)+ (1|Date), data=tast.summary,
                  offset=log(LengthHR), family=truncated_nbinom2(link='log'))

# Model 9 - TAST and Tide and Salmon
model9 <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + scale(Salmon5day) + (1|ID)+ (1|Date), 
                  data=tast.summary,offset=log(LengthHR), 
                  family=truncated_nbinom2(link='log'))

# Model 10 - TAST and Tide and Tag
model10 <- glmmTMB(Occurrence ~ TAST + scale(Tide_Height) + Tag + (1|ID)+ (1|Date), 
                   data=tast.summary,offset=log(LengthHR), 
                   family=truncated_nbinom2(link='log'))

# Model 11 - TAST and Salmon and Tag
model11 <- glmmTMB(Occurrence ~ TAST + scale(Salmon5day) + Tag + (1|ID)+ (1|Date), 
                   data=tast.summary,offset=log(LengthHR),
                   family=truncated_nbinom2(link='log'))

#Compare models to select fixed effects

AIC(model4, model5, model6, model7, model8, model9, model10, model11)
# All models are within 2.0 AIC of each other
# Model 7 is ever so slightly better than the rest
# Run diagnostics to see how the model looks

summary(model7)

#Model 7 diagnostics: TAST and Salmon with ID and Date as random 
testDispersion(model7)
testUniformity(model7)
testOutliers(model7)
simulationOutput7 <- simulateResiduals(fittedModel=model7, plot=T)
plot(simulationOutput7) #significant deviation 
plotResiduals(simulationOutput7, tast.summary$TAST) #issue in TAST on group
plotResiduals(simulationOutput7, tast.summary$Salmon5day)

#Plot model 7
plot_model(model7, type='pred')
ob_d <- emmeans(model7, ~TAST, type='response')
ob_d
plot(ob_d, xlab='Model Estimates: Duration (minute counts/hour)') +
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 13))

#exponentiated estimates and confidence intervals 
exp(confint(model7))

#Look at random effects
plot_model(model7, type='re', show.values = TRUE, value.offset = .6)


tast.summary.on <- tast.summary %>%
  filter(TAST == "ON")
hist(tast.summary.on$Occurrence)
tast.summary.off <- tast.summary %>%
  filter(TAST == "OFF")
hist(tast.summary.off$Occurrence)
