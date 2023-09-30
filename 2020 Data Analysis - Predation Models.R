#Whatcom TAST - 2020 Analysis - Foraging Success/Predation 
#Script written by Kathleen McKeegan
#Last update made 9/30/2023

## Working Directory: 

setwd("~/Desktop/TAST-  Data Cache and R Scripts/Data - csv")

## Load Packages
library('tidyverse')
library('glmmTMB')
library('DHARMa')
library('performance') #for overdispersion and zero inflation
library('sjPlot') #for random effects plots
library('emmeans')


#### Upload and format 2020 TAST data ####
tast.dat <- read.csv("2020_data.csv", header=T)
tast.dat$Fish<- sub("^$", "0", tast.dat$Fish) #fill blanks with zero
tast.dat$Fish[tast.dat$Fish == "Y"] <- "1" #replace 'Y' with one
tast.dat$Fish <- as.numeric(tast.dat$Fish)
tast.dat$Date <- as.Date(tast.dat$Date, "%m/%d/%y")
tast.dat$ID <- str_pad(tast.dat$ID, 4, pad = "0") #make all IDs a 4-digit number
tast.dat$ID[tast.dat$ID == "Sea_lion"]  <- "sea_lion" #fix typo in sea_lion 
tast.dat$ID[tast.dat$ID == "00CF"] <- "0198" #fix incorrect ID
tast.dat$ID[tast.dat$ID == "0265"] <- "0072" #fix incorrect ID
tast.dat$ID <- as.factor(tast.dat$ID)
tast.dat$Salmon3day <- as.character(tast.dat$Salmon3day)
tast.dat$Salmon3day <- as.numeric(tast.dat$Salmon3day)
tast.dat$Salmon5day <- as.character(tast.dat$Salmon5day)
tast.dat$Salmon5day <- as.numeric(tast.dat$Salmon5day)
tast.dat$TAST <- as.factor(tast.dat$TAST)

#Remove any RS or F from the New IDs, keep only LS since we know they are unique individuals
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

# Add presence column for summaries (see GLMM analysis - old models)
tast.summary <- tast.summary %>%
  mutate(Presence = 1)

# Get rid of NAs in Salmon rolling average
tast.summary[is.na(tast.summary)] <- 0

# Add column for observation number across TAST status (potential offset term)
tast.summary <- tast.summary %>%
  mutate(ObsDays = if_else( TAST == 'OFF', 14, 16))

# Add column with sum total number of seals seen per observation period
tast.summary <- tast.summary %>%
  group_by(Date) %>%
  mutate(Total=length(ID))

# Add column with Observation length in hours, not minutes
tast.summary <- tast.summary %>%
  mutate(LengthHR = (Length/60))

#### GLMM Analysis for Predation - New Models ####

#Explore the data a bit - how are these variables related?
# Predation and TAST
plot(tast.summary$TAST,tast.summary$Catches)
# Predation and Tide
plot(tast.summary$Tide_Height,tast.summary$Catches)
# Predation and Salmon
plot(tast.summary$Salmon5day,tast.summary$Catches)
# Predation and ID
plot(tast.summary$ID,tast.summary$Catches)

#What does the catch data look like? Lots of zeros, positive skew
hist(tast.summary$Catches)

# Response Variable is Foraging Success Counts -- need poisson or nb distribution
# Foraging Success counts = number of salmon consumed per ID during observation 

# Model Selection process: 
## Step 1 - optimal combination of random effects with fully populated fixed effect term constant
## Step 2 - keep selected random effects and find optimal fixed effect combination


### STEP 1
#
# Model 1 - TAST and tide as fixed predictors, ID as random, obs effort (Length) as offset
predmodel <- glmmTMB(Catches ~ TAST + Tide_Height + (1|ID), data=tast.summary,
                     offset=log(Length), family=poisson(link='log'))
check_overdispersion(predmodel) #no overdispersion detected, ok with Poisson
check_zeroinflation(predmodel) #model is ok, not zero inflated

#
# Model 2 - Date as Random, obs effort (Length) as offset
predmodel2 <- glmmTMB(Catches ~ TAST + Tide_Height + (1|Date), data=tast.summary,
                      offset=log(Length), family=nbinom2(link='log'))

#
# Model 3 - Number of Cameras (observers) as Random, obs effort (Length) as offset
predmodel3 <- glmmTMB(Catches ~ TAST + Tide_Height + (1|Camera), data=tast.summary,
                      offset=log(Length), family=nbinom2(link='log'))

#
# Model 4 - ID and Date as random
predmodel4 <- glmmTMB(Catches ~ TAST + Tide_Height+  (1|ID) + (1|Date), 
                      data=tast.summary,offset=log(Length), family=poisson(link='log'))

#
# Model 5 - ID and Camera as random
predmodel5 <- glmmTMB(Catches ~ TAST + Tide_Height+ (1|ID) + (1|Camera), 
                      data=tast.summary,offset=log(Length), family=poisson(link='log'))

#
# Model 6 - ID and Date and Camera as random 
predmodel6 <- glmmTMB(Catches ~ TAST + Tide_Height+ (1|ID) + (1|Date) + (1|Camera), 
                      data=tast.summary,offset=log(Length), family=poisson(link='log'))


#Compare models to select random term
AIC(predmodel, predmodel2, predmodel3, predmodel4, predmodel5, predmodel6) 
#First model with just ID as random is best


### STEP 2 - Select combination of Fixed Terms with ID as only random term 

# Model 4 - TAST only, ID as random
predmodel7 <- glmmTMB(Catches ~ TAST + (1|ID), data=tast.summary,
                      offset=log(Length), family=poisson(link='log'))

# Compare model 1 with model 7
AIC(predmodel, predmodel7)
#Model with tide is slightly better, according to AIC

### predmodel: TAST + Tide + (1|ID)
summary(predmodel)

# Model 1 diagnostics:
check_overdispersion(predmodel) #good
check_zeroinflation(predmodel) #good
testDispersion(predmodel)
simulationOutput1 <- simulateResiduals(fittedModel=predmodel)
plot(simulationOutput1) #everything looks good

# Plot Model - quick visualization
plot_model(predmodel, type='pred') 
#Tide is positively associated with catch
#TAST confidence intervals overlap

#Nicer plot
ob_p <- emmeans(predmodel, ~TAST, type='response')
ob_p
plot(ob_p, xlab='Model Estimates: Duration (minute counts/hour)') +
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 13))

#exponentiated estimates and confidence intervals 
exp(confint(predmodel))

# Plot final model random effects, look how catch varies across seal ID 
plot_model(predmodel, type='re', show.values = TRUE, value.offset = .5)



#### GLMM Analysis for Catches - Smaller sample size ####

# Only include seal IDs that were present when TAST was on and off
# Allow for comparison in foraging success across experimental and control conditions

# Seals observed with TAST on and off:
# 0012 0025 0039 0062 0083 0084 0085 00AC 0102 0104 0115 0117 0121 0122 0124 0137 0145 0154 0166 0167 0168 0170 0172 0173 0184 0186 0198 0200 0213 0217 0218 0219 0220 0221 0222 0223 0225 0226 0227 0228 0229 0230 0231 0232 0234 0235 0238 0240 0242 0243 0244 0245 0246 0247 0248 #

#Subset data to the 55 seals seen across TAST status
tast.summary2 <- tast.summary %>%
  filter(ID == "0012" | ID == "0025" | ID == "0039"| ID == "0062"| ID == "0083"| ID == "0084"| ID == "0085"| ID == "00AC"| ID == "0102"| ID == "0104"| ID == "0115"| ID == "0117"| ID == "0121"| ID == "0122"| ID == "0124"| ID == "0137"| ID == "0145"| ID == "0154"| ID == "0166"| ID == "0167"| ID == "0168"| ID == "0170"| ID == "0172"| ID == "0173"| ID == "0184"| ID == "0186"| ID == "0198"| ID == "0200"| ID == "0213"| ID == "0217"| ID == "0218"| ID == "0219"| ID == "0220"| ID == "0221"| ID == "0222"| ID == "0223"| ID == "0225"| ID == "0226"| ID == "0227"| ID == "0228"| ID == "0229"| ID == "0230"| ID == "0231"| ID == "0232"| ID == "0234"| ID == "0235"| ID == "0238"| ID == "0240"| ID == "0242"| ID == "0243"| ID == "0244"| ID == "0245"| ID == "0246"| ID == "0247"| ID == "0248")

# Model 8 - Same model but with smaller sample size:
predmodel8 <- glmmTMB(Catches ~ TAST + Tide_Height + (1|ID), data=tast.summary2,
                      offset=log(LengthHR), family=poisson(link='log'))
summary(predmodel8)

# Model 8 diagnostics:
check_overdispersion(predmodel8) #no overdispersion detected, ok with Poisson
check_zeroinflation(predmodel8) #good
testDispersion(predmodel8)
simulationOutput8 <- simulateResiduals(fittedModel=predmodel8)
plot(simulationOutput8) #good
plotResiduals(simulationOutput8, tast.summary2$TAST) #good
plotResiduals(simulationOutput8, tast.summary2$Tide_Height) #good

# Plot Model8
plot_model(predmodel8, type='pred') #confidence intervals overlap

ob <- emmeans(predmodel8, ~TAST, type='response')
ob
plot(ob, xlab='Model Estimates: Predation Rate (successes/hour)') +
  theme_classic()+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 13))

#exponentiated estimates and confidence intervals 
exp(confint(predmodel8))

# Plot final model random effects, see how catches vary across seal ID
p<-plot_model(predmodel8, type='re', show.values = TRUE, value.offset = .6)
p+theme_sjplot() 
