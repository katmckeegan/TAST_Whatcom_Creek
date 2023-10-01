#Whatcom TAST - Comparison Across Years - Presence 2019 to 2021
#Script written by Kathleen McKeegan
#Last update made 9/30/2023

## Working Directory: 

setwd("~/Desktop/TAST-  Data Cache and R Scripts/Data - csv")

## Load Packages for Script
library('tidyverse')
library('performance') #to check over-dispersion and zero-inflation
library('rsq')
library('MASS') #for glm.nb
library('sjPlot')

#### Compare 2020 presence to 2021 presence based on TAST exposure ####

# Part 1 - Define Exposure to TAST and determine Exposed seals ###

# Load 2020 data 
tast.dat <- read.csv("2020_data.csv", header=T)
tast.dat$Fish<- sub("^$", "0", tast.dat$Fish)
tast.dat$Fish[tast.dat$Fish == "Y"] <- "1"
tast.dat$Fish <- as.numeric(tast.dat$Fish)
tast.dat$Date <- as.Date(tast.dat$Date, "%m/%d/%y")
tast.dat$ID <- str_pad(tast.dat$ID, 4, pad = "0") #make all IDs a 4-digit number
tast.dat$ID[tast.dat$ID == "Sea_lion"]  <- "sea_lion" #fix typo in sea_lion impacting summary
tast.dat$ID[tast.dat$ID == "00CF"] <- "0198" #fix incorrect ID
tast.dat$ID[tast.dat$ID == "0265"] <- "0072" #fix incorrect ID
tast.dat$ID <- as.factor(tast.dat$ID)
tast.dat$Salmon3day <- as.character(tast.dat$Salmon3day) #Add salmon count data
tast.dat$Salmon3day <- as.numeric(tast.dat$Salmon3day)
tast.dat$Salmon5day <- as.character(tast.dat$Salmon5day)
tast.dat$Salmon5day <- as.numeric(tast.dat$Salmon5day)
tast.dat$TAST <- as.factor(tast.dat$TAST)

#Remove any RS from the New IDs, keep only LS since we know they are unique individuals

tast.dat <- subset(tast.dat, ID!= '00AO'& ID!='00BA'& ID!= '00BC'& ID!='00BD'& ID!='00BF'& ID!='00BG'& 
                     ID!='00BR'& ID!='00BT'& ID!='00BW'& ID!='00BY'& ID!='00BZ'& ID!='00CA'& ID!='00CJ'&
                     ID!= '00CL'& ID!='00CT'& ID!='00CU' & ID!= '00CW' &ID!="sea_lion")

# Create dataframe with summary of occurence and catches
## Group by occurrences of each ID per day, narrow it down to the experimental window 10/25-11/25

tast.summary <- tast.dat %>%
  group_by(Date, ID, Length, Camera, Tide_Height, TAST, Salmon5day, Salmon3day) %>%
  summarise(Occurrence = n_distinct(Time),
            Catches = sum(Fish)) %>%
  filter(Date >= '2020-10-25' & Date <= '2020-11-25')


# Create data frame where presence is defined by 5 or more 'surface counts' per observation period

tast.summary2 <- subset(tast.summary, Occurrence>=5)
length(unique(tast.summary$ID)) #98 unique seals identified in 2020 experimental window
length(unique(tast.summary2$ID)) #62 unique seals present for at least 5 surface counts on at least one day

# Create data frame summarizing presence of IDs across TAST 

ID.summary.exp <- tast.summary2 %>%
  group_by(ID, TAST) %>%
  summarise(nDays = length(Date),
            sumCatch = sum(Catches))

ID.summary.exp$TAST <- as.factor(ID.summary.exp$TAST)

# Data frame showing presence for each ID across TAST status (with presence >= 5 surface counts per observation period)
ID.summary.exp <- ID.summary.exp %>%
  group_by(ID) %>%
  complete(TAST, fill = list(nDays = 0))


# Add column for binary yes/no exposed to TAST in 2020

ID.expbi <- ID.summary.exp %>%
  mutate(Exposure2020 = case_when(TAST == 'ON' & nDays > 0 ~ '1',
                                  TAST == 'ON' & nDays == 0 ~ '0',
                                  TAST == 'OFF' ~ '0'))

# Rename nDays column to help with future merging 
ID.expbi <- ID.expbi %>% 
  rename(nDays2020 = nDays)

# Part 2: Determine which seals returned in 2021 ###

yrs.dat <- read.csv("across_years_LabOnly.csv", header=T)
yrs.dat$Date <- as.Date(yrs.dat$Date,format = "%m/%d/%y")
yrs.dat$ID <- str_pad(yrs.dat$ID, 4, pad = "0") #make all IDs a 4-digit number


yrs.dat <- subset(yrs.dat, ID!= '00AO'& ID!='00BA'& ID!= '00BC'& ID!='00BD'& ID!='00BF'& ID!='00BG'& 
                    ID!='00BR'& ID!='00BT'& ID!='00BW'& ID!='00BY'& ID!='00BZ'& ID!='00CA'& ID!='00CJ'&
                    ID!= '00CL'& ID!='00CT'& ID!='00CU' & ID!= '00CW' & ID!='00DI' & ID!="sea_lion")

# Isolate just 2021
yrs.dat$year <- year(ymd(yrs.dat$Date))
yrs.dat$year <- as.factor(yrs.dat$year)
yrs.dat21 <- yrs.dat %>%
  filter(year=="2021")


yrs21.summary <- yrs.dat21 %>%
  group_by(ID) %>%
  summarise(nDays = length(Date))

# Rename nDays column to help with future merging 
yrs21.summary <- yrs21.summary %>% 
  rename(nDays2021 = nDays)

# Part 3: Create Data frame showing Exposure 2020 and Presence 2021 ###

# use ID.expbi dataframe, includes all 62 IDs present in 2020 for more than 5 surface counts
# has the binary yes/no for exposure to TAST, as well a sum Days

full.exp <- merge(ID.expbi, yrs21.summary, by="ID", all.x = TRUE)
full.exp[is.na(full.exp)] <- 0

# Add column for binary yes(1)/no(0) for presence in 2021
full.exp <- full.exp %>%
  mutate(Present2021 = if_else(nDays2021 == 0, "0", "1"))


# Part 4: Run GLM with binomial distribution to assess relationship between exposure and presence ###
full.exp$Exposure2020 <- as.factor(full.exp$Exposure2020)
full.exp$Present2021 <- as.numeric(full.exp$Present2021)


# Model - with count exposure predicting binary presence interacting with treatment
glm <- glm(Present2021 ~ nDays2020*TAST, data=full.exp,
             family=binomial(link='logit'))
summary(glm)
plot(glm)
rsq(glm) 

c <- coef(glm)
exp(c) #exponentiate the coefficients from final model
exp(confint(glm))

plot_model(glm, type='int',
           axis.title=c('Days observed in 2020','Likelihood of Presence in 2021'),
           title='Predicted Probabilities for Individual Presence in 2021')


#### Compare 2019 presence to 2020 presence ####

# Part 1: Isolate 2019 Data ###

# Isolate just 2019
yrs.dat19 <- yrs.dat %>%
  filter(year=="2019")


yrs19.summary <- yrs.dat19 %>%
  group_by(ID) %>%
  summarise(nDays = length(Date))

# Rename nDays column to help with future merging 
yrs19.summary <- yrs19.summary %>% 
  rename(nDays2019 = nDays)

# Part 2 - Determine which seals were present in 2020 ###

# Create data frame summarizing presence of IDs, not across TAST

ID.summary <- tast.summary %>%
  group_by(ID) %>%
  summarise(nDays = length(Date))


# Rename nDays column to help with future merging 
ID.summary <- ID.summary %>% 
  rename(nDays2020 = nDays)

# Part 3: Create Data frame showing days in 2019 as it relates to 2020 ###

#full.exp <- merge(yrs19.summary, ID.summary, by="ID", all.x = TRUE)
full.summary <- merge(yrs19.summary, ID.summary, by="ID", all = TRUE)
full.summary[is.na(full.summary)] <- 0

# Add column for binary yes(1)/no(0) for presence in 2021
full.summary <- full.summary %>%
  mutate(Present2020 = if_else(nDays2020 == 0, "0", "1"))


# Part 3: Run GLM to assess relationship presence across years ###
full.summary$ID <- as.factor(full.summary$ID)
full.summary$Present2020 <- as.numeric(full.summary$Present2020)


#Model 2 - days present in 2019 related to days present in 2020
# over-dispersion detected with Poisson, use neg binomial 
glm2 <- glm.nb(nDays2020 ~ nDays2019, data=full.summary)
summary(glm2)
plot(glm2)
plot_model(glm2, type='pred', 
           title="Predicted number of days observed in 2020",
           axis.title=c('Days observed in 2019', 'Days observed in 2020'))

c <- coef(glm2)
exp(c) #exponentiate the coefficients from final model
exp(confint(glm2))
check_overdispersion(glm2)
