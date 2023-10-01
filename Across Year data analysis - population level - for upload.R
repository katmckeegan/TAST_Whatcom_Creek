#Whatcom TAST - Across Years Analysis - Population Level 
#Script written by Kathleen McKeegan
#Last update made 9/30/2023

## Load Packages 
library('tidyverse')
library('performance')
library('glmmTMB')
library('DHARMa')
library('emmeans')
library('sjPlot')
library('egg')

## Load Data
yrs.pop <- read.csv("pop-acrossyrs.csv", header=T)
yrs.pop$Date <- as.Date(yrs.pop$Date,format = "%m/%d/%y")
yrs.pop$year <- as.factor(yrs.pop$year)
yrs.pop$Tast <- as.factor(yrs.pop$Tast)
yrs.pop$Weather[yrs.pop$Weather == "rainy"]  <- "rain"
yrs.pop$Weather[yrs.pop$Weather == "Sunny"]  <- "sunny"
#Change 10/27/20 to tast on (mistake in data)
yrs.pop$Tast[which(yrs.pop$Date == "2020-10-27")] <- 'ON'

#Create Experimental Windows for Population Level Data
yrs.pop1 <- yrs.pop %>%
  filter(Date >= '2019-10-29' & Date <= '2019-11-29' |
           Date >= '2020-10-29' & Date <= '2020-11-29'|
           Date >= '2021-10-30' & Date <= '2021-11-30')

yrs.pop1 %>% 
  group_by(Tast) %>%
  summarise(num_rows = length(Tast))

yrs.pop1 <- yrs.pop1 %>%
  group_by(Tast) %>%
  mutate(ObsDays = length(Tast))

#### Descriptive Stats ####

yrs.pop19 <- yrs.pop1 %>%
  filter(year == '2019')
yrs.pop20on <- yrs.pop1 %>%
  filter(Tast == 'ON')
yrs.pop20off <- yrs.pop1 %>%
  filter(Tast == 'OFF')
yrs.pop21 <- yrs.pop1 %>%
  filter(year == '2021')

mean(yrs.pop19$Seals) #10.15
sd(yrs.pop19$Seals) #4.54
mean(yrs.pop20on$Seals) #9.2
sd(yrs.pop20on$Seals) #5.31
mean(yrs.pop20off$Seals) #17
sd(yrs.pop20off$Seals) #10.24
mean(yrs.pop21$Seals) #12.37
sd(yrs.pop21$Seals) #7.75

mean(yrs.pop19$Total_Catches) #1.46
sd(yrs.pop19$Total_Catches) #1.45
mean(yrs.pop20on$Total_Catches) #2
sd(yrs.pop20on$Total_Catches) #1.89
mean(yrs.pop20off$Total_Catches) #4.7
sd(yrs.pop20off$Total_Catches) #3.5
mean(yrs.pop21$Total_Catches) #3.3
sd(yrs.pop21$Total_Catches) #4.21

#### Model number of catches: ####

hist(yrs.pop1$Total_Catches) #right skew, counts - Poisson or neg binom

#Step 1 - select random effects, keeping fixed effects constant
# Make it so the reference is TAST OFF in 2020
tast.order <- c('OFF', 'ON', 'before', 'after')
yrs.pop1 = yrs.pop1 %>%
  mutate(Tast = factor(Tast, level=tast.order))

#### Model 1 - camera as random (significant deviation with poisson, select nb)
glm.pop.catch <- glmmTMB(Total_Catches ~ Tast + Roll_5 + Tide_height + 
                           (1|Camera), 
                         data=yrs.pop1, family = nbinom2(link='log'))

### Model 2 -  Year as random
glm.pop.catch2 <- glmmTMB(Total_Catches ~ Tast + Roll_5 + Tide_height + (1|year), 
                          data=yrs.pop1, family = nbinom2(link='log'))

### Model 3 - Year and Camera 

glm.pop.catch3 <- glmmTMB(Total_Catches ~ Tast + Roll_5 + Tide_height +
                            (1|year) + (1|Camera), 
                          data=yrs.pop1, family = nbinom2(link='log'))

AIC(glm.pop.catch, glm.pop.catch2, glm.pop.catch3) #model 1 and 2 are the same, 3 is worse
#try without random effect, see if random is not necessary?

### Model 4 - no random, all fixed (FINAL MODEL)
glm.pop.catch4 <- glmmTMB(Total_Catches ~ Tast + Roll_5 + Tide_height, 
                          data=yrs.pop1, family = nbinom2(link='log'))

AIC(glm.pop.catch, glm.pop.catch2,glm.pop.catch3, glm.pop.catch4) #no random is slightly better
#Move forward with no random effects

# Step 2 - fit fixed effects 

### Model 5 - TAST and salmon
glm.pop.catch5 <- glmmTMB(Total_Catches ~ Tast + Roll_5, 
                          data=yrs.pop1, family = nbinom2(link='log'))

### Model 6 - TAST and tide 
glm.pop.catch6 <- glmmTMB(Total_Catches ~ Tast + Tide_height, 
                          data=yrs.pop1, family = nbinom2(link='log'))

### Model 7 - just TAST
glm.pop.catch7 <- glmmTMB(Total_Catches ~ Tast, 
                          data=yrs.pop1, family = nbinom2(link='log'))


AIC(glm.pop.catch4, glm.pop.catch5, glm.pop.catch6, glm.pop.catch7)
#model 4 with all fixed effects is best

# Final model - model 4
summary(glm.pop.catch4) 
testDispersion(glm.pop.catch4)
simulationOutput4 <- simulateResiduals(fittedModel=glm.pop.catch4)
plot(simulationOutput4) #everything looks good

exp(confint(glm.pop.catch4))
catch.est = emmeans(glm.pop.catch4, ~Tast, type='response')
catch.est
catch.p <- plot(catch.est, xlab='Model Estimates: Mean number of salmon caught per observation', 
                ylab='TAST Status') +
  theme_classic() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 13))
catch.p

plot_model(glm.pop.catch4, type='pred', show.values = TRUE, value.offset = .5)

# Plot with 50% CI
catch.est.5 <- summary(catch.est, level=0.5)
catch.est.5
catch.est.5$lower.CL.50 <- catch.est.5$lower.CL
catch.est.5$upper.CL.50 <- catch.est.5$upper.CL
catch.est.5 <- select(catch.est.5, -c("lower.CL", "upper.CL"))
catch.est.conf <- merge(catch.est, catch.est.5, by = c("Tast","response", "SE", "df"))

tast.order.plot <- c('before', 'ON', 'OFF', 'after')

catch.plot50 <- ggplot(catch.est.conf, aes(x=response, 
                                           y=factor(Tast, level=tast.order.plot)))+
  geom_errorbar(aes(xmin=lower.CL.50, xmax=upper.CL.50),
                color='lightblue', width=0, linewidth=5)+
  geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL),
                color='black', width=0.2)+
  geom_point(aes(x=response))+
  labs(x='Model Estimates - Mean number of catches per observation', 
       y= 'Year - TAST Status')+
  scale_y_discrete(labels=c("OFF" = "2020-OFF", "ON" = "2020-ON",
                            "before" = "2019-before", "after"="2021-after"))
catch.plot50


#### Model number of seals ####

hist(yrs.pop1$Seals) #closer to normal distribution, still slightly skewed. Counts

#Step 1 - select random effects, keeping fixed effects constant

#### Model 1 - camera as random (significant deviation with poisson, select nb)
glm.pop.seal <- glmmTMB(Seals ~ Tast + Roll_5 + Tide_height + (1|Camera), 
                        data=yrs.pop1, family = nbinom2(link='log'))

### Model 2 -  Year as random
glm.pop.seal2 <- glmmTMB(Seals ~ Tast + Roll_5 + Tide_height + (1|year), 
                         data=yrs.pop1, family = nbinom2(link='log'))

AIC(glm.pop.seal, glm.pop.seal2) #same AIC, compare to no random 

### Model 3 -  No random (FINAL MODEL)

glm.pop.seal3 <- glmmTMB(Seals ~ Tast + Roll_5 + Tide_height, 
                         data=yrs.pop1, family = nbinom2(link='log'))

AIC(glm.pop.seal, glm.pop.seal2, glm.pop.seal3) 
#same as catch models, fixed only is slightly better than random. Move forward with no random

# Step 2 - Fit fixed effects

### Model 4 - TAST and salmon 
glm.pop.seal4 <- glmmTMB(Seals ~ Tast + Roll_5, 
                         data=yrs.pop1, family = nbinom2(link='log'))

### Model 5 - TAST and tide
glm.pop.seal5 <- glmmTMB(Seals ~ Tast + Tide_height, 
                         data=yrs.pop1, family = nbinom2(link='log'))

### Model 6 - TAST only 
glm.pop.seal6 <- glmmTMB(Seals ~ Tast, 
                         data=yrs.pop1, family = nbinom2(link='log'))

AIC(glm.pop.seal3, glm.pop.seal4, glm.pop.seal5, glm.pop.seal6)
#Fully populated model has the lowest AIC (model3)

# Final Model - Model 3
summary(glm.pop.seal3) 
testDispersion(glm.pop.seal3)
simulationOutput3 <- simulateResiduals(fittedModel=glm.pop.seal3)
plot(simulationOutput3) #everything looks good
exp(confint(glm.pop.seal3))

seals.est = emmeans(glm.pop.seal3, ~Tast, type='response')
seals.est

seals.p <- plot(seals.est, xlab='Model Estimates: Mean number of seals per observation', 
                ylab='TAST Status') +
  theme_classic() +
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size = 13))
seals.p

ggarrange(seals.p, catch.p)

plot_model(glm.pop.seal3, type='pred', show.values = TRUE, value.offset = .5)

# Plot with 50% CI
seal.est.5 <- summary(seals.est, level=0.5)
seal.est.5
seal.est.5$lower.CL.50 <- seal.est.5$lower.CL
seal.est.5$upper.CL.50 <- seal.est.5$upper.CL
seal.est.5 <- select(seal.est.5, -c("lower.CL", "upper.CL"))
seal.est.conf <- merge(seals.est, seal.est.5, by = c("Tast","response", "SE", "df"))

seal.plot50 <- ggplot(seal.est.conf, aes(x=response, 
                                         y=factor(Tast, level=tast.order.plot)))+
  geom_errorbar(aes(xmin=lower.CL.50, xmax=upper.CL.50),
                color='lightblue', width=0, linewidth=5)+
  geom_errorbar(aes(xmin=lower.CL, xmax=upper.CL),
                color='black', width=0.2)+
  geom_point(aes(x=response))+
  labs(x='Model Estimates - Mean number of seals per observation', 
       y= 'Year - TAST Status')+
  scale_y_discrete(labels=c("OFF" = "2020-OFF", "ON" = "2020-ON",
                            "before" = "2019-before", "after"="2021-after"))
seal.plot50

ggarrange(seal.plot50, catch.plot50)
