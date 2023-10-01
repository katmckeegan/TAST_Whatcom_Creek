#### Short-Term Effects: Fall 2020 - Duration and Presence ###
# TAST Whatcom Creek Manuscript 
# Last update made by Kathleen McKeegan 9/15/2023

## Load Packages for Script
library('tidyverse')


######### 2020 Only - Effects of TAST in ON year ########

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

#Remove any RS from the New IDs, keep only LS since we know they are unique individuals

tast.dat <- subset(tast.dat, ID!= '00AO'& ID!='00BA'& ID!= '00BC'& ID!='00BD'& ID!='00BF'& ID!='00BG'& 
                     ID!='00BR'& ID!='00BT'& ID!='00BW'& ID!='00BY'& ID!='00BZ'& ID!='00CA'& ID!='00CJ'&
                     ID!= '00CL'& ID!='00CT'& ID!='00CU' & ID!= '00CW' &ID!="sea_lion")

## Create dataframe with summary of d and catches
#Group by duration (surface counts) of each ID per day, narrow it down to the experimental window 10/29-11/25

tast.summary <- tast.dat %>%
  group_by(Date, ID, Length, Camera, Tide_Height, TAST, Salmon5day, Salmon3day) %>%
  summarise(Occurrence = n_distinct(Time),
            Catches = sum(Fish)) %>%
  filter(Date >= '2020-10-25' & Date <= '2020-11-25')

###### Data Visualization and Exploration ###### 

#Summarize the information per ID per Day
tast.summary.day <- tast.summary %>%
  group_by(ID, TAST) %>%
  summarise(nDays=length(Date), 
            meanOccurrence = mean(Occurrence),
            sumOccurrence = sum(Occurrence),
            sumCatches = sum(Catches), 
            meanCatches=mean(Catches))

# Look at individuals seen based on TAST
summary.day2 <- tast.summary.day[(1:153), (1:3)]
summary.day2 <- pivot_wider(summary.day2, names_from= 'TAST', values_from="nDays")
summary.day2 <- summary.day2 %>% replace(is.na(.), 0)
sum(summary.day2$ON == 0) #Number of individuals seen only when TAST was OFF
31/98
sum(summary.day2$OFF == 0) #Number of individuals seen only when TAST was ON 
12/98

## quick visualization of duration (surface counts) per ID throughout season
table.occurrence <- tast.dat %>%
  count(ID)
#Summary of duration
#Max
table.occurrence[which.max(table.occurrence$n),] #ID 173 had greatest number of surface counts
#Min
mins.occurrence <- table.occurrence[table.occurrence$n == '1',]

# Number of individuals observed during the experimental window
n_distinct(tast.summary$ID) #98 individuals 

#### Duration across TAST Status #### 

duration.summary <- cbind(tast.summary.day[,c(1:5)])

#Average duration: 
duration.on <- duration.summary %>%
  filter(TAST == "ON")
duration.off <- duration.summary %>%
  filter(TAST == "OFF")
mean(duration.on$meanOccurrence) #6.33 min
sd(duration.on$meanOccurrence)
mean(duration.off$meanOccurrence) #9.67 min

Duration.v2<- ggplot(duration.summary, aes(x=TAST, y=meanOccurrence)) +
  geom_violin(fill='grey')+
  geom_boxplot(color='black', width=0.1, alpha=0.2)+
  labs(y='Mean Duration (min) per Observation')+
  geom_text(data=subset(duration.summary, ID=='0075' | ID=='0236' | ID=='0173'),
            aes(TAST, meanOccurrence, label=ID), nudge_x = 0.1, size=5)+
  theme_classic()+
  theme(legend.position='none',
        axis.text.x = element_text(colour = "black", size = 13, face = "bold"), 
        axis.text.y = element_text(colour = "black", size = 13),
        axis.title = element_text(colour="black", size=15))
Duration.v2
