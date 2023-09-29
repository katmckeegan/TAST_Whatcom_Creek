#### Short term data exploration - for predation 
## Working Directory: 

setwd("~/Desktop/TAST-  Data Cache and R Scripts/Data - csv")

## Load Packages for Script
library('tidyverse')

#Load Data
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

## Create dataframe with summary of occurence and catches
#Group by occurrences of each ID per day, narrow it down to the experimental window 10/29-11/25

tast.summary <- tast.dat %>%
  group_by(Date, ID, Length, Camera, Tide_Height, TAST, Salmon5day, Salmon3day) %>%
  summarise(Occurrence = n_distinct(Time),
            Catches = sum(Fish)) %>%
  filter(Date >= '2020-10-25' & Date <= '2020-11-25')

###### Data Visualization and Exploration ###### 

#Summarise the information per ID per Day
tast.summary.day <- tast.summary %>%
  group_by(ID, TAST) %>%
  summarise(nDays=length(Date), 
            meanOccurrence = mean(Occurrence),
            sumOccurrence = sum(Occurrence),
            sumCatches = sum(Catches), 
            meanCatches=mean(Catches))

## Visualize impact on predation successes
# Look at catches across TAST status per ID graphically 

catches.summary <- tast.summary.day %>%
  select(ID, TAST, sumCatches) %>%
  group_by(ID)%>%
  mutate(totalCatches = sum(sumCatches))

#Keep only seals that caught a fish during season
catches.ID <- catches.summary %>%
  filter(totalCatches > 0)

#Calculate the proportion of catches made by individuals when TAST was on
prop.catches <- catches.ID %>%
  add_column(Proportion = catches.ID$sumCatches/catches.ID$totalCatches) 

#Pull out just proportion data
prop.catches1 <- prop.catches %>%
  select(ID, TAST, Proportion) 

#Make a column for TAST on and TAST off
prop.catches2 <- pivot_wider(prop.catches1, names_from = "TAST",
                             values_from = "Proportion")
#Remove individuals who were not present under both conditions
prop.catches2 <- prop.catches2 %>%
  drop_na()

#What proportion of catches were made when TAST was on vs off? 
#AKA, with TAST on, did seals catch fewer salmon? 
prop.catches2 <- prop.catches2 %>%
  add_column(Difference = prop.catches2$ON - prop.catches2$OFF)

sum(prop.catches2$Difference == 1) #8 seals caught 100% more salmon when TAST was ON
sum(prop.catches2$Difference == -1) #7 seals caught 100% more salmon when TAST was off

#Visualize proportions 
#This plot shows how successful individuals were when TAST was on/off
# Positive values indicate seal was more successful when TAST was on
# Negative values indicate seal was more successful when TAST was off
# +1.0 means the seal only ever caught fish when TAST was On
# -1.0 means seal only ever caught fish when TAST was Off
# 0.0 means seal caught same number of fish when On or Off 
prop.plot <- ggplot(prop.catches2, aes(x=ID, y=Difference, label=ID))+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "darkgray", size=0.5)+
  geom_label(label.padding = unit(0.3, "lines"))+
  labs(x="ID of Seal", y="Proportion of Successes with TAST On")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size=25),
        axis.text.y = element_text(size = 20))
prop.plot
