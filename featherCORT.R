#Social relationships and fCORT in jackdaws
#Author: Luca Hahn
#Last update: 10/07/2026

#(01) DATA CLEANING ---- 

#(1) First cleaning, skip to (2) once these steps are done 

install.packages("ggiraphExtra")
install.packages("mice")
install.packages("miceadds")
library(asnipe)
library(assortnet)
library(bbmle)
library(brms)
library(car)
library(carData)
library(chisq.posthoc.test)
library(ClusterR)
library(corrplot)
library(data.table)
library(DHARMa)
library(dplyr)
library(effects)
library(emmeans)
library(ggeffects)
library(ggiraphExtra)
library(ggplot2)
library(glmmTMB)
library(hms)
library(igraph)
library(interactions)
library(lme4)
library(lmerTest)
library(lubridate)
library(MASS)
library(mice)
library(miceadds)
library(multcomp)
library(MuMIn)
library(loo)
library(performance)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(rptR)
library(sjlabelled)
install.packages('TMB', type = 'source')
install.packages("sjPlot")
library(sjPlot)
library(sjstats)
library(sjmisc)
library(stringi)
library(stringr)
library(svMisc)
library(tidybayes)
library(tidyr)
library(tidyverse)
library(assortnet)

library(tidyverse)
library(emmeans)
library(ggplot2)
library(lme4)
library(lmerTest)
library(Matrix)
library(rptR)
install.packages("datawizard")
library(datawizard)
library(gridGraphics)
library(grid)
library(gridBase)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(ggeffects)
library(rptR)
library(sjlabelled)
library(sjPlot)

library(extrafont)
font_import()
loadfonts()

library(corrr)
library(FactoMineR) 
library(factoextra)
library(psych)

library(sna)

#Read feather_CORT data
fCORT <- read.csv("Data/feather_CORT.csv", header = T, stringsAsFactors = F)

fCORT_perindiv_peryr <- as.data.frame(table(fCORT$feather_year_grown, fCORT$JID))

associations <- read.csv("C:/Users/lh868/Downloads/NeighboursAndKinData.csv", header = T, stringsAsFactors = F)

#Read life history file 
LH <- read.csv("Data/LH20231004.csv", header = T, stringsAsFactors = F)
LH$DATE <- strptime(LH$DATE,format="%d/%m/%Y")
LH$DATE <- as.Date(LH$DATE, format = "%d/%m/%Y") # convert to date

#Add JID_year in LH
LH$year <- as.numeric(format(LH$DATE,"%Y"))   #  #Get year from the date

LH$JID_year_before <- paste(LH$ID, LH$year, sep = "_")
LH$JID_year_grown <- paste(LH$ID, LH$year, sep = "_")
LH$JID_year_sampled <- paste(LH$ID, LH$year, sep = "_")

#LH feather 
LH_feather <- subset(LH, LH$FEATHER == TRUE)
LH_feather <- LH_feather %>%
group_by(JID_year_sampled) %>% 
  slice(1:1)

LH_tarsus <- subset(LH, !LH$TARSUS == "")
LH_tarsus <- subset(LH_tarsus, LH_tarsus$AGE > 3)
LH_tarsus <- LH_tarsus[, c("ID", "TARSUS")]
LH_tarsus$fCORT <- fCORT$CORT_pg.mg[match(LH_tarsus$ID, fCORT$JID)]
LH_tarsus <- subset(LH_tarsus, !is.na(LH_tarsus$fCORT))

tarsus <- as.data.frame(aggregate(x = LH_tarsus$TARSUS, by = list(LH_tarsus$ID), FUN = "mean"))
tarsus$JID <- tarsus$Group.1
tarsus$mean_tarsus <- tarsus$x
tarsus$Group.1 <- NULL
tarsus$x <- NULL

LH_weight <- subset(LH, !LH$WEIGHT == "")
LH_weight <- subset(LH_weight, LH_weight$AGE > 3)
LH_weight <- LH_weight[, c("ID", "WEIGHT")]
LH_weight$fCORT <- fCORT$CORT_pg.mg[match(LH_weight$ID, fCORT$JID)]
LH_weight <- subset(LH_weight, !is.na(LH_weight$fCORT))

weight <- as.data.frame(aggregate(x = LH_weight$WEIGHT, by = list(LH_weight$ID), FUN = "mean"))
weight$JID <- weight$Group.1
weight$mean_weight <- weight$x
weight$Group.1 <- NULL
weight$x <- NULL

LH_sub <- LH

#Separate LH file with (1) pairs and then (2) adding pair ID
#(i) LH file with pairs and boxes
LH_pairs <- subset(LH[c("DATE", "ID", "SEX", "BOX", "PARTNER.ID", "JID_year_before","JID_year_grown", "JID_year_sampled")])
LH_pairs <- LH_pairs[LH_pairs$PARTNER.ID != "",]

#(ii) Adding pair ID
LH_pairs  <- LH_pairs %>% 
  as_tibble() %>% 
  mutate(pair_ID = if_else(LH_pairs$SEX == "F", paste(ID, PARTNER.ID), paste(PARTNER.ID, ID)))

#Only to get most recent entries for pair ID 
#LH_pairs <- LH_pairs %>%
group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

#LH paired and repaired (different from LH pairs)
LH_paired <- subset(LH[c("DATE", "ID", "SEX", "BOX", "CODE", "PARTNER.ID", "JID_year_before","JID_year_grown", "JID_year_sampled")])
LH_paired <- LH_paired[LH_paired$CODE == "PAIRED",]
LH_paired$YEAR <- substr(LH_paired$DATE,1,4)

LH_repaired <- subset(LH[c("DATE", "ID", "SEX", "BOX", "CODE", "PARTNER.ID", "JID_year_before","JID_year_grown", "JID_year_sampled")])
LH_repaired <- LH_repaired[LH_repaired$CODE == "REPAIR",]
LH_repaired$YEAR <- substr(LH_repaired$DATE,1,4)

LH_paired_repaired <- rbind(LH_paired, LH_repaired)
LH_paired_repaired$year <- substr(LH_paired_repaired$YEAR, 3,4)
LH_paired_repaired$pair_ID <- ifelse(LH_paired_repaired$SEX == "F", paste(LH_paired_repaired$ID, LH_paired_repaired$PARTNER.ID, sep = " "), paste(LH_paired_repaired$PARTNER.ID, LH_paired_repaired$ID, sep = " "))

#Read LH_EARLY data
LH_early <- read.csv("Data/LH_EARLY.csv", header = T, stringsAsFactors = F)
LH_early$box_year <- LH_early$BOX_YEAR

#LH only including adults
LH_adults <- subset(LH, LH$AGE == 6)
LH_adults <- subset(LH_adults, !is.na(LH_adults$WEIGHT))
LH_adults <- subset(LH_adults, !is.na(LH_adults$TARSUS))

LH_adults$body_cond <- resid(lm(WEIGHT ~ TARSUS, data = LH_adults))

#Read box densities data
box_densities <- read.csv("Data/box_densities.csv", header = T, stringsAsFactors = F)
box_densities$year <- paste("20",box_densities$Year, sep = "")
box_densities$year<- as.numeric(box_densities$year)
box_densities$box_grown_year <- paste(box_densities$Nest, box_densities$year, sep = "_")

#Read breeding summary data
breeding_summary <- read.csv("Data/BreedingSummary2013-2023.csv", header = T, stringsAsFactors = F)

breeding_summary$pair_ID <- paste(breeding_summary$fem.ID, breeding_summary$male.ID, sep = " ")

breeding_summary2 <- breeding_summary  %>% pivot_longer(
  cols = c("fem.ID", "male.ID"))

breeding_summary2$JID <- breeding_summary2$value

breeding_summary2$JID_year_grown <- paste(breeding_summary2$JID, breeding_summary2$year, sep = "_")
breeding_summary2$JID_year_sampled <- paste(breeding_summary2$JID, breeding_summary2$year, sep = "_")

#Read pair bond strength data
pbvid1 <- read.csv("Data/pbStrength.csv", header = T, stringsAsFactors = F)
pbvid2 <- read.csv("Data/time_budget_summary.csv", header = T, stringsAsFactors = F)
pbvid3 <- read.csv("Data/pairbond.csv", header = T, stringsAsFactors = F)

pbvid_length1 <- read.csv("Data/pbStrength.csv", header = T, stringsAsFactors = F)
pbvid_length2 <- read.csv("Data/time_budget_summary.csv", header = T, stringsAsFactors = F)
pbvid_length3 <- read.csv("Data/pairbond.csv", header = T, stringsAsFactors = F)

sum(pbvid_length1$observation_minus_latency/3600)
sum(pbvid_length2$Vid.duration..s./3600)

pbvid_length1$obs_length <- pbvid_length1$observation_minus_latency/3600
pbvid_length2$obs_length <- pbvid_length2$Vid.duration..s./3600
pbvid_length2 <- subset(pbvid_length2, Behavior == "Both")
pbvid_length3$obs_length <- pbvid_length3$OBSERVATION_LENGTH/3600

obs_length <- c(pbvid_length1$obs_length, pbvid_length2$obs_length, pbvid_length3$obs_length)
mean(obs_length)
sd(obs_length)

mean(pbvid_length1$observation_minus_latency/3600)
mean(pbvid_length2$Vid.duration..s./3600)

#Read pair bond data 22 
pb22_pair <- read.csv("Data/pb22_pair.csv", header = T, stringsAsFactors = F)

#Read pair bond data 22 long 
pb_long_pair <- read.csv("Data/pb_long_pair.csv", header = T, stringsAsFactors = F)
pb22_long_pair <- subset(pb_long_pair, pb_long_pair$year == 2022)

pb22_long_pair$fCORT <- fCORT$CORT_pg.mg[match(pb22_long_pair$JID, fCORT$JID)]

pb22_long_pair$mean_strength_pb <- (pb22_long_pair$strength - pb22_long_pair$pbfeedgmm) / (pb22_long_pair$degree -1)

#Read box visit data 
box_visits_22 <- read.csv("Data/box_visits22_summary.csv", header = T, stringsAsFactors = F)

#Read temperament data 
temperament <- read.csv("Data/temperament.csv", header = T, stringsAsFactors = F)

#Additional information

#Test JID is correct (then reload fCORT again)
fCORT$JID_2 <- LH$ID[match(fCORT$BTO, LH$BTO)]
fCORT$JID_same <- ifelse(fCORT$JID == fCORT$JID_2, "same", "different")
table(fCORT$JID_same)

#Add site
fCORT[substr(fCORT$box,1,1)=="X","site"]<-"X"
fCORT[substr(fCORT$box,1,1)=="Y","site"]<-"Y"
fCORT[substr(fCORT$box,1,1)=="Z","site"]<-"Z"
fCORT[substr(fCORT$box,1,1)=="M","site"]<-"MW"
#two birds with no boxes are presumably resident at Z due to prospecting visits
fCORT[substr(fCORT$box,1,1)=="N","site"]<- "Z"

#Add sex
LH_sex <- LH[LH$SEX != "",]
fCORT$sex <- LH_sex$SEX[match(fCORT$JID, LH_sex$ID)]
fCORT$sex[fCORT$JID == "J5319"] <- "M"

fCORT <- fCORT %>%
  mutate(sex_num = case_when(
    sex =="F" ~ 0,
    sex =="M" ~ 1,    
    is.na(sex) ~ 0.5))

#Add temperament
fCORT$bite <- temperament$bite[match(fCORT$JID, temperament$JID)]

#Add box_year
fCORT$box_year <- paste(fCORT$box, fCORT$feather_year, sep = "_")

#Add year feather was grown
fCORT$feather_year_grown <- fCORT$feather_year_sampled - 1
fCORT$feather_year_grown <- as.numeric(fCORT$feather_year_grown)

#Add year before feather was grown 
fCORT$feather_year_before <- fCORT$feather_year_sampled - 2

#Add age
fCORT$year_diff <- 2022 - fCORT$feather_year_grown
fCORT$age_2022 <- Ringed$min_age[match(fCORT$JID, Ringed$ID)]
fCORT$age <- fCORT$age_2022 - fCORT$year_diff

fCORT$age[fCORT$age == 0] <- 1

#Add JID_year_before
fCORT$JID_year_before <- paste(fCORT$JID, fCORT$feather_year_grown -1, sep = "_")

#Add JID_year_grown
fCORT$JID_year_grown <- paste(fCORT$JID, fCORT$feather_year_grown, sep = "_")

#Add JID_year_sampled
fCORT$JID_year_sampled <- paste(fCORT$JID, fCORT$feather_year_sampled, sep = "_")

#Add tarsus length and weight year before
fCORT$tarsus_year_before <- LH$TARSUS[match(fCORT$JID_year_before, LH$JID_year_before)]
fCORT$weight_year_before <- LH$WEIGHT[match(fCORT$JID_year_before, LH$JID_year_before)]
length(table(fCORT$weight_year_before))

#Add tarsus length and weight year grown
fCORT$tarsus_year_grown <- LH$TARSUS[match(fCORT$JID_year_grown, LH$JID_year_grown)]
fCORT$weight_year_grown <- LH$WEIGHT[match(fCORT$JID_year_grown, LH$JID_year_grown)]
length(table(fCORT$weight_year_grown))

#Add tarsus length and weight year sampled 
fCORT$tarsus_year_sampled <- LH$TARSUS[match(fCORT$JID_year_sampled, LH$JID_year_sampled)]
fCORT$weight_year_sampled <- LH$WEIGHT[match(fCORT$JID_year_sampled, LH$JID_year_sampled)]
length(table(fCORT$weight_year_sampled))

fCORT$tarsus_year_sampled <- LH_feather$TARSUS[match(fCORT$JID_year_sampled, LH_feather$JID_year_sampled)]
fCORT$weight_year_sampled <- LH_feather$WEIGHT[match(fCORT$JID_year_sampled, LH_feather$JID_year_sampled)]
length(table(fCORT$weight_year_sampled))
length(table(fCORT$tarsus_year_sampled))

#Repeatability of weight
weight_rep <- rpt(WEIGHT ~ (1|ID), grname = c("ID"), data = LH, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(weight_rep)
plot(weight_rep)

#Add mean tarsus and weight
fCORT$weight_year_before[c(107, 116)] <- NA #change nonsense values to NA

fCORT <- fCORT %>% 
  mutate(tarsus =rowMeans(.[ , c("tarsus_year_before", "tarsus_year_grown", "tarsus_year_sampled")], na.rm=TRUE))

fCORT$mean_tarsus <- tarsus$mean_tarsus[match(fCORT$JID, tarsus$JID)]

fCORT <- fCORT %>% 
  mutate(weight=rowMeans(.[ , c("weight_year_before", "weight_year_grown", "weight_year_sampled")], na.rm=TRUE))

fCORT$mean_weight <- weight$mean_weight[match(fCORT$JID, weight$JID)]

length(table(fCORT$weight))

#Add body condition 
fCORT_body <- fCORT
fCORT_body <- subset(fCORT, !is.nan(fCORT$weight))

fCORT_body <- subset(fCORT_body, !is.na(fCORT$mean_weight))

fCORT_body$body_cond <- resid(lm(weight ~ tarsus, data = fCORT_body))
fCORT_body$mean_body_cond <- resid(lm(mean_weight ~ mean_tarsus, data = fCORT_body))

fCORT_body$body_cond_binary <- ifelse(fCORT_body$body_cond > 0 , "positive", "negative")

fCORT$body_cond <- fCORT_body$body_cond[match(fCORT$JID_year_grown,fCORT_body$JID_year_grown)]
fCORT$mean_body_cond <- fCORT_body$mean_body_cond[match(fCORT$JID,fCORT_body$JID)]

fCORT$body_cond_binary <- fCORT_body$body_cond_binary[match(fCORT$JID_year_grown,fCORT_body$JID_year_grown)]
fCORT$body_cond_binary_num <- ifelse(fCORT$body_cond_binary == "positive", 1, 0)

fCORT$body_cond_noNA <- ifelse(!is.na(fCORT$body_cond), fCORT$body_cond, fCORT$mean_body_cond)

#Add box year grown 
fCORT$box_grown <- breeding_summary2$box[match(fCORT$JID_year_grown, breeding_summary2$JID_year_grown)]
fCORT$box_grown_year <- paste(fCORT$box_grown, fCORT$feather_year_grown, sep = "_")

breeding_summary$box_grown_year <- paste(breeding_summary$box,breeding_summary$year, sep = "_")

#Add box owner yes or no (year feather grown)
fCORT$box_owner_grown <- ifelse(is.na(fCORT$box_grown), "non-owner", "owner")

#Add partner ID year before feather was grown
fCORT$partner_yr_before <- LH_pairs$PARTNER.ID[match(fCORT$JID_year_before, LH_pairs$JID_year_before)]

#Add partner ID year feather was grown
fCORT$partner_yr_grown <- LH_pairs$PARTNER.ID[match(fCORT$JID_year_grown, LH_pairs$JID_year_grown)]

#Add partner ID year feather was sampled
fCORT$partner_yr_sampled <- LH_pairs$PARTNER.ID[match(fCORT$JID_year_sampled, LH_pairs$JID_year_sampled)]

#Add pair ID year before feather was grown
fCORT$pair_ID_yr_before <- LH_pairs$pair_ID[match(fCORT$JID_year_before, LH_pairs$JID_year_before)]

#Add pair ID year feather was grown
fCORT$pair_ID <- LH_pairs$pair_ID[match(fCORT$JID_year_grown, LH_pairs$JID_year_grown)]

#Add pair ID_year feather was grown combined
fCORT$pair_ID_year <- paste(fCORT$pair_ID, fCORT$feather_year_grown, sep = "_")

#Add pair ID year feather was sampled
fCORT$pair_ID_yr_sampled <-  LH_pairs$pair_ID[match(fCORT$JID_year_sampled, LH_pairs$JID_year_sampled)]

fCORT$pair_ID_yr_sampled <- ifelse(!is.na(fCORT$pair_ID_yr_sampled), fCORT$pair_ID_yr_sampled, ifelse(fCORT$sex == "F", paste(fCORT$JID, "NA", sep = " "), paste("NA", fCORT$JID, sep = " ")))

#Add year individual paired
fCORT$paired_year <- LH_paired$YEAR[match(fCORT$JID, LH_paired$ID)]
fCORT$paired_year <- as.numeric(fCORT$paired_year)

#Add year individual repaired
fCORT$repaired_year <- LH_repaired$YEAR[match(fCORT$JID, LH_repaired$ID)]
fCORT$repaired_year <- as.numeric(fCORT$repaired_year)

#Add difference between year feather grown and paired
fCORT$feather_year_grown_paired <- fCORT$feather_year_grown - fCORT$paired_year

##Add difference between year feather grown and paired
fCORT$feather_year_grown_repaired <- fCORT$feather_year_grown - fCORT$repaired_year

#Add number of partners
fCORT <- fCORT %>%
  mutate(partner_samediff1 = case_when(partner_yr_before == partner_yr_grown ~ 1,
                                       partner_yr_before != partner_yr_grown ~ 2,
                                      is.na(partner_yr_before) | is.na(partner_yr_grown) ~ 1,
                                      is.na(partner_yr_before) & is.na(partner_yr_grown) ~ 1))

fCORT <- fCORT %>%
  mutate(partner_samediff2 = case_when(partner_yr_grown == partner_yr_sampled ~ 1,
                                       partner_yr_grown != partner_yr_sampled ~ 2,
                                       is.na(partner_yr_grown) | is.na(partner_yr_sampled) ~ 1,
                                       is.na(partner_yr_grown) & is.na(partner_yr_sampled) ~ 1))

fCORT <- fCORT %>%
  mutate(partner_samediff3 = case_when(partner_yr_before == partner_yr_sampled ~ 1,
                                       partner_yr_before != partner_yr_sampled ~ 2,
                                       is.na(partner_yr_before) | is.na(partner_yr_sampled) ~ 1,
                                       is.na(partner_yr_before) & is.na(partner_yr_sampled) ~ 1))

fCORT$partner_same <- fCORT$partner_samediff1 + fCORT$partner_samediff2 + fCORT$partner_samediff3

fCORT$partner_fCORT_yr_before <- fCORT$CORT_pg.mg[match(fCORT$partner_yr_before, fCORT$JID)]

#Add pair bond strength from video
pbvid1 <- pbvid1  %>%
  mutate(year =case_when(
    year == 14 ~ 2014,
    year == 15 ~ 2015,
    year == 18 ~ 2018,
    year == 19 ~ 2019))

pbvid1 <- pbvid1  %>% pivot_longer(
  cols = c("male.ID", "fem.ID"))
  
pbvid1$JID <- pbvid1$value
pbvid1$JID_year_grown <- paste(pbvid1$JID, pbvid1$year, sep = "_") 

pbvid1 <- pbvid1[,c("pair.ID", "JID_year_grown", "time_together_total_standardised")]
pbvid1$pbvid <- pbvid1$time_together_total_standardised
pbvid1$pair_ID <- pbvid1$pair.ID
pbvid1 <- pbvid1[,c("pair_ID", "JID_year_grown", "pbvid")]

### 

pbvid2$pair_ID <- paste(pbvid2$FID, pbvid2$MID, sep = "")

pbvid2 <- subset(pbvid2, pbvid2$Behavior == "Both")

pbvid2 <- pbvid2  %>% pivot_longer(
  cols = c("MID", "FID"))

pbvid2$JID <- pbvid2$value
pbvid2$JID_year_grown <- paste(pbvid2$JID, pbvid2$year, sep = "_") 

pbvid2 <- pbvid2[,c("pair_ID", "JID_year_grown", "percent_of_total_length")]
pbvid2$pbvid <- pbvid2$percent_of_total_length
pbvid2 <- pbvid2[,c("pair_ID", "JID_year_grown", "pbvid")]

breeding_summary3 <- subset(breeding_summary, year == 2023)
pbvid3$female_JID <- breeding_summary3$fem.ID[match(pbvid3$BOX, breeding_summary3$box)]
pbvid3$male_JID <- breeding_summary3$male.ID[match(pbvid3$BOX, breeding_summary3$box)]
pbvid3$pair_ID <- paste(pbvid3$female_JID, pbvid3$male_JID, sep = "")

pbvid3 <- pbvid3  %>% pivot_longer(
  cols = c("male_JID", "female_JID"))

pbvid3 <- dplyr::rename(pbvid3, JID = value)

pbvid3$year <- 2023

pbvid3$JID_year_grown <- paste(pbvid3$JID, pbvid3$year, sep = "_") 

pbvid3$pbvid <- pbvid3$BEHAVIOUR_DURATION_SUM_CONTROLLED

pbvid3 <- pbvid3[, c("pair_ID", "JID_year_grown", "pbvid")]

pbvid <- rbind(pbvid1, pbvid2, pbvid3)

#repeatability of pair bond strength (also add 2023?)
pbvid_rep <- rpt(pbvid ~ (1|pair_ID), grname = c("pair_ID"), data = pbvid, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(pbvid_rep)
plot(pbvid_rep)

fCORT$pbvid <- pbvid$pbvid[match(fCORT$JID_year_grown, pbvid$JID_year_grown)]
fCORT$pbvid <- fCORT$pbvid / 100
plot(fCORT$pbvid, fCORT$CORT_pg.mg)
length(table(fCORT$pbvid))

fCORT$pair_ID2 <- gsub(" ", "", fCORT$pair_ID)
fCORT$pair_ID_year2 <- gsub(" ", "", fCORT$pair_ID_year)

fCORT$pair_ID_yr_sampled2 <- gsub(" ", "", fCORT$pair_ID_yr_sampled)
fCORT$pair_ID_yr_sampled2 <- paste(fCORT$pair_ID_yr_sampled2, fCORT$feather_year_sampled, sep = "_")
  
fCORT$pair_ID_yr_before2 <- gsub(" ", "", fCORT$pair_ID_yr_before)
fCORT$pair_ID_yr_before2 <- paste(fCORT$pair_ID_yr_before2, fCORT$feather_year_before, sep = "_")

n <- 4
pbvid$year <- substr(pbvid$JID_year_grown, nchar(pbvid$JID_year_grown) - n + 1, nchar(pbvid$JID_year_grown))
pbvid$pairID_year <- paste(pbvid$pair_ID, pbvid$year, sep = "_")

#mean pair bond strength across all observations
#pbvid_mean <- as.data.frame(aggregate(x = pbvid$pbvid, by = list(pbvid$pair_ID), FUN = "mean"))
#pbvid_mean$pair_ID <- pbvid_mean$Group.1
#pbvid_mean$pbvid <- pbvid_mean$x
#pbvid_mean$Group.1 <- NULL
#pbvid_mean$x <- NULL

#mean pair-bond strength per pair per year
pbvid_mean <- as.data.frame(aggregate(x = pbvid$pbvid, by = list(pbvid$pairID_year), FUN = "mean"))
pbvid_mean$pairID_year <- pbvid_mean$Group.1
pbvid_mean$pbvid <- pbvid_mean$x
pbvid_mean$Group.1 <- NULL
pbvid_mean$x <- NULL

fCORT$pbvid_mean <- pbvid_mean$pbvid[match(fCORT$pair_ID_year2, pbvid_mean$pairID_year)]
fCORT$pbvid_mean <- pbvid_mean$pbvid[match(fCORT$pair_ID_yr_sampled2, pbvid_mean$pairID_year)]

fCORT$pbvid_mean_yr_before <- pbvid_mean$pbvid[match(fCORT$pair_ID_yr_before2, pbvid_mean$pairID_year)]
fCORT$pbvid_mean_yr_grown <- pbvid_mean$pbvid[match(fCORT$pair_ID_year2, pbvid_mean$pairID_year)]
fCORT$pbvid_mean_yr_sampled <- pbvid_mean$pbvid[match(fCORT$pair_ID_yr_sampled2, pbvid_mean$pairID_year)]

fCORT$pbvid_mean <- fCORT$pbvid_mean / 100
fCORT$pbvid_mean_yr_before <- fCORT$pbvid_mean_yr_before / 100
fCORT$pbvid_mean_yr_grown <- fCORT$pbvid_mean_yr_grown / 100
fCORT$pbvid_mean_yr_sampled <- fCORT$pbvid_mean_yr_sampled / 100

fCORT <- fCORT %>% 
  mutate(pbvid_mean =rowMeans(.[ , c("pbvid_mean_yr_before", "pbvid_mean_yr_grown", "pbvid_mean_yr_sampled")], na.rm=TRUE))

#Add early life conditions (not relevant because fCORT not repeatable)
fCORT$hatch_day_early <- LH_early$HATCH_DAY[match(fCORT$JID, LH_early$ID)]
fCORT$brood_size_early <- LH_early$BROOD_SIZE[match(fCORT$JID, LH_early$ID)]
fCORT$body_condition_early <- LH_early$BODY_CONDITION[match(fCORT$JID, LH_early$ID)]
fCORT$growth_rate_early <- LH_early$GROWTH_RATE[match(fCORT$JID, LH_early$ID)]

#Add reproductive success
fCORT$clutch_size_yr_grown <- breeding_summary2$cl.size[match(fCORT$JID_year_grown, breeding_summary2$JID_year_grown)]
fCORT$fledge_nr_yr_grown <- breeding_summary2$n.fledge[match(fCORT$JID_year_grown, breeding_summary2$JID_year_grown)]
fCORT$fledge_weight_yr_grown <- breeding_summary2$mass.ring[match(fCORT$JID_year_grown, breeding_summary2$JID_year_grown)]

fCORT$clutch_size_yr_sampled <- breeding_summary2$cl.size[match(fCORT$JID_year_sampled, breeding_summary2$JID_year_sampled)]
fCORT$fledge_nr_yr_sampled <- breeding_summary2$n.fledge[match(fCORT$JID_year_sampled, breeding_summary2$JID_year_sampled)]
fCORT$fledge_weight_yr_sampled <- breeding_summary2$mass.ring[match(fCORT$JID_year_sampled, breeding_summary2$JID_year_sampled)]

fCORT <- fCORT %>% 
  mutate(fledge_nr =rowMeans(.[ , c("fledge_nr_yr_grown", "fledge_nr_yr_sampled")], na.rm=TRUE))

fCORT <- fCORT %>% 
  mutate(fledge_weight =rowMeans(.[ , c("fledge_weight_yr_grown", "fledge_weight_yr_sampled")], na.rm=TRUE))

#Add box density in year feather grown
fCORT$box_density <- box_densities$Density25m[match(fCORT$box_grown_year, box_densities$box_grown_year)]
fCORT$box_density <- box_densities$Density25m[match(fCORT$box_year, box_densities$box_grown_year)]
fCORT$box_density <- as.numeric(fCORT$box_density)

fCORT$box_density_cat <- box_densities$DensityCat[match(fCORT$box_grown_year, box_densities$box_grown_year)]
fCORT$box_density_cat <- box_densities$DensityCat[match(fCORT$box_year, box_densities$box_grown_year)]

#Add disturbance levels
fCORT$disturbance <- box_densities$Disturbance[match(fCORT$box_grown_year, box_densities$box_grown_year)]
fCORT$disturbance <- box_densities$Disturbance[match(fCORT$box_year, box_densities$box_grown_year)]

#Add early life conditions 
fCORT$body_cond_early <- LH_early$BODY_CONDITION[match(fCORT$JID, LH_early$ID)]
fCORT$brood_size_early <- LH_early$BROOD_SIZE[match(fCORT$JID, LH_early$ID)]
fCORT$site_early <- LH_early$LOCATION[match(fCORT$JID, LH_early$ID)]
fCORT$site_early <- ifelse(fCORT$site_early == "STITHIANS", "Y", "Z")
fCORT$site_change <- ifelse(fCORT$site_early == fCORT$site, "same_site", "diff_site")

#Subset fCORT by sex
fCORT_F <- subset(fCORT, fCORT$sex == "F")
fCORT_M <- subset(fCORT, fCORT$sex == "M")

#Subset fCORT by site
fCORT_Z <- subset(fCORT, fCORT$site == "Z")
fCORT_Y <- subset(fCORT, fCORT$site == "Y")

fCORT_YZ <- subset(fCORT, fCORT$site == "Y" | fCORT$site == "Z")

fCORT_Z_sub <- subset(fCORT_Z, !is.na(fCORT_Z$box_grown))

#Specific box locations
distance_matrix <- read.csv("Data/distance_matrix.csv", header = T, stringsAsFactors = F)
distance_matrix$box_pair <- paste(distance_matrix$InputID, distance_matrix$TargetID, sep = " ")

fCORT_Z_boxes_mean_fCORT <- aggregate(x = fCORT_Z$CORT_pg.mg, by = list(fCORT_Z$box), FUN = "mean")
fCORT_Z_boxes_mean_fCORT$box <- fCORT_Z_boxes_mean_fCORT$Group.1
fCORT_Z_boxes_mean_fCORT$fCORT_mean <- fCORT_Z_boxes_mean_fCORT $x
fCORT_Z_boxes_mean_fCORT$Group.1 <- NULL
fCORT_Z_boxes_mean_fCORT$x <- NULL

fCORT_Z_boxes <- fCORT_Z[, c("box")]

fCORT_Z_boxes <- combn(fCORT_Z_boxes, 2, simplify = FALSE)
fCORT_Z_boxes <- do.call(rbind, fCORT_Z_boxes)

fCORT_Z_boxes <- as.data.frame(fCORT_Z_boxes)
fCORT_Z_boxes$box1 <- fCORT_Z_boxes$V1
fCORT_Z_boxes$box2 <- fCORT_Z_boxes$V2
fCORT_Z_boxes$V1 <- NULL
fCORT_Z_boxes$V2 <- NULL

fCORT_Z_boxes$box1_fCORT <- fCORT_Z_boxes_mean_fCORT$fCORT_mean[match(fCORT_Z_boxes$box1, fCORT_Z_boxes_mean_fCORT$box)]
fCORT_Z_boxes$box2_fCORT <- fCORT_Z_boxes_mean_fCORT$fCORT_mean[match(fCORT_Z_boxes$box2, fCORT_Z_boxes_mean_fCORT$box)]

fCORT_Z_boxes$box_pair <- paste(fCORT_Z_boxes$box1, fCORT_Z_boxes$box2, sep = " ")

fCORT_Z_boxes$box_distance <- distance_matrix$Distance[match(fCORT_Z_boxes$box_pair, distance_matrix$box_pair)]
fCORT_Z_boxes$fCORT_diff <- fCORT_Z_boxes$box1_fCORT - fCORT_Z_boxes$box2_fCORT
fCORT_Z_boxes$fCORT_diff <- sqrt(fCORT_Z_boxes$fCORT_diff^2)

plot(fCORT_Z_boxes$box_distance, fCORT_Z_boxes$fCORT_diff)
cor.test(fCORT_Z_boxes$box_distance, fCORT_Z_boxes$fCORT_diff)

#Subset fCORT_22
fCORT_22 <- subset(fCORT, fCORT$feather_year_sampled == "2023")

#subset fCORT each individual just once
fCORT_individuals <- subset(fCORT, !duplicated(fCORT$JID))
table(fCORT_individuals$sex)

#Add feeder data (visit rates, network positions)
all_individuals$fCORT <- feather_CORT_22$CORT_pg.mg[match(all_individuals$JID,feather_CORT_22$JID)]
all_individuals2 <- all_individuals[,c("JID", "visit_number", "visit_duration", "degree", "eigenvector", "betweenness", "strength", "mean_strength", "closeness", "transitivity")]

fCORT_22 <- merge(x = fCORT_22, y = all_individuals2, by = "JID", all.x = TRUE)

length(table(fCORT_22$visit_number))
fCORT_22 <- subset(fCORT_22, fCORT_22$visit_number > 0)

#Add box data (visit number prospecting and own box)
fCORT_22$visit_nr_own_box <- box_visits_22$visit_nr_own_box[match(fCORT_22$JID, box_visits_22$JID)]
fCORT_22$visit_nr_prospect <- box_visits_22$visit_nr_prospect[match(fCORT_22$JID, box_visits_22$JID)]

#Add feeder and box pair association indices
fCORT_22$pbfeedgmm <- pb22_long_pair$pbfeedgmm[match(fCORT_22$pair_ID, pb22_long_pair$pair_ID)]
fCORT_22$pbprospgmm <- pb22_long_pair$pbprospgmm[match(fCORT_22$pair_ID, pb22_long_pair$pair_ID)]
fCORT_22$pbowngmm <- pb22_long_pair$pbowngmm[match(fCORT_22$pair_ID, pb22_long_pair$pair_ID)]

length(table(fCORT_22$pbfeedgmm))

#Add feeder use/prospecting yes/no
fCORT_22$feeder <- ifelse(is.na(fCORT_22$visit_number), "no", "yes")
fCORT_22$prospect <- ifelse(is.na(fCORT_22$visit_nr_prospect), "no", "yes")

fCORT_22$pbfeeder <- ifelse(is.na(fCORT_22$pbfeedgmm), "no", "yes")
fCORT_22$pbprospect <- ifelse(is.na(fCORT_22$pbprospgmm), "no", "yes")

fCORT_22$pb_feeder_prospect <- ifelse(fCORT_22$pbfeeder == "yes" & fCORT_22$pbprospect =="yes", "yes", "no")

#Co-regulation: fCORT2 to compare F and M within pairs 

#fCORT_coreg <- fCORT[,c("CORT_pg.mg", "sex", "pair_ID", "JID", "feather_year_grown")]
#fCORT_coreg <- fCORT[,c("CORT_pg.mg", "sex", "pair_ID_yr_before", "pair_ID", "pair_ID_yr_sampled", "JID", "feather_year_grown", "feather_year_sampled", "JID_year_grown")]
fCORT_coreg <- fCORT[,c("CORT_pg.mg", "sex", "partner_yr_sampled", "pair_ID", "pair_ID_yr_sampled", "JID", "feather_year_grown", "feather_year_sampled", "JID_year_grown")]

fCORT_mean <- aggregate(x = fCORT$CORT_pg.mg, by = list(fCORT$JID_year_grown), FUN = "mean")
fCORT_mean$JID_year_grown <- fCORT_mean$Group.1
fCORT_mean$fCORT_mean <- fCORT_mean$x
fCORT_mean$Group.1 <- NULL
fCORT_mean$x <- NULL

fCORT_coreg$fCORT_mean <- fCORT_mean$fCORT_mean[match(fCORT_coreg$JID_year_grown, fCORT_mean$JID_year_grown)]

#fCORT_coreg <- fCORT_coreg[-79,] #remove case where 2 feathers per individual in one year
#fCORT_coreg <- subset(fCORT_coreg, !is.na(fCORT_coreg$pair_ID))
#fCORT_coreg <- subset(fCORT_coreg, !is.na(fCORT_coreg$pair_ID_yr_sampled))
fCORT_coreg <- subset(fCORT_coreg, !is.na(fCORT_coreg$partner_yr_sampled))

#fCORT_coreg$JID_year <- paste(fCORT_coreg$JID, fCORT_coreg$feather_year_grown, sep = "_")
#fCORT_coreg$pair_ID_year <- paste(fCORT_coreg$pair_ID, fCORT_coreg$feather_year_grown, sep = "_")
fCORT_coreg$pair_ID_year <- paste(fCORT_coreg$pair_ID_yr_sampled, (fCORT_coreg$feather_year_sampled), sep = "_")

fCORT_coreg <- fCORT_coreg %>% group_by(fCORT_coreg$pair_ID) %>% filter( n() > 1 )
fCORT_coreg <- fCORT_coreg %>% group_by(fCORT_coreg$pair_ID_year) %>% filter( n() > 1 )

fCORT_coreg$pair_ID_year_sex <- paste(fCORT_coreg$pair_ID_year, fCORT_coreg$sex, sep = "_")

fCORT_coreg <- fCORT_coreg %>% filter(!duplicated(pair_ID_year_sex))

fCORT_coreg <- fCORT_coreg %>% pivot_wider(
  id_cols = c("pair_ID_year"), names_from = c("sex"), values_from = c("fCORT_mean"))

#fCORT_coreg <- fCORT_coreg %>% pivot_wider(
#  id_cols = c("pair_ID_year", "pair_ID"), names_from = c("sex"), values_from = c("CORT_pg.mg"))

fCORT_coreg <- na.omit(fCORT_coreg)

plot(fCORT_coreg$M, fCORT_coreg$F)
     
#run next 7 lines only when co-reg mean across years
fCORT_coregF <- aggregate(x = fCORT_F$CORT_pg.mg, by = list(fCORT_F$pair_ID), FUN = "mean")
fCORT_coregM <- aggregate(x = fCORT_M$CORT_pg.mg, by = list(fCORT_M$pair_ID), FUN = "mean")

fCORT_coregF <- aggregate(x = fCORT_F$CORT_pg.mg, by = list(fCORT_F$pair_ID_yr_sampled), FUN = "mean")
fCORT_coregM <- aggregate(x = fCORT_M$CORT_pg.mg, by = list(fCORT_M$pair_ID_yr_sampled), FUN = "mean")

fCORT_coreg <- merge(fCORT_coregF, fCORT_coregM, by = "Group.1", all = TRUE)
fCORT_coreg$CORT_F <- fCORT_coreg$x.x
fCORT_coreg$CORT_M <- fCORT_coreg$x.y
fCORT_coreg$pair_ID <- fCORT_coreg$Group.1
fCORT_coreg <- fCORT_coreg[, c("CORT_F","CORT_M", "pair_ID")]

fCORT_coreg <- na.omit(fCORT_coreg) 

fCORT_coreg$F <- as.numeric(fCORT_coreg$F)
fCORT_coreg$M <- as.numeric(fCORT_coreg$M)
fCORT_coreg$CORT_F <- fCORT_coreg$F
fCORT_coreg$CORT_M <- fCORT_coreg$M
fCORT_coreg$F <- NULL
fCORT_coreg$M <- NULL

write.csv(fCORT_coreg,"fCORT_coreg.csv", row.names = FALSE)

#manually added pairs for which data was not available in the same year
fCORT_coreg <- read.csv("fCORT_coreg_copy.csv")
fCORT_coreg <- subset(fCORT_coreg, !is.na(fCORT_coreg$pair_ID_year))

fCORT_coreg$pair_ID <- str_split(fCORT_coreg$pair_ID_year, "_", simplify = TRUE)
fCORT_coreg$pair_ID <- fCORT_coreg$pair_ID[,1]

fCORT_coreg$female_ID <- str_split(fCORT_coreg$pair_ID, " ", simplify = TRUE)
fCORT_coreg$female_ID <- fCORT_coreg$female_ID[,1]

fCORT_coreg$male_ID <- str_split(fCORT_coreg$pair_ID, " ", simplify = TRUE)
fCORT_coreg$male_ID <- fCORT_coreg$male_ID[,2]

fCORT_coreg$year_paired <- LH_paired_repaired$year[match(fCORT_coreg$pair_ID, LH_paired_repaired$pair_ID)]
fCORT_coreg$year_paired <- as.numeric(fCORT_coreg$year_paired)
fCORT_coreg$years_together <- fCORT_coreg$year - fCORT_coreg$year_paired
fCORT_coreg$pbvid_mean <- fCORT_pbvid$pbvid_mean[match(fCORT_coreg$pair_ID, fCORT_pbvid$pair_ID)]
fCORT_coreg$site <- fCORT$site[match(fCORT_coreg$pair_ID, fCORT$pair_ID_yr_sampled)] 

fCORT_coreg$CORT_diff <- fCORT_coreg$female_CORT - fCORT_coreg$male_CORT

fCORT_coreg$male_body_cond <- fCORT$body_cond_noNA[match(fCORT_coreg$male_ID, fCORT$JID)]
fCORT_coreg$female_body_cond <- fCORT$body_cond_noNA[match(fCORT_coreg$female_ID, fCORT$JID)]

fCORT_coreg_long <- fCORT_coreg  %>% pivot_longer(
  cols = c("female_CORT", "male_CORT"), names_to = "CORT_sex", values_to = "CORT")

fCORT_coreg_long$sex <- ifelse(fCORT_coreg_long$CORT_sex == "female_CORT", "F", "M")

fCORT_coreg_long$sex_numeric <- ifelse(fCORT_coreg_long$sex == "F", 0, 1)

fCORT_coreg_long$JID <- ifelse(fCORT_coreg_long$sex == "F", fCORT_coreg_long$female_ID, fCORT_coreg_long$male_ID)

fCORT_coreg_long$year_sampled <- ifelse(fCORT_coreg_long$sex == "F", fCORT_coreg_long$female_CORT_year, fCORT_coreg_long$male_CORT_year)

fCORT_coreg_long$JID_year_sampled <- paste(fCORT_coreg_long$JID, fCORT_coreg_long$year_sampled, sep = "_")

fCORT_coreg_long$partner_JID <- ifelse(fCORT_coreg_long$sex == "M", fCORT_coreg_long$female_ID, fCORT_coreg_long$male_ID)

fCORT_coreg_long$body_cond <- fCORT$body_cond_noNA[match(fCORT_coreg_long$JID_year_sampled, fCORT$JID_year_sampled)]
fCORT_coreg_long$site <- fCORT$site[match(fCORT_coreg_long$JID_year_sampled, fCORT$JID_year_sampled)]
fCORT_coreg_long$fledge_weight <- fCORT$fledge_weight[match(fCORT_coreg_long$JID_year_sampled, fCORT$JID_year_sampled)]

fCORT_coreg_long_F <- subset(fCORT_coreg_long, fCORT_coreg_long$sex == "F")
fCORT_coreg_long_M <- subset(fCORT_coreg_long, fCORT_coreg_long$sex == "M")

fCORT_coreg_long_F$partner_CORT <- fCORT_coreg_long_M$CORT[match(fCORT_coreg_long_F$pair_ID_year, fCORT_coreg_long_M$pair_ID_year)]
fCORT_coreg_long_M$partner_CORT <- fCORT_coreg_long_F$CORT[match(fCORT_coreg_long_M$pair_ID_year, fCORT_coreg_long_F$pair_ID_year)]

fCORT_coreg_long <- rbind(fCORT_coreg_long_F, fCORT_coreg_long_M)

fCORT_coreg_long <- subset(fCORT_coreg_long, !fCORT_coreg_long$site == "X")

fCORT_coreg_long2 <- subset(fCORT_coreg_long, fCORT_coreg_long$CORT < 7.817757)

fCORT_coreg_long_Y <- subset(fCORT_coreg_long2, fCORT_coreg_long2$site == "Y")
fCORT_coreg_long_Z <- subset(fCORT_coreg_long2, fCORT_coreg_long2$site == "Z")

fCORT_coreg_long$box_grown <- fCORT$box_grown[match(fCORT_coreg_long$JID, fCORT$JID)]
fCORT_coreg_long$box_sampled <- fCORT$box[match(fCORT_coreg_long$JID, fCORT$JID)]

fCORT_coreg2 <- subset(fCORT_coreg, fCORT_coreg$CORT_F < 7)
fCORT_coreg2 <- subset(fCORT_coreg, fCORT_coreg$CORT_F < 6.427079)
fCORT_coreg2 <- subset(fCORT_coreg, fCORT_coreg$CORT_M < 6.427079)

fCORT_individuals <- distinct(fCORT, JID, .keep_all = TRUE)
fCORT_individuals$sample_no <- table(fCORT$JID)
length(table(fCORT$JID))
fCORT_individuals2 <- subset( fCORT_individuals, fCORT_individuals$sample_no > 1)

write.csv(fCORT_coreg_long,"fCORT_coreg_long.csv", row.names = FALSE)

dyad_edge$ind1_CORT <- fCORT_22_sna$CORT_pg.mg[match(dyad_edge$ind1, fCORT_22_sna$JID)]
dyad_edge$ind2_CORT <- fCORT_22_sna$CORT_pg.mg[match(dyad_edge$ind2, fCORT_22_sna$JID)]

dyad_edge$ind1ind2_CORT_diff <- dyad_edge$ind1_CORT - dyad_edge$ind2_CORT
dyad_edge$ind1ind2_CORT_diff_norm <- sqrt(dyad_edge$ind1ind2_CORT_diff^2)

plot(dyad_edge$g_edge_weights, dyad_edge$ind1ind2_CORT_diff_norm)  
cor.test(dyad_edge$g_edge_weights, dyad_edge$ind1ind2_CORT_diff_norm)  

dyad_edge_CORT <- subset(dyad_edge, !is.na(dyad_edge$ind1_CORT))
dyad_edge_CORT <- subset(dyad_edge_CORT, !is.na(dyad_edge_CORT$ind2_CORT))

dyad_edge_CORT_long <- dyad_edge_CORT  %>% pivot_longer(
  cols = c("ind1_CORT", "ind2_CORT"), names_to = "CORT_ind", values_to = "CORT")

dyad_edge_CORT_long$ind <- ifelse(dyad_edge_CORT_long$CORT_ind == "ind1_CORT", "ind1", "ind2")

dyad_edge_CORT_long$JID <- ifelse(dyad_edge_CORT_long$ind == "ind1", dyad_edge_CORT_long$ind1, dyad_edge_CORT_long$ind2)

dyad_edge_CORT_long$associate_JID <- ifelse(dyad_edge_CORT_long$ind == "ind2", dyad_edge_CORT_long$ind1, dyad_edge_CORT_long$ind2)

dyad_edge_CORT_long$JID_year_sampled <- paste(dyad_edge_CORT_long$JID, "2023", sep = "_")

dyad_edge_CORT_long$body_cond <- fCORT$body_cond_noNA[match(dyad_edge_CORT_long$JID_year_sampled, fCORT$JID_year_sampled)]
dyad_edge_CORT_long$site <- fCORT$site[match(dyad_edge_CORT_long$JID_year_sampled, fCORT$JID_year_sampled)]

dyad_edge_CORT_long_ind1 <- subset(dyad_edge_CORT_long, dyad_edge_CORT_long$ind == "ind1")
dyad_edge_CORT_long_ind2 <- subset(dyad_edge_CORT_long, dyad_edge_CORT_long$ind == "ind2")

dyad_edge_CORT_long_ind1$associate_CORT <- dyad_edge_CORT_long_ind2$CORT[match(dyad_edge_CORT_long_ind1$g_dyads, dyad_edge_CORT_long_ind2$g_dyads)]
dyad_edge_CORT_long_ind2$associate_CORT <- dyad_edge_CORT_long_ind1$CORT[match(dyad_edge_CORT_long_ind2$g_dyads, dyad_edge_CORT_long_ind1$g_dyads)]

dyad_edge_CORT_long <- rbind(dyad_edge_CORT_long_ind1, dyad_edge_CORT_long_ind2)

dyad_edge_CORT_long2 <- subset(dyad_edge_CORT_long, dyad_edge_CORT_long$CORT < 7.817757)

dyad_edge_CORT_long$associate_CORT_z <- scale(dyad_edge_CORT_long$associate_CORT)
dyad_edge_CORT_long$g_edge_weights_z <- scale(dyad_edge_CORT_long$g_edge_weights)
dyad_edge_CORT_long$body_cond_z <- scale(dyad_edge_CORT_long$body_cond)

write.csv(dyad_edge_CORT_long,"dyad_edge_CORT_long.csv", row.names = FALSE)

#(02) SUMMARY DESCRIPTIVE STATISTICS ----
#Summary descriptive statistics
length(table(fCORT$JID)) #166 feathers from 124 individuals 
length(table(fCORT$pair_ID)) 
table(fCORT$sex) 
table(fCORT_individuals$sex)
length(table(fCORT$year)) #9 years
table(fCORT$feather_year_grown)
length(fCORT$JID[fCORT$feather_year_grown == 2014])

fCORT_mean <- mean(fCORT$CORT_pg.mg)
fCORT_sd <- sd(fCORT$CORT_pg.mg)
fCORT_mean
fCORT_sd

fCORT_mean + 3 * fCORT_sd

mean(fCORT$extraction_weight_g)
sd(fCORT$extraction_weight_g)

mean(fCORT$feather_length)
sd(fCORT$feather_length)

mean(fCORT$feather_weight)
sd(fCORT$feather_weight)

plot(fCORT$feather_length, fCORT$CORT_pg.mg)
cor.test(fCORT$feather_length, fCORT$CORT_pg.mg)

plot(fCORT$feather_weight, fCORT$CORT_pg.mg)
cor.test(fCORT$feather_weight, fCORT$CORT_pg.mg)

plot(fCORT$weight, fCORT$CORT_pg.mg)
cor.test(fCORT$weight, fCORT$CORT_pg.mg)


#(03) STATISTICAL ANALYSES ----

write.csv(fCORT,"Data for manuscript/fCORT.csv", row.names = FALSE)
write.csv(LH,"Data for manuscript/LH.csv", row.names = FALSE)
write.csv(fCORT_coreg,"Data for manuscript/fCORT_coreg.csv", row.names = FALSE)
write.csv(fCORT_22,"Data for manuscript/fCORT_22.csv", row.names = FALSE)
write.csv(all_individuals_juv,"Data for manuscript/all_individuals_juv.csv", row.names = FALSE)
write.csv(visit_data,"Data for manuscript/visit_data.csv", row.names = FALSE)
write.csv(am,"Data for manuscript/adjacency_matrix.csv", row.names = FALSE)
write.csv(dyad_edge_CORT_long,"dyad_edge_CORT_long.csv", row.names = FALSE)
write.csv(fCORT_coreg_long,"fCORT_coreg_long.csv", row.names = FALSE)

fCORT <- read.csv("Data for manuscript/fCORT.csv")

#Subset fCORT without outliers (3 X sd of fCORT mean)
fCORT_nool <- subset(fCORT, fCORT$CORT_pg.mg < 7.817757)

#1. Repeatability of fCORT ----
fCORT_rep <- rpt(log(CORT_pg.mg) ~ (1|JID), grname = c("JID"), data = fCORT, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(fCORT_rep)

fCORT_rep <- rpt(log(CORT_pg.mg) ~ (1|JID), grname = c("JID"), data = fCORT_nool, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(fCORT_rep)

windowsFonts(A = windowsFont("Garamond"))  

plot(fCORT_rep, main = NA, cex = 2, cex.axis = 1.4, cex.lab = 1.4, family = "A")

fCORT_rep2 <- rpt(CORT_pg.mg ~ (1|pair_ID_yr_sampled), grname = c("pair_ID_yr_sampled"), data = fCORT, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(fCORT_rep2)
plot(fCORT_rep2)

#Repeatability of morphometrics
tarsus_rep <- rpt(TARSUS ~ (1|ID), grname = c("ID"), data = LH_adults, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(tarsus_rep)

weight_rep <- rpt(WEIGHT ~ (1|ID), grname = c("ID"), data = LH_adults, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(weight_rep)

bodycond_rep <- rpt(body_cond ~ (1|ID), grname = c("ID"), data = LH_adults, datatype = "Gaussian", nboot = 1000, npermut = 1000)
summary(bodycond_rep)

#2. Full model of environmental and individual covariates ----
fCORT$feather_year_grown <- as.factor(fCORT$feather_year_grown)
fCORT$site <- as.factor(fCORT$site)

fCORT$body_cond_c <- fCORT$body_cond_noNA - mean(fCORT$body_cond_noNA, na.rm = TRUE)
fCORT$age_c <- fCORT$age - mean(fCORT$age, na.rm = TRUE)

fCORT$body_cond_z <- scale(fCORT$body_cond_noNA)
fCORT$age_z <- scale(fCORT$age)

fCORT$pair_ID_yr_sampled[fCORT$pair_ID_yr_sampled == "J4646 J4645"] <- "J4645 J4646"

fCORT_nool <- subset(fCORT, fCORT$CORT_pg.mg < 7.817757)

#Brms model 

#Model just with prior
default_prior()

fCORT_fullcov_brm1_prior <- brm(CORT_pg.mg ~ site + body_cond_z + sex + 
                         age_z + (1|feather_year_grown) + (1|JID) + 
                         (1|pair_ID_yr_sampled),
                         data = fCORT, 
                         family = lognormal(link = "identity"),
                           prior = c(
                             prior(normal(1.25, 0.3), class = "Intercept"),  
                             prior(normal(0, 0.3), class = "b"),
                             prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                         sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_fullcov_brm1_prior, ndraws = 100) +
  scale_x_log10()

#Model including body condition
fCORT_fullcov_brm1 <- brm(CORT_pg.mg ~ site + sex + body_cond_z +
                                age_z + (1|feather_year_grown) + (1|JID) + 
                                (1|pair_ID_yr_sampled),
                                data = fCORT, 
                                family = lognormal(link = "identity"),
                                prior = c(
                                  prior(normal(1.25, 0.3), class = "Intercept"),  
                                  prior(normal(0, 0.3), class = "b"),
                                  prior(lognormal(log(0.35), 0.3), class = "sigma")),            
)

fCORT_fullcov_nool_brm1 <- brm(CORT_pg.mg ~ site + sex + body_cond_z +
                            age_z + (1|feather_year_grown) + (1|JID) + 
                            (1|pair_ID_yr_sampled),
                          data = fCORT_nool, 
                          family = lognormal(link = "identity"),
                          prior = c(
                            prior(normal(1.25, 0.3), class = "Intercept"),  
                            prior(normal(0, 0.3), class = "b"),
                            prior(lognormal(log(0.35), 0.3), class = "sigma")),            
)

#Model excluding body condition 
fCORT_fullcov_brm1 <- brm(CORT_pg.mg ~ site + sex +
                            age_z + (1|feather_year_grown) + (1|JID) + 
                            (1|pair_ID_yr_sampled),
                          data = fCORT, 
                          family = lognormal(link = "identity"),
                          prior = c(
                            prior(normal(1.25, 0.3), class = "Intercept"),  
                            prior(normal(0, 0.3), class = "b"),
                            prior(lognormal(log(0.35), 0.3), class = "sigma")),            
)

fCORT_fullcov_nool_brm1 <- brm(CORT_pg.mg ~ site + sex +
                                 age_z + (1|feather_year_grown) + (1|JID) + 
                                 (1|pair_ID_yr_sampled),
                               data = fCORT_nool, 
                               family = lognormal(link = "identity"),
                               prior = c(
                                 prior(normal(1.25, 0.3), class = "Intercept"),  
                                 prior(normal(0, 0.3), class = "b"),
                                 prior(lognormal(log(0.35), 0.3), class = "sigma")),            
)

summary(fCORT_fullcov_brm1)
summary(fCORT_fullcov_nool_brm1)

check_collinearity(fCORT_fullcov_brm1)
check_collinearity(fCORT_fullcov_nool_brm1)

#Posterior predictive checks
pp_check(fCORT_fullcov_brm1, type = "dens_overlay", ndraws = 1000)
pp_check(fCORT_fullcov_nool_brm1, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_fullcov_brm1)

launch_shinystan(fCORT_fullcov_brm1)
plot(fCORT_fullcov_nool_brm1)
launch_shinystan(fCORT_fullcov_nool_brm1)

#Posterior distribution
as_draws_df(fCORT_fullcov_brm1)
mcmc_areas(fCORT_fullcov_brm1)
mcmc_intervals(fCORT_fullcov_brm1)

as_draws_df(fCORT_fullcov_nool_brm1)
mcmc_areas(fCORT_fullcov_nool_brm1)
mcmc_intervals(fCORT_fullcov_nool_brm1)

#Evaluation and interpretation
loo(fCORT_fullcov_brm1, moment_match = TRUE)
fitted(fCORT_fullcov_brm1, scale = "response")
conditional_effects(fCORT_fullcov_brm1)

loo(fCORT_fullcov_nool_brm1, moment_match = TRUE)
fitted(fCORT_fullcov_nool_brm1, scale = "response")
conditional_effects(fCORT_fullcov_nool_brm1)

bayes_R2(fCORT_fullcov_brm1)
emmeans(fCORT_fullcov_brm1, ~ site, type = "response") |> pairs()

bayes_R2(fCORT_fullcov_nool_brm1)
emmeans(fCORT_fullcov_nool_brm1, ~ site, type = "response") |> pairs()

#Prior checks 
prior_summary(fCORT_fullcov_nool_brm1)
prior_summary(fCORT_fullcov_nool_brm1)

#Extract predictions
fitted(fCORT_fullcov_brm1)
fitted(fCORT_fullcov_nool_brm1)

hypothesis(fCORT_fullcov_brm1, "body_cond_z < 0") 
hypothesis(fCORT_fullcov_brm1, "age_z < 0") 
hypothesis(fCORT_fullcov_brm1, "sexM < 0") 

hypothesis(fCORT_fullcov_nool_brm1, "body_cond_z < 0") 
hypothesis(fCORT_fullcov_nool_brm1, "age_z < 0") 
hypothesis(fCORT_fullcov_nool_brm1, "sexM < 0") 

summary(fCORT_fullcov_brm1)

posterior_marginal <- fCORT %>%
  add_epred_draws(
    object = fCORT_fullcov_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_fCORT <- posterior_marginal %>%
  dplyr::group_by(.draw, site) %>%
  dplyr::summarise(
    epred = mean(.epred),
    .groups = "drop"
  )
posterior_draw_fCORT

posterior_summary <- posterior_draw_fCORT %>%
  group_by(site) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_fCORT %>%
  compare_levels(epred, by = site)

pairwise_contrasts <- pairwise_contrasts %>%
  dplyr::rename(
    contrast = site,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Body condition plot for manuscript
cond <- conditional_effects(fCORT_fullcov_brm1, effect = "body_cond_z")
p <- plot(cond, points=F)
body_plot <- p[[1]] + 
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  scale_y_continuous(limits = c(0, 9)) +
  geom_point(
    aes(x = body_cond_z, y = CORT_pg.mg), 
    data = fCORT, 
    color = "black",
    size = 2,
    alpha = 0.3,
    shape = 16,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Body condition (residuals)", y = "fCORT (pg/mg)")
body_plot


#Sex plot for manuscript
colours <- c("#FFFFFF", "#CCCCCC")

sex_plot <- ggplot(fCORT, aes(x = sex, y = CORT_pg.mg, fill = sex)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Sex", y = "fCORT (pg/mg)") +
  theme_few(base_size = 14)+
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 2) +
  stat_summary(fun=mean, geom="point", shape=8, size=2, col = "coral", stroke = 1.25) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
sex_plot

#Age plot for manuscript
age_plot <- ggplot(fCORT, aes(x= age, y= CORT_pg.mg)) +
  geom_point(alpha = 0.3, size=2, shape=16) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Age (years)", y = "fCORT (pg/mg)", size =14) +
  theme_few(base_size = 14) + 
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
age_plot

#Site plot for manuscript
colours <- c("grey100", "grey90", "grey80", "grey70")

site_plot <- ggplot(fCORT, aes(x = site, y = CORT_pg.mg, fill = site)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Study site", y = "fCORT (pg/mg)") +
  theme_few(base_size = 14)+
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 2) +
  stat_summary(fun=mean, geom="point", shape=8, size=2, col = "coral", stroke = 1.25) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
site_plot

env_indiv_fCORT_plot <- ggarrange(body_plot, site_plot, age_plot, sex_plot, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"),  widths = c(1, 1), font.label = list(
  family = "Garamond",
  face = "plain",   
  size = 14,
  color = "black"
))

env_indiv_fCORT_plot


#3. Pair bond models ----

fCORT$pbvid_mean_c <- fCORT$pbvid_mean - mean(fCORT$pbvid_mean, na.rm = TRUE)
fCORT$pbvid_mean_z <- scale(fCORT$pbvid_mean)
fCORT$feather_year_grown_z <- scale(fCORT$feather_year_grown)

fCORT_pbvid <- subset(fCORT, !is.na(fCORT$pbvid_mean))
fCORT_pbvid <- subset(fCORT, !is.nan(fCORT$pbvid_mean))

fCORT_pbvid$pair_ID_yr_sampled[fCORT_pbvid$pair_ID_yr_sampled == "J4646 J4645"] <- "J4645 J4646"

table(fCORT_pbvid$site)

fCORT_pbvid_noX <- subset(fCORT_pbvid, !site =="X")
fCORT_pbvid_noX$site <- droplevels(fCORT_pbvid_noX$site)

fCORT_pbvid_nool <- subset(fCORT_pbvid, fCORT_pbvid$CORT_pg.mg < 7.817757) #remove values for fCORT > 3 SD

fCORT_pbvid_nool <- subset(fCORT_pbvid, fCORT_pbvid$pbvid_mean < 0.75) #remove values for pbvid > 3 SD
fCORT_pbvid_nool <- subset(fCORT_pbvid_nool, fCORT_pbvid_nool$pbvid_mean < 0.75) #remove values for pbvid > 3 SD

fCORT_pbvid$pbvid_mean_binary <- ifelse(fCORT_pbvid$pbvid_mean > median(fCORT_pbvid$pbvid_mean), 1, 0)

#Brms model 

#Pair-bond model 1: pair-bond strength + body condition + site

#Model just with prior
default_prior()

fCORT_pbvid_brm1_prior <- brm(CORT_pg.mg ~ pbvid_mean_z + site + body_cond_z + 
                                  (1|feather_year_grown) + (1|JID) + 
                                  (1|pair_ID_yr_sampled),
                                data = fCORT_pbvid, 
                                family = lognormal(link = "identity"),
                                prior = c(
                                  prior(normal(1.25, 0.3), class = "Intercept"),  
                                  prior(normal(0, 0.3), class = "b"),
                                  prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                                sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_pbvid_brm1_prior, ndraws = 100) +
  scale_x_log10()

#model including extreme values
fCORT_pbvid_brm1 <- brm(CORT_pg.mg ~ pbvid_mean_z + site + body_cond_z + 
                                (1|feather_year_grown) + (1|JID) + 
                                (1|pair_ID_yr_sampled),
                              data = fCORT_pbvid, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(1.25, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),
                              save_pars = save_pars(all = TRUE)
                              )

#model excluding extreme values
fCORT_pbvid_nool_brm1 <- brm(CORT_pg.mg ~ pbvid_mean_z + site + body_cond_z + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid_nool, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma")),
                        save_pars = save_pars(all = TRUE)
)

#Model summary
summary(fCORT_pbvid_brm1)
summary(fCORT_pbvid_nool_brm1)

#Check collinearity
check_collinearity(fCORT_pbvid_brm1)
check_collinearity(fCORT_pbvid_nool_brm1)

#Prior checks 
prior_summary(fCORT_pbvid_brm1)
prior_summary(fCORT_pbvid_nool_brm1)

#Posterior predictive checks
pp_check(fCORT_pbvid_brm1, type = "dens_overlay", ndraws = 1000)
pp_check(fCORT_pbvid_nool_brm1, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_pbvid_brm1)
launch_shinystan(fCORT_pbvid_brm1)
plot(fCORT_pbvid_nool_brm1)
launch_shinystan(fCORT_pbvid_nool_brm1)

#Posterior distribution
as_draws_df(fCORT_pbvid_brm1)
mcmc_areas(fCORT_pbvid_brm1)
mcmc_intervals(fCORT_pbvid_brm1)
as_draws_df(fCORT_pbvid_nool_brm1)
mcmc_areas(fCORT_pbvid_nool_brm1)
mcmc_intervals(fCORT_pbvid_nool_brm1)

#Evaluation and interpretation
fCORT_pbvid_brm1_loo <- loo(fCORT_pbvid_brm1, moment_match = TRUE)
fCORT_pbvid_brm1_loo 
fCORT_pbvid_nool_brm1_loo <- loo(fCORT_pbvid_nool_brm1, moment_match = TRUE)
fCORT_pbvid_nool_brm1_loo 

fitted(fCORT_pbvid_brm1, scale = "response")
conditional_effects(fCORT_pbvid_brm1)
fitted(fCORT_pbvid_nool_brm1, scale = "response")
conditional_effects(fCORT_pbvid_nool_brm1)

bayes_R2(fCORT_pbvid_brm1)
emmeans(fCORT_pbvid_brm1, ~ site, type = "response") |> pairs()
bayes_R2(fCORT_pbvid_nool_brm1)
emmeans(fCORT_pbvid_nool_brm1, ~ site, type = "response") |> pairs()

hypothesis(fCORT_pbvid_brm1, "pbvid_mean_z < 0") 
hypothesis(fCORT_pbvid_nool_brm1, "pbvid_mean_z < 0") 

posterior_marginal <- fCORT_pbvid %>%
  add_epred_draws(
    object = fCORT_pbvid_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_fCORT <- posterior_marginal %>%
  group_by(.draw, site) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )
posterior_draw_fCORT

posterior_summary <- posterior_draw_fCORT %>%
  group_by(site) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_relationship %>%
  compare_levels(epred, by = relationship)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = relationship,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

fCORT_pbvid_glmm1 <- glmmTMB(log(CORT_pg.mg) ~ pbvid_mean_z + site + body_cond_z + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid_nool, 
                        family = gaussian(link = "identity"))

Anova(fCORT_pbvid_glmm1)

#Pair-bond model 2: pair-bond strength * body condition + site

#Model just with prior
default_prior()

fCORT_pbvid_brm2_prior <- brm(CORT_pg.mg ~ pbvid_mean_z * body_cond_z + site + 
                                (1|feather_year_grown) + (1|JID) + 
                                (1|pair_ID_yr_sampled),
                              data = fCORT_pbvid, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(1.25, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_pbvid_brm2_prior, ndraws = 100) +
  scale_x_log10()

#model including extreme values
fCORT_pbvid_brm2 <- brm(CORT_pg.mg ~ pbvid_mean_z * body_cond_z + site + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma")),
                        save_pars = save_pars(all = TRUE)
)

#model excluding extreme values
fCORT_pbvid_nool_brm2 <- brm(CORT_pg.mg ~ pbvid_mean_z * body_cond_z + site + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid_nool, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma")),
                        save_pars = save_pars(all = TRUE)
)

#Model summary
summary(fCORT_pbvid_brm2)
summary(fCORT_pbvid_nool_brm2)

#Check collinearity
check_collinearity(fCORT_pbvid_brm2)
check_collinearity(fCORT_pbvid_nool_brm2)

#Prior checks 
prior_summary(fCORT_pbvid_brm2)

#Posterior predictive checks
pp_check(fCORT_pbvid_brm2, type = "dens_overlay", ndraws = 1000)
pp_check(fCORT_pbvid_nool_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_pbvid_brm2)
launch_shinystan(fCORT_pbvid_brm2)
plot(fCORT_pbvid_nool_brm2)
launch_shinystan(fCORT_pbvid_nool_brm2)

#Posterior distribution
as_draws_df(fCORT_pbvid_brm2)
mcmc_areas(fCORT_pbvid_brm2)
mcmc_intervals(fCORT_pbvid_brm2)

as_draws_df(fCORT_pbvid_nool_brm2)
mcmc_areas(fCORT_pbvid_nool_brm2)
mcmc_intervals(fCORT_pbvid_nool_brm2)

#Evaluation and interpretation
fCORT_pbvid_brm2_loo <- loo(fCORT_pbvid_brm2, moment_match = TRUE)
loo_compare(fCORT_pbvid_brm1_loo, fCORT_pbvid_brm2_loo)
fitted(fCORT_pbvid_brm2, scale = "response")
conditional_effects(fCORT_pbvid_nool_brm2)
fCORT_pbvid_nool_brm2_loo <- loo(fCORT_pbvid_nool_brm2, moment_match = TRUE)
loo_compare(fCORT_pbvid_nool_brm1_loo, fCORT_pbvid_nool_brm2_loo)
fitted(fCORT_pbvid_nool_brm2, scale = "response")
conditional_effects(fCORT_pbvid_nool_brm2)

bayes_R2(fCORT_pbvid_brm2)
emmeans(fCORT_pbvid_brm2, ~ site, type = "response") |> pairs()
bayes_R2(fCORT_pbvid_nool_brm2)
emmeans(fCORT_pbvid_nool_brm2, ~ site, type = "response") |> pairs()

hypothesis(fCORT_pbvid_brm2, "pbvid_mean_z < 0") 
hypothesis(fCORT_pbvid_brm2, "pbvid_mean_z:body_cond_z < 0") 
hypothesis(fCORT_pbvid_nool_brm2, "pbvid_mean_z < 0") 
hypothesis(fCORT_pbvid_nool_brm2, "pbvid_mean_z:body_cond_z < 0") 


#Pair-bond model 3: pair-bond strength * site + body condition

#Model just with prior
default_prior()

fCORT_pbvid_brm3_prior <- brm(CORT_pg.mg ~ pbvid_mean_z * site + body_cond_z + 
                                (1|feather_year_grown) + (1|JID) + 
                                (1|pair_ID_yr_sampled),
                              data = fCORT_pbvid, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(1.25, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_pbvid_brm3_prior, ndraws = 100) +
  scale_x_log10()

#model including extreme values
fCORT_pbvid_brm3 <- brm(CORT_pg.mg ~ pbvid_mean_z * site + body_cond_z + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma")),
                        save_pars = save_pars(all = TRUE)
)

#model excluding extreme values
fCORT_pbvid_nool_brm3 <- brm(CORT_pg.mg ~ pbvid_mean_z * site + body_cond_z + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid_nool, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma"))            
)

#Model summary
summary(fCORT_pbvid_brm3)
summary(fCORT_pbvid_nool_brm3)

#Check collinearity
check_collinearity(fCORT_pbvid_brm3)
check_collinearity(fCORT_pbvid_nool_brm3)

#Prior checks
prior_summary(fCORT_pbvid_brm3)
prior_summary(fCORT_pbvid_nool_brm3)

#Posterior predictive checks
pp_check(fCORT_pbvid_brm3, type = "dens_overlay", ndraws = 1000)
pp_check(fCORT_pbvid_nool_brm3, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_pbvid_brm3)
launch_shinystan(fCORT_pbvid_brm3)

plot(fCORT_pbvid_nool_brm3)
launch_shinystan(fCORT_pbvid_nool_brm3)

#Posterior distribution
as_draws_df(fCORT_pbvid_brm3)
mcmc_areas(fCORT_pbvid_brm3)
mcmc_intervals(fCORT_pbvid_brm3)

as_draws_df(fCORT_pbvid_nool_brm3)
mcmc_areas(fCORT_pbvid_nool_brm3)
mcmc_intervals(fCORT_pbvid_nool_brm3)

#Evaluation and interpretation
fCORT_pbvid_brm3_loo <- loo(fCORT_pbvid_brm3)
loo_compare(fCORT_pbvid_brm1_loo, fCORT_pbvid_brm3_loo)
loo_compare(fCORT_pbvid_brm1_loo, fCORT_pbvid_brm2_loo, fCORT_pbvid_brm3_loo)

fitted(fCORT_pbvid_brm3, scale = "response")
conditional_effects(fCORT_pbvid_brm3)

fCORT_pbvid_nool_brm3_loo <- loo(fCORT_pbvid_nool_brm3)
loo_compare(fCORT_pbvid_nool_brm1_loo, fCORT_pbvid_nool_brm3_loo)
fitted(fCORT_pbvid_nool_brm3, scale = "response")
conditional_effects(fCORT_pbvid_nool_brm3)

bayes_R2(fCORT_pbvid_brm3)
emmeans(fCORT_pbvid_brm3, ~ site, type = "response") |> pairs()
bayes_R2(fCORT_pbvid_brm3)
emmeans(fCORT_pbvid_brm3, ~ site, type = "response") |> pairs()

hypothesis(fCORT_pbvid_brm3, "pbvid_mean_z < 0") 
hypothesis(fCORT_pbvid_brm3, "pbvid_mean_z < 0") 


#Pair-bond model 4: pair-bond strength * sex + site + body condition

#Model just with prior
default_prior()

fCORT_pbvid_brm4_prior <- brm(CORT_pg.mg ~ pbvid_mean_z * sex + site + body_cond_z + 
                                (1|feather_year_grown) + (1|JID) + 
                                (1|pair_ID_yr_sampled),
                              data = fCORT_pbvid, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(1.25, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_pbvid_brm4_prior, ndraws = 100) +
  scale_x_log10()

#model including extreme values
fCORT_pbvid_brm4 <- brm(CORT_pg.mg ~ pbvid_mean_z * sex + site + body_cond_z + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma")),
                        save_pars = save_pars(all = TRUE)
)

#model excluding extreme values
fCORT_pbvid_nool_brm4 <- brm(CORT_pg.mg ~ pbvid_mean_z * sex + site + body_cond_z + 
                               (1|feather_year_grown) + (1|JID) + 
                               (1|pair_ID_yr_sampled),
                             data = fCORT_pbvid_nool, 
                             family = lognormal(link = "identity"),
                             prior = c(
                               prior(normal(1.25, 0.3), class = "Intercept"),  
                               prior(normal(0, 0.3), class = "b"),
                               prior(lognormal(log(0.35), 0.3), class = "sigma"))            
)

#Model summary
summary(fCORT_pbvid_brm4)
summary(fCORT_pbvid_nool_brm4)

#Check collinearity
check_collinearity(fCORT_pbvid_brm4)
check_collinearity(fCORT_pbvid_nool_brm4)

#Prior checks
prior_summary(fCORT_pbvid_brm4)
prior_summary(fCORT_pbvid_nool_brm4)

#Posterior predictive checks
pp_check(fCORT_pbvid_brm4, type = "dens_overlay", ndraws = 1000)
pp_check(fCORT_pbvid_nool_brm4, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_pbvid_brm4)
launch_shinystan(fCORT_pbvid_brm4)

plot(fCORT_pbvid_nool_brm4)
launch_shinystan(fCORT_pbvid_nool_brm4)

#Posterior distribution
as_draws_df(fCORT_pbvid_brm4)
mcmc_areas(fCORT_pbvid_brm4)
mcmc_intervals(fCORT_pbvid_brm4)

as_draws_df(fCORT_pbvid_nool_brm4)
mcmc_areas(fCORT_pbvid_nool_brm4)
mcmc_intervals(fCORT_pbvid_nool_brm4)

#Evaluation and interpretation
fCORT_pbvid_brm4_loo <- loo(fCORT_pbvid_brm4)
loo_compare(fCORT_pbvid_brm1_loo, fCORT_pbvid_brm4_loo)
fitted(fCORT_pbvid_brm4, scale = "response")
conditional_effects(fCORT_pbvid_brm4)

fCORT_pbvid_nool_brm4_loo <- loo(fCORT_pbvid_nool_brm4)
loo_compare(fCORT_pbvid_nool_brm1_loo, fCORT_pbvid_nool_brm4_loo)
fitted(fCORT_pbvid_nool_brm4, scale = "response")
conditional_effects(fCORT_pbvid_nool_brm4)

bayes_R2(fCORT_pbvid_brm4)
emmeans(fCORT_pbvid_brm4, ~ site, type = "response") |> pairs()
bayes_R2(fCORT_pbvid_brm4)
emmeans(fCORT_pbvid_nool_brm4, ~ site, type = "response") |> pairs()

hypothesis(fCORT_pbvid_brm4, "pbvid_mean_z < 0") 
hypothesis(fCORT_pbvid_nool_brm4, "pbvid_mean_z < 0") 


#Pair-bond model 5: pair-bond strength * age + site + body condition

#Model just with prior
default_prior()

fCORT_pbvid_brm5_prior <- brm(CORT_pg.mg ~ pbvid_mean_z * age_z + site + body_cond_z + 
                                (1|feather_year_grown) + (1|JID) + 
                                (1|pair_ID_yr_sampled),
                              data = fCORT_pbvid, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(1.25, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_pbvid_brm5_prior, ndraws = 100) +
  scale_x_log10()

#model including extreme values
fCORT_pbvid_brm5 <- brm(CORT_pg.mg ~ pbvid_mean_z * age_z + site + body_cond_z + 
                          (1|feather_year_grown) + (1|JID) + 
                          (1|pair_ID_yr_sampled),
                        data = fCORT_pbvid, 
                        family = lognormal(link = "identity"),
                        prior = c(
                          prior(normal(1.25, 0.3), class = "Intercept"),  
                          prior(normal(0, 0.3), class = "b"),
                          prior(lognormal(log(0.35), 0.3), class = "sigma")),
                        save_pars = save_pars(all = TRUE)
)

#model excluding extreme values
fCORT_pbvid_nool_brm5 <- brm(CORT_pg.mg ~ pbvid_mean_z * age_z + site + body_cond_z + 
                               (1|feather_year_grown) + (1|JID) + 
                               (1|pair_ID_yr_sampled),
                             data = fCORT_pbvid_nool, 
                             family = lognormal(link = "identity"),
                             prior = c(
                               prior(normal(1.25, 0.3), class = "Intercept"),  
                               prior(normal(0, 0.3), class = "b"),
                               prior(lognormal(log(0.35), 0.3), class = "sigma"))            
)

#Model summary
summary(fCORT_pbvid_brm5)
summary(fCORT_pbvid_nool_brm5)

#Check collinearity
check_collinearity(fCORT_pbvid_brm5)
check_collinearity(fCORT_pbvid_nool_brm5)

#Prior checks
prior_summary(fCORT_pbvid_brm5)
prior_summary(fCORT_pbvid_nool_brm5)

#Posterior predictive checks
pp_check(fCORT_pbvid_brm5, type = "dens_overlay", ndraws = 1000)
pp_check(fCORT_pbvid_nool_brm5, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_pbvid_brm5)
launch_shinystan(fCORT_pbvid_brm5)

plot(fCORT_pbvid_nool_brm5)
launch_shinystan(fCORT_pbvid_nool_brm5)

#Posterior distribution
as_draws_df(fCORT_pbvid_brm5)
mcmc_areas(fCORT_pbvid_brm5)
mcmc_intervals(fCORT_pbvid_brm5)

as_draws_df(fCORT_pbvid_nool_brm5)
mcmc_areas(fCORT_pbvid_nool_brm5)
mcmc_intervals(fCORT_pbvid_nool_brm5)

#Evaluation and interpretation
fCORT_pbvid_brm5_loo <- loo(fCORT_pbvid_brm5)
loo_compare(fCORT_pbvid_brm1_loo, fCORT_pbvid_brm5_loo)
fitted(fCORT_pbvid_brm5, scale = "response")
conditional_effects(fCORT_pbvid_brm5)

fCORT_pbvid_nool_brm5_loo <- loo(fCORT_pbvid_nool_brm5)
loo_compare(fCORT_pbvid_nool_brm1_loo, fCORT_pbvid_nool_brm5_loo)
fitted(fCORT_pbvid_nool_brm5, scale = "response")
conditional_effects(fCORT_pbvid_nool_brm5)

bayes_R2(fCORT_pbvid_brm5)
emmeans(fCORT_pbvid_brm5, ~ site, type = "response") |> pairs()
bayes_R2(fCORT_pbvid_brm5)
emmeans(fCORT_pbvid_nool_brm5, ~ site, type = "response") |> pairs()

hypothesis(fCORT_pbvid_brm5, "pbvid_mean_z < 0") 
hypothesis(fCORT_pbvid_nool_brm5, "pbvid_mean_z < 0") 


#Pair bond plot for manuscript 
pairbond_plot <- ggplot(fCORT_pbvid, aes(x= pbvid_mean, y = CORT_pg.mg)) +
  geom_point(alpha = 0.3, size=2, shape= 16) +
  theme_few(base_size = 14) +
  xlab("Pair-bond strength (proportion of time spent together)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,10) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2.5,5, 7.5, 10.0)) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
pairbond_plot

fCORT_pbvid_brm1 %>%
  spread_draws(b_Intercept, b_pbvid_mean_z) %>%
  ggplot(aes(y = , x = b_pbvid_mean_z)) +
  theme_classic(base_size = 35) +
  theme(legend.position="none", text = element_text(size = 28, family = "Garamond")) +
  labs(x = "Pair-bond strength", y = "Density") +
  stat_halfeye()

cond <- conditional_effects(fCORT_pbvid_brm1, effect = "pbvid_mean_z")
p <- plot(cond, points=F)
pairbond_plot <- p[[1]] + 
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  scale_y_continuous(limits = c(0, 10)) +
  geom_point(
    aes(x = pbvid_mean_z, y = CORT_pg.mg), 
    data = fCORT_pbvid, 
    color = "black",
    size = 2,
    alpha = 0.3,
    shape = 16,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Pair-bond strength", y = "fCORT (pg/mg)")
pairbond_plot


#4. Co-regulation ----
mean(fCORT_coreg$F)
mean(fCORT_coreg$M)
sd(fCORT_coreg$F)
sd(fCORT_coreg$M)
summary(fCORT_coreg$F)
summary(fCORT_coreg$M)

#Add year feather was grown
fCORT_coreg$year <- fCORT_coreg$year - 1

#Residualised male and female fCORT
fCORT_coreg$female_CORT_resid <- resid(lm(female_CORT ~ site + female_body_cond, data = fCORT_coreg))
fCORT_coreg$male_CORT_resid <- resid(lm(male_CORT ~ site + male_body_cond, data = fCORT_coreg))

#Log-transformed male and female fCORT
fCORT_coreg$female_CORT_log <- log(fCORT_coreg$female_CORT)
fCORT_coreg$male_CORT_log <- log(fCORT_coreg$male_CORT)

#Scaled variables
fCORT_coreg$female_CORT_z <- scale(fCORT_coreg$female_CORT)
fCORT_coreg$male_CORT_z <- scale(fCORT_coreg$male_CORT)
fCORT_coreg$female_body_cond_z <- scale(fCORT_coreg$female_body_cond)
fCORT_coreg$male_body_cond_z <- scale(fCORT_coreg$male_body_cond)

#Drop site MW with 0 observations
fCORT_coreg$site <- droplevels(fCORT_coreg$site)

#Correct one pair ID
fCORT_coreg$pair_ID[fCORT_coreg$pair_ID == "J4646 J4645"] <- "J4645 J4646"

#Subset by site 
fCORT_coreg_Y <- subset(fCORT_coreg, site == "Y")
fCORT_coreg_Z <- subset(fCORT_coreg, site == "Z")

#Remove outliers with high fCORT levels
fCORT_coreg_nool <- subset(fCORT_coreg, fCORT_coreg$female_CORT < 8)
fCORT_coreg_nool <- subset(fCORT_coreg_nool, fCORT_coreg_nool$male_CORT < 8)

#Only feathers sampled in same year
fCORT_coreg2 <- subset(fCORT_coreg, fCORT_coreg$year_diff ==0)
fCORT_coreg2 <- subset(fCORT_coreg2, fCORT_coreg2$year_diff ==0)

#Remove site X
fCORT_coreg2 <- subset(fCORT_coreg, !fCORT_coreg$site == "X")
fCORT_coreg2$site <- droplevels(fCORT_coreg2$site)

#Model 

# Separate formulas for each trait
bf_fCORT_coreg_f <- bf(female_CORT_log ~ factor(site) + female_body_cond_z + (1|pair_ID) + (1|female_CORT_year))
bf_fCORT_coreg_m   <- bf(male_CORT_log   ~ factor(site) + male_body_cond_z + (1|pair_ID) + (1|male_CORT_year))

bf_fCORT_coreg_f <- bf(female_CORT_log ~ factor(site) + (1|pair_ID) + (1|female_CORT_year))
bf_fCORT_coreg_m   <- bf(male_CORT_log   ~ factor(site) + (1|pair_ID) + (1|male_CORT_year))

#Model just with prior
default_prior()

fCORT_coreg_brm1_prior <- brm(bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE), 
                              family = gaussian(link = "identity"),
                              data = fCORT_coreg,
                              prior = c(
                                prior(normal(0, 0.3), class = "b", resp = "femaleCORTlog"),
                                prior(normal(0, 0.3), class = "b", resp = "maleCORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "femaleCORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "maleCORTlog"),
                                prior(student_t(3, 0, 0.35), class = "sigma", resp = "femaleCORTlog"),
                                prior(student_t(3, 0, 0.35), class = "sigma", resp = "maleCORTlog"),
                                prior(lkj(2), class = "rescor")),                              
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_coreg_brm1_prior, resp = "femaleCORTlog", ndraws = 100) 
pp_check(fCORT_coreg_brm1_prior, resp = "maleCORTlog", ndraws = 100) 

#model including extreme values
#model for manuscript
fCORT_coreg_brm1 <- brm(bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE), 
                              family = gaussian(link = "identity"),
                              data = fCORT_coreg,
                              prior = c(
                                prior(normal(0, 0.3), class = "b", resp = "femaleCORTlog"),
                                prior(normal(0, 0.3), class = "b", resp = "maleCORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "femaleCORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "maleCORTlog"),
                                prior(student_t(3, 0, 0.35), class = "sigma", resp = "femaleCORTlog"),
                                prior(student_t(3, 0, 0.35), class = "sigma", resp = "maleCORTlog"),
                                prior(lkj(2), class = "rescor")),
                                control = list(adapt_delta = 0.99)
)

fCORT_coreg_brm1 <- brm(bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE), 
                        family = gaussian(link = "identity"),
                        data = fCORT_coreg,
                        prior = c(
                          prior(normal(0, 0.3), class = "b", resp = "femaleCORTlog"),
                          prior(normal(0, 0.3), class = "b", resp = "maleCORTlog"),
                          prior(normal(1.25, 0.3), class = "Intercept", resp = "femaleCORTlog"),
                          prior(normal(1.25, 0.3), class = "Intercept", resp = "maleCORTlog"),
                          prior(student_t(7, 0, 0.35), class = "sigma", resp = "femaleCORTlog"),
                          prior(student_t(7, 0, 0.35), class = "sigma", resp = "maleCORTlog"),
                          prior(lkj(2), class = "rescor")),
                        control = list(adapt_delta = 0.99)
)

fCORT_coreg_brm1 <- brm(bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE), 
                        family = gaussian(link = "identity"),
                        data = fCORT_coreg,
                        prior = c(
                          prior(normal(0, 0.3), class = "b", resp = "femaleCORTlog"),
                          prior(normal(0, 0.3), class = "b", resp = "maleCORTlog"),
                          prior(normal(1.25, 0.3), class = "Intercept", resp = "femaleCORTlog"),
                          prior(normal(1.25, 0.3), class = "Intercept", resp = "maleCORTlog"),
                          prior(normal(0, 0.35), class = "sigma", resp = "femaleCORTlog"),
                          prior(normal(0, 0.35), class = "sigma", resp = "maleCORTlog"),
                          prior(lkj(2), class = "rescor")),
                        control = list(adapt_delta = 0.99)
)

#model excluding extreme values
fCORT_coreg_nool_brm1 <- brm(bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE), 
                        family = gaussian(link = "identity"),
                        data = fCORT_coreg_nool,
                        prior = c(
                          prior(normal(0, 0.3), class = "b", resp = "femaleCORTlog"),
                          prior(normal(0, 0.3), class = "b", resp = "maleCORTlog"),
                          prior(normal(1.25, 0.3), class = "Intercept", resp = "femaleCORTlog"),
                          prior(normal(1.25, 0.3), class = "Intercept", resp = "maleCORTlog"),
                          prior(student_t(3, 0, 0.35), class = "sigma", resp = "femaleCORTlog"),
                          prior(student_t(3, 0, 0.35), class = "sigma", resp = "maleCORTlog"),
                          prior(lkj(2), class = "rescor")),
                        control = list(adapt_delta = 0.99)
)

#Model summary
summary(fCORT_coreg_brm1)
summary(fCORT_coreg_nool_brm1)

#Check collinearity
check_collinearity(fCORT_coreg_brm1)
check_collinearity(fCORT_coreg_nool_brm1)

#Prior checks
prior_summary(fCORT_coreg_brm1)
prior_summary(fCORT_coreg_nool_brm1)

#Posterior predictive checks
pp_check(fCORT_coreg_brm1, type = "dens_overlay", resp = "femaleCORTlog", ndraws = 1000)
pp_check(fCORT_coreg_brm1, type = "dens_overlay", resp = "maleCORTlog", ndraws = 1000)

pp_check(fCORT_coreg_nool_brm1, type = "dens_overlay", resp = "femaleCORTlog", ndraws = 1000)
pp_check(fCORT_coreg_nool_brm1, type = "dens_overlay", resp = "maleCORTlog", ndraws = 1000)

#Plot model
plot(fCORT_coreg_brm1)
launch_shinystan(fCORT_coreg_brm1)

plot(fCORT_coreg_brm1)
launch_shinystan(fCORT_coreg_nool_brm1)

#Posterior distribution
as_draws_df(fCORT_coreg_brm1)
mcmc_areas(fCORT_coreg_brm1)
mcmc_intervals(fCORT_coreg_brm1)

as_draws_df(fCORT_coreg_nool_brm1)
mcmc_areas(fCORT_coreg_nool_brm1)
mcmc_intervals(fCORT_coreg_nool_brm1)

#Evaluation and interpretation
fCORT_coreg_brm1 <- loo(fCORT_coreg_brm1)
fitted(fCORT_coreg_brm1, scale = "response")
conditional_effects(fCORT_coreg_brm1)

fCORT_coreg_nool_brm1_loo <- loo(fCORT_coreg_nool_brm1)
fitted(fCORT_coreg_nool_brm1, scale = "response")
conditional_effects(fCORT_coreg_nool_brm1)

bayes_R2(fCORT_coreg_brm1)
bayes_R2(fCORT_coreg_brm1)

posterior <- as_draws_df(fCORT_coreg_brm1)
posterior$rescor__femaleCORTlog__maleCORTlog
rescor_samples <- posterior$rescor__femaleCORTlog__maleCORTlog
prob_positive <- mean(rescor_samples > 0)
prob_positive

posterior <- as_draws_df(fCORT_coreg_nool_brm1)
posterior$rescor__femaleCORTlog__maleCORTlog
rescor_samples <- posterior$rescor__femaleCORTlog__maleCORTlog
prob_positive <- mean(rescor_samples > 0)
prob_positive

#Coregulation plot for manuscript
coreg_plot <- ggplot(fCORT_coreg, aes(x= male_CORT, y = female_CORT)) +
  geom_point(alpha = 0.3, size=2, shape= 16) +
  theme_few(base_size = 14) +
  xlab("Male fCORT (pg/mg)") +
  ylab("Female fCORT (pg/mg)") +
  xlim(0,10) +
  ylim(0,10) +
  scale_x_continuous(limits = c(0, 10), breaks = c(0,2.5,5, 7.5, 10)) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2.5,5, 7.5, 10)) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
coreg_plot

pair_fCORT_plot <- ggarrange(pairbond_plot, coreg_plot, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),  
                             font.label = list(family = "Garamond", face = "plain", size = 14, color = "black"), 
                             heights = c(1, 1), widths = c(1, 1))
pair_fCORT_plot


fCORT_coreg_brm1 %>%
  spread_draws(rescor__femaleCORTlog__maleCORTlog) %>%
  ggplot(aes(x = rescor__femaleCORTlog__maleCORTlog, y = 0)) +  # y=0 just puts it on one line
  stat_halfeye() +
  theme_classic(base_size = 35) +
  theme(legend.position = "none", text = element_text(size = 28, family = "Garamond")) +
  labs(x = "Residual correlation of fCORT", y = "Density")

library(tidybayes)

# Create a grid of female CORT values
grid <- data.frame(
  female_body_cond_z = mean(fCORT_coreg$female_body_cond_z),
  male_body_cond_z = mean(fCORT_coreg$male_body_cond_z),
  site = "X" # choose a site
)

# Get posterior predictions for male CORT given female values
pp <- add_predicted_draws(fCORT_coreg, fCORT_coreg_brm1, n = 500)

# Plot male vs female with posterior uncertainty
ggplot(pp, aes(x = female_CORT_log, y = .prediction)) +
  stat_halfeye(aes(y = .prediction), .width = c(0.5, 0.9)) +
  geom_point(data = fCORT_coreg, aes(y = male_CORT_log), color = "black") +
  geom_smooth(aes(y = .prediction), method = "lm", color = "blue")


post <- as_draws_df(fCORT_coreg_brm1)
# Fitted values (posterior mean) for female
fitted_f <- fitted(fCORT_coreg_brm1, resp = "femaleCORTlog")[, "Estimate"]
# Fitted values (posterior mean) for male
fitted_m <- fitted(fCORT_coreg_brm1, resp = "maleCORTlog")[, "Estimate"]

resid_f <- fCORT_coreg$female_CORT_log - fitted_f
resid_m <- fCORT_coreg$male_CORT_log - fitted_m

ggplot(data.frame(resid_f, resid_m), aes(x = resid_f, y = resid_m)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  theme_minimal() +
  labs(x = "Female residual CORT", y = "Male residual CORT")

 ggplot(data.frame(resid_f, resid_m), aes(x = resid_f, y = resid_m)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()

 
#5. 2022 social network ----

fCORT_22_sna <- subset(fCORT_22, !is.na(fCORT_22$degree))
fCORT_22_sna$social_diff <- am_data$social_diff[match(fCORT_22_sna$JID, am_data$JID)]
fCORT_22_sna$site <- ifelse(fCORT_22_sna$site == "Y", "Y", "Z")

fCORT_22_sna$pair_ID_yr_sampled[fCORT_22_sna$pair_ID_yr_sampled == "J4646 J4645"] <- "J4645 J4646"

fCORT_22_sna$degree_z <- scale(fCORT_22_sna$degree)
fCORT_22_sna$strength_z <- scale(fCORT_22_sna$strength)
fCORT_22_sna$social_diff_z <- scale(fCORT_22_sna$social_diff)
fCORT_22_sna$transitivity_z <- scale(fCORT_22_sna$transitivity)
fCORT_22_sna$visit_number_z <- scale(fCORT_22_sna$visit_number)
fCORT_22_sna$body_cond_z <- scale(fCORT_22_sna$body_cond)

#Model 1 Degree 

#Model just with prior
default_prior()

fCORT_sna_brm1_prior <- brm(CORT_pg.mg ~ degree_z  + visit_number_z +
                              (1|pair_ID_yr_sampled),
                              data = fCORT_22_sna, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(1.25, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_sna_brm1_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_sna_brm1 <- brm(CORT_pg.mg ~ degree_z  + visit_number_z +
                            site + body_cond_z +                        
                            (1|pair_ID_yr_sampled),
                            data = fCORT_22_sna, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(1.25, 0.3), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),
                      control = list(adapt_delta = 0.99)
                      
)
                            
#Model summary
summary(fCORT_sna_brm1)

#Check collinearity
check_collinearity(fCORT_sna_brm1)

#Prior checks 
prior_summary(fCORT_sna_brm1)

#Posterior predictive checks
pp_check(fCORT_sna_brm1, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_sna_brm1)
launch_shinystan(fCORT_sna_brm1)

#Posterior distribution
as_draws_df(fCORT_sna_brm1)
mcmc_areas(fCORT_sna_brm1)
mcmc_intervals(fCORT_sna_brm1)

#Evaluation and interpretation
fCORT_sna_brm1_loo <- loo(fCORT_sna_brm1, moment_match = TRUE)
fCORT_sna_brm1_loo 

fitted(fCORT_sna_brm1, scale = "response")
conditional_effects(fCORT_sna_brm1)

bayes_R2(fCORT_sna_brm1)
emmeans(fCORT_sna_brm1, ~ site, type = "response") |> pairs()

hypothesis(fCORT_sna_brm1, "degree_z < 0") 


#Model 2 Strength 

#Model just with prior
default_prior()

fCORT_sna_brm2_prior <- brm(CORT_pg.mg ~ strength_z  + visit_number_z +
                              (1|pair_ID_yr_sampled),
                            data = fCORT_22_sna, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(1.25, 0.3), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                            sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_sna_brm2_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_sna_brm2 <- brm(CORT_pg.mg ~ strength_z  + visit_number_z +
                        site + body_cond_z +                        
                        (1|pair_ID_yr_sampled),
                      data = fCORT_22_sna, 
                      family = lognormal(link = "identity"),
                      prior = c(
                        prior(normal(1.25, 0.3), class = "Intercept"),  
                        prior(normal(0, 0.3), class = "b"),
                        prior(lognormal(log(0.35), 0.3), class = "sigma")),
                      control = list(adapt_delta = 0.99)
                      
)

#Model summary
summary(fCORT_sna_brm2)

#Check collinearity
check_collinearity(fCORT_sna_brm2)

#Prior checks 
prior_summary(fCORT_sna_brm2)

#Posterior predictive checks
pp_check(fCORT_sna_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_sna_brm2)
launch_shinystan(fCORT_sna_brm2)

#Posterior distribution
as_draws_df(fCORT_sna_brm2)
mcmc_areas(fCORT_sna_brm2)
mcmc_intervals(fCORT_sna_brm2)

#Evaluation and interpretation
fCORT_sna_brm2_loo <- loo(fCORT_sna_brm2, moment_match = TRUE)
fCORT_sna_brm2_loo 

fitted(fCORT_sna_brm2, scale = "response")
conditional_effects(fCORT_sna_brm2)

bayes_R2(fCORT_sna_brm2)

hypothesis(fCORT_sna_brm2, "strength_z < 0") 


#Model 3 Social differentiation

#Model just with prior
default_prior()

fCORT_sna_brm3_prior <- brm(CORT_pg.mg ~ social_diff_z  + visit_number_z +
                              (1|pair_ID_yr_sampled),
                            data = fCORT_22_sna, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(1.25, 0.3), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                            sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_sna_brm3_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_sna_brm3 <- brm(CORT_pg.mg ~ social_diff_z  + visit_number_z +
                        site + body_cond_z +                        
                        (1|pair_ID_yr_sampled),
                      data = fCORT_22_sna, 
                      family = lognormal(link = "identity"),
                      prior = c(
                        prior(normal(1.25, 0.3), class = "Intercept"),  
                        prior(normal(0, 0.3), class = "b"),
                        prior(lognormal(log(0.35), 0.3), class = "sigma")),
                      control = list(adapt_delta = 0.99)
)

#Model summary
summary(fCORT_sna_brm3)

#Check collinearity
check_collinearity(fCORT_sna_brm3)

#Prior checks 
prior_summary(fCORT_sna_brm3)

#Posterior predictive checks
pp_check(fCORT_sna_brm3, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_sna_brm3)
launch_shinystan(fCORT_sna_brm3)

#Posterior distribution
as_draws_df(fCORT_sna_brm3)
mcmc_areas(fCORT_sna_brm3)
mcmc_intervals(fCORT_sna_brm3)

#Evaluation and interpretation
fCORT_sna_brm3_loo <- loo(fCORT_sna_brm3, moment_match = TRUE)
fCORT_sna_brm3_loo 

fitted(fCORT_sna_brm3, scale = "response")
conditional_effects(fCORT_sna_brm3)

bayes_R2(fCORT_sna_brm3)

hypothesis(fCORT_sna_brm3, "social_diff_z < 0") 


#Model 4 Local transitivity

#Model just with prior
default_prior()

fCORT_sna_brm4_prior <- brm(CORT_pg.mg ~ transitivity_z  + visit_number_z +
                              (1|pair_ID_yr_sampled),
                            data = fCORT_22_sna, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(1.25, 0.3), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                            sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_sna_brm4_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_sna_brm4 <- brm(CORT_pg.mg ~ transitivity_z  + visit_number_z +
                        site + body_cond_z +                        
                        (1|pair_ID_yr_sampled),
                      data = fCORT_22_sna, 
                      family = lognormal(link = "identity"),
                      prior = c(
                        prior(normal(1.25, 0.3), class = "Intercept"),  
                        prior(normal(0, 0.3), class = "b"),
                        prior(lognormal(log(0.35), 0.3), class = "sigma")),
                      control = list(adapt_delta = 0.99)
)

#Model summary
summary(fCORT_sna_brm4)

#Check collinearity
check_collinearity(fCORT_sna_brm4)

#Prior checks 
prior_summary(fCORT_sna_brm4)

#Posterior predictive checks
pp_check(fCORT_sna_brm4, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_sna_brm4)
launch_shinystan(fCORT_sna_brm4)

#Posterior distribution
as_draws_df(fCORT_sna_brm4)
mcmc_areas(fCORT_sna_brm4)
mcmc_intervals(fCORT_sna_brm4)

#Evaluation and interpretation
fCORT_sna_brm4_loo <- loo(fCORT_sna_brm4, moment_match = TRUE)
fCORT_sna_brm4_loo 

fitted(fCORT_sna_brm4, scale = "response")
conditional_effects(fCORT_sna_brm4)

bayes_R2(fCORT_sna_brm4)

hypothesis(fCORT_sna_brm4, "transitivity_z < 0") 

loo_compare(fCORT_sna_brm1_loo, fCORT_sna_brm2_loo, fCORT_sna_brm3_loo, fCORT_sna_brm4_loo)


#Coregulation in social network ----
dyad_edge_CORT$ind1_site <- fCORT$site[match(dyad_edge_CORT$ind1, fCORT$JID)]
dyad_edge_CORT$ind2_site <- fCORT$site[match(dyad_edge_CORT$ind2, fCORT$JID)]
dyad_edge_CORT$ind1_site <- droplevels(dyad_edge_CORT$ind1_site)
dyad_edge_CORT$ind2_site <- droplevels(dyad_edge_CORT$ind2_site)

dyad_edge_CORT$ind1_body_cond <- fCORT$body_cond[match(dyad_edge_CORT$ind1, fCORT$JID)]
dyad_edge_CORT$ind2_body_cond <- fCORT$body_cond[match(dyad_edge_CORT$ind2, fCORT$JID)]
dyad_edge_CORT$ind1_body_cond_z <- scale(dyad_edge_CORT$ind1_body_cond)
dyad_edge_CORT$ind2_body_cond_z <- scale(dyad_edge_CORT$ind2_body_cond)

dyad_edge_CORT$ind1 <- as.factor(dyad_edge_CORT$ind1)
dyad_edge_CORT$ind2 <- as.factor(dyad_edge_CORT$ind2)

dyad_edge_CORT$ind1_CORT_log <- log(dyad_edge_CORT$ind1_CORT)
dyad_edge_CORT$ind2_CORT_log <- log(dyad_edge_CORT$ind2_CORT)

dyad_edge_CORT$g_edge_weights_z <- scale(dyad_edge_CORT$g_edge_weights)

fCORT_coreg_brm2_prior <- brm(ind1ind2_CORT_diff_norm ~ g_edge_weights_z + 
                              (1|mm(ind1,ind2)),
                              data = dyad_edge_CORT, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(0.4, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma")),            
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_coreg_brm2_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_coreg_brm2 <- brm(ind1ind2_CORT_diff_norm ~ g_edge_weights_z + 
                                (1|mm(ind1,ind2)),
                              data = dyad_edge_CORT, 
                              family = lognormal(link = "identity"),
                              prior = c(
                                prior(normal(0.4, 0.3), class = "Intercept"),  
                                prior(normal(0, 0.3), class = "b"),
                                prior(lognormal(log(0.35), 0.3), class = "sigma"))
)


#Model summary
summary(fCORT_coreg_brm2)

#Check collinearity
check_collinearity(fCORT_coreg_brm2)

#Posterior predictive checks
pp_check(fCORT_coreg_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_coreg_brm2)
launch_shinystan(fCORT_coreg_brm2)

#Posterior distribution
as_draws_df(fCORT_coreg_brm2)
mcmc_areas(fCORT_coreg_brm2)
mcmc_intervals(fCORT_coreg_brm2)

#Evaluation and interpretation
fCORT_coreg_brm2_loo <- loo(fCORT_coreg_brm2, moment_match = TRUE)
fitted(fCORT_coreg_brm2, scale = "response")
conditional_effects(fCORT_coreg_brm2)

bayes_R2(fCORT_coreg_brm2)

hypothesis(fCORT_coreg_brm2, "g_edge_weights_z > 0") 


#Bivariate model? Convergence issues.

# Separate formulas for each trait
bf_fCORT_coreg_ind1 <- bf(ind1_CORT_log ~ factor(ind1_site) + ind1_body_cond_z + (1|mm(ind1, ind2)))
bf_fCORT_coreg_ind2   <- bf(ind2_CORT_log ~ factor(ind2_site) + ind2_body_cond_z + (1|mm(ind1, ind2)))

bf_fCORT_coreg_ind1 <- bf(ind1_CORT_log ~ factor(ind1_site) + ind1_body_cond_z + (1|ind1))
bf_fCORT_coreg_ind2 <- bf(ind2_CORT_log ~ factor(ind2_site) + ind2_body_cond_z + (1|ind2))

#Model just with prior
default_prior()

fCORT_coreg_brm2_prior <- brm(bf_fCORT_coreg_ind1 + bf_fCORT_coreg_ind2 + set_rescor(TRUE), 
                              family = gaussian(link = "identity"),
                              data = dyad_edge_CORT,
                              prior = c(
                                prior(normal(0, 0.3), class = "b", resp = "ind1CORTlog"),
                                prior(normal(0, 0.3), class = "b", resp = "ind2CORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "ind1CORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "ind2CORTlog"),
                                prior(student_t(3, 0, 0.35), class = "sigma", resp = "ind1CORTlog"),
                                prior(student_t(3, 0, 0.35), class = "sigma", resp = "ind2CORTlog"),
                                prior(lkj(2), class = "rescor")),                              
                              sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_coreg_brm2_prior, resp = "ind1CORTlog", ndraws = 100) 
pp_check(fCORT_coreg_brm2_prior, resp = "ind2CORTlog", ndraws = 100) 

#Model
fCORT_coreg_brm2 <- brm(bf_fCORT_coreg_ind1 + bf_fCORT_coreg_ind2 + set_rescor(TRUE), 
                              family = gaussian(link = "identity"),
                              data = dyad_edge_CORT,
                              prior = c(
                                prior(normal(0, 0.3), class = "b", resp = "ind1CORTlog"),
                                prior(normal(0, 0.3), class = "b", resp = "ind2CORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "ind1CORTlog"),
                                prior(normal(1.25, 0.3), class = "Intercept", resp = "ind2CORTlog"),
                                prior(student_t(3, 0, 0.25), class = "sigma", resp = "ind1CORTlog"),
                                prior(student_t(3, 0, 0.25), class = "sigma", resp = "ind2CORTlog"),
                                prior(lkj(6), class = "rescor")),
                                prior(student_t(3, 0, 0.15), class = "sd"),
                        iter = 6000, warmup = 2000,
                        control = list(adapt_delta = 0.90, max_treedepth = 12)
)

fCORT_coreg_brm2 <- brm(bf_fCORT_coreg_ind1 + bf_fCORT_coreg_ind2 + set_rescor(TRUE),
                            family = gaussian(),
                            data = dyad_edge_CORT,
                            prior = c(
                              prior(normal(0, 0.3), class = "b", resp = "ind1CORTlog"),
                              prior(normal(0, 0.3), class = "b", resp = "ind2CORTlog"),
                              prior(normal(1.25, 0.3), class = "Intercept", resp = "ind1CORTlog"),
                              prior(normal(1.25, 0.3), class = "Intercept", resp = "ind2CORTlog"),
                              prior(student_t(3, 0, 0.25), class = "sigma", resp = "ind1CORTlog"),
                              prior(student_t(3, 0, 0.25), class = "sigma", resp = "ind2CORTlog"),
                              prior(lkj(2), class = "rescor")),
                            iter = 6000, warmup = 2000,
                            control = list(adapt_delta = 0.95, max_treedepth = 15)
)


#Model summary
summary(fCORT_coreg_brm2)

#Check collinearity
check_collinearity(fCORT_coreg_brm2)

#Prior checks
prior_summary(fCORT_coreg_brm2)

#Posterior predictive checks
pp_check(fCORT_coreg_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_coreg_brm2)
launch_shinystan(fCORT_coreg_brm2)

#Posterior distribution
as_draws_df(fCORT_coreg_brm2)
mcmc_areas(fCORT_coreg_brm2)
mcmc_intervals(fCORT_coreg_brm2)

#Evaluation and interpretation
fCORT_coreg_brm2 <- loo(fCORT_coreg_brm1)
fitted(fCORT_coreg_brm2, scale = "response")
conditional_effects(fCORT_coreg_brm2)

bayes_R2(fCORT_coreg_brm2)

posterior <- as_draws_df(fCORT_coreg_brm2)
posterior$rescor__ind1CORTlog__ind2CORTlog
rescor_samples <- posterior$rescor__ind1CORTlog__ind2CORTlog
prob_positive <- mean(rescor_samples > 0)
prob_positive









#PCA on social network metrics
all_individuals3 <- all_individuals2
all_individuals3$JID <- NULL
all_individuals3$visit_number <- NULL
all_individuals3$visit_duration <- NULL

all_individuals3 <- na.omit(all_individuals3)

all_individuals3$mean_strength <- NULL

#normalise data
norm_data <- scale(all_individuals3)

#correlation matrix
corr_matrix <- cor(norm_data)
corr_matrix

#correlation plot
ggcorrplot(corr_matrix)

#PCA diagnostics
eigen(corr_matrix)
cortest.bartlett(corr_matrix, n= 309) #correlation
cortest.mat(corr_matrix, n1= 309)
KMO(corr_matrix)

factors(corr_matrix, n.obs = 310)
fa.parallel(corr_matrix, n.obs = 310)

#PCA
data.pca <- princomp(norm_data)
data.pca <- prcomp(norm_data)

#Create dataframe with PCs
data.pca <- prcomp(norm_data, scale = FALSE)$x
data.pca <- as.data.frame(data.pca)

summary(data.pca)
print(data.pca)
data.pca$scores

#PCA loadings
data.pca$loadings[, 1:2]

all_individuals2$pc1 <- data.pca$PC1[match(row.names(all_individuals2), row.names(data.pca))]
all_individuals2$pc2 <- data.pca$PC2[match(row.names(all_individuals2), row.names(data.pca))]

fCORT_22_sna$pc1 <- all_individuals2$pc1[match(fCORT_22_sna$JID, all_individuals2$JID)]
fCORT_22_sna$pc2 <- all_individuals2$pc2[match(fCORT_22_sna$JID, all_individuals2$JID)]

plot(fCORT_22_sna$age_2022, fCORT_22_sna$pc1)
plot(fCORT_22_sna$age_2022, fCORT_22_sna$pc2)
plot(fCORT_22_sna$body_cond, fCORT_22_sna$pc1)
plot(fCORT_22_sna$body_cond, fCORT_22_sna$pc2)
plot(fCORT_22_sna$sex_num, fCORT_22_sna$pc1)
plot(fCORT_22_sna$sex_num, fCORT_22_sna$pc2)

fCORT_sna_m0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m1 <- lmer(log(CORT_pg.mg) ~ degree + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m2 <- lmer(log(CORT_pg.mg) ~ strength + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
#fCORT_sna_m3 <- lmer(log(CORT_pg.mg) ~ mean_strength + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
#fCORT_sna_m4 <- lmer(log(CORT_pg.mg) ~ eigenvector + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
#fCORT_sna_m5 <- lmer(log(CORT_pg.mg) ~ betweenness + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m6 <- lmer(log(CORT_pg.mg) ~ social_diff + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
#fCORT_sna_m7 <- lmer(log(CORT_pg.mg) ~ closeness + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m8 <- lmer(log(CORT_pg.mg) ~ transitivity + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m9 <- lmer(log(CORT_pg.mg) ~ visit_number + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m10 <- lmer(log(CORT_pg.mg) ~ body_cond_noNA + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
fCORT_sna_m11 <- lmer(log(CORT_pg.mg) ~ site + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
#fCORT_sna_m12 <- lmer(log(CORT_pg.mg) ~ pc1 + (1|pair_ID_yr_sampled), data = fCORT_22_sna)
#fCORT_sna_m13 <- lmer(log(CORT_pg.mg) ~ pc2 + (1|pair_ID_yr_sampled), data = fCORT_22_sna)

fCORT_sna_m0 <- lm(log(CORT_pg.mg) ~ 1, data = fCORT_22_sna)
#fCORT_sna_m1 <- lm(log(CORT_pg.mg) ~ degree, data = fCORT_22_sna)
fCORT_sna_m2 <- lm(log(CORT_pg.mg) ~ strength + visit_number, data = fCORT_22_sna)
#fCORT_sna_m3 <- lm(log(CORT_pg.mg) ~ mean_strength, data = fCORT_22_sna)
#fCORT_sna_m4 <- lm(log(CORT_pg.mg) ~ eigenvector, data = fCORT_22_sna)
#fCORT_sna_m5 <- lm(log(CORT_pg.mg) ~ betweenness, data = fCORT_22_sna)
fCORT_sna_m6 <- lm(log(CORT_pg.mg) ~ social_diff + + visit_number, data = fCORT_22_sna)
#fCORT_sna_m7 <- lm(log(CORT_pg.mg) ~ closeness, data = fCORT_22_sna)
fCORT_sna_m8 <- lm(log(CORT_pg.mg) ~ transitivity + visit_number, data = fCORT_22_sna)
#fCORT_sna_m9 <- lm(log(CORT_pg.mg) ~ visit_number, data = fCORT_22_sna)
fCORT_sna_m10 <- lm(log(CORT_pg.mg) ~ body_cond_noNA + site, data = fCORT_22_sna)
fCORT_sna_m11 <- lm(log(CORT_pg.mg) ~ site, data = fCORT_22_sna)
#fCORT_sna_m12 <- lm(log(CORT_pg.mg) ~ pc1, data = fCORT_22_sna)
#fCORT_sna_m13 <- lm(log(CORT_pg.mg) ~ pc2, data = fCORT_22_sna)
fCORT_sna_m14 <- lm(log(CORT_pg.mg) ~ group_size, data = fCORT_22_sna)

summary(fCORT_sna_m14)
Anova(fCORT_sna_m14, type = "II")
Anova(fCORT_sna_m14, type = "III")
confint(fCORT_sna_m3)

AICc(fCORT_sna_m0, fCORT_sna_m1, fCORT_sna_m2, fCORT_sna_m3, fCORT_sna_m4, fCORT_sna_m5, fCORT_sna_m6, fCORT_sna_m7, fCORT_sna_m8, fCORT_sna_m9, fCORT_sna_m10, fCORT_sna_m11, fCORT_sna_m12, fCORT_sna_m13)
AICctab(fCORT_sna_m2, fCORT_sna_m6, fCORT_sna_m8, fCORT_sna_m10)

AICc(fCORT_sna_m1, fCORT_sna_m2, fCORT_sna_m3, fCORT_sna_m4, fCORT_sna_m5, fCORT_sna_m6, fCORT_sna_m9, fCORT_sna_m10, fCORT_sna_m11, fCORT_sna_m12)
AICctab(fCORT_sna_m1, fCORT_sna_m2, fCORT_sna_m3, fCORT_sna_m4, fCORT_sna_m5, fCORT_sna_m6, fCORT_sna_m9, fCORT_sna_m10, fCORT_sna_m11, fCORT_sna_m12)

plot(fCORT_22_sna$mean_strength, fCORT_22_sna$CORT_pg.mg)
cor.test(fCORT_22_sna$mean_strength, fCORT_22_sna$CORT_pg.mg)

#Assortment and fCORT?
library(sna)

r = assortment.continuous(amfCORT, V(gfCORT), weighted = TRUE, SE = TRUE)
r

t= 10000
permutation=vector(length=t)
for (i in 1:t){
  s=sample(V(gfCORT), length(V(gfCORT)),replace=F) 
  permutation[i]=assortment.continuous(amfCORT, s)$r 
  }

hist(permutation)
abline(v= r$r, lty=2, col="red", lwd=3)

p =(length(which(permutation > r$r))+1)/(t+1) 
p

#Edge weight (association strength) and similarity in fCORT?

#One observation per dyad
fCORT_sna_m14 <- lmer(log(ind1ind2_CORT_diff_norm) ~ g_edge_weights + (1|ind1) + (1|ind2), data = dyad_edge)

#Two observations per dyad
  
fCORT_sna_m15 <- lmer(log(ind1ind2_CORT_diff_norm) ~ g_edge_weights + (1|ind1) + (1|ind2), data = dyad_edge_CORT_long)

fCORT_sna_m15 <- lmer(log(CORT) ~ associate_CORT * g_edge_weights + body_cond + site + (1|JID) + (1|associate_JID), data = dyad_edge_CORT_long)

#Scaling (z-transforming) numeric variables - can help with model convergence
fCORT_sna_m15 <- lmer(log(CORT) ~ associate_CORT_z * g_edge_weights_z + body_cond_z + site + (1|JID) + (1|associate_JID), data = dyad_edge_CORT_long)

summary(fCORT_sna_m14)
Anova(fCORT_sna_m14, type = "II")
Anova(fCORT_sna_m14, type = "III")
vif(fCORT_sna_m15)
confint(fCORT_sna_m15)

simulationOutput1 <- simulateResiduals(fittedModel = fCORT_sna_m15, plot = F)
plot(simulationOutput1)

testDispersion(simulationOutput1)
testZeroInflation(simulationOutput1)

#6. Parent fCORT and offspring social behaviour ----

all_individuals_juv <- subset(all_individuals, all_individuals$age == 1)

all_individuals_juv$mother_ID <- LH_early$MOTHER.ID[match(all_individuals_juv$JID, LH_early$ID)]
all_individuals_juv$father_ID <- LH_early$FATHER.ID[match(all_individuals_juv$JID, LH_early$ID)]
all_individuals_juv$parent_ID <- paste(all_individuals_juv$mother_ID, all_individuals_juv$father_ID, sep = "_")

all_individuals_juv$mother_fCORT <- fCORT_22$CORT_pg.mg[match(all_individuals_juv$mother_ID, fCORT_22$JID)]
all_individuals_juv$father_fCORT <- fCORT_22$CORT_pg.mg[match(all_individuals_juv$father_ID, fCORT_22$JID)]

all_individuals_juv$dyad_ID_mother1 <- paste(all_individuals_juv$JID, all_individuals_juv$mother_ID, sep = " ")
all_individuals_juv$dyad_ID_mother2 <- paste(all_individuals_juv$mother_ID, all_individuals_juv$JID, sep = " ")

all_individuals_juv$JID <- as.character(all_individuals_juv$JID)
all_individuals_juv$mother_ID <- as.character(all_individuals_juv$mother_ID)

all_individuals_juv$mother_dyad_ID <- ifelse(all_individuals_juv$JID < all_individuals_juv$mother_ID, paste(all_individuals_juv$mother_ID, all_individuals_juv$JID, sep = " "), paste(all_individuals_juv$mother_ID, all_individuals_juv$JID, sep = " "))
all_individuals_juv$mother_edge <- dyad_edge$g_edge_weights[match(all_individuals_juv$mother_dyad_ID, dyad_edge$g_dyads)]

all_individuals_juv$father_dyad_ID <- ifelse(all_individuals_juv$JID < all_individuals_juv$father_ID, paste(all_individuals_juv$father_ID, all_individuals_juv$JID, sep = " "), paste(all_individuals_juv$father_ID, all_individuals_juv$JID, sep = " "))
all_individuals_juv$father_edge <- dyad_edge$g_edge_weights[match(all_individuals_juv$father_dyad_ID, dyad_edge$g_dyads)]

all_individuals_juv <- all_individuals_juv %>% 
  rowwise() %>% 
  mutate(parent_fCORT = mean(c_across(c(mother_fCORT, father_fCORT)), na.rm = TRUE))

all_individuals_juv <- all_individuals_juv %>% 
  rowwise() %>% 
  mutate(parent_edge = mean(c_across(c(mother_edge, father_edge)), na.rm = TRUE))

all_individuals_juv$parent_edge_binary <- ifelse(is.nan(all_individuals_juv$parent_edge), 0, 1)

all_individuals_juv$growth_rate <- LH_early$GROWTH_RATE[match(all_individuals_juv$JID, LH_early$ID)]
all_individuals_juv$body_condition <- LH_early$BODY_CONDITION[match(all_individuals_juv$JID, LH_early$ID)]
all_individuals_juv$brood_size <- LH_early$BROOD_SIZE[match(all_individuals_juv$JID, LH_early$ID)]
all_individuals_juv$hatch_day <- LH_early$HATCH_DAY[match(all_individuals_juv$JID, LH_early$ID)]

all_individuals_juv$mother_visits <- fCORT_22$visit_nr_own_box[match(all_individuals_juv$mother_ID, fCORT_22$JID)]
all_individuals_juv$father_visits <- fCORT_22$visit_nr_own_box[match(all_individuals_juv$father_ID, fCORT_22$JID)]

all_individuals_juv <- all_individuals_juv %>% 
  rowwise() %>% 
  mutate(parent_visits = mean(c_across(c(mother_visits, father_visits)), na.rm = TRUE))

all_individuals_juv$social_diff <- am_data$social_diff[match(all_individuals_juv$JID, am_data$JID)]

all_individuals_juv$closeness <- all_individuals$closeness[match(all_individuals_juv$JID, all_individuals$JID)]
all_individuals_juv$transitivity <- all_individuals$transitivity[match(all_individuals_juv$JID, all_individuals$JID)]
all_individuals_juv$pc1 <- all_individuals2$pc1[match(all_individuals_juv$JID, all_individuals2$JID)]
all_individuals_juv$pc2 <- all_individuals2$pc2[match(all_individuals_juv$JID, all_individuals2$JID)]
all_individuals_juv$betweenness_binary <-ifelse(all_individuals_juv$betweenness == 0, 0, 1)

all_individuals_juv$strength_z <- scale(all_individuals_juv$strength)
all_individuals_juv$visit_number_z <- scale(all_individuals_juv$visit_number)
all_individuals_juv$parent_fCORT_z <- scale(all_individuals_juv$parent_fCORT)

all_individuals_juv <- subset(all_individuals_juv, !is.nan(all_individuals_juv$parent_fCORT))

all_individuals_juv$obs <- seq_len(nrow(all_individuals_juv))

#Model with degree
#Model just with prior
default_prior()

fCORT_juv_brm1_prior <- brm(degree ~ parent_fCORT_z + visit_number_z +(1|parent_ID) + (1|obs),
                              data = all_individuals_juv,
                              family = poisson(link = "log"),
                              prior = c(
                                prior(normal(3.5, 0.5), class = "Intercept"),
                                prior(normal(0, 0.15), class = "b"),
                                prior(exponential(3), class = "sd")),
                              control = list(adapt_delta = 0.95),
                            sample_prior = "only")

#Prior predictive checks 
pp_check(fCORT_juv_brm1_prior, ndraws = 100) +
  scale_x_log10()

#Model 
fCORT_juv_brm1 <- brm(degree ~ parent_fCORT_z + visit_number_z + (1|parent_ID) + (1|obs),
                            data = all_individuals_juv,
                            family = poisson(link = "log"),
                            prior = c(
                              prior(normal(3.5, 0.5), class = "Intercept"),
                              prior(normal(0, 0.15), class = "b"),
                              prior(exponential(3), class = "sd")),
                            control = list(adapt_delta = 0.99)
)

summary(fCORT_juv_brm1)

check_collinearity(fCORT_juv_brm1)

#Posterior predictive checks
pp_check(fCORT_juv_brm1, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_juv_brm1)
launch_shinystan(fCORT_juv_brm1)

#Posterior distribution
as_draws_df(fCORT_juv_brm1)
mcmc_areas(fCORT_juv_brm1)
mcmc_intervals(fCORT_juv_brm1)

#Evaluation and interpretation
loo(fCORT_juv_brm1, moment_match = TRUE)
fitted(fCORT_juv_brm1, scale = "response")
conditional_effects(fCORT_juv_brm1)

bayes_R2(fCORT_juv_brm1)

#Sensitivity analysis and prior checks 
prior_summary(fCORT_juv_brm1)

#Extract predictions
fitted(fCORT_juv_brm1)

hypothesis(fCORT_juv_brm1, "parent_fCORT_z < 0") 

cond <- conditional_effects(fCORT_juv_brm1, effect = "parent_fCORT_z")
p <- plot(cond, points=F)
offspring_degree_plot <- p[[1]] + 
  theme_few(base_size = 14) +
  ylim(0,100) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = parent_fCORT_z, y = degree), 
    data = all_individuals_juv, 
    color = "black",
    size = 2,
    alpha = 0.3,
    shape = 16,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Parent fCORT (pg/mg)", y = "Degree centrality")
offspring_degree_plot


#Model with strength

#Model just with prior
default_prior()

fCORT_juv_brm2_prior <- brm(strength ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                           data = all_individuals_juv, 
                           family = gaussian,
                           prior = c(
                             set_prior("normal(0.5, 0.1)", class = "Intercept"),  
                             set_prior("normal(0, 0.1)", class = "b"),
                             set_prior("exponential(1)", class = "sigma")), 
                           sample_prior = "only"
)

fCORT_juv_brm2_prior <- brm(strength ~ parent_fCORT_z + visit_number_z +(1|parent_ID), 
                      data = all_individuals_juv, 
                      family = lognormal(link = "identity"),
                      prior = c(
                        prior(normal(0.09, 0.75), class = "Intercept"),  
                        prior(normal(0, 0.3), class = "b"),
                        prior(lognormal(log(0.35), 0.3), class = "sigma")),
                      sample_prior = "only"
)

#Prior predictive checks 
pp_check(fCORT_juv_brm2_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_juv_brm2 <- brm(strength ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                            data = all_individuals_juv, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(0.09, 0.75), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),            
)

summary(fCORT_juv_brm2)

check_collinearity(fCORT_juv_brm2)

#Posterior predictive checks
pp_check(fCORT_juv_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_juv_brm2)
launch_shinystan(fCORT_juv_brm2)

#Posterior distribution
as_draws_df(fCORT_juv_brm2)
mcmc_areas(fCORT_juv_brm2)
mcmc_intervals(fCORT_juv_brm2)

#Evaluation and interpretation
loo(fCORT_juv_brm2, moment_match = TRUE)
fitted(fCORT_juv_brm2, scale = "response")
conditional_effects(fCORT_juv_brm2)

bayes_R2(fCORT_juv_brm2)

#Sensitivity analysis and prior checks 
prior_summary(fCORT_juv_brm2)

#Extract predictions
fitted(fCORT_juv_brm2)

hypothesis(fCORT_juv_brm2, "parent_fCORT_z < 0") 

cond <- conditional_effects(fCORT_juv_brm2, effect = "parent_fCORT_z")
p <- plot(cond, points=F)
offspring_strength_plot <- p[[1]] + 
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = parent_fCORT_z, y = strength), 
    data = all_individuals_juv, 
    color = "black",
    size = 2,
    alpha = 0.3,
    shape = 16,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Parent fCORT (pg/mg)", y = "Strength")
offspring_strength_plot

offspring_plot <- ggarrange(offspring_degree_plot, offspring_strength_plot, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),  widths = c(1, 1), font.label = list(
  family = "Garamond",
  face = "plain",   # or "bold", "italic"
  size = 14,
  color = "black"
))
offspring_plot

#Model with eigenvector

#Model just with prior
default_prior()

fCORT_juv_brm3_prior <- brm(eigenvector ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                            data = all_individuals_juv, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(-1.65, 0.75), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),
                            sample_prior = "only"
)

#Prior predictive checks 
pp_check(fCORT_juv_brm3_prior, ndraws = 100) +
  scale_x_log10()

#Model
fCORT_juv_brm3 <- brm(eigenvector ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                      data = all_individuals_juv, 
                      family = lognormal(link = "identity"),
                      prior = c(
                        prior(normal(-1.65, 0.75), class = "Intercept"),  
                        prior(normal(0, 0.3), class = "b"),
                        prior(lognormal(log(0.35), 0.3), class = "sigma")),            
)

summary(fCORT_juv_brm3)

check_collinearity(fCORT_juv_brm3)

#Posterior predictive checks
pp_check(fCORT_juv_brm3, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_juv_brm3)
launch_shinystan(fCORT_juv_brm3)

#Posterior distribution
as_draws_df(fCORT_juv_brm3)
mcmc_areas(fCORT_juv_brm3)
mcmc_intervals(fCORT_juv_brm3)

#Evaluation and interpretation
loo(fCORT_juv_brm3, moment_match = TRUE)
fitted(fCORT_juv_brm3, scale = "response")
conditional_effects(fCORT_juv_brm3)

bayes_R2(fCORT_juv_brm3)

#Sensitivity analysis and prior checks 
prior_summary(fCORT_juv_brm3)

#Extract predictions
fitted(fCORT_juv_brm3)

hypothesis(fCORT_juv_brm3, "parent_fCORT_z < 0") 


#Model with betweenness

#Model just with prior
default_prior()

fCORT_juv_brm4_prior <- brm(betweenness_binary ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                            data = all_individuals_juv, 
                            family = bernoulli(),
                            prior = c(
                              prior(normal(-0.7, 0.7), class = "Intercept"),  
                              prior(normal(0, 0.5), class = "b"),
                              prior(student_t(3, 0, 2.5), class = "sd")),
                            sample_prior = "only"
)

#Prior predictive checks 
pp_check(fCORT_juv_brm4_prior, ndraws = 100) 

fCORT_juv_brm4 <- brm(betweenness_binary ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                            data = all_individuals_juv, 
                            family = bernoulli(),
                            prior = c(
                              prior(normal(-0.7, 0.7), class = "Intercept"),  
                              prior(normal(0, 0.5), class = "b"),
                              prior(student_t(3, 0, 2.5), class = "sd"))
)

summary(fCORT_juv_brm4)

check_collinearity(fCORT_juv_brm4)

#Posterior predictive checks
pp_check(fCORT_juv_brm4, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_juv_brm4)
launch_shinystan(fCORT_juv_brm4)

#Posterior distribution
as_draws_df(fCORT_juv_brm4)
mcmc_areas(fCORT_juv_brm4)
mcmc_intervals(fCORT_juv_brm4)

#Evaluation and interpretation
loo(fCORT_juv_brm4, moment_match = TRUE)
fitted(fCORT_juv_brm4, scale = "response")
conditional_effects(fCORT_juv_brm4)

bayes_R2(fCORT_juv_brm4)

#Sensitivity analysis and prior checks 
prior_summary(fCORT_juv_brm4)

#Extract predictions
fitted(fCORT_juv_brm4)

hypothesis(fCORT_juv_brm4, "parent_fCORT_z < 0") 


#Model with social differentiation

#Model just with prior
default_prior()

fCORT_juv_brm5_prior <- brm(social_diff ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                            data = all_individuals_juv, 
                            family = lognormal(link = "identity"),
                            prior = c(
                              prior(normal(1.65, 0.75), class = "Intercept"),  
                              prior(normal(0, 0.3), class = "b"),
                              prior(lognormal(log(0.35), 0.3), class = "sigma")),
                            sample_prior = "only"
)

#Prior predictive checks 
pp_check(fCORT_juv_brm5_prior, ndraws = 100) +
  scale_x_log10()

fCORT_juv_brm5 <- brm(social_diff ~ parent_fCORT_z + visit_number_z + (1|parent_ID), 
                      data = all_individuals_juv, 
                      family = lognormal(link = "identity"),
                      prior = c(
                        prior(normal(1.65, 0.75), class = "Intercept"),  
                        prior(normal(0, 0.3), class = "b"),
                        prior(lognormal(log(0.35), 0.3), class = "sigma")),
                      control = list(adapt_delta = 0.90)
)

summary(fCORT_juv_brm5)

check_collinearity(fCORT_juv_brm5)

#Posterior predictive checks
pp_check(fCORT_juv_brm5, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_juv_brm5)
launch_shinystan(fCORT_juv_brm5)

#Posterior distribution
as_draws_df(fCORT_juv_brm5)
mcmc_areas(fCORT_juv_brm5)
mcmc_intervals(fCORT_juv_brm5)

#Evaluation and interpretation
loo(fCORT_juv_brm5, moment_match = TRUE)
fitted(fCORT_juv_brm5, scale = "response")
conditional_effects(fCORT_juv_brm5)

bayes_R2(fCORT_juv_brm5)

#Sensitivity analysis and prior checks 
prior_summary(fCORT_juv_brm5)

#Extract predictions
fitted(fCORT_juv_brm5)

hypothesis(fCORT_juv_brm5, "parent_fCORT_z > 0") 


#Model with parent foraging associations

#Model just with prior
default_prior()

fCORT_juv_brm6_prior <- brm(parent_edge_binary ~ parent_fCORT_z + (1|parent_ID), 
                            data = all_individuals_juv, 
                            family = bernoulli(),
                            prior = c(
                              prior(normal(-0.4, 0.5), class = "Intercept"),  
                              prior(normal(0, 0.5), class = "b"),
                              prior(student_t(3, 0, 2.5), class = "sd")),
                            sample_prior = "only"
)

#Prior predictive checks 
pp_check(fCORT_juv_brm6_prior, ndraws = 100) 

fCORT_juv_brm6 <- brm(parent_edge_binary ~ parent_fCORT_z + (1|parent_ID), 
                      data = all_individuals_juv, 
                      family = bernoulli(),
                      prior = c(
                        prior(normal(-0.7, 0.7), class = "Intercept"),  
                        prior(normal(0, 0.5), class = "b"),
                        prior(student_t(3, 0, 2.5), class = "sd"))
)

summary(fCORT_juv_brm6)

check_collinearity(fCORT_juv_brm6)

#Posterior predictive checks
pp_check(fCORT_juv_brm6, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(fCORT_juv_brm6)
launch_shinystan(fCORT_juv_brm6)

#Posterior distribution
as_draws_df(fCORT_juv_brm6)
mcmc_areas(fCORT_juv_brm6)
mcmc_intervals(fCORT_juv_brm6)

#Evaluation and interpretation
loo(fCORT_juv_brm6, moment_match = TRUE)
fitted(fCORT_juv_brm6, scale = "response")
conditional_effects(fCORT_juv_brm6)

bayes_R2(fCORT_juv_brm6)

#Sensitivity analysis and prior checks 
prior_summary(fCORT_juv_brm6)

#Extract predictions
fitted(fCORT_juv_brm6)

hypothesis(fCORT_juv_brm6, "parent_fCORT_z < 0") 


#(04) ARCHIVE ----

#Archive
fCORT_22_F <- subset(fCORT_22, fCORT_22$sex == "F")
feather_CORT_22_F$partner_id <- pairbond$MID[match( feather_CORT_22_F$JID, pairbond$FID)]

feather_CORT_22_M <- subset(feather_CORT_22, feather_CORT_22$sex == "M")
feather_CORT_22_M$partner_id <- pairbond$FID[match( feather_CORT_22_M$JID, pairbond$MID)]

feather_CORT_22 <- rbind(feather_CORT_22_F, feather_CORT_22_M)

feather_CORT_22$partner_fCORT <- feather_CORT_22$CORT_pg.mg[match(feather_CORT_22$JID, feather_CORT_22$partner_id)]

feather_CORT_222 <- feather_CORT_22 %>% pivot_wider(
  id_cols = c("box"), names_from = c("sex"), values_from = c("CORT_pg.mg"))

feather_CORT_222 <- na.omit(feather_CORT_222)
feather_CORT_222 <- subset(feather_CORT_222, !feather_CORT_222$box == "N")

plot(feather_CORT_222$F, feather_CORT_222$M)
cor.test(feather_CORT_222$F, feather_CORT_222$M)

co_reg <- lm(F ~ M, data = feather_CORT_222)
co_reg1 <- lm(F ~ 1, data = feather_CORT_222)
anova(co_reg1, co_reg)

feathers23$partner_id <- LH_pairs$PARTNER.ID[match(feathers23$JID, LH_pairs$ID)]
feathers23$pair_id <- LH_pairs$pair_ID[match(feathers23$JID, LH_pairs$ID)]

#Adding information about site: X, Y, Z
feathers23[substr(feathers23$box,1,1)=="X","site"]<-"X"
feathers23[substr(feathers23$box,1,1)=="Y","site"]<-"Y"
feathers23[substr(feathers23$box,1,1)=="Z","site"]<-"Z"
feathers23[substr(feathers23$box,1,1)=="M","site"]<-"MW"
feathers23[substr(feathers23$box,1,1)=="N","site"]<-"NN"

table(feathers23$site)

#Get pair bond strength data
pairbond <- read.csv("Data/time_budget_summary.csv", header = T, stringsAsFactors = F)
pairbond$bond_strength <- pairbond$percent_of_total_length

pairbond  <- pairbond %>% 
  as_tibble() %>% 
  mutate(pair_id = paste(FID, MID))

pairbond$box <- pairbond$Box

pairbond <- subset(pairbond, pairbond$Behavior == "Both")
pairbond <- pairbond[, c("bond_strength", "pair_id", "box", "FID", "MID")]

#Merge feather and pair bond datasets
feathers23 <- merge(x = feathers23, y = pair_bond, by = "pair_id", all.x = TRUE)

#Write feather data set as CSV
write.csv(feathers23,"feathers23.csv", row.names = FALSE)

#read feather 2023 data
feathers23 <- read.csv("feathers23.csv", header = T, stringsAsFactors = F)

feathers23 <- read.csv("feathers2023.csv", header = T, stringsAsFactors = F)

#Lab sheet
feather_CORT_lab <- read.csv("Data/Feather_CORT_1123.csv", header = T, stringsAsFactors = F)
summary(feather_CORT_lab)

#Fill in sex
LH_sex <- LH[LH$SEX != "",]
feathers23$sex <- LH_sex$SEX[match(feathers23$JID, LH_sex$ID)]
table(feathers23$sex)

#Add female year grown 
fCORT$female_ID <- breeding_summary$fem.ID[match(fCORT$box_grown_year, breeding_summary$box_grown_year)]

#Add male year grown
fCORT$male_ID <- breeding_summary$male.ID[match(fCORT$box_grown_year, breeding_summary$box_grown_year)]

#Add pair_ID
fCORT$pair_ID1 <- paste(fCORT$female_ID, fCORT$male_ID, sep = " ")
fCORT$pair_ID2 <- paste(fCORT$male_ID, fCORT$female_ID, sep = " ")

fCORT$pair_ID2 <- LH_pairs$pair_ID[match(fCORT$JID_year_grown, LH_pairs$JID_year_grown)]

#Add feeder edge weight
fCORT$g_edge_weights1 <- dyad_edge$g_edge_weights[match(fCORT$pair_ID1, dyad_edge$pair_ID1)]
fCORT$g_edge_weights2 <- dyad_edge$g_edge_weights[match(fCORT$pair_ID1, dyad_edge$pair_ID2)]

fCORT$g_edge_weights1[is.na(fCORT$g_edge_weights1)] <- 0
fCORT$g_edge_weights2[is.na(fCORT$g_edge_weights2)] <- 0

fCORT$feeder_edge_weight <- fCORT$g_edge_weights1 + fCORT$g_edge_weights2

#Model archive

pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~ body_cond_binary * pbvid + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT)
pairbond_lmm2 <- lmer(log(CORT_pg.mg) ~ site + pbvid + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~ body_cond_binary * pbvid + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_body)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  scale(pbvid) * site +  scale(weight) + sex + scale(age) + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_test)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  scale(pbvid) * scale(weight) +  site +  scale(age) + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_test)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  scale(weight) + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_test)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  sex  + site +  body_cond_binary + age + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_body)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  body_cond_binary + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_body)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  scale(pbvid) * site + scale(pbvid) * body_cond_binary + scale(pbvid) * sex + scale(age) + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_test)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  scale(pbvid) * sex  + site + scale(pbvid) * body_cond_binary + scale(age) + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_test)
pairbond_lmm1 <- lmer(log(CORT_pg.mg) ~  feather_year_grown + (1|JID) + (1|pair_ID), data = fCORT_test2)

#Individual characteristics 
ind_fCORT_m <- lmer(log(CORT_pg.mg) ~ sex + age + weight + site +  (1|pair_ID), data = fCORT)
ind_fCORT_m <- lmer(log(CORT_pg.mg) ~ bite +  (1|pair_ID), data = fCORT)

summary(ind_fCORT_m)
Anova(ind_fCORT_m)

#Social characteristics
social_fCORT_m  <- lmer(log(CORT_pg.mg) ~ scale(time_together) * scale(brood_weight) + (1|JID), data = fCORT)
social_fCORT_m  <- lmer(log(CORT_pg.mg) ~ degree + (1|JID), data = fCORT)

social_fCORT_m  <- lm(log(CORT_pg.mg) ~ feeder_edge_weight, data = fCORT_22)
social_fCORT_m  <- lm(log(CORT_pg.mg) ~ mean_strength, data = fCORT_22)

summary(social_fCORT_m)
Anova(social_fCORT_m)

#Environmental characteristics

env_fCORT_m <- lmer(log(CORT_pg.mg) ~ site + feather_year_grown + (1|JID), data = fCORT)
summary(env_fCORT_m)
Anova(env_fCORT_m)

AIC(ind_fCORT_m, env_fCORT_m)

#Global model 
fCORT_m <- lmer(log(CORT_pg.mg) ~ site + weight + (1|JID) + (1|feather_year_grown), data = fCORT)
summary(fCORT_m)
Anova(fCORT_m)


fCORT_int_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid * age + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_pbvid)

fCORT_int_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ pbvid * site + age + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_pbvid)

fCORT_int_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ pbvid + body_cond_binary + age + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_pbvid)
fCORT_int_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ pbvid * body_cond_binary + age + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_pbvid)

fCORT_int_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid + sex + age + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_pbvid)
fCORT_int_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid * sex + age + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID), data = fCORT_pbvid)

fCORT_int_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid + feather_year_grown + age + site + (1|JID) + (1|pair_ID), data = fCORT_pbvid)
fCORT_int_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid * feather_year_grown + age + site + (1|JID) + (1|pair_ID), data = fCORT_pbvid)

#(3) Full model pair bond (video) and covariates, without interactions
fCORT_full_pbvid_lmm <- lmer(log(CORT_pg.mg) ~ pbvid + site + body_cond (1|feather_year_grown), data = fCORT)

summary(fCORT_full_pbvid_lmm)
Anova(fCORT_full_pbvid_lmm)
anova(pairbond_lmm1, pairbond_lmm)

vif(lmer(log(CORT_pg.mg) ~ pbvid + site + sex + age + body_cond + (1|feather_year_grown), data = fCORT_body), type = "predictor")

simulationOutput1 <- simulateResiduals(fittedModel = fCORT_full_pbvid_lmm, plot = F)
plot(simulationOutput1)

plot(fCORT$age, fCORT$body_cond)
plot(fCORT$age, fCORT$feather_year_grown)

cor.test(fCORT$age, fCORT$body_cond)

fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond  + sex + age + feather_year_grown + (1|JID), data = fCORT)
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|feather_year_grown) + (1|JID), data = fCORT)
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ age + (1|feather_year_grown) + (1|JID), data = fCORT)

plot(fCORT_full_cov_brm)
conditional_effects(fCORT_full_cov_brm, effects = "body_cond")

emmeans(fCORT_full_cov_lmm, list(pairwise ~ feather_year_grown), adjust = "bonferroni")

fCORT_bite_cov_lmm <- lmer(log(CORT_pg.mg) ~ bite + (1|JID), data = fCORT)

fCORT$feather_year_grown <- as.numeric(fCORT$feather_year_grown)

coreg_m2 <- lmer(log(CORT_F) ~ CORT_M + body_cond + (1|pair_ID), data = fCORT_coreg)

coreg_m1 <- lm(log(CORT_F) ~ CORT_M + body_cond, data = fCORT_coreg)
coreg_m2 <- lm(log(CORT_F) ~ CORT_M, data = fCORT_coreg)

plot(fCORT_coreg$CORT_M, fCORT_coreg$CORT_F)
cor.test(fCORT_coreg$CORT_F, fCORT_coreg$CORT_M)
plot(fCORT_coreg2$CORT_F, fCORT_coreg2$CORT_M)
cor.test(fCORT_coreg2$CORT_F, fCORT_coreg2$CORT_M)
plot(fCORT_coreg_test$CORT_F, fCORT_coreg_test$CORT_M)

plot(fCORT$pbvid, fCORT$CORT_pg.mg)
plot(fCORT_22$pbfeedgmm, fCORT_22$CORT_pg.mg)
plot(fCORT_22$pbprospgmm, fCORT_22$CORT_pg.mg)
plot(fCORT_22$pbowngmm, fCORT_22$CORT_pg.mg)
plot(fCORT_22$degree, fCORT_22$CORT_pg.mg)
plot(fCORT_22$strength, fCORT_22$CORT_pg.mg)
plot(fCORT_22$mean_strength, fCORT_22$CORT_pg.mg)
plot(fCORT_22$eigenvector, fCORT_22$CORT_pg.mg)
plot(fCORT_22$betweenness, fCORT_22$CORT_pg.mg)

boxplot(fCORT$pbvid ~  fCORT$site)
boxplot(fCORT$pbvid ~  fCORT$age)
boxplot(fCORT$pbvid ~  fCORT$body_cond_bindary)

p <- sjp.lmer(coreg_m1, type = "fe.prob", vars = "CORT_M", 
              show.scatter = TRUE, geom.colors = "Black", 
              facet.grid = FALSE, show.ci = TRUE)

fCORT_coreg$predF = predict(coreg_m1, newdata = fCORT_coreg, re.form = NA)
coreg_ci = model.matrix(formula(coreg_m1)[-2], fCORT_coreg)

df <- ggpredict(coreg_m1, terms = c("CORT_M"))

ggplot(fCORT_coreg, aes(x = CORT_M, y = CORT_F)) + 
  geom_line(color= "black") +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"))

ggPredict(coreg_m1,se=TRUE)

effects_coreg <- effect(term= "CORT_M", mod= coreg_m1)
summary(effects_coreg) 
effects_coreg <- as.data.frame(exp(effects_coreg))

coreg_plot <- ggplot() + 
  geom_point(data= fCORT_coreg, aes(CORT_M, CORT_F)) + 
  geom_line(data= effects_coreg, aes(x= CORT_M, y=fit), color="black") +
  geom_ribbon(data= effects_coreg, aes(x= CORT_M, ymin=lower, ymax=upper), alpha= 0.3, fill="grey") +
  labs(x="Male fCORT (pg/mg)", y="Feale fCORT (pg/mg)")

coreg_plot

fCORT_pbvid <- subset(fCORT_pbvid, !fCORT_pbvid$site == "X")
fCORT_pbvid <- subset(fCORT_pbvid, !fCORT_pbvid$site == "Y")
fCORT_pbvid <- subset(fCORT_pbvid, !fCORT_pbvid$sex == "F")

fCORT_pbvid_glmm0 <- glmmTMB(log(CORT_pg.mg) ~ body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm1 <- glmmTMB(log(CORT_pg.mg) ~ pbvid + body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm2 <- glmmTMB(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm3 <- glmmTMB(log(CORT_pg.mg) ~ pbvid * site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm4 <- glmmTMB(log(CORT_pg.mg) ~ pbvid * body_cond + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm5 <- glmmTMB(log(CORT_pg.mg) ~ pbvid * age + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm6 <- glmmTMB(log(CORT_pg.mg) ~ pbvid * sex + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))
fCORT_pbvid_glmm7 <- glmmTMB(log(CORT_pg.mg) ~ site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))

fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid * body_cond + pbvid * site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ pbvid + body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ pbvid * body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ pbvid * site + body_cond + (1|feather_year_grown)+ (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid * sex  + body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid * feather_year_grown + body_cond + site + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid * age + body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)

fCORT_pbvid_lmm <- lmer(log(CORT_pg.mg) ~ pbvid * sex * body_cond + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_lmm <- lmer(log(CORT_pg.mg) ~ pbvid * body_cond_binary_num + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)
fCORT_pbvid_glmm1 <- glmmTMB(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid, family = Gamma(link = "log"))
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid)

fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_glmm2 <- glmmTMB(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA, family = Gamma(link = "log"))

anova(fCORT_pbvid_lmm0, fCORT_int_pbvid_lmm2)

AICctab(fCORT_pbvid_glmm0, fCORT_pbvid_glmm1, fCORT_pbvid_glmm2, fCORT_pbvid_glmm3, fCORT_pbvid_glmm4, fCORT_pbvid_glmm5, fCORT_pbvid_glmm6, fCORT_pbvid_glmm7)
AICctab(fCORT_pbvid_lmm2, fCORT_pbvid_lmm3, fCORT_pbvid_lmm5, fCORT_pbvid_lmm6)
AICctab(fCORT_pbvid_lmm2, fCORT_pbvid_glmm2)

plot(fCORT_22$visit_number, fCORT_22$CORT_pg.mg)
plot(fCORT_22$visit_duration, fCORT_22$CORT_pg.mg)
plot(fCORT_22$degree, fCORT_22$CORT_pg.mg)
plot(fCORT_22$eigenvector, fCORT_22$CORT_pg.mg)
cor.test(fCORT_22$eigenvector, fCORT_22$CORT_pg.mg)
plot(fCORT_22$betweenness, fCORT_22$CORT_pg.mg)
cor.test(fCORT_22$betweenness, fCORT_22$CORT_pg.mg)
plot(fCORT_22$strength, fCORT_22$CORT_pg.mg)
plot(fCORT_22$strength, fCORT_22$visit_number)
plot(fCORT_22_sna$social_diff, fCORT_22_sna$CORT_pg.mg)

cor.test(fCORT_22$strength, fCORT_22$CORT_pg.mg)
plot(fCORT_22$mean_strength, fCORT_22$CORT_pg.mg)
cor.test(fCORT_22$mean_strength, fCORT_22$CORT_pg.mg)
plot(fCORT_22$feeder_edge_weight, fCORT_22$CORT_pg.mg)
plot(fCORT_22$feeder_edge_weight, fCORT_22$time_together)
plot(fCORT_22.2$feeder_edge_weight, fCORT_22.2$CORT_pg.mg)
cor.test(fCORT_22.2$feeder_edge_weight, fCORT_22.2$CORT_pg.mg)
plot(fCORT_22$visit_nr_own_box, fCORT_22$CORT_pg.mg)
plot(fCORT_22$visit_nr_prospect, fCORT_22$CORT_pg.mg)
plot(fCORT_22$pbfeedgmm, fCORT_22$CORT_pg.mg)
cor.test(fCORT_22$pbfeedgmm, log(fCORT_22$CORT_pg.mg))
plot(fCORT_22$pbprospgmm, log(fCORT_22$CORT_pg.mg))
plot(fCORT_22$pbowngmm, log(fCORT_22$CORT_pg.mg))

parent_fCORT_m4 <- lmer(eigenvector ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
parent_fCORT_m5 <- lmer(betweenness ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)

fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ body_cond + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid + body_cond + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid + body_cond + site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid * site + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm9 <- lmer(log(CORT_pg.mg) ~ pbvid * body_cond + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm10 <- lmer(log(CORT_pg.mg) ~ pbvid * age + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm11 <- lmer(log(CORT_pg.mg) ~ pbvid * sex + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA)

fCORT_pbvid_lmm0.1 <- lmer(log(CORT_pg.mg) ~ pbvid * site + pbvid * age + pbvid * sex + pbvid * body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA, na.action = "na.fail")
fCORT_pbvid_lmm0.2 <- dredge(fCORT_pbvid_lmm0.1)

model.avg(fCORT_pbvid_lmm0.2)
confset.95p <- get.models(fCORT_pbvid_lmm0.2, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)

fCORT_pbvid_glmm1 <- glmmTMB(CORT_pg.mg ~ pbvid_mean + (1|feather_year_grown) + (1|JID), data = fCORT, family = Gamma(link = "log"))

fCORT_pbvid_lmm12 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * box_density + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm13 <- lmer(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA_Y)
fCORT_pbvid_lmm14 <- lmer(log(CORT_pg.mg) ~ pbvid + (1|feather_year_grown) + (1|JID), data = fCORT_pbvid_noNA_Z)

fCORT_sna_m0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m1 <- lmer(log(CORT_pg.mg) ~ scale(degree) + scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m2 <- lmer(log(CORT_pg.mg) ~ scale(strength) + scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m3 <- lmer(log(CORT_pg.mg) ~ scale(mean_strength) + scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m4 <- lmer(log(CORT_pg.mg) ~ scale(eigenvector) + scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m5 <- lmer(log(CORT_pg.mg) ~ scale(betweenness) + scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m6 <- lmer(log(CORT_pg.mg) ~ scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m7 <- lmer(log(CORT_pg.mg) ~ scale(social_diff) + scale(visit_number) + (1|pair_ID), data = fCORT_22_sna)

fCORT_sna_m0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m1 <- lmer(log(CORT_pg.mg) ~ degree + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m2 <- lmer(log(CORT_pg.mg) ~ strength + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m3 <- lmer(log(CORT_pg.mg) ~ mean_strength + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m4 <- lmer(log(CORT_pg.mg) ~ eigenvector + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m5 <- lmer(log(CORT_pg.mg) ~ betweenness + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m6 <- lmer(log(CORT_pg.mg) ~ visit_number + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m7 <- lmer(log(CORT_pg.mg) ~ social_diff + (1|pair_ID), data = fCORT_22_sna)

fCORT_sna_m7 <- lmer(log(CORT_pg.mg) ~ pbfeedgmm + (1|pair_ID), data = fCORT_22_sna)
fCORT_sna_m8 <- lmer(log(CORT_pg.mg) ~ pbprospgmm + (1|pair_ID), data = fCORT_22_sna)

fCORT_m1 <- lmer(log(CORT_pg.mg) ~ strength + site, data = fCORT_22)
fCORT_m1 <- lmer(log(CORT_pg.mg) ~ scale(mean_strength) * scale(body_cond) + (1|pair_ID), data = fCORT_22)
fCORT_m1 <- lmer(log(CORT_pg.mg) ~ strength * weight + visit_number + (1|pair_ID), data = fCORT_22)
fCORT_m2 <- lmer(log(CORT_pg.mg) ~ strength + weight + visit_number + (1|pair_ID), data = fCORT_22)
fCORT_m3 <- lmer(log(CORT_pg.mg) ~ weight + visit_number + (1|pair_ID), data = fCORT_22)
fCORT_m4 <- lmer(log(CORT_pg.mg) ~ scale(pbfeedgmm) * scale(weight) + scale(visit_number) + (1|pair_ID), data = fCORT_22)
fCORT_m5 <- lmer(log(CORT_pg.mg) ~ feeder + (1|pair_ID), data = fCORT_22)

pairbond_glmm <- lmer(log(CORT_pg.mg) ~ pbprospgmm * site + (1|pair_ID), data = fCORT22_test)

pairbond_glmm <- lmer(log(CORT_pg.mg) ~ pbfeedgmm * weight + (1|pair_ID), data = fCORT22_test)
pairbond_glmm1 <- lmer(log(CORT_pg.mg) ~ weight + (1|pair_ID), data = fCORT22_test)

pairbond_glmm1 <- lmer(log(CORT_pg.mg) ~  site + (1|pair_ID), data = fCORT22_test)

interact_plot(pairbond_glmm, pred = pbprospgmm, modx = site)
interact_plot(pairbond_glmm, pred = pbprospgmm, modx = weight)

interact_plot(pairbond_glmm, pred = pbfeedgmm, modx = site)
interact_plot(pairbond_glmm, pred = pbfeedgmm, modx = weight)

summary(fCORT_m1)
Anova(fCORT_m1)
anova(fCORT_m3, fCORT_m4)
interact_plot(fCORT_m4 , pred = pbfeedgmm, modx = weight)

#Social characteristics pb_22_long_pair
plot(pb22_long_pair$pb_vid, pb22_long_pair$fCORT)
plot(pb22_long_pair$pbfeedgmm, pb22_long_pair$fCORT)
plot(pb22_long_pair$pbprospgmm, pb22_long_pair$fCORT)
cor.test(pb22_long_pair$pbprospgmm, pb22_long_pair$fCORT)
plot(pb22_long_pair$pbowngmm, pb22_long_pair$fCORT)
plot(pb22_long_pair$degree, pb22_long_pair$fCORT)
plot(pb22_long_pair$strength, pb22_long_pair$fCORT)
plot(pb22_long_pair$mean_strength, pb22_long_pair$fCORT)
plot(pb22_long_pair$mean_strength_pb, pb22_long_pair$fCORT)
plot(pb22_long_pair$mean_strength_pb, pb22_long_pair$pbfeedgmm)
cor.test(pb22_long_pair$mean_strength_pb, pb22_long_pair$pbfeedgmm)

social_fCORT_m1 <- lmer(log(fCORT) ~  scale(visit_number) + scale(strength) * scale(pbfeedgmm)  + (1|pair_ID), data = pb22_long_pair)
social_fCORT_m1 <- lmer(log(fCORT) ~  scale(visit_number) + scale(strength) + (1|pair_ID), data = pb22_long_pair)

social_fCORT_m1 <- lmer(log(fCORT) ~  scale(visit_number) + scale(mean_strength_pb) * scale(pbfeedgmm)  + (1|pair_ID), data = pb22_long_pair)

social_fCORT_m2 <- lmer(log(fCORT) ~  scale(visit_number) + scale(pbfeedgmm) + (1|pair_ID), data = pb22_long_pair)

social_fCORT_m1 <- lmer(log(fCORT) ~  pbprospgmm * pbfeedgmm + (1|pair_ID), data = pb22_long_pair)

summary(social_fCORT_m1)
Anova(social_fCORT_m1)
anova(social_fCORT_m1, social_fCORT_m2)

interact_plot(social_fCORT_m, pred = pbfeedgmm, modx = strength)
interact_plot(social_fCORT_m, pred = strength, modx = pbfeedgmm)

cor.test(pb22_long_pair$pbfeedgmm, pb22_long_pair$strength)

fCORT_coreg_m1 <- lmer(log(CORT_F) ~ CORT_M  + (1|pair_ID_year), data = fCORT_coreg)

fCORT_pbvid_lmm1 <- glmmTMB(log(CORT_pg.mg) ~ pbvid_mean + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA, family = Gamma)

fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond + sex + age + (1|feather_year_grown) + (1|JID) + (1|box), data = fCORT)

fCORT_full_cov_lmm <- with(fCORT_mi, lmer(log(CORT_pg.mg) ~ site + body_cond + age + sex + pbvid_mean + (1|feather_year_grown) + (1|JID)))
fCORT_full_cov_lmm_pooled <- pool(fCORT_full_cov_lmm)
summary(fCORT_full_cov_lmm_pooled)

m <- 10
fCORT_mi <- mice(fCORT, m = m, print = FALSE)

fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + scale(body_cond) + sex + scale(age) + + scale(pbvid_mean) + (1|feather_year_grown) + (1|JID), data = fCORT)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ site + body_cond + age + sex + (1|feather_year_grown) + (1|JID), data = fCORT)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ site + body_cond + age + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm2 <- lmer(log(CORT_pg.mg) ~ site + body_cond + age + sex + pbvid_mean + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm3 <- lmer(log(CORT_pg.mg) ~ site * pbvid_mean + body_cond + age + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm4 <- lmer(log(CORT_pg.mg) ~ body_cond * pbvid_mean + site + age + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm5 <- lmer(log(CORT_pg.mg) ~ sex * pbvid_mean + site + body_cond + age + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm6 <- lmer(log(CORT_pg.mg) ~ age * pbvid_mean + body_cond + site + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm7 <- lmer(log(CORT_pg.mg) ~ site + body_cond + age + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm8 <- lmer(log(CORT_pg.mg) ~ site + body_cond + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm9 <- lmer(log(CORT_pg.mg) ~ site + age + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))
fCORT_full_cov_lmm10 <- lmer(log(CORT_pg.mg) ~ body_cond + age + sex + (1|feather_year_grown) + (1|JID), data = subset(fCORT, !is.nan(pbvid_mean)))

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * (site + body_cond + sex + age) + (1|feather_year_grown) + (1|JID), data = fCORT, REML = FALSE)
fCORT_full_cov_lmm2 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + site + body_cond + sex + age + (1|feather_year_grown) + (1|JID), data = fCORT, REML = FALSE)
anova(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2)

model_list <- list(fCORT_pbvid_lmm1, fCORT_pbvid_lmm2, fCORT_pbvid_lmm3, fCORT_pbvid_lmm4, fCORT_pbvid_lmm5, fCORT_pbvid_lmm6, fCORT_pbvid_lmm7, fCORT_pbvid_lmm8,  fCORT_pbvid_lmm9,  fCORT_pbvid_lmm10,  fCORT_pbvid_lmm11)

model_sel <- model.sel(model_list)
print(model_sel)

avg_model <- model.avg(model_sel, subset = delta < 2)
summary(avg_model)

coef_table <- summary(avg_model)$coefmat.full
print(coef_table)
Anova(avg_model)

fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ pbvid + (1|JID), data = fCORT_pbvid)

fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ poly(pbvid_mean,2) + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA)

fCORT_full_cov_brm <- brm(
  bf(CORT_pg.mg ~ mi(pbvid_mean) + site + body_cond + age + sex + (1|JID) + (1|feather_year_grown)) +
    bf(pbvid_mean | mi() ~ 1),
  data = fCORT,
  chains = 4, cores = 4
)

fCORT_full_cov_brm <- brm_multiple(CORT_pg.mg ~ pbvid_mean + site + body_cond + age + sex + (1|JID) + (1|feather_year_grown),
                                   data = fCORT_mi,
                                   chains = 4, cores = 4)

summary(fCORT_full_cov_brm)
loo(fCORT_full_cov_brm)

fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA)

fCORT_pbvid <- subset(fCORT, !is.na(fCORT$pbvid_mean))
fCORT_pbvid <- subset(fCORT, !is.nan(fCORT$pbvid_mean))

fCORT_pbvid_noNA <- subset(fCORT_pbvid, !is.na(fCORT_pbvid$body_cond))
#fCORT_pbvid_noNA <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$pbvid < 0.75)
fCORT_pbvid_noNA_Y <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$site == "Y")
fCORT_pbvid_noNA_Z <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$site == "Z")

#NA removed, N = 78
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + (1|JID) + (1|pair_ID_yr_sampled) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm9 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm10 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * age + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm11 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * sex + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm12 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * box_density + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)

#NA removed, one outlier removed, N = 77
fCORT_pbvid_noNA_nool <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$CORT_pg.mg < 7.817757)
#fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ poly(pbvid_mean,2) + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm9 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm10 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * age + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm11 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * sex + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm12 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * box_density + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)

#NA not removed, N = 89
#fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ site + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ body_cond + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + site + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + site + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * site + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm9 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * body_cond + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm10 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * age + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm11 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * sex + (1|JID) + (1|pair_ID), data = fCORT)
fCORT_pbvid_lmm12 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * box_density + (1|JID) + (1|pair_ID), data = fCORT)

summary(fCORT_pbvid_lmm11)
Anova(fCORT_pbvid_lmm11, type = "II")
Anova(fCORT_pbvid_lmm11, type = "III")

summary(fCORT_pbvid_lmm8, dispersion = 1) 

AICc(fCORT_pbvid_lmm1, fCORT_pbvid_lmm2, fCORT_pbvid_lmm3, fCORT_pbvid_lmm4, fCORT_pbvid_lmm5, fCORT_pbvid_lmm6, fCORT_pbvid_lmm7, fCORT_pbvid_lmm8,  fCORT_pbvid_lmm9,  fCORT_pbvid_lmm10,  fCORT_pbvid_lmm11)
AICctab(fCORT_pbvid_lmm1, fCORT_pbvid_lmm2, fCORT_pbvid_lmm3, fCORT_pbvid_lmm4, fCORT_pbvid_lmm5, fCORT_pbvid_lmm6, fCORT_pbvid_lmm7, fCORT_pbvid_lmm8,  fCORT_pbvid_lmm9,  fCORT_pbvid_lmm10,  fCORT_pbvid_lmm11)

vif(fCORT_pbvid_lmm1)

simulationOutput1 <- simulateResiduals(fittedModel = fCORT_pbvid_lmm2, plot = F)
plot(simulationOutput1)

interact_plot(fCORT_pbvid_lmm, pred = pbvid, modx = sex, mod2 = body_cond)

#NA removed, N = 78
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + (1|JID) + (1|pair_ID_yr_sampled) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * site + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm9 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * body_cond + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm10 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * age + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm11 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * sex + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm12 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * box_density + (1|JID) + (1|pair_ID) + (1|feather_year_grown), data = fCORT_pbvid_noNA)

#fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ poly(pbvid_mean,2) + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ site + body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean + body_cond + site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm8 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * site + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm9 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * body_cond + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm10 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * age + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm11 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * sex + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)
fCORT_pbvid_lmm12 <- lmer(log(CORT_pg.mg) ~ pbvid_mean * box_density + (1|JID) + (1|pair_ID), data = fCORT_pbvid_noNA_nool)

#fCORT_pbvid_lmm0 <- lmer(log(CORT_pg.mg) ~ 1 + (1|JID) + (1|pair_ID), data = fCORT)

AICc(fCORT_pbvid_lmm1, fCORT_pbvid_lmm2, fCORT_pbvid_lmm3, fCORT_pbvid_lmm4, fCORT_pbvid_lmm5, fCORT_pbvid_lmm6, fCORT_pbvid_lmm7, fCORT_pbvid_lmm8,  fCORT_pbvid_lmm9,  fCORT_pbvid_lmm10,  fCORT_pbvid_lmm11)
AICctab(fCORT_pbvid_lmm1, fCORT_pbvid_lmm2, fCORT_pbvid_lmm3, fCORT_pbvid_lmm4, fCORT_pbvid_lmm5, fCORT_pbvid_lmm6, fCORT_pbvid_lmm7, fCORT_pbvid_lmm8,  fCORT_pbvid_lmm9,  fCORT_pbvid_lmm10,  fCORT_pbvid_lmm11)

summary(fCORT_pbvid_lmm8, dispersion = 1) 

fCORT_coreg$pbvid <- fCORT$pbvid[match(fCORT_coreg$pair_ID, fCORT$pair_ID)]
fCORT_coreg$pbfeedgmm <- fCORT_22$pbfeedgmm[match(fCORT_coreg$pair_ID, fCORT_22$pair_ID)]
fCORT_coreg$pbprospgmm <- fCORT_22$pbprospgmm[match(fCORT_coreg$pair_ID, fCORT_22$pair_ID)]

fCORT_coreg$pair_ID1 <- fCORT_coreg$pair_ID
fCORT_coreg$pair_ID2 <- fCORT_coreg$pair_ID

fCORT_coreg$JID <- str_split(fCORT_coreg$pair_ID1, " ", simplify = TRUE)
fCORT_coreg$JID <- fCORT_coreg$JID[,1]

fCORT_coreg$body_cond <- fCORT$body_cond[match(fCORT_coreg$JID, fCORT$JID)]

interact_plot(fCORT_coreg_m1, pred = partner_CORT, modx = site)

site_plot <- ggplot(fCORT, aes(x = site, y = CORT_pg.mg, fill = site)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 8)) +
  labs(x="Study site", y = "fCORT (pg/mg)") +
  theme_classic(base_size = 14)+
  geom_jitter(alpha = 0.5, shape=16, position=position_jitter(0.1), size = 2) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "black")  
site_plot

set_theme(base = theme_classic(base_size = 14), theme.font = "Garamond")
body_plot <- plot_model(fCORT_pbvid_lmm1, type="pred", title = "", colors= "bw", show.data = FALSE, terms = "body_cond_z", jitter = 0.001) + scale_y_continuous(limits = c(0, 9))
body_plot

body_plot2 <- body_plot + ylab("fCORT (pg/mg)") + xlab("Body condition (residuals)")
body_plot2

body_plot3 <- body_plot2 + geom_point(alpha = 0.5, size = 2, shape = 16, data = fCORT_pbvid, aes(x=body_cond_z, y=CORT_pg.mg, colour = "black")) + theme(text=element_text(size=14))
body_plot3

# Prior: weakly informative for 2x2 residual covariance
prior <- list(R = list(V = diag(2), nu = 3), G = list())

prior <- list(
  R = list(V = diag(2), nu = 3),
  G = list(G1 = list(V = 1, nu = 1))   # univariate random intercept for pair_ID
)

# Model: adjust covariates to your variable names
model <- MCMCglmm(
  cbind(female_CORT_log, male_CORT_log) ~ trait - 1 +
    trait:(factor(site) + female_body_cond_z + male_body_cond_z),
  random = ~ pair_ID,        # univariate random intercept, NOT us(trait):pair_ID
  rcov = ~ us(trait):units,
  family = c("gaussian","gaussian"),
  data = fCORT_coreg2,
  prior = prior,
  nitt = 260000, burnin = 60000, thin = 200
)

model <- MCMCglmm(
  cbind(female_CORT_log, male_CORT_log) ~ trait - 1,
  random = ~ pair_ID,        # univariate random intercept, NOT us(trait):pair_ID
  rcov = ~ us(trait):units,
  family = c("gaussian","gaussian"),
  data = fCORT_coreg2,
  prior = prior,
  nitt = 260000, burnin = 60000, thin = 200
)

summary(model)

# Extract posterior correlation
VCV <- model$VCV
cov_fm <- VCV[ , "traitfemale_CORT_log:traitmale_CORT_log.units"]
var_f  <- VCV[ , "traitfemale_CORT_log:traitfemale_CORT_log.units"]
var_m  <- VCV[ , "traitmale_CORT_log:traitmale_CORT_log.units"]

cor_post <- cov_fm / sqrt(var_f * var_m)

# posterior mean and 95% credible interval
mean(cor_post)
HPDinterval(as.mcmc(cor_post))

# approximate pMCMC (two-sided)
mean(abs(cor_post) >= abs(mean(cor_post)))  # not the ideal pMCMC; better to use fraction of posterior >0
mean(cor_post > 0)   # posterior probability correlation > 0

prior <- list(
  R = list(V = diag(2), nu = 3),
  G = list(
    G1 = list(V = 1, nu = 1)   # random intercept per pair
  )
)

model <- MCMCglmm(
  cbind(female_CORT, male_CORT) ~ trait - 1,
  random = ~ pair_ID,             # univariate random intercept
  rcov = ~ us(trait):units,       # still estimates the partner covariance
  family = c("gaussian", "gaussian"),
  data = fCORT_coreg2,
  prior = prior,
  nitt = 130000,
  burnin = 30000,
  thin = 100
)

VCV <- model$VCV

cov_fm <- VCV[ , "traitfemale_CORT:traitmale_CORT.units"]
var_f  <- VCV[ , "traitfemale_CORT:traitfemale_CORT.units"]
var_m  <- VCV[ , "traitmale_CORT:traitmale_CORT.units"]

cor_fm <- cov_fm / sqrt(var_f * var_m)

mean(cor_fm)
HPDinterval(as.mcmc(cor_fm))

summary(model)

fCORT_coreg_long <- as.data.frame(fCORT_coreg_long)
fCORT_coreg_long$logCORT1 <- log(fCORT_coreg_long$CORT)
fCORT_coreg_long$logCORT2 <- log(fCORT_coreg_long$CORT)

model <- MCMCglmm(
  cbind(logCORT1, logCORT2) ~ trait - 1,      # trait-specific intercepts only
  random = ~ us(trait):pair_ID,               # dyad-level covariance
  rcov = ~ us(trait):units,                   # residual covariance
  family = c("gaussian", "gaussian"),
  data = fCORT_coreg_long,
  prior = prior,
  nitt = 130000,
  burnin = 30000,
  thin = 100
)

fCORT_coreg_long2 <- subset(fCORT_coreg_long, !fCORT_coreg_long$site == "X")

model_site <- MCMCglmm(
  cbind(logCORT1, logCORT2) ~ trait - 1 + trait:site + trait:body_cond,
  random = ~ us(trait):pair_ID,
  rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"),
  data = fCORT_coreg_long2,
  prior = list(
    R = list(V = diag(2) * 0.5, nu = 2),  # slightly stronger prior
    G = list(G1 = list(V = diag(2) * 0.5, nu = 2))
  ),
  nitt = 130000,
  burnin = 30000,
  thin = 100,
  singular.ok = TRUE
)


library(MCMCglmm)

model <- MCMCglmm(
  cbind(log(female_CORT), log(male_CORT)) ~ trait - 1,     # site can have separate effects per trait
  random = ~ us(trait):pair_ID,
  rcov = ~ us(trait):units,
  family = c("gaussian", "gaussian"),
  data = fCORT_coreg2,
  prior = prior,
  nitt = 130000,
  burnin = 30000,
  thin = 100,
  singular.ok = TRUE
)

summary(model)

summary(model)
plot(model$VCV)        # check mixing
plot(model$Sol)

# dyad-level posterior correlation:
var_f <- model$VCV[,"traitfemale_CORT:traitfemale_CORT.pair_ID"]
var_m <- model$VCV[,"traitmale_CORT:traitmale_CORT.pair_ID"]
cov_fm <- model$VCV[,"traitfemale_CORT:traitmale_CORT.pair_ID"]
cor_dyad <- cov_fm / sqrt(var_f * var_m)
summary(cor_dyad)
plot(density(cor_dyad)); abline(v=0, col="red", lty=2)


model_site$VCV[, c("traitlogCORT1:traitlogCORT1.pair_ID",
                   "traitlogCORT2:traitlogCORT2.pair_ID",
                   "traitlogCORT1:traitlogCORT2.pair_ID")]

var1 <- model_site$VCV[,"traitlogCORT1:traitlogCORT1.pair_ID"]
var2 <- model_site$VCV[,"traitlogCORT2:traitlogCORT2.pair_ID"]
cov12 <- model_site$VCV[,"traitlogCORT1:traitlogCORT2.pair_ID"]

cor_dyad <- cov12 / sqrt(var1 * var2)
summary(cor_dyad)

plot(model_site$VCV)       # dyad-level covariance traces
plot(model_site$Sol)

install.packages("MCMCglmm")

library(MCMCglmm)


prior <- list(
  R = list(V = diag(2), nu = 3),           # residual 2x2 covariance
  G = list(
    G1 = list(V = diag(2), nu = 3)         # pair_ID 2x2 covariance
  )
)


model <- MCMCglmm(
  cbind(female_CORT_log, male_CORT_log) ~ trait - 1 +
    trait:(factor(site) + female_body_cond_z + male_body_cond_z),
  
  random = ~ us(trait):pair_ID,      # multivariate random intercept for pair_ID
  rcov   = ~ us(trait):units,        # unstructured residual covariance
  family = c("gaussian", "gaussian"),
  data   = fCORT_coreg2,
  prior  = prior,
  
  nitt   = 520000,   # total iterations
  burnin = 120000,   # burn-in
  thin   = 100       # thinning interval
)

prior <- list(
  R = list(V = diag(2), nu = 3)  # weakly informative prior for residual covariance
)

model_simple <- MCMCglmm(
  cbind(female_CORT_log, male_CORT_log) ~ trait - 1 +
    trait:(factor(site) + female_body_cond_z + male_body_cond_z),
  
  rcov = ~ us(trait):units,       # unstructured residual covariance
  family = c("gaussian", "gaussian"),
  data = fCORT_coreg2,
  prior = prior,
  
  nitt = 520000, 
  burnin = 120000, 
  thin = 100
)

model_simple <- MCMCglmm(
  cbind(female_CORT_log, male_CORT_log) ~ trait - 1,
  rcov = ~ us(trait):units,       # unstructured residual covariance
  family = c("gaussian", "gaussian"),
  data = fCORT_coreg2,
  prior = prior,
  
  nitt = 520000, 
  burnin = 120000, 
  thin = 100
)

summary(model)
summary(model_simple)

#Plots ----
colours <- c("#FFFFFF", "#CCCCCC")
colours <- c("#CCCCCC", "#333333")
colours <- c("#FFFFFF", "#CCCCCC", "#999999", "#333333")
colours <- c("#FFFFFF", "#CCCCCC", "#999999", "#333333", "#FFFFFF", "#CCCCCC", "#999999", "#333333", "#FFFFFF")

#Plot 
plot(fCORT$feather_weight, fCORT$extraction_weight_g)

fCORT_hist <-ggplot(fCORT, aes(x= CORT_pg.mg)) + 
  geom_histogram(color="white", fill="black") +
  theme_bw(base_size = 14) +
  labs(x = "fCORT (pg/mg)", y = "Frequency")
fCORT_hist

fCORT_hist <-ggplot(fCORT, aes(x= CORT_pg.mg)) + 
  geom_histogram(color="white", fill="black") +
  theme_base(base_size = 14) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond")) +
  labs(x = "fCORT (pg/mg)", y = "Frequency")
fCORT_hist

#Individual 
plot(fCORT$bite, fCORT$CORT_pg.mg)
plot(fCORT$age, fCORT$CORT_pg.mg)
plot(fCORT$tarsus, fCORT$CORT_pg.mg)
plot(fCORT$weight, fCORT$CORT_pg.mg)
plot(fCORT_body$body_cond, fCORT_body$CORT_pg.mg)

boxplot(fCORT$CORT_pg.mg ~ fCORT$box_owner_grown)
boxplot(fCORT_body$CORT_pg.mg ~ fCORT_body$body_cond_binary)

#Body condition continuous for manuscript
pred <- ggpredict(fCORT_full_cov_lmm, terms = "body_cond_z", type = "fixed")

body_plot <- ggplot(pred, aes(x, predicted)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey70", alpha = 0.3) +
  geom_point(data = fCORT,
             aes(x = body_cond_z, y = CORT_pg.mg),
             alpha = 0.3, size = 2, shape = 16) +
  labs(x = "Body condition (residuals)",
       y = "fCORT (pg/mg)") +
  scale_y_continuous(limits = c(0, 9)) +
  theme_few(base_size = 14) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond")) 
body_plot

body_plot <- ggplot(fCORT_body, aes(x= body_cond, y= CORT_pg.mg)) +
  geom_point(size=3, shape=20) + geom_rug() +
  labs(x="Body condition", y = "fCORT (pg/mg)", size =14) +
  theme_bw(base_size = 14)
body_plot

set_theme(base = theme_few(base_size = 14), theme.font = "Garamond")
body_plot <- plot_model(fCORT_full_cov_lmm, type="pred", title = "", colors= "bw", show.data = FALSE, terms = "body_cond_z", jitter = 0.001) + scale_y_continuous(limits = c(0, 9))
body_plot

body_plot2 <- body_plot + ylab("fCORT (pg/mg)") + xlab("Body condition (residuals)")
body_plot2

body_plot3 <- body_plot2 + geom_point(alpha = 0.3, size = 2, shape = 16, data = fCORT, aes(x=body_cond_z, y=CORT_pg.mg, colour = "black")) + theme(text=element_text(size=14))
body_plot3


#Body condition binary
body_plot <- ggplot(fCORT_body, aes(x = body_cond_binary, y = CORT_pg.mg, fill = body_cond_binary)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Body condition", y = "fCORT (pg/mg)") +
  theme_bw(base_size = 14)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "tomato")  
body_plot

#Sex plot for manuscript
colours <- c("#FFFFFF", "#CCCCCC")

sex_plot <- ggplot(fCORT, aes(x = sex, y = CORT_pg.mg, fill = sex)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(x="Sex", y = "fCORT (pg/mg)") +
  theme_classic(base_size = 14)+
  geom_jitter(alpha = 0.5, shape=16, position=position_jitter(0.1), size = 2) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "black")  
sex_plot

sex_plot <- ggplot(fCORT, aes(x = sex, y = CORT_pg.mg, fill = sex)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Sex", y = "fCORT (pg/mg)") +
  theme_few(base_size = 14)+
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 2) +
  stat_summary(fun=mean, geom="point", shape=8, size=2, col = "coral", stroke = 1.25) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
sex_plot

site_plot <- ggplot(fCORT, aes(x = site, y = CORT_pg.mg, fill = site)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Study site", y = "fCORT (pg/mg)") +
  theme_classic(base_size = 14)+
  geom_jitter(alpha = 0.5, shape=16, position=position_jitter(0.1), size = 2) +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral") +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
site_plot

#Age
age_plot <- ggplot(fCORT, aes(x= age, y= CORT_pg.mg)) +
  geom_point(alpha = 0.5, size=2, shape=16) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(x="Age (Years)", y = "fCORT (pg/mg)", size =14) +
  theme_classic(base_size = 14)
age_plot

age_plot <- ggplot(fCORT, aes(x= age, y= CORT_pg.mg)) +
  geom_point(alpha = 0.3, size=2, shape=16) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Age (years)", y = "fCORT (pg/mg)", size =14) +
  theme_few(base_size = 14) + 
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
age_plot

set_theme(base = theme_classic(base_size = 14))
age_plot <- plot_model(fCORT_full_cov_lmm, title = "", type="pred", colors= "bw", show.data = FALSE, terms = "age", jitter = 0.001) + scale_y_continuous(limits = c(0, 10))
age_plot

age_plot2 <- age_plot + ylab("fCORT (pg/mg)") + xlab("Age (years)")
age_plot2

age_plot3 <- age_plot2 + geom_point(alpha = 0.5, size = 2, shape = 16, data= fCORT, aes(x=age, y=CORT_pg.mg, colour = "black")) + theme(text=element_text(size=14))
age_plot3

#Social 

#Pair bond 
#see for fCORT_pbvid below

library(reformulas)

set_theme(base = theme_classic(base_size = 14), theme.font = "Garamond")
pairbond_plot <- plot_model(fCORT_pbvid_lmm1, type="pred", title = "", colors= "bw", show.data = FALSE, terms = "pbvid_mean_z") 
pairbond_plot

pairbond_plot2 <- pairbond_plot + ylab("fCORT (pg/mg)") + xlab("Pair-bond strength")
pairbond_plot2

pairbond_plot3 <- pairbond_plot2 + geom_point(alpha = 0.1, size = 2, shape = 16, data = fCORT_pbvid, aes(x=pbvid_mean_z, y= CORT_pg.mg, colour = "black")) + theme(text=element_text(size=14))
pairbond_plot3

pred <- ggpredict(fCORT_pbvid_lmm1, terms = "pbvid_mean_z", condition = c(site = c("Y","Z")), type = "fixed")

pairbond_plot <- ggplot(pred, aes(x, predicted)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "grey70", alpha = 0.3) +
  geom_point(data = fCORT_pbvid,
             aes(x = pbvid_mean_z, y = CORT_pg.mg),
             alpha = 0.3, size = 2) +
  labs(x = "Pair-bond strength",
       y = "fCORT (pg/mg)") +
  scale_y_continuous(limits = c(0, 10)) +
  theme_few(base_size = 14) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond")) 
pairbond_plot

#Pair bond plot for manuscript 
pairbond_plot <- ggplot(fCORT_pbvid, aes(x= pbvid_mean, y = CORT_pg.mg)) +
  geom_point(alpha = 0.3, size=2, shape= 16) +
  theme_few(base_size = 14) +
  xlab("Pair-bond strength (proportion of time spent together)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,9) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0,2.5,5, 7.5)) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
pairbond_plot

pairbond_plot <- ggplot(fCORT, aes(x= pbvid_mean, y = CORT_pg.mg)) +
  geom_point(alpha = 0.5, size=2, shape= 16) +
  theme_classic(base_size = 14) +
  xlab("Pair-bond strength (proportion of time spent together)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,10) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2,4,6,8,10, 12)) 
pairbond_plot

ggplot(fCORT_pbvid, aes(x= pbvid, y = CORT_pg.mg, col = sex)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Pair-bond strength (proportion of time spent together)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,8) 

ggplot(fCORT_test, aes(x= pbvid, y = CORT_pg.mg, col = body_cond_binary)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Pair-bond strength (time spent together, %)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,100) +
  ylim(0,10) 

ggplot(fCORT_pbvid, aes(x= pbvid, y = CORT_pg.mg, col = site)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Pair-bond strength (time spent together, %)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,10) 

ggplot(fCORT_pbvid, aes(x= pbvid, y = CORT_pg.mg, col = feather_year_grown)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Pair-bond strength (time spent together, %)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,10) 

ggplot(fCORT_test2, aes(x= pbvid, y = CORT_pg.mg, col = body_cond_binary)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Pair-bond strength (time spent together, %)") +
  ylab("fCORT (pg/mg)") +
  xlim(0,50) +
  ylim(0,10) 

#Coregulation plot for manuscript
coreg_plot <- ggplot(fCORT_coreg, aes(x= male_CORT, y = female_CORT)) +
  geom_point(alpha = 0.3, size=2, shape= 16) +
  theme_few(base_size = 14) +
  xlab("Male fCORT (pg/mg)") +
  ylab("Female fCORT (pg/mg)") +
  xlim(0,9) +
  ylim(0,9) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0,2.5,5, 7.5)) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
coreg_plot

ggplot(fCORT_coreg, aes(x= CORT_M, y= CORT_F, col = body_cond)) +
  geom_point(size=4, shape=20) +
  geom_smooth(method = lm, se = TRUE, col = "black") +
  theme_bw(base_size = 14) +
  xlab("Male fCORT (pg/mg)") +
  ylab("Female fCORT (pg/mg)") +
  xlim(0,8) +
  ylim(0,8)

ggplot(fCORT_coreg, aes(x= CORT_M, y= CORT_F)) +
  geom_point(size=4, shape=20) +
  geom_smooth(method = lm, se = TRUE) +
  theme_bw(base_size = 14) +
  xlab("Male fCORT (pg/mg)") +
  ylab("Female fCORT (pg/mg)") +
  xlim(0,8) +
  ylim(0,8)

coreg_plot <- ggplot(fCORT_coreg, aes(x= CORT_M, y= CORT_F, col = year_same_diff)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Male fCORT (pg/mg)") +
  ylab("Female fCORT (pg/mg)") +
  labs(fill = "x") +
  xlim(0,8) +
  ylim(0,8) 

coreg_plot
coreg_plot + guides(col=guide_legend(title="Year"))

plot_model(coreg_m1)
set_theme(base = theme_bw())

p <- plot_model(coreg_m1, title = "", axis.lim = c(0,10), type = "pred", colors= "bw", terms = "CORT_M", transform = "exp", show.data = TRUE)

#Pair bond plot
pairbond_plot <- plot_model(fCORT_pbvid_glmm1, title = "", type = "pred", colors= "bw", terms = "pbvid_mean_c", transform = "exp", show.data = FALSE) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2,4,6,8,10, 12)) 
pairbond_plot 

pairbond_plot2 <- pairbond_plot + ylab("Pair-bond strength") + xlab("fCORT (pg/mg)") 
pairbond_plot2

pairbond_plot3 <- pairbond_plot2 + theme_classic(base_size = 14) + geom_point(alpha = 0.5, size = 2, shape = 16, data= fCORT, aes(x= pbvid, y= CORT_pg.mg)) + theme(text=element_text(size=14))
pairbond_plot3

#Co-regulation plot for manuscript
coreg_plot <- plot_model(fCORT_coreg_m1, title = "", type = "pred", colors= "bw", terms = "partner_CORT", transform = "exp", show.data = FALSE) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0,2.5,5,7.5)) +
  scale_x_continuous(limits = c(0, 9), breaks = c(0,2.5,5,7.5))
coreg_plot

coreg_plot2 <- coreg_plot + ylab("fCORT (pg/mg)") + xlab("Partner fCORT (pg/mg)") 
coreg_plot2

coreg_plot3 <- coreg_plot2 + theme_classic(base_size = 14) + geom_point(alpha = 0.5, size = 2, shape = 16, data= fCORT_coreg_long, aes(x=partner_CORT, y=CORT)) + theme(text=element_text(size=14, family = "Garamond"))
coreg_plot3

coreg_plot <- plot_model(fCORT_coreg_m1, title = "", type = "pred", colors= "bw", terms = "CORT_M", transform = "exp", show.data = FALSE) +
  scale_y_continuous(limits = c(0, 9), breaks = c(0,2.5,5,7.5)) +
  scale_x_continuous(limits = c(0, 9), breaks = c(0,2.5,5,7.5))
coreg_plot

coreg_plot2 <- coreg_plot + ylab("Female fCORT (pg/mg)") + xlab("Male fCORT (pg/mg)") 
coreg_plot2

coreg_plot3 <- coreg_plot2 + theme_classic(base_size = 14) + geom_point(alpha = 0.5, size = 2, shape = 16, data= fCORT_coreg, aes(x=CORT_M, y=CORT_F)) + theme(text=element_text(size=14, family = "Garamond"))
coreg_plot3

#Social network 
#Foraging associates plot for manuscript
assoc_plot <- plot_model(fCORT_sna_m14, title = "", type = "pred", colors= "bw", terms = "associate_CORT", transform = "exp", show.data = FALSE) +
  scale_y_continuous(limits = c(0, 8), breaks = c(0,2.5,5,7.5)) +
  scale_x_continuous(limits = c(0, 8), breaks = c(0,2.5,5,7.5))
assoc_plot

assoc_plot2 <- assoc_plot + ylab("fCORT (pg/mg)") + xlab("Associate fCORT (pg/mg)") 
assoc_plot2

assoc_plot3 <- assoc_plot2 + theme_classic(base_size = 14) + geom_point(alpha = 0.5, size = 2, shape = 16, data= fCORT_coreg_long, aes(x=partner_CORT, y=CORT)) + theme(text=element_text(size=14, family = "Garamond"))
assoc_plot3

#Mean strength
mean_strength_plot <- plot_model(fCORT_sna_m3, title = "", type = "pred", colors= "bw", terms = "mean_strength", transform = "exp", show.data = FALSE) +
  scale_y_continuous(limits = c(0, 10), breaks = c(0,2,4,6,8,10, 12)) +
  scale_x_continuous(limits = c(0, 0.05), breaks = c(0.01, 0.02, 0.03, 0.04, 0.05)) 
mean_strength_plot

mean_strength_plot2 <- mean_strength_plot + ylab("fCORT (pg/mg)") + xlab("Mean strength") 
mean_strength_plot2

mean_strength_plot3 <- mean_strength_plot2 + theme_classic(base_size = 14) + geom_point(alpha = 0.5, size = 2, shape = 16, data= fCORT_22_sna, aes(x= mean_strength, y= CORT_pg.mg)) + theme(text=element_text(size=14))
mean_strength_plot3

ggplot(fCORT_22_sna, aes(x= pc1, y = CORT_pg.mg, col = body_cond)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Pair-bond strength (proportion of time spent together)") +
  ylab("fCORT (pg/mg)") 

ggplot(fCORT_22, aes(x= strength, y = CORT_pg.mg, col = weight)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Weighted degree") +
  ylab("fCORT (pg/mg)") +
  xlim(0,5) +
  ylim(0,8) 

ggplot(fCORT_22, aes(x= pbfeedgmm, y = CORT_pg.mg, col = weight)) +
  geom_point(size=4, shape=20) +
  theme_bw(base_size = 14) +
  xlab("Weighted degree") +
  ylab("fCORT (pg/mg)") +
  xlim(0,1) +
  ylim(0,8) 

#Parental effects on offspring
parent_plot <- plot_model(parent_fCORT_m2, title = "", type = "pred", colors= "bw", terms = "parent_fCORT", show.data = FALSE) 
parent_plot

parent_plot2 <- parent_plot + ylab("Juvenile weighted degree") + xlab("Parent fCORT (pg/mg)") 
parent_plot2

parent_plot3 <- parent_plot2 + theme_classic(base_size = 14) + geom_point(alpha = 0.5, size = 2, shape = 16, data= all_individuals_juv, aes(x= parent_fCORT, y= strength)) + theme(text=element_text(size=14, family = "Garamond"))
parent_plot3

#Box density 
box_density_plot <- ggplot(fCORT, aes(x = box_density, y = CORT_pg.mg, col = site)) + 
  labs(x="Box density", y = "fCORT (pg/mg)") +
  theme_bw(base_size = 14)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5) + theme(legend.position="bottom") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "tomato")  
box_density_plot

box_density_plot <- ggplot(fCORT, aes(x = box_density_cat, y = CORT_pg.mg, col = site)) + 
  labs(x="Box density", y = "fCORT (pg/mg)") +
  geom_violin() + 
  theme_bw(base_size = 14)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5) + theme(legend.position="bottom") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "tomato")  
box_density_plot

plot(fCORT$box_density, fCORT$CORT_pg.mg)
plot(fCORT_Z$box_density, fCORT_Z$CORT_pg.mg)

cor.test(fCORT$box_density, fCORT$CORT_pg.mg)

#Ecological
boxplot(fCORT$CORT_pg.mg ~ fCORT$site)
boxplot(fCORT$CORT_pg.mg ~ fCORT$feather_year_grown)

#Reproductive 
plot(fCORT$clutch_size, fCORT$CORT_pg.mg)
plot(fCORT$fledge_nr, fCORT$CORT_pg.mg)
plot(fCORT$fledge_weight, fCORT$CORT_pg.mg)

colours <- c("#CCCCCC", "#333333")
colours <- c("#FFFFFF", "#CCCCCC", "#999999", "#333333")
colours <- c("#FFFFFF", "#CCCCCC", "#999999", "#666666")
colours <- c("#FFFFFF", "#CCCCCC", "#999999", "#333333", "#FFFFFF", "#CCCCCC", "#999999", "#333333", "#FFFFFF")


ggplot(fCORT22_test, aes(x= pbprospgmm, y= CORT_pg.mg, col = site)) +
  geom_point(size=3, shape=20) + geom_rug() +
  labs(x="Degree", y = "Feeding duration (s)", size =14) +
  theme_bw(base_size = 18)

ggplot(fCORT22_test, aes(x= pbprospgmm, y= CORT_pg.mg, col = weight)) +
  geom_point(size=3, shape=20) + geom_rug() +
  labs(x="Degree", y = "Feeding duration (s)", size =14) +
  theme_bw(base_size = 18)

ggplot(fCORT22_test, aes(x= pbfeedgmm, y= CORT_pg.mg, col = weight)) +
  geom_point(size=3, shape=20) + geom_rug() +
  labs(x="Degree", y = "Feeding duration (s)", size =14) +
  theme_bw(base_size = 18)

#Feeder use
feeder_plot <- ggplot(fCORT_22, aes(x = pbfeeder, y = CORT_pg.mg, fill = pbfeeder)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Feeder use", y = "fCORT (pg/mg)") +
  theme_bw(base_size = 14)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "tomato")  
feeder_plot

#Prospecting
prospect_plot <- ggplot(fCORT_22, aes(x = pbprospect, y = CORT_pg.mg, fill = pbprospect)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Prospecting", y = "fCORT (pg/mg)") +
  theme_bw(base_size = 14)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "tomato")  
prospect_plot

#Site plot for manuscript
#colours <- c("#FFFFFF", "#CCCCCC", "#999999", "#666666")
colours <- c("grey100", "grey90", "grey80", "grey70")

site_plot <- ggplot(fCORT, aes(x = site, y = CORT_pg.mg, fill = site)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.1, fill='white', color="black", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(limits = c(0, 9)) +
  labs(x="Study site", y = "fCORT (pg/mg)") +
  theme_few(base_size = 14)+
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 2) +
  stat_summary(fun=mean, geom="point", shape=8, size=2, col = "coral", stroke = 1.25) +
  theme(legend.position = "none", text = element_text(size = 14, family = "Garamond"))
site_plot

#Years
fCORT$feather_year_grown <- as.factor(fCORT$feather_year_grown)

year_plot <- ggplot(fCORT, aes(x = feather_year_grown, y = CORT_pg.mg, fill = feather_year_grown)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Year", y = "fCORT (pg/mg)") +
  theme_bw(base_size = 14)+
  geom_jitter(shape=16, position=position_jitter(0.2), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "tomato")  
year_plot

year_plot <- ggplot(data = fCORT, aes(x = feather_year_grown, y = CORT_pg.mg, group = JID)) + 
  geom_line(color = "grey")+
  geom_point(size=3, shape=20) +
  theme_bw(base_size = 14) + 
  labs(x="Year", y = "fCORT (pg/mg)")
year_plot

library(patchwork)
site_plot + sex_plot + body_plot + age_plot

env_indiv_fCORT_plot <- ggarrange(site_plot, sex_plot, body_plot, age_plot, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"),  widths = c(1, 1), font.label = list(
  family = "Garamond",
  face = "plain",   # or "bold", "italic"
  size = 14,
  color = "black"
))

env_indiv_fCORT_plot

pair_fCORT_plot <- ggarrange(pairbond_plot, coreg_plot, ncol = 2, nrow = 1, labels = c("(a)", "(b)"),  font.label = list(family = "Garamond"), heights = c(1, 1), widths = c(1, 1))
pair_fCORT_plot


#Full model of environmental and individual covariates ----

#Full model with N = 153 (N = 13 missing for body cond)

#Are independent variables associated with each other? 
body_cond_lmm <- lmer(body_cond ~ age + sex + site + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT)
summary(body_cond_lmm)
Anova(body_cond_lmm)

#Older birds and males have relatively better body condition

#scaled variables
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + scale(body_cond) + sex + scale(age) + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT)
summary(fCORT_full_cov_lmm)
Anova(fCORT_full_cov_lmm)

#Centred variables
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond_c + sex + age_c + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT)
summary(fCORT_full_cov_lmm)
Anova(fCORT_full_cov_lmm)

#Non-scaled and non-centred variables
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond + sex + age + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT)
summary(fCORT_full_cov_lmm)
Anova(fCORT_full_cov_lmm)
confint(fCORT_full_cov_lmm)
vif(fCORT_full_cov_lmm)

#Full model with N = 166 (either remaining 13 datapoints for body condition added or body condition removed)
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond_z + sex + age_z + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT)

summary(fCORT_full_cov_lmm)
Anova(fCORT_full_cov_lmm)

confint(fCORT_full_cov_lmm)
vif(fCORT_full_cov_lmm)

simulationOutput1 <- simulateResiduals(fittedModel = fCORT_full_cov_lmm, plot = F)
plot(simulationOutput1)

emmeans(fCORT_full_cov_lmm, list(pairwise ~ site), adjust = "Tukey")

#subset without 2 outliers (more than 3 SD)
fCORT_full_cov_lmm <- lmer(log(CORT_pg.mg) ~ site + body_cond_noNA + sex + age + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_nool)
summary(fCORT_full_cov_lmm)
Anova(fCORT_full_cov_lmm)

plot_model(fCORT_full_cov_lmm)

interact_plot(fCORT_full_cov_lmm1, pred = box_density, modx = site, data = fCORT_YZ)

#Level of disturbance and box density (just an exploration)
fCORT_YZ_density_noNA <- subset(fCORT_YZ, !is.na(fCORT_YZ$box_density))

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ site * (box_density * disturbance) + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_YZ)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ site * box_density + disturbance + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_YZ)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ box_density * disturbance + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_YZ)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ box_density + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_YZ)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ disturbance + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_YZ)

fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ site + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_YZ)

interaction.plot(fCORT_YZ$box_density, fCORT_YZ$site, fCORT_YZ$CORT_pg.mg, fun = mean)

#Stithians
fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ box_density + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Y)
fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ box_density_cat + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Y)

fCORT_full_cov_lmm2 <- lmer(log(CORT_pg.mg) ~  disturbance + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Y)

fCORT_full_cov_lmm2 <- lmer(log(CORT_pg.mg) ~  box_density + disturbance + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Y)

#Pencoose
fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ box_density + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Z)
fCORT_full_cov_lmm1 <- lmer(log(CORT_pg.mg) ~ box_density_cat + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Z)

fCORT_full_cov_lmm2 <- lmer(log(CORT_pg.mg) ~  disturbance + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Z)

fCORT_full_cov_lmm2 <- lmer(log(CORT_pg.mg) ~  box_density + disturbance + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_Z)

summary(fCORT_full_cov_lmm1)
Anova(fCORT_full_cov_lmm1)
confint(fCORT_full_cov_lmm1)
vif(fCORT_full_cov_lmm1)

summary(fCORT_full_cov_lmm2)
Anova(fCORT_full_cov_lmm2)
confint(fCORT_full_cov_lmm2)
vif(fCORT_full_cov_lmm2)
emmeans(fCORT_full_cov_lmm2, list(pairwise ~ disturbance), adjust = "Tukey")

plot(fCORT$age, fCORT$body_cond)
plot(fCORT$sex_num, fCORT$body_cond)
boxplot(fCORT_Y$CORT_pg.mg ~ fCORT_Y$box_density_cat)

AICc(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2)
AICctab(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2)

AICc(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2, fCORT_full_cov_lmm3, fCORT_full_cov_lmm4, fCORT_full_cov_lmm5, fCORT_full_cov_lmm6)
AICctab(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2, fCORT_full_cov_lmm3, fCORT_full_cov_lmm4, fCORT_full_cov_lmm5, fCORT_full_cov_lmm6)

AICc(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2, fCORT_full_cov_lmm3, fCORT_full_cov_lmm4, fCORT_full_cov_lmm5, fCORT_full_cov_lmm6, fCORT_full_cov_lmm7, fCORT_full_cov_lmm8, fCORT_full_cov_lmm9, fCORT_full_cov_lmm10)
AICctab(fCORT_full_cov_lmm1, fCORT_full_cov_lmm2, fCORT_full_cov_lmm3, fCORT_full_cov_lmm4, fCORT_full_cov_lmm5, fCORT_full_cov_lmm6, fCORT_full_cov_lmm7, fCORT_full_cov_lmm8, fCORT_full_cov_lmm9, fCORT_full_cov_lmm10)



#Pair bonds ----

#fCORT_pbvid_noNA <- subset(fCORT_pbvid, !is.na(fCORT_pbvid$body_cond))
#fCORT_pbvid_noNA <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$pbvid < 0.75)
fCORT_pbvid_noNA_Y <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$site == "Y")
fCORT_pbvid_noNA_Z <- subset(fCORT_pbvid_noNA, fCORT_pbvid_noNA$site == "Z")

#N = 97
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z + site + body_cond_z + (1|JID) +  (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)
fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z + site + body_cond_z + (1|JID) +  (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)

fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * site + body_cond_z + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * body_cond_z + site + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)
#fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * age_z + (1|JID) + site + body_cond_z + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)
#fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * sex + site + body_cond_z + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)
#fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * feather_year_grown_z + site + body_cond_z + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)
#fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * box_density + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid)

#Gamma GLMM? does not seem to be a better fit than LMM
fCORT_pbvid_glmm1 <- glmmTMB(CORT_pg.mg ~ pbvid_mean_z + body_cond_z + site + (1|feather_year_grown) + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_pbvid, family = Gamma(link = "log"))

#One outlier removed, N = 95
fCORT_pbvid_nool <- subset(fCORT_pbvid, fCORT_pbvid$CORT_pg.mg < 7.817757)
#fCORT_pbvid_nool <- subset(fCORT_pbvid, fCORT_pbvid$pbvid_mean_c < 0.5)

fCORT_pbvid_lmm1 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z + site + body_cond_z + (1|JID) +  (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)
fCORT_pbvid_lmm2 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * site + body_cond_z + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)
fCORT_pbvid_lmm3 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * body_cond_z + site + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)
fCORT_pbvid_lmm4 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * age_z + (1|JID) + site + body_cond_z + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)
fCORT_pbvid_lmm5 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * sex + site + body_cond_z + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)
fCORT_pbvid_lmm6 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * feather_year_grown_z + site + body_cond_z + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)
#fCORT_pbvid_lmm7 <- lmer(log(CORT_pg.mg) ~ pbvid_mean_z * box_density + (1|JID) + (1|feather_year_grown) + (1|pair_ID_yr_sampled), data = fCORT_pbvid_nool)

summary(fCORT_pbvid_lmm1)
Anova(fCORT_pbvid_lmm1, type = "II")
Anova(fCORT_pbvid_lmm1, type = "III")

summary(fCORT_pbvid_lmm2)
Anova(fCORT_pbvid_lmm2, type = "II")
Anova(fCORT_pbvid_lmm2, type = "III")

summary(fCORT_pbvid_lmm3)
Anova(fCORT_pbvid_lmm3, type = "II")
Anova(fCORT_pbvid_lmm3, type = "III")

AICc(fCORT_pbvid_lmm1, fCORT_pbvid_glmm1)

vif(fCORT_pbvid_lmm1)
confint(fCORT_pbvid_lmm1)

simulationOutput1 <- simulateResiduals(fittedModel = fCORT_pbvid_lmm1, plot = F)
plot(simulationOutput1)

interact_plot(fCORT_pbvid_lmm2, pred = pbvid, modx = site)
interact_plot(fCORT_pbvid_lmm, pred = pbvid, modx = sex)

#Coregulation ----
#Coregulation? With or without outliers/subset

cor.test(fCORT_coreg$male_CORT_resid, fCORT_coreg$female_CORT_resid)
cor.test(fCORT_coreg2$male_CORT_resid, fCORT_coreg2$female_CORT_resid)

fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + (1|pair_ID), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + (1|pair_ID), data = fCORT_coreg2)

fCORT_coreg_m1 <- lm(log(female_CORT) ~ male_CORT, data = fCORT_coreg_Y)
fCORT_coreg_m1 <- lm(log(female_CORT) ~ male_CORT, data = fCORT_coreg_Z)
fCORT_coreg_m1 <- lm(log(male_CORT) ~ female_CORT, data = fCORT_coreg_Y)
fCORT_coreg_m1 <- lm(log(male_CORT) ~ female_CORT, data = fCORT_coreg_Z)

fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + female_body_cond + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + female_body_cond + (1|pair_ID) + (1|year), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + site + male_body_cond + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + site + male_body_cond + (1|pair_ID) + (1|year), data = fCORT_coreg2)

fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + (1|pair_ID), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + (1|pair_ID), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + site + (1|pair_ID), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + site + (1|pair_ID), data = fCORT_coreg2)

fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + (1|pair_ID), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + (1|pair_ID), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + site + (1|pair_ID), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + site + (1|pair_ID), data = fCORT_coreg2)

summary(fCORT_coreg_m1)
Anova(fCORT_coreg_m1, type = "II")
Anova(fCORT_coreg_m1, type = "III")
vif(fCORT_coreg_m1)
confint(fCORT_coreg_m1)

fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + site + (1|pair_ID) + (1|year), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + (1|pair_ID) + (1|year), data = fCORT_coreg_Y)
fCORT_coreg_m1 <- lmer(log(female_CORT) ~ male_CORT + (1|pair_ID) + (1|year), data = fCORT_coreg_Z)
fCORT_coreg_m1 <- lmer(female_CORT_resid ~ male_CORT_resid + female_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(female_CORT_resid ~ male_CORT_resid + female_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg2)

fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + male_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + male_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + male_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + male_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + male_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + male_body_cond + site + (1|pair_ID) + (1|year), data = fCORT_coreg2)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + (1|pair_ID) + (1|year), data = fCORT_coreg_Y)
fCORT_coreg_m1 <- lmer(log(male_CORT) ~ female_CORT + (1|pair_ID) + (1|year), data = fCORT_coreg_Z)
fCORT_coreg_m1 <- lmer(male_CORT_resid ~ female_CORT_resid + (1|pair_ID) + (1|year), data = fCORT_coreg)
fCORT_coreg_m1 <- lmer(male_CORT_resid ~ female_CORT_resid + (1|pair_ID) + (1|year), data = fCORT_coreg2)

#Does difference in fCORT depend on years or time spent together?
fCORT_coreg_m1 <- lm(log(CORT_diff) ~ years_together, data = fCORT_coreg)
fCORT_coreg_m1 <- lm(log(CORT_F) ~ CORT_M * years_together, data = fCORT_coreg)
fCORT_coreg_m1 <- lm(log(CORT_F) ~ CORT_M * pbvid_mean, data = fCORT_coreg)

summary(fCORT_coreg_m1)
Anova(fCORT_coreg_m1, type = "II")
Anova(fCORT_coreg_m1, type = "III")
vif(fCORT_coreg_m1)
confint(fCORT_coreg_m1)

#fCORT_coreg_long: each indivdiual per dyad once as response and once as predictor
fCORT_coreg_long2 <- subset(fCORT_coreg_long2, fCORT_coreg_long2$CORT < 7.817757) #remove values > 3 SD
fCORT_coreg_long2 <- subset(fCORT_coreg_long2, fCORT_coreg_long2$partner_CORT < 7.817757) #remove values > 3 SD

fCORT_coreg_long2 <- subset(fCORT_coreg_long, fCORT_coreg_long$year_diff == 0)
fCORT_coreg_long2 <- subset(fCORT_coreg_long2, fCORT_coreg_long2$year_diff == 0)
fCORT_coreg_long2 <- subset(fCORT_coreg_long, !fCORT_coreg_long$site == "X")
fCORT_coreg_long2 <- subset(fCORT_coreg_long2, !fCORT_coreg_long2$site == "X")

#decompose within and between dyad variation
fCORT_coreg_long2 <- fCORT_coreg_long %>%
  group_by(pair_ID) %>%
  mutate(partner_between = mean(partner_CORT),
         partner_within  = partner_CORT - partner_between) %>%
  ungroup()

fCORT_coreg_long2$partner_between_z <- scale(fCORT_coreg_long2$partner_between)
fCORT_coreg_long2$partner_within_z <- scale(fCORT_coreg_long2$partner_within)
fCORT_coreg_long2$body_cond_z <- scale(fCORT_coreg_long2$body_cond)

fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT + body_cond + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long)

#With study site
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT * site + body_cond + (1|pair_ID) +  (1|year_sampled), data = fCORT_coreg_long2)
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT * site + body_cond +  (1|year_sampled), data = fCORT_coreg_long)
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT + body_cond + (1|pair_ID) +  (1|year_sampled), data = fCORT_coreg_long2)
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT + body_cond +  (1|year_sampled), data = fCORT_coreg_long2)
summary(fCORT_coreg_m2)

#With years together
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT * years_together + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long)

#Stithians
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT + body_cond + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long_Y)

#Pencoose
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT + body_cond + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long_Z)

#within and between dyad variation 
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_between_z + site + partner_within_z + body_cond_z + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long2)
fCORT_coreg_m2 <- lmer(log(CORT) ~ partner_CORT + site + body_cond + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long2)

fCORT_coreg_m2 <- lmer(log(CORT) ~ site + body_cond + (1|pair_ID) + (1|year_sampled), data = fCORT_coreg_long2)

summary(fCORT_coreg_m2)
Anova(fCORT_coreg_m2, type = "II")
Anova(fCORT_coreg_m2, type = "III")
vif(fCORT_coreg_m2)
confint(fCORT_coreg_m2)

simulationOutput1 <- simulateResiduals(fittedModel = fCORT_coreg_m2, plot = F)
plot(simulationOutput1)

plotResiduals(simulationOutput1, form = fCORT_coreg$M)
plotResiduals(simulationOutput1,  fCORT_coreg$M)

testUniformity(simulationOutput1)
testDispersion(simulationOutput1)
testOutliers(simulationOutput1)
outliers(simulationOutput1)
testSimulatedResiduals(simulationOutput1)
testQuantiles(simulationOutput1)

interact_plot(fCORT_coreg_m2, pred = partner_CORT, modx = site)

#Fit bivariate model
fCORT_coreg_brm1 <- brm(
  bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE),  # residual correlation estimated
  data = fCORT_coreg,
  family = gaussian())

fCORT_coreg_brm1 <- brm(
  bf_fCORT_coreg_f + bf_fCORT_coreg_m + set_rescor(TRUE),  # residual correlation estimated
  data = fCORT_coreg2,
  family = gaussian())

summary(fCORT_coreg_brm1)

#Parent-offspring ----


plot(all_individuals_juv$parent_fCORT, all_individuals_juv$degree)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$strength)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$growth_rate)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$eigenvector)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$betweenness)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$closeness)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$transitivity)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$mean_strength)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$pc1)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$pc2)

plot(all_individuals_juv$parent_fCORT, all_individuals_juv$parent_edge)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$parent_edge_binary)

plot(all_individuals_juv$growth_rate, all_individuals_juv$degree)
plot(all_individuals_juv$growth_rate, all_individuals_juv$strength)
plot(all_individuals_juv$parent_fCORT, all_individuals_juv$visit_number)
plot(all_individuals_juv$visit_number, all_individuals_juv$degree)
plot(all_individuals_juv$visit_number, all_individuals_juv$strength)

cor.test(all_individuals_juv$parent_fCORT, all_individuals_juv$degree)
cor.test(all_individuals_juv$parent_fCORT, all_individuals_juv$strength)
cor.test(all_individuals_juv$parent_fCORT, all_individuals_juv$growth_rate)
cor.test(all_individuals_juv$growth_rate, all_individuals_juv$degree)
cor.test(all_individuals_juv$growth_rate, all_individuals_juv$strength)
cor.test(all_individuals_juv$parent_fCORT, all_individuals_juv$visit_number)











parent_fCORT_m1 <- glmmTMB(degree ~  parent_fCORT + visit_number + (1|parent_ID), data = all_individuals_juv, family = poisson())

parent_fCORT_m2.1 <- lmer(strength ~   parent_fCORT_z + visit_number_z + (1|parent_ID), data = all_individuals_juv[!is.na(all_individuals_juv$parent_fCORT_z),])
parent_fCORT_m2.2 <- lmer(strength ~   parent_fCORT_z + (1|parent_ID), data = all_individuals_juv[!is.na(all_individuals_juv$parent_fCORT_z),])
parent_fCORT_m2.3 <- lmer(strength ~   visit_number_z + (1|parent_ID), data = all_individuals_juv[!is.na(all_individuals_juv$parent_fCORT_z),])

#parent_fCORT_m3 <- lmer(log(mean_strength) ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
parent_fCORT_m4 <- lmer(eigenvector ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
parent_fCORT_m5 <- lmer(betweenness ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
parent_fCORT_m6 <- lmer(social_diff ~  parent_fCORT + visit_number + (1|parent_ID), data = all_individuals_juv)
#parent_fCORT_m7 <- lmer(closeness ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
#parent_fCORT_m8 <- lmer(transitivity ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
#parent_fCORT_m9 <- lmer(pc1 ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
#parent_fCORT_m10 <- lmer(pc2 ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
parent_fCORT_m11 <- glmmTMB(visit_number ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv, family = poisson)
#parent_fCORT_m12 <- lmer(body_condition ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
#parent_fCORT_m13 <- lmer(parent_edge ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
parent_fCORT_m14 <- glmmTMB(parent_edge_binary ~  parent_fCORT + visit_number + (1|parent_ID), data = all_individuals_juv, family = binomial)

#Role of parental visits etc. 
parent_fCORT_m15 <- lmer(log(CORT_pg.mg) ~  visit_nr_own_box + body_cond_noNA + (1|pair_ID_yr_sampled), data = fCORT_22)
parent_fCORT_m16 <- lmer(log(CORT_pg.mg) ~  visit_nr_own_box + (1|pair_ID_yr_sampled), data = fCORT_22)
parent_fCORT_m17 <- lmer(log(CORT_pg.mg) ~  body_cond_noNA + (1|pair_ID_yr_sampled), data = fCORT_22)
parent_fCORT_m18 <- lmer(body_cond_noNA ~ visit_nr_own_box + (1|pair_ID_yr_sampled), data = fCORT_22)
parent_fCORT_m19 <- glmmTMB(visit_nr_own_box ~  body_cond_noNA + (1|pair_ID_yr_sampled), data = fCORT_22, family = "negative.binomial")

parent_fCORT_m20 <- lmer(body_condition ~  parent_visits + (1|parent_ID), data = all_individuals_juv)

#parent_fCORT_m8 <- lmer(growth_rate ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv)
#parent_fCORT_m10 <- glmmTMB(hatch_day ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv, family = poisson)
#parent_fCORT_m11 <- glmmTMB(brood_size ~  parent_fCORT + (1|parent_ID), data = all_individuals_juv, family = poisson)

summary(parent_fCORT_m1)
Anova(parent_fCORT_m1)

summary(parent_fCORT_m2.2)
Anova(parent_fCORT_m2.2)

summary(parent_fCORT_m3)
Anova(parent_fCORT_m3)

summary(parent_fCORT_m4)
Anova(parent_fCORT_m4)

summary(parent_fCORT_m5)
Anova(parent_fCORT_m5)

summary(parent_fCORT_m6)
Anova(parent_fCORT_m6)

summary(parent_fCORT_m7)
Anova(parent_fCORT_m7)

summary(parent_fCORT_m8)
Anova(parent_fCORT_m8)

summary(parent_fCORT_m9)
Anova(parent_fCORT_m9)

summary(parent_fCORT_m10)
Anova(parent_fCORT_m10)

summary(parent_fCORT_m11)
Anova(parent_fCORT_m11)

summary(parent_fCORT_m12)
Anova(parent_fCORT_m12)

summary(parent_fCORT_m13)
Anova(parent_fCORT_m13)

summary(parent_fCORT_m14)
Anova(parent_fCORT_m14)

summary(parent_fCORT_m15)
Anova(parent_fCORT_m15)

summary(parent_fCORT_m16)
Anova(parent_fCORT_m16)

summary(parent_fCORT_m17)
Anova(parent_fCORT_m17)

summary(parent_fCORT_m18)
Anova(parent_fCORT_m18)

summary(parent_fCORT_m19)
Anova(parent_fCORT_m19)

summary(parent_fCORT_m20)
Anova(parent_fCORT_m20)

simulationOutput1 <- simulateResiduals(fittedModel = parent_fCORT_m12, plot = F)
plot(simulationOutput1)

AICc(parent_fCORT_m1, parent_fCORT_m2, parent_fCORT_m3, parent_fCORT_m4, parent_fCORT_m5, parent_fCORT_m6, parent_fCORT_m7, parent_fCORT_m8, parent_fCORT_m9, parent_fCORT_m10, parent_fCORT_m11, parent_fCORT_m12)
AICctab(parent_fCORT_m1, parent_fCORT_m2, parent_fCORT_m3, parent_fCORT_m4, parent_fCORT_m5, parent_fCORT_m6, parent_fCORT_m7, parent_fCORT_m8, parent_fCORT_m9, parent_fCORT_m10, parent_fCORT_m11, parent_fCORT_m12)

#Fitness
fledge_nr_glmm <- glmmTMB(fledge_nr_yr_sampled ~ CORT_pg.mg + body_cond_noNA + (1|JID) + (1|feather_year_sampled) + (1|pair_ID_yr_sampled), data = fCORT, family = compois())

fledge_nr_lmm <- lmer(fledge_nr ~ CORT_pg.mg + body_cond_noNA + (1|JID) + (1|feather_year_sampled) + (1|pair_ID_yr_sampled), data = fCORT)

fledge_nr_glmm <- glmmTMB(fledge_nr_yr_sampled ~ CORT_pg.mg * pbvid_mean + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT, family = compois())

fledge_nr_glmm <- glmmTMB(fledge_nr_yr_grown ~ CORT_pg.mg + site + body_cond_noNA + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT, family = compois())

summary(fledge_nr_lmm)
Anova(fledge_nr_lmm, type = "II")
confint(fledge_nr_lmm)

simulationOutput1 <- simulateResiduals(fittedModel = fledge_nr_glmm, plot = F)
plot(simulationOutput1)

fledge_weight_lmm <- lmer(fledge_weight_yr_sampled ~ CORT_pg.mg + body_cond_noNA + (1|JID) + (1|feather_year_sampled) + (1|pair_ID_yr_sampled), data = fCORT)
fledge_weight_lmm <- lmer(fledge_weight_yr_grown ~ CORT_pg.mg + site + body_cond_noNA + (1|JID) + (1|pair_ID_yr_sampled), data = fCORT)
fledge_weight_lmm <- lmer(fledge_weight ~ CORT_pg.mg + body_cond_noNA + (1|JID) + (1|feather_year_sampled) + (1|pair_ID_yr_sampled), data = fCORT)

summary(fledge_weight_lmm)
Anova(fledge_weight_lmm, type = "II")
confint(fledge_weight_lmm)

simulationOutput1 <- simulateResiduals(fittedModel = fledge_weight_lmm, plot = F)
plot(simulationOutput1)