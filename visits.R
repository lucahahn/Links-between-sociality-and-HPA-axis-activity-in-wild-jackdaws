#RFID Feeder Visit Data from CJP (2022)
#Author: Luca Hahn
#Last update: 11/08/2025

#(1) IMPORT, CLEANING ----

#Load packages
install.packages("asnipe")
install.packages("car")
install.packages("carData")
install.packages("chisq.posthoc.test")
install.packages("ClusterR")
install.packages("corrplot")
install.packages("data.table")
install.packages("DHARMa")
install.packages("dplyr")
install.packages("emmeans2")
install.packages("ggplot2")
install.packages("glmmTMB")
install.packages("hms")
install.packages("igraph")
install.packages("lme4")
install.packages("lmerTest")
install.packages("lubridate")
install.packages("MASS")
install.packages("multcomp")
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("rptR")
install.packages("scales")
install.packages("stringi")
install.packages("stringr")
install.packages("svMisc")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("assortnet")

library(asnipe)
library(car)
library(carData)
library(chisq.posthoc.test)
library(ClusterR)
library(corrplot)
library(data.table)
library(DHARMa)
library(dplyr)
library(emmeans)
library(ggplot2)
library(glmmTMB)
library(hms)
library(igraph)
library(lme4)
library(lmerTest)
library(lubridate)
library(MASS)
library(multcomp)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(rptR)
library(scales)
library(stringi)
library(stringr)
library(svMisc)
library(tidyr)
library(tidyverse)
library(assortnet)

#Shows milliseconds
op <- options(digits.secs=1)

#Use directory where you want to look for RT files, concatenate paths
RT_files <- list.files("passive2", pattern = "RT", recursive = TRUE)
RT_paths <- paste("passive2",RT_files, sep = "/")

#Load saved life history csv file 
LH <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 2/Data/Ch 2 Social Bonds, Social Support, Social Foraging/Data/LH20231004.csv", header = T, stringsAsFactors = F)
LH$DATE <- strptime(LH$DATE,format="%d/%m/%Y")
LH$DATE <- as.Date(LH$DATE, format = "%d/%m/%Y") # convert to date

#Separate LH file with (1) pairs and then (2) adding pair ID
#(i) LH file with pairs and boxes
LH_pairs <- subset(LH[c("DATE", "ID", "SEX", "BOX", "PARTNER.ID")])
LH_pairs <- LH_pairs[LH_pairs$PARTNER.ID != "",]
length(table(LH_pairs$ID)) #570 individuals with pair data

#get most recent entries for pair ID 
LH_pairs <- LH_pairs %>%
  group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

#(ii) Adding pair ID
LH_pairs  <- LH_pairs %>% 
  as_tibble() %>% 
  mutate(pair_ID = if_else(LH_pairs$SEX == "F", paste(ID, PARTNER.ID), paste(PARTNER.ID, ID)))

#LH RFID
LH_RFID <- subset(LH, !LH$RFID == "")

#Create sub-strings that contain feeder ID (e.g. "Y1.1") and date
day_arrays <- unique(substr(RT_files,5,14))

#Create empty list to place data into as we go
collapsed_list <- list()

#Run through RT files in list 
for(i in 1:length(day_arrays)) {
  progress(i, max.value = length(day_arrays))
  
  day_array_list <- list()
  
  for (j in 1:length(RT_files[which(substr(RT_files,5,14) == day_arrays[i])])) {
    temp_day_file <- (read.delim(RT_paths[which(substr(RT_files,5,14) == day_arrays[i])][j], header = T, stringsAsFactors = F))[,1:13]
    temp_day_file$feeder <- substr(RT_files[which(substr(RT_files,5,14) == day_arrays[i])][j],5,7)
    day_array_list[[j]] <- temp_day_file
  }
  
  temp_RT <- do.call(rbind, day_array_list)
  
  temp_RT$Time <- strptime(paste(temp_RT$Date,stri_sub(temp_RT$Hmsec/1024,2,5), sep = ""), "%Y-%m-%d %H:%M:%OS")  # Add in miliseconds (1024 in a second) and format time
  
  temp_RT %>% filter(nchar(TagID_hex) == 10) -> temp_tags  #remove times when no tag
  
  
  if(dim(temp_tags)[1] >0){  
    
    visits <- data.frame(Event = temp_tags$Event, Start = temp_tags$Time, End = temp_tags$Time+(0.5*(temp_tags$Reps -1)), tag = temp_tags$TagID_hex, feeder = temp_tags$feeder)  # adds reps to visit length (0.5 seconds for every extra detection as that was resampling speed)
    visits <- arrange(visits, Start)
    visits <- arrange(visits, feeder)
    within_errors <- which(visits$End < lag(visits$End) & visits$tag == lag(visits$tag) & visits$feeder == lag(visits$feeder))
    if(length(within_errors) > 0){
      visits <- visits[-which(visits$End < lag(visits$End) & visits$tag == lag(visits$tag) & visits$feeder == lag(visits$feeder)),]  ## Get rid of reads within bouts - almost always erroneous single reads that are repeated later in the datastream
    }
    visits <- arrange(visits, Start)
    visits$count <- sapply(1:nrow(visits),function(x)sum(visits$tag[x]==visits$tag[1:x]))  ## add individual visit counter - if two bouts don't have sequential counts then bird seen elsewhere in between
    visits <- arrange(visits, feeder)
    
    visits$collapse <- "0"  # Temp column to tell me if this bout is to be collapsed
    
    visit_time <- 10 ## How long between detections before a new bout if classed
    
    # Is the last visit ending within 'visit_time' seconds of this one starting, and with the same tag & in sequence?
    visits$collapse[ which (visits$Start-lag(visits$End) < visit_time & lead(visits$Start) - visits$End < visit_time & visits$feeder == lag(visits$feeder) & visits$feeder == lead(visits$feeder) & visits$tag == lag(visits$tag) & visits$tag == lead(visits$tag) & (visits$count-lag(visits$count)) == 1 & (lead(visits$count)-visits$count) == 1 )] <- "1"  ## What about times they hop in between?!
    
    visits$collapse[which(visits$collapse == 0 & lag(visits$collapse == 1))] <- "End"  ## Can work out the end based on the 0s and 1s
    visits$collapse[which(visits$collapse == 0 & lead(visits$collapse == 1))] <- "Start" ## Can work out the start based on the 0s and 1s
    
    ## Mop up those that only have two potential detections (so no middle values to get assigned 1 above)
    visits$collapse[which(visits$collapse == 0 & lead(visits$Start) - visits$End < visit_time & lead(visits$Start) - visits$End > -2 & visits$tag == lead(visits$tag) & lead(visits$count) - visits$count == 1)] <- "Start"
    visits$collapse[which(visits$collapse == 0 & visits$Start - lag(visits$End) < visit_time & lag(visits$feeder) == visits$feeder & visits$tag == lag(visits$tag) & lag(visits$count) - visits$count == -1)] <- "End"
    
    # Add this file's data to the list
    collapsed_list[[i]] <- data.frame(start = visits$Start[which(visits$collapse %in% c("Start","0"))], end = visits$End[which(visits$collapse %in% c("End","0"))], tag = visits$tag[which(visits$collapse %in% c("Start","0"))], event = visits$Event[which(visits$collapse %in% c("Start","0"))], feeder =  visits$feeder[which(visits$collapse %in% c("Start","0"))], array = rep(stri_sub(day_arrays[i], 8,11), length(which(visits$collapse %in% c("Start","0")))))  }
}

#Turn the list into a data frame
visit_data <- do.call(rbind,collapsed_list)

#Find instances in which time is still NA
visit_data_NA <- subset(visit_data, is.na(visit_data$start))

#Filter instances in which time is not NA
visit_data <- subset(visit_data, !is.na(visit_data$start))

#Remove instances where JID = NA and test tags
visit_data$JID      <- LH$ID[match(visit_data$tag,LH$RFID)]
visit_data <- subset(visit_data, JID != "NA")



#For importing and concatenating raw RT data
RT_files <- list.files("Data/passive2", pattern = "RT", recursive = TRUE)
RT_paths <- paste("Data/passive2",RT_files, sep = "/")

day_arrays <- unique(substr(RT_files,5,14))

collapsed_list <- list()

for(i in 1:length(day_arrays)) {
  progress(i, max.value = length(day_arrays))
  
  day_array_list <- list()
  
  for (j in 1:length(RT_files[which(substr(RT_files,5,14) == day_arrays[i])])) {
    temp_day_file <- (read.delim(RT_paths[which(substr(RT_files,5,14) == day_arrays[i])][j], header = T, stringsAsFactors = F))[,1:13]
    temp_day_file$feeder <- substr(RT_files[which(substr(RT_files,5,14) == day_arrays[i])][j],5,7)
    day_array_list[[j]] <- temp_day_file
  }
  
  temp_RT <- do.call(rbind, day_array_list)
  
  temp_RT$Time <- strptime(paste(temp_RT$Date,stri_sub(temp_RT$Hmsec/1024,2,5), sep = ""), "%Y-%m-%d %H:%M:%OS")  # Add in miliseconds (1024 in a second) and format time
  
  temp_RT %>% filter(nchar(TagID_hex) == 10) -> temp_tags
  
  collapsed_list[[i]] <- data.frame(tag = temp_RT$TagID_hex, feeder =  temp_RT$feeder, time = temp_RT$Time)
}

visit_data_RT <- do.call(rbind,collapsed_list)

visit_data_RT <- subset(visit_data_RT, !visit_data_RT$tag == "")

visit_data_RT$JID <- LH_RFID$ID[match(visit_data_RT$tag, LH_RFID$RFID)]

sum(is.na(visit_data_RT$JID))

visit_data_RT <- subset(visit_data_RT, !is.na(visit_data_RT$JID))

visit_data_RT$day <- yday(visit_data_RT$time) #day of year

# (2) ADDITIONAL INFORMATION ---- 

#Adding information about visit duration
visit_data$interval <- interval(visit_data$start,visit_data$end)
visit_data$visit_duration <- as.duration(visit_data$interval)

#Average visit duration per individual data set
visit_duration <- as.data.frame(aggregate(visit_data$visit_duration, by = list(visit_data$JID), FUN = "mean", na.rm = TRUE))
visit_duration$JID <- visit_duration$Group.1
visit_duration$visit_duration <- visit_duration$x
visit_duration <- subset(visit_duration, select = -c(Group.1, x))

#Visit number per individual data set
visits <- as.data.frame(table(visit_data$JID))
visits$JID <- visits$Var1
visits$visit_number <- visits$Freq
visits <- subset(visits, select = -c(Var1, Freq))

#visit number and duration per individual data set
visits_per_indiv <- merge(x = visits, y = visit_duration, by = "JID", all = TRUE)

#Adding information about site: Y, Z
visit_data[substr(visit_data$feeder,1,1)=="Y","site"]<-"Y"
visit_data[substr(visit_data$feeder,1,1)=="Z","site"]<-"Z"

#Adding information about feeder position 
visit_data[substr(visit_data$feeder,1,1)=="1","position"]<-"1"
visit_data[substr(visit_data$feeder,1,1)=="2","position"]<-"2"
visit_data[substr(visit_data$feeder,1,1)=="3","position"]<-"3"
visit_data[substr(visit_data$feeder,1,1)=="4","position"]<-"4"

#Adding information about feeder position 
visit_data_RT[substr(visit_data_RT$feeder,1,1)=="1","position"]<-"1"
visit_data_RT[substr(visit_data_RT$feeder,1,1)=="2","position"]<-"2"
visit_data_RT[substr(visit_data_RT$feeder,1,1)=="3","position"]<-"3"
visit_data_RT[substr(visit_data_RT$feeder,1,1)=="4","position"]<-"4"

#Adding day of year and day of study period 
visit_data$day <- yday(visit_data$start) #day of year

#Remove "array" column
visit_data <- subset(visit_data, select = -c(array))

#Adding information about time of the day 
visit_data$time <- as_hms(visit_data$start)
visit_data$hour <- hour(visit_data$start)

#Adding individual jackdaw information 
visit_data$JID      <- LH$ID[match(visit_data$tag,LH$RFID)]
visit_data$rings    <- LH$COMBINATION[match(visit_data$tag,LH$RFID)]
LH_sex <- LH[LH$SEX != "",]
visit_data$sex      <- LH_sex$SEX[match(visit_data$tag,LH_sex$RFID)]   
LH_partner <- LH_pairs[LH_pairs$PARTNER.ID != "",]
LH_mother <- LH[LH$MOTHER.ID != "",]
LH_father <- LH[LH$FATHER.ID != "",]
LH_box <- LH[LH$BOX != "",]
LH_box$year <- as.numeric(format(LH_box$DATE,"%Y"))   
LH_box <- LH_box[LH_box$year == "2022",]
visit_data$partnerID <- LH_partner$PARTNER.ID[match(visit_data$JID,LH_partner$ID)]
visit_data$pairID <- LH_pairs$pair_ID[match(visit_data$JID, LH_pairs$ID)]
visit_data$box <- LH_box$BOX[match(visit_data$JID, LH_box$ID)] 
visit_data$motherID <- LH_mother$MOTHER.ID[match(visit_data$tag,LH$RFID)]
visit_data$fatherID <- LH_father$FATHER.ID[match(visit_data$tag,LH$RFID)]

#Adding individual age
Current_year <- 2022  #Set reference year

LH$DATE = as.Date(LH$DATE, format="%d-%m-%Y")
LH$year <- as.numeric(format(LH$DATE,"%Y"))   #  #Get year from the date
Ringed <- LH %>% filter(CODE == "RINGED")   #Get only records of when birds were ringed for the first time

Ringed$known_age <- 0   #binary 0/1 do we know the exact age (e.g. birds ringed as a 6 are 0)
Ringed$min_age <- 0  #Either actual age (if known_age = 1), or minimum age (if known_age = 0) - currently as number of new years crossed.

for (i in  1:nrow(Ringed)) {
  if(Ringed[i,11] == "4"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2)
  }
  else if(Ringed[i,11] %in% c("1","1J","3","3J")){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i]+1) 
  }
  else if(Ringed[i,11] == "5"){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2) 
  }
  else if (Ringed[i,11] == "6"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +3) 
  }
  else { Ringed$known_age[i] = 0
  Ringed$min_age[i] = NA }
}

sum(Ringed$known_age)
length(Ringed$known_age)

unique_ID = unique(visits$JID)
individuals = data.frame(JID = unique_ID)
individuals$known_age = Ringed$known_age[match(individuals$JID,Ringed$ID)]
individuals$min_age = Ringed$min_age[match(individuals$JID,Ringed$ID)]
individuals$sex = Ringed$SEX[match(individuals$JID,Ringed$ID)]
individuals$box = visits$visitor_owned_box[match(individuals$JID,visits$JID)]
visit_data$age <- individuals$min_age[match(visit_data$JID,individuals$JID)]

#Subset data set by site 
visit_data_Y <- subset(visit_data, site == 'Y')
visit_data_Z <- subset(visit_data, site == 'Z')

#Subset data set by location 
visit_data_1 <- subset(visit_data, position == '1')
visit_data_2 <- subset(visit_data, position == '2')
visit_data_3 <- subset(visit_data, position == '3')
visit_data_4 <- subset(visit_data, position == '4')

#Subset data set by sex 
visit_data_F <- subset(visit_data, sex == 'F')
visit_data_M<- subset(visit_data, sex == 'M')
visit_data_N <- subset(visit_data,is.na(visit_data$sex))
visit_data_J <- subset(visit_data,age == 1)
visit_data_J2 <- subset(visit_data_J,!is.na(visit_data_J$box))
visit_data_J2 <- visit_data_J2[, c("JID", "rings")]
write.csv(visit_data_J2,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 1/Data/PhD Ch 1 Social Bonds and Stress/visit_data_J2.csv", row.names = FALSE)

#Data from other years
visit_data_J2020 <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/CJP RFID feeders/visit_data_20_J.csv", header = T, stringsAsFactors = F)
visit_data_J2019 <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/CJP RFID feeders/visit_data_19_J.csv", header = T, stringsAsFactors = F)

visit_data_test <- visit_data %>% distinct(JID, .keep_all=TRUE)
visit_data_test <- merge(x = visit_data_test, y = visit_data_J2019, by = "JID", all.x = TRUE)
visit_data_test <- subset(visit_data_test, !is.na(visit_data_test$age.y))
length(table(visit_data_test$JID))

#(3) SUMMARY AND EXPLORATION ----

#Number of individuals in total
length(table(visit_data$JID)) #315 individuals
table(visit_data$JID)
length(table(visit_data$tag))
visit_data2 <- subset(visit_data[, c("tag", "JID")])
visit_data2_unique <- unique(visit_data2$JID) 

#individuals per site
length(table(visit_data_Y$JID)) 
length(table(visit_data_Z$JID)) 

#Number of individuals per feeder location
length(table(visit_data_1$JID)) #84
length(table(visit_data_2$JID)) #156
length(table(visit_data_3$JID)) #108
length(table(visit_data_4$JID)) #116

#Number of individuals per sex 
length(table(visit_data_F$JID)) #122
length(table(visit_data_M$JID)) #170
length(table(visit_data_N$JID)) #23

#Number of visits per site, feeder, position, perch, sex
table(visit_data$position)
table(visit_data$feeder)
table(visit_data$sex)#F = 12430, M = 31715

#Visit duration 
mean(visit_data$visit_duration, na.rm =TRUE) #9.40 s
median(visit_data$visit_duration, na.rm =TRUE) #4.27 s
sd(visit_data$visit_duration, na.rm =TRUE) #14.50 s

#Number of visits per individual per site
table(visit_data$JID, visit_data$feeder)
visit_data_sum <- data.frame(table(visit_data$feeder, visit_data$JID))

#Visits per individual per site and per feeder, preferred sites
perindivpersite <- visit_data  %>%  count(JID, site)
perindivpersite2 <- reshape(perindivpersite, idvar = "JID", timevar = "site", direction = "wide")
perindivpersite3 <- perindivpersite2
perindivpersite3[is.na(perindivpersite3)] <- 0.99
perindivpersite3$pref <- perindivpersite3$n.Y / perindivpersite3$n.Z
perindivpersite3 <- perindivpersite3 %>% 
  as_tibble() %>% 
  mutate(pref.site = if_else(pref > 0.5,"Y", "Z"))
perindivpersite3  <- subset(perindivpersite3 [,c(1,6)])

#Individuals that were detected at different sites
perindivperpos <- visit_data %>% count(JID, position)
indiv_1site <- perindivpersite %>% distinct(JID, .keep_all=TRUE)
indiv_2sites <- perindivpersite[perindivpersite$JID %in% perindivpersite$JID[duplicated(perindivpersite$JID)],]
length(table(indiv_2sites$JID)) #91 individuals detected at Y and Z

#Individuals' preferred position
perindivperpos2 <- dcast(setDT(perindivperpos), JID ~ position, value.var = "n")
perindivperpos2$pref_pos <- colnames(perindivperpos2)[apply(perindivperpos2,1,which.max)]

#Individual attribute data set 
n <- 315 #number of individuals in study
all_individuals <- visits_per_indiv
all_individuals$age <- individuals$min_age[match(all_individuals$JID,individuals$JID)]
all_individuals$sex <- LH_sex$SEX[match(all_individuals$JID,LH_sex$ID)]
all_individuals$pair_ID <- LH_pairs$pair_ID[match(all_individuals$JID,LH_pairs$ID)]

all_individuals_fCORT <- subset(all_individuals, !is.na(all_individuals$fCORT))

# (4) SOCIAL NETWORK ----

#Gaussian Mixture Models 

#Prepare variables 
## (i) Time stamp: seconds since start of study period (1 s before first visit)
int <- interval(ymd_hms("2022-07-11 16:45:00 UTC"), ymd_hms(visit_data$start))
visit_data$time <- time_length(int, "second")

int <- interval(ymd_hms("2022-07-11 16:45:00 UTC"), ymd_hms(visit_data_RT$time))
visit_data_RT$time <- time_length(int, "second")

## (ii) Identity 
visit_data$JID
global_ids <- visits$JID

visit_data_RT$JID
global_ids <- sort(unique(visit_data_RT$JID))

## (iii) Location 
visit_data$position

visit_data_RT$position

## (iv) Location in time
visit_data$loc_date <-
  paste(visit_data$position,
        visit_data$day,sep="_")

visit_data_RT$loc_date <-
  paste(visit_data_RT$position,
        visit_data_RT$day,sep="_")

## Data set 
visit_data_gmm <- visit_data[, c("time", "JID", "position", "loc_date")]
visit_data_gmm <- na.omit(visit_data_gmm)

visit_data_RT_gmm <- visit_data_RT[, c("time", "JID", "position", "loc_date")]
visit_data_RTgmm <- na.omit(visit_data_RT_gmm)

# Generate GMM data
gmm_data <- gmmevents(time= visit_data_gmm$time,
                      identity=visit_data_gmm$JID,
                      location=visit_data_gmm$loc_date,
                      global_ids=global_ids)

gmm_data_RT <- gmmevents(time= visit_data_RT_gmm$time,
                      identity=visit_data_RT_gmm$JID,
                      location=visit_data_RT_gmm$loc_date,
                      global_ids=global_ids)

# Extract output
gbi <- gmm_data$gbi
events <- gmm_data$metadata
observations_per_event <- gmm_data$B

gbi <- gmm_data_RT$gbi
events <- gmm_data_RT$metadata
observations_per_event <- gmm_data_RT$B

# Can also subset gbi to only individuals observed
# in the dataset to give same answer as if
# global_ids had not been provided
gbi <- gbi[,which(colSums(gbi)>0)]

gbi_data <- as.data.frame(gbi)
gbi_data$group_size <- rowSums(gbi_data)

gbi_individuals <- setdiff(colnames(gbi_data), "group_size")

gbi_indiv_group_sizes <- sapply(gbi_individuals, function(id) {
  mean(gbi_data$group_size[gbi_data[[id]] == 1], na.rm = TRUE)
})

gbi_indiv_group_sizes

all_individuals$group_size <- gbi_indiv_group_sizes

# Split up location and date data
tmp <- strsplit(events$Location,"_")
tmp <- do.call("rbind",tmp)
events$Location <- tmp[,1]
events$Date <- tmp[,2]

#Get the adjacency matrix from group by individual matrix 
am <- get_network(gbi, data_format = "GBI",
                  association_index = "SRI", identities = NULL,
                  which_identities = NULL, times = NULL, occurrences = NULL,
                  locations = NULL, which_locations = NULL, start_time = NULL,
                  end_time = NULL, classes = NULL, which_classes = NULL,
                  enter_time = NULL, exit_time = NULL)

am_data <- as.data.frame(am)
am_data$mean_strength <- rowMeans(am_data)
am_data$sd_strength <- apply(am_data, 1, sd, na.rm=TRUE)
am_data$social_diff <- am_data$sd_strength / am_data$mean_strength
am_data$JID <- row.names(am_data)

#Get graph from the adjacency matrix
g <- graph_from_adjacency_matrix(am, mode= "undirected",weighted=TRUE,diag=FALSE)

g_dyads <- as_ids(E(g))
g_edge_weights <- E(g)$weight

dyad_edge <- data.frame(matrix(NA, nrow = 6725, ncol = 2))
dyad_edge <- data.frame(matrix(NA, nrow = 6300, ncol = 2))

dyad_edge$X1 <- g_dyads
dyad_edge$X2 <- g_edge_weights

dyad_edge$g_dyads <- chartr("|", " ", dyad_edge$X1)
dyad_edge$g_edge_weights <- dyad_edge$X2

dyad_edge <- dyad_edge[, c("g_dyads", "g_edge_weights")]

dyad_edge$pair_ID1 <- dyad_edge$g_dyads

dyad_edge$pair_ID2 <- dyad_edge$pair_ID1

dyad_edge$ind1 <- str_split(dyad_edge$pair_ID1, " ", simplify = TRUE)
dyad_edge$ind1 <- dyad_edge$ind1[,1]

dyad_edge$ind2 <- str_split(dyad_edge$pair_ID2, " ", simplify = TRUE)
dyad_edge$ind2 <- dyad_edge$ind2[,2]

dyad_edge$pair_ID2 <- paste(dyad_edge$ind2,dyad_edge$ind1,sep=" ")

#Get sub-graph
edge_weight <- E(g)$weight
summary(edge_weight)
quantile(edge_weight, probs = seq(0, 1, 1/100))
mean(E(g)$weight) * 2

gs <- subgraph.edges(g, E(g)[E(g)$weight >  0.019607843], del=F) #top 50%
gs <- subgraph.edges(g, E(g)[E(g)$weight >  0.035714286], del=F) #top 25%
gs <- subgraph_from_edges(g, E(g)[E(g)$weight >  0.055555556], del=F) #top 10%

gs2 <- delete_vertices(gs, degree(gs)==0)

gfCORT <- subgraph(g, V(g)[!is.na(V(g)$fCORT)])
amfCORT <- as_adjacency_matrix(gfCORT, attr = "weight") 

gbonded <- subgraph_from_edges(g, E(g)[!E(g)$relationship == "other"])

?subgraph

#vertex attributes
g <- set_vertex_attr(g, "sex", value = all_individuals$sex)
g <- set_vertex_attr(g, "age", value = all_individuals$age)
g <- set_vertex_attr(g, "pref_pos", value = all_individuals$pref_pos)
g <- set_vertex_attr(g, "pref_site", value = all_individuals$pref_site)
g <- set_vertex_attr(g, "fCORT", value = all_individuals$fCORT)
g <- set_vertex_attr(g, "strength", value = all_individuals$strength)

#Edge attributes 
g <- set_edge_attr(g , "kin", value = dyad_edge$parent_offspring_kin)
g <- set_edge_attr(g , "relationship", value = dyad_edge$relationship)

E(g)$kin

#Plot social network
V(g)$colour <- ifelse(V(g)$sex == "M", "lightblue", "orange")
V(g)$colour <- ifelse(V(g)$pref_pos == "Y4", "lightblue", "orange")
V(g)$colour <- ifelse(V(g)$pref_site == "Y", "lightblue", "orange")

V(g)$colour <- ifelse(V(g)$sex == "M",  "#009E73", "#E69F00")
V(gs)$colour <- ifelse(V(gs)$sex == "M",  "#009E73", "#E69F00")
V(gs2)$colour <- ifelse(V(gs2)$sex == "M",  "#009E73", "#E69F00")
V(gfCORT)$colour <- ifelse(V(gfCORT)$sex == "M",  "#009E73", "#E69F00")

V(g)$colour <- ifelse(V(g)$fCORT < 2,"#FFFFFF", ifelse(V(g)$fCORT > 2 & V(g)$fCORT < 3, "#CCCCCC", ifelse(V(g)$fCORT > 3 & V(g)$fCORT < 4, "#999999", ifelse(V(g)$fCORT > 4 & V(g)$fCORT < 5, "#666666", ifelse(V(g)$fCORT > 5 & V(g)$fCORT < 6,  "#333333", "#000000")))))
V(gs)$colour <- ifelse(V(gs)$fCORT < 2,"#FFFFFF", ifelse(V(gs)$fCORT > 2 & V(gs)$fCORT < 3, "#CCCCCC", ifelse(V(gs)$fCORT > 3 & V(gs)$fCORT < 4, "#999999", ifelse(V(gs)$fCORT > 4 & V(gs)$fCORT < 5, "#666666", ifelse(V(gs)$fCORT > 5 & V(gs)$fCORT < 6,  "#333333", "#000000")))))

V(g)$colour <- ifelse(is.na(V(g)$fCORT),"#FFFFFF", ifelse(V(g)$fCORT < 2,"#EBEBEB", ifelse(V(g)$fCORT > 2 & V(g)$fCORT < 3, "#D4D4D4", ifelse(V(g)$fCORT > 3 & V(g)$fCORT < 4, "#BABABA", ifelse(V(g)$fCORT > 4 & V(g)$fCORT < 5, "#9B9B9B", ifelse(V(g)$fCORT > 5 & V(g)$fCORT < 6,  "#717171", "#000000"))))))
V(gs)$colour <- ifelse(is.na(V(gs)$fCORT),"#FFFFFF", ifelse(V(gs)$fCORT < 2,"#EBEBEB", ifelse(V(gs)$fCORT > 2 & V(gs)$fCORT < 3, "#D4D4D4", ifelse(V(gs)$fCORT > 3 & V(gs)$fCORT < 4, "#BABABA", ifelse(V(gs)$fCORT > 4 & V(gs)$fCORT < 5, "#9B9B9B", ifelse(V(gs)$fCORT > 5 & V(gs)$fCORT < 6,  "#717171", "#000000"))))))
V(gfCORT)$colour <- ifelse(is.na(V(gfCORT)$fCORT),"#FFFFFF", ifelse(V(gfCORT)$fCORT < 2,"#EBEBEB", ifelse(V(gfCORT)$fCORT > 2 & V(gfCORT)$fCORT < 3, "#D4D4D4", ifelse(V(gfCORT)$fCORT > 3 & V(gfCORT)$fCORT < 4, "#BABABA", ifelse(V(gfCORT)$fCORT > 4 & V(gfCORT)$fCORT < 5, "#9B9B9B", ifelse(V(gfCORT)$fCORT > 5 & V(gfCORT)$fCORT < 6,  "#717171", "#000000"))))))

E(g)$colour <- ifelse(E(g)$kin == "kin", "#332288", "grey")
E(gs)$colour <- ifelse(E(gs)$kin == "kin", "#332288", "grey")

E(g)$colour <- ifelse(E(g)$relationship == "pair", "#661100", ifelse(E(g)$relationship == "kin", "#332288", "#EBEBEB"))
E(gs)$colour <- ifelse(E(gs)$relationship == "pair", "#661100", ifelse(E(gs)$relationship == "kin", "#332288", "#EBEBEB"))

curve_multiple(g, start = 1.0)

coords <- layout_(g, nicely())
coords <- layout_(gs, nicely())
coords <- layout_(gs2, nicely())
coords <- layout_(gfCORT, nicely())
coords <- layout_(gbonded, nicely())

plot(g, layout = coords, vertex.size= 4, vertex.label=NA, vertex.color= V(g)$colour, edge.width = E(g)$weight * 10, edge.color = "grey", edge.curved = 0.35)

plot(gs, layout = coords, vertex.size= V(gs)$strength + 1, vertex.label=NA, vertex.color= V(gs)$colour, edge.width = E(gs)$weight * 10, edge.color = "grey", edge.curved = 0.35)

plot(gs, layout = coords, vertex.size= V(gs)$strength + 1, vertex.label=NA, vertex.color= V(gs)$colour, edge.width = E(gs)$weight * 15, edge.color = E(gs)$colour, edge.curved = 0.35)

plot(gbonded, layout = coords, vertex.size= V(gbonded)$strength + 1, vertex.label=NA, vertex.color= V(gbonded)$colour, edge.width = E(gbonded)$weight * 20, edge.color = E(gbonded)$colour, edge.curved = 0.35)

plot(gfCORT, layout = coords, vertex.size= V(gfCORT)$strength + 2, vertex.label=NA, vertex.color= V(gfCORT)$colour, edge.width = E(gfCORT)$weight * 20, edge.color = "grey", edge.curved = 0.35)

plot(gs2, layout = coords, vertex.size= 4, vertex.label=NA, vertex.color= V(gs2)$colour, edge.width = E(gs2)$weight * 10, edge.color = "grey", edge.curved = 0.35)

plot(gfCORT, vertex.size = V(gfCORT)$fCORT * 2, vertex.label = NA, vertex.color = V(gfCORT)$colour, edge.width = E(gfCORT)$weight * 15, edge.curved = 0.35)

#Edge weight and pair bond strength
pairbond2$pair_ID1 <- paste(pairbond2$FID, pairbond2$MID, sep = " ")
pairbond2$pair_ID2 <- paste(pairbond2$MID, pairbond2$FID, sep = " ")

g_dyads <- as_ids(E(g))
g_edge_weights <- E(g)$weight

dyad_edge <- data.frame(matrix(NA, nrow = 6725, ncol = 2))
dyad_edge <- data.frame(matrix(NA, nrow = 6300, ncol = 2))

dyad_edge$X1 <- g_dyads
dyad_edge$X2 <- g_edge_weights

dyad_edge$g_dyads <- chartr("|", " ", dyad_edge$X1)
dyad_edge$g_edge_weights <- dyad_edge$X2

dyad_edge <- dyad_edge[, c("g_dyads", "g_edge_weights")]

dyad_edge$pair_ID1 <- dyad_edge$g_dyads

dyad_edge$pair_ID2 <- dyad_edge$pair_ID1

dyad_edge$ind1 <- str_split(dyad_edge$pair_ID1, " ", simplify = TRUE)
dyad_edge$ind1 <- dyad_edge$ind1[,1]

dyad_edge$ind2 <- str_split(dyad_edge$pair_ID2, " ", simplify = TRUE)
dyad_edge$ind2 <- dyad_edge$ind2[,2]

dyad_edge$pair_ID2 <- paste(dyad_edge$ind2,dyad_edge$ind1,sep=" ")

write.csv(dyad_edge,"C:/Users/lh868/OneDrive - University of Exeter/feedgmm22.csv", row.names = FALSE)

pairbond2$g_edge_weights1 <- dyad_edge$g_edge_weights[match(pairbond2$pair_ID1, dyad_edge$pair_ID1)]
pairbond2$g_edge_weights2 <- dyad_edge$g_edge_weights[match(pairbond2$pair_ID1, dyad_edge$pair_ID2)]

pairbond2$g_edge_weights1[is.na(pairbond2$g_edge_weights1)] <- 0
pairbond2$g_edge_weights2[is.na(pairbond2$g_edge_weights2)] <- 0

pairbond2$g_edge_weights3 <- pairbond2$g_edge_weights1 + pairbond2$g_edge_weights2

pairbond2 <- pairbond2[pairbond2$g_edge_weights3 > 0,]

pairbond$ind1_visits <- all_individuals$visit_number[match(pairbond$JID1, all_individuals$JID)]
pairbond$ind2_visits <- all_individuals$visit_number[match(pairbond$JID2, all_individuals$JID)]

pairbond$ind1_visits[is.na(pairbond$ind1_visits)] <- 0
pairbond$ind2_visits[is.na(pairbond$ind2_visits)] <- 0

pairbond$visits <- pairbond$ind1_visits + pairbond$ind2_visits

plot(pairbond$g_edge_weights3, pairbond$bond_strength)
plot(pairbond2$g_edge_weights3, pairbond2$bond_strength)

#Add mother ID
LH_mother <- LH[LH$MOTHER.ID != "",]

dyad_edge$ind1_motherID <- LH_mother$MOTHER.ID[match(dyad_edge$ind1,LH_mother$ID)]
dyad_edge$ind2_motherID <- LH_mother$MOTHER.ID[match(dyad_edge$ind2,LH_mother$ID)]

#Add initiator and joiner father ID
LH_father <- LH[LH$FATHER.ID != "",]

dyad_edge$ind1_fatherID <- LH_father$FATHER.ID[match(dyad_edge$ind1,LH_mother$ID)]
dyad_edge$ind2_fatherID <- LH_father$FATHER.ID[match(dyad_edge$ind2,LH_mother$ID)]

#Add offspring-parent
dyad_edge <- dyad_edge %>% 
  as_tibble() %>% 
  mutate(parent_offspring = if_else(ind1 == ind2_motherID | ind1 == ind2_fatherID,"parent_offspring", "no"))

#Add parent-offspring
dyad_edge <- dyad_edge %>% 
  as_tibble() %>% 
  mutate(offspring_parent = if_else(ind2 == ind1_motherID | ind2 == ind1_fatherID,"offspring_parent", "no"))

#Add kin (parent and offspring only, no siblings)
dyad_edge <- dyad_edge %>% 
  as_tibble() %>% 
  mutate(parent_offspring_kin = if_else(parent_offspring == "parent_offspring" | offspring_parent == "offspring_parent","kin", "non-kin"))

dyad_edge$parent_offspring_kin <- ifelse(is.na(dyad_edge$parent_offspring_kin), "non-kin", ifelse(dyad_edge$parent_offspring_kin == "kin", "kin", "non-kin"))

#Add pairs 
dyad_edge$ind1_partner <- LH_partner$PARTNER.ID[match(dyad_edge$ind1, LH_partner$ID)]

dyad_edge <- dyad_edge %>% 
  as_tibble() %>% 
  mutate(pair = if_else(ind1_partner == ind2,"pair", "non-pair"))

dyad_edge  <- dyad_edge  %>%
  mutate(relationship =case_when(
    pair =="pair" & parent_offspring_kin == "non-kin" ~ "pair",
    pair == "non-pair" & parent_offspring_kin == "kin" ~ "kin",
    pair == "pair" & is.na(parent_offspring_kin) ~ "pair",
    is.na(pair) & parent_offspring_kin == "kin" ~ "kin",
    pair == "non-pair" & parent_offspring_kin == "non-kin" ~ "other",
    is.na(pair) & is.na(parent_offspring_kin) ~ "other",
    is.na(pair) & parent_offspring_kin == "non-kin"~ "other",
    pair == "non-pair" & is.na(parent_offspring_kin) ~ "other",
    pair == "pair" & parent_offspring_kin == "kin" ~ "pair"))

dyad_edge$relationship <- factor(dyad_edge$relationship, levels=c("kin","pair","other"))

#Centrality measures
degree <- as.data.frame(degree(g))
degree$JID <- rownames(degree)
all_individuals$degree <- degree$`degree(g)`[match(all_individuals$JID, degree$JID)] 

eigenvector <- as.data.frame(eigen_centrality(g))
eigenvector <- eigenvector[, c("vector", "value")]
eigenvector$JID <- rownames(eigen)
eigenvector$eigenvector <- eigenvector$vector
all_individuals$eigenvector <- eigenvector$eigenvector[match(all_individuals$JID, eigenvector$JID)] 

betweenness <- as.data.frame(betweenness(g))
betweenness$JID <- rownames(betweenness)
betweenness$betweenness <- betweenness$`betweenness(g)`
all_individuals$betweenness <- betweenness$betweenness[match(all_individuals$JID, betweenness$JID)] 

strength <- as.data.frame(strength(g))
strength$JID <- rownames(strength)
strength$strength <- strength$`strength(g)` 
all_individuals$strength <- strength$strength[match(all_individuals$JID, strength$JID)] 

closeness <- as.data.frame(closeness(g))
closeness$JID <- rownames(closeness)
closeness$closeness <- closeness$`closeness(g)` 
all_individuals$closeness <- closeness$closeness[match(all_individuals$JID, closeness$JID)] 

transitivity <- as.data.frame(transitivity(g, type = "local"))
transitivity$JID <- rownames(transitivity)
transitivity$transitivity <- transitivity$`transitivity(g)` 
all_individuals$transitivity <- transitivity$transitivity[match(all_individuals$JID, transitivity$JID)] 

all_individuals$mean_strength <- all_individuals$strength /all_individuals$degree

write.csv(all_individuals,"all_individuals_feed22.csv", row.names = FALSE)

#Box visits: box owners and prospectors 


#(5) ARCHIVE ----

plot(gs2, vertex.size =2, vertex.label = NA, vertex.color = "green", edge.width = E(gs2)$weight)
plot(gs2, layout = coords, vertex.size =3, vertex.label = NA, vertex.color = V(gs2)$colour, edge.color = "black", edge.width = E(gs2)$weight * 15, edge.curved = 0.35)

gs <- subgraph.edges(g, E(g)[E(g)$weight > 0.0573237], del=F)
