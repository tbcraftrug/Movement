######################Combined godwit location data##########################################################################################################
# https://gist.github.com/ 

library(move) 
library(tidyverse)

#Movebank login
username = "t.b.craft"
password = "GodwitSnl24!!"
login<-movebankLogin(username,password)

######################reference data######################
ib19 <- as.data.frame(getMovebankReferenceTable(study=652989041,login=login))
microwave2021 <- as.data.frame(getMovebankReferenceTable(study=1498143083,login=login))
extremadura2022 <- as.data.frame(getMovebankReferenceTable(study=1923591036,login=login))
southholland2021 <- as.data.frame(getMovebankReferenceTable(study=1145538280,login=login))
BtgTagus2021 <- as.data.frame(getMovebankReferenceTable(study=1693518103,login=login))
ad_dum2018 <- as.data.frame(getMovebankReferenceTable(study=484019425,login=login))
ad_dum2019 <- as.data.frame(getMovebankReferenceTable(study=831990025,login=login))
ad_dum2020 <- as.data.frame(getMovebankReferenceTable(study=1105026166,login=login))
ad_dum2021 <- as.data.frame(getMovebankReferenceTable(study=1482506572,login=login))
ch_dum2018 <- as.data.frame(getMovebankReferenceTable(study=500187586,login=login))
ch_dum2019 <- as.data.frame(getMovebankReferenceTable(study=878914763,login=login))
ch_dum2020 <- as.data.frame(getMovebankReferenceTable(study=1183466126,login=login))
ch_dum2021 <- as.data.frame(getMovebankReferenceTable(study=1482505185,login=login))
ad_dum2022 <- as.data.frame(getMovebankReferenceTable(study=1751337831,login=login))
ch_dum2022 <- as.data.frame(getMovebankReferenceTable(study=2098519852,login=login))
hmadults <- as.data.frame(getMovebankReferenceTable(study=69402287,login=login))
hrjuv <- as.data.frame(getMovebankReferenceTable(study=76429224,login=login))
hrjuv2016 <- as.data.frame(getMovebankReferenceTable(study=175328223,login=login))
hrjuv2017 <- as.data.frame(getMovebankReferenceTable(study=293970900,login=login))
HQXS_Black_tailed_godwits <- as.data.frame(getMovebankReferenceTable(study=1658294759,login=login))
iberiaBlackwits <- as.data.frame(getMovebankReferenceTable(study= 49547785,login=login))
icarus <- as.data.frame(getMovebankReferenceTable(study=1487044886,login=login))
polish <- as.data.frame(getMovebankReferenceTable(study=163516781,login=login))
teamPiersmaHQXS <- as.data.frame(getMovebankReferenceTable(study=1563249841,login=login))
teamPiersmaHQXS2022 <- as.data.frame(getMovebankReferenceTable(study=2083443328,login=login))
wildjuv <- as.data.frame(getMovebankReferenceTable(study=75360602,login=login))
wildjuv2016 <- as.data.frame(getMovebankReferenceTable(study=170829089,login=login))
wildjuv2017 <- as.data.frame(getMovebankReferenceTable(study=282596404,login=login))

######################location data######################
ib19data <- as.data.frame(getMovebankData(study=652989041,login=login,removeDuplicatedTimestamps=TRUE))
microwave2021data <- as.data.frame(getMovebankData(study=1498143083,login=login,removeDuplicatedTimestamps=TRUE))
extremadura2022data <- as.data.frame(getMovebankData(study=1923591036,login=login,removeDuplicatedTimestamps=TRUE))
southholland2021data <- as.data.frame(getMovebankData(study=1145538280,login=login,removeDuplicatedTimestamps=TRUE))
BtgTagus2021data <- as.data.frame(getMovebankData(study=1693518103,login=login,removeDuplicatedTimestamps=TRUE))
ad_dum2018data <- as.data.frame(getMovebankData(study=484019425,login=login,removeDuplicatedTimestamps=TRUE))
ad_dum2019data <- as.data.frame(getMovebankData(study=831990025,login=login,removeDuplicatedTimestamps=TRUE))
ad_dum2020data <- as.data.frame(getMovebankData(study=1105026166,login=login,removeDuplicatedTimestamps=TRUE))
ad_dum2021data <- as.data.frame(getMovebankData(study=1482506572,login=login,removeDuplicatedTimestamps=TRUE))
ch_dum2018data <- as.data.frame(getMovebankData(study=500187586,login=login,removeDuplicatedTimestamps=TRUE))
ch_dum2019data <- as.data.frame(getMovebankData(study=878914763,login=login,removeDuplicatedTimestamps=TRUE))
ch_dum2020data <- as.data.frame(getMovebankData(study=1183466126,login=login,removeDuplicatedTimestamps=TRUE))
ch_dum2021data <- as.data.frame(getMovebankData(study=1482505185,login=login,removeDuplicatedTimestamps=TRUE))
ad_dum2022data <- as.data.frame(getMovebankData(study=1751337831,login=login,removeDuplicatedTimestamps=TRUE))
ch_dum2022data <- as.data.frame(getMovebankData(study=2098519852,login=login,removeDuplicatedTimestamps=TRUE))
hmadultsData <- as.data.frame(getMovebankData(study=69402287,login=login,removeDuplicatedTimestamps=TRUE))
hrjuvData <- as.data.frame(getMovebankData(study=76429224,login=login,removeDuplicatedTimestamps=TRUE))
hrjuv2016data <- as.data.frame(getMovebankData(study=175328223,login=login,removeDuplicatedTimestamps=TRUE))
hrjuv2017data <- as.data.frame(getMovebankData(study=293970900,login=login,removeDuplicatedTimestamps=TRUE))
HQXS_Black_tailed_godwitsData <- as.data.frame(getMovebankData(study=1658294759,login=login,removeDuplicatedTimestamps=TRUE))
iberiaBlackwitsData <- as.data.frame(getMovebankData(study= 49547785,login=login,removeDuplicatedTimestamps=TRUE))
icarusData <- as.data.frame(getMovebankData(study=1487044886,login=login,removeDuplicatedTimestamps=TRUE))
polishData <- as.data.frame(getMovebankData(study=163516781,login=login,removeDuplicatedTimestamps=TRUE))
teamPiersmaHQXSdata <- as.data.frame(getMovebankData(study=1563249841,login=login,removeDuplicatedTimestamps=TRUE))
teamPiersmaHQXS2022data <- as.data.frame(getMovebankData(study=2083443328,login=login,removeDuplicatedTimestamps=TRUE))
wildjuvData <- as.data.frame(getMovebankData(study=75360602,login=login,removeDuplicatedTimestamps=TRUE))
wildjuv2016data <- as.data.frame(getMovebankData(study=170829089,login=login,removeDuplicatedTimestamps=TRUE))
wildjuv2017data <- as.data.frame(getMovebankData(study=282596404,login=login,removeDuplicatedTimestamps=TRUE))

library(dplyr)
library(plyr)
combinedReferenceData <- rbind.fill(ib19,microwave2021,extremadura2022,southholland2021,BtgTagus2021,ad_dum2018,ad_dum2019,ad_dum2020,ad_dum2021,ch_dum2018,ch_dum2019,ch_dum2020,ch_dum2021,ad_dum2022,ch_dum2022,hmadults,hrjuv,hrjuv2016,
                                    hrjuv2017,polish,HQXS_Black_tailed_godwits,iberiaBlackwits,icarus,teamPiersmaHQXS,teamPiersmaHQXS2022,wildjuv,wildjuv2016,wildjuv2017)
combinedLocationData <- rbind.fill(ib19data,microwave2021data,extremadura2022data,southholland2021data,BtgTagus2021data,ad_dum2018data,ad_dum2019data,ad_dum2020data,ad_dum2021data,ch_dum2018data,ch_dum2019data,ch_dum2020data,ch_dum2021data,ad_dum2022data,ch_dum2022data,hmadultsData,hrjuvData,hrjuv2016data,
                                             hrjuv2017data,polishData,HQXS_Black_tailed_godwitsData,iberiaBlackwitsData,icarusData,teamPiersmaHQXSdata,teamPiersmaHQXS2022data,wildjuvData,wildjuv2016data,wildjuv2017data)

######################add study name column######################
combinedReferenceData$study_name <- combinedReferenceData$study_id
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 652989041] <- "black-tailed godwit iberia - NL 2019"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 69402287] <- "Haanmeer Adults"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 163516781] <- "Polish Adults"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1498143083] <- "black-tailed godwit nl 2021 microwave"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1145538280] <- "black-tailed godwits nl south holland"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 76429224] <- "hand-raised juveniles"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 175328223] <- "hand-raised juveniles 2016"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 293970900] <- "hand-raised juveniles 2017"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 49547785] <- "Iberia Blackwits"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 75360602] <- "wild juveniles"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 170829089] <- "wild juveniles 2"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 282596404] <- "wild juveniles 3"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 484019425] <- "dummer adults 2018"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 831990025] <- "dummer adults 2019"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1105026166] <- "dummer adults 2020"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1482506572] <- "dummer adults 2021"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1751337831] <- "dummer adults 2022"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 500187586] <- "dummer chicks 2018"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 878914763] <- "dummer chicks 2019"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1183466126] <- "dummer chicks 2020"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1482505185] <- "dummer chicks 2021"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 2098519852] <- "dummer chicks 2022"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1563249841] <- "teampiersma hqxs"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 2083443328] <- "teampiersma hqxs 2022"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1693518103] <- "BtgTagus2021"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1658294759] <- "HQXS Black-tailed Godwits"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1487044886] <- "icarus black-tailed godwits theunis piersma"
combinedReferenceData["study_name"][combinedReferenceData["study_name"] == 1923591036] <- "Black-tailed Godwits Extremadura (Lotek)"

############################combine reference and location data############################
library('data.table') # for setnames function
setnames(combinedReferenceData, "animal_local_identifier", "local_identifier") #change name to match with location data table
allLocations <-  merge(combinedReferenceData,combinedLocationData,by="local_identifier",all=FALSE)
############################filter############################
allLocations <- subset(allLocations, animal_taxon_detail == "Limosa limosa limosa") #remove bar-tailed/icelandic
allLocations <- subset(allLocations, argos_lc!="A" & argos_lc!="B" & argos_lc!="C" & argos_lc!="Z" | is.na(argos_lc)) #remove low quality locations
allLocations <- allLocations[,colSums(is.na(allLocations))<nrow(allLocations)] #remove NA columns
allLocations <- distinct(allLocations)#remove duplicate rows
allLocations <- subset(allLocations,select = c(local_identifier,study_name,timestamp,location_lat,location_long,tag_manufacturer_name,study_site,argos_lc,sensor))
############################export############################################################################################################################
write.csv(allLocations,"allLocations.csv", row.names = FALSE)
