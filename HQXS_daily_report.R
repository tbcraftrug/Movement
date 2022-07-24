library('move')
library('tidyverse')
library('plyr')

#Movebank login Dutch godwits
username = "t.b.craft"
password = "GodwitSnl24!!"
login<-movebankLogin(username,password)

#HQXS reference data
teamPiersmaHQXS <- as.data.frame(getMovebankReferenceTable(study=1563249841,login=login))
teamPiersmaHQXS2022 <- as.data.frame(getMovebankReferenceTable(study=	2083443328,login=login))
BtgTagus2021 <- as.data.frame(getMovebankReferenceTable(study=1693518103, login=login))
HQXS_Black_tailed_godwits <- as.data.frame(getMovebankReferenceTable(study=1658294759, login=login))
southholland2021 <- as.data.frame(getMovebankReferenceTable(study=1145538280,login=login))
southholland2021HQXS <- subset(southholland2021,tag_manufacturer_name == "HQXS" ) #remove microwave tags
#combine HQXS reference data
combinedHQXS <- rbind.fill(teamPiersmaHQXS,teamPiersmaHQXS2022,BtgTagus2021, HQXS_Black_tailed_godwits,southholland2021HQXS)
combinedHQXS <- combinedHQXS[!duplicated(combinedHQXS$animal_local_identifier), ]#remove duplicate rows
#add study name column
combinedHQXS$study_name <- combinedHQXS$study_id
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1563249841] <- "teampiersma hqxs"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 	2083443328] <- "teampiersma hqxs 2022"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1693518103] <- "BtgTagus2021"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1658294759] <- "HQXS Black-tailed Godwits"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1145538280] <- "black-tailed godwits nl south holland"

#HQXS location data
teamPiersmaHQXSlocations <- as.data.frame(getMovebankData(study=1563249841,login=login, removeDuplicatedTimestamps=T))
teamPiersmaHQXS2022locations <- as.data.frame(getMovebankData(study=	2083443328,login=login, removeDuplicatedTimestamps=T))
BtgTagus2021locations <- as.data.frame(getMovebankData(study=1693518103, login=login, removeDuplicatedTimestamps=T))
HQXS_Black_tailed_godwitsLocations <- as.data.frame(getMovebankData(study=1658294759, login=login, removeDuplicatedTimestamps=T))
southholland2021locations <- as.data.frame(getMovebankData(study=1145538280,login=login, removeDuplicatedTimestamps=T))
southholland2021HQXSlocations <- subset(southholland2021locations,sensor_type == "GPS")
#combined HQXS location data
combinedHQXSlocations <- rbind.fill(teamPiersmaHQXSlocations,teamPiersmaHQXS2022locations,BtgTagus2021locations,HQXS_Black_tailed_godwitsLocations,
                               southholland2021HQXSlocations)

#combine reference and location table
library('data.table') # for setnames function
setnames(combinedHQXS, "animal_local_identifier", "local_identifier") #change name to match with location data table
combinedHQXSdata <- merge(combinedHQXS,combinedHQXSlocations,by="local_identifier",all=FALSE)
combinedHQXSdata$timestamp <- as.POSIXct(combinedHQXSdata$timestamp)
combinedHQXSdata$animal_timestamp_end <- as.POSIXct(combinedHQXSdata$animal_timestamp_end)

#last timestamps of all HQXS tags
HQXSlastTimestamps <- combinedHQXSdata %>% 
  group_by(local_identifier) %>%
  slice(which.max(timestamp))

#remove unnecessary columns
HQXSlastTimestamps <- HQXSlastTimestamps[,c("local_identifier", "tag_serial_no", "animal_ring_id","study_name","tag_voltage", "external_temperature", 
                                                                 "timestamp", "location_lat", "location_long")]

#active in last 30 days
library(lubridate)
HQXSactive30days <- HQXSlastTimestamps %>%
  filter(timestamp  >= today() - days(30))
#write.csv(HQXSactive30days,"combinedHQXSactive30days.csv", row.names = FALSE)

#active in last 3 days
HQXSactive3days <- HQXSlastTimestamps %>%
  filter(timestamp  >= today() - days(3))
#write.csv(HQXSactive3days,"combinedHQXSactive3days.csv", row.names = FALSE)

#voltage <3.65 (active in past 30 days)
reduceFixes <- filter(HQXSactive30days, tag_voltage < 3650)
#write.csv(reduceFixes,"reduceFixes.csv", row.names = FALSE)

#inactive for over 60 days
inactive60daysHQXS <- HQXSlastTimestamps %>%
  filter(timestamp  <= today() - days(60)) %>%
  filter(local_identifier != "Gellehuister" & 
           local_identifier != "Swarte_Walde" & 
           local_identifier != "Sudhoeke")
write.csv(inactive60daysHQXS,"inactive60daysHQXS.csv", row.names = FALSE)

#active in last 60 days
active60daysHQXS <- HQXSlastTimestamps %>%
  filter(timestamp  >= today() - days(60))
write.csv(active60daysHQXS,"active60daysHQXS.csv", row.names = FALSE)

#fixes per day Teampiersma HQXS 2022
teamPiersmaHQXS2022$ActiveDays <- difftime(as.Date(teamPiersmaHQXS2022$animal_timestamp_end),as.Date(teamPiersmaHQXS2022$animal_timestamp_start), units = c("days"))
teamPiersmaHQXS2022$fixesPerDay <- teamPiersmaHQXS2022$number_of_location_events/as.numeric(teamPiersmaHQXS2022$ActiveDays)

#fixes per day combined
combinedHQXS$ActiveDays <- difftime(as.Date(combinedHQXS$animal_timestamp_end),as.Date(combinedHQXS$animal_timestamp_start), units = c("days"))
combinedHQXS$fixesPerDay <- combinedHQXS$number_of_location_events/as.numeric(combinedHQXS$ActiveDays)
