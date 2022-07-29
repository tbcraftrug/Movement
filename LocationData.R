######################Combined godwit location data##########################################################################################################
library(move) 
#https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
library(tidyverse)

#Movebank login
username = "t.b.craft"
password = "GodwitSnl24!!"
login<-movebankLogin(username,password)

############################Microwave######################################################################################################################
#Microwave reference data
ib19 <- as.data.frame(getMovebankReferenceTable(study=652989041,login=login))
hmadults <- as.data.frame(getMovebankReferenceTable(study=69402287,login=login))
polish <- as.data.frame(getMovebankReferenceTable(study=163516781,login=login))
microwave2021 <- as.data.frame(getMovebankReferenceTable(study=1498143083,login=login))
southholland2021microwave <- as.data.frame(getMovebankReferenceTable(study=1145538280,login=login))
hrjuv <- as.data.frame(getMovebankReferenceTable(study=76429224,login=login))
hrjuv2016 <- as.data.frame(getMovebankReferenceTable(study=175328223,login=login))
hrjuv2017 <- as.data.frame(getMovebankReferenceTable(study= 293970900,login=login))
iberiaBlackwits <- as.data.frame(getMovebankReferenceTable(study= 49547785,login=login))
wildjuv <- as.data.frame(getMovebankReferenceTable(study=75360602,login=login))
wildjuv2016 <- as.data.frame(getMovebankReferenceTable(study=170829089,login=login))
wildjuv2017 <- as.data.frame(getMovebankReferenceTable(study=282596404,login=login))
ad_dum2018 <- as.data.frame(getMovebankReferenceTable(study=484019425,login=login))
ad_dum2019 <- as.data.frame(getMovebankReferenceTable(study=831990025,login=login))
ad_dum2020microwave <- as.data.frame(getMovebankReferenceTable(study=1105026166,login=login))
ch_dum2018 <- as.data.frame(getMovebankReferenceTable(study=500187586,login=login))
ch_dum2019 <- as.data.frame(getMovebankReferenceTable(study=878914763,login=login))
ch_dum2020 <- as.data.frame(getMovebankReferenceTable(study=1183466126,login=login))
ch_dum2021microwave <- as.data.frame(getMovebankReferenceTable(study=1482505185,login=login))

#Microwave location data
ib19data <- as.data.frame(getMovebankData(study=652989041,login=login, removeDuplicatedTimestamps=T))
hmadultsData <- as.data.frame(getMovebankData(study=69402287,login=login, removeDuplicatedTimestamps=T))
polishData <- as.data.frame(getMovebankData(study=163516781,login=login, removeDuplicatedTimestamps=T))
microwave2021data <- as.data.frame(getMovebankData(study=1498143083,login=login, removeDuplicatedTimestamps=T))
southholland2021microwaveData <- as.data.frame(getMovebankData(study=1145538280,login=login, removeDuplicatedTimestamps=T))
hrjuvData <- as.data.frame(getMovebankData(study=76429224,login=login, removeDuplicatedTimestamps=T))
hrjuv2016data <- as.data.frame(getMovebankData(study=175328223,login=login, removeDuplicatedTimestamps=T))
hrjuv2017data <- as.data.frame(getMovebankData(study= 293970900,login=login, removeDuplicatedTimestamps=T))
iberiaBlackwitsData <- as.data.frame(getMovebankData(study=49547785,login=login, removeDuplicatedTimestamps=T))
wildjuvData <- as.data.frame(getMovebankData(study=75360602,login=login, removeDuplicatedTimestamps=T))
wildjuv2016data <- as.data.frame(getMovebankData(study=170829089,login=login, removeDuplicatedTimestamps=T))
wildjuv2017data <- as.data.frame(getMovebankData(study=282596404,login=login, removeDuplicatedTimestamps=T))
ad_dum2018data <- as.data.frame(getMovebankData(study=484019425,login=login, removeDuplicatedTimestamps=T))
ad_dum2019data <- as.data.frame(getMovebankData(study=831990025,login=login, removeDuplicatedTimestamps=T))
ad_dum2020microwaveData <- as.data.frame(getMovebankData(study=1105026166,login=login, removeDuplicatedTimestamps=T))
ch_dum2018data <- as.data.frame(getMovebankData(study=500187586,login=login, removeDuplicatedTimestamps=T))
ch_dum2019data <- as.data.frame(getMovebankData(study=878914763,login=login, removeDuplicatedTimestamps=T))
ch_dum2020data <- as.data.frame(getMovebankData(study=1183466126,login=login, removeDuplicatedTimestamps=T))
ch_dum2021microwaveData <- as.data.frame(getMovebankData(study=1482505185,login=login, removeDuplicatedTimestamps=T))

#combine Microwave reference data
library(dplyr)
library(plyr) # for fill function
combinedMicrowave <- rbind.fill(ib19,hmadults,microwave2021,southholland2021microwave,hrjuv,
                                hrjuv2016,hrjuv2017,iberiaBlackwits, wildjuv,wildjuv2016,wildjuv2017,ad_dum2018,
                                ad_dum2019,ad_dum2020microwave, ch_dum2018,ch_dum2019,ch_dum2020,ch_dum2021microwave,polish)

combinedMicrowave <- combinedMicrowave[!duplicated(combinedMicrowave$animal_local_identifier), ]#remove duplicate rows
combinedMicrowave <- subset(combinedMicrowave, animal_taxon_detail == "Limosa limosa limosa") #remove bar-tailed/icelandic
combinedMicrowave <- subset(combinedMicrowave, tag_manufacturer_name == "Microwave") #remove lotek/HQXS

#add study name column
combinedMicrowave$study_name <- combinedMicrowave$study_id
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 652989041] <- "black-tailed godwit iberia - NL 2019"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 69402287] <- "Haanmeer Adults"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 163516781] <- "Polish Adults"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 1498143083] <- "black-tailed godwit nl 2021 microwave"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 1145538280] <- "black-tailed godwits nl south holland"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 76429224] <- "hand-raised juveniles"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 175328223] <- "hand-raised juveniles 2016"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 293970900] <- "hand-raised juveniles 2017"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 49547785] <- "Iberia Blackwits"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 75360602] <- "wild juveniles"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 170829089] <- "wild juveniles 2"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 282596404] <- "wild juveniles 3"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 484019425] <- "dummer adults 2018"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 831990025] <- "dummer adults 2019"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 1105026166] <- "dummer adults 2020"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 500187586] <- "dummer chicks 2018"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 878914763] <- "dummer chicks 2019"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 1183466126] <- "dummer chicks 2020"
combinedMicrowave["study_name"][combinedMicrowave["study_name"] == 1482505185] <- "dummer chicks 2021"

#combine Microwave location data
combinedMicrowaveData <- rbind.fill(ib19data,hmadultsData,microwave2021data,southholland2021microwaveData,hrjuvData,
                                    hrjuv2016data,hrjuv2017data,iberiaBlackwitsData,wildjuvData,wildjuv2016data,wildjuv2017data,ad_dum2018data,
                                    ad_dum2019data,ad_dum2020microwaveData,ch_dum2018data,ch_dum2019data,ch_dum2020data,ch_dum2021microwaveData,polishData)
#remove HQXS/Lotek/Icarus tags
combinedMicrowaveData <- subset(combinedMicrowaveData,sensor == "Argos Doppler Shift")

#remove islandica/lapponica
combinedMicrowaveData <- subset(combinedMicrowaveData, trackId != "Iwik2" &
                                  trackId !="Iwik1" &
                                  trackId !="Iwik3" &
                                  trackId != "Iwik4" &
                                  trackId !="C1WPPW - NL" &
                                  trackId !="Montenovo" &
                                  trackId !="Kaldadarnes" &
                                  trackId !="Vorsaber")

#filter low quality locations
combinedMicrowaveData <- subset(combinedMicrowaveData, argos_lc=="3"|argos_lc=="2"|argos_lc=="1")

#remove NA columns
combinedMicrowaveData <- combinedMicrowaveData[,colSums(is.na(combinedMicrowaveData))<nrow(combinedMicrowaveData)]


############################HQXS##########################################################################################################################################################
#HQXS reference data
teamPiersmaHQXS <- as.data.frame(getMovebankReferenceTable(study=1563249841,login=login))
teamPiersmaHQXS2022 <- as.data.frame(getMovebankReferenceTable(study=2083443328,login=login))
BtgTagus2021 <- as.data.frame(getMovebankReferenceTable(study=1693518103,login=login))
HQXS_Black_tailed_godwits <- as.data.frame(getMovebankReferenceTable(study=1658294759,login=login))
southholland2021HQXS <- as.data.frame(getMovebankReferenceTable(study=1145538280,login=login))
ch_dum2021HQXS <- as.data.frame(getMovebankReferenceTable(study=1482505185,login=login))

#HQXS location data
teamPiersmaHQXSdata <- as.data.frame(getMovebankData(study=1563249841,login=login, removeDuplicatedTimestamps=T))
teamPiersmaHQXS2022data <- as.data.frame(getMovebankData(study=2083443328,login=login, removeDuplicatedTimestamps=T))
BtgTagus2021data <- as.data.frame(getMovebankData(study=1693518103, login=login, removeDuplicatedTimestamps=T))
HQXS_Black_tailed_godwitsData <- as.data.frame(getMovebankData(study=1658294759, login=login, removeDuplicatedTimestamps=T))
southholland2021HQXSdata <- as.data.frame(getMovebankData(study=1145538280,login=login, removeDuplicatedTimestamps=T))
ch_dum2021HQXSdata <- as.data.frame(getMovebankData(study=1482505185,login=login, removeDuplicatedTimestamps=T))

#combine HQXS reference data
combinedHQXS <- rbind.fill(teamPiersmaHQXS,teamPiersmaHQXS2022,BtgTagus2021, HQXS_Black_tailed_godwits,southholland2021HQXS,ch_dum2021HQXS)
combinedHQXS <- subset(combinedHQXS, animal_taxon_detail == "Limosa limosa limosa") #remove bar-tailed/icelandic
combinedHQXS <- subset(combinedHQXS, tag_manufacturer_name == "HQXS") #remove lotek/microwave/icarus
combinedHQXS <- combinedHQXS[!duplicated(combinedHQXS$animal_local_identifier), ]#remove duplicate rows

#add study name column
combinedHQXS$study_name <- combinedHQXS$study_id
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1563249841] <- "teampiersma hqxs"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 2083443328] <- "teampiersma hqxs 2022"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1693518103] <- "BtgTagus2021"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1658294759] <- "HQXS Black-tailed Godwits"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1145538280] <- "black-tailed godwits nl south holland"
combinedHQXS["study_name"][combinedHQXS["study_name"] == 1482505185] <- "dummer chicks 2021"

#combined HQXS location data
combinedHQXSdata <- rbind.fill(teamPiersmaHQXSdata,teamPiersmaHQXS2022data,BtgTagus2021data,HQXS_Black_tailed_godwitsData,
                               southholland2021HQXSdata,ch_dum2021HQXSdata)
combinedHQXSdata <- subset(combinedHQXSdata, trackId != "Mouchao") #remove islandica
combinedHQXSdata <- subset(combinedHQXSdata, sensor_type == "GPS") #remove Argos locations
combinedHQXSdata <- combinedHQXSdata[,colSums(is.na(combinedHQXSdata))<nrow(combinedHQXSdata)] #remove NA columns


############################ICARUS##########################################################################################################################################################
#Icarus reference data
icarus <- as.data.frame(getMovebankReferenceTable(study=1487044886,login=login))
icarus <- icarus[!duplicated(icarus$animal_local_identifier), ] #remove duplicate rows
icarus <- subset(icarus, animal_taxon_detail == "Limosa limosa limosa") #remove bar-tailed/icelandic

#add study name column
icarus$study_name <- icarus$study_id
icarus["study_name"][icarus["study_name"] == 1487044886] <- "icarus black-tailed godwits theunis piersma"

#Icarus location data
icarusData <-  as.data.frame(getMovebankData(study=1487044886,login=login, removeDuplicatedTimestamps=T))
icarusData <- subset(icarusData, trackId != "Benavente") #remove islandica
icarusData <- icarusData[,colSums(is.na(icarusData))<nrow(icarusData)] #remove NA columns

############################Lotek##########################################################################################################################################################
#Lotek reference data
ad_dum2020lotek <- as.data.frame(getMovebankReferenceTable(study=1105026166,login=login))
ad_dum2021 <- as.data.frame(getMovebankReferenceTable(study=1482506572,login=login)) #all tags are LOTEK in ad_dum2021
ad_dum2021 <- ad_dum2021[!duplicated(ad_dum2021$animal_local_identifier), ]#remove duplicate rows
extremadura2022 <- as.data.frame(getMovebankReferenceTable(study=1923591036,login=login))
extremadura2022 <- extremadura2022[!duplicated(extremadura2022$animal_local_identifier), ] #remove duplicate rows

#Lotek location data
ad_dum2020lotekData <- as.data.frame(getMovebankData(study=1105026166,login=login, removeDuplicatedTimestamps=T))
ad_dum2021data <- as.data.frame(getMovebankData(study=1482506572,login=login, removeDuplicatedTimestamps=T))
extremadura2022data <- as.data.frame(getMovebankData(study=1923591036,login=login, removeDuplicatedTimestamps=T))

#combine Lotek reference data
combinedLotek <- rbind.fill(ad_dum2020lotek,ad_dum2021,extremadura2022)
combinedLotek <- combinedLotek[!duplicated(combinedLotek$animal_local_identifier), ]#remove duplicate rows
combinedLotek <- subset(combinedLotek, animal_taxon_detail == "Limosa limosa limosa") #remove bar-tailed/icelandic
combinedLotek <- subset(combinedLotek, tag_manufacturer_name == "Lotek") #remove hqxs/microwave/icarus

#add study name column
combinedLotek$study_name <- combinedLotek$study_id
combinedLotek["study_name"][combinedLotek["study_name"] == 1105026166] <- "dummer adults 2020"
combinedLotek["study_name"][combinedLotek["study_name"] == 1482506572] <- "dummer adults 2021"
combinedLotek["study_name"][combinedLotek["study_name"] == 1923591036] <- "Black-tailed Godwits Extremadura (Lotek)"

#combine Lotek location data 
combinedLotekData <- rbind.fill(ad_dum2020lotekData,ad_dum2021data, extremadura2022)
combinedLotekData <- subset(combinedLotekData, trackId != "WWN-WY(V)") #remove islandica
combinedLotekData <- combinedLotekData[,colSums(is.na(combinedLotekData))<nrow(combinedLotekData)] #remove NA columns


############################Combine############################################################################################################################
#combine all reference data
combinedReferenceData <- rbind.fill(combinedMicrowave,combinedHQXS,icarus,combinedLotek)

#combine all location data
combinedLocationData <- rbind.fill(combinedMicrowaveData,combinedHQXSdata,icarusData,combinedLotekData)
combinedLocationData$timestamp <- as.character(combinedLocationData$timestamp) #convert timestamps
combinedLocationData$timestamps <- as.character(combinedLocationData$timestamps) #convert timestamps
combinedLocationData <- combinedLocationData %>% mutate_all(na_if,"") #convert blank columns to NA
combinedLocationData <- combinedLocationData[,colSums(is.na(combinedLocationData))<nrow(combinedLocationData)] #remove NA columns
combinedLocationData$timestamp <- as.Date(combinedLocationData$timestamp) #reconvert timestamps
combinedLocationData$timestamps <- as.Date(combinedLocationData$timestamps) #reconvert timestamps

#combine reference and location data
library('data.table') # for setnames function
setnames(combinedReferenceData, "animal_local_identifier", "local_identifier") #change name to match with location data table
allLocations <-  merge(combinedReferenceData,combinedLocationData,by="local_identifier",all=FALSE)
allLocations <- subset(allLocations,select = c(local_identifier,study_name,timestamp,location_lat,location_long))
                      
############################Export######################################################################################################################
# #All data
#write.csv(allLocations,"allLocations.csv", row.names = FALSE)
#
# #All data since xxxx-xx-xx
# allLocationsSince <- allLocations%>%filter(timestamp>"2022-04-15" & timestamp < "2022-05-10")%>%
#   write.csv("allLocationsSpring2022.csv", row.names = FALSE)
# allLocationsPreBreeding <- allLocations%>%filter(timestamp>"2022-01-01" & timestamp < "2022-03-01")%>%
#   write.csv("allLocationsSpring2022.csv", row.names = FALSE)
# 
# #last timestamps of all tags
# allLastTimestamps <- allLocations %>% 
#    group_by(local_identifier) %>%
#    slice(which.max(timestamp))
# 
# # #Last timestamps since xxxx-xx-xx
# allLocationsSince <- allLastTimestamps%>%filter(timestamp>"2022-04-01" & timestamp < "2022-05-10")%>%
#   write.csv("allLocationsSpring2022.csv", row.names = FALSE)
# allLocationsPreBreeding <- allLastTimestamps%>%filter(timestamp>"2022-01-01" & timestamp < "2022-03-01")%>%
#   write.csv("allLocationsSpring2022.csv", row.names = FALSE)

# 
# #Dutch birds only
# locations%>%filter(!study_name %like% "dummer",)%>%
#   write.csv("DutchBirds.csv", row.names = FALSE)
# 
# #German birds only
# locations%>%filter(study_name %like% "dummer",)%>%
#   write.csv("GermanBirds.csv", row.names = FALSE)
# 
# #Microwave only
# locations%>%filter(tag_manufacturer_name == "Microwave")%>%
#   write.csv("MicrowaveBirds.csv", row.names = FALSE)
# 
# #HQXS only
# locations%>%filter(tag_manufacturer_name == "HQXS")%>%
#   write.csv("HQXSbirds.csv", row.names = FALSE)
# 
# #Icarus only
# locations%>%filter(tag_manufacturer_name == "Icarus")%>%
#   write.csv("IcarusBirds.csv", row.names = FALSE)
# 
# #Lotek only
# locations%>%filter(tag_manufacturer_name == "Lotek")%>%
#   write.csv("/LotekBirds.csv", row.names = FALSE)
# 
################################post-fledge period############################################################
# locations%>%filter(study_name == "wild juveniles")%>%
#   filter(timestamp > "2015-04-01" & timestamp < "2015-08-01")%>%
#   write.csv("postFledge2015.csv", row.names = FALSE)
# 
# locations%>%filter(study_name == "wild juveniles 2")%>%
#   filter(timestamp > "2016-04-01" & timestamp < "2016-08-01")%>%
#   write.csv("/postFledge2016.csv", row.names = FALSE)
# 
# locations%>%filter(study_name == "wild juveniles 3")%>%
#   filter(timestamp > "2017-04-01" & timestamp < "2017-08-01")%>%
#   write.csv("/postFledge2017.csv", row.names = FALSE)
# 
# 
##############################################kmz############################################################################################
#   # create points layer and sliding time bar
#   library("GISTools")
#   library("spacetime")
#   
#   # make it spatial
#   sp <-SpatialPoints(locations[,c("location_long","location_lat")])
#   
#   # add the projection WG84
#   proj4string(sp) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#   locations$timestamp <- as.POSIXct(locations$timestamp)
#   
#   # define time element
#   sp.st<-STIDF(sp,time=timestamp,data=locations[,c("local_identifier","timestamp","animal_timestamp_start","animal_timestamp_end","study_name")])
#   sp.st<-sp.st[order(sp.st$local_identifier, sp.st$timestamp),] # order by id and timestamp
#   
#   # creating lines
#   locations2<-locations[,c("location_long","location_lat","local_identifier")]
#   locations3<-locations2[order(locations2$local_identifier, locations$timestamp),] # order by id and timestamp
#   coordinates(locations3) <- ~location_long+location_lat
#   ## list of Lines per id, each with one Line in a list
#   x <- lapply(split(locations3, locations3$local_identifier), function(x) Lines(list(Line(coordinates(x))), x$local_identifier[1L]))
#   lines <- SpatialLines(x)
#   data <- data.frame(id = unique(locations3$local_identifier))
#   rownames(data) <- data$id
#   l <- SpatialLinesDataFrame(lines, data)
#   # add the projection WG84
#   proj4string(l) <- CRS("+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#   
#   todaydate <- Sys.Date()
#   todaydate
#   
#   #plot as kmz
#   library("plotKML")
#   shape="http://maps.google.com/mapfiles/kml/pal2/icon18.png"
#   kml_open(paste("allbtg",todaydate,".kmz", sep=""))
#   # kml_open("query.kmz")
#   kml_layer(lapply(split(sp.st,sp.st$local_identifier.x), function(kmz.pts) kml_layer(kmz.pts,shape=shape,colour=local_identifier.x, points_names=kmz.pts$argos_lc,size=0.5,balloon=TRUE,subfolder.name=unique(kmz.pts$local_identifier.x))))
#   kml_layer(l,colour=l$id, lwd=2)
#   kml_close(paste("allbtg",todaydate,".kmz", sep=""))
#   # kml_close("query.kmz")


