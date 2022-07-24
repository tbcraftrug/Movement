######################Combined godwit location data##########################################################################################################
library('move')
library('tidyverse')

#Movebank login Dutch godwits
username = "t.b.craft"
password = "GodwitSnl24!!"
login<-movebankLogin(username,password)

#Movebank login German godwits
username2 = "Dummer_Life"
password2 = "Uferschnepfe_2018"
login2<-movebankLogin(username2,password2)

############################Microwave######################################################################################################################
#Microwave reference data
ib19 <- as.data.frame(getMovebankReferenceTable(study=652989041,login=login))
ib19 <- subset(ib19, animal_local_identifier != "Sarilhos")
hmadults <- as.data.frame(getMovebankReferenceTable(study=69402287,login=login))
polish <- as.data.frame(getMovebankReferenceTable(study=163516781,login=login))
microwave2021 <- as.data.frame(getMovebankReferenceTable(study=1498143083,login=login))
southholland2021 <- as.data.frame(getMovebankReferenceTable(study=1145538280,login=login))
southholland2021microwave <- subset(southholland2021,tag_manufacturer_name == "Microwave") #remove 5 HQXS tags
hrjuv <- as.data.frame(getMovebankReferenceTable(study=76429224,login=login))
hrjuv2016 <- as.data.frame(getMovebankReferenceTable(study=175328223,login=login))
hrjuv2017 <- as.data.frame(getMovebankReferenceTable(study= 293970900,login=login))
iberiaBlackwits <- as.data.frame(getMovebankReferenceTable(study= 49547785,login=login))
wildjuv <- as.data.frame(getMovebankReferenceTable(study=75360602,login=login))
wildjuv2016 <- as.data.frame(getMovebankReferenceTable(study=170829089,login=login))
wildjuv2017 <- as.data.frame(getMovebankReferenceTable(study=282596404,login=login))
ad_dum2018 <- as.data.frame(getMovebankReferenceTable(study=484019425,login=login2))
ad_dum2019 <- as.data.frame(getMovebankReferenceTable(study=831990025,login=login2))
ad_dum2020 <- as.data.frame(getMovebankReferenceTable(study=1105026166,login=login2))
ad_dum2020microwave <- subset(ad_dum2020,tag_manufacturer_name == "Microwave" ) #remove 9 Lotec tags
ch_dum2018 <- as.data.frame(getMovebankReferenceTable(study=500187586,login=login2))
ch_dum2019 <- as.data.frame(getMovebankReferenceTable(study=878914763,login=login2))
ch_dum2020 <- as.data.frame(getMovebankReferenceTable(study=1183466126,login=login2))
ch_dum2021 <- as.data.frame(getMovebankReferenceTable(study=1482505185,login=login2))
ch_dum2021microwave <- subset(ch_dum2021,tag_manufacturer_name == "Microwave" ) #remove 1 HQXS tag

#Microwave location data
ib19data <- as.data.frame(getMovebankData(study=652989041,login=login, removeDuplicatedTimestamps=T))
ib19data <- filter(ib19data, trackId != "Sarilhos") #remove islandicas
hmadultsData <- as.data.frame(getMovebankData(study=69402287,login=login, removeDuplicatedTimestamps=T))
polishData <- as.data.frame(getMovebankData(study=163516781,login=login, removeDuplicatedTimestamps=T))
polishData$study_name <- "polish adults"
microwave2021data <- as.data.frame(getMovebankData(study=1498143083,login=login, removeDuplicatedTimestamps=T))
southholland2021data <- as.data.frame(getMovebankData(study=1145538280,login=login, removeDuplicatedTimestamps=T))
southholland2021microwaveData <- subset(southholland2021data,sensor == "Argos Doppler Shift") #remove HQXS tag
hrjuvData <- as.data.frame(getMovebankData(study=76429224,login=login, removeDuplicatedTimestamps=T))
hrjuvData <- filter(hrjuvData, individual_id != "80444306"&individual_id != "80444217"&individual_id != "90924588"&individual_id !=	"90944036") # remove bar-tailed
hrjuv2016data <- as.data.frame(getMovebankData(study=175328223,login=login, removeDuplicatedTimestamps=T))
hrjuv2017data <- as.data.frame(getMovebankData(study= 293970900,login=login, removeDuplicatedTimestamps=T))
hrjuv2017data <- filter(hrjuv2017data, individual_id != "293973960") # remove bar-tailed
iberiaBlackwitsData <- as.data.frame(getMovebankData(study=49547785,login=login, removeDuplicatedTimestamps=T))
iberiaBlackwitsData <- filter(iberiaBlackwitsData, individual_id != "54790555"&individual_id !="79307795"&individual_id !="82322234") #remove islandicas
wildjuvData <- as.data.frame(getMovebankData(study=75360602,login=login, removeDuplicatedTimestamps=T))
wildjuv2016data <- as.data.frame(getMovebankData(study=170829089,login=login, removeDuplicatedTimestamps=T))
wildjuv2017data <- as.data.frame(getMovebankData(study=282596404,login=login, removeDuplicatedTimestamps=T))
ad_dum2018data <- as.data.frame(getMovebankData(study=484019425,login=login2, removeDuplicatedTimestamps=T))
ad_dum2019data <- as.data.frame(getMovebankData(study=831990025,login=login2, removeDuplicatedTimestamps=T))
ad_dum2020data <- as.data.frame(getMovebankData(study=1105026166,login=login2, removeDuplicatedTimestamps=T))
ad_dum2020microwaveData <- subset(ad_dum2020data,sensor == "Argos Doppler Shift" ) #remove Lotec tags
ch_dum2018data <- as.data.frame(getMovebankData(study=500187586,login=login2, removeDuplicatedTimestamps=T))
ch_dum2019data <- as.data.frame(getMovebankData(study=878914763,login=login2, removeDuplicatedTimestamps=T))
ch_dum2020data <- as.data.frame(getMovebankData(study=1183466126,login=login2, removeDuplicatedTimestamps=T))
ch_dum2021data <- as.data.frame(getMovebankData(study=1482505185,login=login2, removeDuplicatedTimestamps=T))
ch_dum2021microwaveData <- subset(ch_dum2021data,sensor == "Argos Doppler Shift" ) #remove HQXS tag

#combine Microwave reference data
library('plyr') # for .fill function
combinedMicrowave <- rbind.fill(ib19,hmadults,microwave2021,southholland2021microwave,hrjuv,
                                hrjuv2016,hrjuv2017,iberiaBlackwits, wildjuv,wildjuv2016,wildjuv2017,ad_dum2018,
                                ad_dum2019,ad_dum2020microwave, ch_dum2018,ch_dum2019,ch_dum2020,ch_dum2021microwave,polish)

combinedMicrowave <- combinedMicrowave[!duplicated(combinedMicrowave$animal_local_identifier), ]#remove duplicate rows

combinedMicrowave <- subset(combinedMicrowave, animal_taxon_detail == "Limosa limosa limosa") #remove bar-tailed/icelandic

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
                                    ad_dum2019data,ad_dum2020microwaveData, ch_dum2018data,ch_dum2019data,ch_dum2020data,ch_dum2021microwaveData, polishData)

#filter low quality locations
combinedMicrowaveData <- subset(combinedMicrowaveData, argos_lc=="3"|argos_lc=="2"|argos_lc=="1")

############################HQXS##########################################################################################################################################################
#HQXS reference data
teamPiersmaHQXS <- as.data.frame(getMovebankReferenceTable(study=1563249841,login=login))
teamPiersmaHQXS2022 <- as.data.frame(getMovebankReferenceTable(study=2083443328,login=login))
BtgTagus2021 <- as.data.frame(getMovebankReferenceTable(study=1693518103, login=login))
BtgTagus2021 <- filter(BtgTagus2021, animal_local_identifier != "Mouchao") #remove icelandic
HQXS_Black_tailed_godwits <- as.data.frame(getMovebankReferenceTable(study=1658294759, login=login))
southholland2021HQXS <- subset(southholland2021,tag_manufacturer_name == "HQXS" ) #remove microwave tags
ch_dum2021HQXS <- subset(ch_dum2021,tag_manufacturer_name == "HQXS" )  #remove microwave tags

#HQXS location data
teamPiersmaHQXSdata <- as.data.frame(getMovebankData(study=1563249841,login=login, removeDuplicatedTimestamps=T))
teamPiersmaHQXS2022data <- as.data.frame(getMovebankData(study=2083443328,login=login, removeDuplicatedTimestamps=T))
BtgTagus2021data <- as.data.frame(getMovebankData(study=1693518103, login=login, removeDuplicatedTimestamps=T))
BtgTagus2021data <- filter(BtgTagus2021data, local_identifier != "Mouchao") #remove icelandic
HQXS_Black_tailed_godwitsData <- as.data.frame(getMovebankData(study=1658294759, login=login, removeDuplicatedTimestamps=T))
southholland2021HQXSdata <- subset(southholland2021data,sensor_type == "GPS")
ch_dum2021HQXSdata <- subset(ch_dum2021data,sensor_type == "GPS")

#combine HQXS reference data
combinedHQXS <- rbind.fill(teamPiersmaHQXS,teamPiersmaHQXS2022,BtgTagus2021, HQXS_Black_tailed_godwits,southholland2021HQXS,ch_dum2021HQXS)
combinedHQXS <- filter(combinedHQXS, animal_id != "NA")

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
combinedHQXSdata <- rbind.fill(teamPiersmaHQXSdata,teamPiersmaHQXS2022,data,BtgTagus2021data,HQXS_Black_tailed_godwitsData,
                               southholland2021HQXSdata,ch_dum2021HQXSdata)

############################ICARUS##########################################################################################################################################################
#Icarus reference data
icarus <- as.data.frame(getMovebankReferenceTable(study=1487044886,login=login))
icarus <- icarus[!duplicated(icarus$animal_local_identifier), ] #remove duplicate rows

#add study name column
icarus$study_name <- icarus$study_id
icarus["study_name"][icarus["study_name"] == 1487044886] <- "icarus black-tailed godwits theunis piersma"

#Icarus location data
icarusData <-  as.data.frame(getMovebankData(study=1487044886,login=login, removeDuplicatedTimestamps=T))

############################Lotek##########################################################################################################################################################
#Lotek reference data
ad_dum2020lotek <- subset(ad_dum2020,tag_manufacturer_name =="Lotek") #remove 2 microwave tags
ad_dum2021 <- as.data.frame(getMovebankReferenceTable(study=1482506572,login=login2)) #all tags are LOTEK in ad_dum2021
ad_dum2021 <- ad_dum2021[!duplicated(ad_dum2021$animal_local_identifier), ]#remove duplicate rows
extremadura2022 <- as.data.frame(getMovebankReferenceTable(study=1923591036,login=login))
extremadura2022 <- filter(extremadura2022,)
extremadura2022 <- extremadura2022[!duplicated(extremadura2022$animal_local_identifier), ] #remove duplicate rows

#Lotek location data
ad_dum2020lotekData <- subset(ad_dum2020data,sensor_type == "GPS")
ad_dum2021data <- as.data.frame(getMovebankData(study=1482506572,login=login2, removeDuplicatedTimestamps=T))
extremadura2022data <- as.data.frame(getMovebankData(study=1923591036,login=login, removeDuplicatedTimestamps=T))

#combine Lotek reference data
combinedLotek <- rbind.fill(ad_dum2020lotek,ad_dum2021,extremadura2022)

combinedLotek <- combinedLotek[!duplicated(combinedLotek$animal_local_identifier), ]#remove duplicate rows

#add study name column
combinedLotek$study_name <- combinedLotek$study_id
combinedLotek["study_name"][combinedLotek["study_name"] == 1105026166] <- "dummer adults 2020"
combinedLotek["study_name"][combinedLotek["study_name"] == 1482506572] <- "dummer adults 2021"
combinedLotek["study_name"][combinedLotek["study_name"] == 1923591036] <- "Black-tailed Godwits Extremadura (Lotek)"

#combine Lotek location data 
combinedLotekData <- rbind.fill(ad_dum2020lotekData,ad_dum2021data, extremadura2022)

############################Active tags############################################################################################################################

activeMicrowave <- filter(combinedMicrowave, animal_timestamp_end > "2022-03-01")
activeMicrowaveSpring2021 <- filter(combinedMicrowave, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01")
activeHQXS <- filter(combinedHQXS, animal_timestamp_end > "2022-03-01")
activeIcarus <- filter(icarus, animal_timestamp_end > "2022-03-01")
activeLotek <- filter(combinedLotek, animal_timestamp_end > "2022-03-01")
activeLotekSpring2021 <- filter(combinedLotek, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01")
activeTags <- rbind.fill(activeMicrowave,activeHQXS,activeIcarus,activeLotek)

activeMicrowaveL <- nrow(filter(combinedMicrowave, animal_timestamp_end > "2022-03-01"))
activeMicrowaveSpring2021L <- nrow(filter(combinedMicrowave, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01"))
activeHQXSL <- nrow(filter(combinedHQXS, animal_timestamp_end > "2022-03-01"))
activeIcarusL <- nrow(filter(icarus, animal_timestamp_end > "2022-03-01"))
activeLotekL <- nrow(filter(combinedLotek, animal_timestamp_end > "2022-03-01"))
activeLotekSpring2021L <- nrow(filter(combinedLotek, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01"))

############################Combine############################################################################################################################
#combine all reference data
combinedReferenceData <- rbind.fill(combinedMicrowave,combinedHQXS,icarus,combinedLotek)

#combine all location data
combinedLocationData <- rbind.fill(combinedMicrowaveData,combinedHQXSdata,icarusData,combinedLotekData)

#combine reference and location data
library('data.table') # for setnames function
setnames(combinedReferenceData, "animal_local_identifier", "local_identifier") #change name to match with location data table
allLocations <-  merge(combinedReferenceData,combinedLocationData,by="local_identifier",all=FALSE)
allLocations <- subset(allLocations, select = - c(tag_local_identifier.x, sensor_type_id.x, number_of_location_events.x,
                                            animal_number_of_deployments.x, animal_sensor_type_ids.x, tag_id.x, 
                                            tag_comments, deployment_id.x, deployment_comments, deploy_on_person, 
                                            deployment_local_identifier, study_id.x, animal_death_comments, 
                                            deploy_off_person, animal_exact_date_of_birth, 
                                            animal_reproductive_condition, tag_production_date, attachment_type, 
                                            tag_readout_method, tag_id.y, sensor_type_id.y, gps_fix_type_raw, 
                                            import_marked_outlier, lotek_crc_status, lotek_crc_status_text, 
                                            update_ts, algorithm_marked_outlier, argos_altitude, argos_best_level, 
                                            argos_calcul_freq, argos_iq, argos_lat1, argos_lat2, argos_lon1, 
                                            argos_lon2, argos_nb_mes, argos_nb_mes_120, argos_nopc, 
                                            argos_pass_duration,argos_sensor_1, argos_sensor_2, argos_sensor_3, 
                                            argos_sensor_4, visible, deployment_id.y, event_id, sensor_type, 
                                            tag_local_identifier.y, location_long.1, location_lat.1, optional, 
                                            sensor, timestamps, trackId, comments, death_comments, earliest_date_born, 
                                            exact_date_of_birth, individual_id, latest_date_born, nick_name, 
                                            sex, taxon_canonical_name, timestamp_start, timestamp_end, 
                                            number_of_events, number_of_deployments, sensor_type_ids, taxon_detail, 
                                            argos_gdop, argos_orientation, argos_semi_major, argos_semi_minor, 
                                            argos_sat_id, internal_temperature, tag_battery_voltage, gt_activity_count, 
                                            gt_sys_week, gt_tx_count, manually_marked_outlier, height_above_ellipsoid, 
                                            mortality_status, manually_marked_valid, data_decoding_software,
                                            gps_horizontal_accuracy_estimate, gps_speed_accuracy_estimate, 
                                            gps_time_to_fix, icarus_ecef_vx, icarus_ecef_vy, icarus_ecef_vz, 
                                            icarus_ecef_x, icarus_ecef_y, icarus_ecef_z, icarus_reset_counter, 
                                            icarus_timestamp_accuracy, icarus_timestamp_source, icarus_uplink_counter,
                                            transmission_protocol, acceleration_raw_x, acceleration_raw_y, 
                                            acceleration_raw_z, barometric_pressure, icarus_barometric_pressure_raw, 
                                            icarus_battery_voltage_raw, icarus_internal_reference_raw, 
                                            icarus_relative_humidity, icarus_relative_humidity_raw, icarus_resistance_raw,
                                            icarus_solar_cell_current, icarus_solar_cell_current_raw, icarus_solar_cell_voltage_raw, 
                                            icarus_temperature0_raw, icarus_temperature1_raw, icarus_temperature2_raw, 
                                            magnetic_field_raw_x, magnetic_field_raw_y, magnetic_field_raw_z, magnetic_field_x,
                                            magnetic_field_y, magnetic_field_z))
allLocations <- subset(allLocations,select = c(local_identifier,study_name.x,timestamp,location_lat,location_long))
                      
############################Export######################################################################################################################
# #All data
#write.csv(allLocations,"allLocations.csv", row.names = FALSE)
#
# #All data since xxxx-xx-xx
allLocationsSince <- allLocations%>%filter(timestamp>"2022-04-15" & timestamp < "2022-05-10")%>%
  write.csv("allLocationsSpring2022.csv", row.names = FALSE)
allLocationsPreBreeding <- allLocations%>%filter(timestamp>"2022-01-01" & timestamp < "2022-03-01")%>%
  write.csv("allLocationsSpring2022.csv", row.names = FALSE)

#last timestamps of all tags
allLastTimestamps <- allLocations %>% 
   group_by(local_identifier) %>%
   slice(which.max(timestamp))

# #Last timestamps since xxxx-xx-xx
allLocationsSince <- allLastTimestamps%>%filter(timestamp>"2022-04-01" & timestamp < "2022-05-10")%>%
  write.csv("allLocationsSpring2022.csv", row.names = FALSE)
allLocationsPreBreeding <- allLastTimestamps%>%filter(timestamp>"2022-01-01" & timestamp < "2022-03-01")%>%
  write.csv("allLocationsSpring2022.csv", row.names = FALSE)

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
# ################################post-fledge period############################################################
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
# ##############################################kmz############################################################################################
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

#################################Transmitter Performance##############################################

# #combine spring 2021 data
# locationsSince2021 <- filter(locations, timestamp_start > "2021-01-01")
# 
# #get number of total tags used
# totalTagsMicrowave <- nrow(combinedMicrowave)
# totalTagsHQXS <- nrow(combinedHQXS)
# totalTagsIcarus <- nrow(icarus)
# totalTagsLotek <- nrow(combinedLotek)
# 
# #tidy timestamps
combinedMicrowave$animal_timestamp_start <- as.Date(combinedMicrowave$animal_timestamp_start)
combinedMicrowave$animal_timestamp_end <- as.Date(combinedMicrowave$animal_timestamp_end)
combinedHQXS$animal_timestamp_start <- as.Date(combinedHQXS$animal_timestamp_start)
combinedHQXS$animal_timestamp_end <- as.Date(combinedHQXS$animal_timestamp_end)
icarus$animal_timestamp_start <- as.Date(icarus$animal_timestamp_start)
icarus$animal_timestamp_end <- as.Date(icarus$animal_timestamp_end)
combinedLotek$animal_timestamp_start <- as.Date(combinedLotek$animal_timestamp_start)
combinedLotek$animal_timestamp_end <- as.Date(combinedLotek$animal_timestamp_end)
# 
# #get number of active tags
# activeMicrowave <- filter(combinedMicrowave, animal_timestamp_end > "2022-03-01")
# activeMicrowaveSpring2021 <- filter(combinedMicrowave, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01")
# activeHQXS <- filter(combinedHQXS, animal_timestamp_end > "2022-03-01")
# activeIcarus <- filter(icarus, animal_timestamp_end > "2022-03-01")
# activeLotek <- filter(combinedLotek, animal_timestamp_end > "2022-03-01")
# activeLotekSpring2021 <- filter(combinedLotek, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01")
# activeTags <- rbind.fill(activeMicrowave,activeHQXS,activeHQXS,activeLotek)
# 
# activeMicrowaveL <- nrow(filter(combinedMicrowave, animal_timestamp_end > "2022-03-01"))
# activeMicrowaveSpring2021L <- nrow(filter(combinedMicrowave, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01"))
# activeHQXSL <- nrow(filter(combinedHQXS, animal_timestamp_end > "2022-03-01"))
# activeIcarusL <- nrow(filter(icarus, animal_timestamp_end > "2022-03-01"))
# activeLotekL <- nrow(filter(combinedLotek, animal_timestamp_end > "2022-03-01"))
# activeLotekSpring2021L <- nrow(filter(combinedLotek, animal_timestamp_start > "2021-01-01" & animal_timestamp_end > "2022-03-01"))
# 
# #add manufacturer column
combinedMicrowave['Manufacturer'] <- 'Microwave'
combinedHQXS['Manufacturer'] <- 'HQXS'
icarus['Manufacturer'] <- 'Icarus'
combinedLotek['Manufacturer'] <- 'Lotek'
# 
# #mean number of days in use
# daysInUseMicrowave <- as.numeric(mean(na.rm = TRUE, Sys.Date()-combinedMicrowave$timestamp_start))
# daysInUseHQXS <- as.numeric(mean(na.rm = TRUE,Sys.Date()-combinedHQXS$timestamp_start))
# daysInUseIcarus <- as.numeric(mean(na.rm = TRUE,Sys.Date()-icarus$timestamp_start))
# daysInUseLotek <- as.numeric(mean(na.rm = TRUE,Sys.Date()-combinedLotek$timestamp_start))
# daysInUseMicrowaveSpring2021 <- as.numeric(mean(na.rm = TRUE, Sys.Date()-combinedMicrowaveSpring2021$timestamp_start))
# 
# #locations per day
locationsPerDayMicrowave<- mean((combinedMicrowave$number_of_location_events)/as.numeric(Sys.Date()-combinedMicrowave$animal_timestamp_start),na.rm=TRUE)
locationsPerDayHQXS <- mean((combinedHQXS$number_of_location_events)/as.numeric(Sys.Date()-combinedHQXS$animal_timestamp_start),na.rm=TRUE)
locationsPerDayIcarus <- mean((icarus$number_of_location_events)/as.numeric(Sys.Date()-icarus$animal_timestamp_start),na.rm=TRUE)
locationsPerDayLotek <- mean((combinedLotek$number_of_location_events)/as.numeric(Sys.Date()-combinedLotek$animal_timestamp_start),na.rm=TRUE)
# locationsPerDayMicrowaveSpring2021<- mean((combinedMicrowaveSpring2021$number_of_events)/as.numeric(Sys.Date()-combinedMicrowaveSpring2021$timestamp_start),na.rm=TRUE)
# 
# #total locations
# totalLocationsMicrowave <- as.numeric(sum(combinedMicrowave[18]))
# totalLocationsHQXS <- as.numeric(sum(combinedHQXS[18]))
# totalLocationsIcarus <- as.numeric(sum(icarus[18]))
# totalLocationsLotek <- as.numeric(sum(combinedLotek[18]))
# totalLocationsMicrowaveSpring2021 <- as.numeric(sum(combinedMicrowaveSpring2021[18]))
# 
# #remove inactive tags
# #combinedDataActiveTags <- filter(combinedData, timestamp_end > "2021-09-01")
# #activeMicrowaveTable <- subset(combinedDataActiveTags,Manufacturer=="Microwave")
# #activeMicrowave <- nrow(activeMicrowaveTable)
# #activeMicrowaveSpring2021 <- filter(activeMicrowaveTable,timestamp_start > "2021-01-01")
# #activeHQXSTable <- subset(combinedDataActiveTags,Manufacturer=="HQXS")
# #activeHQXS <- nrow(activeHQXSTable)
# #activeIcarusTable <- subset(combinedDataActiveTags,Manufacturer=="Icarus")
# #activeIcarus <- nrow(activeIcarusTable)
# #activeLotekTable <- subset(combinedDataActiveTags,Manufacturer=="Lotek")
# #activeLotek <-nrow(activeLotekTable)
# #totalActive <- nrow(rbind(activeMicrowaveTable,activeHQXSTable,activeIcarusTable))
# 
# #mean transmitter duration of use
# timediffMicrowave <- na.omit(as.data.frame(combinedMicrowave$timestamp_end-combinedMicrowave$timestamp_start))
# meanTransmitterLifespanMicrowave <- as.numeric(mean(timediffMicrowave[,1]))
# timediffHQXS <- na.omit(as.data.frame(combinedHQXS$timestamp_end-combinedHQXS$timestamp_start))
# meanTransmitterLifespanHQXS <- as.numeric(mean(timediffHQXS[,1]))
# timediffIcarus <- na.omit(as.data.frame(icarus$timestamp_end-icarus$timestamp_start))
# meanTransmitterLifespanIcarus <- as.numeric(mean(timediffIcarus[,1]))
# timediffLotek <- na.omit(as.data.frame(combinedLotek$timestamp_end-combinedLotek$timestamp_start))
# meanTransmitterLifespanLotek <- as.numeric(mean(timediffLotek[,1]))
# 
# #combine studies for plotting
locationsPerDay <- c(locationsPerDayMicrowave,locationsPerDayHQXS,locationsPerDayIcarus,locationsPerDayLotek)
# locationsPerDaySpring2021 <- c(locationsPerDayMicrowaveSpring2021,locationsPerDayHQXS,locationsPerDayIcarus,locationsPerDayLotek)
# totalLocations <- c(totalLocationsMicrowave,totalLocationsHQXS,totalLocationsIcarus,totalLocationsLotek)
# totalLocationsSpring2021 <- c(totalLocationsMicrowaveSpring2021,totalLocationsHQXS,totalLocationsIcarus,totalLocationsLotek)
# daysInUse <- c(daysInUseMicrowave,daysInUseHQXS,daysInUseIcarus,daysInUseLotek)
# daysInUseSpring2021 <- c(daysInUseMicrowaveSpring2021,daysInUseHQXS,daysInUseIcarus,daysInUseLotek)
# active <- c(activeMicrowave,activeHQXS,activeIcarus,activeLotek)
# activeSpring2021 <- c(activeMicrowaveSpring2021,activeHQXS,activeIcarus,activeLotek)
# meanTransmitterLifespan <- c(meanTransmitterLifespanMicrowave,meanTransmitterLifespanHQXS,meanTransmitterLifespanIcarus,meanTransmitterLifespanLotek)
# totalTags <- c(nrow(combinedMicrowave),nrow(combinedHQXS),nrow(icarus),nrow(combinedLotek))
# totalTagsSpring2021 <- c(nrow(combinedMicrowaveSpring2021),nrow(combinedHQXS),nrow(icarus),nrow(combinedLotek))
manufacturer <- c("Microwave","HQXS","Icarus","Lotek")
# combined <- data.frame(manufacturer,locationsPerDay,totalLocations,daysInUse,active,meanTransmitterLifespan)
# combinedSpring2021 <- data.frame(manufacturer,locationsPerDaySpring2021,totalLocationsSpring2021,daysInUseSpring2021,activeSpring2021,meanTransmitterLifespan)
# 
# #location quality Microwave/Lotek
# table(combinedMicrowaveData$argos_lc)
# table(combinedLotekData$argos_lc)
# lotekGPS <- subset(combinedLotekData, sensor_type=="GPS")
# lotekArgos <- subset(combinedLotekData, sensor_type=="Argos Doppler Shift")
# nrow(lotekArgos)
# nrow(lotekGPS)
# nrow(combinedLotekData)
# 
# ############################Plots##############################################
# 
# 
# #combined locations/day plot
 pl1 <- ggplot(data=combinedLocationsPerDay,aes(x=manufacturer,y=locationsPerDay))
 plot1 <- pl1+geom_bar(stat='identity', color='steelblue4',fill='steelblue4',width = 0.4)+
   ylab("Locations/Day") +xlab("")+
   geom_text(aes(label=locationsPerDay),vjust=-0.2, color="black",
             size=6.5, label = round(locationsPerDay,2)) + ylim(0,3) + theme(axis.title = element_text(size = 15))+ theme(axis.text = element_text(size = 15))             
 plot1
 
# 
# pl2 <- ggplot(data=combined,aes(x=manufacturer,y=active))
# plot2 <- pl2+geom_bar(stat='identity', color='steelblue4',fill='steelblue4',width = 0.4)+
#   ylab("Active tags")+xlab("")+
#   geom_text(aes(label=active),vjust=-0.3, color="black",
#             size=3.5) +ylim(0,80)
# 
# pl3 <- ggplot(data=combined,aes(x=manufacturer,y=totalTags))
# plot3 <- pl3 +geom_bar(stat='identity', color='steelblue4',fill='steelblue4',width = 0.4)+
#   ylim(0,450)+
#   ylab("Total tags deployed")+xlab("")+
#   geom_text(aes(label=totalTags),vjust=-0.3, color="black",
#             size=3.5)
# 
# ggarrange(plot1, plot2, plot3, ncol = 2, nrow = 2)


############################all inactive over 30 days########################################################################################################################

# locations$timestamp <- as.POSIXct(locations$timestamp)
# locations$animal_timestamp_end <- as.POSIXct(locations$animal_timestamp_end)
# 
# #last timestamps of all tags
# allLastTimestamps <- locations %>% 
#   group_by(local_identifier) %>%
#   slice(which.max(timestamp))
# 
# #remove unnecessary columns
# allLastTimestamps <- allLastTimestamps[,c("local_identifier", "animal_ring_id.x",
#                                             "timestamp", "location_lat", "location_long")]
# 
# #inactive for over 30 days
# library(lubridate)
# allNonActive <- allLastTimestamps %>%
#   filter(timestamp  <= today() - days(30))
# write.csv(allNonActive,"allNonActive.csv", row.names = FALSE)
# 
# #voltage <3.65 (active in past 30 days)
# reduceFixes <- filter(HQXSactive, tag_voltage < 3650)
# 


