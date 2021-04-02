#This script will download current NEON files, calculate gas concentrations and measurement error,
#and flag or remove NEON quality flagged data.

library(neonUtilities)
library(ggplot2)
library(dplyr)
library(padr)
library(tidyr)
library(lubridate)
library(gridExtra)
library(glue)
library(stringr)
library(httr)
library(jsonlite)
library(readr)
library(feather)
library(tidyverse)


#######USING NEONUTILITIES - USUALLY GLITCHES WITH ALL SITES###########
#load the files
gases<-loadByProduct(dpID = "DP1.20097.001", package = "basic", #expanded brings in all the readmes
                     startdate = "2015-01", enddate = "2019-12", site='all', #starting with these experimental 3
                     check.size = F)

list2env(gases, .GlobalEnv)



#Kaelin's scripts state that the default water volume is 40 mL, default gas is 20
#replace with default values where NA
##### Default values #####
volH2O <- 40 #mL
volGas <- 20 #mL
sdg_fieldDataProc$waterVolumeSyringe[is.na(sdg_fieldDataProc$waterVolumeSyringe)] <- volH2O
sdg_fieldDataProc$gasVolumeSyringe[is.na(sdg_fieldDataProc$gasVolumeSyringe)] <- volGas


#at this point export the lab data and check that the QF actually represents all QF issues
#then re-import it and run the next part of the script
#filtering will include flagging samples with any lack of clarity on results - will not worry as much
#about the dates yet
#use the external remarks field and the sample condition fields
#on Mar 15, 2021 - <20 samples had issues but were not flagged

write.csv(sdg_externalLabData, 'Lab data_unchecked.csv')
sdg_externalLabData<-read.csv('Lab data_checked.csv', header=TRUE)

write.csv(sdg_fieldDataAir, 'air_unchecked.csv')
sdg_fieldDataAir<-read.csv('air_checked.csv', header=TRUE)
#approximately 10-20 samples changed to flagged based on remarks

#rename the fields to be used for merging for ease (air samples and reference samples all 
#listed on the same lab output - so these will pull appropriate records)

sdg_fieldDataProc<-sdg_fieldDataProc %>%
  rename(sampleID = equilibratedAirSampleID)

sdg_fieldDataAir<-sdg_fieldDataAir %>%
  rename(sampleID = referenceAirSampleID)


concs1<-merge(sdg_externalLabData, sdg_fieldDataProc, by='sampleID')
names(concs1)
#the collect time for the field and lab data differ - going with the field data entry
concs1<-concs1 %>%
  rename(fieldcollectdate = collectDate.y)

concs2<-merge(sdg_externalLabData, sdg_fieldDataAir, by='sampleID')

concs3<-merge(concs1, concs2, by.x = 'referenceAirSampleID', by.y = 'sampleID', all=FALSE)
#this has now merged datasets horizontally so the reference air samples can be linked to
#the headspace samples.  x is the water sample info 

concs<-concs3
names(concs)[2:7]<-c('sampleID', 'uid', 'domainID', 'siteID', 'namedLocation', 'labcollectDate')
names(concs)

concs$date<-as.Date(concs$fieldcollectdate, format = '%Y-%m-%d %H:%M:%S', tz="")
#timezone already standardized to GMT, same as UTC, default
concs$month<-month(concs$date)
#air sample concentrations are .y, water headspace .x

#replace the water temps with the stream temps where missing 
#If run within three hours then stream temp, if not then storage - 
#refer to Kaelin email

#create index for what temp
concs$tempsource[is.na(concs$storageWaterTemp)] <-
  'stream'
concs$tempsource[is.na(concs$storageWaterTemp)==FALSE] <-
  'storage'

#Kaelin's language for replacing the water temp
for(m in 1:length(concs$waterSampleID)){
  try(concs$storageWaterTemp[m] <- sdg_fieldSuperParent$waterTemp[sdg_fieldSuperParent$parentSampleID == concs$waterSampleID[m]],silent = T)
  if(is.na(concs$storageWaterTemp[m])){
    try(
      concs$storageWaterTemp[m] <- fieldSuperParent$waterTemp[fieldSuperParent$parentSampleID == concs$waterSampleID[m]],
      silent = T)
  }
}
#should have replaced the missing SWT values with the water temp of the stream

#also bring in the stream water temperatures
concs$fieldtemp<-NA
for(m in 1:length(concs$waterSampleID)){
  try(concs$fieldtemp[m]<-sdg_fieldSuperParent$waterTemp[sdg_fieldSuperParent$parentSampleID == concs$waterSampleID[m]], silent=T) }

#address quality flagging for low volume injections 
#concs$lowAirVolumeQF #just 1,-1 flag, same for concs$lowgasvolume
#concs$lowGasVolumeQF
#Remove any nonzero (error) values

concsfiltered<-concs[concs$lowAirVolumeQF==0,]
concsfiltered<-concsfiltered[concsfiltered$lowGasVolumeQF==0,]
#these should be removed bc completely unknown error

#and for no standard check
#gasCheckStandardQF:Quality flag indicating when the gas check standard percent deviation 
#is above 2 %; 0 indicates percent deviation less than 2 %; 1 indicates percent deviation 
#2 % or above; and -1 indicates that the test could not be performed

concs$gasCheckStandardQF.x #so many of these have error - 
#need to make this a separate choice to remove
unique(concs$gasCheckStandardQF.x)
concsfiltered_laberr<-concs[is.na(concs$gasCheckStandardQF.x)==TRUE,]
#this lost way too many observations - will choose not to remove but to symbolize differently
#if choice to remove, just set concs as concsfiltered_laberr

concs<-concsfiltered

#next step to flag when BDL - so maybe flags will include the high SD standards or BDL
names(concs)
which(concs3$concentrationCO2.x<concs3$runDetectionLimitCO2.x)

#these are the negative value issues that pop up later?
which(concs3$concentrationCO2.x<concs3$concentrationCO2.y)

#and where are concentrations NA?
which(is.na(concs3$concentrationCO2.x))
which(is.na(concs3$concentrationCO2.y))

hist(concs$runDetectionLimitCO2.x)
#why is this so spread out?
hist(concs$runDetectionLimitCH4.x)
hist(concs$runDetectionLimitN2O.x)

####need to check air concentrations
hist(concs$concentrationCO2.y)
max(concs$concentrationCO2.y, na.rm = TRUE)
qtilemax<-quantile(concs$concentrationCO2.y, 0.90, na.rm=TRUE)
#this also looks very flawed - will remove observations with air concentrations >90% quantile ppm?
which(concs$concentrationCO2.y>qtilemax)
concs$concentrationCO2.y[concs$concentrationCO2.y>qtilemax]<-NA

qtileCH4<-quantile(concs$concentrationCH4.y, 0.90, na.rm=TRUE)
concs$concentrationCH4.y[concs$concentrationCH4.y>qtilemax]<-NA

qtileN2O<-quantile(concs$concentrationN2O.y, 0.90, na.rm=TRUE)
concs$concentrationN2O.y[concs$concentrationN2O.y>qtilemax]<-NA

#remove samples BDL
concs$concentrationCO2.x[concs$concentrationCO2.x<concs$runDetectionLimitCO2.x]<-NA
concs$concentrationCH4.x[concs$concentrationCH4.x<concs$runDetectionLimitCH4.x]<-NA
concs$concentrationN2O.x[concs$concentrationN2O.x<concs$runDetectionLimitN2O.x]<-NA
concs$concentrationCO2.y[concs$concentrationCO2.y<concs$runDetectionLimitCO2.y]<-NA
concs$concentrationCH4.y[concs$concentrationCH4.y<concs$runDetectionLimitCH4.y]<-NA
concs$concentrationN2O.y[concs$concentrationN2O.y<concs$runDetectionLimitN2O.y]<-NA


concs$flagging<-ifelse(is.na(concs$gasCheckStandardQF.x)==TRUE, '>2% SD stds', 
                       ifelse(concs$concentrationCO2.x<concs$concentrationCO2.y, 'air>CO2HS',
                       ifelse(concs$concentrationCH4.x<concs$concentrationCH4.y, 'air>CH4HS',
                      ifelse(concs$concentrationN2O.x<concs$concentrationN2O.y, 'air>N2OHS','valid'))))
#flagging also where the air concentrations used for equilibration were higher than the measured headspace 
#concentrations
                     
           
#NA are resulting anytime the quantitative expressions have a NA
#also replace the flagging NA with 'err' - only run this if using specific error flags
concs$flagging[is.na(concs$flagging)==TRUE]<-'missing data'
#so error also where any data is missing on lab issues


#Next to calculate concentrations and error
#L kPa K-1 mol -1
R<-8.3122598

#H constants in mol m-3 Pa-1
HCO2<-0.00033 #Sander 2015, NEON guide
HN2O<-0.00024 #Sander 2015
HCH4<-0.000014 #Sander 2015


#calculate Henry's law constants and error (mol m-3 Pa-1)
#no longer have error in this because the temperatures dont have error
concs$HCO2_adj<-HCO2*exp(2400*((1/(concs$storageWaterTemp+273.15))-(1/298.15)))
mean(concs$HCO2_adj, na.rm=TRUE)
concs$HCO2_field<-HCO2*exp(2400*((1/(concs$fieldtemp+273.15))-(1/298.15)))

concs$HN2O_adj<-HN2O*exp(2700*((1/(concs$storageWaterTemp+273.15))-(1/298.15)))
concs$HN2O_field<-HN2O*exp(2700*((1/(concs$fieldtemp+273.15))-(1/298.15)))

concs$HCH4_adj<-HCH4*exp(1900*((1/(concs$storageWaterTemp+273.15))-(1/298.15)))
concs$HCH4_field<-HCH4*exp(1900*((1/(concs$fieldtemp+273.15))-(1/298.15)))

#These equations match Kaelin's

concs$CO2aq<-concs$ptBarometricPressure*(10^-6)*
  (((concs$gasVolumeSyringe/1000)*(concs$concentrationCO2.x-concs$concentrationCO2.y))/
                                                   (R*(concs$storageWaterTemp+273.15)*
                                                      (concs$waterVolumeSyringe/1000)) +
                                                   (concs$HCO2_adj*concs$concentrationCO2.x))*
  10^6 #to umol

concs$CO2aq[concs$CO2aq<0]<-NA

######also add in saturation concentrations and differences
hist(concs$concentrationCO2.y)
mean(concs$concentrationCO2.y, na.rm = TRUE) #660?!
concs$CO2_satconc<-(concs$ptBarometricPressure*1000)*409.8*(concs$HCO2_field/1000) #used global avg
hist(concs$CO2_satconc)
concs$CO2_excess<-concs$CO2aq - concs$CO2_satconc #excess umol/L
concs$CO2_pctsat<-concs$CO2aq/concs$CO2_satconc *100#pct saturation
hist(concs$CO2_pctsat, xlim=c(0,200), breaks=10000)

hist(concs$concentrationCH4.y, breaks=500, xlim=c(0,10))

###also do it my way to check
#L/mol at elevation
concs$molperL<-(concs$ptBarometricPressure*0.00986923)/(0.08206*(concs$storageWaterTemp+273)) #in atm
concs$HSconcCO2<-concs$concentrationCO2.x/10^6*concs$molperL #in mol/L
#dont forget the Kh is in mol/m3 pascal
concs$molHS<-concs$HSconcCO2*concs$gasVolumeSyringe/1000 #moles in HS
concs$partpressCO2<-concs$HSconcCO2*R*(concs$storageWaterTemp+273) #used R in kPa here
concs$waterconc<-concs$HCO2_adj*concs$partpressCO2 #mol/L bc pressure in kPa (so x 1000) but convert m3 to m by divie by 1000
concs$molwater<-concs$waterconc*concs$waterVolumeSyringe/1000
concs$molinitair<-concs$concentrationCO2.y/10^6*concs$molperL*concs$gasVolumeSyringe/1000 #mol/L * L
concs$initCO2aq<-(concs$molHS+concs$molwater-concs$molinitair)/(concs$waterVolumeSyringe/1000)*10^6 #umol/L
plot(initCO2aq~CO2aq, data=concs)
#exact - good exercise

#propagate error from measurement precision
hist(concs$precisionCO2.x)#up to ~2.5% (this is already in percent)
concs$CO2aq_err<-sqrt((sqrt((concs$precisionCO2.x/100*concs$concentrationCO2.x)^2+(concs$precisionCO2.y/100*concs$concentrationCO2.y)^2)/
                         (concs$concentrationCO2.x + concs$concentrationCO2.y))^2 + (concs$precisionCO2.x/100)^2)
mean(concs$CO2aq_err, na.rm=TRUE) #proportional error
hist(concs$CO2aq_err)




#10^-6 is the conversion ppm to mol
concs$N2Oaq<-concs$ptBarometricPressure*(10^-6)*(((concs$gasVolumeSyringe/1000)*(concs$concentrationN2O.x-concs$concentrationN2O.y))/
                                                   (R*(concs$storageWaterTemp+273.15)*(concs$waterVolumeSyringe/1000)) +
                                                   (concs$HN2O_adj*concs$concentrationN2O.x))*10^6 #to umol
concs$N2Oaq_err<-sqrt((sqrt((concs$precisionN2O.x/100*concs$concentrationN2O.x)^2+(concs$precisionN2O.y/100*concs$concentrationN2O.y)^2)/
                         (concs$concentrationN2O.x + concs$concentrationN2O.y))^2 + (concs$precisionN2O.x/100)^2)
mean(concs$N2Oaq_err, na.rm=TRUE) #proportional error
concs$n2oaq[concs$N2Oaq<0]<-NA

concs$N2O_satconc<-(concs$ptBarometricPressure*1000)*0.328*(concs$HN2O_field/1000) #used global avg
hist(concs$N2O_satconc)
concs$N2O_excess<-concs$N2Oaq - concs$N2O_satconc #excess umol/L
concs$N2O_pctsat<-concs$N2Oaq/concs$N2O_satconc *100#pct saturation
hist(concs$N2O_pctsat, xlim=c(0,300), breaks=60000)

#10^-6 is the conversion ppm to mol
concs$CH4aq<-concs$ptBarometricPressure*(10^-6)*(((concs$gasVolumeSyringe/1000)*(concs$concentrationCH4.x-concs$concentrationCH4.y))/
                                                   (R*(concs$storageWaterTemp+273.15)*(concs$waterVolumeSyringe/1000)) +
                                                   (concs$HCH4_adj*concs$concentrationCH4.x))*10^6 #to umol
concs$CH4aq_err<-sqrt((sqrt((concs$precisionCH4.x/100*concs$concentrationCH4.x/100)^2+(concs$precisionCH4.y/100*concs$concentrationCH4.y/100)^2)/
                         (concs$concentrationCH4.x/100 + concs$concentrationCH4.y/100))^2 + (concs$precisionCH4.x/100)^2)
mean(concs$CH4aq_err, na.rm=TRUE) #proportional error
concs$CH4aq[concs$CH4aq<0]<-NA

concs$CH4_satconc<-(concs$ptBarometricPressure*1000)*1.851*(concs$HCH4_field/1000) #used global avg
hist(concs$CH4_satconc)
concs$CH4_excess<-concs$CH4aq - concs$CH4_satconc #excess umol/L
concs$CH4_pctsat<-concs$CH4aq/concs$CH4_satconc *100#pct saturation
hist(concs$CH4_pctsat, xlim=c(0,300), breaks=60000)

#format date time
concs$date<-as.Date(concs$fieldcollectdate, format = '%Y-%m-%d %H:%M:%S', tz="")
concs$month<-month(concs$date)

ids<-unique(concs$siteID)

unique(concs$flagging)

write.csv(concs, 'currentGHGdataset.csv')


#next step summarize by date and flag samples where the cv >10 pct?

concs$datetime<-as_datetime(concs$fieldcollectdate, tz='UTC')
concs$genflag<-ifelse(concs$flagging=='valid', 'valid', '>2 error or air>headspace')
concs$genflag_binary<-as.numeric(ifelse(concs$flagging=='valid', '0', '1'))

concs_sum<- concs %>%
  dplyr::select(siteID, date, datetime, storageWaterTemp, fieldtemp, ptBarometricPressure,
                CO2aq, CO2aq_err, CO2_excess, CO2_pctsat, CH4aq, CH4aq_err, 
                CH4_excess, CH4_pctsat, N2Oaq, N2Oaq_err, N2O_excess,
                N2O_pctsat, genflag_binary)

concs_sumdate<- concs_sum %>%
  group_by(siteID, date) %>%
  summarize(CO2aq1 = mean(CO2aq, na.rm=TRUE),
            CO2_excess = mean(CO2_excess, na.rm=TRUE),
            CO2_pctsat = mean(CO2_pctsat, na.rm=TRUE),
            CO2aq_cv = sd(CO2aq, na.rm=TRUE)/mean(CO2aq, na.rm=TRUE),
            CH4aq1 = mean(CH4aq, na.rm=TRUE),
            CH4_excess = mean(CH4_excess, na.rm=TRUE),
            CH4_pctsat = mean(CH4_pctsat, na.rm=TRUE),
            CH4aq_cv = sd(CH4aq, na.rm=TRUE)/mean(CH4aq, na.rm=TRUE),
            N2Oaq1 = mean(N2Oaq, na.rm=TRUE),
            N2O_excess = mean(N2O_excess, na.rm=TRUE),
            N2O_pctsat = mean(N2O_pctsat, na.rm=TRUE),
            N2Oaq_cv = sd(N2Oaq, na.rm=TRUE)/mean(N2Oaq, na.rm=TRUE),
            storageWaterTemp = mean(storageWaterTemp, na.rm=TRUE),
            fieldtemp = mean(fieldtemp, na.rm=TRUE),
            ptBarometricPressure = mean(ptBarometricPressure, na.rm=TRUE),
            flagging = sum(genflag_binary, na.rm=TRUE))
  
#note that the way I've summarized this will not yield times to match, just dates          
hist(concs_sumdate$CO2aq_cv)           
hist(concs_sumdate$CH4aq_cv)
hist(concs_sumdate$N2Oaq_cv)

mean(concs_sumdate$CO2aq_cv, na.rm = TRUE)           
mean(concs_sumdate$CH4aq_cv, na.rm = TRUE)
mean(concs_sumdate$N2Oaq_cv, na.rm = TRUE)

#should remove the ones with CV>50%?


concs_sum1<- concs %>%
  dplyr::select(siteID, date, datetime, storageWaterTemp, fieldtemp, ptBarometricPressure,
                CO2aq, CO2aq_err, CO2_excess, CO2_pctsat, CH4aq, CH4aq_err, 
                CH4_excess, CH4_pctsat, N2Oaq, N2Oaq_err, N2O_excess,
                N2O_pctsat, genflag_binary) %>%
  group_by(siteID, date) %>%
  summarize(CO2aq1 = mean(CO2aq, na.rm=TRUE),
            CO2_excess = mean(CO2_excess, na.rm=TRUE),
            CO2_pctsat = mean(CO2_pctsat, na.rm=TRUE),
            CO2aq_cv = sd(CO2aq, na.rm=TRUE)/mean(CO2aq, na.rm=TRUE),
            CH4aq1 = mean(CH4aq, na.rm=TRUE),
            CH4_excess = mean(CH4_excess, na.rm=TRUE),
            CH4_pctsat = mean(CH4_pctsat, na.rm=TRUE),
            CH4aq_cv = sd(CH4aq, na.rm=TRUE)/mean(CH4aq, na.rm=TRUE),
            N2Oaq1 = mean(N2Oaq, na.rm=TRUE),
            N2O_excess = mean(N2O_excess, na.rm=TRUE),
            N2O_pctsat = mean(N2O_pctsat, na.rm=TRUE),
            N2Oaq_cv = sd(N2Oaq, na.rm=TRUE)/mean(N2Oaq, na.rm=TRUE),
            storageWaterTemp = mean(storageWaterTemp, na.rm=TRUE),
            fieldtemp = mean(fieldtemp, na.rm=TRUE),
            ptBarometricPressure = mean(ptBarometricPressure, na.rm=TRUE),
            flagging = sum(genflag_binary, na.rm=TRUE)) %>%
  filter(CO2aq_cv<=0.5 | is.na(CO2aq_cv) ==TRUE, 
         CH4aq_cv<=0.5 | is.na(CH4aq_cv) ==TRUE,
         N2Oaq_cv<=0.5 | is.na(N2Oaq_cv) == TRUE)

#filtering lost 74 observations
write.csv(concs_sum1, 'Daily mean GHG stats.csv')











