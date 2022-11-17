setwd('C://Users/adelv/Dropbox/NEON/CSVs_R/')
methDB<-read.csv('./External data/MethDB_Concentrations.csv')
sitesDB<-read.csv('./External data/MethDB_Sites.csv')

names(methDB)

methdb<-merge(methDB, sitesDB, by=c('Site_Nid', 'Site_Name'))
hist(methdb$lat_decimal)
north<-methdb[(methdb$lat_decimal>=23 & methdb$lat_decimal<70) | (methdb$lat_decimal<=-23 & methdb$lat_decimal>-70),]
#as in not in tropics

ggplot(methdb) +
  geom_bar(aes(x=Sample_Season))

seasons<- c('Spring', 'Summer', 'Fall', 'Winter', 'all')
north2<-north[north$Sample_Season %in% seasons == TRUE,]

methdb_site<- methDB %>%
  group_by(Site_Nid) %>%
  summarise(meanmeth = mean(concmean_uM, na.rm=TRUE),
            meantempa = mean(WaterTemp_actual, na.rm=TRUE),
            meantempe = mean(WaterTemp_est, na.rm=TRUE),
            meanDOC = mean(DOC_actual_uM/1000, na.rm=TRUE),
            meanTN = mean(TN_actual_uM/1000, na.rm=TRUE))


mdbtempmin<-min(c(min(methdb_site$meantempa, na.rm=TRUE), min(methdb_site$meantempe, na.rm=TRUE)))
mdbtempmax<-max(c(max(methdb_site$meantempa[methdb_site$meantempa<60], na.rm=TRUE),
                  max(methdb_site$meantempe[methdb_site$meantempe<60], na.rm=TRUE)))


range(methdb_site$meanmeth, na.rm=TRUE)

methdb_site$Site_Nid[which(methdb_site$meanmeth == max(methdb_site$meanmeth, na.rm=TRUE))]
methDB$Site_Name[methDB$Site_Nid==36180]
#Uatuma River seems to be approximately sea level
#lowest is 0
HCH4<-0.000014 #Sander 2015
HCH4_MDB<-HCH4*exp(2400*((1/(methdb_site$meantempa[279]+273.15))-(1/298.15)))
#no temperature either actual or estimated
HCH4_MDB<-HCH4*exp(2400*((1/(26+273.15))-(1/298.15))) #right now wikipedia mean 26 to 17 degree mean annually
ppm_max<- 77/((1*10^6)*(HCH4_MDB/1000))
log10(ppm_max)
#approximate max ppm for CH4 is 3.75 then? range 0 to 3.75



methdbs<-merge(methdb_site, sitesDB, by='Site_Nid')
#966

methdbs<-methdbs[is.na(methdbs$meanmeth)==FALSE,]
#933 obs

methdbs<-methdbs[is.na(methdbs$Landuse)==FALSE,]

#want to convert all other land uses to 'other'
included<-c('Wetland', 'Agriculture', 'Urban')

methdbs<-methdbs[is.na(methdbs$Landuse)==FALSE,]
methdbs<-methdbs[methdbs$Landuse!='',]
methdbs$newlanduse<-as.character(methdbs$Landuse)
methdbs$newlanduse[!(methdbs$Landuse %in% included)]<-'Other'

ws$NEON<-rep('NEON', length(ws$site))
#ws comes from the NLCD table load_changed categories script where the
#landuse was manually reclassified the same way

