setwd('C://Users/adelv/Dropbox/NEON/CSVs_R/')
methDB<-read.csv('./External data/MethDB_Concentrations.csv')
sitesDB<-read.csv('./External data/MethDB_Sites.csv')

names(methDB)

methdb<-merge(methDB, sitesDB, by='Site_Nid')
hist(methdb$lat_decimal)
north<-methdb[(methdb$lat_decimal>=23 & methdb$lat_decimal<70) | (methdb$lat_decimal<=-23 & methdb$lat_decimal>-70),]
#as in not in tropics

ggplot(methdb) +
  geom_bar(aes(x=Sample_Season))

seasons<- c('Spring', 'Summer', 'Fall', 'Winter', 'all')
north2<-north[north$Sample_Season %in% seasons == TRUE,]

methdb<-merge(methDB, sitesDB, by='Site_Nid')


methdb_site<- methDB %>%
  group_by(Site_Nid) %>%
  summarize(meanmeth = mean(concmean_uM, na.rm=TRUE),
            meantempa = mean(WaterTemp_actual, na.rm=TRUE),
            meantempe = mean(WaterTemp_est, na.rm=TRUE),
            meanDOC = mean(DOC_actual_uM/1000, na.rm=TRUE),
            meanTN = mean(TN_actual_uM/1000, na.rm=TRUE))

mdbtempmin<-min(c(min(methdb_site$meantempa, na.rm=TRUE), min(methdb_site$meantempe, na.rm=TRUE)))
mdbtempmax<-max(c(max(methdb_site$meantempa[methdb_site$meantempa<60], na.rm=TRUE),
                  max(methdb_site$meantempe[methdb_site$meantempe<60], na.rm=TRUE)))

log10(range(methdb_site$meanDOC, na.rm=TRUE))
log10(range(methdb_site$meanTN, na.rm=TRUE))


range(methdb_site$meanmeth, na.rm=TRUE)
methdb_site$Site_Name[which(methdb_site$meanmeth == max(methdb_site$meanmeth, na.rm=TRUE))]
methDB$Site_Name[methDB$Site_Nid==36180]
#Uatuma River seems to be approximately sea level
#lowest is 0
HCH4<-0.000014 #Sander 2015
HCH4_MDB<-HCH4*exp(2400*((1/(methdb_site$meantemp[279]+273.15))-(1/298.15)))
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


#get the values from methdbs
names(methdbs)
mdbmin<-log10(min(methdbs$meanmeth[methdbs$meanmeth>0]))
mdbmax<-log10(max(methdbs$meanmeth))
mdbtempmin<-min(c(min(methdbs$meantempa, na.rm=TRUE), min(methdbs$meantempe, na.rm=TRUE)))
mdbtempmax<-max(c(max(methdbs$meantempa[methdbs$meantempa<60], na.rm=TRUE),
                  max(methdbs$meantempe[methdbs$meantempe<60], na.rm=TRUE)))

tiff('ghgprecip.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
  geom_rect(aes(ymin=1.8, ymax=3.7, xmin=0, xmax=3.2), fill=NA, alpha=0.3, colour='grey20', 
            linetype='solid', size=0.3)+
  geom_point(aes(x=meanprecip_mm/1000, y=log10(meanpCO2)), size=2, colour = 'grey30') +
  geom_errorbar(aes(x=meanprecip_mm/1000, ymin=log10(meanpCO2-sepCO2), ymax=log10(meanpCO2+sepCO2))) +
  geom_errorbar(aes(x=meanprecip_mm/1000, ymin=log10(meanpCH4-sepCH4), 
                    ymax=log10(meanpCH4+sepCH4)), alpha=0.8, colour='coral2') +
  geom_point(aes(x=meanprecip_mm/1000, y=log10(meanpCH4)), size=2, colour = 'coral2') +
  geom_errorbar(aes(x=meanprecip_mm/1000, ymin=(meanpN2O-sepN2O), 
                    ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
  geom_point(aes(x=meanprecip_mm/1000, y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
  theme_classic() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)) +
  scale_x_continuous(limits=c(0, 3.5), breaks = seq(0, 3.5, 0.5)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0), 
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='P [m a^-1]')
dev.off()
#no significant relationships so no lines

hist(sites$meanpCO2)
co2mod<-lm(log10(meanpCO2)~I(meanprecip_mm/1000), data=sites)
summary(co2mod)

hist(sites$meanpCH4, breaks=20)
CHmod<-lm(log10(meanpCH4)~I(meanprecip_mm/1000), data=sites)
summary(CHmod)

hist(sites$meanpN2O, breaks=20)
N2Omod<-lm(meanpN2O~I(meanprecip_mm/1000), data=sites)
summary(N2Omod)

N2Omin<-0.2
N2Omax<-37.4

#have edited this so that the methane units are concentrations and water temperature
tiff('gases_temp.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
  geom_rect(aes(ymin=1.8, ymax=3.7, xmin=-3, xmax=23), fill=NA, alpha=0.3, colour='grey20', 
            linetype='solid', size=.3)+
  geom_vline(aes(xintercept=(mdbtempmin)), colour='coral2', size=0.9, alpha=0.6, linetype='dashed')+
  geom_vline(aes(xintercept=(mdbtempmax)), colour='coral2', size=0.9, alpha=0.6, linetype='dashed')+
  geom_rect(aes(ymin=0, ymax=3.75, #from the Historgrams of GHG compilations file
                xmin=mdbtempmin, xmax=mdbtempmax), fill=NA, alpha=0.1, colour='coral2', 
            linetype='solid', size=0.3)+
  geom_point(aes(x=meanairtemp, y=log10(meanpCO2)), size=2, colour = 'grey30') +
  geom_abline(aes(slope=0.019954, intercept=3.051690), color='grey30', linetype='solid') +
  geom_errorbar(aes(x=meanairtemp, ymin=log10(meanpCO2-sepCO2), ymax=log10(meanpCO2+sepCO2))) +
  geom_errorbar(aes(x=meanairtemp, ymin=log10(meanpCH4-sepCH4), 
                    ymax=log10(meanpCH4+sepCH4)), alpha=0.8, colour='coral2') +
  geom_point(aes(x=meanairtemp, y=log10(meanpCH4)), size=2, colour = 'coral2') +
  geom_abline(aes(slope=0.05313, intercept=1.16371), color='coral2', linetype='solid') +
 # geom_rect(aes(ymin=N2Omin, ymax=N2Omax, xmin=9.9, xmax=27.1), fill=NA, alpha=0.3, colour='grey30', 
  #          linetype='dotted', size=1.2)+ #using the temp values from Beaulieau but the concentratiosn from Yao
 #the full rect wont fit with the actual max values from Yao, so will need to just plot the three lines
 # geom_vline(aes(xintercept=(9.9)), colour='grey60', size=1.2, alpha=0.7, linetype='dotted')+
  #geom_vline(aes(xintercept=(27.1)), colour='grey60', size=1.2, alpha=0.7, linetype='dotted')+
   geom_errorbar(aes(x=meanairtemp, ymin=(meanpN2O-sepN2O), 
                    ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
  geom_point(aes(x=meanairtemp, y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
  geom_abline(aes(slope=0.017224, intercept=0.546253), color='grey30', linetype='longdash') +
  theme_classic() +
    scale_y_continuous(limits = c(0, 4), breaks = seq(0,4, 1)) +
  scale_x_continuous(limits=c(-4, 36), breaks = seq(-4, 36, 4)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0), 
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='Temperature [deg C]')
dev.off()
#methDB is water temperature and approximate ppm.  will need to do a separate plot?

hist(sites$meanpCO2)
co2mod<-lm(log10(meanpCO2)~meanairtemp, data=sites)
summary(co2mod)
#significant

hist(sites$meanpCH4, breaks=20)
CHmod<-lm(log10(meanpCH4)~meanairtemp, data=sites)
summary(CHmod)
#also significant

hist(sites$meanpN2O, breaks=20)
N2Omod<-lm(meanpN2O~meanairtemp, data=sites)
summary(N2Omod)
#and also significant

#and now slope
tiff('gases_slope.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
  geom_rect(aes(ymin=1.8, ymax=3.7, xmin=exp(-1.1), xmax=exp(3.4)), fill=NA, alpha=0.3, colour='grey30', 
            linetype='solid', size=0.3)+
  geom_point(aes(x=slope, y=log10(meanpCO2)), size=2, colour = 'grey30') +
  geom_abline(aes(slope=-0.016132, intercept=3.446441), color='grey30', linetype='solid') +
  geom_errorbar(aes(x=slope, ymin=log10(meanpCO2-sepCO2), ymax=log10(meanpCO2+sepCO2))) +
  geom_errorbar(aes(x=slope, ymin=log10(meanpCH4-sepCH4), 
                    ymax=log10(meanpCH4+sepCH4)), alpha=0.8, colour='coral2') +
  geom_point(aes(x=slope, y=log10(meanpCH4)), size=2, colour = 'coral2') +
  geom_abline(aes(slope=-0.06633, intercept=2.49430), color='coral2', linetype='solid') +
  geom_errorbar(aes(x=slope, ymin=(meanpN2O-sepN2O), 
                    ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
  geom_point(aes(x=slope, y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
 # geom_abline(aes(slope=0.017224, intercept=0.546253), color='grey30', linetype='longdash') +
  theme_classic() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4, 1)) +
  scale_x_continuous(limits=c(0, 30), breaks = seq(0, 30, 5)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0), 
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='Slope [degrees]')
dev.off()

hist(sites$meanpCO2)
co2mod<-lm(log10(meanpCO2)~slope, data=sites)
summary(co2mod)
#significant

hist(sites$meanpCH4, breaks=20)
CHmod<-lm(log10(meanpCH4)~slope, data=sites)
summary(CHmod)
#also significant

hist(sites$meanpN2O, breaks=20)
N2Omod<-lm(meanpN2O~slope, data=sites)
summary(N2Omod)
#and also significant

#need to summarize the TN data and add that too
#run the TN script

tiff('tn_ghg_sites.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
   geom_rect(aes(ymin=0, ymax=3.75, #already log10 transformed
                xmin=-1.659, xmax=2.364), fill=NA, alpha=0.1, colour='coral2',
             linetype='solid', size=0.3)+
  # geom_vline(aes(xintercept=log10(mdbTNmin)), colour='coral2', size=0.9, alpha=0.6, linetype='dashed')+
  #geom_vline(aes(xintercept=log10(mdbTNmax)), colour='coral2', size=0.9, alpha=0.6, linetype='dashed')+
  #  geom_rect(aes(ymin=1.8, ymax=3.7, xmin=-1.33, xmax=0.31), fill=NA, alpha=0.3, colour='grey70',
  #           linetype='dotted', size=1.2)+
  geom_point(aes(x=log10(meanTN), y=log10(meanpCO2)), size=1, colour = 'grey30') +
  geom_abline(aes(slope=0.39750, intercept=3.47001), color='grey30', linetype='solid') +
  geom_errorbar(aes(x=log10(meanTN), ymin=log10(meanpCO2-sepCO2), ymax=log10(meanpCO2+sepCO2))) +
  geom_errorbarh(aes(xmin=log10(meanTN-seTN), xmax=log10(meanTN+seTN), y=log10(meanpCO2))) +
  geom_errorbar(aes(x=log10(meanTN), ymin=log10(meanpCH4-sepCH4),
                    ymax=log10(meanpCH4+sepCH4)), alpha=0.8, colour='coral2') +
  geom_errorbarh(aes(xmin=log10(meanTN-seTN), xmax=log10(meanTN+seTN), y=log10(meanpCH4)),
                 alpha=0.8, colour='coral2') +
  geom_point(aes(x=log10(meanTN), y=log10(meanpCH4)), size=1, colour = 'coral2') +
  geom_abline(aes(slope=1.2673, intercept=2.3908), color='coral2', linetype='solid', size=1) +
  geom_errorbar(aes(x=log10(meanTN), ymin=(meanpN2O-sepN2O),
                    ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
  geom_errorbarh(aes(xmin=log10(meanTN-seTN), xmax=log10(meanTN+seTN), y=log10(meanpCO2)),
                 alpha=0.8, colour='grey30') +
  geom_point(aes(x=log10(meanTN), y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
  geom_abline(aes(slope=0.34153, intercept=0.91531), color='grey30', linetype='longdash') +
  theme_classic() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4, 1)) +
  scale_x_continuous(limits=c(-3.5, 2.5), breaks = seq(-3, 2, 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='log10(TN) [mg/L]')
dev.off()

tiff('doc_ghg_sites.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
   geom_rect(aes(ymin=0, ymax=3.75, #already log10 transformed
                xmin=-1.617, xmax=0.751), fill=NA, alpha=0.1, colour='coral2',
             linetype='solid', size=0.3) +
  # geom_vline(aes(xintercept=log10(mdbDOCmin)), colour='coral2', size=0.9, alpha=0.6, linetype='dashed')+
  #geom_vline(aes(xintercept=log10(mdbDOCmax)), colour='coral2', size=0.9, alpha=0.6, linetype='dashed')+
  geom_rect(aes(ymin=1.8, ymax=3.7, xmin=-1.33, xmax=0.31), fill=NA, alpha=0.3, colour='grey30',
            linetype='solid', size=0.3)+
  geom_point(aes(x=log10(meanDOC), y=log10(meanpCO2)), size=1, colour = 'grey30') +
  # geom_abline(aes(slope=0.3475, intercept=3.1206), color='grey30', linetype='solid') +
  geom_errorbar(aes(x=log10(meanDOC), ymin=log10(meanpCO2-sepCO2), ymax=log10(meanpCO2+sepCO2))) +
  geom_errorbarh(aes(xmin=log10(meanDOC-seDOC), xmax=log10(meanDOC+seDOC), y=log10(meanpCO2))) +
  geom_errorbar(aes(x=log10(meanDOC), ymin=log10(meanpCH4-sepCH4),
                    ymax=log10(meanpCH4+sepCH4)), alpha=0.8, colour='coral2') +
  geom_errorbarh(aes(xmin=log10(meanDOC-seDOC), xmax=log10(meanDOC+seDOC), y=log10(meanpCH4)),
                 alpha=0.8, colour='coral2') +
  geom_point(aes(x=log10(meanDOC), y=log10(meanpCH4)), size=1, colour = 'coral2') +
  geom_abline(aes(slope=1.9114, intercept=0.9705), color='coral2', linetype='solid', size=1) +
  geom_errorbar(aes(x=log10(meanDOC), ymin=(meanpN2O-sepN2O),
                    ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
  geom_errorbarh(aes(xmin=log10(meanDOC-seDOC), xmax=log10(meanDOC+seDOC), y=meanpN2O),
                 alpha=0.8, colour='grey30') +
  #  geom_vline(aes(xintercept=log10(0.4)), colour='grey60', size=1.2, alpha=0.7, linetype='dotted')+
  # geom_vline(aes(xintercept=log10(25.56)), colour='grey60', size=1.2, alpha=0.7,
  #   linetype='dotted')+
  geom_point(aes(x=log10(meanDOC), y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
  #  geom_abline(aes(slope=0.017224, intercept=0.546253), color='grey30', linetype='longdash') +
  theme_classic() +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4, 1)) +
  scale_x_continuous(limits=c(-1.7, 1.6), breaks = seq(-2, 2, 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='log10(DOC) [mg/L]')
dev.off()
