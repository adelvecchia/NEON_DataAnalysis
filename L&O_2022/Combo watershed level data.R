#======================
#Bring in and summarize the concentration data, site level data, imagery data
#=====================

#first the concentration data
setwd('C://Users/adelv/Dropbox/NEON')
concs<-read.csv('concs_sum1.csv')
names(concs)

concs$date<-as.Date(concs$date, format = '%Y-%m-%d')
names(concs)

se<-function(x) {sd(x)/sqrt(length(x))}

#want to summarize boxplots of the concentrations by site as a start

bpframe<-data.frame(site=concs$siteID, date=concs$date, month=month(concs$date),
                    CH4=concs$CH4aq1, N2O=concs$N2Oaq1, CO2=concs$CO2aq1, pCO2 = concs$pCO21,
                    pCH4=concs$pCH41, pN2O = concs$pN2O1, pCO2_err = concs$pCO2_err1,
                    pCH4_err=concs$pCH4_err1, pN2O_err = concs$pN2O_err1, fieldtemp = concs$fieldtemp)

bulkstats<-bpframe %>%
  group_by(site) %>%
  summarise(meanCO2 = mean(CO2, na.rm=TRUE),
            meanpCO2 = mean(pCO2, na.rm=TRUE),
            meanCH4 = mean(CH4, na.rm = TRUE),
            maxCH4 = max(CH4, na.rm = TRUE),
            meanpCH4 = mean(pCH4, na.rm = TRUE),
            meanN2O = mean(N2O, na.rm = TRUE),
            maxN2O = max(N2O, na.rm = TRUE),
            meanpN2O = mean(pN2O, na.rm = TRUE),
            seCO2 = se(CO2),
            seCH4 = se(CH4),
            seN2O = se(N2O),
            sepCO2 = se(pCO2),
            sepCH4 = se(pCH4),
            sepN2O = se(pN2O),
            medCO2 = median(CO2, na.rm=TRUE),
            medCH4 = median(CH4, na.rm = TRUE),
            medN2O = median(N2O, na.rm = TRUE),
            meanfieldtemp = mean(fieldtemp), na.rm=TRUE)

#write.csv(bulkstats, './CSVs_R/Gases/summary_table.csv')

####Run the above, then the imagery script, then load in sites
setwd('C://Users/adelv/Dropbox/NEON/Mapping/neon_ws')
sites_slopes<-read.csv('strslopes.csv', header=TRUE)
names(sites_slopes)
sites_slopes1<-merge(sites_slopes, bulkstats, by='site')

#load the NEON provided site-level data that has already been cropped to columns of interest
setwd('C://Users/adelv/Dropbox/NEON')
siteinfo<-read.csv('./CSVs_R/Site level/sitecsv.csv', header=TRUE)
names(siteinfo)

sites_slopes_metad<-merge(siteinfo, sites_slopes1, by='site')

names(sites_slopes_metad)

siteinfo2<-read.csv('./CSVs_R/Site level/metadata_cropped.csv', header=TRUE)
names(siteinfo2)

sites_slopes_metad2<-merge(siteinfo2, sites_slopes_metad, by='site')

names(sites_slopes_metad2)
#now has the basic metadata and the gas concentrations and slopes

#load the imagery-summarized slopes and stream slopes
#load the NLCD data, also summarized to columns of interest already
setwd('C://Users/adelv/Dropbox/NEON/CSVs_R/Site level')
#nlcd<-read.csv('watershed NLCD_cropped.csv', header=TRUE)
nlcd<-read.csv('NLCD_summed_newcats.csv', header=TRUE)

library(reshape2)

#nlcd<-melt(nlcd, id.vars='siteID')
names(nlcd)
names(sites_slopes_metad2)

sites<-merge(sites_slopes_metad2, nlcd, by='site')
names(sites)
#write.csv(sites, 'site_summary.csv')

#bring in the tn and doc
names(tndoc)[1]<-'site'
sites1<-merge(sites, tndoc, by='site')
sites<-sites1

#write.csv(sites, 'site summary.csv')

#okay ready watershed level dataset for correlations now


#is Pearson possible?

hist(sites$meanpCO2)
#needs transformation
#

hist(log10(sites$meanpCO2))
#maybe possible to use this - close to normal
sites$log10pCO2<-log10(sites$meanpCO2)
names(sites)
hist(sites$meanprecip_mm, xlim=c(0,1600))

#now this table is good.  time to bring in stream order and average DO
order<-read.csv('neon_stream_order.csv')
sites1<-merge(sites, order, by='site', all=TRUE)   ##########changed the all=TRUE here
sites<-sites1

#and bring in the DO
#run the DO_grab script first

names(aggregated_DOxSite)[1]<-'site'
#set_DO<-set[,c(2,4,5)][is.na(set$dailymean)==FALSE,]
#set_DO_site<-aggregate(set_DO$dailymean, by=list(site = set_DO$site), FUN=mean)
#names(set_DO_site)[2]<-'DOmean'
sites1<-merge(sites, aggregated_DOxSite, by='site')
sites<-sites1
