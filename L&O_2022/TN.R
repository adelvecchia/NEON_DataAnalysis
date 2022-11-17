library(neonUtilities)
library(glue)
library(httr)
library(jsonlite)
library(stringr)
library(readr)
library(feather)
library(tidyverse)
library(lubridate)
library(gridExtra)

###Downloading data from Macrosheds folder
#####Macrosheds functions
read_macrosheds_all <- function(file_path){

  neon_files <- list.files(file_path, full.names = TRUE, recursive = TRUE)
  file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
  file_names <- str_split_fixed(file_names, '[/]', n = Inf)

  final_list <- list()
  for(i in 1:length(site_name)){

    all_files <- tibble()

    path <- glue('{f}/{i}.feather',
                 f = file_path,
                 i = site_name[i])

    one_site <- read_feather(path)

    all_files <- rbind(all_files, one_site)
  }

  all_site_files <- list(all_files)

  names(all_site_files) <- data_files[i]

  final_list <- c(final_list, all_site_files)


  #  meta_data_files <-  inv_file_names[grepl('categoricalCodes|readme|validation|variables', inv_file_names)]

  for(i in 1:length(meta_data_files)){

    path <- glue('{f}/{s}/{i}.feather',
                 f = file_path,
                 i = meta_data_files[i],
                 s = site_names[1])

    meta_file <- read_feather(path)

    meta_list <- list(meta_file)
    names(meta_list) <- meta_data_files[i]

    final_list <- c(final_list, meta_list)
  }
  return(final_list)
}
#####

#####
setwd('C://Users/adelv/Dropbox/NEON')
#setwd('C:/Users/adelv/Dropbox/NEON/Macrosheds drive download/all_stream_chem/neon/neon/derived/stream_chemistry__ms001')

filepath<-'C:/Users/adelv/Dropbox/NEON/Macrosheds drive download/all_stream_chem/neon/neon/derived/stream_chemistry__ms001'

waterchem<-read_neon_feathers(filepath, by_site=FALSE)

list.files(filepath, full.names = TRUE, recursive=TRUE)
neon_files <- list.files(filepath, full.names = TRUE, recursive = TRUE)
file_names <- str_split_fixed(neon_files, '.feather', n = Inf)[,1]
file_names <- str_split_fixed(file_names, '[/]', n = Inf)
#specify index starting site names
site_name <- file_names[628:684]
#file location copied: C:\Users\adelv\Dropbox\NEON\Macrosheds drive download\all_stream_chem\neon\neon\derived\stream_chemistry__ms001

all_files<-tibble()

for(i in 1:length(site_name)){

  path <- glue('{f}/{i}.feather',
               f = filepath,
               i = site_name[i])

  one_site <- read_feather(path)

  all_files <- rbind(all_files, one_site)
}


allfiles<-data.frame(all_files)
#write.csv(allfiles, 'allfiles.csv')  #all macrosheds data now
unique(allfiles[,3])
#temp<-allfiles[allfiles$var=='IS_temp',]
tn<-allfiles[allfiles$var=='GN_TN',]
names(tn)[4]<-'TN'
names(tn)[5]<-'TNflag'
doc<-allfiles[ allfiles$var=='GN_DOC',]
names(doc)[4]<-'DOC'
names(doc)[5]<-'DOCflag'
unique(allfiles$site_name)
tn1<-merge(tn, doc, by=c('datetime', 'site_name'))

#####
tn1$date<-as.Date(tn1$datetime, format = '%Y-%m-%d %H:%M:%S', tz="")
tn<-tn1
#tn$datetime1<-as_datetime(temp$datetime, tz="UTC")
#tn$hour<-as.numeric(format(temp$datetime1, '%H'))
#tn$year<-year(temp$date)
tn<-tn[tn$DOCflag==0,] #lost 200 of 900 observations
tn<-tn[tn$TNflag==0,]
names(tn)[2]<-'siteID'

lakes=c('LIRO', 'TOOK', 'SUGG', 'CRAM', 'BARC', 'PRPO', 'PRLA')
tn<-tn[(tn$site %in% lakes) == FALSE,]

#pause here and make a summary by site

tndoc<-tn %>%
  group_by(siteID) %>%
  summarise(meanDOC = mean(DOC, na.rm=TRUE),
            meanTN = mean(TN, na.rm=TRUE),
            seDOC = se(DOC),
            seTN = se(TN))

#write.csv(tndoc, 'tndoc.csv')


#okay can use that for other things now

set_tn<-merge(tn, concs, by=c('siteID', 'date'), all=FALSE)
#why are there only 1four3four
#write.csv(set_tn, 'set_tn.csv')
#set_tn<-read.csv('set_tn.csv')

names(set_tn)

plot(N2Oaq1~TN, data=set_tn)
tnmod<-lm(N2Oaq1~val, data=set_tn)
summary(tnmod)
abline(tnmod)
#the one outlier has to come out

set_tn<-set_tn[set_tn$N2Oaq1<=1,]

plot(N2Oaq1~TN, data=set_tn)
tnmod<-lm(N2Oaq1~TN, data=set_tn)
summary(tnmod)
abline(tnmod)
#significant but explains barely any variation

hist(set_tn$N2Oaq1, breaks=100)
#mostly normal

#Use the individual measurements or summarize by catchment?
#I think TN and N2O, and then DOC and CO2, CHfour
#so two panels on those?

#run the Lauerwald style plots script first
mdbDOCmin<-min(methDB$DOC_actual_uM/1000, na.rm=TRUE)
mdbDOCmax<-max(methDB$DOC_actual_uM/1000, na.rm=TRUE)

tiff('doc_ghg.tif', width=4, height=2.2, res=390, units='in')
ggplot(set_tn) +
    geom_rect(aes(ymin=mdbmin, ymax=mdbmax,
                 xmin=log10(mdbDOCmin), xmax=log10(mdbDOCmax)), fill=NA, alpha=0.1, colour='coral2',
      linetype='solid', size=0.3)+
  geom_rect(aes(ymin=1.8, ymax=3.7, xmin=-1.33, xmax=0.31), fill=NA, alpha=0.3, colour='grey70',
            linetype='dotted', size=1.2)+
  geom_point(aes(x=log10(DOC), y=log10(pCO21)), size=1, colour = 'grey30') +
  geom_abline(aes(slope=0.07957, intercept=3.21849), color='grey30', linetype='solid') +
#  geom_errorbar(aes(x=DOC, ymin=log10(pCO21-pCO2_err1), ymax=log10(pCO21+pCO2_err1))) +
 # geom_errorbar(aes(x=DOC, ymin=log10(pCH41-pCH4_err1),
  #                  ymax=log10(pCH41+pCH4_err1)), alpha=0.8, colour='coral2') +
  geom_point(aes(x=log10(DOC), y=log10(pCH41)), size=1, colour = 'coral2') +
  geom_abline(aes(slope=0.95802, intercept=1.26505), color='coral2', linetype='solid') +
 # geom_errorbar(aes(x=meanairtemp, ymin=(meanpN2O-sepN2O),
  #                  ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
#  geom_point(aes(x=meanairtemp, y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
#  geom_abline(aes(slope=0.017224, intercept=0.546253), color='grey30', linetype='longdash') +
  theme_classic() +
#  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4, 1)) +
#  scale_x_continuous(limits=c(-4, 24), breaks = seq(-4, 24, 4)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='log10(DOC) [mg/L]')
dev.off()

docmod<-lm(log10(pCO21)~log10(DOC), data=set_tn)
summary(docmod)

docmod1<-lm(log10(pCH41)~log10(DOC), data=set_tn)
summary(docmod1)
#significant for methane

docmod2<-lm(pN2O1~log10(DOC), data=set_tn)
summary(docmod2)

#and maybe redo these plots with central tendency/per site?
#whats our argument for only looking at these?

hist(set_tn$TN, breaks=100)
tiff('tn_N2O.tif', width=4, height=2.2, res=390, units='in')
ggplot(set_tn_filt) +
 # geom_rect(aes(ymin=mdbmin, ymax=mdbmax,
  #              xmin=log10(mdbDOCmin), xmax=log10(mdbDOCmax)), fill=NA, alpha=0.1, colour='coral2',
   #         linetype='solid', size=0.7)+
 # geom_rect(aes(ymin=1.8, ymax=3.7, xmin=-1.33, xmax=0.31), fill=NA, alpha=0.3, colour='grey70',
  #          linetype='dotted', size=1.2)+
  geom_point(aes(x=log10(TN), y=pN2O1), size=2, shape=1, colour = 'grey30') +
#  geom_abline(aes(slope=0.07957, intercept=3.21849), color='grey30', linetype='solid') +
  #  geom_errorbar(aes(x=DOC, ymin=log10(pCO21-pCO2_err1), ymax=log10(pCO21+pCO2_err1))) +
  # geom_errorbar(aes(x=DOC, ymin=log10(pCH41-pCH4_err1),
  #                  ymax=log10(pCH41+pCH4_err1)), alpha=0.8, colour='coral2') +
 # geom_point(aes(x=log10(DOC), y=log10(pCH41)), size=1, colour = 'coral2') +
#  geom_abline(aes(slope=0.95802, intercept=1.26505), color='coral2', linetype='solid') +
  # geom_errorbar(aes(x=meanairtemp, ymin=(meanpN2O-sepN2O),
  #                  ymax=(meanpN2O+sepN2O)), alpha=0.8, colour='grey30') +
  #  geom_point(aes(x=meanairtemp, y=meanpN2O), size=2, shape=1,  colour = 'grey30') +
    geom_abline(aes(slope=0.11884, intercept=0.78262), color='grey30', linetype='longdash') +
  theme_classic() +
  #  scale_y_continuous(limits = c(0, 4), breaks = seq(0,4, 1)) +
  #  scale_x_continuous(limits=c(-4, 24), breaks = seq(-4, 24, 4)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='log10(TN) [mg/L]')
dev.off()

set_tn_filt<-set_tn[set_tn$TN>=0.0001,]
hist(set_tn_filt$TN, breaks=100)
hist(set_tn_filt$pN2O1, breaks=100)
N2O1tnmod2<-lm(pN2O1~log10(TN), data=set_tn_filt)
summary(tnmod2)

names(set_tn)

tiff('doc_ghg_sites.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
 # geom_rect(aes(ymin=mdbmin, ymax=mdbmax, #already log10 transformed
  #              xmin=log10(mdbDOCmin), xmax=log10(mdbDOCmax)), fill=NA, alpha=0.1, colour='coral2',
 #           linetype='solid', size=0.3)+
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
    scale_x_continuous(limits=c(-1.6, 1.6), breaks = seq(-2, 2, 1)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=0, hjust = 0),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_blank()) +
  labs(x='log10(DOC) [mg/L]')
dev.off()

docmod<-lm(log10(meanpCO2)~log10(meanDOC), data=sites)
summary(docmod)

docmod1<-lm(log10(meanpCH4)~log10(meanDOC), data=sites)
summary(docmod1)
#significant for methane

docmod2<-lm(log10(meanpN2O)~log10(meanDOC), data=sites)
summary(docmod2)
#NS for N2O

mdbTNmin<-min(methDB$TN_actual_uM/1000, na.rm=TRUE)
mdbTNmax<-max(methDB$TN_actual_uM/1000, na.rm=TRUE)


tiff('tn_ghg_sites.tif', width=4, height=2.2, res=390, units='in')
ggplot(sites) +
  # geom_rect(aes(ymin=mdbmin, ymax=mdbmax, #already log10 transformed
  #              xmin=log10(mdbTNmin), xmax=log10(mdbTNmax)), fill=NA, alpha=0.1, colour='coral2',
  #           linetype='solid', size=0.3)+
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

tnmod<-lm(log10(meanpCO2)~log10(meanTN), data=sites)
summary(tnmod)

tnmod1<-lm(log10(meanpCH4)~log10(meanTN), data=sites)
summary(tnmod1)
#significant for methane

tnmod2<-lm(meanpN2O~log10(meanTN), data=sites)
summary(tnmod2)
#NS for N2O
