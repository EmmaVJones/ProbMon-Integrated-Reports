suppressPackageStartupMessages(library(tidyverse))#1.2.1
suppressPackageStartupMessages(library(sf))#0.6-1
extrafont::loadfonts(device="win") # run this once each session, see https://cran.r-project.org/web/packages/extrafont/README.html for more info on the font_import() function required for "	font family not found in Windows font database" warning to go away
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(rgdal))
suppressPackageStartupMessages(library(micromap))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(docxtools))


# Data needed for maps (needs to be in report for figure generation vs in global file)
basinssmooth <- rgdal::readOGR('originalData','VAbasins_smoothNoChesPeeDee')# basin shapefile for micromaps
VAoutline <- rgdal::readOGR('originalData','Va_otlne') # for basic site maps
allProb <- rgdal::readOGR('processedData','allProb2018')# for basic site maps
# a little data management action
pre2010 <- subset(allProb, allProb@data$Year <= '2010')
post2010 <- subset(allProb, allProb@data$Year > '2010')

# Benthic Priorities layer
priorities <- rgdal:: readOGR('NMC', 'benthicPriorities') %>%
  spTransform(CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'))

# Vafrm stream layer
vafrm99 <- rgdal::readOGR('F:/GIS/ProbMonGIS','vafrm_99_05Albers') %>%
  spTransform(CRS('+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'))

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(VAoutline,border='wheat3',col='wheat1',ylim=c(37.4,38.5))
plot(vafrm99, col= 'blue',lwd = 0.25,add=T)
legend('left',legend=c('National Hydrography \nDataset 1:100k (1999)')
       ,bty='n',inset=0.01,lty = 1, lwd = 0.25,col=c('blue'))


par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(VAoutline,border='wheat3',col='wheat1',ylim=c(37.4,38.5))
plot(priorities, col= 'red',lwd = 0.25,add=T)
legend('left',legend=c('Benthic Impairments \n(2016 - 2022 Prioritization Effort) ')
       ,bty='n',inset=0.01,lty = 1, lwd = 0.25,col=c('red'))



par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(VAoutline,border='wheat3',col='wheat1',ylim=c(37.4,38.5))
plot(allProb,col='gray21',pch=19,cex=0.5,add=T)
legend('left',legend=c('Wadeable Sites 2001 - 2016 (n = 735) ')
       ,title='Legend',bty='n',inset=0.01,pch=c(19,19),cex=0.8
       ,col=c('grey21','red'))




par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(VAoutline,border='wheat3',col='wheat1',ylim=c(37.4,38.5))
plot(pre2010,col='gray21',pch=19,cex=0.5,add=T)
legend('left',legend=c('Wadeable Sites (2001 - 2010)')
       ,title='Legend',bty='n',inset=0.01,pch=c(19,19),cex=0.8
       ,col=c('grey21','red'))


par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(VAoutline,border='wheat3',col='wheat1',ylim=c(37.4,38.5))
plot(pre2010,col='gray21',pch=19,cex=0.5,add=T)
plot(post2010,col='red',pch=19,cex=0.5,add=T)
legend('left',legend=c('Wadeable Sites (2001 - 2010)', 'Wadeable Sites (2010 - 2016)')
       ,title='Legend',bty='n',inset=0.01,pch=c(19,19),cex=0.8
       ,col=c('grey21','red'))

# Bring in relative risk data
rr <- read.csv('processedData/relriskIR2018.csv') %>%
  mutate(Stressor=dplyr::recode(Stressor,"TotHabstatus"="Habitat Disturbance",
                                "TDSstatus"='Ionic Strength',
                                "TNstatus"='Total Nitrogen',
                                "MetalCCUstatus"='Cumulative Dissolved Metals',
                                "LRBSstatus"='Streambed Sedimentation',
                                "TPstatus"='Total Phosphorus'))# %>%

# RELATIVE RISK PLOT
rr$Stressor <- as.factor(c('Habitat Disturbance', 'Total Phosphorus','Ionic Strength','Cumulative Dissolved Metals','Streambed Sedimentation','Total Nitrogen'))
rr$Stressor <-  factor(rr$Stressor,
                       levels=c('Total Nitrogen','Total Phosphorus','Cumulative Dissolved Metals','Streambed Sedimentation','Ionic Strength','Habitat Disturbance'))




#rr$Stressor <- factor(rr$Stressor,levels=c('Total Phosphorus','Streambed Sedimentation',
#                                           'Cumulative Dissolved Metals','Total Nitrogen',
#                                           'Ionic Strength','Habitat Disturbance'))

ggplot(rr,aes(Stressor,Estimate,fill=Subpopulation,label=Estimate))+
  geom_bar(stat='identity',width=.5)+
  theme(aspect.ratio=4/6)+
  labs(x="Stressor", y = "Relative Risk")+ylim(-0.1,5)+
  scale_x_discrete(breaks=unique(rr$Stressor)#,labels=addline_format(rr$Stressor)
                   )+
  coord_flip()+ geom_hline(yintercept=1,linetype = "longdash")+
  geom_text(aes(label=paste(format(Estimate,digits=1,nsmall=1)),vjust=-0.5,hjust=-0.5),size=3,family='Arial')+
  geom_errorbar(aes(ymin=LCB95Pct, ymax=UCB95Pct), width=0.2)+ # error bars
  theme_minimal()+scale_fill_manual(values=c("#D55E00"))+
  theme(legend.position="none")+ # no legend, center title
  theme( text=element_text(family="Arial"),
         panel.grid.major.x = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),
         axis.line = element_line(colour = "black"))+ # no background grid
  theme(plot.margin=unit(c(0,1,0,0),'cm'))





# stream order

s1 <- rgdal::readOGR('NMC','streams11')
s2 <- rgdal::readOGR('NMC','streams2')
s3 <- rgdal::readOGR('NMC','streams3')


par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(s1,col='blue',pch=19,cex=0.5)
plot(s2,col='blue',pch=19,lwd=4,add=T)
plot(s3,col='blue',pch=19,lwd=9,add=T)

legend('left',legend=c('Wadeable Sites (2001 - 2010)')
       ,title='Legend',bty='n',inset=0.01,pch=c(19,19),cex=0.8
       ,col=c('grey21','red'))


