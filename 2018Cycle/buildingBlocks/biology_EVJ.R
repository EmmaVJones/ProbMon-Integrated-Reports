# File: Biology.R
# Purpose: Demonstrate creation of a design status file, 
# estimation of stream status, and estimation of biology stream condition CDFs for combined 1999 and 2005 survey
# Programmer: Tony Olson, modified for Virginia by Jason Hill to work with new R upgrade
# Modified by Emma Jones to work R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
# Date: April 10, 2018 
# This code can subset by basin, year, and ecoregion - then make CDF and percentiles for each subset

# Load all libraries
library(spsurvey) #v3.3
library(tidyverse) #v1.2.1
#library("scales")
#library("reshape2")
#library ("plyr")

#Load Rfuctions
# Purpose: R functions to be initiated at R startup
# Programmer: Tony Olsen
# Date: February 25, 2003
# To use copy and paste into R.

plot.map <- function (x, y) {
  rx <- range (x[!is.na(y)], na.rm=TRUE)
  ry <- range (y[!is.na(x)], na.rm=TRUE)
  plot.new ()
  plot.window (rx, ry, asp=1)
}

head <- function (x, n=10) {
  if (is.null(dim(x))) x[1:n] else x[1:n,]
}

tail <- function (x, n=10) {
  if (is.null(dim(x))) x[(length(x)-n+1):length(x)]
  else x[(nrow(x)-n+1):nrow(x),]
}

#### Set up sample frame stream length summary information
# Current example uses statewide Strahler Order 
# List stream kilometer by km

sframe <- c('1st'=51210, '2nd'=13680, '3rd'=7781.08, '4th'=4448.257, 
            '5th'=1731.302, '6th'=163.901, '7th'=14.7099 )

### Read design and status file and do a few summaries
statusdata <- read.delim('analysis/biology.tab',comment.char = "")
names(statusdata)

## Extract data for design status file. keep only sites evaluated
dsgnstatus <- select(statusdata,'sampleID','strahler.order','Longitude.DD','Latitude.DD','design.weight',
                     'weight.category', 'station', 'state','status', 'set','Basin','SubBasin','BayShed','BayPanel',
                     'EcoRegion','BioRegion','Year','Panel','BioPanel','Order','BasinSize','IR2008','IR2010',
                     'IR2012','IR2014','IR2016','StreamSizeCat','StreamSizeCatPhase') 
names(dsgnstatus) <- c('siteID','StrahlerOrder','LongDD','LatDD',
                       'InitialWeight','MDCaty','SiteNum','State','Status','Set','Basin','SubBasin',
                       'BayShed','BayPanel','EcoRegion','BioRegion','Year','Panel','BioPanel','Order',
                       'BasinSize','IR2008','IR2010','IR2012','IR2014','IR2016','StreamSizeCat','StreamSizeCatPhase')

# add marinus projection coordinates
tmp <- marinus(dsgnstatus$LatDD,-dsgnstatus$LongDD)
dsgnstatus$xmarinus <- tmp[,'x']
dsgnstatus$ymarinus <- tmp[,'y']

# Plot site locations
plot(dsgnstatus$xmarinus,dsgnstatus$ymarinus ,type='n')
box()
points(dsgnstatus$xmarinus,dsgnstatus$ymarinus,pch=16)
table(dsgnstatus$Status)
tst <- dsgnstatus$Status == 'TS'
points(dsgnstatus$xmarinus[tst],dsgnstatus$ymarinus[tst],pch=16,col='Green')

## Create weights for population estimation
# recode StrahlerOrder to match sframe 
dsgnstatus$StrahlerOrder <- as.factor(dsgnstatus$StrahlerOrder)
levels(dsgnstatus$StrahlerOrder) <- c('1st','2nd','3rd','4th','5th','6th','7th')
table(dsgnstatus$StrahlerOrder,dsgnstatus$MDCaty)
wgt <- sframe/table(dsgnstatus$StrahlerOrder)

dsgnstatus$finalwgt <- adjwgt(rep(TRUE,nrow(dsgnstatus)), 
                              dsgnstatus$InitialWeight,
                              dsgnstatus$StrahlerOrder, sframe)

write.table(dsgnstatus,'dsgnstatus.tab',sep='\t',row.names=FALSE)



####### Stream Extent and Status estimation
sites.ext <- data.frame(siteID=dsgnstatus$siteID, Use=rep(TRUE,nrow(dsgnstatus)) )

subpop.ext <- data.frame(siteID=dsgnstatus$siteID,Region=rep('Virginia',nrow(dsgnstatus)) )

design.ext <- data.frame(siteID=dsgnstatus$siteID,stratum=rep(1,nrow(dsgnstatus)),
                         wgt=dsgnstatus$finalwgt, xcoord=dsgnstatus$xmarinus,
                         ycoord=dsgnstatus$ymarinus)

StatusTNT <- dsgnstatus$Status
levels(StatusTNT)
levels(StatusTNT) <- list(T=c('TS','PD','OT'), NT=c('NT') )

data.cat.ext <- data.frame(dsgnstatus[,c('siteID','Status')],StatusTNT)

popstatus.est <- cat.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cat = data.cat.ext, 
                              conf=95, vartype="Local")

write.table(popstatus.est,'popstatus2.tab',sep='\t',row.names=FALSE)


#### Read in master field data
masterdata <- read.delim('analysis/masterdata.tab')
head(masterdata)
tail(masterdata)
names(masterdata)
# Note: all wadeable sampled target sites from 2001-2014
svydata <- masterdata[1:646,]

### Match dsgnstatus with masterdata
indx <- match(svydata$Station,dsgnstatus$siteID)
# does site status agree with having data available?
table(dsgnstatus$Status[indx])
cbind(dsgnstatus$siteID[indx],dsgnstatus$Status[indx])


### This section of code produces statewide estimates  
# merge two data sets
finaldata <- data.frame(dsgnstatus[indx,],svydata)
names(finaldata)
head(finaldata)
tail(finaldata)

write.table(finaldata,'finaldata.tab',sep='\t',row.names=FALSE)


## biology estimates
sites.ext <- data.frame(siteID=finaldata$siteID, Use=rep(TRUE,nrow(finaldata)) )

subpop.ext <- data.frame(siteID=finaldata$siteID,
                         Region=rep('Virginia',nrow(finaldata)) )

design.ext <- data.frame(siteID=finaldata$siteID,stratum=rep('1',nrow(finaldata)),
                         wgt=finaldata$finalwgt, xcoord=finaldata$xmarinus,
                         ycoord=finaldata$ymarinus)

data.cont.ext <- finaldata[,c('siteID','VSCIAll')]

biostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
                           pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
                           conf=95, vartype="Local")

write.table(biostatus$CDF,'biostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(biostatus$Pct,'biostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(biostatus$Tot,'biostatus.Tot.tab',sep='\t',row.names=FALSE)

# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(biostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='biostatusCDF.pdf',biostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='biostatusCDFkm.pdf',biostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')



# Now list all the ways Jason has to look at the data in biology.R
# Will need to repeat all of these (and maybe more) for each of: VSCI, DO, dissolved Metals, LRBS, 
# MetalsCCU, pH, SfClNaK, SpCond, TDS, TN, TotHab, TP
# - Statewide
# - Basin (Potomac-Shenandoah James, Rappahannock-York, Roanoke, Chowan, Tennessee, New) finaldata$Basin
# - SubBasin (Holston, Big Sandy, Clinch-Powell, Potomac, Shenandoah, Rappahannock, York) finaldata$SubBasin 
#   less the basins that were already run (James, Roanoke, Chowan, New)
#   Cheasapeake Bay not included bc too low n still
# - Ecoregion ( Piedmont, Northern Piedmont, Central Appalachian Ridges and Valleys, 
#   Blue Ridge Mountains, Southeastern Plains, Central Appalachians) finaldata$EcoRegion
#   Middle Atlantic Coastal Plain not included bc too low n still
# - Bioregion (Mountain, Piedmont, Coast) finaldata$BioRegion
# - Order (1:5) finaldata$Order *Note* 6 not run bc n too low
# - Basin Size (1 = less than 1 sq mile, 2 = 1 to 10 sq mile, 
#   3 = 10 to 200 sq mile, 4 = >200 sq mile) finaldata$BasinSize
# - Year (2001:2014/2016) finaldata$Year
# - Panel (Phase1 = 2001-2007, Phase2 = 2008-2014) finaldata$Panel
# - Bay vs NonBay (Bay, NonBay) finaldata$BayShed
# - Bay/NonBay by Panel (BayPhase1, BayPhase2, NonBayPhase1, NonBayPhase2) finaldata$BayPanel
# - Biophase (Phase1 = 2001-2003, Phase2 = 2004-2006, Phase3 = 2007-2010,
#   Phase4= 2011-2014, Phase5 = 2015-2016?) finaldata$BioPanel
# - IR window (IR2008 = 2001-2006, IR2010 = 2003-2008, IR2012 = 2005-2010, IR2014 = 2007-2012,
#   IR2016 = 2009-2014, IR2018 = 2011-2016) finaldata$IR2008, etc
# - StreamSize (Small = orders 1,2,3 & <10 sq miles, Medium = orders 1,2,3,4 & >10 but <50 sq miles,
#   Large = orders 3,4,5,6 & >50 sq miles) finaldata$StreamSizeCat
# - StreamSize and Years (Phase1Small = (orders 1,2,3 & <10 sq miles and phase 1 -2001-2007,
#   Phase2Small = orders 1,2,3 & <10 sq miles and phase 2 -2008-2014,
#   Phase1Medium = medium stream and phase 1 -2001-2007,
#   Phase2Medium = medium stream and phase 2 -2008-2014,
#   Phase1Large = large stream and phase 2 -2001-2007,
#   Phase2Large = large stream and phase 2 -2008-2014) finaldata$StreamSizeCatPhase

# So taking the biostatus$CDF from each subpopulation run and stacking them together (individually first
# for each parameter) then all parameters together, we can get the datasets that the IR report .Rmd is 
# expecting to look similar to C:/HardDriveBackup/R/IR2016/data/probmon_trend_status03152017.csv