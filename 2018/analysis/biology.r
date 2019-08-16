# File: Biology.R
# Purpose: Demonstrate creation of a design status file, 
# estimation of stream status, and estimation of biology stream condition CDFs for combined 1999 and 2005 survey
# Programmer: Tony Olson, modified for Virginia by Jason Hill to work with new R upgrade
# Date: February 20, 2017 to work in 'R' 2.15.3 and Psurvey Package 2.5 (spsurvey)
# This code can subset by basin, year, and ecoregion - then make CDF and percentiles for each subset

# Load all libraries
library("spsurvey")
library("sp")
library("ggplot2")
library("scales")
library("reshape2")
library ("plyr")

setwd("C:/Users/ktq89598/Desktop/ProbMon2006-2010/2014/IR2016ProbMonChapter/IRFinal/CDF/VSCI_QA_MAR2016")

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
statusdata <- read.delim('biology.tab',comment.char = "")
names(statusdata)

## Extract data for design status file. keep only sites evaluated
dsgnstatus <- statusdata[1:948,c('sampleID','strahler.order',
				'Longitude.DD','Latitude.DD','design.weight',
				'weight.category', 'station', 'state',
				'status', 'set','Basin','SubBasin','BayShed','BayPanel',
        'EcoRegion','BioRegion','Year','Panel',
        'BioPanel','Order','BasinSize','IR2008','IR2010','IR2012','IR2014','IR2016','StreamSizeCat','StreamSizeCatPhase')]
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
levels(dsgnstatus$StrahlerOrder) <- 
			c('1st','2nd','3rd','4th','5th','6th','7th')
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
masterdata <- read.delim('masterdata.tab')
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

#### Subset by major basin, allow VDEQ to evualuate conditions by basin
# subset by basin (roanoke)
finaldata2<-finaldata[finaldata$Basin=='Roanoke',]; #subset for roanoke basin;
names(finaldata2)
head(finaldata2)
tail(finaldata2)

write.table(finaldata2,'finaldata2.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata2$siteID, Use=rep(TRUE,nrow(finaldata2)) )

subpop.ext <- data.frame(siteID=finaldata2$siteID,
			Region=rep('Roanoke Basin',nrow(finaldata2)) )

design.ext <- data.frame(siteID=finaldata2$siteID,stratum=rep('1',nrow(finaldata2)),
			wgt=finaldata2$finalwgt, xcoord=finaldata2$xmarinus,
			ycoord=finaldata2$ymarinus)

data.cont.ext <- finaldata2[,c('siteID','VSCIAll')]

roanokestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(roanokestatus$CDF,'roanokestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(roanokestatus$Pct,'roanokestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(roanokestatus$Tot,'roanokestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(roanokestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='roanokestatusCDF.pdf',roanokestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='roanokestatusCDFkm.pdf',roanokestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')
                             
                 
# subset by basin (james)
finaldata3<-finaldata[finaldata$Basin=='James',]; #subset for james basin;
names(finaldata3)
head(finaldata3)
tail(finaldata3)

write.table(finaldata3,'finaldata3.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata3$siteID, Use=rep(TRUE,nrow(finaldata3)) )

subpop.ext <- data.frame(siteID=finaldata3$siteID,
			Region=rep('James Basin',nrow(finaldata3)) )

design.ext <- data.frame(siteID=finaldata3$siteID,stratum=rep('1',nrow(finaldata3)),
			wgt=finaldata3$finalwgt, xcoord=finaldata3$xmarinus,
			ycoord=finaldata3$ymarinus)

data.cont.ext <- finaldata3[,c('siteID','VSCIAll')]

jamesstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(jamesstatus$CDF,'jamesstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(jamesstatus$Pct,'jamesstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(jamesstatus$Tot,'jamesstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(jamesstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='jamesstatusCDF.pdf',jamesstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='jamesstatusCDFkm.pdf',jamesstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  

# subset by basin (potomac-shenandoah)
finaldata4<-finaldata[finaldata$Basin=='Potomac-Shenandoah',]; #subset for potomac-shenandoah basin;
names(finaldata4)
head(finaldata4)
tail(finaldata4)

write.table(finaldata4,'finaldata4.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata4$siteID, Use=rep(TRUE,nrow(finaldata4)) )

subpop.ext <- data.frame(siteID=finaldata4$siteID,
			Region=rep('Potomac-Shenandoah',nrow(finaldata4)) )

design.ext <- data.frame(siteID=finaldata4$siteID,stratum=rep('1',nrow(finaldata4)),
			wgt=finaldata4$finalwgt, xcoord=finaldata4$xmarinus,
			ycoord=finaldata4$ymarinus)

data.cont.ext <- finaldata4[,c('siteID','VSCIAll')]

potomacshenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(potomacshenstatus$CDF,'potomacshenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(potomacshenstatus$Pct,'potomacshenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(potomacshenstatus$Tot,'potomacshenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(potomacshenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='potomacshenstatusCDF.pdf',potomacshenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='potomacshenstatusCDFkm.pdf',potomacshenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  
                 
# subset by basin (rappahannock-york)
finaldata5<-finaldata[finaldata$Basin=='Rappahannock-York',]; #subset for rappahannock-york basin;
names(finaldata5)
head(finaldata5)
tail(finaldata5)

write.table(finaldata5,'finaldata5.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata5$siteID, Use=rep(TRUE,nrow(finaldata5)) )

subpop.ext <- data.frame(siteID=finaldata5$siteID,
			Region=rep('Rappahannock-York',nrow(finaldata5)) )

design.ext <- data.frame(siteID=finaldata5$siteID,stratum=rep('1',nrow(finaldata5)),
			wgt=finaldata5$finalwgt, xcoord=finaldata5$xmarinus,
			ycoord=finaldata5$ymarinus)

data.cont.ext <- finaldata5[,c('siteID','VSCIAll')]

rappyorkstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(rappyorkstatus$CDF,'rappyorkstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(rappyorkstatus$Pct,'rappyorkstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(rappyorkstatus$Tot,'rappyorkstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(rappyorkstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='rappyorkstatusCDF.pdf',rappyorkstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='rappyorkstatusCDFkm.pdf',rappyorkstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                   
                                                                                                  
# subset by basin (new)
finaldata6<-finaldata[finaldata$Basin=='New',]; #subset for new basin;
names(finaldata6)
head(finaldata6)
tail(finaldata6)

write.table(finaldata6,'finaldata6.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata6$siteID, Use=rep(TRUE,nrow(finaldata6)) )

subpop.ext <- data.frame(siteID=finaldata6$siteID,
			Region=rep('New',nrow(finaldata6)) )

design.ext <- data.frame(siteID=finaldata6$siteID,stratum=rep('1',nrow(finaldata6)),
			wgt=finaldata6$finalwgt, xcoord=finaldata6$xmarinus,
			ycoord=finaldata6$ymarinus)

data.cont.ext <- finaldata6[,c('siteID','VSCIAll')]

newstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(newstatus$CDF,'newstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(newstatus$Pct,'newstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(newstatus$Tot,'newstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(newstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='newstatusCDF.pdf',newstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='newstatusCDFkm.pdf',newstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                   


# subset by basin (chowan)
finaldata7<-finaldata[finaldata$Basin=='Chowan',]; #subset for chowan basin;
names(finaldata7)
head(finaldata7)
tail(finaldata7)

write.table(finaldata7,'finaldata7.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata7$siteID, Use=rep(TRUE,nrow(finaldata7)) )

subpop.ext <- data.frame(siteID=finaldata7$siteID,
			Region=rep('Chowan',nrow(finaldata7)) )

design.ext <- data.frame(siteID=finaldata7$siteID,stratum=rep('1',nrow(finaldata7)),
			wgt=finaldata7$finalwgt, xcoord=finaldata7$xmarinus,
			ycoord=finaldata7$ymarinus)

data.cont.ext <- finaldata7[,c('siteID','VSCIAll')]

chowanstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(chowanstatus$CDF,'chowanstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(chowanstatus$Pct,'chowanstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(chowanstatus$Tot,'chowanstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(chowanstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='chowanstatusCDF.pdf',chowanstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='chowanstatusCDFkm.pdf',chowanstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                   

# subset by basin (tennessee)
finaldata8<-finaldata[finaldata$Basin=='Tennessee',]; #subset for tennessee basin;
names(finaldata8)
head(finaldata8)
tail(finaldata8)

write.table(finaldata8,'finaldata8.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata8$siteID, Use=rep(TRUE,nrow(finaldata8)) )

subpop.ext <- data.frame(siteID=finaldata8$siteID,
			Region=rep('Tennessee',nrow(finaldata8)) )

design.ext <- data.frame(siteID=finaldata8$siteID,stratum=rep('1',nrow(finaldata8)),
			wgt=finaldata8$finalwgt, xcoord=finaldata8$xmarinus,
			ycoord=finaldata8$ymarinus)

data.cont.ext <- finaldata8[,c('siteID','VSCIAll')]

tennstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(tennstatus$CDF,'tennstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(tennstatus$Pct,'tennstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(tennstatus$Tot,'tennstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(tennstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='tennstatusCDF.pdf',tennstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='tennstatusCDFkm.pdf',tennstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                   


#### Subset by subbasin, the subbasin has lower sample n, but show interesting patterns within basin
# subset by subbasin (tennessee - big sandy, holston, clinch-powell)
finaldata9<-finaldata[finaldata$SubBasin=='Holston',]; #subset for tennessee subbasin;
names(finaldata9)
head(finaldata9)
tail(finaldata9)

write.table(finaldata9,'finaldata9.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata9$siteID, Use=rep(TRUE,nrow(finaldata9)) )

subpop.ext <- data.frame(siteID=finaldata9$siteID,
			Region=rep('Holston',nrow(finaldata9)) )

design.ext <- data.frame(siteID=finaldata9$siteID,stratum=rep('1',nrow(finaldata9)),
			wgt=finaldata9$finalwgt, xcoord=finaldata9$xmarinus,
			ycoord=finaldata9$ymarinus)

data.cont.ext <- finaldata9[,c('siteID','VSCIAll')]

holstonstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(holstonstatus$CDF,'holstonstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(holstonstatus$Pct,'holstonstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(holstonstatus$Tot,'holstonstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(holstonstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='holstonstatusCDF.pdf',holstonstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='holstonstatusCDFkm.pdf',holstonstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  
                 
# subset by subbasin (tennessee - big sandy, holston, clinch-powell)
finaldata10<-finaldata[finaldata$SubBasin=='Big Sandy',]; #subset for big sandy subbasin;
names(finaldata10)
head(finaldata10)
tail(finaldata10)

write.table(finaldata10,'finaldata10.tab',sep='\t',row.names=FALSE)

## biology estimate
sites.ext <- data.frame(siteID=finaldata10$siteID, Use=rep(TRUE,nrow(finaldata10)) )

subpop.ext <- data.frame(siteID=finaldata10$siteID,
			Region=rep('Big Sandy',nrow(finaldata10)) )

design.ext <- data.frame(siteID=finaldata10$siteID,stratum=rep('1',nrow(finaldata10)),
			wgt=finaldata10$finalwgt, xcoord=finaldata10$xmarinus,
			ycoord=finaldata10$ymarinus)

data.cont.ext <- finaldata10[,c('siteID','VSCIAll')]

bigsandystatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(bigsandystatus$CDF,'bigsandystatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(bigsandystatus$Pct,'bigsandystatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(bigsandystatus$Tot,'bigsandystatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(bigsandystatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='bigsandystatusCDF.pdf',bigsandystatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='bigsandystatusCDFkm.pdf',bigsandystatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)') 
                 
# subset by subbasin (tennessee - big sandy, holston, clinch-powell)
finaldata11<-finaldata[finaldata$SubBasin=='Clinch-Powell',]; #subset for clinch-powell subbasin;
names(finaldata11)
head(finaldata11)
tail(finaldata11)

write.table(finaldata11,'finaldata11.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata11$siteID, Use=rep(TRUE,nrow(finaldata11)) )

subpop.ext <- data.frame(siteID=finaldata11$siteID,
			Region=rep('Clinch-Powell',nrow(finaldata11)) )

design.ext <- data.frame(siteID=finaldata11$siteID,stratum=rep('1',nrow(finaldata11)),
			wgt=finaldata11$finalwgt, xcoord=finaldata11$xmarinus,
			ycoord=finaldata11$ymarinus)

data.cont.ext <- finaldata11[,c('siteID','VSCIAll')]

clinchstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(clinchstatus$CDF,'clinchstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(clinchstatus$Pct,'clinchstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(clinchstatus$Tot,'clinchstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(clinchstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='clinchstatusCDF.pdf',clinchstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='clinchstatusCDFkm.pdf',clinchstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')   
                 

# subset by subbasin (potomac and shenandoah)
finaldata12<-finaldata[finaldata$SubBasin=='Potomac',]; #subset for potomac subbasin;
names(finaldata12)
head(finaldata12)
tail(finaldata12)

write.table(finaldata12,'finaldata12.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata12$siteID, Use=rep(TRUE,nrow(finaldata12)) )

subpop.ext <- data.frame(siteID=finaldata12$siteID,
			Region=rep('Potomac',nrow(finaldata12)) )

design.ext <- data.frame(siteID=finaldata12$siteID,stratum=rep('1',nrow(finaldata12)),
			wgt=finaldata12$finalwgt, xcoord=finaldata12$xmarinus,
			ycoord=finaldata12$ymarinus)

data.cont.ext <- finaldata12[,c('siteID','VSCIAll')]

potomacstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(potomacstatus$CDF,'potomacstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(potomacstatus$Pct,'potomacstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(potomacstatus$Tot,'potomacstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(potomacstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='potomacstatusCDF.pdf',potomacstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='potomacstatusCDFkm.pdf',potomacstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')     
                 
# subset by subbasin (potomac and shenandoah)
finaldata13<-finaldata[finaldata$SubBasin=='Shenandoah',]; #subset for shenandoah subbasin;
names(finaldata13)
head(finaldata13)
tail(finaldata13)

write.table(finaldata13,'finaldata13.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata13$siteID, Use=rep(TRUE,nrow(finaldata13)) )

subpop.ext <- data.frame(siteID=finaldata13$siteID,
			Region=rep('Shenandoah',nrow(finaldata13)) )

design.ext <- data.frame(siteID=finaldata13$siteID,stratum=rep('1',nrow(finaldata13)),
			wgt=finaldata13$finalwgt, xcoord=finaldata13$xmarinus,
			ycoord=finaldata13$ymarinus)

data.cont.ext <- finaldata13[,c('siteID','VSCIAll')]

shenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(shenstatus$CDF,'shenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(shenstatus$Pct,'shenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(shenstatus$Tot,'shenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(shenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='shenstatusCDF.pdf',shenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='shenstatusCDFkm.pdf',shenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')   
                 
# subset by subbasin (rappahanock and york)
finaldata14<-finaldata[finaldata$SubBasin=='Rappahannock',]; #subset for rappahannock subbasin;
names(finaldata14)
head(finaldata14)
tail(finaldata14)

write.table(finaldata14,'finaldata14.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata14$siteID, Use=rep(TRUE,nrow(finaldata14)) )

subpop.ext <- data.frame(siteID=finaldata14$siteID,
			Region=rep('Rappahanock',nrow(finaldata14)) )

design.ext <- data.frame(siteID=finaldata14$siteID,stratum=rep('1',nrow(finaldata14)),
			wgt=finaldata14$finalwgt, xcoord=finaldata14$xmarinus,
			ycoord=finaldata14$ymarinus)

data.cont.ext <- finaldata14[,c('siteID','VSCIAll')]

rappstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(rappstatus$CDF,'rappstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(rappstatus$Pct,'rappstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(rappstatus$Tot,'rappstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(rappstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='rappstatusCDF.pdf',rappstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='rappstatusCDFkm.pdf',rappstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  
                 
# subset by subbasin (rappahanock and york)
finaldata15<-finaldata[finaldata$SubBasin=='York',]; #subset for york subbasin;
names(finaldata15)
head(finaldata15)
tail(finaldata15)

write.table(finaldata15,'finaldata15.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata15$siteID, Use=rep(TRUE,nrow(finaldata15)) )

subpop.ext <- data.frame(siteID=finaldata15$siteID,
			Region=rep('York',nrow(finaldata15)) )

design.ext <- data.frame(siteID=finaldata15$siteID,stratum=rep('1',nrow(finaldata15)),
			wgt=finaldata15$finalwgt, xcoord=finaldata15$xmarinus,
			ycoord=finaldata15$ymarinus)

data.cont.ext <- finaldata15[,c('siteID','VSCIAll')]

yorkstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yorkstatus$CDF,'yorkstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yorkstatus$Pct,'yorkstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yorkstatus$Tot,'yorkstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yorkstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yorkstatusCDF.pdf',yorkstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yorkstatusCDFkm.pdf',yorkstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                      
                 
# subset by ecoregion (piedmont)
finaldata16<-finaldata[finaldata$EcoRegion=='Piedmont',]; #subset for piedmont ecoregion;
names(finaldata16)
head(finaldata16)
tail(finaldata16)

write.table(finaldata16,'finaldata16.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata16$siteID, Use=rep(TRUE,nrow(finaldata16)) )

subpop.ext <- data.frame(siteID=finaldata16$siteID,
			Region=rep('Piedmont',nrow(finaldata16)) )

design.ext <- data.frame(siteID=finaldata16$siteID,stratum=rep('1',nrow(finaldata16)),
			wgt=finaldata16$finalwgt, xcoord=finaldata16$xmarinus,
			ycoord=finaldata16$ymarinus)

data.cont.ext <- finaldata16[,c('siteID','VSCIAll')]

piedmontstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(piedmontstatus$CDF,'piedmontstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(piedmontstatus$Pct,'piedmontstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(piedmontstatus$Tot,'piedmontstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(piedmontstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='piedmontstatusCDF.pdf',piedmontstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='piedmontstatusCDFkm.pdf',piedmontstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                    


# subset by ecoregion (northern piedmont)
finaldata17<-finaldata[finaldata$EcoRegion=='Northern Piedmont',]; #subset for northern piedmont ecoregion;
names(finaldata17)
head(finaldata17)
tail(finaldata17)

write.table(finaldata17,'finaldata17.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata17$siteID, Use=rep(TRUE,nrow(finaldata17)) )

subpop.ext <- data.frame(siteID=finaldata17$siteID,
			Region=rep('Northern Piedmont',nrow(finaldata17)) )

design.ext <- data.frame(siteID=finaldata17$siteID,stratum=rep('1',nrow(finaldata17)),
			wgt=finaldata17$finalwgt, xcoord=finaldata17$xmarinus,
			ycoord=finaldata17$ymarinus)

data.cont.ext <- finaldata17[,c('siteID','VSCIAll')]

npiedmontstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(npiedmontstatus$CDF,'npiedmontstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(npiedmontstatus$Pct,'npiedmontstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(npiedmontstatus$Tot,'npiedmontstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(npiedmontstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='npiedmontstatusCDF.pdf',npiedmontstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='npiedmontstatusCDFkm.pdf',npiedmontstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)') 
                 
# subset by ecoregion (central appalachian ridge and valleys)
finaldata18<-finaldata[finaldata$EcoRegion=='Central Appalachian Ridges and Valleys',]; #subset for central appalachian ridge and valleys ecoregion;
names(finaldata18)
head(finaldata18)
tail(finaldata18)

write.table(finaldata18,'finaldata18.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata18$siteID, Use=rep(TRUE,nrow(finaldata18)) )

subpop.ext <- data.frame(siteID=finaldata18$siteID,
			Region=rep('Central Appalachian Ridges and Valleys',nrow(finaldata18)) )

design.ext <- data.frame(siteID=finaldata18$siteID,stratum=rep('1',nrow(finaldata18)),
			wgt=finaldata18$finalwgt, xcoord=finaldata18$xmarinus,
			ycoord=finaldata18$ymarinus)

data.cont.ext <- finaldata18[,c('siteID','VSCIAll')]

carvstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(carvstatus$CDF,'carvstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(carvstatus$Pct,'carvstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(carvstatus$Tot,'carvstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(carvstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='carvstatusCDF.pdf',carvstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='carvstatusCDFkm.pdf',carvstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                  
                 

# subset by ecoregion (Southeastern Plains)
finaldata19<-finaldata[finaldata$EcoRegion=='Southeastern Plains',]; #subset for Southeastern Plains ecoregion;
names(finaldata19)
head(finaldata19)
tail(finaldata19)

write.table(finaldata19,'finaldata19.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata19$siteID, Use=rep(TRUE,nrow(finaldata19)) )

subpop.ext <- data.frame(siteID=finaldata19$siteID,
			Region=rep('Southeastern Plains',nrow(finaldata19)) )

design.ext <- data.frame(siteID=finaldata19$siteID,stratum=rep('1',nrow(finaldata19)),
			wgt=finaldata19$finalwgt, xcoord=finaldata19$xmarinus,
			ycoord=finaldata19$ymarinus)

data.cont.ext <- finaldata19[,c('siteID','VSCIAll')]

splainsstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(splainsstatus$CDF,'splainsstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(splainsstatus$Pct,'splainsstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(splainsstatus$Tot,'splainsstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(splainsstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='splainsstatusCDF.pdf',splainsstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='splainsstatusCDFkm.pdf',splainsstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')          
                 

# subset by ecoregion (Blue Ridge Mountains)
finaldata20<-finaldata[finaldata$EcoRegion=='Blue Ridge Mountains',]; #subset for Blue Ridge Mountains ecoregion;
names(finaldata20)
head(finaldata20)
tail(finaldata20)

write.table(finaldata20,'finaldata20.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata20$siteID, Use=rep(TRUE,nrow(finaldata20)) )

subpop.ext <- data.frame(siteID=finaldata20$siteID,
			Region=rep('Blue Ridge Mountains',nrow(finaldata20)) )

design.ext <- data.frame(siteID=finaldata20$siteID,stratum=rep('1',nrow(finaldata20)),
			wgt=finaldata20$finalwgt, xcoord=finaldata20$xmarinus,
			ycoord=finaldata20$ymarinus)

data.cont.ext <- finaldata20[,c('siteID','VSCIAll')]

brmstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(brmstatus$CDF,'brmstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(brmstatus$Pct,'brmstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(brmstatus$Tot,'brmstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(brmstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='brmstatusCDF.pdf',brmstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='brmstatusCDFkm.pdf',brmstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  
                 

# subset by ecoregion (Central Appalachians)
finaldata21<-finaldata[finaldata$EcoRegion=='Central Appalachians',]; #subset for Central Appalachians ecoregion;
names(finaldata21)
head(finaldata21)
tail(finaldata21)

write.table(finaldata21,'finaldata21.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata21$siteID, Use=rep(TRUE,nrow(finaldata21)) )

subpop.ext <- data.frame(siteID=finaldata21$siteID,
			Region=rep('Central Appalachians',nrow(finaldata21)) )

design.ext <- data.frame(siteID=finaldata21$siteID,stratum=rep('1',nrow(finaldata21)),
			wgt=finaldata21$finalwgt, xcoord=finaldata21$xmarinus,
			ycoord=finaldata21$ymarinus)

data.cont.ext <- finaldata21[,c('siteID','VSCIAll')]

centappstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(centappstatus$CDF,'centappstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(centappstatus$Pct,'centappstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(centappstatus$Tot,'centappstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(centappstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='centappstatusCDF.pdf',centappstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='centappstatusCDFkm.pdf',centappstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                   
                 


# subset by bioregion (Mountain)
finaldata22<-finaldata[finaldata$BioRegion=='Mountain',]; #subset for Mountain bioregion;
names(finaldata22)
head(finaldata22)
tail(finaldata22)

write.table(finaldata22,'finaldata22.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata22$siteID, Use=rep(TRUE,nrow(finaldata22)) )

subpop.ext <- data.frame(siteID=finaldata22$siteID,
			Region=rep('Mountain Bioregion',nrow(finaldata22)) )

design.ext <- data.frame(siteID=finaldata22$siteID,stratum=rep('1',nrow(finaldata22)),
			wgt=finaldata22$finalwgt, xcoord=finaldata22$xmarinus,
			ycoord=finaldata22$ymarinus)

data.cont.ext <- finaldata22[,c('siteID','VSCIAll')]

mountainstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(mountainstatus$CDF,'mountainstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(mountainstatus$Pct,'mountainstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(mountainstatus$Tot,'mountainstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(mountainstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='mountainstatusCDF.pdf',mountainstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='mountainstatusCDFkm.pdf',mountainstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)') 


# subset by bioregion (Piedmont)
finaldata23<-finaldata[finaldata$BioRegion=='Piedmont',]; #subset for Piedmont bioregion;
names(finaldata23)
head(finaldata23)
tail(finaldata23)

write.table(finaldata23,'finaldata23.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata23$siteID, Use=rep(TRUE,nrow(finaldata23)) )

subpop.ext <- data.frame(siteID=finaldata23$siteID,
			Region=rep('Piedmont Bioregion',nrow(finaldata23)) )

design.ext <- data.frame(siteID=finaldata23$siteID,stratum=rep('1',nrow(finaldata23)),
			wgt=finaldata23$finalwgt, xcoord=finaldata23$xmarinus,
			ycoord=finaldata23$ymarinus)

data.cont.ext <- finaldata23[,c('siteID','VSCIAll')]

peidstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(peidstatus$CDF,'peidstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(peidstatus$Pct,'peidstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(peidstatus$Tot,'peidstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(peidstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='peidstatusCDF.pdf',peidstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='peidstatusCDFkm.pdf',peidstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')    


# subset by bioregion (Coast)
finaldata24<-finaldata[finaldata$BioRegion=='Coast',]; #subset for Coast bioregion;
names(finaldata24)
head(finaldata24)
tail(finaldata24)

write.table(finaldata24,'finaldata24.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata24$siteID, Use=rep(TRUE,nrow(finaldata24)) )

subpop.ext <- data.frame(siteID=finaldata24$siteID,
			Region=rep('Coastal Bioregion',nrow(finaldata24)) )

design.ext <- data.frame(siteID=finaldata24$siteID,stratum=rep('1',nrow(finaldata24)),
			wgt=finaldata24$finalwgt, xcoord=finaldata24$xmarinus,
			ycoord=finaldata24$ymarinus)

data.cont.ext <- finaldata24[,c('siteID','VSCIAll')]

coaststatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(coaststatus$CDF,'coaststatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(coaststatus$Pct,'coaststatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(coaststatus$Tot,'coaststatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(coaststatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='coaststatusCDF.pdf',coaststatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='coaststatusCDFkm.pdf',coaststatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by order (First Order)
finaldata25<-finaldata[finaldata$Order=='1',]; #subset for First Order;
names(finaldata25)
head(finaldata25)
tail(finaldata25)

write.table(finaldata25,'finaldata25.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata25$siteID, Use=rep(TRUE,nrow(finaldata25)) )

subpop.ext <- data.frame(siteID=finaldata25$siteID,
			Region=rep('First Order',nrow(finaldata25)) )

design.ext <- data.frame(siteID=finaldata25$siteID,stratum=rep('1',nrow(finaldata25)),
			wgt=finaldata25$finalwgt, xcoord=finaldata25$xmarinus,
			ycoord=finaldata25$ymarinus)

data.cont.ext <- finaldata25[,c('siteID','VSCIAll')]

firstorderstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(firstorderstatus$CDF,'firstorderstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(firstorderstatus$Pct,'firstorderstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(firstorderstatus$Tot,'firstorderstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(firstorderstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='firstorderstatusCDF.pdf',firstorderstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='firstorderstatusCDFkm.pdf',firstorderstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by order (Second Order)
finaldata26<-finaldata[finaldata$Order=='2',]; #subset for Second Order;
names(finaldata26)
head(finaldata26)
tail(finaldata26)

write.table(finaldata26,'finaldata26.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata26$siteID, Use=rep(TRUE,nrow(finaldata26)) )

subpop.ext <- data.frame(siteID=finaldata26$siteID,
			Region=rep('Second Order',nrow(finaldata26)) )

design.ext <- data.frame(siteID=finaldata26$siteID,stratum=rep('1',nrow(finaldata26)),
			wgt=finaldata26$finalwgt, xcoord=finaldata26$xmarinus,
			ycoord=finaldata26$ymarinus)

data.cont.ext <- finaldata26[,c('siteID','VSCIAll')]

secondorderstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(secondorderstatus$CDF,'secondorderstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(secondorderstatus$Pct,'secondorderstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(secondorderstatus$Tot,'secondorderstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(secondorderstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='secondorderstatusCDF.pdf',secondorderstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='secondorderstatusCDFkm.pdf',secondorderstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by order (Third Order)
finaldata27<-finaldata[finaldata$Order=='3',]; #subset for Third Order;
names(finaldata27)
head(finaldata27)
tail(finaldata27)

write.table(finaldata27,'finaldata27.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata27$siteID, Use=rep(TRUE,nrow(finaldata27)) )

subpop.ext <- data.frame(siteID=finaldata27$siteID,
			Region=rep('Third Order',nrow(finaldata27)) )

design.ext <- data.frame(siteID=finaldata27$siteID,stratum=rep('1',nrow(finaldata27)),
			wgt=finaldata27$finalwgt, xcoord=finaldata27$xmarinus,
			ycoord=finaldata27$ymarinus)

data.cont.ext <- finaldata27[,c('siteID','VSCIAll')]

thirdorderstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(thirdorderstatus$CDF,'thirdorderstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(thirdorderstatus$Pct,'thirdorderstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(thirdorderstatus$Tot,'thirdorderstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(thirdorderstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='thirdorderstatusCDF.pdf',thirdorderstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='thirdorderstatusCDFkm.pdf',thirdorderstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by order (Fourth Order)
finaldata28<-finaldata[finaldata$Order=='4',]; #subset for Fourth Order;
names(finaldata28)
head(finaldata28)
tail(finaldata28)

write.table(finaldata28,'finaldata28.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata28$siteID, Use=rep(TRUE,nrow(finaldata28)) )

subpop.ext <- data.frame(siteID=finaldata28$siteID,
			Region=rep('Fourth Order',nrow(finaldata28)) )

design.ext <- data.frame(siteID=finaldata28$siteID,stratum=rep('1',nrow(finaldata28)),
			wgt=finaldata28$finalwgt, xcoord=finaldata28$xmarinus,
			ycoord=finaldata28$ymarinus)

data.cont.ext <- finaldata28[,c('siteID','VSCIAll')]

fourthorderstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(fourthorderstatus$CDF,'fourthorderstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(fourthorderstatus$Pct,'fourthorderstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(fourthorderstatus$Tot,'fourthorderstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(fourthorderstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='fourthorderstatusCDF.pdf',fourthorderstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='fourthorderstatusCDFkm.pdf',fourthorderstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           


# subset by order (Fifth Order)
finaldata29<-finaldata[finaldata$Order=='5',]; #subset for Fifth Order;
names(finaldata29)
head(finaldata29)
tail(finaldata29)

write.table(finaldata29,'finaldata29.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata29$siteID, Use=rep(TRUE,nrow(finaldata29)) )

subpop.ext <- data.frame(siteID=finaldata29$siteID,
			Region=rep('Fifth Order',nrow(finaldata29)) )

design.ext <- data.frame(siteID=finaldata29$siteID,stratum=rep('1',nrow(finaldata29)),
			wgt=finaldata29$finalwgt, xcoord=finaldata29$xmarinus,
			ycoord=finaldata29$ymarinus)

data.cont.ext <- finaldata29[,c('siteID','VSCIAll')]

fifthorderstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(fifthorderstatus$CDF,'fifthorderstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(fifthorderstatus$Pct,'fifthorderstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(fifthorderstatus$Tot,'fifthorderstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(fifthorderstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='fifthorderstatusCDF.pdf',fifthorderstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='fifthorderstatusCDFkm.pdf',fifthorderstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by basin size (less than 1 sq mile)
finaldata30<-finaldata[finaldata$BasinSize=='1',]; #subset for less than 1 sq mile;
names(finaldata30)
head(finaldata30)
tail(finaldata30)

write.table(finaldata30,'finaldata30.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata30$siteID, Use=rep(TRUE,nrow(finaldata30)) )

subpop.ext <- data.frame(siteID=finaldata30$siteID,
			Region=rep('<1 square mile',nrow(finaldata30)) )

design.ext <- data.frame(siteID=finaldata30$siteID,stratum=rep('1',nrow(finaldata30)),
			wgt=finaldata30$finalwgt, xcoord=finaldata30$xmarinus,
			ycoord=finaldata30$ymarinus)

data.cont.ext <- finaldata30[,c('siteID','VSCIAll')]

basinonestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(basinonestatus$CDF,'basinonestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(basinonestatus$Pct,'basinonestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(basinonestatus$Tot,'basinonestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(basinonestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='basinonestatusCDF.pdf',basinonestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='basinonestatusCDFkm.pdf',basinonestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by basin size (1 to 10 sq mile)
finaldata31<-finaldata[finaldata$BasinSize=='2',]; #subset for 1 to 10 sq mile;
names(finaldata31)
head(finaldata31)
tail(finaldata31)

write.table(finaldata31,'finaldata31.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata31$siteID, Use=rep(TRUE,nrow(finaldata31)) )

subpop.ext <- data.frame(siteID=finaldata31$siteID,
			Region=rep('1 to 10 square mile',nrow(finaldata31)) )

design.ext <- data.frame(siteID=finaldata31$siteID,stratum=rep('1',nrow(finaldata31)),
			wgt=finaldata31$finalwgt, xcoord=finaldata31$xmarinus,
			ycoord=finaldata31$ymarinus)

data.cont.ext <- finaldata31[,c('siteID','VSCIAll')]

basintwostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(basintwostatus$CDF,'basintwostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(basintwostatus$Pct,'basintwostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(basintwostatus$Tot,'basintwostatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(basintwostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='basintwostatusCDF.pdf',basintwostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='basintwostatusCDFkm.pdf',basintwostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by basin size (10 to 200 sq mile)
finaldata32<-finaldata[finaldata$BasinSize=='3',]; #subset for 10 to 200 sq mile;
names(finaldata32)
head(finaldata32)
tail(finaldata32)

write.table(finaldata32,'finaldata32.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata32$siteID, Use=rep(TRUE,nrow(finaldata32)) )

subpop.ext <- data.frame(siteID=finaldata32$siteID,
			Region=rep('10 to 200 square mile',nrow(finaldata32)) )

design.ext <- data.frame(siteID=finaldata32$siteID,stratum=rep('1',nrow(finaldata32)),
			wgt=finaldata32$finalwgt, xcoord=finaldata32$xmarinus,
			ycoord=finaldata32$ymarinus)

data.cont.ext <- finaldata32[,c('siteID','VSCIAll')]

basinthreestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(basinthreestatus$CDF,'basinthreestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(basinthreestatus$Pct,'basinthreestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(basinthreestatus$Tot,'basinthreestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(basinthreestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='basinthreestatusCDF.pdf',basinthreestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='basinthreestatusCDFkm.pdf',basinthreestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by basin size (>200 sq mile)
finaldata33<-finaldata[finaldata$BasinSize=='4',]; #subset for >200 sq mile;
names(finaldata33)
head(finaldata33)
tail(finaldata33)

write.table(finaldata33,'finaldata33.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata33$siteID, Use=rep(TRUE,nrow(finaldata33)) )

subpop.ext <- data.frame(siteID=finaldata33$siteID,
			Region=rep('>200 square mile',nrow(finaldata33)) )

design.ext <- data.frame(siteID=finaldata33$siteID,stratum=rep('1',nrow(finaldata33)),
			wgt=finaldata33$finalwgt, xcoord=finaldata33$xmarinus,
			ycoord=finaldata33$ymarinus)

data.cont.ext <- finaldata33[,c('siteID','VSCIAll')]

basinfourstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(basinfourstatus$CDF,'basinfourstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(basinfourstatus$Pct,'basinfourstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(basinfourstatus$Tot,'basinfourstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(basinfourstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='basinfourstatusCDF.pdf',basinfourstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='basinfourstatusCDFkm.pdf',basinfourstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           


# subset by year (2001)
finaldata34<-finaldata[finaldata$Year=='2001',]; #subset for year 2001;
names(finaldata34)
head(finaldata34)
tail(finaldata34)

write.table(finaldata34,'finaldata34.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata34$siteID, Use=rep(TRUE,nrow(finaldata34)) )

subpop.ext <- data.frame(siteID=finaldata34$siteID,
			Region=rep('Year 2001',nrow(finaldata34)) )

design.ext <- data.frame(siteID=finaldata34$siteID,stratum=rep('1',nrow(finaldata34)),
			wgt=finaldata34$finalwgt, xcoord=finaldata34$xmarinus,
			ycoord=finaldata34$ymarinus)

data.cont.ext <- finaldata34[,c('siteID','VSCIAll')]

yearonestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearonestatus$CDF,'yearonestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearonestatus$Pct,'yearonestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearonestatus$Tot,'yearonestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearonestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearonestatusCDF.pdf',yearonestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearonestatusCDFkm.pdf',yearonestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by year (2002)
finaldata35<-finaldata[finaldata$Year=='2002',]; #subset for year 2002;
names(finaldata35)
head(finaldata35)
tail(finaldata35)

write.table(finaldata35,'finaldata35.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata35$siteID, Use=rep(TRUE,nrow(finaldata35)) )

subpop.ext <- data.frame(siteID=finaldata35$siteID,
			Region=rep('Year 2002',nrow(finaldata35)) )

design.ext <- data.frame(siteID=finaldata35$siteID,stratum=rep('1',nrow(finaldata35)),
			wgt=finaldata35$finalwgt, xcoord=finaldata35$xmarinus,
			ycoord=finaldata35$ymarinus)

data.cont.ext <- finaldata35[,c('siteID','VSCIAll')]

yeartwostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yeartwostatus$CDF,'yeartwostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yeartwostatus$Pct,'yeartwostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yeartwostatus$Tot,'yeartwostatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yeartwostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yeartwostatusCDF.pdf',yeartwostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yeartwostatusCDFkm.pdf',yeartwostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by year (2003)
finaldata36<-finaldata[finaldata$Year=='2003',]; #subset for year 2003;
names(finaldata36)
head(finaldata36)
tail(finaldata36)

write.table(finaldata36,'finaldata36.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata36$siteID, Use=rep(TRUE,nrow(finaldata36)) )

subpop.ext <- data.frame(siteID=finaldata36$siteID,
			Region=rep('Year 2003',nrow(finaldata36)) )

design.ext <- data.frame(siteID=finaldata36$siteID,stratum=rep('1',nrow(finaldata36)),
			wgt=finaldata36$finalwgt, xcoord=finaldata36$xmarinus,
			ycoord=finaldata36$ymarinus)

data.cont.ext <- finaldata36[,c('siteID','VSCIAll')]

yearthreestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearthreestatus$CDF,'yearthreestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearthreestatus$Pct,'yearthreestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearthreestatus$Tot,'yearthreestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearthreestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearthreestatusCDF.pdf',yearthreestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearthreestatusCDFkm.pdf',yearthreestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by year (2004)
finaldata37<-finaldata[finaldata$Year=='2004',]; #subset for year 2004;
names(finaldata37)
head(finaldata37)
tail(finaldata37)

write.table(finaldata37,'finaldata37.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata37$siteID, Use=rep(TRUE,nrow(finaldata37)) )

subpop.ext <- data.frame(siteID=finaldata37$siteID,
			Region=rep('Year 2004',nrow(finaldata37)) )

design.ext <- data.frame(siteID=finaldata37$siteID,stratum=rep('1',nrow(finaldata37)),
			wgt=finaldata37$finalwgt, xcoord=finaldata37$xmarinus,
			ycoord=finaldata37$ymarinus)

data.cont.ext <- finaldata37[,c('siteID','VSCIAll')]

yearfourstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearfourstatus$CDF,'yearfourstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearfourstatus$Pct,'yearfourstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearfourstatus$Tot,'yearfourstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearfourstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearfourstatusCDF.pdf',yearfourstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearfourstatusCDFkm.pdf',yearfourstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by year (2005)
finaldata38<-finaldata[finaldata$Year=='2005',]; #subset for year 2005;
names(finaldata38)
head(finaldata38)
tail(finaldata38)

write.table(finaldata38,'finaldata38.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata38$siteID, Use=rep(TRUE,nrow(finaldata38)) )

subpop.ext <- data.frame(siteID=finaldata38$siteID,
			Region=rep('Year 2005',nrow(finaldata38)) )

design.ext <- data.frame(siteID=finaldata38$siteID,stratum=rep('1',nrow(finaldata38)),
			wgt=finaldata38$finalwgt, xcoord=finaldata38$xmarinus,
			ycoord=finaldata38$ymarinus)

data.cont.ext <- finaldata38[,c('siteID','VSCIAll')]

yearfivestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearfivestatus$CDF,'yearfivestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearfivestatus$Pct,'yearfivestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearfivestatus$Tot,'yearfivestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearfivestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearfivestatusCDF.pdf',yearfivestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearfivestatusCDFkm.pdf',yearfivestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by year (2006)
finaldata39<-finaldata[finaldata$Year=='2006',]; #subset for year 2006;
names(finaldata39)
head(finaldata39)
tail(finaldata39)

write.table(finaldata39,'finaldata39.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata39$siteID, Use=rep(TRUE,nrow(finaldata39)) )

subpop.ext <- data.frame(siteID=finaldata39$siteID,
			Region=rep('Year 2006',nrow(finaldata39)) )

design.ext <- data.frame(siteID=finaldata39$siteID,stratum=rep('1',nrow(finaldata39)),
			wgt=finaldata39$finalwgt, xcoord=finaldata39$xmarinus,
			ycoord=finaldata39$ymarinus)

data.cont.ext <- finaldata39[,c('siteID','VSCIAll')]

yearsixstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearsixstatus$CDF,'yearsixstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearsixstatus$Pct,'yearsixstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearsixstatus$Tot,'yearsixstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearsixstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearsixstatusCDF.pdf',yearsixstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearsixstatusCDFkm.pdf',yearsixstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           



# subset by year (2007)
finaldata40<-finaldata[finaldata$Year=='2007',]; #subset for year 2007;
names(finaldata40)
head(finaldata40)
tail(finaldata40)

write.table(finaldata40,'finaldata40.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata40$siteID, Use=rep(TRUE,nrow(finaldata40)) )

subpop.ext <- data.frame(siteID=finaldata40$siteID,
			Region=rep('Year 2007',nrow(finaldata40)) )

design.ext <- data.frame(siteID=finaldata40$siteID,stratum=rep('1',nrow(finaldata40)),
			wgt=finaldata40$finalwgt, xcoord=finaldata40$xmarinus,
			ycoord=finaldata40$ymarinus)

data.cont.ext <- finaldata40[,c('siteID','VSCIAll')]

yearsevenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearsevenstatus$CDF,'yearsevenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearsevenstatus$Pct,'yearsevenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearsevenstatus$Tot,'yearsevenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearsevenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearsevenstatusCDF.pdf',yearsevenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearsevenstatusCDFkm.pdf',yearsevenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           
                                                                   


# subset by year (2008)
finaldata41<-finaldata[finaldata$Year=='2008',]; #subset for year 2008;
names(finaldata41)
head(finaldata41)
tail(finaldata41)

write.table(finaldata41,'finaldata41.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata41$siteID, Use=rep(TRUE,nrow(finaldata41)) )

subpop.ext <- data.frame(siteID=finaldata41$siteID,
			Region=rep('Year 2008',nrow(finaldata41)) )

design.ext <- data.frame(siteID=finaldata41$siteID,stratum=rep('1',nrow(finaldata41)),
			wgt=finaldata41$finalwgt, xcoord=finaldata41$xmarinus,
			ycoord=finaldata41$ymarinus)

data.cont.ext <- finaldata41[,c('siteID','VSCIAll')]

yeareightstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yeareightstatus$CDF,'yeareightstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yeareightstatus$Pct,'yeareightstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yeareightstatus$Tot,'yeareightstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yeareightstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yeareightstatusCDF.pdf',yeareightstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yeareightstatusCDFkm.pdf',yeareightstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           
                                                                   



# subset by year (2009)
finaldata42<-finaldata[finaldata$Year=='2009',]; #subset for year 2009;
names(finaldata42)
head(finaldata42)
tail(finaldata42)

write.table(finaldata42,'finaldata42.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata42$siteID, Use=rep(TRUE,nrow(finaldata42)) )

subpop.ext <- data.frame(siteID=finaldata42$siteID,
			Region=rep('Year 2009',nrow(finaldata42)) )

design.ext <- data.frame(siteID=finaldata42$siteID,stratum=rep('1',nrow(finaldata42)),
			wgt=finaldata42$finalwgt, xcoord=finaldata42$xmarinus,
			ycoord=finaldata42$ymarinus)

data.cont.ext <- finaldata42[,c('siteID','VSCIAll')]

yearninestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearninestatus$CDF,'yearninestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearninestatus$Pct,'yearninestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearninestatus$Tot,'yearninestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearninestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearninestatusCDF.pdf',yearninestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearninestatusCDFkm.pdf',yearninestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           
                                                                   


# subset by year (2010)
finaldata43<-finaldata[finaldata$Year=='2010',]; #subset for year 2010;
names(finaldata43)
head(finaldata43)
tail(finaldata43)

write.table(finaldata43,'finaldata43.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata43$siteID, Use=rep(TRUE,nrow(finaldata43)) )

subpop.ext <- data.frame(siteID=finaldata43$siteID,
			Region=rep('Year 2010',nrow(finaldata43)) )

design.ext <- data.frame(siteID=finaldata43$siteID,stratum=rep('1',nrow(finaldata43)),
			wgt=finaldata43$finalwgt, xcoord=finaldata43$xmarinus,
			ycoord=finaldata43$ymarinus)

data.cont.ext <- finaldata43[,c('siteID','VSCIAll')]

yeartenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yeartenstatus$CDF,'yeartenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yeartenstatus$Pct,'yeartenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yeartenstatus$Tot,'yeartenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yeartenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yeartenstatusCDF.pdf',yeartenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yeartenstatusCDFkm.pdf',yeartenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           
                                                         
# subset by year (2011)
finaldata44<-finaldata[finaldata$Year=='2011',]; #subset for year 2011;
names(finaldata44)
head(finaldata44)
tail(finaldata44)

write.table(finaldata44,'finaldata44.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata44$siteID, Use=rep(TRUE,nrow(finaldata44)) )

subpop.ext <- data.frame(siteID=finaldata44$siteID,
			Region=rep('Year 2011',nrow(finaldata44)) )

design.ext <- data.frame(siteID=finaldata44$siteID,stratum=rep('1',nrow(finaldata44)),
			wgt=finaldata44$finalwgt, xcoord=finaldata44$xmarinus,
			ycoord=finaldata44$ymarinus)

data.cont.ext <- finaldata44[,c('siteID','VSCIAll')]

yearelevenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearelevenstatus$CDF,'yearelevenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearelevenstatus$Pct,'yearelevenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearelevenstatus$Tot,'yearelevenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearelevenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearelevenstatusCDF.pdf',yearelevenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearelevenstatusCDFkm.pdf',yearelevenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                     
                                                         
# subset by year (2012)
finaldata45<-finaldata[finaldata$Year=='2012',]; #subset for year 2012;
names(finaldata45)
head(finaldata45)
tail(finaldata45)

write.table(finaldata45,'finaldata45.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata45$siteID, Use=rep(TRUE,nrow(finaldata45)) )

subpop.ext <- data.frame(siteID=finaldata45$siteID,
			Region=rep('Year 2012',nrow(finaldata45)) )

design.ext <- data.frame(siteID=finaldata45$siteID,stratum=rep('1',nrow(finaldata45)),
			wgt=finaldata45$finalwgt, xcoord=finaldata45$xmarinus,
			ycoord=finaldata45$ymarinus)

data.cont.ext <- finaldata45[,c('siteID','VSCIAll')]

yeartwelvestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yeartwelvestatus$CDF,'yeartwelvestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yeartwelvestatus$Pct,'yeartwelvestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yeartwelvestatus$Tot,'yeartwelvestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yeartwelvestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yeartwelvestatusCDF.pdf',yeartwelvestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yeartwelvestatusCDFkm.pdf',yeartwelvestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                      


# subset by year (2013)
# note final data is 56 - out of order!
finaldata56<-finaldata[finaldata$Year=='2013',]; #subset for year 2013;
names(finaldata56)
head(finaldata56)
tail(finaldata56)

write.table(finaldata56,'finaldata56.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata56$siteID, Use=rep(TRUE,nrow(finaldata56)) )

subpop.ext <- data.frame(siteID=finaldata56$siteID,
			Region=rep('Year 2013',nrow(finaldata56)) )

design.ext <- data.frame(siteID=finaldata56$siteID,stratum=rep('1',nrow(finaldata56)),
			wgt=finaldata56$finalwgt, xcoord=finaldata56$xmarinus,
			ycoord=finaldata56$ymarinus)

data.cont.ext <- finaldata56[,c('siteID','VSCIAll')]

yearthirteenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearthirteenstatus$CDF,'yearthirteenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearthirteenstatus$Pct,'yearthirteenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearthirteenstatus$Tot,'yearthirteenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearthirteenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearthirteenstatusCDF.pdf',yearthirteenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearthirteenstatusCDFkm.pdf',yearthirteenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')     
                 
# subset by year (2014)
# note final data is 57 - out of order!
finaldata57<-finaldata[finaldata$Year=='2014',]; #subset for year 2014;
names(finaldata57)
head(finaldata57)
tail(finaldata57)

write.table(finaldata57,'finaldata57.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata57$siteID, Use=rep(TRUE,nrow(finaldata57)) )

subpop.ext <- data.frame(siteID=finaldata57$siteID,
			Region=rep('Year 2014',nrow(finaldata57)) )

design.ext <- data.frame(siteID=finaldata57$siteID,stratum=rep('1',nrow(finaldata57)),
			wgt=finaldata57$finalwgt, xcoord=finaldata57$xmarinus,
			ycoord=finaldata57$ymarinus)

data.cont.ext <- finaldata57[,c('siteID','VSCIAll')]

yearfourteenstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(yearfourteenstatus$CDF,'yearfourteenstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(yearfourteenstatus$Pct,'yearfourteenstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(yearfourteenstatus$Tot,'yearfourteenstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(yearfourteenstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='yearfourteenstatusCDF.pdf',yearfourteenstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='yearfourteenstatusCDFkm.pdf',yearfourteenstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                            




# subset by panel (2001-2007)
finaldata46<-finaldata[finaldata$Panel=='Phase1',]; #subset for phase one 2001-2007;
names(finaldata46)
head(finaldata46)
tail(finaldata46)

write.table(finaldata46,'finaldata46.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata46$siteID, Use=rep(TRUE,nrow(finaldata46)) )

subpop.ext <- data.frame(siteID=finaldata46$siteID,
			Region=rep('Phase One 2001-2007',nrow(finaldata46)) )

design.ext <- data.frame(siteID=finaldata46$siteID,stratum=rep('1',nrow(finaldata46)),
			wgt=finaldata46$finalwgt, xcoord=finaldata46$xmarinus,
			ycoord=finaldata46$ymarinus)

data.cont.ext <- finaldata46[,c('siteID','VSCIAll')]

phaseonestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(phaseonestatus$CDF,'phaseonestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(phaseonestatus$Pct,'phaseonestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(phaseonestatus$Tot,'phaseonestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(phaseonestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='phaseonestatusCDF.pdf',phaseonestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='phaseonestatusCDFkm.pdf',phaseonestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           
                                                                   


# subset by panel (2008-2014)
finaldata47<-finaldata[finaldata$Panel=='Phase2',]; #subset for phase one 2008-2014;
names(finaldata47)
head(finaldata47)
tail(finaldata47)

write.table(finaldata47,'finaldata47.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata47$siteID, Use=rep(TRUE,nrow(finaldata47)) )

subpop.ext <- data.frame(siteID=finaldata47$siteID,
			Region=rep('Phase Two 2008-2014',nrow(finaldata47)) )

design.ext <- data.frame(siteID=finaldata47$siteID,stratum=rep('1',nrow(finaldata47)),
			wgt=finaldata47$finalwgt, xcoord=finaldata47$xmarinus,
			ycoord=finaldata47$ymarinus)

data.cont.ext <- finaldata47[,c('siteID','VSCIAll')]

phasetwostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(phasetwostatus$CDF,'phasetwostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(phasetwostatus$Pct,'phasetwostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(phasetwostatus$Tot,'phasetwostatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(phasetwostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='phasetwostatusCDF.pdf',phasetwostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='phasetwostatusCDFkm.pdf',phasetwostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                           
                                                                 
# subset by bay (All 2001-2014)
finaldata48<-finaldata[finaldata$BayShed=='Bay',]; #subset for bay 2001-2014;
names(finaldata48)
head(finaldata48)
tail(finaldata48)

write.table(finaldata48,'finaldata48.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata48$siteID, Use=rep(TRUE,nrow(finaldata48)) )

subpop.ext <- data.frame(siteID=finaldata48$siteID,
			Region=rep('Bay Watersheds 2001-2014',nrow(finaldata48)) )

design.ext <- data.frame(siteID=finaldata48$siteID,stratum=rep('1',nrow(finaldata48)),
			wgt=finaldata48$finalwgt, xcoord=finaldata48$xmarinus,
			ycoord=finaldata48$ymarinus)

data.cont.ext <- finaldata48[,c('siteID','VSCIAll')]

baystatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(baystatus$CDF,'baystatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(baystatus$Pct,'baystatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(baystatus$Tot,'baystatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(baystatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='baystatusCDF.pdf',baystatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='baystatusCDFkm.pdf',baystatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                                                      
                                                                              
# subset by non-bay (All 2001-2014)
finaldata49<-finaldata[finaldata$BayShed=='NonBay',]; #subset for non-bay 2001-2014;
names(finaldata49)
head(finaldata49)
tail(finaldata49)

write.table(finaldata49,'finaldata49.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata49$siteID, Use=rep(TRUE,nrow(finaldata49)) )

subpop.ext <- data.frame(siteID=finaldata49$siteID,
			Region=rep('Non-Bay Watersheds 2001-2014',nrow(finaldata49)) )

design.ext <- data.frame(siteID=finaldata49$siteID,stratum=rep('1',nrow(finaldata49)),
			wgt=finaldata49$finalwgt, xcoord=finaldata49$xmarinus,
			ycoord=finaldata49$ymarinus)

data.cont.ext <- finaldata49[,c('siteID','VSCIAll')]

nonbaystatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(nonbaystatus$CDF,'nonbaystatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(nonbaystatus$Pct,'nonbaystatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(nonbaystatus$Tot,'nonbaystatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(nonbaystatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='nonbaystatusCDF.pdf',nonbaystatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='nonbaystatusCDFkm.pdf',nonbaystatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')   
                 
                 
                 
# subset by bay phase 1 (All bay 2001-2007)
finaldata50<-finaldata[finaldata$BayPanel=='BayPhase1',]; #subset for bay phase 1 - 2001-2007;
names(finaldata50)
head(finaldata50)
tail(finaldata50)

write.table(finaldata50,'finaldata50.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata50$siteID, Use=rep(TRUE,nrow(finaldata50)) )

subpop.ext <- data.frame(siteID=finaldata50$siteID,
			Region=rep('Bay Watersheds 2001-2007',nrow(finaldata50)) )

design.ext <- data.frame(siteID=finaldata50$siteID,stratum=rep('1',nrow(finaldata50)),
			wgt=finaldata50$finalwgt, xcoord=finaldata50$xmarinus,
			ycoord=finaldata50$ymarinus)

data.cont.ext <- finaldata50[,c('siteID','VSCIAll')]

bayphaseonestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(bayphaseonestatus$CDF,'bayphaseonestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(bayphaseonestatus$Pct,'bayphaseonestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(bayphaseonestatus$Tot,'bayphaseonestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(bayphaseonestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='bayphaseonestatusCDF.pdf',bayphaseonestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='bayphaseonestatusCDFkm.pdf',bayphaseonestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                                                      
                                                                              

# subset by bay phase 2 (All bay 2008-2014)
finaldata51<-finaldata[finaldata$BayPanel=='BayPhase2',]; #subset for bay phase 2 - 2008-2014;
names(finaldata51)
head(finaldata51)
tail(finaldata51)

write.table(finaldata51,'finaldata51.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata51$siteID, Use=rep(TRUE,nrow(finaldata51)) )

subpop.ext <- data.frame(siteID=finaldata51$siteID,
			Region=rep('Bay Watersheds 2008-2014',nrow(finaldata51)) )

design.ext <- data.frame(siteID=finaldata51$siteID,stratum=rep('1',nrow(finaldata51)),
			wgt=finaldata51$finalwgt, xcoord=finaldata51$xmarinus,
			ycoord=finaldata51$ymarinus)

data.cont.ext <- finaldata51[,c('siteID','VSCIAll')]

bayphasetwostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(bayphasetwostatus$CDF,'bayphasetwostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(bayphasetwostatus$Pct,'bayphasetwostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(bayphasetwostatus$Tot,'bayphasetwostatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(bayphasetwostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='bayphasetwostatusCDF.pdf',bayphasetwostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='bayphasetwostatusCDFkm.pdf',bayphasetwostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')   
                 
                 
                  
 
 # subset by non-bay phase 1 (All non-bay 2001-2007)
# note final data 58 out of order
finaldata58<-finaldata[finaldata$BayPanel=='NonBayPhase1',]; #subset for non-bay phase 1 - 2001-2007;
names(finaldata58)
head(finaldata58)
tail(finaldata58)

write.table(finaldata58,'finaldata58.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata58$siteID, Use=rep(TRUE,nrow(finaldata58)) )

subpop.ext <- data.frame(siteID=finaldata58$siteID,
			Region=rep('Non-Bay Watersheds 2001-2007',nrow(finaldata58)) )

design.ext <- data.frame(siteID=finaldata58$siteID,stratum=rep('1',nrow(finaldata58)),
			wgt=finaldata58$finalwgt, xcoord=finaldata58$xmarinus,
			ycoord=finaldata58$ymarinus)

data.cont.ext <- finaldata58[,c('siteID','VSCIAll')]

nonbayphaseonestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(nonbayphaseonestatus$CDF,'nonbayphaseonestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(nonbayphaseonestatus$Pct,'nonbayphaseonestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(nonbayphaseonestatus$Tot,'nonbayphaseonestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(nonbayphaseonestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='nonbayphaseonestatusCDF.pdf',nonbayphaseonestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='nonbayphaseonestatusCDFkm.pdf',nonbayphaseonestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  
                 
                 
# subset by non-bay phase 2 (All non-bay 2008-2014)
# note final data 59 out of order
finaldata59<-finaldata[finaldata$BayPanel=='NonBayPhase2',]; #subset for non-bay phase 2 - 2008-2014;
names(finaldata59)
head(finaldata59)
tail(finaldata59)

write.table(finaldata59,'finaldata59.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata59$siteID, Use=rep(TRUE,nrow(finaldata59)) )

subpop.ext <- data.frame(siteID=finaldata59$siteID,
			Region=rep('Non-Bay Watersheds 2008-2014',nrow(finaldata59)) )

design.ext <- data.frame(siteID=finaldata59$siteID,stratum=rep('1',nrow(finaldata59)),
			wgt=finaldata59$finalwgt, xcoord=finaldata59$xmarinus,
			ycoord=finaldata59$ymarinus)

data.cont.ext <- finaldata59[,c('siteID','VSCIAll')]

nonbayphasetwostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(nonbayphasetwostatus$CDF,'nonbayphasetwostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(nonbayphasetwostatus$Pct,'nonbayphasetwostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(nonbayphasetwostatus$Tot,'nonbayphasetwostatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(nonbayphasetwostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='nonbayphasetwostatusCDF.pdf',nonbayphasetwostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='nonbayphasetwostatusCDFkm.pdf',nonbayphasetwostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                   
                                                       
                                                                                                                    
                                                                              
# subset by bio phase 1 (All 2001-2003)
finaldata52<-finaldata[finaldata$BioPanel=='Phase1',]; #subset for bio phase 1 - 2001-2003;
names(finaldata52)
head(finaldata52)
tail(finaldata52)

write.table(finaldata52,'finaldata52.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata52$siteID, Use=rep(TRUE,nrow(finaldata52)) )

subpop.ext <- data.frame(siteID=finaldata52$siteID,
			Region=rep('VSCI Scores 2001-2003',nrow(finaldata52)) )

design.ext <- data.frame(siteID=finaldata52$siteID,stratum=rep('1',nrow(finaldata52)),
			wgt=finaldata52$finalwgt, xcoord=finaldata52$xmarinus,
			ycoord=finaldata52$ymarinus)

data.cont.ext <- finaldata52[,c('siteID','VSCIAll')]

biophaseonestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(biophaseonestatus$CDF,'biophaseonestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(biophaseonestatus$Pct,'biophaseonestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(biophaseonestatus$Tot,'biophaseonestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(biophaseonestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='biophaseonestatusCDF.pdf',biophaseonestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='biophaseonestatusCDFkm.pdf',biophaseonestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')              
                 
# subset by bio phase 2 (All 2004-2006)
finaldata53<-finaldata[finaldata$BioPanel=='Phase2',]; #subset for bio phase 2 - 2004-2006;
names(finaldata53)
head(finaldata53)
tail(finaldata53)

write.table(finaldata53,'finaldata53.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata53$siteID, Use=rep(TRUE,nrow(finaldata53)) )

subpop.ext <- data.frame(siteID=finaldata53$siteID,
			Region=rep('VSCI Scores 2004-2006',nrow(finaldata53)) )

design.ext <- data.frame(siteID=finaldata53$siteID,stratum=rep('1',nrow(finaldata53)),
			wgt=finaldata53$finalwgt, xcoord=finaldata53$xmarinus,
			ycoord=finaldata53$ymarinus)

data.cont.ext <- finaldata53[,c('siteID','VSCIAll')]

biophasetwostatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(biophasetwostatus$CDF,'biophasetwostatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(biophasetwostatus$Pct,'biophasetwostatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(biophasetwostatus$Tot,'biophasetwostatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(biophasetwostatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='biophasetwostatusCDF.pdf',biophasetwostatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='biophasetwostatusCDFkm.pdf',biophasetwostatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                                                                       


# subset by bio phase 3 (All 2007-2010)
finaldata54<-finaldata[finaldata$BioPanel=='Phase3',]; #subset for bio phase 3 - 2007-2010;
names(finaldata54)
head(finaldata54)
tail(finaldata54)

write.table(finaldata54,'finaldata54.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata54$siteID, Use=rep(TRUE,nrow(finaldata54)) )

subpop.ext <- data.frame(siteID=finaldata54$siteID,
			Region=rep('VSCI Scores 2007-2010',nrow(finaldata54)) )

design.ext <- data.frame(siteID=finaldata54$siteID,stratum=rep('1',nrow(finaldata54)),
			wgt=finaldata54$finalwgt, xcoord=finaldata54$xmarinus,
			ycoord=finaldata54$ymarinus)

data.cont.ext <- finaldata54[,c('siteID','VSCIAll')]

biophasethreestatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(biophasethreestatus$CDF,'biophasethreestatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(biophasethreestatus$Pct,'biophasethreestatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(biophasethreestatus$Tot,'biophasethreestatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(biophasethreestatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='biophasethreestatusCDF.pdf',biophasethreestatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='biophasethreestatusCDFkm.pdf',biophasethreestatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')    
                 
                                                                                                                      

# subset by bio phase 4 (All 2011-2014)
finaldata55<-finaldata[finaldata$BioPanel=='Phase4',]; #subset for bio phase 4 - 2011-2014;
names(finaldata55)
head(finaldata55)
tail(finaldata55)

write.table(finaldata55,'finaldata55.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata55$siteID, Use=rep(TRUE,nrow(finaldata55)) )

subpop.ext <- data.frame(siteID=finaldata55$siteID,
			Region=rep('VSCI Scores 2011-2014',nrow(finaldata55)) )

design.ext <- data.frame(siteID=finaldata55$siteID,stratum=rep('1',nrow(finaldata55)),
			wgt=finaldata55$finalwgt, xcoord=finaldata55$xmarinus,
			ycoord=finaldata55$ymarinus)

data.cont.ext <- finaldata55[,c('siteID','VSCIAll')]

biophasefourstatus <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(biophasefourstatus$CDF,'biophasefourstatus.CDF.tab',sep='\t',row.names=FALSE)
write.table(biophasefourstatus$Pct,'biophasefourstatus.Pct.tab',sep='\t',row.names=FALSE)
write.table(biophasefourstatus$Tot,'biophasefourstatus.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(biophasefourstatus$CDF$Indicator))
cont.cdfplot.fcn(pdffile='biophasefourstatusCDF.pdf',biophasefourstatus$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='biophasefourstatusCDFkm.pdf',biophasefourstatus$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                                                         
                                                      
                                                      
                                                      
#######                                                      
# subset by IR2008 (All 2001-2006 VSCI Data)
# note final data 60 out of order
finaldata60<-finaldata[finaldata$IR2008=='2008',]; #subset for IR2008 - 2001-2006;
names(finaldata60)
head(finaldata60)
tail(finaldata60)

write.table(finaldata60,'finaldata60.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata60$siteID, Use=rep(TRUE,nrow(finaldata60)) )

subpop.ext <- data.frame(siteID=finaldata60$siteID,
			Region=rep('IR2008',nrow(finaldata60)) )

design.ext <- data.frame(siteID=finaldata60$siteID,stratum=rep('1',nrow(finaldata60)),
			wgt=finaldata60$finalwgt, xcoord=finaldata60$xmarinus,
			ycoord=finaldata60$ymarinus)

data.cont.ext <- finaldata60[,c('siteID','VSCIAll')]

IR2008 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(IR2008$CDF,'IR2008.CDF.tab',sep='\t',row.names=FALSE)
write.table(IR2008$Pct,'IR2008.Pct.tab',sep='\t',row.names=FALSE)
write.table(IR2008$Tot,'IR2008.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(IR2008$CDF$Indicator))
cont.cdfplot.fcn(pdffile='IR2008CDF.pdf',IR2008$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='IR2008CDFkm.pdf',IR2008$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                   
                                                     

# subset by IR2010 (All 2003-2008 VSCI Data)
finaldata61<-finaldata[finaldata$IR2010=='2010',]; #subset for IR2010 - 2003-2008;
names(finaldata61)
head(finaldata61)
tail(finaldata61)

write.table(finaldata61,'finaldata61.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata61$siteID, Use=rep(TRUE,nrow(finaldata61)) )

subpop.ext <- data.frame(siteID=finaldata61$siteID,
			Region=rep('IR2010',nrow(finaldata61)) )

design.ext <- data.frame(siteID=finaldata61$siteID,stratum=rep('1',nrow(finaldata61)),
			wgt=finaldata61$finalwgt, xcoord=finaldata61$xmarinus,
			ycoord=finaldata61$ymarinus)

data.cont.ext <- finaldata61[,c('siteID','VSCIAll')]

IR2010 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(IR2010$CDF,'IR2010.CDF.tab',sep='\t',row.names=FALSE)
write.table(IR2010$Pct,'IR2010.Pct.tab',sep='\t',row.names=FALSE)
write.table(IR2010$Tot,'IR2010.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(IR2010$CDF$Indicator))
cont.cdfplot.fcn(pdffile='IR2010CDF.pdf',IR2010$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='IR2010CDFkm.pdf',IR2010$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')  


# subset by IR2012 (All 2005-2010 VSCI Data)
finaldata62<-finaldata[finaldata$IR2012=='2012',]; #subset for IR2012 - 2005-2010;
names(finaldata62)
head(finaldata62)
tail(finaldata62)

write.table(finaldata62,'finaldata62.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata62$siteID, Use=rep(TRUE,nrow(finaldata62)) )

subpop.ext <- data.frame(siteID=finaldata62$siteID,
			Region=rep('IR2012',nrow(finaldata62)) )

design.ext <- data.frame(siteID=finaldata62$siteID,stratum=rep('1',nrow(finaldata62)),
			wgt=finaldata62$finalwgt, xcoord=finaldata62$xmarinus,
			ycoord=finaldata62$ymarinus)

data.cont.ext <- finaldata62[,c('siteID','VSCIAll')]

IR2012 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(IR2012$CDF,'IR2012.CDF.tab',sep='\t',row.names=FALSE)
write.table(IR2012$Pct,'IR2012.Pct.tab',sep='\t',row.names=FALSE)
write.table(IR2012$Tot,'IR2012.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(IR2012$CDF$Indicator))
cont.cdfplot.fcn(pdffile='IR2012CDF.pdf',IR2012$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='IR2012CDFkm.pdf',IR2012$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')         
                 


# subset by IR2014 (All 2007-2012 VSCI Data)
finaldata63<-finaldata[finaldata$IR2014=='2014',]; #subset for IR2014 - 2007-2012;
names(finaldata63)
head(finaldata63)
tail(finaldata63)

write.table(finaldata63,'finaldata63.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata63$siteID, Use=rep(TRUE,nrow(finaldata63)) )

subpop.ext <- data.frame(siteID=finaldata63$siteID,
			Region=rep('IR2014',nrow(finaldata63)) )

design.ext <- data.frame(siteID=finaldata63$siteID,stratum=rep('1',nrow(finaldata63)),
			wgt=finaldata63$finalwgt, xcoord=finaldata63$xmarinus,
			ycoord=finaldata63$ymarinus)

data.cont.ext <- finaldata63[,c('siteID','VSCIAll')]

IR2014 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(IR2014$CDF,'IR2014.CDF.tab',sep='\t',row.names=FALSE)
write.table(IR2014$Pct,'IR2014.Pct.tab',sep='\t',row.names=FALSE)
write.table(IR2014$Tot,'IR2014.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(IR2014$CDF$Indicator))
cont.cdfplot.fcn(pdffile='IR2014CDF.pdf',IR2014$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='IR2014CDFkm.pdf',IR2014$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')     
                 

# subset by IR2016 (All 2009-2014 VSCI Data)
finaldata64<-finaldata[finaldata$IR2016=='2016',]; #subset for IR2016 - 2009-2014;
names(finaldata64)
head(finaldata64)
tail(finaldata64)

write.table(finaldata64,'finaldata64.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata64$siteID, Use=rep(TRUE,nrow(finaldata64)) )

subpop.ext <- data.frame(siteID=finaldata64$siteID,
			Region=rep('IR2016',nrow(finaldata64)) )

design.ext <- data.frame(siteID=finaldata64$siteID,stratum=rep('1',nrow(finaldata64)),
			wgt=finaldata64$finalwgt, xcoord=finaldata64$xmarinus,
			ycoord=finaldata64$ymarinus)

data.cont.ext <- finaldata64[,c('siteID','VSCIAll')]

IR2016 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(IR2016$CDF,'IR2016.CDF.tab',sep='\t',row.names=FALSE)
write.table(IR2016$Pct,'IR2016.Pct.tab',sep='\t',row.names=FALSE)
write.table(IR2016$Tot,'IR2016.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(IR2016$CDF$Indicator))
cont.cdfplot.fcn(pdffile='IR2016CDF.pdf',IR2016$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='IR2016CDFkm.pdf',IR2016$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')       
                 
                 
# subset by stream size (Using order and sq miles)
finaldata65<-finaldata[finaldata$StreamSizeCat=='Small',]; #subset for small stream (orders 1,2,3 & <10 sq miles;
names(finaldata65)
head(finaldata65)
tail(finaldata65)

write.table(finaldata65,'finaldata65.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata65$siteID, Use=rep(TRUE,nrow(finaldata65)) )

subpop.ext <- data.frame(siteID=finaldata65$siteID,
			Region=rep('Small',nrow(finaldata65)) )

design.ext <- data.frame(siteID=finaldata65$siteID,stratum=rep('1',nrow(finaldata65)),
			wgt=finaldata65$finalwgt, xcoord=finaldata65$xmarinus,
			ycoord=finaldata65$ymarinus)

data.cont.ext <- finaldata65[,c('siteID','VSCIAll')]

Small <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Small$CDF,'Small.CDF.tab',sep='\t',row.names=FALSE)
write.table(Small$Pct,'Small.Pct.tab',sep='\t',row.names=FALSE)
write.table(Small$Tot,'Small.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Small$CDF$Indicator))
cont.cdfplot.fcn(pdffile='SmallCDF.pdf',Small$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='SmallCDFkm.pdf',Small$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')      
                 

# subset by stream size (Using order and sq miles)
finaldata66<-finaldata[finaldata$StreamSizeCat=='Medium',]; #subset for medium stream (orders 1,2,3,4 & >10 but <50 sq miles;
names(finaldata66)
head(finaldata66)
tail(finaldata66)

write.table(finaldata66,'finaldata66.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata66$siteID, Use=rep(TRUE,nrow(finaldata66)) )

subpop.ext <- data.frame(siteID=finaldata66$siteID,
			Region=rep('Medium',nrow(finaldata66)) )

design.ext <- data.frame(siteID=finaldata66$siteID,stratum=rep('1',nrow(finaldata66)),
			wgt=finaldata66$finalwgt, xcoord=finaldata66$xmarinus,
			ycoord=finaldata66$ymarinus)

data.cont.ext <- finaldata66[,c('siteID','VSCIAll')]

Medium <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Medium$CDF,'Medium.CDF.tab',sep='\t',row.names=FALSE)
write.table(Medium$Pct,'Medium.Pct.tab',sep='\t',row.names=FALSE)
write.table(Medium$Tot,'Medium.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Medium$CDF$Indicator))
cont.cdfplot.fcn(pdffile='MediumCDF.pdf',Medium$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='MediumCDFkm.pdf',Medium$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                     
                 

# subset by stream size (Using order and sq miles)
finaldata67<-finaldata[finaldata$StreamSizeCat=='Large',]; #subset for medium stream (orders 3,4,5,6 & >50 sq miles;
names(finaldata67)
head(finaldata67)
tail(finaldata67)

write.table(finaldata67,'finaldata67.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata67$siteID, Use=rep(TRUE,nrow(finaldata67)) )

subpop.ext <- data.frame(siteID=finaldata67$siteID,
			Region=rep('Large',nrow(finaldata67)) )

design.ext <- data.frame(siteID=finaldata67$siteID,stratum=rep('1',nrow(finaldata67)),
			wgt=finaldata67$finalwgt, xcoord=finaldata67$xmarinus,
			ycoord=finaldata67$ymarinus)

data.cont.ext <- finaldata67[,c('siteID','VSCIAll')]

Large <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Large$CDF,'Large.CDF.tab',sep='\t',row.names=FALSE)
write.table(Large$Pct,'Large.Pct.tab',sep='\t',row.names=FALSE)
write.table(Large$Tot,'Large.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Large$CDF$Indicator))
cont.cdfplot.fcn(pdffile='LargeCDF.pdf',Large$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='LargeCDFkm.pdf',Large$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)') 
                 
                 
# subset by stream size and years (Using order and sq miles)
finaldata68<-finaldata[finaldata$StreamSizeCatPhase=='Phase1Small',]; #subset for small stream (orders 1,2,3 & <10 sq miles and phase 1 -2001-2007;
names(finaldata68)
head(finaldata68)
tail(finaldata68)

write.table(finaldata68,'finaldata68.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata68$siteID, Use=rep(TRUE,nrow(finaldata68)) )

subpop.ext <- data.frame(siteID=finaldata68$siteID,
			Region=rep('Phase1Small',nrow(finaldata68)) )

design.ext <- data.frame(siteID=finaldata68$siteID,stratum=rep('1',nrow(finaldata68)),
			wgt=finaldata68$finalwgt, xcoord=finaldata68$xmarinus,
			ycoord=finaldata68$ymarinus)

data.cont.ext <- finaldata68[,c('siteID','VSCIAll')]

Small1 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Small1$CDF,'Small1.CDF.tab',sep='\t',row.names=FALSE)
write.table(Small1$Pct,'Small1.Pct.tab',sep='\t',row.names=FALSE)
write.table(Small1$Tot,'Small1.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Small1$CDF$Indicator))
cont.cdfplot.fcn(pdffile='Small1CDF.pdf',Small1$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='Small1CDFkm.pdf',Small1$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')   
                 
# subset by stream size and years (Using order and sq miles)
finaldata69<-finaldata[finaldata$StreamSizeCatPhase=='Phase2Small',]; #subset for small stream (orders 1,2,3 & <10 sq miles and phase 2 -2008-2014;
names(finaldata69)
head(finaldata69)
tail(finaldata69)

write.table(finaldata69,'finaldata69.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata69$siteID, Use=rep(TRUE,nrow(finaldata69)) )

subpop.ext <- data.frame(siteID=finaldata69$siteID,
			Region=rep('Phase2Small',nrow(finaldata69)) )

design.ext <- data.frame(siteID=finaldata69$siteID,stratum=rep('1',nrow(finaldata69)),
			wgt=finaldata69$finalwgt, xcoord=finaldata69$xmarinus,
			ycoord=finaldata69$ymarinus)

data.cont.ext <- finaldata69[,c('siteID','VSCIAll')]

Small2 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Small2$CDF,'Small2.CDF.tab',sep='\t',row.names=FALSE)
write.table(Small2$Pct,'Small2.Pct.tab',sep='\t',row.names=FALSE)
write.table(Small2$Tot,'Small2.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Small2$CDF$Indicator))
cont.cdfplot.fcn(pdffile='Small2CDF.pdf',Small2$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='Small2CDFkm.pdf',Small2$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                 
                 
                 
# subset by stream size and years (Using order and sq miles)
finaldata70<-finaldata[finaldata$StreamSizeCatPhase=='Phase1Medium',]; #subset for medium stream and phase 1 -2001-2007;
names(finaldata70)
head(finaldata70)
tail(finaldata70)

write.table(finaldata70,'finaldata70.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata70$siteID, Use=rep(TRUE,nrow(finaldata70)) )

subpop.ext <- data.frame(siteID=finaldata70$siteID,
			Region=rep('Phase1Medium',nrow(finaldata70)) )

design.ext <- data.frame(siteID=finaldata70$siteID,stratum=rep('1',nrow(finaldata70)),
			wgt=finaldata70$finalwgt, xcoord=finaldata70$xmarinus,
			ycoord=finaldata70$ymarinus)

data.cont.ext <- finaldata70[,c('siteID','VSCIAll')]

Medium1 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Medium1$CDF,'Medium1.CDF.tab',sep='\t',row.names=FALSE)
write.table(Medium1$Pct,'Medium1.Pct.tab',sep='\t',row.names=FALSE)
write.table(Medium1$Tot,'Medium1.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Medium1$CDF$Indicator))
cont.cdfplot.fcn(pdffile='Medium1CDF.pdf',Medium1$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='Medium1CDFkm.pdf',Medium1$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')           
                 
# subset by stream size and years (Using order and sq miles)
finaldata71<-finaldata[finaldata$StreamSizeCatPhase=='Phase2Medium',]; #subset for medium stream and phase 2 -2008-2014;
names(finaldata71)
head(finaldata71)
tail(finaldata71)

write.table(finaldata71,'finaldata71.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata71$siteID, Use=rep(TRUE,nrow(finaldata71)) )

subpop.ext <- data.frame(siteID=finaldata71$siteID,
			Region=rep('Phase2Medium',nrow(finaldata71)) )

design.ext <- data.frame(siteID=finaldata71$siteID,stratum=rep('1',nrow(finaldata71)),
			wgt=finaldata71$finalwgt, xcoord=finaldata71$xmarinus,
			ycoord=finaldata71$ymarinus)

data.cont.ext <- finaldata71[,c('siteID','VSCIAll')]

Medium2 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Medium2$CDF,'Medium2.CDF.tab',sep='\t',row.names=FALSE)
write.table(Medium2$Pct,'Medium2.Pct.tab',sep='\t',row.names=FALSE)
write.table(Medium2$Tot,'Medium2.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Medium2$CDF$Indicator))
cont.cdfplot.fcn(pdffile='Medium2CDF.pdf',Medium2$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='Medium2CDFkm.pdf',Medium2$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                        
                 
                 
# subset by stream size and years (Using order and sq miles)
finaldata72<-finaldata[finaldata$StreamSizeCatPhase=='Phase1Large',]; #subset for large stream and phase 2 -2001-2007;
names(finaldata72)
head(finaldata72)
tail(finaldata72)

write.table(finaldata72,'finaldata72.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata72$siteID, Use=rep(TRUE,nrow(finaldata72)) )

subpop.ext <- data.frame(siteID=finaldata72$siteID,
			Region=rep('Phase1Large',nrow(finaldata72)) )

design.ext <- data.frame(siteID=finaldata72$siteID,stratum=rep('1',nrow(finaldata72)),
			wgt=finaldata72$finalwgt, xcoord=finaldata72$xmarinus,
			ycoord=finaldata72$ymarinus)

data.cont.ext <- finaldata72[,c('siteID','VSCIAll')]

Large1 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Large1$CDF,'Large1.CDF.tab',sep='\t',row.names=FALSE)
write.table(Large1$Pct,'Large1.Pct.tab',sep='\t',row.names=FALSE)
write.table(Large1$Tot,'Large1.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Large1$CDF$Indicator))
cont.cdfplot.fcn(pdffile='Large1CDF.pdf',Large1$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='Large1CDFkm.pdf',Large1$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')              
                 

# subset by stream size and years (Using order and sq miles)
finaldata73<-finaldata[finaldata$StreamSizeCatPhase=='Phase2Large',]; #subset for large stream and phase 2 -2008-2014;
names(finaldata73)
head(finaldata73)
tail(finaldata73)

write.table(finaldata73,'finaldata73.tab',sep='\t',row.names=FALSE)

## biology estimates
sites.ext <- data.frame(siteID=finaldata73$siteID, Use=rep(TRUE,nrow(finaldata73)) )

subpop.ext <- data.frame(siteID=finaldata73$siteID,
			Region=rep('Phase2Large',nrow(finaldata73)) )

design.ext <- data.frame(siteID=finaldata73$siteID,stratum=rep('1',nrow(finaldata73)),
			wgt=finaldata73$finalwgt, xcoord=finaldata73$xmarinus,
			ycoord=finaldata73$ymarinus)

data.cont.ext <- finaldata73[,c('siteID','VSCIAll')]

Large2 <- cont.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cont = data.cont.ext, 
		pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
		conf=95, vartype="Local")

write.table(Large2$CDF,'Large2.CDF.tab',sep='\t',row.names=FALSE)
write.table(Large2$Pct,'Large2.Pct.tab',sep='\t',row.names=FALSE)
write.table(Large2$Tot,'Large2.Tot.tab',sep='\t',row.names=FALSE)


# plot estimates
source('cdfplot.fcn.r')
source('cont.plot.fcn.r')

xlab <- as.character(unique(Large2$CDF$Indicator))
cont.cdfplot.fcn(pdffile='Large2CDF.pdf',Large2$CDF,cdftype='Percent',
                 xlab=xlab,ylab='Percent')
cont.cdfplot.fcn(pdffile='Large2CDFkm.pdf',Large2$CDF,cdftype='km',
                 xlab=xlab,ylab='Stream Length(km)')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   