# ------------------------
#Van Sickle Code modified March 3rd, 2017
#Relative Risk Estimate
#R3.3.3 64 Bit
#SPsurvey 3.3




### Install packages
#install.packages('plyr')
#install.packages('dplyr')
#install.packages('ggplot2')
#install.packages('spsurvey')

# Loading packages 
library(plyr)
library(dplyr)
library(spsurvey)
# ALWAYS LOAD plyr BEFORE dplyr

### Load data
# As .csv 
prob <- readxl::read_excel('originalData/ProbMonWadeable2001-2016.xlsx',sheet='ProbMonData2001-2016_EVJ') 

# Pipes, these are streamlining operations allowing you to flow from one command to the next
prob3 <- select(prob,StationID:MetalCCU) # lose dissolved metals and GIS data

####Need to remove the fair....need to run code at home to check results

prob3$VSCIstatus <- cut(prob3$VSCIVCPMI,c(0,50,60,100),labels=c('Poor','Fair','Good'))  
prob3$TotHabstatus <- cut(prob3$TotHab, c(0,120,150,200),labels=c('Poor','Fair','Good'))
prob3$TNstatus <- cut(prob3$TN, c(0,1,2,100),labels=c('Good','Fair','Poor'))
prob3$TPstatus <- cut(prob3$TP, c(0,0.02,0.05,100),labels=c('Good','Fair','Poor'))
prob3$TDSstatus <- cut(prob3$TDS, c(0,100,350,20000),labels=c('Good','Fair','Poor'))
prob3$MetalCCUstatus <- cut(prob3$MetalCCU, c(0,1,2,100),labels=c('Good','Fair','Poor'))
prob3$LRBSstatus <- cut(prob3$LRBS, c(-10,-1,-0.5,0.5,10),labels=c('Poor','Fair','Good','Fair2'))
names(prob3)
prob3$VSCIstatus[ prob3$VSCIstatus == "Fair" ] = NA
prob3$TotHabstatus[ prob3$TotHabstatus == "Fair" ] = NA
prob3$TNstatus[ prob3$TNstatus == "Fair" ] = NA
prob3$TPstatus[ prob3$TPstatus == "Fair" ] = NA
prob3$TDSstatus[ prob3$TDSstatus == "Fair" ] = NA
prob3$MetalCCUstatus[ prob3$MetalCCUstatus == "Fair" ] = NA
prob3$LRBSstatus[ prob3$LRBSstatus == "Fair" ] = NA
prob3$LRBSstatus[ prob3$LRBSstatus == "Fair2" ] = NA
###This code puts fair back into data frame
#prob3$VSCIstatus[is.na(prob3$VSCIstatus)] <- "Fair"

###JRH Attempt at RR in new SP survey


# Create a variable that contains only the name of the response variable.
resp.var<-"VSCIstatus";
resp.var

#########.
# For stressor variables, we will initially select a few stressor variables 
# All of these must be condition class variables.

stres.vars<- c("TotHabstatus","TNstatus","TPstatus","TDSstatus","MetalCCUstatus","LRBSstatus");
stres.vars

# Create a vector containing the names of all selected stressor variables
#stres.vars<-c("PTL_COND","NTL_COND","TURB_COND","ANC_COND","SALINITY_COND",
#              "RDIS_COND", "RVEG_COND", "LITCVR_COND","LITRIPCVR_COND");


# First, set up the 4 data frames, using columns in ProbMetrics

# 4.1) The sites data frame:
sites.va.rr<-data.frame(siteID=prob3$StationID_Trend, Use=rep(TRUE, nrow(prob3)));

# ii) The subpopulation data frame. RR and AR estimates have high uncertainties
# and require large sample sizes. Thus, we will estimate RR and AR only for our largest sample size, 
# which is the whole state basis.
subpop.va.rr <- data.frame(siteID=prob3$StationID_Trend, all.virginia=rep("All_of_VA", nrow(prob3)));


# iii) The design data frame is the same as for extent estimation. However, to be safe,  
#      let's rebuild this data frame, since we are now working from texas.dat.rr;

# add marinus projection coordinates
#tmp <- marinus(prob3$LatitudeDD,-prob3$LongitudeDD)
#tmp$xmarinus <- tmp[,'x']
#tmp$ymarinus <- tmp[,'y']
####Noticed USEPA using albers ver marinus projection now (using geodalbers function spsurvey)


albers.cord.rr<-geodalbers(lon=prob3$LongitudeDD,lat=prob3$LatitudeDD);
design.va.rr<-data.frame(siteID=prob3$StationID_Trend, 
                         xcoord=albers.cord.rr$xcoord,ycoord=albers.cord.rr$ycoord, 
                         wgt=prob3$filwgt_trend);


# iv) The data.cat data frame should contain siteID, plus all stressor
# and response variables;
data.cat.va.rr<-subset(prob3, select=c("StationID_Trend",resp.var, stres.vars),drop=T);
names(data.cat.va.rr)[1]<-"StationID_Trend";


#Van Sickle RR example
# First, set up the 4 data frames, using columns in texas.dat.rr.

# 4.1) The sites data frame:
#sites.tx.rr<-data.frame(siteID=texas.dat.rr$SITE_ID, Use=rep(TRUE, nrow(texas.dat.rr)));

# ii) The subpopulation data frame. RR and AR estimates have high uncertainties
# and require large sample sizes. Thus, we will estimate RR and AR only for our largest sample size, 
# which is the whole stateon a whole-state basis.
#subpop.tx.rr <- data.frame(siteID=texas.dat.rr$SITE_ID,
#                           all.texas=rep("All_of_Texas", nrow(texas.dat.rr)));

# iii) The design data frame is the same as for extent estimation. However, to be safe,  
#      let's rebuild this data frame, since we are now working from texas.dat.rr;

#albers.cord.rr<-geodalbers(lon=texas.dat.rr$LON_DD,lat=texas.dat.rr$LAT_DD);
#design.tx.rr<-data.frame(siteID=texas.dat.rr$SITE_ID, 
#                         xcoord=albers.cord.rr$xcoord,ycoord=albers.cord.rr$ycoord, 
#                         wgt=texas.dat.rr$WGT_NLA);

# iv) The data.cat data frame should contain siteID, plus all stressor
# and response variables;
#data.cat.tx.rr<-subset(texas.dat.rr, select=c("SITE_ID",resp.var, stres.vars),drop=T);
#names(data.cat.tx.rr)[1]<-"siteID";
#
# The 4 input data frames are complete.
##################.

#Finally, we implement relrisk.analysis().

relrisk.estimates<-relrisk.analysis(sites=sites.va.rr,subpop=subpop.va.rr,
                                    design=design.va.rr,data.rr=data.cat.va.rr,
                                    response.var=rep(resp.var,length(stres.vars)),
                                    stressor.var=stres.vars,
                                    response.levels=rep(list(c("Poor","Good")),length(stres.vars)), 
                                    stressor.levels=rep(list(c("Poor","Good")),length(stres.vars)));


write.csv(relrisk.estimates, file = "relriskIR2016.csv", row.names = FALSE)
write.csv(prob3, file = "prob3.csv", row.names = FALSE)
