# Analysis starting point

library(tidyverse)
library(spsurvey)

# Bring in sites sampled and one parameter, for this example we will use biology
masterdata <- read.delim('analysis/masterdata.tab') %>%
  mutate(siteID=Station)

# Bring in design status data, this includes all information about sites that weren't sampled, etc.
# the nrow of survey data should be more than the sampled data
statusdata <- read.delim('analysis/biology.tab',comment.char = "")

# rename columns
dsgnstatus <- select(statusdata,'sampleID','strahler.order','Longitude.DD','Latitude.DD','design.weight',
                     'weight.category', 'station', 'state','status', 'set','Basin','SubBasin','BayShed','BayPanel',
                     'EcoRegion','BioRegion','Year','Panel','BioPanel','Order','BasinSize','IR2008','IR2010',
                     'IR2012','IR2014','IR2016','StreamSizeCat','StreamSizeCatPhase') 
names(dsgnstatus) <- c('siteID','StrahlerOrder','LongDD','LatDD',
                       'InitialWeight','MDCaty','SiteNum','State','Status','Set','Basin','SubBasin',
                       'BayShed','BayPanel','EcoRegion','BioRegion','Year','Panel','BioPanel','Order',
                       'BasinSize','IR2008','IR2010','IR2012','IR2014','IR2016','StreamSizeCat','StreamSizeCatPhase')

# add marinus projection coordinates
dsgnstatus <- cbind(dsgnstatus,marinus(dsgnstatus$LatDD,-dsgnstatus$LongDD))%>%
  rename(xmarinus=x,ymarinus=y)

# Current example uses statewide Strahler Order 
# List stream kilometer by km

sframe <- c('1st'=51210, '2nd'=13680, '3rd'=7781.08, '4th'=4448.257, 
            '5th'=1731.302, '6th'=163.901, '7th'=14.7099 )



## Create weights for population estimation
# recode StrahlerOrder to match sframe 
dsgnstatus$StrahlerOrder <- as.factor(dsgnstatus$StrahlerOrder)
levels(dsgnstatus$StrahlerOrder) <- c('1st','2nd','3rd','4th','5th','6th','7th')
table(dsgnstatus$StrahlerOrder,dsgnstatus$MDCaty)
wgt <- sframe/table(dsgnstatus$StrahlerOrder)

dsgnstatus$finalweight <- adjwgt(rep(TRUE,nrow(dsgnstatus)), 
                              dsgnstatus$InitialWeight,
                              dsgnstatus$StrahlerOrder, sframe)

####### Stream Extent and Status estimation
# need to make sure no siteID names are duplicated, e.g. use the StationID_Trend
#sites.ext <- data.frame(siteID=dsgnstatus$siteID, Use=rep(TRUE,nrow(dsgnstatus)) )

#subpop.ext <- data.frame(siteID=dsgnstatus$siteID,Region=rep('Virginia',nrow(dsgnstatus)) )

#design.ext <- data.frame(siteID=dsgnstatus$siteID,stratum=rep(1,nrow(dsgnstatus)),
#                         wgt=dsgnstatus$finalweight, xcoord=dsgnstatus$xmarinus,
#                         ycoord=dsgnstatus$ymarinus)

#StatusTNT <- dsgnstatus$Status
#levels(StatusTNT)
#levels(StatusTNT) <- list(T=c('TS','PD','OT'), NT=c('NT') )

#data.cat.ext <- data.frame(dsgnstatus[,c('siteID','Status')],StatusTNT)

#popstatus.est <- cat.analysis(sites = sites.ext, subpop = subpop.ext, design = design.ext, data.cat = data.cat.ext, 
#                              conf=95, vartype="Local")


### Match dsgnstatus with masterdata
indx <- match(masterdata$Station,dsgnstatus$siteID)
# does site status agree with having data available?
table(dsgnstatus$Status[indx])
cbind(dsgnstatus$siteID[indx],dsgnstatus$Status[indx])

finalData <- left_join(masterdata,dsgnstatus,by='siteID') %>%
  select(siteID,everything(),-c(Station))




## Subpopulation estimates by parameter
subpopEstimate <- function(finalData,parameter, subpopulationCategory, subpopulation, altName, specialWeight){
  # Build in catch in case whole population needed
  if(is.na(subpopulationCategory) & subpopulation == "Virginia"){
    subpopData <- finalData
    sites.ext <- select(subpopData, siteID) %>% mutate(Use = TRUE)
    subpop.ext <- select(subpopData,siteID) %>% mutate(Region = 'Virginia')
  }else{
    subpopData <- finalData[ finalData[[subpopulationCategory]] == subpopulation, ] 
    if(nrow(subpopData) == nrow(finalData)){ # special catch for IR windows not filtering correctly
      subpopData <- finalData[ !is.na(finalData[[subpopulationCategory]] ), ] 
    }}
  
  # If no data in subpopulation, keep moving
  if(nrow(subpopData) !=0){
    sites.ext <- select(subpopData, siteID) %>% mutate(Use = TRUE)
    # Special Cases to match existing terminology for each Subpopulation
    if(is.na(altName)){
      subpop.ext <- select(subpopData,siteID) %>% mutate(Region = subpopulation)
    }else{
      subpop.ext <- select(subpopData,siteID) %>% mutate(Region = altName)
    }
    
   
    
  # Choose correct final weight and filter to stratum = 1
    if(specialWeight==FALSE){
      finalweight <- subpopData$finalweight_all
    }else{
      finalweight <- as.numeric(as.matrix(subpopData[,specialWeight]))
    }

  
  design.ext <- mutate(subpopData,siteID = siteID, stratum = "1", wgt = finalweight, xcoord = xmarinus, ycoord = ymarinus) %>%
    select(siteID, stratum, wgt, xcoord, ycoord)
  
  data.cont.ext <- select(finalData, siteID, parameter)
  
  subpopStatus <- cont.analysis(sites = data.frame(sites.ext), subpop = data.frame(subpop.ext), design = data.frame(design.ext), 
                                data.cont = data.frame(data.cont.ext),
                                pctval=c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 99), 
                                conf=95, vartype="Local")
  return(subpopStatus)}else{return(list())}
}



# Now build it all into a single function that returns a list of lists (one list object per estimate, 
# each estimate is a list of of 4 dataframes where the important ones ($CDF and $Pct) can easily be 
# reached with list calls noted at bottom of page)
listOfResults <- function(popstatus.est, finalData, parameterName){
  list(
    popstatus.est = popstatus.est,
    dataAnalyzed = finalData,
    # All Virginia
    estimateVirginia = subpopEstimate(finalData,parameterName,NA, 'Virginia',NA, specialWeight=FALSE),
    # By Basin
    estimateRoanoke = subpopEstimate(finalData,parameterName, 'Basin', 'Roanoke','Roanoke Basin', specialWeight=FALSE),
    estimateJames = subpopEstimate(finalData,parameterName, 'Basin', 'James','James Basin', specialWeight=FALSE) ,
    estimatePotomacShenandoah = subpopEstimate(finalData,parameterName, 'Basin', 'Potomac-Shenandoah',NA, specialWeight=FALSE) ,
    estimateRappahannockYork = subpopEstimate(finalData,parameterName, 'Basin', 'Rappahannock-York',NA, specialWeight=FALSE),
    estimateNew = subpopEstimate(finalData,parameterName, 'Basin', 'New',NA, specialWeight=FALSE) ,
    estimateChowan = subpopEstimate(finalData,parameterName, 'Basin', 'Chowan',NA, specialWeight=FALSE), 
    estimateTennessee = subpopEstimate(finalData,parameterName, 'Basin', 'Tennessee',NA, specialWeight=FALSE), 
    estimateHolston = subpopEstimate(finalData,parameterName, 'SubBasin', 'Holston',NA, specialWeight=FALSE) ,
    estimateBigSandy = subpopEstimate(finalData,parameterName, 'SubBasin', 'Big Sandy',NA, specialWeight=FALSE) ,
    estimateClinchPowell = subpopEstimate(finalData,parameterName, 'SubBasin', 'Clinch-Powell',NA, specialWeight=FALSE) ,
    estimatePotomac = subpopEstimate(finalData,parameterName, 'SubBasin', 'Potomac',NA, specialWeight=FALSE) ,
    estimateShenandoah = subpopEstimate(finalData,parameterName, 'SubBasin', 'Shenandoah',NA, specialWeight=FALSE) ,
    estimateRappahannock = subpopEstimate(finalData,parameterName, 'SubBasin', 'Rappahannock',NA, specialWeight=FALSE) ,
    estimateYork = subpopEstimate(finalData,parameterName, 'SubBasin', 'York',NA, specialWeight=FALSE) ,
    # By Ecoregion
    estimatePiedmont = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Piedmont',NA, specialWeight=FALSE) ,
    estimateNorthernPiedmont = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Northern Piedmont',NA, specialWeight=FALSE) ,
    estimateCARV = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Central Appalachian Ridges and Valleys',NA, specialWeight=FALSE) ,
    estimateSEplains = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Southeastern Plains',NA, specialWeight=FALSE) ,
    estimateBRM = subpopEstimate(finalData,parameterName, 'EcoRegion', 'Blue Ridge Mountains',NA, specialWeight=FALSE) ,
    estimateMountainBioregion = subpopEstimate(finalData,parameterName, 'BioRegion', 'Mountain','Mountain Bioregion', specialWeight=FALSE),
    estimatePiedmontBioregion = subpopEstimate(finalData,parameterName, 'BioRegion', 'Piedmont','Piedmont Bioregion', specialWeight=FALSE) ,
    estimateCoastBioregion = subpopEstimate(finalData,parameterName, 'BioRegion', 'Coast', 'Coast Bioregion', specialWeight=FALSE),
    # By Order
    estimateFirstOrder = subpopEstimate(finalData,parameterName, 'Order', '1','First Order', specialWeight=FALSE) ,
    estimateSecondOrder = subpopEstimate(finalData,parameterName, 'Order', '2','Second Order', specialWeight=FALSE) ,
    estimateThirdOrder = subpopEstimate(finalData,parameterName, 'Order', '3','Third Order', specialWeight=FALSE) ,
    estimateFourthOrder = subpopEstimate(finalData,parameterName, 'Order', '4','Fourth Order', specialWeight=FALSE) ,
    estimateFifthOrder = subpopEstimate(finalData,parameterName, 'Order', '5','Fifth Order', specialWeight=FALSE) ,
    # By Basin Size
    estimateBasin1 = subpopEstimate(finalData,parameterName, 'BasinSize', '1','<1 square mile', specialWeight=FALSE) ,
    estimateBasin2 = subpopEstimate(finalData,parameterName, 'BasinSize', '2','1 to 10 square mile', specialWeight=FALSE) ,
    estimateBasin3 = subpopEstimate(finalData,parameterName, 'BasinSize', '3','10 to 200 square mile', specialWeight=FALSE) ,
    estimateBasin4 = subpopEstimate(finalData,parameterName, 'BasinSize', '4','>200 square mile', specialWeight=FALSE) ,
    # By Year
    estimate2001 = subpopEstimate(finalData,parameterName, 'Year', '2001','Year 2001', specialWeight='finalweight_Year') ,
    estimate2002 = subpopEstimate(finalData,parameterName, 'Year', '2002','Year 2002', specialWeight='finalweight_Year'), 
    estimate2003 = subpopEstimate(finalData,parameterName, 'Year', '2003','Year 2003', specialWeight='finalweight_Year') ,
    estimate2004 = subpopEstimate(finalData,parameterName, 'Year', '2004','Year 2004', specialWeight='finalweight_Year') ,
    estimate2005 = subpopEstimate(finalData,parameterName, 'Year', '2005','Year 2005', specialWeight='finalweight_Year') ,
    estimate2006 = subpopEstimate(finalData,parameterName, 'Year', '2006','Year 2006', specialWeight='finalweight_Year') ,
    estimate2007 = subpopEstimate(finalData,parameterName, 'Year', '2007','Year 2007', specialWeight='finalweight_Year') ,
    estimate2008 = subpopEstimate(finalData,parameterName, 'Year', '2008','Year 2008', specialWeight='finalweight_Year') ,
    estimate2009 = subpopEstimate(finalData,parameterName, 'Year', '2009','Year 2009', specialWeight='finalweight_Year') ,
    estimate2010 = subpopEstimate(finalData,parameterName, 'Year', '2010','Year 2010', specialWeight='finalweight_Year') ,
    estimate2011 = subpopEstimate(finalData,parameterName, 'Year', '2011','Year 2011', specialWeight='finalweight_Year') ,
    estimate2012 = subpopEstimate(finalData,parameterName, 'Year', '2012','Year 2012', specialWeight='finalweight_Year') ,
    estimate2013 = subpopEstimate(finalData,parameterName, 'Year', '2013','Year 2013', specialWeight='finalweight_Year') ,
    estimate2014 = subpopEstimate(finalData,parameterName, 'Year', '2014','Year 2014', specialWeight='finalweight_Year') ,
    estimate2015 = subpopEstimate(finalData,parameterName, 'Year', '2015','Year 2015', specialWeight='finalweight_Year') ,
    estimate2016 = subpopEstimate(finalData,parameterName, 'Year', '2016','Year 2016', specialWeight='finalweight_Year') ,
    # By Phase
    estimatePhase1 = subpopEstimate(finalData,parameterName, 'Panel', 'Phase1','Phase One 2001-2009', specialWeight='finalweight_Panel1') ,
    estimatePhase2 = subpopEstimate(finalData,parameterName, 'Panel', 'Phase2','Phase Two 2009-2016', specialWeight='finalweight_Panel2') ,
    # Bay/ NonBay
    estimateBay = subpopEstimate(finalData,parameterName, 'BayShed', 'Bay','Bay Watersheds 2001-2016', specialWeight=FALSE),
    estimateNonBay = subpopEstimate(finalData,parameterName, 'BayShed', 'NonBay','Non-Bay Watersheds 2001-2016', specialWeight=FALSE),
    estimateBayPhase1 = subpopEstimate(finalData,parameterName, 'BayPanel', 'BayPhase1','Bay Watersheds 2001-2009', specialWeight='finalweight_Panel1'),
    estimateBayPhase2 = subpopEstimate(finalData,parameterName, 'BayPanel', 'BayPhase2','Bay Watersheds 2009-2016', specialWeight='finalweight_Panel2'),
    estimateNonBayPhase1 = subpopEstimate(finalData,parameterName, 'BayPanel', 'NonBayPhase1','Non-Bay Watersheds 2001-2009', specialWeight='finalweight_Panel1'),
    estimateNonBayPhase2 = subpopEstimate(finalData,parameterName, 'BayPanel', 'NonBayPhase2','Non-Bay Watersheds 2009-2016', specialWeight='finalweight_Panel2'),
    # Sampling Phase
    estimateBioPanelPhase1 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase1','VSCI Scores 2001-2004', specialWeight='finalweight_BioPanel1'),
    estimateBioPanelPhase2 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase2','VSCI Scores 2005-2008', specialWeight='finalweight_BioPanel2'),
    estimateBioPanelPhase3 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase3','VSCI Scores 2009-2012', specialWeight='finalweight_BioPanel3'),
    estimateBioPanelPhase4 = subpopEstimate(finalData,parameterName, 'BioPanel', 'Phase4','VSCI Scores 2013-2016', specialWeight='finalweight_BioPanel4'),
    # IR window
    estimateIR2008 = subpopEstimate(finalData,parameterName, 'IR2008', '2008','IR2008', specialWeight='finalweight_IR2008'),
    estimateIR2010 = subpopEstimate(finalData,parameterName, 'IR2010', '2010','IR2010', specialWeight='finalweight_IR2010'),
    estimateIR2012= subpopEstimate(finalData,parameterName, 'IR2012', '2012','IR2012', specialWeight='finalweight_IR2012'),
    estimateIR2014 = subpopEstimate(finalData,parameterName, 'IR2014', '2014','IR2014', specialWeight='finalweight_IR2014'),
    estimateIR2016 = subpopEstimate(finalData,parameterName, 'IR2016', '2016','IR2016', specialWeight='finalweight_IR2016'),
    estimateIR2018 = subpopEstimate(finalData,parameterName, 'IR2018', '2018','IR2018', specialWeight='finalweight_IR2018'),
    
    # Stream size
    estimateSmall = subpopEstimate(finalData,parameterName, 'StreamSizeCat','Small',NA, specialWeight=FALSE),
    estimateMedium = subpopEstimate(finalData,parameterName, 'StreamSizeCat','Medium',NA, specialWeight=FALSE),
    estimateLarge = subpopEstimate(finalData,parameterName, 'StreamSizeCat','Large',NA, specialWeight=FALSE),
    # Stream size and sample phase
    estimatePhase1Small = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase1Small',NA, specialWeight='finalweight_Panel1'),
    estimatePhase2Small = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase2Small',NA, specialWeight='finalweight_Panel2'),
    estimatePhase1Medium = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase1Medium',NA, specialWeight='finalweight_Panel1'),
    estimatePhase2Medium = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase2Medium',NA, specialWeight='finalweight_Panel2'),
    estimatePhase1Large = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase1Large',NA, specialWeight='finalweight_Panel1'),
    estimatePhase2Large = subpopEstimate(finalData,parameterName, 'StreamSizeCatPhase','Phase2Large',NA, specialWeight='finalweight_Panel2')
  )
}

VSCI <- listOfSubpopulations(finalData, 'VSCIAll')
VSCI_CDFdf <- suppressWarnings(VSCI %>% map_df(1))
VSCI_PCTdf <-  suppressWarnings(VSCI %>% map_df(2))


# How to view individual parts of lists
names(VSCI)
View(VSCI[['estimateVirginia']])
View(VSCI[['estimateVirginia']]$CDF)
View(VSCI[['estimateVirginia']]$Pct)

# lets imagine for next parameter
allOtherParameters <- read_csv('C:/HardDriveBackup/R/IR_2016/data/ProbMetrics_2001-2014_Final_Web_March_7_2017.csv') %>%
  select(StationID_Trend, Year, DO:TDS,MetalCCU) %>%
  rename(Station=StationID_Trend)
masterdata <- left_join(masterdata,allOtherParameters,by='Station')


finalDataDO <- select(masterdata,Station, DO, siteID) %>%
  left_join(dsgnstatus,by='siteID') %>%
  select(siteID,everything(),-c(Station))

DO <- listOfSubpopulations(finalDataDO, 'DO')
DO_CDFdf <- suppressWarnings(DO %>% map_df(1))
DO_PCTdf <-  suppressWarnings(DO %>% map_df(2))

# Cool, now how to run through all desired parameters
parameters <- c('VSCIAll','DO','pH')
allStuff <- list()

for(i in parameters){
  print(i)
  finalData <- select(masterdata,Station, i, siteID) %>%
    left_join(dsgnstatus,by='siteID') %>%
    select(siteID,everything(),-c(Station))
  parameter <- listOfSubpopulations(finalData, i)
  CDFdf <- suppressWarnings(parameter %>% map_df(1))
  PCTdf <-  suppressWarnings(parameter %>% map_df(2))
  
  allStuff[[i]] <- list(estimates = parameter, CDF = CDFdf, PCT = PCTdf)

}

# everything stored in allStuff neatly
# how to get each list element stacked into single dataframe
allCDF <- suppressWarnings(allStuff %>% map_df(2))
allPCT <- suppressWarnings(allStuff %>% map_df(3))



### Code ends###################################################################################################

# Building block steps


##### Virginia can be used in this same function
#finalData <- finaldata
#parameter <- 'VSCIAll'
#subpopulationCategory <- "Basin" 
#subpopulation <- "Roanoke"
#altName <- 'Roanoke Basin'

#subpopulationCategory <- NA
#subpopulation <- 'Virginia'
#altName <- NA


# Statewide
estimateVirginia <- subpopEstimate(finalData,'VSCIAll',NA, 'Virginia',NA)
# By Basin
estimateRoanoke <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'Roanoke','Roanoke Basin')
estimateJames <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'James','James Basin') 
estimatePotomacShenandoah <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'Potomac-Shenandoah',NA) 
estimateRappahannockYork <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'Rappahannock-York',NA) 
estimateNew <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'New',NA) 
estimateChowan <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'Chowan',NA) 
estimateTennessee <- subpopEstimate(finalData,'VSCIAll', 'Basin', 'Tennessee',NA) 
estimateHolston <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'Holston',NA) 
estimateBigSandy <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'Big Sandy',NA) 
estimateClinchPowell <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'Clinch-Powell',NA) 
estimatePotomac <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'Potomac',NA) 
estimateShenandoah <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'Shenandoah',NA) 
estimateRappahannock <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'Rappahannock',NA) 
estimateYork <- subpopEstimate(finalData,'VSCIAll', 'SubBasin', 'York',NA) 
# By Ecoregion
estimatePiedmont <- subpopEstimate(finalData,'VSCIAll', 'EcoRegion', 'Piedmont',NA) 
estimateNorthernPiedmont <- subpopEstimate(finalData,'VSCIAll', 'EcoRegion', 'Northern Piedmont',NA) 
estimateCARV <- subpopEstimate(finalData,'VSCIAll', 'EcoRegion', 'Central Appalachian Ridges and Valleys',NA) 
estimateSEplains <- subpopEstimate(finalData,'VSCIAll', 'EcoRegion', 'Southeastern Plains',NA) 
estimateBRM <- subpopEstimate(finalData,'VSCIAll', 'EcoRegion', 'Blue Ridge Mountains',NA) 
estimateMountainBioregion <- subpopEstimate(finalData,'VSCIAll', 'BioRegion', 'Mountain','Mountain Bioregion')
estimatePiedmontBioregion <- subpopEstimate(finalData,'VSCIAll', 'BioRegion', 'Piedmont','Piedmont Bioregion') 
estimateCoastBioregion <- subpopEstimate(finalData,'VSCIAll', 'BioRegion', 'Coast', 'Coast Bioregion')
# By Order
estimateFirstOrder <- subpopEstimate(finalData,'VSCIAll', 'Order', '1','First Order') 
estimateSecondOrder <- subpopEstimate(finalData,'VSCIAll', 'Order', '2','Second Order') 
estimateThirdOrder <- subpopEstimate(finalData,'VSCIAll', 'Order', '3','Third Order') 
estimateFourthOrder <- subpopEstimate(finalData,'VSCIAll', 'Order', '4','Fourth Order') 
estimateFifthOrder <- subpopEstimate(finalData,'VSCIAll', 'Order', '5','Fifth Order') 
# By Basin Size
estimateBasin1 <- subpopEstimate(finalData,'VSCIAll', 'BasinSize', '1','<1 square mile') 
estimateBasin2 <- subpopEstimate(finalData,'VSCIAll', 'BasinSize', '2','1 to 10 square mile') 
estimateBasin3 <- subpopEstimate(finalData,'VSCIAll', 'BasinSize', '3','10 to 200 square mile') 
estimateBasin4 <- subpopEstimate(finalData,'VSCIAll', 'BasinSize', '4','>200 square mile') 
# By Year
estimate2001 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2001','Year 2001') 
estimate2002 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2002','Year 2002') 
estimate2003 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2003','Year 2003') 
estimate2004 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2004','Year 2004') 
estimate2005 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2005','Year 2005') 
estimate2006 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2006','Year 2006') 
estimate2007 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2007','Year 2007') 
estimate2008 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2008','Year 2008') 
estimate2009 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2009','Year 2009') 
estimate2010 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2010','Year 2010') 
estimate2011 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2011','Year 2011') 
estimate2012 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2012','Year 2012') 
estimate2013 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2013','Year 2013') 
estimate2014 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2014','Year 2014') 
#estimate2015 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2015','Year 2015') 
#estimate2016 <- subpopEstimate(finalData,'VSCIAll', 'Year', '2016','Year 2016') 
# By Phase
estimatePhase1 <- subpopEstimate(finalData,'VSCIAll', 'Panel', 'Phase1','Phase One 2001-2007') 
estimatePhase2 <- subpopEstimate(finalData,'VSCIAll', 'Panel', 'Phase2','Phase Two 2008-2014') 
estimateBay <- subpopEstimate(finalData,'VSCIAll', 'BayShed', 'Bay','Bay Watersheds 2001-2014')
# Bay/ NonBay
estimateNonBay <- subpopEstimate(finalData,'VSCIAll', 'BayShed', 'NonBay','Non-Bay Watersheds 2001-2014')
estimateBayPhase1 <- subpopEstimate(finalData,'VSCIAll', 'BayPanel', 'BayPhase1','Bay Watersheds 2001-2007')
estimateBayPhase2 <- subpopEstimate(finalData,'VSCIAll', 'BayPanel', 'BayPhase2','Bay Watersheds 2008-2014')
estimateNonBayPhase1 <- subpopEstimate(finalData,'VSCIAll', 'BayPanel', 'NonBayPhase1','Non-Bay Watersheds 2001-2007')
estimateNonBayPhase2 <- subpopEstimate(finalData,'VSCIAll', 'BayPanel', 'NonBayPhase2','Non-Bay Watersheds 2008-2014')
# Sampling Phase
estimateBioPanelPhase1 <- subpopEstimate(finalData,'VSCIAll', 'BioPanel', 'Phase1','VSCI Scores 2001-2003')
estimateBioPanelPhase2 <- subpopEstimate(finalData,'VSCIAll', 'BioPanel', 'Phase2','VSCI Scores 2004-2006')
estimateBioPanelPhase3 <- subpopEstimate(finalData,'VSCIAll', 'BioPanel', 'Phase3','VSCI Scores 2007-2010')
estimateBioPanelPhase4 <- subpopEstimate(finalData,'VSCIAll', 'BioPanel', 'Phase4','VSCI Scores 2011-2014')
# IR window
estimateIR2008 <- subpopEstimate(finalData,'VSCIAll', 'IR2008', '2008','IR2008')
estimateIR2010 <- subpopEstimate(finalData,'VSCIAll', 'IR2010', '2010','IR2010')
estimateIR2012<- subpopEstimate(finalData,'VSCIAll', 'IR2012', '2012','IR2012')
estimateIR2014 <- subpopEstimate(finalData,'VSCIAll', 'IR2014', '2014','IR2014')
estimateIR2016 <- subpopEstimate(finalData,'VSCIAll', 'IR2016', '2016','IR2016')
#estimateIR2018 <- subpopEstimate(finalData,'VSCIAll', 'IR2018', '2018','IR2018')
# Stream size
estimateSmall <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCat','Small',NA)
estimateMedium <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCat','Medium',NA)
estimateLarge <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCat','Large',NA)
# Stream size and sample phase
estimatePhase1Small <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCatPhase','Phase1Small',NA)
estimatePhase2Small <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCatPhase','Phase2Small',NA)
estimatePhase1Medium <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCatPhase','Phase1Medium',NA)
estimatePhase2Medium <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCatPhase','Phase2Medium',NA)
estimatePhase1Large <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCatPhase','Phase1Large',NA)
estimatePhase2Large <- subpopEstimate(finalData,'VSCIAll', 'StreamSizeCatPhase','Phase2Large',NA)


# now pull together all the data from each list
allTogether <- list(estimateVirginia,estimateRoanoke,estimateJames)
allTogether1 <- list(grep("estimate",names(.GlobalEnv),value=TRUE))
allTogether2 <- map(get(.GlobalEnv), grep("estimate",names(.GlobalEnv),value=TRUE))
allTogether3 <- do.call("list",mget(grep("estimate",names(.GlobalEnv),value=TRUE)))
allTogether4 <- invoke_map(list(grep("estimate",names(.GlobalEnv),value=TRUE)))
VSCI <- list(
  estimateVirginia = subpopEstimate(finalData,'VSCIAll',NA, 'Virginia',NA),
  # By Basin
  estimateRoanoke = subpopEstimate(finalData,'VSCIAll', 'Basin', 'Roanoke','Roanoke Basin'),
  estimateJames = subpopEstimate(finalData,'VSCIAll', 'Basin', 'James','James Basin') ,
  estimatePotomacShenandoah = subpopEstimate(finalData,'VSCIAll', 'Basin', 'Potomac-Shenandoah',NA) ,
  estimateRappahannockYork = subpopEstimate(finalData,'VSCIAll', 'Basin', 'Rappahannock-York',NA) 
  
)

CDFlist <- VSCI %>% map(1)
CDFdf <- suppressWarnings(allTogether5 %>% map_df(1))


CDFlist <- allTogether %>% map(1) 
CDFdf <- suppressWarnings(allTogether %>% map_df(1))
PCTlist <- allTogether %>% map(2)
PCTcombined <-  suppressWarnings(allTogether %>% map_df(2))



# Test to make sure everything matches
all.equal(statewideStatus[["CDF"]],estimateVirginia[["CDF"]])
all.equal(roanokestatus[["CDF"]],estimateRoanoke[["CDF"]])
all.equal(jamesstatus[["CDF"]],estimateJames[["CDF"]])
all.equal(potomacshenstatus[["CDF"]],estimatePotomacShenandoah[["CDF"]])
all.equal(rappyorkstatus[["CDF"]],estimateRappahannockYork[["CDF"]])
all.equal(newstatus[["CDF"]],estimateNew[["CDF"]])
all.equal(chowanstatus[["CDF"]],estimateChowan[["CDF"]])
all.equal(tennstatus[["CDF"]],estimateTennessee[["CDF"]])
all.equal(holstonstatus[["CDF"]],estimateHolston[["CDF"]])
all.equal(bigsandystatus[["CDF"]],estimateBigSandy[["CDF"]])
all.equal(clinchstatus[["CDF"]],estimateClinchPowell[["CDF"]])
all.equal(potomacstatus[["CDF"]],estimatePotomac[["CDF"]])
all.equal(shenstatus[["CDF"]],estimateShenandoah[["CDF"]])
all.equal(rappstatus[["CDF"]],estimateRappahannock[["CDF"]]) # doesnt match bc Rappahannock spelled wrong in biology.R code in rappstatus[["CDF"]]$Subpopulation
all.equal(yorkstatus[["CDF"]],estimateYork[["CDF"]])
all.equal(piedmontstatus[["CDF"]],estimatePiedmont[["CDF"]])
all.equal(npiedmontstatus[["CDF"]],estimateNorthernPiedmont[["CDF"]])
all.equal(carvstatus[["CDF"]],estimateCARV[["CDF"]])
all.equal(splainsstatus[["CDF"]],estimateSEplains[["CDF"]])
all.equal(brmstatus[["CDF"]],estimateBRM[["CDF"]])
all.equal(centappstatus[["CDF"]],estimateCentralAppalachians[["CDF"]])
all.equal(mountainstatus[["CDF"]],estimateMountainBioregion[["CDF"]])
all.equal(peidstatus[["CDF"]],estimatePiedmontBioregion[["CDF"]])
all.equal(coaststatus[["CDF"]],estimateCoastBioregion[["CDF"]])
all.equal(firstorderstatus[["CDF"]],estimateFirstOrder[["CDF"]])
all.equal(secondorderstatus[["CDF"]],estimateSecondOrder[["CDF"]])
all.equal(thirdorderstatus[["CDF"]],estimateThirdOrder[["CDF"]])
all.equal(fourthorderstatus[["CDF"]],estimateFourthOrder[["CDF"]])
all.equal(fifthorderstatus[["CDF"]],estimateFifthOrder[["CDF"]])
all.equal(basinonestatus[["CDF"]],estimateBasin1[["CDF"]])
all.equal(basintwostatus[["CDF"]],estimateBasin2[["CDF"]])
all.equal(basinthreestatus[["CDF"]],estimateBasin3[["CDF"]])
all.equal(basinfourstatus[["CDF"]],estimateBasin4[["CDF"]])
all.equal(yearonestatus[["CDF"]],estimate2001[["CDF"]])
all.equal(yeartwostatus[["CDF"]],estimate2002[["CDF"]])
all.equal(yearthreestatus[["CDF"]],estimate2003[["CDF"]])
all.equal(yearfourstatus[["CDF"]],estimate2004[["CDF"]])
all.equal(yearfivestatus[["CDF"]],estimate2005[["CDF"]])
all.equal(yearsixstatus[["CDF"]],estimate2006[["CDF"]])
all.equal(yearsevenstatus[["CDF"]],estimate2007[["CDF"]])
all.equal(yeareightstatus[["CDF"]],estimate2008[["CDF"]])
all.equal(yearninestatus[["CDF"]],estimate2009[["CDF"]])
all.equal(yeartenstatus[["CDF"]],estimate2010[["CDF"]])
all.equal(yearelevenstatus[["CDF"]],estimate2011[["CDF"]])
all.equal(yeartwelvestatus[["CDF"]],estimate2012[["CDF"]])
all.equal(yearthirteenstatus[["CDF"]],estimate2013[["CDF"]])
all.equal(yearfourteenstatus[["CDF"]],estimate2014[["CDF"]])
all.equal(phaseonestatus[["CDF"]],estimatePhase1[["CDF"]])
all.equal(phasetwostatus[["CDF"]],estimatePhase2[["CDF"]])
all.equal(baystatus[["CDF"]],estimateBay[["CDF"]])
all.equal(nonbaystatus[["CDF"]],estimateNonBay[["CDF"]])
all.equal(bayphaseonestatus[["CDF"]],estimateBayPhase1[["CDF"]])
all.equal(bayphasetwostatus[["CDF"]],estimateBayPhase2[["CDF"]])
all.equal(nonbayphaseonestatus[["CDF"]],estimateNonBayPhase1[["CDF"]])
all.equal(nonbayphasetwostatus[["CDF"]],estimateNonBayPhase2[["CDF"]])
all.equal(biophaseonestatus[["CDF"]],estimateBioPanelPhase1[["CDF"]])
all.equal(biophasetwostatus[["CDF"]],estimateBioPanelPhase2[["CDF"]])
all.equal(biophasethreestatus[["CDF"]],estimateBioPanelPhase3[["CDF"]])
all.equal(biophasefourstatus[["CDF"]],estimateBioPanelPhase4[["CDF"]])
all.equal(IR2008[["CDF"]],estimateIR2008[["CDF"]])
all.equal(IR2010[["CDF"]],estimateIR2010[["CDF"]])
all.equal(IR2012[["CDF"]],estimateIR2012[["CDF"]])
all.equal(IR2014[["CDF"]],estimateIR2014[["CDF"]])
all.equal(IR2016[["CDF"]],estimateIR2016[["CDF"]])
all.equal(Small[["CDF"]],estimateSmall[["CDF"]])
all.equal(Medium[["CDF"]],estimateMedium[["CDF"]])
all.equal(Large[["CDF"]],estimateLarge[["CDF"]])
all.equal(Small1[["CDF"]],estimatePhase1Small[["CDF"]])
all.equal(Small2[["CDF"]],estimatePhase2Small[["CDF"]])
all.equal(Medium1[["CDF"]],estimatePhase1Medium[["CDF"]])
all.equal(Medium2[["CDF"]],estimatePhase2Medium[["CDF"]])
all.equal(Large1[["CDF"]],estimatePhase1Large[["CDF"]])
all.equal(Large2[["CDF"]],estimatePhase2Large[["CDF"]])

