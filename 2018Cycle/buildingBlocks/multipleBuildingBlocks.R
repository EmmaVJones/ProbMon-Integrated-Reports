suppressPackageStartupMessages(library(tidyverse))#1.2.1
suppressPackageStartupMessages(library(sf))#0.6-1
suppressPackageStartupMessages(library(spsurvey))#3.3
suppressPackageStartupMessages(library(DT))#0.4


surveyData <- readxl::read_excel('originalData/ProbMonWadeable2001-2016.xlsx',sheet='ProbMonData2001-2016_EVJ') %>%
  select(StationID_Trend,Year,DO:MetalCCU,CALCIUM:MERCURY,wshdImpPCT) %>% 
  dplyr::rename(siteID=StationID_Trend,Ortho_P= "Ortho-P",
                X70331VFine=`70331VFine`)

## Match new parameter data with designStatus data
designStatus <- readxl::read_excel('originalData/biology.xlsx',sheet='biology2018') %>%
  mutate(#BioPanel= dplyr::recode(BioPanel,'Phase1'='BioPhase1','Phase2'='BioPhase2',
    #                         'Phase3'='BioPhase3','Phase4'='BioPhase4'),   # recode biophase
    IR2008 = replace(IR2008, IR2008==1, NA),
    IR2010 = replace(IR2010, IR2010==1, NA),
    IR2012 = replace(IR2012, IR2012==1, NA),
    IR2014 = replace(IR2014, IR2014==1, NA),
    IR2016 = replace(IR2016, IR2016==1, NA),
    IR2018 = replace(IR2018, IR2018==1, NA),
    Panel1 = ifelse(Panel == "Phase1", 1, NA),
    Panel2 = ifelse(Panel == "Phase2", 1, NA),
    BioPanel1 = ifelse(BioPanel == "Phase1", 1, NA),
    BioPanel2 = ifelse(BioPanel == "Phase2", 1, NA),
    BioPanel3 = ifelse(BioPanel == "Phase3", 1, NA),
    BioPanel4 = ifelse(BioPanel == "Phase4", 1, NA)) %>% # housekeeping: recode biophase and change 1's to NA's to indicate it wasnt sampled in that window, break up Panel and BioPanel windows to separate columns for weight adjustment purposes
  dplyr::rename(siteID=sampleID ) # start playing nicely with spsurvey)


## Start with VSCI then move through list

allCDFresults <- function(designStatus, surveyData, parameter){
  
  # Organize new parameter data
  parameterData <- data.frame(select(surveyData,siteID,parameter))
  parameterData[,2] <- as.numeric(as.character(parameterData[,2])) # change to numeric in case it was saved as anything else
  
  
  # Change status if it wasnt sampled
  initialWeightData <- left_join(designStatus,parameterData,by='siteID') %>%
    mutate(parameter_status=ifelse(status == "TS" & is.na(!!as.name(parameter)),"OT", # if TS but not sampled, change to OT
                                   ifelse(status == "OT" & is.na(!!as.name(parameter)),"OT", # if OT but not sampled, keep OT
                                          ifelse(status == "PD", "PD", # if PD, keep it
                                                 ifelse(status == 'NT','NT', # if NT, keep it
                                                        ifelse(status == "TS",'TS', NA))))), # if TS AND sampled, keep as TS
           designweightoriginal = as.factor(`strahler order`), # easier to change as factor
           designweightoriginal = dplyr::recode(designweightoriginal,`1`="3790.5165999999999",
                                                `2`="947.62919999999997",
                                                `3`="541.50239999999997",
                                                `4`="315.87639999999999", 
                                                `5`="140.3895", 
                                                `6`="140.3895")) # overwrite all design weights back to original
  
  # Adjust initial weights for trend stations across various data windows
  # First calculate number of times sampled in each window
  
  trendWeightAdjustments <- mutate(initialWeightData,
                                   designweightoriginal = as.numeric(as.character(designweightoriginal)), #factor to numeric
                                   siteIDoriginal=gsub("_.*$", "", siteID)) %>% # get rid of concatenated year for trends to make calculations easier
    # Full window
    group_by(siteIDoriginal) %>% 
    mutate(nYearsSampled = ifelse(is.na(siteIDoriginal),NA,n())) %>%
    ungroup()%>%
    # 2018 IR
    group_by(siteIDoriginal, IR2018) %>%
    mutate(nYearsSampled_IR2018 = ifelse(is.na(IR2018),NA,n())) %>%
    ungroup() %>%
    # 2016 IR
    group_by(siteIDoriginal, IR2016) %>%
    mutate(nYearsSampled_IR2016 = ifelse(is.na(IR2016),NA,n())) %>%
    ungroup()%>%
    # 2014 IR
    group_by(siteIDoriginal, IR2014) %>%
    mutate(nYearsSampled_IR2014 = ifelse(is.na(IR2014),NA,n()))%>%
    ungroup()%>%
    # 2012 IR
    group_by(siteIDoriginal, IR2012) %>%
    mutate(nYearsSampled_IR2012 = ifelse(is.na(IR2012),NA,n())) %>%
    ungroup()%>%
    # 2010 IR
    group_by(siteIDoriginal, IR2010) %>%
    mutate(nYearsSampled_IR2010 = ifelse(is.na(IR2010),NA,n())) %>%
    ungroup()%>%
    # 2008 IR
    group_by(siteIDoriginal, IR2008) %>%
    mutate(nYearsSampled_IR2008 = ifelse(is.na(IR2008),NA,n())) %>%
    # Panels
    group_by(siteIDoriginal, Panel1) %>%
    mutate(nYearsSampled_Panel1 = ifelse(is.na(Panel1),NA,n())) %>%
    ungroup()%>%
    group_by(siteIDoriginal, Panel2) %>%
    mutate(nYearsSampled_Panel2 = ifelse(is.na(Panel2),NA,n())) %>%
    ungroup()%>%
    # Biopanels
    group_by(siteIDoriginal, BioPanel1) %>%
    mutate(nYearsSampled_BioPanel1 = ifelse(is.na(BioPanel1),NA,n())) %>%
    ungroup() %>% 
    group_by(siteIDoriginal, BioPanel2) %>%
    mutate(nYearsSampled_BioPanel2 = ifelse(is.na(BioPanel2),NA,n())) %>%
    ungroup() %>% 
    group_by(siteIDoriginal, BioPanel3) %>%
    mutate(nYearsSampled_BioPanel3 = ifelse(is.na(BioPanel3),NA,n())) %>%
    ungroup() %>% 
    group_by(siteIDoriginal, BioPanel4) %>%
    mutate(nYearsSampled_BioPanel4 = ifelse(is.na(BioPanel4),NA,n())) %>%
    ungroup() %>% 
    # Divide out weight by years sampled in each window
    # if the site wasnt sampled in a window, give it the original weight but it will be adjusted according to an OT status later
    mutate(designweight_all= designweightoriginal/nYearsSampled,
           designweight_IR2018= ifelse(is.na(nYearsSampled_IR2018),designweightoriginal,designweightoriginal/nYearsSampled_IR2018),
           designweight_IR2016= ifelse(is.na(nYearsSampled_IR2016),designweightoriginal,designweightoriginal/nYearsSampled_IR2016),
           designweight_IR2014= ifelse(is.na(nYearsSampled_IR2014),designweightoriginal,designweightoriginal/nYearsSampled_IR2014),
           designweight_IR2012= ifelse(is.na(nYearsSampled_IR2012),designweightoriginal,designweightoriginal/nYearsSampled_IR2012),
           designweight_IR2010= ifelse(is.na(nYearsSampled_IR2010),designweightoriginal,designweightoriginal/nYearsSampled_IR2010),
           designweight_IR2008= ifelse(is.na(nYearsSampled_IR2008),designweightoriginal,designweightoriginal/nYearsSampled_IR2008),
           designweight_Panel1= ifelse(is.na(nYearsSampled_Panel1),designweightoriginal,designweightoriginal/nYearsSampled_Panel1),
           designweight_Panel2= ifelse(is.na(nYearsSampled_Panel2),designweightoriginal,designweightoriginal/nYearsSampled_Panel2),
           designweight_BioPanel1= ifelse(is.na(nYearsSampled_BioPanel1),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel1),
           designweight_BioPanel2= ifelse(is.na(nYearsSampled_BioPanel2),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel2),
           designweight_BioPanel3= ifelse(is.na(nYearsSampled_BioPanel3),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel3),
           designweight_BioPanel4= ifelse(is.na(nYearsSampled_BioPanel4),designweightoriginal,designweightoriginal/nYearsSampled_BioPanel4))
  
  
  ## Adjust design weights to get final weights
  
  # Initial sample frame inputs
  # List stream order by kilometer it represents
  sframe <- c('1st'=51210, '2nd'=13680, '3rd'=7781.08, '4th'=4448.257, 
              '5th'=1731.302, '6th'=163.901, '7th'=14.7099 )
  
  # recode to factor to make sframe match up to stream order
  trendWeightAdjustments$`strahler order` <- as.factor(trendWeightAdjustments$`strahler order`)
  levels(trendWeightAdjustments$`strahler order`) <- c('1st','2nd','3rd','4th','5th','6th','7th')
  
  finalWeights <- select(trendWeightAdjustments,siteID:IR2018,siteIDoriginal,designweightoriginal,designweight_all:designweight_BioPanel4,
                         parameter_status,!!as.name(parameter))
  finalWeights$finalweight_Year <- adjwgt(rep(TRUE,length(finalWeights$designweightoriginal)),finalWeights$designweightoriginal,
                                          finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_all <- adjwgt(rep(TRUE,length(finalWeights$designweight_all)),finalWeights$designweight_all,
                                         finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2018 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2018)),finalWeights$designweight_IR2018,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2016 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2016)),finalWeights$designweight_IR2016,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2014 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2014)),finalWeights$designweight_IR2014,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2012 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2012)),finalWeights$designweight_IR2012,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2010 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2010)),finalWeights$designweight_IR2010,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_IR2008 <- adjwgt(rep(TRUE,length(finalWeights$designweight_IR2008)),finalWeights$designweight_IR2008,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_Panel1 <- adjwgt(rep(TRUE,length(finalWeights$designweight_Panel1)),finalWeights$designweight_Panel1,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_Panel2 <- adjwgt(rep(TRUE,length(finalWeights$designweight_Panel2)),finalWeights$designweight_Panel2,
                                            finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel1 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel1)),finalWeights$designweight_BioPanel1,
                                               finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel2 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel2)),finalWeights$designweight_BioPanel2,
                                               finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel3 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel3)),finalWeights$designweight_BioPanel3,
                                               finalWeights$`strahler order`, sframe)
  finalWeights$finalweight_BioPanel4 <- adjwgt(rep(TRUE,length(finalWeights$designweight_BioPanel4)),finalWeights$designweight_BioPanel4,
                                               finalWeights$`strahler order`, sframe)
  
  
  
  # Change decimal degree coordinates to marinus for equal area projection (needed for spsurvey::cat.analysis())
  marinus <- spsurvey::marinus(finalWeights$`Latitude-DD`,finalWeights$`Longitude-DD`)
  finalWeights <- cbind(finalWeights,data.frame(xmarinus=marinus$x,ymarinus=marinus$y))
  rm(marinus)
  
  ####Stream Extent and Status Estimate
  
  siteExtent <- data.frame(siteID=finalWeights$siteID, Use=rep(TRUE,nrow(finalWeights)) )
  subpopExtent <- data.frame(siteID=finalWeights$siteID,Region=rep('Virginia',nrow(finalWeights)) )
  designExtent <- data.frame(siteID=finalWeights$siteID,stratum=rep(1,nrow(finalWeights)),
                             wgt=finalWeights$finalweight_all, xcoord=finalWeights$xmarinus,ycoord=finalWeights$ymarinus)
  
  StatusTNT <- finalWeights$status
  levels(StatusTNT) <- list(T=c('TS','PD','OT'), NT=c('NT') )
  
  data.cat.ext <- data.frame(finalWeights[,c('siteID','status')],StatusTNT)
  
  popstatus.est <- cat.analysis(sites = siteExtent, subpop = subpopExtent, design = designExtent,
                                data.cat = data.cat.ext, conf=95, vartype="Local")
  
  dataAnalyzed <- filter(finalWeights, parameter_status=='TS') %>%
    select(siteID,!!as.name(parameter),parameter_status,everything())
  
  output <- listOfResults(popstatus.est, dataAnalyzed, parameter)
  
  return(output)
}

LRBS <- allCDFresults(designStatus, surveyData,'LRBS')
LRBS_CDFdf <- suppressWarnings(LRBS[3:77] %>% map_df(1))
LRBS_PCTdf <-  suppressWarnings(LRBS[3:77] %>% map_df(2))


VSCI <- allCDFresults(designStatus, surveyData,'VSCIVCPMI')

DO <- allCDFresults(designStatus, surveyData,'DO')


# first go through each parameter and analyze all subpopulations
#allResults <- list()

for(i in 3:length(surveyData)){
  #allResults[[i]] <- suppressWarnings(allCDFresults(designStatus, surveyData,names(surveyData)[i]))
  #names(allResults)[i] <- names(surveyData)[i]
  # or to save all as separate objects
  assign(names(surveyData)[i], suppressWarnings(allCDFresults(designStatus, surveyData,names(surveyData)[i])))
}





# Then unpack each parameter to get a long df of CDF and PCT data
# Needs to be in separate loop so names of each object to be called will be in the environment
allCDF <- data.frame(Type=NA,Subpopulation=NA,Indicator=NA,Value=NA,NResp=NA,Estimate.P=NA,StdError.P=NA,LCB95Pct.P=NA,UCB95Pct.P=NA,   
                     Estimate.U=NA,StdError.U=NA,LCB95Pct.U=NA,UCB95Pct.U=NA)
allPCT <- data.frame(Type=NA,Subpopulation=NA,Indicator=NA,Statistic=NA,NResp=NA,Estimate=NA,StdError=NA,LCB95Pct=NA,UCB95Pct=NA)

for(i in 3:length(surveyData)){
  assign(paste(names(surveyData)[i],"CDFdf",sep="_"),suppressWarnings(get(names(surveyData)[i])[3:77] %>% map_df(1)))
  assign(paste(names(surveyData)[i],"PCTdf",sep="_"),suppressWarnings(get(names(surveyData)[i])[3:77] %>% map_df(2)))
  allCDF <- rbind(allCDF,get(paste(names(surveyData)[i],"CDFdf",sep="_")))
  allPCT <- rbind(allPCT,get(paste(names(surveyData)[i],"PCTdf",sep="_")))
  
}

write.csv(allCDF,'processedData/allCDF.csv',row.names = F)
write.csv(allPCT,'processedData/allPCT.csv',row.names = F)

