#Run in R3.4.4


# Presets for easier reports year to year
IRyear <- '2018'
IRyearWindowEnd <- '2016'
IRrange <- '2011 - 2016'





dat <- read.csv('processedData/allCDF.csv')# CDF results
surveyData <- readxl::read_excel('processedData/Wadeable_ProbMon_2001-2016_EVJ.xlsx',sheet='Wadeable_ProbMon_2001-2016_EVJ') %>%
  select(StationID_Trend,Year,DO:MetalCCU,CALCIUM:MERCURY,wshdImpPCT) %>% 
  dplyr::rename(siteID=StationID_Trend,Ortho_P= "Ortho-P",
                X70331VFine=`70331VFine`)
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




# Bring in relative risk data
rr <- read.csv('processedData/relriskIR2018.csv') %>%
  mutate(Stressor=dplyr::recode(Stressor,"TotHabstatus"="Habitat Disturbance",
                                "TDSstatus"='Ionic Strength',
                                "TNstatus"='Total Nitrogen',
                                "MetalCCUstatus"='Cumulative Dissolved Metals',
                                "LRBSstatus"='Streambed Sedimentation',
                                "TPstatus"='Total Phosphorus'))# %>%
  #arrange(desc(Estimate))
# rearrange factor level order based on descending Estimate, have to list backwards to get plot to
# show up in correct order
#rr$Stressor <- factor(rr$Stressor,levels=c('Total Phosphorus','Streambed Sedimentation',
#                                         'Cumulative Dissolved Metals','Total Nitrogen',
#                                         'Ionic Strength','Habitat Disturbance'))
#rr$Stressor = factor(rr$Stressor, levels=unique(rr$Stressor[order(rr$Estimate)]), ordered=TRUE)


# Bring in basin rank data
rankBasins <- read.csv('processedData/rankBasinsbyMicromap.csv') 
rankBasins1.4 <- read.csv('processedData/rankBasinsbyMicromap1-4.csv') 


# VLOOKUP (Excel function hack) by Julin Maloof
vlookup <- function(ref, #the value or values that you want to look for
                    table, #the table where you want to look for it; will look in first column
                    column, #the column that you want the return data to come from,
                    range=FALSE, #if there is not an exact match, return the closest?
                    larger=FALSE) #if doing a range lookup, should the smaller or larger key be used?)
{
  if(!is.numeric(column) & !column %in% colnames(table)) {
    stop(paste("can't find column",column,"in table"))
  }
  if(range) {
    if(!is.numeric(table[,1])) {
      stop(paste("The first column of table must be numeric when using range lookup"))
    }
    table <- table[order(table[,1]),] 
    index <- findInterval(ref,table[,1])
    if(larger) {
      index <- ifelse(ref %in% table[,1],index,index+1)
    }
    output <- table[index,column]
    output[!index <= dim(table)[1]] <- NA
    
  } else {
    output <- table[match(ref,table[,1]),column]
    output[!ref %in% table[,1]] <- NA #not needed?
  }
  dim(output) <- dim(ref)
  output
}
# Add line to ggplot function
addline_format <- function(x,...){
  gsub('\\s','\n',x)}

# Get data in correct format for micromap plots function
statslookup <- function(indicator,measure,category,revOrder){
  chowan <-  filter(dat,Subpopulation=='Chowan'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  rappahannock <-  filter(dat,Subpopulation=='Rappahannock'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  york <-  filter(dat,Subpopulation=='York'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  potomac <-  filter(dat,Subpopulation=='Potomac'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  shenandoah <-  filter(dat,Subpopulation=='Shenandoah'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  roanoke <-  filter(dat,Subpopulation=='Roanoke Basin'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  james <-  filter(dat,Subpopulation=='James Basin'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  new <-  filter(dat,Subpopulation=='New'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  bigsandy <-  filter(dat,Subpopulation=='Big Sandy'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  cp <-  filter(dat,Subpopulation=='Clinch-Powell'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  holston <-  filter(dat,Subpopulation=='Holston'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  Virginia <-  filter(dat,Subpopulation=='Virginia'&Indicator==indicator)%>%
    select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
  
  x <- data.frame(Subpopulation=c('Chowan','Rappahannock','York','Potomac','Shenandoah','Roanoke','James',
                                  'New','Big Sandy','Clinch-Powell','Holston','Virginia'),
                  Category=category,NResp=c(max(chowan$NResp),max(rappahannock$NResp),max(york$NResp),
                                            max(potomac$NResp),max(shenandoah$NResp),max(roanoke$NResp),
                                            max(james$NResp),max(new$NResp),max(bigsandy$NResp),max(cp$NResp),
                                            max(holston$NResp),max(Virginia$NResp)),
                  matchingdfname=c('chowan','rappahannock','york','potomac','shenandoah','roanoke','james','new','bigsandy','cp','holston','Virginia'))
  y=data.frame(Estimate.P=NA,LCB95Pct.P=NA, UCB95Pct.P=NA)
  for(i in 1:nrow(x)){
    y[i,1] <- vlookup(measure,get(as.character(x[i,4])),2,TRUE)
    y[i,2] <- vlookup(measure,get(as.character(x[i,4])),3,TRUE)
    y[i,3] <- vlookup(measure,get(as.character(x[i,4])),4,TRUE)
    y[is.na(y)] <- 0
    y[y>100] <- 100
  }
  y2 <- mutate(y,Estimate.P2=100-Estimate.P,LCB95Pct.P2=Estimate.P2-(Estimate.P-LCB95Pct.P),
               UCB95Pct.P2=Estimate.P2+(UCB95Pct.P-Estimate.P))%>%select(Estimate.P2,LCB95Pct.P2,UCB95Pct.P2)
  y2[y2>100] <- 100
  names(y) <- c(paste(indicator,"Estimate.P",sep=""),paste(indicator,"LCB95Pct.P",sep=""),paste(indicator,"UCB95Pct.P",sep=""))
  names(y2) <- c(paste(indicator,"Estimate.P",sep=""),paste(indicator,"LCB95Pct.P",sep=""),paste(indicator,"UCB95Pct.P",sep=""))
  if(revOrder==FALSE){z <- cbind(x,y)}else{z<-cbind(x,y2)}
  return(z)
}




# DO DATA
DO <- filter(dat,Subpopulation=='IR2018'&Indicator=='DO') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,Estimate.U,StdError.U,LCB95Pct.U,UCB95Pct.U,StdError.P)
totalsDO <- data.frame(Condition=c('Suboptimal'),
                       pct=vlookup(4,DO,2,TRUE),
                       ErrorLCB95=c(vlookup(4,DO,2,TRUE)-vlookup(4,DO,3,TRUE)),
                       ErrorUCB95=c(vlookup(4,DO,4,TRUE)-vlookup(4,DO,2,TRUE)))%>%
  mutate(Parameter='Dissolved Oxygen')%>%select(Parameter,everything())
DOsummary <- mutate(totalsDO,BelowStd=paste(formatC(pct,digits=1),"% ( +/- ",formatC(ErrorLCB95,digits=1),"% )", sep=""))%>%select(Parameter,BelowStd)
colnames(DOsummary)<-c('Parameter','Below Standard ( 4 mg/L )')



# PH DATA
pH <- filter(dat,Subpopulation=='IR2018'&Indicator=='pH') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,Estimate.U,StdError.U,LCB95Pct.U,UCB95Pct.U,StdError.P)
totalspH <- data.frame(Condition=c('Suboptimal'),
                       pct=vlookup(6,pH,2,TRUE),
                       ErrorLCB95=c(vlookup(6,pH,2,TRUE)-vlookup(6,pH,3,TRUE)),
                       ErrorUCB95=c(vlookup(6,pH,4,TRUE)-vlookup(6,pH,2,TRUE)))%>%
  mutate(Parameter='pH (Below 6)')%>%select(Parameter,everything())

pHsummary <- data.frame(Parameter='pH',Below=vlookup(6,pH,2,TRUE),Above=100-vlookup(9,pH,2,TRUE),
                        BelowError=vlookup(6,pH,2,TRUE)-vlookup(6,pH,3,TRUE),
                        AboveError=vlookup(9,pH,2,TRUE)-vlookup(9,pH,3,TRUE))
pHsummary <- mutate(pHsummary,BelowStd=paste(formatC(Below,digits=2),"% ( +/- ",formatC(BelowError,digits=2),"% )", sep=""),AboveStd=paste(formatC(Above,digits=2),"% ( +/- ",formatC(AboveError,digits=2),"% )", sep=""))%>%select(Parameter,BelowStd,AboveStd)
colnames(pHsummary)<-c('Parameter','Below Standard ( pH 6 )','Above Standard ( pH 9)')


# HABITAT DATA
hab <- filter(dat,Subpopulation=='IR2018'&Indicator=='TotHab') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P, NResp)
totalshab <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                        pct=c(vlookup(120,hab,2,TRUE), #suboptimal
                              vlookup(150,hab,2,TRUE)-vlookup(120,hab,2,TRUE), # fair
                              100-vlookup(150,hab,2,TRUE)), #optimal
                        # Error for middle ranges calculated by using midpoint of range error estimates
                        ErrorLCB95=c(vlookup(120,hab,2,TRUE)-vlookup(120,hab,3,TRUE),
                                     vlookup(135,hab,2,TRUE)-vlookup(135,hab,3,TRUE),
                                     vlookup(150,hab,2,TRUE)-vlookup(150,hab,3,TRUE)),
                        ErrorUCB95=c(vlookup(120,hab,4,TRUE)-vlookup(120,hab,2,TRUE),
                                     vlookup(135,hab,4,TRUE)-vlookup(135,hab,2,TRUE),
                                     vlookup(150,hab,4,TRUE)-vlookup(150,hab,2,TRUE))) %>%
  mutate(Parameter='Habitat Disturbance',
         n=c(vlookup(120,hab,5,TRUE),
             vlookup(150,hab,5,TRUE)-vlookup(120,hab,5,TRUE),
             max(hab$NResp)-vlookup(150,hab,5,TRUE)))%>%
  select(Parameter,everything())
totalshab$Condition <- factor(totalshab$Condition,levels = unique(totalshab$Condition))

totalshabsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                  pct=c(vlookup(120,hab,2,TRUE)),
                                  ErrorLCB95=c(vlookup(120,hab,2,TRUE)-vlookup(120,hab,3,TRUE)),
                                  ErrorUCB95=c(vlookup(120,hab,4,TRUE)-vlookup(120,hab,2,TRUE)))%>%
  mutate(Parameter='Habitat Disturbance')%>%select(Parameter,everything())


# LRBS DATA
LRBS <- filter(dat,Subpopulation=='IR2018'&Indicator=='LRBS') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
totalsLRBS <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                         pct=c(vlookup(-1,LRBS,2,TRUE), #suboptimal
                               #have to be creative here since two fair categories
                               sum((vlookup(-0.5,LRBS,2,TRUE)-vlookup(-1,LRBS,2,TRUE)), # fair soft
                                   (100- vlookup(0.5,LRBS,2,TRUE))), #fair hardening
                               vlookup(0.5,LRBS,2,TRUE)-vlookup(-0.5,LRBS,2,TRUE)),# optimal
                         # Error for middle ranges calculated by using midpoint of range error estimates
                         ErrorLCB95=c(vlookup(-1,LRBS,2,TRUE)-vlookup(-1,LRBS,3,TRUE),
                                      vlookup(-0.75,LRBS,2,TRUE)-vlookup(-0.75,LRBS,3,TRUE),
                                      vlookup(-0.5,LRBS,2,TRUE)-vlookup(-0.5,LRBS,3,TRUE)),
                         ErrorUCB95=c(vlookup(-1,LRBS,4,TRUE)-vlookup(-1,LRBS,2,TRUE),
                                      vlookup(-0.75,LRBS,4,TRUE)-vlookup(-0.75,LRBS,2,TRUE),
                                      vlookup(-0.5,LRBS,4,TRUE)-vlookup(-0.5,LRBS,2,TRUE))) %>%
  mutate(Parameter='Streambed Sedimentation',
         n=c(vlookup(-1,LRBS,5,TRUE),
             # again, had to be creative bc 2 fair ranges
             sum((vlookup(-0.5,LRBS,5,TRUE)-vlookup(-1,LRBS,5,TRUE)), # fair soft
                 (max(LRBS$NResp)-vlookup(0.5,LRBS,5,TRUE))), # fair hardening
             vlookup(0.5,LRBS,5,TRUE)-vlookup(-0.5,LRBS,5,TRUE)))%>%
  select(Parameter,everything())
totalsLRBS$Condition <- factor(totalsLRBS$Condition,levels = unique(totalsLRBS$Condition))

totalsLRBSsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                   pct=c(vlookup(-1,LRBS,2,TRUE)),
                                   ErrorLCB95=c(vlookup(-1,LRBS,2,TRUE)-vlookup(-1,LRBS,3,TRUE)),
                                   ErrorUCB95=c(vlookup(-1,LRBS,4,TRUE)-vlookup(-1,LRBS,2,TRUE)))%>%
  mutate(Parameter='Streambed Sedimentation')%>%select(Parameter,everything())


#TN DATA
TN <- filter(dat,Subpopulation=='IR2018'&Indicator=='TN') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,NResp)
totalsTN <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                       pct=c(100-vlookup(2,TN,2,TRUE),#suboptimal
                             vlookup(2,TN,2,TRUE)-vlookup(1,TN,2,TRUE), # fair
                             vlookup(1,TN,2,TRUE)),# optimal
                       # Error for middle ranges calculated by using midpoint of range error estimates
                       ErrorLCB95=c(vlookup(2,TN,2,TRUE)-vlookup(2,TN,3,TRUE),
                                    vlookup(1.5,TN,2,TRUE)-vlookup(1.5,TN,3,TRUE),
                                    vlookup(1,TN,2,TRUE)-vlookup(1,TN,3,TRUE)),
                       ErrorUCB95=c(vlookup(2,TN,4,TRUE)-vlookup(2,TN,2,TRUE),
                                    vlookup(1.5,TN,4,TRUE)-vlookup(1.5,TN,2,TRUE),
                                    vlookup(1,TN,4,TRUE)-vlookup(1,TN,2,TRUE))) %>%
  mutate(Parameter='Total Nitrogen',
         n=c(max(TN$NResp)-vlookup(2,TN,5,TRUE),
             vlookup(2,TN,5,TRUE) - vlookup(1,TN,5,TRUE),
             vlookup(1,TN,5,TRUE))) %>%
  select(Parameter,everything())
totalsTN$Condition <- factor(totalsTN$Condition,levels = unique(totalsTN$Condition))

totalsTNsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                 pct=c(100-vlookup(2,TN,2,TRUE)),
                                 ErrorLCB95=c(vlookup(2,TN,2,TRUE)-vlookup(2,TN,3,TRUE)),
                                 ErrorUCB95=c(vlookup(2,TN,4,TRUE)-vlookup(2,TN,2,TRUE)))%>%
  mutate(Parameter='Total Nitrogen')%>%select(Parameter,everything())

#TP DATA
TP <- filter(dat,Subpopulation=='IR2018'&Indicator=='TP') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P, NResp)
totalsTP <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                       pct=c(100-vlookup(0.05,TP,2,TRUE),#suboptimal
                             vlookup(0.05,TP,2,TRUE)-vlookup(0.02,TP,2,TRUE), # fair
                             vlookup(0.02,TP,2,TRUE)), #optimal
                       # Error for middle ranges calculated by using midpoint of range error estimates
                       ErrorLCB95=c(vlookup(0.05,TP,2,TRUE)-vlookup(0.05,TP,3,TRUE),
                                    vlookup(0.035,TP,2,TRUE)-vlookup(0.035,TP,3,TRUE),
                                    vlookup(0.02,TP,2,TRUE)-vlookup(0.02,TP,3,TRUE)),
                       ErrorUCB95=c(vlookup(0.05,TP,4,TRUE)-vlookup(0.05,TP,2,TRUE),
                                    vlookup(0.035,TP,4,TRUE)-vlookup(0.035,TP,2,TRUE),
                                    vlookup(0.02,TP,4,TRUE)-vlookup(0.02,TP,2,TRUE))) %>%
  mutate(Parameter='Total Phosphorus',
         n=c(max(TP$NResp)-vlookup(0.05,TP,5,TRUE),
             vlookup(0.05,TP,5,TRUE) - vlookup(0.02,TP,5,TRUE),
             vlookup(0.02,TP,5,TRUE)))%>%
  select(Parameter,everything())
totalsTP$Condition <- factor(totalsTP$Condition,levels = unique(totalsTP$Condition))

totalsTPsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                 pct=c(100-vlookup(0.05,TP,2,TRUE)),
                                 ErrorLCB95=c(vlookup(0.05,TP,2,TRUE)-vlookup(0.05,TP,3,TRUE)),
                                 ErrorUCB95=c(vlookup(0.05,TP,4,TRUE)-vlookup(0.05,TP,2,TRUE)))%>%
  mutate(Parameter='Total Phosphorus')%>%select(Parameter,everything())

# TDS DATA
TDS <- filter(dat,Subpopulation=='IR2018'&Indicator=='TDS') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P, NResp)
totalsTDS <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                        pct=c(100-vlookup(350,TDS,2,TRUE),#suboptimal
                              vlookup(350,TDS,2,TRUE)-vlookup(100,TDS,2,TRUE), # fair
                              vlookup(100,TDS,2,TRUE)), #optimal
                        # Error for middle ranges calculated by using midpoint of range error estimates
                        ErrorLCB95=c(vlookup(350,TDS,2,TRUE)-vlookup(350,TDS,3,TRUE),
                                     vlookup(225,TDS,2,TRUE)-vlookup(225,TDS,3,TRUE),
                                     vlookup(100,TDS,2,TRUE)-vlookup(100,TDS,3,TRUE)),
                        ErrorUCB95=c(vlookup(350,TDS,4,TRUE)-vlookup(350,TDS,2,TRUE),
                                     vlookup(225,TDS,4,TRUE)-vlookup(225,TDS,2,TRUE),
                                     vlookup(100,TDS,4,TRUE)-vlookup(100,TDS,2,TRUE))) %>%
  mutate(Parameter='Ionic Strength',
         n=c(max(TDS$NResp)-vlookup(350,TDS,5,TRUE),
             vlookup(350,TDS,5,TRUE) - vlookup(100,TDS,5,TRUE),
             vlookup(100,TDS,5,TRUE))) %>%
  select(Parameter,everything())
totalsTDS$Condition <- factor(totalsTDS$Condition,levels = unique(totalsTDS$Condition))

totalsTDSsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                  pct=c(100-vlookup(350,TDS,2,TRUE)),
                                  ErrorLCB95=c(vlookup(350,TDS,2,TRUE)-vlookup(350,TDS,3,TRUE)),
                                  ErrorUCB95=c(vlookup(350,TDS,4,TRUE)-vlookup(350,TDS,2,TRUE)))%>%
  mutate(Parameter='Ionic Strength')%>%select(Parameter,everything())

# METALS CCU DATA
mCCU <- filter(dat,Subpopulation=='IR2018'&Indicator=='MetalCCU') %>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P, NResp)
totalsmCCU <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                         pct=c(100-vlookup(2,mCCU,2,TRUE),#suboptimal
                               vlookup(2,mCCU,2,TRUE)-vlookup(1,mCCU,2,TRUE), # fair
                               vlookup(1,mCCU,2,TRUE)), #optimal
                         # Error for middle ranges calculated by using midpoint of range error estimates
                         ErrorLCB95=c(vlookup(2,mCCU,2,TRUE)-vlookup(2,mCCU,3,TRUE),
                                      vlookup(1.5,mCCU,2,TRUE)-vlookup(1.5,mCCU,3,TRUE),
                                      vlookup(1,mCCU,2,TRUE)-vlookup(1,mCCU,3,TRUE)),
                         ErrorUCB95=c(vlookup(2,mCCU,4,TRUE)-vlookup(2,mCCU,2,TRUE),
                                      vlookup(1.5,mCCU,4,TRUE)-vlookup(1.5,mCCU,2,TRUE),
                                      vlookup(1,mCCU,4,TRUE)-vlookup(1,mCCU,2,TRUE))) %>%
  mutate(Parameter='Cumulative Dissolved Metals',
         n=c(max(mCCU$NResp)-vlookup(2,mCCU,5,TRUE),
             vlookup(2,mCCU,5,TRUE) - vlookup(1,mCCU,5,TRUE),
             vlookup(1,mCCU,5,TRUE))) %>%
  select(Parameter,everything())
totalsmCCU$Condition <- factor(totalsmCCU$Condition,levels = unique(totalsmCCU$Condition))

totalsmCCUsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                   pct=c(100-vlookup(2,mCCU,2,TRUE)),
                                   ErrorLCB95=c(vlookup(2,mCCU,2,TRUE)-vlookup(2,mCCU,3,TRUE)),
                                   ErrorUCB95=c(vlookup(2,mCCU,4,TRUE)-vlookup(2,mCCU,2,TRUE)))%>%
  mutate(Parameter='Cumulative Dissolved Metals')%>%select(Parameter,everything())


# VSCI DATA
VSCI <- filter(dat,Subpopulation=='IR2018'&Indicator=='VSCIVCPMI')%>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,Estimate.U,StdError.U,LCB95Pct.U,UCB95Pct.U,StdError.P)
totalsVSCI <- data.frame(Condition=c('Suboptimal'),
                         pct=vlookup(60,VSCI,2,TRUE),
                         ErrorLCB95=c(vlookup(60,VSCI,2,TRUE)-vlookup(60,VSCI,3,TRUE)),
                         ErrorUCB95=c(vlookup(60,VSCI,4,TRUE)-vlookup(60,VSCI,2,TRUE)))%>%
  mutate(Parameter='VSCI/VCPMI (Biomonitoring)')%>%select(Parameter,everything())
#VSCIsummary <- mutate(totalsVSCI,BelowStd=paste(formatC(pct,digits=2),"% ( +/- ",formatC(ErrorLCB95,digits=2),"% )", sep=""))%>%select(Parameter,BelowStd)
#colnames(VSCIsummary)<-c('Parameter','Below Standard ( 4 mg/L )')


# DNickel DATA
Nickel <- filter(dat,Subpopulation=='IR2018'&Indicator=='NICKEL')%>%
  select(Value,Estimate.P,LCB95Pct.P,UCB95Pct.P,Estimate.U,StdError.U,LCB95Pct.U,UCB95Pct.U,StdError.P)
totalsNickel <- data.frame(Condition=c('Suboptimal'),
                           pct=vlookup(0.03,Nickel,2,TRUE),
                           ErrorLCB95=c(vlookup(0.03,Nickel,2,TRUE)-vlookup(0.03,Nickel,3,TRUE)),
                           ErrorUCB95=c(vlookup(0.03,Nickel,4,TRUE)-vlookup(0.03,Nickel,2,TRUE)))%>%
  mutate(Parameter='Dissolved Nickel')%>%select(Parameter,everything())


# initial Graph, may need to change the order of factor levels each IR window based on %
paramsummary <- rbind(totalsDO,totalshabsuboptimal,totalsLRBSsuboptimal,totalsTNsuboptimal,totalspH,totalsTPsuboptimal,
                      totalsTDSsuboptimal,totalsmCCUsuboptimal,totalsVSCI)%>%
  mutate(Standard=c('yes','no','no','no','yes','no','no','no','yes'))
paramsummary$Parameter <- c('Dissolved Oxygen','Habitat\nDisturbance','Streambed\nSedimentation','Total Nitrogen',
                            'pH\n(Below 6)','Total\nPhosphorus','Ionic\nStrength','Cumulative \nDissolved Metals','VSCI/VCPMI\n(Biomonitoring)')
paramsummary$Parameter <- factor(paramsummary$Parameter)
paramsummary$Parameter <- factor(paramsummary$Parameter,
                                 levels=c('Dissolved Oxygen','Cumulative \nDissolved Metals',
                                          'Total Nitrogen','pH\n(Below 6)','Ionic\nStrength','Habitat\nDisturbance',
                                          'Total\nPhosphorus','Streambed\nSedimentation','VSCI/VCPMI\n(Biomonitoring)'))

# Stressor extent
stressorext <- rbind(totalsLRBSsuboptimal,totalsTNsuboptimal,totalsTPsuboptimal,totalsTDSsuboptimal,totalsmCCUsuboptimal,totalshabsuboptimal) %>%
  arrange(desc(pct))
# rearrange factor level order based on descending pct, have to list backwards to get plot to
# show up in correct order
stressorext$Parameter <- factor(stressorext$Parameter,
                                levels=c('Cumulative Dissolved Metals','Total Nitrogen','Ionic Strength',
                                         'Habitat Disturbance','Total Phosphorus','Streambed Sedimentation'))
                                 

## Play with stacked bar charts
together <- rbind(totalshab,totalsLRBS,totalsTN,totalsTP,totalsTDS,totalsmCCU)

togetherstacked <- together
togetherstacked$Parameter <- factor(togetherstacked$Parameter,
                                    levels=c('Streambed Sedimentation','Total Phosphorus',
                                             'Habitat Disturbance','Total Nitrogen',
                                             'Cumulative Dissolved Metals','Ionic Strength'))

togetherstackedflip <- together
togetherstackedflip$Parameter <- factor(togetherstackedflip$Parameter,
                                        levels=c('Ionic Strength','Cumulative Dissolved Metals',
                                                 'Total Nitrogen','Habitat Disturbance',
                                                 'Total Phosphorus','Streambed Sedimentation'))


# RELATIVE RISK PLOT
rr$Stressor <- as.factor(c('Habitat Disturbance', 'Total Nitrogen','Total Phosphorus','Ionic Strength','Cumulative Dissolved Metals','Streambed Sedimentation'))
rr$Stressor <-  factor(rr$Stressor,
                       levels=c('Total Phosphorus','Cumulative Dissolved Metals','Total Nitrogen','Streambed Sedimentation','Ionic Strength','Habitat Disturbance'))

