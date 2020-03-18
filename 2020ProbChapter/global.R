#Run in R3.6.1


# Presets for easier reports year to year
IRyear <- '2020'
IRyearWindowBegin <- '2013'
IRyearWindowEnd <- '2018'
IRrange <- '2013 - 2018'





dat <- read.csv('processedData/allCDF.csv')# CDF results
surveyData <- read_csv('processedData/Wadeable_ProbMon_2001-2018.csv') %>%
  select(StationID_Trend,Year,DO:MetalCCU,CALCIUM:MERCURY,wshdImpPCT) %>% 
  dplyr::rename(siteID=StationID_Trend,Ortho_P= "Ortho-P",
                X70331VFine=`70331VFine`)
designStatus <- readxl::read_excel('processedData/biology20012018final.xlsx',sheet='biology2020') %>%
  mutate(
    IR2008 = replace(IR2008, IR2008==1, NA),
    IR2010 = replace(IR2010, IR2010==1, NA),
    IR2012 = replace(IR2012, IR2012==1, NA),
    IR2014 = replace(IR2014, IR2014==1, NA),
    IR2016 = replace(IR2016, IR2016==1, NA),
    IR2018 = replace(IR2018, IR2018==1, NA),
    IR2020 = replace(IR2020, IR2020==1, NA),
    Panel1 = ifelse(Panel == "Phase1", 1, NA),
    Panel2 = ifelse(Panel == "Phase2", 1, NA),
    BioPanel1 = ifelse(BioPanel == "Phase1", 1, NA),
    BioPanel2 = ifelse(BioPanel == "Phase2", 1, NA),
    BioPanel3 = ifelse(BioPanel == "Phase3", 1, NA),
    BioPanel4 = ifelse(BioPanel == "Phase4", 1, NA)) %>% # housekeeping: recode biophase and change 1's to NA's to indicate it wasnt sampled in that window, break up Panel and BioPanel windows to separate columns for weight adjustment purposes
  dplyr::rename(siteID=sampleID ) # start playing nicely with spsurvey)



# Bring in relative risk data
rr <- read.csv('processedData/relriskIR2020.csv') %>%
  mutate(Stressor=dplyr::recode(Stressor,"TotHabstatus"="Habitat Disturbance",
                                "TDSstatus"='Ionic Strength',
                                "TNstatus"='Total Nitrogen',
                                "MetalCCUstatus"='Cumulative Dissolved Metals',
                                "LRBSstatus"='Streambed Sedimentation',
                                "TPstatus"='Total Phosphorus')) %>%
  mutate(MoE = StdError.log *1.96)


# VLOOKUP (Excel function hack) by Julin Maloof
vlookup <- function(ref, #the value or values that you want to look for
                    table, #the table where you want to look for it; will look in first column
                    column, #the column that you want the return data to come from,
                    range=FALSE, #if there is not an exact match, return the closest?
                    larger=FALSE) #if doing a range lookup, should the smaller or larger key be used?)
{
  # 2020 addition, make tibbles dataframes
  table <- as.data.frame(table)
  
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

#indicator <- 'TotHab'
#measure <- 120
#category <- 'SubOptimal'
#revOrder <- FALSE

statslookup <- function(indicator,measure,category,revOrder){
  chowan <-  filter(dat,Subpopulation=='Chowan'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  rappahannock <-  filter(dat,Subpopulation=='Rappahannock'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  york <-  filter(dat,Subpopulation=='York'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  potomac <-  filter(dat,Subpopulation=='Potomac'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  shenandoah <-  filter(dat,Subpopulation=='Shenandoah'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  roanoke <-  filter(dat,Subpopulation=='Roanoke Basin'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  james <-  filter(dat,Subpopulation=='James Basin'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  new <-  filter(dat,Subpopulation=='New'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  bigsandy <-  filter(dat,Subpopulation=='Big Sandy'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  cp <-  filter(dat,Subpopulation=='Clinch-Powell'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  holston <-  filter(dat,Subpopulation=='Holston'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  Virginia <-  filter(dat,Subpopulation=='Virginia'&Indicator==indicator)%>%
    select(Value,Estimate.P,StdError.P,NResp) %>%
    mutate(MoE = StdError.P * 1.96)
  
  x <- data.frame(Subpopulation=c('Chowan','Rappahannock','York','Potomac','Shenandoah','Roanoke','James',
                                  'New','Big Sandy','Clinch-Powell','Holston','Virginia'),
                  Category=category,NResp=c(max(chowan$NResp),max(rappahannock$NResp),max(york$NResp),
                                            max(potomac$NResp),max(shenandoah$NResp),max(roanoke$NResp),
                                            max(james$NResp),max(new$NResp),max(bigsandy$NResp),max(cp$NResp),
                                            max(holston$NResp),max(Virginia$NResp)),
                  matchingdfname=c('chowan','rappahannock','york','potomac','shenandoah','roanoke','james','new','bigsandy','cp','holston','Virginia'))
  y=data.frame(Estimate.P=NA, MoE = NA)
  for(i in 1:nrow(x)){
    y[i,1] <- vlookup(measure,get(as.character(x[i,4])),2,TRUE) # get percentile
    y[i,2] <- vlookup(measure,get(as.character(x[i,4])),5,TRUE) # get MoE
    y[is.na(y)] <- 0
    y[y>100] <- 100
  }
  
  y <- mutate(y, LCB95Pct.P = Estimate.P-(MoE),
              UCB95Pct.P=Estimate.P+(MoE))
  
  y2 <- mutate(y, Estimate.P2=100-Estimate.P, 
               LCB95Pct.P2=Estimate.P2-(MoE),
               UCB95Pct.P2=Estimate.P2+(MoE)) %>%
    select(Estimate.P2,LCB95Pct.P2,UCB95Pct.P2)
  y2[y2>100] <- 100
  y <- select(y, -MoE)
  names(y) <- c(paste(indicator,"Estimate.P",sep=""),paste(indicator,"LCB95Pct.P",sep=""),paste(indicator,"UCB95Pct.P",sep=""))
  names(y2) <- c(paste(indicator,"Estimate.P",sep=""),paste(indicator,"LCB95Pct.P",sep=""),paste(indicator,"UCB95Pct.P",sep=""))
  if(revOrder==FALSE){z <- cbind(x,y)}else{z<-cbind(x,y2)}
  return(z)
}




# DO DATA
DO <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'DO') %>%
  select(Value, Estimate.P, StdError.P)
totalsDO <- data.frame(Condition=c('Suboptimal'),
                       pct = vlookup(4,DO,2,TRUE),
                       MoE = vlookup(4, DO, 3, TRUE) * 1.96) %>%
  replace_na(list(pct = 0, MoE = 0)) %>% # change NA to 0 if necessary
  mutate(Parameter='Dissolved Oxygen') %>%
  select(Parameter,everything())
DOsummary <- mutate(totalsDO, BelowStd= paste(formatC(pct,digits=1),
                                              "% ( +/- ",formatC(MoE,digits=1),"% )", sep="")) %>%
  select(Parameter,BelowStd)
colnames(DOsummary)<-c('Parameter','Below Standard ( 4 mg/L )')



# PH DATA
pH <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'pH') %>%
  select(Value, Estimate.P, StdError.P)
totalspH <- data.frame(Condition=c('Suboptimal'),
                       pct = vlookup(6,pH,2,TRUE),
                       MoE = vlookup(6, pH, 3, TRUE) * 1.96) %>%
  replace_na(list(pct = 0, MoE = 0)) %>% # change NA to 0 if necessary
  mutate(Parameter='pH\n(Below 6)') %>%
  select(Parameter,everything())

pHsummary <- data.frame(Parameter='pH', Below = vlookup(6,pH,2,TRUE),
                        Above = 100-vlookup(9,pH,2,TRUE),
                        BelowError = vlookup(6,pH,3,TRUE)* 1.96,
                        AboveError = vlookup(9,pH,3,TRUE)* 1.96)
pHsummary <- mutate(pHsummary, BelowStd = paste(formatC(Below,digits=2),
                                                "% ( +/- ",formatC(BelowError,digits=2),"% )", sep=""),
                    AboveStd = paste(formatC(Above,digits=2),
                                     "% ( +/- ",formatC(AboveError,digits=2),"% )", sep="")) %>%
  select(Parameter,BelowStd,AboveStd)
colnames(pHsummary)<-c('Parameter','Below Standard ( pH 6 )','Above Standard ( pH 9)')


# HABITAT DATA
hab <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'TotHab') %>%
  select(Value,Estimate.P,StdError.P, NResp)
totalshab <- data.frame(Condition = c('Suboptimal','Fair','Optimal'),
                        pct = c(vlookup(120,hab,2,TRUE), #suboptimal
                                vlookup(150,hab,2,TRUE)-vlookup(120,hab,2,TRUE), # fair
                                100-vlookup(150,hab,2,TRUE)), #optimal
                        # Error for middle ranges calculated by using midpoint of range error estimates
                        MoE = c(vlookup(120,hab,3,TRUE) * 1.96,
                                vlookup(135,hab,3,TRUE) * 1.96,
                                vlookup(150,hab,3,TRUE) * 1.96)) %>%
  mutate(Parameter='Habitat Disturbance',
         n = c(vlookup(120,hab,4,TRUE),
             vlookup(150,hab,4,TRUE)-vlookup(120,hab,4,TRUE),
             max(hab$NResp)-vlookup(150,hab,4,TRUE))) %>%
  select(Parameter,everything())
totalshab$Condition <- factor(totalshab$Condition,levels = unique(totalshab$Condition))

totalshabsuboptimal <- data.frame(Condition=c('Suboptimal'),
                                  pct=c(vlookup(120,hab,2,TRUE)),
                                  MoE = c(vlookup(120, hab, 3, TRUE) * 1.96)) %>%
  mutate(Parameter='Habitat\nDisturbance') %>%
  select(Parameter,everything())


# LRBS DATA
LRBS <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'LRBS') %>%
  select(Value, Estimate.P, StdError.P, NResp)
totalsLRBS <- data.frame(Condition=c('Suboptimal','Fair','Optimal'),
                         pct= c(vlookup(-1,LRBS,2,TRUE), #suboptimal
                                #have to be creative here since two fair categories
                                sum((vlookup(-0.5,LRBS,2,TRUE)-vlookup(-1,LRBS,2,TRUE)), # fair soft
                                    (100- vlookup(0.5,LRBS,2,TRUE))), #fair hardening
                                vlookup(0.5,LRBS,2,TRUE)-vlookup(-0.5,LRBS,2,TRUE)),# optimal
                         # Error for middle ranges calculated by using midpoint of range error estimates
                         MoE = c(vlookup(-1,LRBS,3,TRUE) * 1.96,
                                 vlookup(-0.75,LRBS,3,TRUE) * 1.96,
                                 vlookup(-0.5,LRBS,3,TRUE) * 1.96)) %>%
  mutate(Parameter = 'Streambed\nSedimentation',
         n = c(vlookup(-1,LRBS,4,TRUE),
               # again, had to be creative bc 2 fair ranges
               sum((vlookup(-0.5,LRBS,4,TRUE)-vlookup(-1,LRBS,4,TRUE)), # fair soft
                   (max(LRBS$NResp)-vlookup(0.5,LRBS,4,TRUE))), # fair hardening
               vlookup(0.5,LRBS,4,TRUE)-vlookup(-0.5,LRBS,4,TRUE))) %>%
  select(Parameter,everything())
totalsLRBS$Condition <- factor(totalsLRBS$Condition,levels = unique(totalsLRBS$Condition))

totalsLRBSsuboptimal <- data.frame(Condition = c('Suboptimal'),
                                   pct = c(vlookup(-1,LRBS,2,TRUE)),
                                   MoE = c(vlookup(-1,LRBS,3,TRUE)*1.96)) %>%
  mutate(Parameter='Streambed Sedimentation') %>%
  select(Parameter,everything())


#TN DATA
TN <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'TN') %>%
  select(Value,Estimate.P, StdError.P, NResp)
totalsTN <- data.frame(Condition = c('Suboptimal','Fair','Optimal'),
                       pct = c(100-vlookup(2,TN,2,TRUE),#suboptimal
                               vlookup(2,TN,2,TRUE)-vlookup(1,TN,2,TRUE), # fair
                               vlookup(1,TN,2,TRUE)),# optimal
                       # Error for middle ranges calculated by using midpoint of range error estimates
                       MoE = c(vlookup(2,TN,3,TRUE) * 1.96, 
                               vlookup(1.5,TN,3,TRUE)* 1.96,
                               vlookup(1,TN,3,TRUE)* 1.96)) %>%
  mutate(Parameter = 'Total Nitrogen',
         n=c(max(TN$NResp)-vlookup(2,TN,4,TRUE),
             vlookup(2,TN,4,TRUE) - vlookup(1,TN,4,TRUE),
             vlookup(1,TN,4,TRUE))) %>%
  select(Parameter,everything())
totalsTN$Condition <- factor(totalsTN$Condition,levels = unique(totalsTN$Condition))

totalsTNsuboptimal <- data.frame(Condition = c('Suboptimal'),
                                 pct = c(100-vlookup(2,TN,2,TRUE)),
                                 MoE = c(vlookup(2,TN,3,TRUE)*1.96)) %>%
  mutate(Parameter='Total Nitrogen')%>%select(Parameter,everything())

#TP DATA
TP <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'TP') %>%
  select(Value, Estimate.P, StdError.P, NResp)
totalsTP <- data.frame(Condition = c('Suboptimal','Fair','Optimal'),
                       pct = c(100-vlookup(0.05,TP,2,TRUE),#suboptimal
                               vlookup(0.05,TP,2,TRUE)-vlookup(0.02,TP,2,TRUE), # fair
                               vlookup(0.02,TP,2,TRUE)), #optimal
                       # Error for middle ranges calculated by using midpoint of range error estimates
                       MoE = c(vlookup(0.05,TP,3,TRUE) * 1.96, 
                               vlookup(0.035,TP,3,TRUE)* 1.96,
                               vlookup(0.02,TP,3,TRUE)* 1.96)) %>%
  mutate(Parameter = 'Total\nPhosphorus',
         n = c(max(TP$NResp)-vlookup(0.05,TP,4,TRUE),
               vlookup(0.05,TP,4,TRUE) - vlookup(0.02,TP,4,TRUE),
               vlookup(0.02,TP,4,TRUE))) %>%
  select(Parameter,everything())
totalsTP$Condition <- factor(totalsTP$Condition,levels = unique(totalsTP$Condition))

totalsTPsuboptimal <- data.frame(Condition = c('Suboptimal'),
                                 pct = c(100-vlookup(0.05,TP,2,TRUE)),
                                 MoE = c(vlookup(0.05,TP,3,TRUE)*1.96)) %>%
  mutate(Parameter='Total Phosphorus')%>%select(Parameter,everything())

# TDS DATA
TDS <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'TDS') %>%
  select(Value, Estimate.P, StdError.P, NResp)
totalsTDS <- data.frame(Condition = c('Suboptimal','Fair','Optimal'),
                        pct = c(100-vlookup(350,TDS,2,TRUE),#suboptimal
                                vlookup(350,TDS,2,TRUE)-vlookup(100,TDS,2,TRUE), # fair
                                vlookup(100,TDS,2,TRUE)), #optimal
                        # Error for middle ranges calculated by using midpoint of range error estimates
                        MoE = c(vlookup(350,TDS,3,TRUE) * 1.96, 
                                vlookup(225,TDS,3,TRUE)* 1.96,
                                vlookup(100,TDS,3,TRUE)* 1.96)) %>%
  mutate(Parameter = 'Ionic\nStrength',
         n = c(max(TDS$NResp)-vlookup(350,TDS,4,TRUE),
               vlookup(350,TDS,4,TRUE) - vlookup(100,TDS,4,TRUE),
               vlookup(100,TDS,4,TRUE))) %>%
  select(Parameter, everything())
totalsTDS$Condition <- factor(totalsTDS$Condition,levels = unique(totalsTDS$Condition))

totalsTDSsuboptimal <- data.frame(Condition = c('Suboptimal'),
                                  pct = c(100-vlookup(350,TDS,2,TRUE)),
                                  MoE = c(vlookup(350,TDS,3,TRUE) * 1.96)) %>%
  mutate(Parameter = 'Ionic Strength') %>% 
  select(Parameter, everything())

# METALS CCU DATA
mCCU <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'MetalCCU') %>%
  select(Value, Estimate.P, StdError.P, NResp)
totalsmCCU <- data.frame(Condition = c('Suboptimal','Fair','Optimal'),
                         pct = c(100-vlookup(2,mCCU,2,TRUE),#suboptimal
                                 vlookup(2,mCCU,2,TRUE)-vlookup(1,mCCU,2,TRUE), # fair
                                 vlookup(1,mCCU,2,TRUE)), #optimal
                         # Error for middle ranges calculated by using midpoint of range error estimates
                         MoE = c(vlookup(2,mCCU,3,TRUE) * 1.96, 
                                 vlookup(1.5,mCCU,3,TRUE)* 1.96,
                                 vlookup(1,mCCU,3,TRUE)* 1.96)) %>%
  mutate(Parameter = 'Cumulative \nDissolved Metals',
         n = c(max(mCCU$NResp)-vlookup(2,mCCU,4,TRUE),
               vlookup(2,mCCU,4,TRUE) - vlookup(1,mCCU,4,TRUE),
               vlookup(1,mCCU,4,TRUE))) %>%
  select(Parameter,everything())
totalsmCCU$Condition <- factor(totalsmCCU$Condition,levels = unique(totalsmCCU$Condition))

totalsmCCUsuboptimal <- data.frame(Condition = c('Suboptimal'),
                                   pct = c(100-vlookup(2,mCCU,2,TRUE)),
                                   MoE = c(vlookup(2,mCCU,3,TRUE) * 1.96)) %>%
  mutate(Parameter='Cumulative Dissolved Metals') %>%
  select(Parameter,everything())


# VSCI DATA
VSCI <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'VSCIVCPMI')%>%
  select(Value, Estimate.P, StdError.P, NResp)
totalsVSCI <- data.frame(Condition = c('Suboptimal'),
                         pct = vlookup(60,VSCI,2,TRUE),
                         MoE = vlookup(60,VSCI,3,TRUE) *1.96) %>%
  mutate(Parameter='VSCI/VCPMI\n(Biomonitoring)') %>%
  select(Parameter,everything())
#VSCIsummary <- mutate(totalsVSCI,BelowStd=paste(formatC(pct,digits=2),"% ( +/- ",formatC(ErrorLCB95,digits=2),"% )", sep=""))%>%select(Parameter,BelowStd)
#colnames(VSCIsummary)<-c('Parameter','Below Standard ( 4 mg/L )')


# DNickel DATA
Nickel <- filter(dat, Subpopulation == 'IR2020' & Indicator == 'NICKEL') %>%
  select(Value, Estimate.P, StdError.P, NResp)
totalsNickel <- data.frame(Condition = c('Suboptimal'),
                           pct = vlookup(0.03,Nickel,2,TRUE),
                           MoE = vlookup(0.03,Nickel,3,TRUE)*1.96) %>%
  mutate(Parameter='Dissolved Nickel') %>%
  select(Parameter,everything())


# initial Graph
paramsummary <- rbind(totalsDO,totalshabsuboptimal,totalsLRBSsuboptimal,totalsTNsuboptimal,totalspH,totalsTPsuboptimal,
                      totalsTDSsuboptimal,totalsmCCUsuboptimal,totalsVSCI)%>%
  mutate(Standard = c('yes','no','no','no','yes','no','no','no','yes')) %>%
  arrange(desc(pct))
paramsummary$Parameter <- factor(paramsummary$Parameter)
paramsummary$Parameter <- reorder(paramsummary$Parameter, paramsummary$pct) # reorder for ggplot based on pct results

# Stressor extent
stressorext <- rbind(totalsLRBSsuboptimal,totalsTNsuboptimal,totalsTPsuboptimal,totalsTDSsuboptimal,totalsmCCUsuboptimal,totalshabsuboptimal) %>%
  arrange(desc(pct))
stressorext$Parameter <- factor(stressorext$Parameter) %>% # rearrange factor level order based on descending pct
  reorder(stressorext$pct)


## Play with stacked bar charts
together <- rbind(totalshab,totalsLRBS,totalsTN,totalsTP,totalsTDS,totalsmCCU)

togetherstacked <- together
togetherstacked$Parameter <- factor(togetherstacked$Parameter,
                                    levels=c('Streambed Sedimentation','Total Phosphorus',
                                             'Habitat Disturbance','Total Nitrogen',
                                             'Cumulative Dissolved Metals','Ionic Strength'))

togetherstackedflip <- together
togetherstackedflip$Parameter <- factor(togetherstackedflip$Parameter,
                                        levels=c('Cumulative \nDissolved Metals','Total Nitrogen',
                                                 'Ionic\nStrength','Total\nPhosphorus',
                                                 'Streambed\nSedimentation','Habitat Disturbance'))


# RELATIVE RISK PLOT
rr <- read.csv('processedData/relriskIR2020.csv') %>%
  mutate(Stressor=dplyr::recode(Stressor,"TotHabstatus"="Habitat Disturbance",
                                "TDSstatus"='Ionic Strength',
                                "TNstatus"='Total Nitrogen',
                                "MetalCCUstatus"='Cumulative Dissolved Metals',
                                "LRBSstatus"='Streambed Sedimentation',
                                "TPstatus"='Total Phosphorus')) %>%
  mutate(MoE = StdError.log *1.96)



rr$Stressor <- factor(rr$Stressor)
rr$Stressor <- reorder(rr$Stressor, rr$Estimate)

