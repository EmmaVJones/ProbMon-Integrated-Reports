listOfResults <- function(popstatus.est, finalData, parameterName){
  list(
    popstatus.est = popstatus.est,
    dataAnalyzed = finalData,
    # All Virginia
    estimateVirginia = subpopEstimate(finalData,parameterName,NA, 'Virginia',NA, specialWeight=FALSE),
    # By Basin
    estimateJames = subpopEstimate(finalData,parameterName, 'Basin', 'James','James Basin', specialWeight=FALSE) ,
    estimateBigSandy = subpopEstimate(finalData,parameterName, 'SubBasin', 'Big Sandy',NA, specialWeight=FALSE) ,
    # IR window
    estimateIR2008 = subpopEstimate(finalData,parameterName, 'IR2008', '2008','IR2008', specialWeight='finalweight_IR2008'),
    estimateIR2010 = subpopEstimate(finalData,parameterName, 'IR2010', '2010','IR2010', specialWeight='finalweight_IR2010'),
    estimateIR2012= subpopEstimate(finalData,parameterName, 'IR2012', '2012','IR2012', specialWeight='finalweight_IR2012'),
    estimateIR2014 = subpopEstimate(finalData,parameterName, 'IR2014', '2014','IR2014', specialWeight='finalweight_IR2014'),
    estimateIR2016 = subpopEstimate(finalData,parameterName, 'IR2016', '2016','IR2016', specialWeight='finalweight_IR2016'),
    estimateIR2018 = subpopEstimate(finalData,parameterName, 'IR2018', '2018','IR2018', specialWeight='finalweight_IR2018')
  )
}


subpopulationCategory <- 'IR2008'; subpopulation<- '2008';  altName<- 'IR2008'; specialWeight='finalweight_IR2008')
subpopulation<- 1





surveyDataJBS <- filter(surveyData, SubBasin %in% c("Big Sandy", "James")) # 30 Big Sandy
surveyDataJBS_2008 <- filter(surveyDataJBS, Year <= 2006) # 10 Big Sandy
surveyDataJBS_2010 <- filter(surveyDataJBS, Year <= 2008)# 12 Big Sandy
surveyDataJBS_2012 <- filter(surveyDataJBS, Year <= 2010)# 17 Big Sandy
surveyDataJBS_2014 <- filter(surveyDataJBS, Year <= 2012)# 22 Big Sandy
surveyDataJBS_2016 <- filter(surveyDataJBS, Year <= 2014)# 27 Big Sandy


nrow(filter(surveyDataJBS_2016, SubBasin == 'Big Sandy'))


VSCIJBS_2008 <- allCDFresults(designStatus, surveyDataJBS_2008,'VSCIVCPMI','finalweight_IR2008')
VSCIJBS_2010 <- allCDFresults(designStatus, surveyDataJBS_2010,'VSCIVCPMI','finalweight_IR2010')
VSCIJBS_2012 <- allCDFresults(designStatus, surveyDataJBS_2012,'VSCIVCPMI','finalweight_IR2012')
VSCIJBS_2014 <- allCDFresults(designStatus, surveyDataJBS_2014,'VSCIVCPMI','finalweight_IR2014')
VSCIJBS_2016 <- allCDFresults(designStatus, surveyDataJBS_2016,'VSCIVCPMI','finalweight_IR2016')
VSCIJBS_2018 <- allCDFresults(designStatus, surveyDataJBS,'VSCIVCPMI','finalweight_IR2018')


# find CDF results at 60 for desired subpopulations
VSCI60JBS_2008_CDF <- VSCIJBS_2008[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P'))  # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
VSCI60JBS_2008 <- suppressWarnings(map(VSCI60JBS_2008_CDF, vlookupEVJ2, 60, 2:3, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
                                     tibble(
                                       Subpopulation = names(.), 
                                       IRwindow = 2008,
                                       Estimate.P = map_dbl(., 'Estimate.P'), 
                                       StdError.P =  map_dbl(., 'StdError.P')) %>%
                                     select(-.))

VSCI60JBS_2010_CDF <- VSCIJBS_2010[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P')) # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
VSCI60JBS_2010 <- suppressWarnings(map(VSCI60JBS_2010_CDF, vlookupEVJ2, 60, 2:3, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
                                     tibble(
    Subpopulation = names(.), 
    IRwindow = 2010,
    Estimate.P = map_dbl(., 'Estimate.P'), 
    StdError.P =  map_dbl(., 'StdError.P')) %>%
  select(-.))

VSCI60JBS_2012_CDF <- VSCIJBS_2012[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P'))  # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
VSCI60JBS_2012 <- suppressWarnings(map(VSCI60JBS_2012_CDF, vlookupEVJ2, 60, 2:3, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
                                     tibble(
                                       Subpopulation = names(.), 
                                       IRwindow = 2012,
                                       Estimate.P = map_dbl(., 'Estimate.P'), 
                                       StdError.P =  map_dbl(., 'StdError.P')) %>%
                                     select(-.))

VSCI60JBS_2014_CDF <- VSCIJBS_2014[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P'))  # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
VSCI60JBS_2014 <- suppressWarnings(map(VSCI60JBS_2014_CDF, vlookupEVJ2, 60, 2:3, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
                                     tibble(
                                       Subpopulation = names(.), 
                                       IRwindow = 2014,
                                       Estimate.P = map_dbl(., 'Estimate.P'), 
                                       StdError.P =  map_dbl(., 'StdError.P')) %>%
                                     select(-.))


VSCI60JBS_2016_CDF <- VSCIJBS_2016[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P')) # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
VSCI60JBS_2016 <- suppressWarnings(map(VSCI60JBS_2016_CDF, vlookupEVJ2, 60, 2:3, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
                                     tibble(
                                       Subpopulation = names(.), 
                                       IRwindow = 2016,
                                       Estimate.P = map_dbl(., 'Estimate.P'), 
                                       StdError.P =  map_dbl(., 'StdError.P')) %>%
                                     select(-.))

VSCI60JBS_2018_CDF <- VSCIJBS_2018[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P'))  # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
VSCI60JBS_2018 <- suppressWarnings(map(VSCI60JBS_2018_CDF, vlookupEVJ2, 60, 2:3, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
                                     tibble(
                                       Subpopulation = names(.), 
                                       IRwindow = 2018,
                                       Estimate.P = map_dbl(., 'Estimate.P'), 
                                       StdError.P =  map_dbl(., 'StdError.P')) %>%
                                     select(-.))

VSCIoverWindows <- bind_rows(VSCI60JBS_2008,VSCI60JBS_2010,VSCI60JBS_2012,VSCI60JBS_2014,VSCI60JBS_2016,VSCI60JBS_2018)










table <- VSCIJBS_2014[c("estimateJames","estimateBigSandy")] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P','StdError.P'))
table <- table[[2]]
ref <- 60
column <- 2
range=TRUE

str(vlookupEVJ2(table, 60 , 2:3, TRUE ))

vlookupEVJ2 <- function(table, #the table where you want to look for it; will look in first column
                       ref, #the value or values that you want to look for
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
    names(output) <- names(table[column])
    
  } else {
    output <- table[match(ref,table[,1]),column]
    output[!ref %in% table[,1]] <- NA #not needed?
    names(output) <- names(table[column])
  }
  #dim(output) <- dim(ref)
  output
}











#### CDF Curves Over Time (Scatterplot)

```{r CDF plotly}

VSCI_CDFdf <- rbind(VSCI60JBS_2008_CDF[1] %>% map(`[`, c('Value','Estimate.P','StdError.P')), VSCI60JBS_2008_CDF[2] %>% map_df(1) %>% rename('estimate'))

VSCI60JBS_2008_CDF[1] %>% map(`[`, c('Value','Estimate.P','StdError.P')) %>% mutate(Subpopulation = names(1))#mutate(Subpopulation = strsplit(names(.),'estimate')[[1]][2])

VSCI60JBS_2008_CDF[1] %>% map_df(1) %>% mutate(Subpopulation = strsplit(names(.),'estimate')[[1]][2]) %>% rename(1 = ')
```



```{r}
# add margin of error to plotly plots for each subpopulation efficiently
addMoE <- function(p, dataset, subpopulation){
  add_ribbons(p, data = filter(dataset, Subpopulation== subpopulation),
              x = ~Value, ymin = ~ymin, ymax = ~ymax, line = list(color = 'rgba(7, 164, 181, 0.05)'),
              fillcolor = 'rgba(7, 164, 181, 0.2)', name = paste(subpopulation," Margin of Error",sep=""), visible = 'legendonly')
}

plot_ly(cdfData()) %>%
        add_trace(data = cdfData(), x = ~Value, y = ~Estimate.P, mode = 'line', color = ~Subpopulation, 
                  hoverinfo = 'text', text = ~paste(sep = '<br>',
                                                    paste("Subpopulation: ", Subpopulation),
                                                    paste(unique(Indicator),":", Value), # add units!!!!!!!!!! see benthic stressor tool data manipulation??
                                                    paste('Percentile:',format(Estimate.P,digits=2), '+/-', format(MoE, digits = 2) ))) %>%
        addMoE(cdfData(),'Coast Bioregion') %>%
        addMoE(cdfData(),'Piedmont Bioregion') %>%
        addMoE(cdfData(),'Mountain Bioregion') %>%
        layout(#showlegend=FALSE,
          yaxis=list(title="Percentile"),
          xaxis=list(title="Subpopulation"))
```
