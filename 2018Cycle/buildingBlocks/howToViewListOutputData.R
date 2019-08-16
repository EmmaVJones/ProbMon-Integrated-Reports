# remember this one?

# VLOOKUP (Excel function hack) by Julin Maloof, edited for list application by Emma Jones
vlookupEVJ <- function(table, #the table where you want to look for it; will look in first column
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
    
  } else {
    output <- table[match(ref,table[,1]),column]
    output[!ref %in% table[,1]] <- NA #not needed?
  }
  dim(output) <- dim(ref)
  output
}


# Step1: bring in VSCI data

listsOnListsOnLists <- readRDS('processedData/VSCI2018.RDS')

# Step2: How to view VSCI data by subpopulations

# all stuff we ran
names(listsOnListsOnLists) #75 ways and counting!

# Just virginia
View(listsOnListsOnLists[['estimateVirginia']]) # look at what is inside of the estimateVirginia list 
# (best if done in more recent verion of RStudio)

# CDF data for VA
View(listsOnListsOnLists[['estimateVirginia']]$CDF)
#PCT data for VA
View(listsOnListsOnLists[['estimateVirginia']]$Pct)


# same format for any subpopulation
View(listsOnListsOnLists[['estimateRoanoke']]) 
View(listsOnListsOnLists[['estimateRoanoke']]$CDF)
View(listsOnListsOnLists[['estimateRoanoke']]$Pct)




# Step 3: pull out VSCI at 60 across all subpopulations
library(purrr)
listsOnListsOnLists %>% map('CDF')  # extract all the CDF datasets from each list item, cool but we want just where each
# subpopulation has a VSCI of 60


### the purrr package is AMAZING
VSCI60 <- listsOnListsOnLists %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P')) %>% # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
  map( vlookupEVJ, 60  , 2, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
  tibble::enframe(value = "CDFestimateAt60") # change the results from that list into a tibble (kinda like dataframe)


### Compare to Trend adjusted data with NEW method
listsOnListsOnLists_trend <- readRDS('processedData/VSCI2018.RDS')

VSCI60_trend <- listsOnListsOnLists_trend %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P')) %>% # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
  map( vlookupEVJ, 60  , 2, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
  tibble::enframe(value = "CDFestimateAt60") # change the results from that list into a tibble (kinda like dataframe)

VSCI60df <- data.frame(name=VSCI60$name, oldpopPCT_VSCIat60=(as.matrix(unlist(VSCI60$CDFestimateAt60))))
VSCI60_trenddf <- data.frame(name=VSCI60_trend$name, newpopPCT_VSCIat60=(as.matrix(unlist(VSCI60_trend$CDFestimateAt60))))

overall <- full_join(data.frame(VSCI60df),data.frame(VSCI60_trenddf),by='name') %>%
  mutate(absdiff=abs(oldpopPCT_VSCIat60-newpopPCT_VSCIat60))


write.csv(overall,'processedData/comparisonVSCIoldvsnewTrendMethod.csv', row.names = F)






### LRBS
listsOnListsOnLists <- readRDS('processedData/LRBS.RDS')

LRBS0 <- listsOnListsOnLists[3:77] %>% # start with our big 'ol list of lists
  map('CDF') %>% #extract just the 'CDF' list from each list item (subpopulation)
  map(`[`, c('Value','Estimate.P')) %>% # kinda like dplyr::select here, extract the Value and Estimate.P columns from each list item (CDF)
  map( vlookupEVJ, 0  , 2, TRUE) %>% # use the vlookup function on each of the tables to find the Estimate.P where value = 60
  tibble::enframe(value = "CDFestimateAt0") # change the results from that list into a tibble (kinda like dataframe)

