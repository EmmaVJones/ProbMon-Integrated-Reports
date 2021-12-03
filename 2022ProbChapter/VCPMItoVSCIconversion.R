# This script converts average VCPMI scores to a 1-100 scale VSCI equivalent for IR reporting purposes
# Method originally developed by Jason Hill and converted to R by Emma Jones

library(tidyverse)
probRatio <- tribble(~`0 to 100 convesion`,	~`VSCI Score`,	~`ratio`,	~`category`,
                     0,	10,	10,	'severe stress',
                     2,	10,	10,		'severe stress',
                     4,	10,	10,	'severe stress',
                     6,	10,	10,	'severe stress',
                     8,	10,	1, 'severe stress',
                     10,	16.5,	1.375,		'severe stress',
                     12,	23,	1.642857143,		'severe stress',
                     14,	29.5,	1.84375,		'severe stress',
                     16,	36,	2,		'severe stress',
                     18,	42,	2.1,		'severe stress',
                     20,	43.8,	1.990909091,		'stress',
                     22,	45.6,	1.9,		'stress',
                     24,	47.4,	1.823076923,		'stress',
                     26,	49.2,	1.757142857,		'stress',
                     28,	51,	1.7,		'stress',
                     30,	52.8,	1.65,		'stress',
                     32,	54.6,	1.605882353,		'stress',
                     34,	56.4,	1.566666667,		'stress',
                     36,	58.2,	1.531578947,		'stress',
                     38,	60,	1.5,		'good',
                     40,	61.1,	1.454761905,		'good',
                     42,	62.2,	1.413636364,		'good',
                     44,	63.3,	1.376086957,		'good',
                     46,	64.4,	1.341666667,		'good',
                     48,	65.5,	1.31,		'good',
                     50,	66.6,	1.280769231,		'good',
                     52,	67.7,	1.253703704,		'good',
                     54,	68.8,	1.228571429,		'good',
                     56,	69.9,	1.205172414,		'good',
                     58,	72,	1.2,		'excellent',
                     60,	73.125,	1.179435484,		'excellent',
                     62,	74.25,	1.16015625,		'excellent',
                     64,	75.375,	1.142045455,		'excellent',
                     66,	76.5,	1.125,		'excellent',
                     68,	77.625,	1.108928571,		'excellent',
                     70,	78.75,	1.09375,		'excellent',
                     72,	79.875,	1.079391892,		'excellent',
                     74,	81,	1.065789474,		'excellent',
                     76,	82.125,	1.052884615,		'excellent',
                     78,	83.25,	1.040625,		'excellent',
                     80,	84.375,	1.028963415,		'excellent',
                     82,	85.5,	1.017857143,		'excellent',
                     84,	86.625,	1.007267442,		'excellent',
                     86,	87.75,	0.997159091,		'excellent',
                     88,	90,	1,	'excellent',
                     90,	90,	90,	'excellent',
                     92,	90,	90,	'excellent',
                     94,	90,	90,		'excellent',
                     96,	90,	90,'excellent',
                     98,	90,	90,		'excellent' )

VCPMItoVSCIconversion <- function(averageVCPMI, probRatio){
  averageVCPMI * vlookup(averageVCPMI, probRatio, 3, range = T, larger = F)
}

# How to use function
#madeUpVCPMI <- c(19.79, 29.92, 33.09, 16.99, 26.63, 45.15, 26.13, 37.66, 62.31, 29.03, 60.37, 77.33)
#VCPMItoVSCIconversion(madeUpVCPMI, probRatio)





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
