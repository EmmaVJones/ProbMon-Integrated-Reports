library(raster)
library(rgdal)
library(maptools)
library(rgeos)
library(reshape)
library(reshape2)
library(plyr)
library(tidyverse)
library(sf)

landcover <- raster("E:/evjones/GIS/ProbMonGIS/GISdata/nlcd2011va.TIF")
#VAecoregionsL3 <- rgdal::readOGR('originalData', 'VA_level3ecoregion')# for basic site maps
VAecoregionsL3 <- st_read('originalData/VA_level3ecoregion.shp') %>%# for basic site maps
  #filter(US_L3NAME != "Middle Atlantic Coastal Plain") %>%
  filter(US_L3NAME %in% c('Central Appalachians','Blue Ridge')) %>%
  group_by(US_L3NAME) %>%
  summarize() %>%
  #st_combine() %>%
  as_Spatial()
wshdList <- as.character(VAecoregionsL3$US_L3NAME)
wshdPolys <- VAecoregionsL3

#wshdPolys <- readOGR('E:/evjones/GIS/ProbMonGIS/DelineatedWatersheds/YearlyAnalyses','2015_final_EVJ')
#wshdSites <- readOGR('E:/evjones/GIS/ProbMonGIS/DelineatedWatersheds/YearlyAnalyses','2015_finalsites_EVJ')
#wshdPolys@data$StationID <- sub("\r\n" ,"",wshdPolys@data$SiteID)
#wshdList <- as.character(wshdPolys$SiteID)
#siteList <- as.character(wshdSites$SiteID)


df <- data.frame(matrix(NA, ncol = 17, nrow = 30))
names(df) <- c('StationID','VALUE_11','VALUE_12','VALUE_21','VALUE_22','VALUE_23','VALUE_24'
               ,'VALUE_31','VALUE_41','VALUE_42','VALUE_43','VALUE_52','VALUE_71','VALUE_81'
               ,'VALUE_82','VALUE_90','VALUE_95')
template <- data.frame(VALUE_11=0,VALUE_12=0,VALUE_21=0,VALUE_22=0, VALUE_23=0,VALUE_24=0
                       ,VALUE_31=0,VALUE_41=0,VALUE_42=0,VALUE_43=0,VALUE_52=0,VALUE_71=0
                       ,VALUE_81=0,VALUE_82=0,VALUE_90=0,VALUE_95=0)



# Loop through polygons in shapefile to get landuse information
for(i in 1:length(wshdList)){ 
  print(i)
  e = raster::extract(landcover, wshdPolys[i,],small=T, na.rm=F)
  et = lapply(e,table)
  if(!length(et[[1]])==0){ 
    t <- melt(et)
    t.cast <- cast(t, L1 ~ Var1, sum)
    names(t.cast)[1] <- "StationID"
    print(t.cast[1,1] <- wshdList[i])
    colnames(t.cast)[2:length(names(t.cast))] <- paste("VALUE_",colnames(t.cast)
                                                       [2:length(names(t.cast))], sep = "")
  }
  if(length(et[[1]])==0){ 
    t.cast <- data.frame(wshdPolys[i,]@data[1])
    print(head(t.cast))
    rownames(t.cast) <- 1
    t.cast$VALUE_11= 0
  }
  zeros = template[is.na(match(names(template), names(t.cast)))]
  results = data.frame(t.cast, zeros)
  results <- results[,order(names(results))]
  df[i,] <- results
}

# Combine all subwatersheds by StationID
df$VALUE_12 <- NULL
df1 <- melt(df, 'StationID')
df2 <- aggregate(. ~ StationID + variable, data=df1,sum)
landuse <- mutate(df2, acre=900*value*0.0002471053814672,sqMile=acre*0.0015625
                  ,hectare=sqMile*258.9988110336)
landuse2 <- aggregate(cbind(sqMile) ~ StationID, data=landuse, FUN='sum')
names(landuse2) <- c('StationID','totalArea_sqMile')
landuse3 <- merge(landuse, landuse2, by='StationID')
landuse4 <- ddply(landuse3,c('StationID','variable'),mutate
                  ,Water_sqMile=sum(sqMile[(variable=='VALUE_11')])
                  ,PWater=(Water_sqMile/totalArea_sqMile)*100
                  ,N_sqMile=sum(sqMile[(variable=='VALUE_31')|(variable=='VALUE_41')
                                       |(variable=='VALUE_42')|(variable=='VALUE_43')
                                       |(variable=='VALUE_51')|(variable=='VALUE_52')
                                       |(variable=='VALUE_71')|(variable=='VALUE_72')
                                       |(variable=='VALUE_73')|(variable=='VALUE_74')
                                       |(variable=='VALUE_91')|(variable=='VALUE_95')])
                  ,N_INDEX=(N_sqMile/totalArea_sqMile)*100
                  ,Forest_sqMile=sum(sqMile[(variable=='VALUE_41')|(variable=='VALUE_42')
                                            |(variable=='VALUE_43')])
                  ,PFOR=(sum(Forest_sqMile)/totalArea_sqMile)*100
                  ,Wetland_sqMile=sum(sqMile[(variable=='VALUE_90')|(variable=='VALUE_95')])
                  ,PWETL=(Wetland_sqMile/totalArea_sqMile)*100
                  ,Shrub_sqMile=sum(sqMile[(variable=='VALUE_51')|(variable=='VALUE_52')])
                  ,PSHRB=(Shrub_sqMile/totalArea_sqMile)*100
                  ,Ngrasslands_sqMile=sum(sqMile[(variable=='VALUE_71')|(variable=='VALUE_72')
                                                 |(variable=='VALUE_73')|(variable=='VALUE_74')])
                  ,PNG=(Ngrasslands_sqMile/totalArea_sqMile)*100
                  ,Barren_sqMile=sum(sqMile[(variable=='VALUE_31')])
                  ,PBAR=(Barren_sqMile/totalArea_sqMile)*100
                  ,TotBarren_sqMile=sum(sqMile[(variable=='VALUE_21')|(variable=='VALUE_31')])
                  ,PTotBAR=(TotBarren_sqMile/totalArea_sqMile)*100
                  ,U_sqMile=sum(sqMile[(variable=='VALUE_21')|(variable=='VALUE_22')
                                       |(variable=='VALUE_23')|(variable=='VALUE_24')
                                       |(variable=='VALUE_81')|(variable=='VALUE_82')])
                  ,U_INDEX=(U_sqMile/totalArea_sqMile)*100
                  ,Urban_sqMile=sum(sqMile[(variable=='VALUE_21')|(variable=='VALUE_22')
                                           |(variable=='VALUE_23')|(variable=='VALUE_24')])
                  ,PURB=(Urban_sqMile/totalArea_sqMile)*100
                  ,MBAR_sqMile=sum(sqMile[(variable=='VALUE_21')])
                  ,PMBAR=(MBAR_sqMile/totalArea_sqMile)*100
                  ,AGT_sqMile=sum(sqMile[(variable=='VALUE_81')|(variable=='VALUE_82')])
                  ,PAGT=(AGT_sqMile/totalArea_sqMile)*100
                  ,AGP_sqMile=sum(sqMile[(variable=='VALUE_81')])
                  ,PAGP=(AGP_sqMile/totalArea_sqMile)*100
                  ,AGC_sqMile=sum(sqMile[(variable=='VALUE_82')])
                  ,PAGC=(AGC_sqMile/totalArea_sqMile)*100)
landuselong <- melt(landuse4,id.vars=c('StationID','totalArea_sqMile')
                    ,measure.vars=c('PWater','N_INDEX','PFOR','PWETL','PSHRB','PNG','PBAR','PTotBAR','U_INDEX'
                                    ,'PURB','PMBAR','PAGT','PAGP','PAGC'))
landusewide<- dcast(landuselong,StationID+totalArea_sqMile~variable,value.var='value',sum)
saveRDS(landusewide, 'preliminaryLanduse.RDS')
