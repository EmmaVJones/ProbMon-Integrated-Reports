# This script walks users through making a boatable site spatial file for 2020 IR Prob Chapter
library(tidyverse)
library(sf)
library(readxl)

# Bring in 2017, 2018 sites
x2017 <- read_excel('C:/HardDriveBackup/ProbMon/2018/EmmaGIS20172018.xlsx', sheet = 'GIScrossWalk2017') %>%
  filter(str_detect(`Sample Code`,'boatable')) %>%
  mutate(StationID = DEQSITEID, Year = 2017) %>%
  st_as_sf(coords = c("LONG_DD", "LAT_DD"),  # make spatial layer using these columns
                      remove = T, # don't remove these lat/lon cols from df
                      crs = 4326) %>%
  dplyr::select(StationID, Year)
x2018 <- read_excel('C:/HardDriveBackup/ProbMon/2018/EmmaGIS20172018.xlsx', sheet = 'GIScrossWalk2018') %>%
  filter(str_detect(`Sample Code`,'boatable')) %>%
  mutate(StationID = DEQSITEID, Year = 2018) %>%
  st_as_sf(coords = c("LONG_DD", "LAT_DD"),  # make spatial layer using these columns
           remove = T, # don't remove these lat/lon cols from df
           crs = 4326) %>%
  dplyr::select(StationID, Year)

# bring in previous boatable sites information (2008-2016)
boatableSites <- st_read("originalData/boatableSites2018.shp") %>%
  st_transform(4326) %>%
  mutate(StationID = StatnID) %>%
  dplyr::select(StationID, Year)


# Combine
boatableSites <- rbind(boatableSites, x2017, x2018)
st_write(boatableSites, 'processedData/boatableSites2008_2018.shp')
