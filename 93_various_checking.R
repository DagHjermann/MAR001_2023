

library(dplyr)
dat_a <- readRDS("Data/01_Combined_data_with_duplicates_a.rds")

#
# Why are there no data from east Baltic? ----
# Answer: no blue mussels or oyster species
#
dat_a %>%
  filter(Country %in% c("Finland", "Estonia", "Latvia")) %>% 
  distinct(Country, WoRMS_scientificName, SD_StationCode) %>% 
  xtabs(~WoRMS_scientificName + Country, .)

#
# Checking number of time series by country/parameter ----
#

# a. number of measurements, all species

d1a <- dat_a %>% 
  filter(PARAM %in% c("CD","CU","PB","HG","HCB","HCHG","DDEPP","BAP", "CB118"),
         Year >= 2010)
d1a %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# b. number of measurements, blue mussels and oysters

d1b <- d1a %>% 
  filter(grepl("Mytilus", WoRMS_scientificName) | WoRMS_scientificName %in% c("Crassostrea gigas", "Ostrea edulis"))

d1b %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# c. number of measurements designated to a station (time series)  

d2 <- d1b %>% 
  filter(!is.na(SD_StationCode))

d2 %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# d. number of stations    

d3 <- d2 %>%  
  distinct(Country, SD_StationCode, SD_StationName, Year, PARAM) %>% 
  count(Country, SD_StationCode, SD_StationName, PARAM, name = "n_year")

d3 %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# e. number of stations with at least 3 measurements after 2010       

d4 <- d3 %>%  
  filter(n_year >= 3)

d4 %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

