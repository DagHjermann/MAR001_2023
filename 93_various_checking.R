

# load package(s) ----
library(dplyr)
library(ggplot2)
library(ggeasy)
library(forcats)
library(glue)

# load data (from script 01)
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

# . a. number of measurements, all species ----

d1a <- dat_a %>% 
  filter(PARAM %in% c("CD","CU","PB","HG","HCB","HCHG","DDEPP","BAP", "CB118"),
         Year >= 2010)
d1a %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# . b. number of measurements, blue mussels and oysters ----

d1b <- d1a %>% 
  filter(grepl("Mytilus", WoRMS_scientificName) | WoRMS_scientificName %in% c("Crassostrea gigas", "Ostrea edulis"))

d1b %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# . c. number of measurements designated to a station (time series)  ----

d2 <- d1b %>% 
  filter(!is.na(SD_StationCode))

d2 %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# . d. number of stations    ----

d3 <- d2 %>%  
  distinct(Country, SD_StationCode, SD_StationName, Year, PARAM) %>% 
  count(Country, SD_StationCode, SD_StationName, PARAM, name = "n_year")

d3 %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)

# . e. number of stations with at least 3 measurements after 2010 ----      

d4 <- d3 %>%  
  filter(n_year >= 3)

d4 %>% 
  count(Country, PARAM) %>% 
  pivot_wider(names_from = PARAM, values_from = n)


#
# Check meta-analysis of trends (Cd in particular) ----
#


# . functions and data (from script 06) ----  

# Functions
source("06_MAR001_results_functions.R")

# Trend + status

fn <- "Data/07_dat_status_trend_2024.rds"
dat_status_trend <- readRDS(fn)

# Delete one PCB118 outlier
dat_status_trend <- dat_status_trend %>%
  filter(!(id %in% "56.02833333_20.83333333" & PARAM %in% "CB118"))
# nrow(dat_status_trend)

cat("File modified at: \n")
file.info(fn)$mtime

# Trend including slope, will be used for meta-analysis
dat_series_trend <- readRDS("Data/05_dat_series_trend_2024.rds")

table_trends_region("North-East Atlantic Ocean")
table_trends_region("Baltic")  

# . check one parameter / area  ----

param <- "CD"
region <- "Baltic"
# region <- "North-East Atlantic Ocean"

# Trend classification only
# Note: only 3 or more years 

# .. tables ----

# dat1 = Trend classification only
dat1 <- dat_status_trend %>%
  filter(Region == region & PARAM == param) 

xtabs(~addNA(Trend), dat1)
xtabs(~n_years + addNA(Trend), dat1)

# dat2 = Trend including slope
# Note: onlyalso including 1-2 year time series   
dat2 <- dat_series_trend %>%
  mutate(Region = sub(".", "", Region, fixed = TRUE)) %>% #xtabs(~Region, .)
  filter(Region == region & PARAM == param) 

xtabs(~n_years, dat2)
xtabs(~(n_years <= 2), dat2)
xtabs(~(n_years <= 2) + trend, dat2)

# .. plots ----

# histogram of slope estimates
ggplot(dat2 %>% 
         filter(n_years >= 3 & trend != "Trend unknown"), 
       aes(x = Slope)) +
  geom_histogram(fill = "blue") +
  geom_vline(xintercept = 0) +
  # easy_rotate_x_labels(angle = -45) +
  labs(title = glue("{param}, {region}"),
       x = "Trend estimate")

# point and line plot
ggplot(dat2 %>% 
         filter(n_years >= 3 & trend != "Trend unknown"), 
       aes(x = reorder(id, Slope), y = Slope, ymin = CI_lower, ymax = CI_upper)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(color = "blue") +
  easy_rotate_x_labels(angle = -45) +
  labs(title = glue("{param}, {region}"),
       x = "Station ID",
       y = "Trend estimate")

ggplot(dat2 %>% 
         mutate(
           Slope2 = case_when(
             n_years < 3 ~ 0,
             trend %in% "Trend unknown" ~ 0,
             TRUE ~ Slope),
           pointcolour = case_when(
             n_years < 3 ~ "(b) 1 or 2 years data",
             trend %in% "Trend unknown" ~ "(c) Estimation failed",
             TRUE ~ "(a) Trend estimated"),
           order = case_when(
             n_years < 3 ~ 2000,
             trend %in% "Trend unknown" ~ 3000,
             TRUE ~ Slope)),
       aes(x = reorder(id, order), y = Slope2, ymin = CI_lower, ymax = CI_upper, colour = pointcolour)) +
  geom_hline(yintercept = 0) +
  geom_pointrange() +
  scale_color_manual(values = c("blue", "grey", "pink2")) +
  easy_rotate_x_labels(angle = -45) +
  labs(title = glue("{param}, {region}"),
       x = "Station ID",
       y = "Trend estimate")  

