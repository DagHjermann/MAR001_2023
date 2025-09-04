

library(maps)
library(ggplot2)
library(ggeasy)
library(dplyr)
library(purrr)
library(mapdata)
library(ggeasy)
library(ggimage)    # geom_subview(), theme_transparent()

source("92_Various_for_metadata_delivery_functions.R")

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# For dataset "Spatial extent" ----
# :
# Upper-, lower-, left- and right coordinates for the bounding box

# Run first the start of script 01, until both EMODnet (df1a) and ICES data have been read
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

# EMODnet
range(df1a$Longitude)
range(df1a$Latitude)
# [1] -9.213 39.625
# [1] 35.820 65.818
# Upper = 65.818, Lower = 35.820 , Right = 39.625, Left = -9.213, 


# ICES
range(df2_orig$Longitude)
range(df2_orig$Latitude)
# [1] -41.75  54.45
# [1] 36.56667 80.16667

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# For map example ----
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
# . Map with points per station ----
#

# Station data
fn <- "Data/07_dat_status_trend_2025.rds"
dat_status_trend <- readRDS(fn)

very_simple_map <- map_data("world")


gg <- dat_status_trend %>%
  mutate(Concentrations = case_when(
    Status == 1 ~ "Low", 
    Status == 2 ~ "Moderate", 
    Status == 3 ~ "High")) %>% 
  ggplot(aes(Longitude, Latitude, color = Concentrations)) +
  annotation_map(very_simple_map, fill = "navajowhite2") +
  geom_point() +
  coord_map("lambert", parameters = c(2, 50)) +
  facet_wrap(vars(PARAM)) +
  easy_remove_axes() +
  theme(
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
    ) 
#  easy_text_size(which = "strip.text", size = 10)
# ?ggeasy::.all_theme_els

ggsave("Figures/92_example_map_all_pointsonly_2025.png", gg, width = 10, height = 9, dpi = 200)


#
# Region data ----
#

fn <- "Submitted/MAR001_status_trend_by_region.xlsx"
readxl::excel_sheets(fn)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
# Simple categorical trend + status (main data)
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

fn <- "Data/07_dat_status_trend_2025.rds"
dat_status_trend <- readRDS(fn)

# Delete one PCB118 outlier
dat_status_trend <- dat_status_trend %>%
  filter(!(id %in% "56.02833333_20.83333333" & PARAM %in% "CB118"))

dat_region_status <- dat_status_trend %>%
  count(PARAM, Region, Status)
  # mutate(
  # Status = case_when(
  #   Status %in% 1 ~ "Low",
  #   Status %in% 2 ~ "Moderate",
  #   Status %in% 3 ~ "High"),
  # Status = factor(Status, levels = c("Low", "Moderate", "High"))
  # )


#
# Results of trend metaanalysis  
#
dat_region_trend <- readRDS("Data/06_meta_trends_region.rds") %>%
  select(-data)

# check
table(meta_trends_region$PARAM)
# table(meta_trends_region$PARAM, meta_trends_region$Region)

# Fix region name
sel <- dat_region_trend$Region %in% "North-East Atlantic Ocean." 
dat_region_trend$Region[sel] <- "North-East Atlantic Ocean"

# check
table(dat_region_trend$Region)

# check
table(dat_region_trend$`Overall trend`)

# Change trend text
# - and select only the substances we want  
dat_region_trend <- dat_region_trend %>%
  filter(PARAM %in% unique(dat_status_trend$PARAM)) %>%
  mutate(
    `Overall trend` = case_when(
      `Overall trend` %in% "No sign. pattern" ~ "No change",
      `Overall trend` %in% "Downward" ~ "Improvement",
      `Overall trend` %in% "Upward" ~ "Detoriation"))
      
# check
table(dat_region_trend$`Overall trend`)

#
# . Save excel data for EEA ----
#

data_map_points <- dat_status_trend %>%
  select(PARAM, Region, Longitude, Latitude, Trend, Status) %>%
  mutate(
    `Point shape (trend)` = case_when(
      Trend %in% "Decrease" ~ "triangle pointing down",
      Trend %in% "Increase" ~ "triangle pointing up",
      TRUE ~ "round")
    ) %>%
  left_join(
    data.frame(Status = 1:3, Point_color = c("green", "orange", "red")), 
    relationship = "many-to-one") %>%
  rename(`Point color (status)` = Point_color)

data_map_pies <- dat_region_status_trend <- dat_region_status %>%
  left_join(
    data.frame(Status = 1:3, Pie_color = c("green", "orange", "red")), 
    relationship = "many-to-one")

data_map_pietext <- dat_region_trend %>% 
  select(PARAM, Region, `Overall trend`) 

writexl::write_xlsx(
  list(
    points = data_map_points,
    pies = data_map_pies,
    pie_text = data_map_pietext,
    info = data.frame(info = "One map per 'PARAM'")),
  "Figures/2025/Maps_by_parameter_plotdata_2025.xlsx"
)

#
# . Plot pies ----
#

# Based on seksjon 318/Elveovervakning:
# 'make_pie_from_colours' in 03_map_pie_functions.R


make_pie <- function(data, xvar, yvar, cols = c("1"="lightgreen", "2"="orange", "3"="red2")){
  data$x <- factor(data[[xvar]])
  data$y <- data[[yvar]]
  ggplot(data, aes(x=1, y, fill = x)) + 
    geom_bar(stat="identity", color = "black", size = 1) + coord_polar(theta="y") +
    scale_fill_manual(values = cols) +
    theme_void() + theme(legend.position="none") + theme_transparent()
}
test <- dat_region_status %>%
  filter(Region == "Baltic" & PARAM == "HG")
make_pie(test, xvar = "Status", yvar = "n")

#
# . Map with pie charts ----
# HG only
#

pies <- dat_region_status %>%
  filter(PARAM == "HG") %>%
  split(~Region) %>%
  purrr::map(make_pie, xvar = "Status", yvar = "n")
length(pies)

# Make data set for pie center coordinates
df_pies <- c(
  28.92, 56.97,     # Baltic
  8, 36,            # Mediterranean
  0.1, 66.1) %>%    # `North-East Atlantic Ocean` 
  matrix(ncol = 2, byrow = TRUE) %>%
  as.data.frame() %>%
  set_names(c("x", "y"))

# Add pies as a variable
# df_pies$pie <- pies[c(1,3,4)]   # skip numer 2, Black Sea
df_pies$pie <- pies   # no Black Sea here

# Sets size of pies
df_pies$width = 12  
df_pies$height = 12   

# Plot pies only
df_range <- data.frame(x = c(-9, 35), y = c(28, 73))
ggplot(data=df_range, aes(x, y)) + 
  geom_blank() +
  annotation_map(very_simple_map, fill = "sienna") +
  geom_subview(data=df_pies, aes(x=x, y=y, subview=pie, width=width, height=height))
# gg


# Plot points + pies for HG
gg <- dat_status_trend %>%
  filter(PARAM == "HG") %>%
  mutate(Concentrations = case_when(
    Status == 1 ~ "Low", 
    Status == 2 ~ "Moderate", 
    Status == 3 ~ "High")) %>% 
  ggplot(aes(Longitude, Latitude)) +
  annotation_map(very_simple_map, fill = "navajowhite2") +
  geom_point(aes(color = Concentrations)) +
  annotate("text", x = -Inf, y = Inf, label = "Mercury", hjust = -0.3, vjust = 1.3, size = 7) +  # 
  geom_subview(data=df_pies, aes(x=x, y=y, subview=pie, width=width, height=height)) +
  easy_remove_axes() +
  theme(panel.background = element_rect(fill = "azure"))
gg

ggsave("Figures/92_example_map_mercury_noarrows_2025.png", width = 6, height = 4.5, dpi = 200)

dat_region_trend %>%
  filter(PARAM == "HG")  
# PARAM Region                     Change_perc_10yr Change_perc_10yr_lo Change_perc_10yr_up `Overall trend` 
# 1 HG    Baltic                                -24.3              -47.2                 8.66 No sign. pattern
# 2 HG    North-East Atlantic Ocean.             10.8                1.07               21.4  Upward  


#
# . All pie charts + "arrow" maps  ----
#
# GAve up arrows - just made text instaed
#


# debugonce(plot_contaminant)


plot_contaminant("HG")
plot_contaminant("BAP")
# debugonce(plot_contaminant)
plot_contaminant("BDE47")

# Make all plots
plots <- unique(dat_region_status$PARAM) %>% map(plot_contaminant)

# Combine plots in a big plot (plot to screen)
# cowplot::plot_grid(plotlist = plots, nrow = 3)

# Combine plots in a big plot (plot to file)
comb_plot <- cowplot::plot_grid(plotlist = plots)
cowplot::save_plot(
  filename = "Figures/92_map_all_contaminants_25.png",
  plot = comb_plot, 
  nrow = 3, base_asp = 2.4)


