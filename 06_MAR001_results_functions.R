
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for relative class ----
#
##-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

#
# Create function which creates relative class, based on two thresholds  
# - by default, works on log scale  
# - 'lowest' (= relclass 1.0)  could be taken to be just below the lowest value observed 
#
create_class_function <- function(thresh1, thresh2,  
                                  lowest = 0, log = TRUE,
                                  from_class = 1, to_class = 100){
  
  # Make class data
  df_classes_all <- data.frame(
    conc = c(lowest, thresh1, thresh2),
    class = 1:3
  )
  
  if (to_class > 3){
    df_classes_all <- bind_rows(
      df_classes_all,
      data.frame(
        conc = seq(2,(to_class-2))*thresh2,
        class = seq(4,to_class))
    )
  }
  
  df_classes <- df_classes_all %>%
    filter(class >= from_class & class <= to_class)
  
  if (log){
    fun <- approxfun(log(df_classes$conc), df_classes$class)
    result <- function(x){ fun(log(x))  }
  } else {
    result <- approxfun(df_classes$conc, df_classes$class)
  }
  
  result

  }

# TEST
if (FALSE){
  f1 <- create_class_function(thresh1 = 20, thresh2 = 500, log = FALSE,
                             from_class = 1, to_class = 100)
  f1(c(8, 20, 30, 300, 500, 700))
  
  f2 <- create_class_function(thresh1 = 20, thresh2 = 500, log = TRUE,
                             from_class = 1, to_class = 100)
  f2(c(8, 20, 30, 300, 500, 700))
}


#
# Get Relclass when the two thresholds are on the same basis  
# - to be used on a 'preselected' dataset (only one PARAM and one Group)
# - returns a vector of values (same sequence as data given, can be used in mutate)
# - also works when only the lowest threshold is given
#
get_relclass_same_basis <- function(dataset, data_thresholds = df_thresholds){
  
  dataset_summ <- dataset %>%
    distinct(PARAM, Group)

  # Thresholds for test
  # - get parameter + Group from the 1st obs of the data
  dfth <- data_thresholds %>% 
    filter(PARAM %in% dataset_summ$PARAM & Group %in% dataset_summ$Group)
  
  basis <- dfth$Thresh1_basis

  if (basis == "W"){
    data_values <- dataset$perc90_value_ww
  } else if (basis == "D"){
    data_values <- dataset$perc90_value_dw
  } else if (basis == "L"){
    data_values <- dataset$perc90_value_fw
  }
  
  # If only threshold 1 is given, we set threshold 0 (lowest value),
  #   threshold 1 and threshold 2 at equal distance on log scaøe
  # Unless threshold 2 then is lower than the max value in the data   
  if (is.na(dfth$Thresh2)){
    dfth$Thresh2 <- exp(
      log(dfth$Thresh1) + 
        (log(dfth$Thresh1) - min(log(data_values), na.rm = TRUE))
    )
    # Check if threshold 2 then is lower than the max value in the data
    if (dfth$Thresh2 < max(data_values, na.rm = TRUE)){
      dfth$Thresh2 <- max(data_values, na.rm = TRUE)*1.1
    }
    # BAsis is the same
    dfth$Thresh2_basis <- dfth$Thresh1_basis
  }
    
    if (dfth$Thresh2_basis != basis){
    stop("The two bases are not the same")
  }
  
  f_interpol <- create_class_function(dfth$Thresh1, dfth$Thresh2,
                                      lowest = min(data_values)*0.99, 
                                      log = TRUE)
  f_interpol(data_values)

}


#
# Get Relclass when lower threshold is DW and higher is WW
#
get_relclass_dw_ww <- function(dataset, data_thresholds = df_thresholds){

  dataset_summ <- dataset %>%
    distinct(PARAM, Group)
  
  if (nrow(dataset_summ) > 1){
    stop("Dataset has >1 PARAM x Group")
  }
  
  # Add log-transformed data (for transforming thresholds)
  dataset <- dataset  %>% 
    mutate(log_perc90_value_dw = log(perc90_value_dw),
           log_perc90_value_ww = log(perc90_value_ww))
  
  # Thresholds for test
  # - get parameter + Group from the 1st obs of the data
  dfth <- data_thresholds %>% 
    filter(PARAM %in% dataset_summ$PARAM & Group %in% dataset_summ$Group)
  
  if (dfth$Thresh1_basis != "D")
    stop("Thresh1 basis must be D")
  if (dfth$Thresh2_basis != "W")
    stop("Thresh2 basis must be W")
  
  # 1. Make class data based on d.w. basis

  # 1a.   Thresh2 is transfomed to dw
  mod_ww_to_dw <- lm(log_perc90_value_dw ~ log_perc90_value_ww, data = dataset)
  Thresh2_dw_l <- predict(mod_ww_to_dw, 
                          newdata = data.frame(log_perc90_value_ww = log(dfth$Thresh2)))
  Thresh2_dw <- exp(Thresh2_dw_l)
  
  # 1b. Make approximation function
  f_dw <- create_class_function(dfth$Thresh1, Thresh2_dw,
                                lowest = min(dataset$perc90_value_dw)*0.99, 
                                log = TRUE)
  
  # 2. Make class data based on w.w. basis

  # 2a. Thresh1 is transfomed to ww
  mod_dw_to_ww <- lm(log_perc90_value_ww ~ log_perc90_value_dw, data = dataset)
  Thresh1_ww_l <- predict(mod_dw_to_ww, 
                          newdata = data.frame(log_perc90_value_dw = log(dfth$Thresh1)))
  Thresh1_ww <- exp(Thresh1_ww_l)
  
  # 2b. Make approximation function
  f_ww <- create_class_function(Thresh1_ww, dfth$Thresh2,
                                lowest = min(dataset$perc90_value_ww)*0.99, 
                                log = TRUE)
  
  dataset <- dataset %>%
    mutate(
      Status = factor(Status),
      Relclass_dw = f_dw(perc90_value_dw),
      Relclass_ww = f_ww(perc90_value_ww),
      Relclass = case_when(
        Status %in% 1 ~ Relclass_dw,
        Status %in% 2 ~ (Relclass_dw + Relclass_ww)/2,
        Status %in% 3 ~ Relclass_ww),
      Relclass = case_when(
        Status %in% 1 & Relclass >= 2 ~ 1.99,
        Status %in% 2 & Relclass < 2 ~ 2,
        Status %in% 2 & Relclass >= 3 ~ 2.99,
        TRUE ~ Relclass)
    )           
  
  dataset$Relclass
  
}

#
# Get Relclass when lower threshold is DW and higher is FW (fat weightm = lipid weight LW)
#
get_relclass_dw_fw <- function(dataset, data_thresholds = df_thresholds){
  
  dataset_summ <- dataset %>%
    distinct(PARAM, Group)
  
  if (nrow(dataset_summ) > 1){
    stop("Dataset has >1 PARAM x Group")
  }
  
  # Add log-transformed data (for transforming thresholds)
  dataset <- dataset  %>% 
    mutate(log_perc90_value_dw = log(perc90_value_dw),
           log_perc90_value_fw = log(perc90_value_fw))
  
  # Thresholds for test
  # - get parameter + Group from the 1st obs of the data
  dfth <- data_thresholds %>% 
    filter(PARAM %in% dataset_summ$PARAM & Group %in% dataset_summ$Group)
  
  if (dfth$Thresh1_basis != "D")
    stop("Thresh1 basis must be D")
  if (dfth$Thresh2_basis != "L")
    stop("Thresh2 basis must be L")
  
  # 1. Make class data based on d.w. basis
  
  # 1a.   Thresh2 is transfomed to dw
  mod_fw_to_dw <- lm(log_perc90_value_dw ~ log_perc90_value_fw, data = dataset)
  Thresh2_dw_l <- predict(mod_fw_to_dw, 
                          newdata = data.frame(log_perc90_value_fw = log(dfth$Thresh2)))
  Thresh2_dw <- exp(Thresh2_dw_l)
  
  # 1b. Make approximation function
  f_dw <- create_class_function(dfth$Thresh1, Thresh2_dw,
                                lowest = min(dataset$perc90_value_dw)*0.99, 
                                log = TRUE)
  
  # 2. Make class data based on w.w. basis
  
  # 2a. Thresh1 is transfomed to fw
  mod_dw_to_fw <- lm(log_perc90_value_fw ~ log_perc90_value_dw, data = dataset)
  Thresh1_fw_l <- predict(mod_dw_to_fw, 
                          newdata = data.frame(log_perc90_value_dw = log(dfth$Thresh1)))
  Thresh1_fw <- exp(Thresh1_fw_l)
  
  # 2b. Make approximation function
  f_fw <- create_class_function(Thresh1_fw, dfth$Thresh2,
                                lowest = min(dataset$perc90_value_fw)*0.99, 
                                log = TRUE)
  
  dataset <- dataset %>%
    mutate(
      Status = factor(Status),
      Relclass_dw = f_dw(perc90_value_dw),
      Relclass_fw = f_fw(perc90_value_fw),
      Relclass = case_when(
        Status %in% 1 ~ Relclass_dw,
        Status %in% 2 ~ (Relclass_dw + Relclass_fw)/2,
        Status %in% 3 ~ Relclass_fw),
      Relclass = case_when(
        Status %in% 1 & Relclass >= 2 ~ 1.99,
        Status %in% 2 & Relclass < 2 ~ 2,
        Status %in% 2 & Relclass >= 3 ~ 2.99,
        TRUE ~ Relclass)
    )           
  
  dataset$Relclass
  
}




#
# Function for adding relative class to data  
#
# The 'dataset' needs to contain the 'threshold_comb' variable
#

get_data_with_relclass <- function(param, group, threshold_comb, dataset) {
  if (threshold_comb %in% c("W-W","W-NA","D-NA")){
    result <- dataset %>%
      filter(PARAM %in% param & Group == group) %>%
      mutate(Status = factor(Status),
             Relclass = get_relclass_same_basis(.))
  } else if (threshold_comb %in% "D-W"){
    result <- dataset %>%
      filter(PARAM %in% param & Group == group) %>%
      mutate(Status = factor(Status),
             Relclass = get_relclass_dw_ww(.))
  } else if (threshold_comb %in% "D-L"){
    result <- dataset %>%
      filter(PARAM %in% param & Group == group) %>%
      mutate(Status = factor(Status),
             Relclass = get_relclass_dw_fw(.))
  }
  result
}



#
## Relative class by parameter - FUNCTION    
#
# As above (plot_relclass_region, just switching PARAM and region)
#

plot_relclass_param <- function(param, 
                                msfd_regions = FALSE,
                                exclude_norway_regions = TRUE){
  
  data_for_plot <- dat_status_trend_relclass %>%
    filter(PARAM %in% param) %>%
    as.data.frame()
  
  if (msfd_regions)
    data_for_plot <- data_for_plot %>%
      mutate(Region = MSFD_region)
  
  if (exclude_norway_regions)
    data_for_plot <- data_for_plot %>%
      filter(!Region %in% c("Norwegian Sea", "Barents Sea", "Icelandic Ocean"))
  
  above_plot_border <- data_for_plot %>%
    filter(Relclass > 10) %>%
    count(Region) %>%
    mutate(Relclass = 10.1, labeltext = paste("Over 10:\n", n, "stations"))
  
  gg <- ggplot(data_for_plot %>% filter(Relclass <= 10), 
               aes(x = Region, y = Relclass)) +
    geom_hline(aes(yintercept = 2), linetype = "dashed", color = "yellow3", size = rel(1))  +
    geom_hline(aes(yintercept = 3), linetype = "dashed", color = "red3", size = rel(1))  +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(fill = Status), alpha = 0.3, width = 0.25, shape = 21) +
    scale_fill_manual(values = c(`1` = "green3", `2` = "yellow2", `3` = "red")) +
    scale_y_continuous(breaks = 0:10, minor_breaks = NULL, limits = c(0,10.8)) + 
    geom_text(data = above_plot_border, aes(label = labeltext), 
              vjust = -0.3, size = rel(3), lineheight = 0.9) +
    labs(y = "Relative status class") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.title.y = element_text(margin = margin(r = 20))
    )
  
  if (msfd_regions)
    gg <- gg +
    theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
    labs(x = "MSFD region")
  
  gg
  
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for status tables ----
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

table_status_region <- function(region){
  dat_status_trend %>%
    filter(Region == region) %>%    # filter by region here
    mutate(PARAM = factor(
      PARAM, 
      levels = c("CD","CU","PB","HG","HCB","HCHG","DDEPP","BAP", "CB118"))
    ) %>%
    tabyl(Status, PARAM) %>%
    # adorn_totals("row") %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns() %>%
    kable(caption = region) %>%       # add region as caption  
    kable_classic(full_width = F, html_font = "Arial Narrow")
}

table_status_msfdregion <- function(region){
  dat_status_trend %>%
    rename(MSFD_region2 = MSFD_region) %>%
    # Add Overall trend
    mutate(
      MSFD_region = case_when(
        MSFD_region2 %in% "MWE" ~ "Western Mediterranean Sea",
        MSFD_region2 %in% "MAD" ~ "Adriatic Sea",
        MSFD_region2 %in% "MIC" ~ "Ionia Sea and central Med.",
        MSFD_region2 %in% "BLK" ~ "Black Sea",
        MSFD_region2 %in% "ABI" ~ "Bay of Biscay / Iberian coast",
        MSFD_region2 %in% "ACS" ~ "Celtic Seas",
        MSFD_region2 %in% "ANS" ~ "Greater North Sea",
        MSFD_region2 %in% "BAL" ~ "Baltic Sea",
        MSFD_region2 %in% "NOR" ~ "Norwegian Sea",
        MSFD_region2 %in% "BAR" ~ "Barents Sea",
        MSFD_region2 %in% "ICE" ~ "Icelandic Ocean")
    ) %>%
    mutate(MSFD_region = factor(MSFD_region, levels = msfd_regions)) %>%
    filter(MSFD_region %in% region) %>%    # filter by region here
    mutate(PARAM = factor(
      PARAM, 
      levels = c("CD","CU","PB","HG","HCB","HCHG","DDEPP","BAP", "CB118"))
    ) %>%
    tabyl(Status, PARAM) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns() %>%
    kable(caption = region) %>%       # add region as caption  
    kable_classic(full_width = F, html_font = "Arial Narrow")
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for plotting relative class ----
# 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


#
## Relative class by region - FUNCTION    
#

plot_relclass_region <- function(region){
  
  data_for_plot <- dat_status_trend_relclass %>%
    filter(Region %in% region)
  
  above_plot_border <- data_for_plot %>%
    filter(Relclass > 10) %>%
    count(PARAM) %>%
    mutate(Relclass = 10.1, labeltext = paste("> 10:\nn =", n))
  
  ggplot(data_for_plot %>% filter(Relclass <= 10), 
         aes(x = PARAM, y = Relclass)) +
    geom_boxplot() +
    geom_jitter(aes(fill = Status), alpha = 0.3, width = 0.25, shape = 21) +
    scale_fill_manual(values = c(`1` = "green3", `2` = "yellow2", `3` = "red")) +
    scale_y_continuous(breaks = 0:10, minor_breaks = NULL, limits = c(0,10.8)) + 
    geom_text(data = above_plot_border, aes(label = labeltext), 
              vjust = -0.3, size = rel(3.5), lineheight = 0.9) +
    theme_bw() 
  
}

plot_relclass_msfdregion <- function(msfdregion){
  
  data_for_plot <- dat_status_trend_relclass %>%
    filter(MSFD_region %in% msfdregion)
  
  above_plot_border <- data_for_plot %>%
    filter(Relclass > 10) %>%
    count(PARAM) %>%
    mutate(Relclass = 10.1, labeltext = paste("> 10:\nn =", n))
  
  ggplot(data_for_plot %>% filter(Relclass <= 10), 
         aes(x = PARAM, y = Relclass)) +
    geom_boxplot() +
    geom_jitter(aes(fill = Status), alpha = 0.3, width = 0.25, shape = 21) +
    scale_fill_manual(values = c(`1` = "green3", `2` = "yellow2", `3` = "red")) +
    scale_y_continuous(breaks = 0:10, minor_breaks = NULL, limits = c(0,10.8)) + 
    geom_text(data = above_plot_border, aes(label = labeltext), 
              vjust = -0.3, size = rel(3.5), lineheight = 0.9) +
    theme_bw() 
  
}

# table(dat_status_trend_relclass$Region) 
# table(dat_status_trend_relclass$MSFD_region) %>% names()   



#
## Relative class by parameter - FUNCTIONs    
#
# As above (plot_relclass_region, just switching PARAM and region)
#

plot_relclass_param <- function(param, 
                                msfd_regions = FALSE,
                                exclude_norway_regions = TRUE,
                                get_boxplot_only = FALSE,
                                get_points_only = FALSE,
                                max_relclass = 10,
                                y_axis_text = TRUE){
  
  data_for_plot <- dat_status_trend_relclass %>%
    filter(PARAM %in% param) %>%
    as.data.frame()
  
  if (msfd_regions)
    data_for_plot <- data_for_plot %>%
      mutate(Region = MSFD_region)
  
  if (exclude_norway_regions)
    data_for_plot <- data_for_plot %>%
      filter(!Region %in% c("Norwegian Sea", "Barents Sea", "Icelandic Ocean"))
  
  above_plot_border <- data_for_plot %>%
    filter(Relclass > max_relclass) %>%
    count(Region) %>%
    mutate(Relclass = max_relclass + 0.1, 
           labeltext = paste0("Over ", max_relclass-2, "x EQS:\n", n, " stations"))
  
  gg <- ggplot(data_for_plot %>% filter(Relclass <= max_relclass), 
               aes(x = Region, y = Relclass)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(fill = Status), alpha = 0.3, width = 0.25, shape = 21) +
    scale_fill_manual(values = c(`1` = "green3", `2` = "yellow2", `3` = "red")) +
    scale_y_continuous(
      breaks = 0:max_relclass, 
      minor_breaks = NULL, 
      limits = c(0, max_relclass + 0.8),
      labels = c("", "0", "BAC", "EQS", paste0(seq(2, max_relclass-2), "x EQS"))) + 
    geom_text(data = above_plot_border, aes(label = labeltext), 
              vjust = -0.3, size = rel(3), lineheight = 0.9) +
    labs(y = "Relative status class") +
    theme_bw() 
  
  if (y_axis_text){
    gg <- gg +
      scale_y_continuous(
        breaks = 1:max_relclass, 
        minor_breaks = NULL, 
        limits = c(1, max_relclass + 0.8),
        labels = c("0", "BAC", "EQS", paste0(seq(2, max_relclass-2), "x EQS")))
  } else {
    gg <- gg +
      scale_y_continuous(
        breaks = 0:max_relclass, 
        minor_breaks = NULL, 
        limits = c(0, max_relclass + 0.8))
  }
  
  if (msfd_regions){      
    gg <- gg +
      theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
      labs(x = "MSFD region")
  }
  
  # get_boxplot_only and get_points_only are only for extracting plot data
  if (get_boxplot_only){
    gg <- ggplot(data_for_plot %>% filter(Relclass <= max_relclass), 
                 aes(x = Region, y = Relclass)) +
      geom_boxplot(outlier.shape = NA)
  }
  if (get_points_only){
    gg <- ggplot(data_for_plot %>% filter(Relclass <= max_relclass), 
                 aes(x = Region, y = Relclass)) +
      geom_jitter(aes(fill = Status), alpha = 0.3, width = 0.25, shape = 21) +
      scale_fill_manual(values = c(`1` = "green3", `2` = "yellow2", `3` = "red"))
  }  
  gg +
    theme(panel.grid.minor.y = element_blank())
  
}



# table(dat_status_trend_relclass$Region) 
# table(dat_status_trend_relclass$MSFD_region) %>% names()   


#
# plot_relclass_param2   
# - as plot_relclass_param, but using geom_dots instead of geom_jitter
# - note new parameters binwidth and overflow
#   (see https://mjskay.github.io/ggdist/articles/dotsinterval.html)

plot_relclass_param2 <- function(param, 
                                msfd_regions = FALSE,
                                exclude_norway_regions = TRUE,
                                get_boxplot_only = FALSE,
                                get_points_only = FALSE,
                                binwidth = 0.15, overflow = "compress",
                                max_relclass = 10,
                                y_axis_text = TRUE){
  
  data_for_plot <- dat_status_trend_relclass %>%
    filter(PARAM %in% param) %>%
    mutate(
      Measurement = factor(
        ifelse(Below_LOQ, "Under LOQ", "Quantified"),
        levels = c("Quantified", "Under LOQ"))
      ) %>%
    as.data.frame()
  
  if (msfd_regions)
    data_for_plot <- data_for_plot %>%
      mutate(Region = MSFD_region)
  
  if (exclude_norway_regions)
    data_for_plot <- data_for_plot %>%
      filter(!Region %in% c("Norwegian Sea", "Barents Sea", "Icelandic Ocean"))
  
  above_plot_border <- data_for_plot %>%
    filter(Relclass > max_relclass) %>%
    count(Region) %>%
    mutate(Relclass = max_relclass + 0.1, 
           labeltext = paste0("Over ", max_relclass-2, "x EQS:\n", n, " stations"))
  
  gg1 <- ggplot(data_for_plot %>% filter(Relclass <= (max_relclass-2)), 
               aes(x = Region, y = Relclass))
  
  gg_add_box <- function(gg) {
    gg + geom_boxplot(outlier.shape = NA)
  }

  
  gg_add_points <- function(gg) {
    gg + 
      ggdist::geom_dots(aes(fill = Status, shape = Measurement), color = "black", linewidth = 0.35, side = "both", 
                        binwidth = binwidth, overflow = overflow) +
      # jitter version was: geom_jitter(aes(fill = Status), alpha = 0.3, width = 0.25, shape = 21)
      scale_shape_manual(values = c(21,25)) +
      scale_fill_manual(values = c(`1` = "green3", `2` = "yellow2", `3` = "red"))
  }
  
  # get_boxplot_only and get_points_only are mostly used for extracting plot data
  if (get_boxplot_only){
    gg2 <- gg1 |>
      gg_add_box()
  } else if (get_points_only){
    gg2 <- gg1 |>
      gg_add_points()
  } else {
    gg2 <- gg1 |>
      gg_add_box() |>
      gg_add_points()
  }
  gg <- gg2 +
    geom_text(data = above_plot_border, aes(label = labeltext), 
              vjust = -0.3, size = rel(3), lineheight = 0.9) +
    # guides(shape = "none") + 
    labs(y = "Relative status class") +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank()) 
  
  if (y_axis_text){
    gg <- gg +
      scale_y_continuous(
        breaks = 1:max_relclass, 
        minor_breaks = NULL, 
        limits = c(1, max_relclass + 0.8),
        labels = c("0", "BAC", "EQS", paste0(seq(2, max_relclass-2), "x EQS")))
  } else {
    gg <- gg +
      scale_y_continuous(
        breaks = 0:max_relclass, 
        minor_breaks = NULL, 
        limits = c(0, max_relclass + 0.8))
  }
  
  if (msfd_regions){      
    gg <- gg +
      theme(axis.text.x = element_text(angle = -35, hjust = 0)) +
      labs(x = "MSFD region")
  }
   
  gg
  
}


# table(dat_status_trend_relclass$Region) 
# table(dat_status_trend_relclass$MSFD_region) %>% names()   



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for trend tables ----
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


table_trends_region <- function(region){
  
  dat_status_trend %>%
    filter(Region == region) %>%
    mutate(PARAM = factor(
      PARAM, 
      levels = c("CD","CU","PB","HG","HCB","HCHG","DDEPP","BAP", "CB118"))
    ) %>%
    tabyl(Trend, PARAM) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns() %>%
    kable(caption = region) %>%       # add region as caption  
    kable_classic(full_width = F, html_font = "Arial Narrow") 
  
}

table_trends_msfdregion <- function(region){
  
  dat_status_trend %>%
    rename(MSFD_region2 = MSFD_region) %>%
    # Add Overall trend
    mutate(
      MSFD_region = case_when(
        MSFD_region2 %in% "MWE" ~ "Western Mediterranean Sea",
        MSFD_region2 %in% "MAD" ~ "Adriatic Sea",
        MSFD_region2 %in% "MIC" ~ "Ionia Sea and central Med.",
        MSFD_region2 %in% "BLK" ~ "Black Sea",
        MSFD_region2 %in% "ABI" ~ "Bay of Biscay / Iberian coast",
        MSFD_region2 %in% "ACS" ~ "Celtic Seas",
        MSFD_region2 %in% "ANS" ~ "Greater North Sea",
        MSFD_region2 %in% "BAL" ~ "Baltic Sea",
        MSFD_region2 %in% "NOR" ~ "Norwegian Sea",
        MSFD_region2 %in% "BAR" ~ "Barents Sea",
        MSFD_region2 %in% "ICE" ~ "Icelandic Ocean")
    ) %>%
    mutate(MSFD_region = factor(MSFD_region, levels = msfd_regions)) %>%
    filter(MSFD_region == region) %>%
    mutate(PARAM = factor(
      PARAM, 
      levels = c("CD","CU","PB","HG","HCB","HCHG","DDEPP","BAP", "CB118"))
    ) %>%
    tabyl(Trend, PARAM) %>%
    adorn_percentages("col") %>%
    adorn_pct_formatting(digits = 1) %>%
    adorn_ns() %>%
    kable(caption = region) %>%       # add region as caption  
    kable_classic(full_width = F, html_font = "Arial Narrow") 
  
}



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for trend meta analysis ----
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# Function - get estimates of slope   

# xtabs(~PARAM + Region, result_loq)

# Get estimates of slope, regardless of whether there is only on Region (no Baltic data) or not  
# - regardless  
get_slope_estimates_grouped <- function(df, group_variable){
  # Filter for useful data
  df <- df %>% 
    filter(!is.na(Slope))
  # Add varable Group (copy of 'group_variable' variable) 
  sel <- names(df) %in% group_variable
  if (sum(sel) != 1)
    stop(group_variable, " not found")
  df[["Group"]] <- df[[which(sel)]]
  # Find number of groups
  tab <- xtabs(~Group, df)
  if (length(tab) > 1){
    m <- lm(Slope ~ Group, data = df)     # include Group effect
  } else {
    m <- lm(Slope ~ 1, data = df)          # exclude Group effect
  }
  df_for_predict <- data.frame(Group = names(tab))
  result_list <- predict(m, df_for_predict, se.fit = TRUE)
  result <- df_for_predict %>%
    mutate(Estimate = result_list$fit,
           SE = result_list$se.fit,
           CI_lower = Estimate - 2*SE,
           CI_upper = Estimate + 2*SE
    )
  # Change name of Group variable to 'group_variable'
  names(result)[names(result) == "Group"] <- group_variable
  result
}

if (FALSE){
  # test
  # debugonce(get_slope_estimates_grouped)
  dat_series_trend %>% filter(PARAM == "BAP") %>% get_slope_estimates_grouped(group_variable = "Region")
  dat_series_trend %>% filter(PARAM == "BAP") %>% get_slope_estimates_grouped(group_variable = "MSFD_region")
}




#
# Function 'plot_slope_region'   
#

plot_slope_region <- function(param){
  meta_trends_region %>%
    filter(PARAM %in% param) %>%
    ggplot(
      aes(Region, Change_perc_10yr, color = `Overall trend`)) +
    geom_pointrange(
      aes(ymin = Change_perc_10yr_lo, ymax = Change_perc_10yr_up)) +
    scale_color_manual(values = col_scale) +
    geom_hline(yintercept = 0) +
    labs(title = param,
         x = "Region",
         y = "Percentage change / 10 year") +
    coord_flip() +
    theme_bw() +
    # Increase distance to axis labels
    theme(axis.title.x = element_text(margin = margin(t = 30)),
          axis.title.y = element_text(margin = margin(r = 30)))
}

#
# Function 'plot_slope_msfdregion'   
#

plot_slope_msfdregion <- function(param, exclude_norwegian_areas = TRUE){
  if (exclude_norwegian_areas){ 
    data_plot <- meta_trends_msfdregion %>%
      filter(!MSFD_region %in% c("Norwegian Sea", "Barents Sea", "Icelandic Ocean"))
  } else {
    data_plot <- meta_trends_msfdregion
  }
  data_plot %>%
    filter(PARAM %in% param) %>%
    ggplot(
      aes(MSFD_region, Change_perc_10yr, color = `Overall trend`)) +
    geom_pointrange(
      aes(ymin = Change_perc_10yr_lo, ymax = Change_perc_10yr_up)) +
    scale_color_manual(values = col_scale) +
    geom_hline(yintercept = 0) +
    labs(title = param,
         x = "Region",
         y = "Percentage change / 10 year") +
    coord_flip() +
    theme_bw() +
    # Increase distance to axis labels
    theme(axis.title.x = element_text(margin = margin(t = 20)),
          axis.title.y = element_text(margin = margin(r = 30)))
  
}

# plot_slope_msfdregion("CD")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for plotting level ----
# 
# Plots 90th percentile
#   perc90_value_ww, _dw or _fw
# With level as colours
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


plot_level_vs_species <- function(param, log = TRUE, threshold = 1, 
                                  data = dat_status_trend, df_thresh = df_thresholds){
  
  if (threshold == 1){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh1_basis) %>% unique()
  } else if (threshold == 2){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh2_basis) %>% unique()
  } else {
    stop("'threshold' must be either 1 or 2 (lower/upper")
  }
  
  if (length(basis) == 1 & basis == "D"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(WoRMS_scientificName, perc90_value_dw))
    
  } else if (length(basis) == 1 & basis == "W"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(WoRMS_scientificName, perc90_value_ww))
    
  } else if (length(basis) == 1 & basis == "L"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(WoRMS_scientificName, perc90_value_fw))
    
  }
  
  
  if (length(basis) == 1){
    
    if (threshold == 1){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh1) 
    } else if (threshold == 2){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh2)
    }

    gg <- gg +
      geom_jitter(aes(color = factor(Status)), width = 0.25, height = 0) +
      scale_color_manual(values = c(`1` = "blue", `2` = "green", `3` = "red")) +
      labs(
        title = param,
        subtitle = paste0("Threshold level ", threshold, 
                          "value = ", df_threshold$Threshold, ", basis = ", basis)) +
      geom_hline(
        data = df_threshold, aes(yintercept = Threshold), linetype = "dashed")

    if (log){
      gg <- gg + scale_y_log10() + coord_flip() 
    } else {
      gg <- gg + coord_flip() 
    }
    
    print(gg)
    
  } else if (length(basis) > 1){
    
    stop("More than one basis for this threshold 1 of this parameter! Levels =", paste(basis, collapse = ", "))
    
  } else if (length(basis) == 0){
    
    stop("This parameter not found, or basis not given for threshold 1")
    
  }
}



plot_level_vs_region <- function(param, log = TRUE, threshold = 1, 
                                 data = dat_status_trend, df_thresh = df_thresholds){
  
  if (threshold == 1){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh1_basis) %>% unique()
  } else if (threshold == 2){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh2_basis) %>% unique()
  } else {
    stop("'threshold' must be either 1 or 2 (lower/upper")
  }
  
  if (length(basis) == 1 & basis == "D"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Region, perc90_value_dw))
    
  } else if (length(basis) == 1 & basis == "W"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Region, perc90_value_ww))
    
  } else if (length(basis) == 1 & basis == "L"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Region, perc90_value_fw))
    
  }
  
  
  if (length(basis) == 1){
    
    if (threshold == 1){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh1) 
    } else if (threshold == 2){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh2)
    }
    
    gg <- gg +
      geom_jitter(aes(color = factor(Status)), width = 0.25, height = 0) +
      scale_color_manual(values = c(`1` = "blue", `2` = "green", `3` = "red")) +
      labs(
        title = param,
        subtitle = paste0("Threshold level ", threshold, 
                          "value = ", df_threshold$Threshold, ", basis = ", basis)) +
      geom_hline(
        data = df_threshold, aes(yintercept = Threshold), linetype = "dashed")
    
    if (log){
      gg <- gg + scale_y_log10() + coord_flip() 
    } else {
      gg <- gg + coord_flip() 
    }
    
    print(gg)
    
  } else if (length(basis) > 1){
    
    stop("More than one basis for this threshold 1 of this parameter! Levels =", paste(basis, collapse = ", "))
    
  } else if (length(basis) == 0){
    
    stop("This parameter not found, or basis not given for threshold 1")
    
  }
}



plot_level_vs_country <- function(param, log = TRUE, threshold = 1, 
                                  data = dat_status_trend, df_thresh = df_thresholds){
  
  if (threshold == 1){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh1_basis) %>% unique()
  } else if (threshold == 2){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh2_basis) %>% unique()
  } else {
    stop("'threshold' must be either 1 or 2 (lower/upper")
  }
  
  if (length(basis) == 1 & basis == "D"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Country, perc90_value_dw))
    
  } else if (length(basis) == 1 & basis == "W"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Country, perc90_value_ww))
    
  } else if (length(basis) == 1 & basis == "L"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Country, perc90_value_fw))
    
  }
  
  
  if (length(basis) == 1){
    
    if (threshold == 1){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh1) 
    } else if (threshold == 2){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh2)
    }
    
    gg <- gg +
      geom_jitter(aes(color = factor(Status)), width = 0.25, height = 0) +
      scale_color_manual(values = c(`1` = "blue", `2` = "green", `3` = "red")) +
      labs(
        title = param,
        subtitle = paste0("Threshold level ", threshold, 
                          "value = ", df_threshold$Threshold, ", basis = ", basis)) +
      geom_hline(
        data = df_threshold, aes(yintercept = Threshold), linetype = "dashed")

    if (log){
      gg <- gg + scale_y_log10() + coord_flip() 
    } else {
      gg <- gg + coord_flip() 
    }
    
    print(gg)
    
  } else if (length(basis) > 1){
    
    stop("More than one basis for this threshold 1 of this parameter! Levels =", paste(basis, collapse = ", "))
    
  } else if (length(basis) == 0){
    
    stop("This parameter not found, or basis not given for threshold 1")
    
  }
}



plot_level_vs_country_trend <- function(param, log = TRUE, threshold = 1, 
                                        data = dat_status_trend, df_thresh = df_thresholds){
  
  if (threshold == 1){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh1_basis) %>% unique()
  } else if (threshold == 2){
    basis <- df_thresholds %>% filter(PARAM == param) %>% pull(Thresh2_basis) %>% unique()
  } else {
    stop("'threshold' must be either 1 or 2 (lower/upper")
  }
  
  if (length(basis) == 1 & basis == "D"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      arrange(Status) %>%
      ggplot(aes(Country, perc90_value_dw))
    
  } else if (length(basis) == 1 & basis == "W"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Country, perc90_value_ww))
    
  } else if (length(basis) == 1 & basis == "L"){
    
    gg <- data %>%
      filter(PARAM == param & !is.na(Status)) %>%    #  & Country == "Norway"
      ggplot(aes(Country, perc90_value_fw))
    
  }
  
  
  if (length(basis) == 1){
    
    if (threshold == 1){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh1) 
    } else if (threshold == 2){
      df_threshold <- df_thresholds %>% 
        filter(PARAM == param) %>%
        mutate(Threshold = Thresh2)
    }
    
    gg <- gg +
      geom_jitter(aes(color = factor(Status)), width = 0.25, height = 0) +
      scale_color_manual(values = c(`1` = "blue", `2` = "green", `3` = "red")) +
      labs(
        title = param,
        subtitle = paste0("Threshold level ", threshold, 
                          "value = ", df_threshold$Threshold, ", basis = ", basis)) +
      geom_hline(
        data = df_threshold, aes(yintercept = Threshold), linetype = "dashed")

    if (log){
      gg <- gg + scale_y_log10() + coord_flip() 
    } else {
      gg <- gg + coord_flip() 
    }
    
    print(gg)
    
  } else if (length(basis) > 1){
    
    stop("More than one basis for this threshold 1 of this parameter! Levels =", paste(basis, collapse = ", "))
    
  } else if (length(basis) == 0){
    
    stop("This parameter not found, or basis not given for threshold 1")
    
  }
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#
# Functions for maps of levels + trends ----
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-


# Function for making map
# GLOBAL: 'dat_status_trend', colors_status 

# Note: the 'map_param' uses 3 other functions that also are used by 'map_param_for_combined_plot' below

# NOTE: 'labeltype' is set to either
# - 'low,medium,high' (or similar) - for using these labels for all substances  
# - 'custom'          (or similar) - for using separate labels per substance (Below EQS, Below MPC, etc.)  

map_param <- function(param_sel, labeltype = "custom", data_thresholds = df_thresholds){
  
  df_plot <- dat_status_trend %>%
    filter(PARAM == param_sel)    #  & Country == "Norway"
  
  lab <- get_labels(labeltype = labeltype, param = param_sel, data_thresholds = data_thresholds)
  
  # Set Status2 (used for labeling fill colors in map)
  df_plot$Status2 <- make_status2_column(df_plot, "Status", lab)
  
  # Set colors for the fill  
  lookup_color <- get_colors(lab)
  
  gg <- ggplot(df_plot, aes(Longitude, Latitude, fill = Status2, shape = Trend)) +
    annotation_map(my_map, fill = "grey70") +
    geom_point() +
    scale_fill_manual("Status", values = lookup_color) +
    scale_shape_manual(values = c(Decrease = 25, Increase = 24,
                                  `No trend` = 21, `Trend unknown` = 22))
  
  result <- gg +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 3)),
      shape = guide_legend(override.aes = list(size = 3))) +
    coord_map("ortho", orientation = c(52, 10, 0),
              xlim = c(-10.10816667, 30.36166667),  # fix limits of plot so all are equal
              ylim = c(36.541, 70.13833333)) +
    labs(title = param_sel) +
    theme_bw()  
  
  result
}

# get_labels
# Function used by 'map_param' + 'map_param_for_combined_plot' (below)

get_labels <- function(labeltype, param = NULL, data_thresholds){
  
  # Set labels to use for fill color legend depending on 'labeltype'
  if (grepl("low", labeltype, ignore.case = TRUE)){
    lab1 <- "Low"
    lab2 <- "Moderate"
    lab3 <- "High"
  } else if (grepl("custom", labeltype, ignore.case = TRUE)){
    # Get custom labels for map 
    df_thresholds_sel <- data_thresholds %>%
      filter(PARAM %in% param)
    lab1 <- df_thresholds_sel %>% pull(`Label below Th1`) %>% .[1]
    lab2 <- df_thresholds_sel %>% pull(`Label between Th1 and Th2`) %>% .[1]
    lab3 <- df_thresholds_sel %>% pull(`Label above Th2`) %>% .[1]
  }
  c(lab1, lab2, lab3)
}

# get_labels("low,medium,high")
# get_labels("custom")   # should fail
# get_labels("custom", "HG", df_thresholds)


# make_status2_column
# Function used by 'map_param' + 'map_param_for_combined_plot' (below)

# data = dataset (already filtered for parameter)
# status_column = name of status column, has values 1,2,3   
# labels = labels to use 
make_status2_column <- function(data, status_column, labels){
  
  result <- ""
  result[data[[status_column]] %in% 1] <- labels[1]
  result[data[[status_column]] %in% 2] <- labels[2]
  result[data[[status_column]] %in% 3] <- labels[3]
  
  # Set order (factor)
  if (is.na(labels[3])){
    result <- factor(result, levels = c(labels[1:2]))
  } else {
    result <- factor(result, levels = c(labels[1:3]))
  }
  
  result
  
}

# get_colors
# Function used by 'map_param' + 'map_param_for_combined_plot' (below)

# GLOBAL: colors_status

# df <- dat_status_trend %>% filter(PARAM == "CD")
# make_status2_column(df, "Status", c("Low", "Moderate", "High"))

get_colors <- function(labels){
  
  lookup_color <- c("","","")
  
  if (labels[1] == "Low"){
    lookup_color[1] <- colors_status["blue"]
  } else {
    stop("Please specify color for label no. 1 in this case")
  }
  if (labels[2] %in% c("Above EQS")){
    lookup_color[2] <- colors_status["orange"]
  } else if (labels[2] %in% c("Below EQS", "Below MPC", "Above background", "Moderate")){
    lookup_color[2] <- colors_status["yellow"]
  } else {
    stop("Please specify color for label no. 2 in this case")
  }
  if (is.na(labels[3])){
    lookup_color <- lookup_color[1:2]
  } else if (labels[3] %in% "Above EQS") {
    lookup_color[3] <- colors_status["orange"]
  } else if (labels[3] %in% c("Above MPC", "High")) {
    lookup_color[3] <- colors_status["red"]
  } else {
    stop("Please specify color for label no. 3 in this case")
  } 
  
  lookup_color
  
}

# get_colors(c("Low", "Moderate", "High"))
# get_labels("low,medium,high") %>% get_colors()
# get_labels("custom", "HG", df_thresholds) %>% get_colors()



# As map_param, but some changes for making combined plot:  
# - remove trend legend and put status legend inside plot
# - marked by # CHANGE

# As map_param, the 'labeltype' decides whether to use custom labels per parameter, 
#   or "low,medium,high"

map_param_for_combined_plot <- function(param_sel, 
                                        labeltype = "low,medium,high", 
                                        legends = "fill",
                                        legend_in_map = TRUE,
                                        data_thresholds = df_thresholds){
  
  df_plot <- dat_status_trend %>%
    filter(PARAM == param_sel)    #  & Country == "Norway"
  
  lab <- get_labels(labeltype = labeltype, param = param_sel, data_thresholds = data_thresholds)
  
  # Set Status2 (used for labeling fill colors in map)
  df_plot$Status2 <- make_status2_column(df_plot, "Status", lab)
  
  # Set colors for the fill  
  lookup_color <- get_colors(lab)
  
  gg <- ggplot(df_plot, aes(Longitude, Latitude, fill = Status2, shape = Trend)) +
    annotation_map(my_map, fill = "grey70") +
    geom_point() +
    scale_fill_manual("Status", values = lookup_color) +
    scale_shape_manual(values = c(Decrease = 25, Increase = 24,
                                  `No trend` = 21, `Trend unknown` = 22))
  
  if (legends == "shape"){
    # We want to plot a single shape legend below the graphs
    # So we need to make a plot with the shape legend only (no fill legend)
    result <- gg + 
      guides(fill = "none")
  } else if (legends == "fill"){
    # We want to plot a single shape legend below the graphs
    # So we need to make a plot with the shape legend only (no fill legend)
    result <- gg +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)),
        shape = "none")                          
  } else if (legends == "both"){
    result <- gg +
      guides(
        fill = guide_legend(override.aes = list(shape = 21)))                          
  } else if (legends == "none"){
    # We want to plot a single shape legend below the graphs
    # So we need to make a plot with the shape legend only (no fill legend)
    result <- gg +
      guides(
        fill = "none",
        shape = "none")                          
  }
  
  # The ordinary plot is plotted without the shape legend (with fill legend only)
  result <- result +
    coord_map("ortho", orientation = c(52, 10, 0),
              xlim = c(-10.10816667, 30.36166667),  # CHANGE: fix limits of plot so all are equal
              ylim = c(36.541, 70.13833333)) +
    labs(title = param_sel) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  if (legend_in_map){
    result <- result +
      theme(
        legend.position = c(0.01, .98),        # CHANGE (put status legend inside plot)  
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.margin = margin(1,1,1,1),
        legend.key.height = unit(0.8,"line"),
        legend.text = element_text(size = 9),
      )
  }
  
  
  result
  
}

