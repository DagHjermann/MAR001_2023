
#
# From "K:\Avdeling\214-Oseanografi\DHJ\R\R selfmade functions.R"  
#

# 23. distance in kilometers between two long/lat positions (from "fossil" package)
earth.dist <- function (long1, lat1, long2, lat2) 
{
  rad <- pi/180
  a1 <- as.numeric(lat1) * rad
  a2 <- as.numeric(long1) * rad
  b1 <- as.numeric(lat2) * rad
  b2 <- as.numeric(long2) * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}

# Distance from position number i in data frame 'data' to all other positions  
# 'var_id' is the name of the variable that uniquely identifies each row in the data
earth.dist.fromone <- function(i, data, var_id, var_long = "Longitude", var_lat = "Latitude"){
  dist <- earth.dist(data[[var_long]][i], data[[var_lat]][i], data[[var_long]], data[[var_lat]])
  data.frame(id1 = data[[var_id]][i],
             id2 = data[[var_id]],
             dist = dist)
}

# Distance between all pairs of rows in data frame 'data'    
# 'var_id' is the name of the variable that uniquely identifies each row in the data
earth.dist.pairwise <- function(data, var_id, var_long = var_long, var_lat = var_lat){
  if (length(unique(data[[var_id]])) != nrow(data)){
    stop("Variable ", sQuote(var_id), " doesn't uniquely identify each row in the data")
  }
  result <- map_dfr(
    1:nrow(data),
    ~earth.dist.fromone(i = ., data = data, var_id = var_id)
  )
  # Distance to itself should be set to NA    
  result[["dist"]][result[["id1"]] == result[["id2"]]] <- NA
  result
}

if (FALSE){
  
  test_data <- structure(list(
    Country = c("Finland_ICES", "Finland_ICES", "Finland_ICES"), 
    Latitude = c(58.75, 58.76, 59.533), 
    Longitude = c(20.50, 20.51, 22.001), Dataset = c("ICES", "ICES", "ICES"), 
    Station = c("A", "B", "C")), 
    row.names = c(NA, -3L), 
    class = c("tbl_df", "tbl", "data.frame"))
  
  test1 <- earth.dist.fromone(1, test_data, var_id = "Station")
  # debugonce(earth.dist.pairwise)
  test2 <- earth.dist.pairwise(test_data, var_id = "Station")
  test2 <- earth.dist.pairwise(test_data, var_id = "Country")

}

#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o
#
# Checking duplicates ----  
#
# Main workhorse function: 'check_pair_matrix_param_year'  
# - make changes here to add/remove variables  
#
# Main function to use: 'check_pair'  
#
# - which calls 'check_pair_species', which calls 'check_pair_matrix', etc. back to
#   'check_pair_matrix_param_year'
#
#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o#o

#
# Which basis to use?
#
# METALS assessment  
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-
#        |  2 DW  | 1 DW  | 0 DW |    
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-
#  2 WW  |   DW   |  DW   |  WW
#  1 WW  |   DW   |  DW   |  WW
#  0 WW  |   DW   |  DW   |  -
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-
#
# METALS comparison  
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-
#        |  2 DW  | 1 DW  | 0 DW |    
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-
#  2 WW  |   DW   |  WW   |  WW
#  1 WW  |   DW   |  -    |  -
#  0 WW  |   DW   |  -    |  -
# -o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-


check_pair_matrix_param_year <- function(id1 = "59.684_23.683", id2 = "59.68333333_23.63333333", 
                                         var_id = "LatLong", 
                                         dataset, 
                                         species = "Clupea harengus",
                                         matrix = "LI", param = "CD", year = 2000){
  
  #
  # Pick data (for both 'stations')
  # - calculate dry-weight and wet-weight value  
  #
  # Value_dw = Value_ww/Drywt = (Value_fw*Fatwt)/Drywt 
  # Value_ww = Value_dw*Drywt = Value_fw*Fatwt 
  #
  # We use species-specific fat weight in any case  
  #
  df_both <- dataset[dataset[[var_id]] %in% c(id1,id2),] %>%
    filter(WoRMS_scientificName == species & MATRX == matrix & PARAM == param & Year == year) %>%
    mutate(
      # Calculate dry and wet-weight basis (but DRYWT often lacking in one of the pair)
      Value_dw = case_when(
        BASIS %in% "W" & DRYWT_MUNIT %in% "%" ~ Value/(DRYWT*0.01),
        BASIS %in% "L" & DRYWT_MUNIT %in% "%" ~ (Value*FATWT_species*0.01)/(DRYWT*0.01),
        BASIS %in% "D" ~ Value),
      Value_ww = case_when(
        BASIS %in% "W" ~ Value,
        BASIS %in% "L" ~ Value*(FATWT_species*0.01),
        BASIS %in% "D" & DRYWT_MUNIT %in% "%" ~ Value*(DRYWT*0.01))
    )

  # Pick data for each station  
  df1_orig <- df_both[df_both[[var_id]] == id1,]
  df2_orig <- df_both[df_both[[var_id]] == id2,]
  
  # Sampling time   
  if (nrow(df1_orig) > 0){
    t1_all <- ymd_hms(paste0(as.character(df1_orig$Samplingtime), "000000"))
    t1 <- min(t1_all) + diff(range(t1_all))/2   # midpoint of sampling times
  } else {
    t1 <- as.POSIXct(NA)
  }
  if (nrow(df2_orig) > 0){
    t2_all <- ymd_hms(paste0(as.character(df2_orig$Samplingtime), "000000"))
    t2 <- min(t2_all) + diff(range(t2_all))/2   # midpoint of sampling times
  } else {
    t2 <- as.POSIXct(NA)
  }
  days_diff <- (t1-t2)/ddays(1)  
  
  # Original values, basis and sample size  
  val1 <- median(df1_orig$Value, na.rm = TRUE)
  val2 <- median(df2_orig$Value, na.rm = TRUE)
  n1 <- sum(!is.na(df1_orig$Value))
  n2 <- sum(!is.na(df2_orig$Value))
  if (n1>0){
    val1_range <- range(df1_orig$Value, na.rm = TRUE) %>% paste(collapse = "-")
  } else {
    val1_range <- as.character(NA)
  }
  if (n2>0){
    val2_range <- range(df2_orig$Value, na.rm = TRUE) %>% paste(collapse = "-")
  } else {
    val2_range <- as.character(NA)
  }
  bas1 <- paste(unique(df1_orig$BASIS), collapse = ",")
  bas2 <- paste(unique(df2_orig$BASIS), collapse = ",")

  # LOQ
  under_loq1 = sum(!is.na(df1_orig$QFLAG))
  under_loq2 = sum(!is.na(df2_orig$QFLAG))
  
  # same_length
  
  # Count number of dry-weight and wet-weight values
  n1_dw <- sum(!is.na(df1_orig$Value_dw))
  n2_dw <- sum(!is.na(df2_orig$Value_dw))
  n1_ww <- sum(!is.na(df1_orig$Value_ww))
  n2_ww <- sum(!is.na(df2_orig$Value_ww))  
  
  metals <- c("ZN", "CU", "CD", "HG", "CR", "NI", "PB", "AS")
  metal <- param %in% metals
  
  # For comparison
  if (metal){
    # Prefer dry weight, but more important to have two values  
    if (n1_dw > 0 & n2_dw > 0){
      compare_basis = "D"
    } else if (n1_ww > 0 & n2_ww > 0){
      compare_basis = "W" 
    } else {
      compare_basis = "D"
    }
  } else {
    # Prefer wet weight, but more important to have two values  
    if (n1_ww > 0 & n2_ww > 0){
      compare_basis = "W"
    } else if (n1_dw > 0 & n2_dw > 0){
      compare_basis = "D" 
    } else {
      compare_basis = "W"
    }
  }
  
  # For assessment
  if (metal){
    # Prefer dry weight, even when there are two wet-weight values  
    if (n1_dw > 0 | n2_dw > 0){
      assess_basis = "D"
    } else {
      assess_basis = "W"
    }
  } else {
    # Prefer wet weight, even when there are two dry-weight values  
    if (n1_ww > 0 | n2_ww > 0){
      assess_basis = "W"
    } else {
      assess_basis = "D"
    }
  }

  # If both data sets have at least one dry-weight value, we compare Value_dw
  if (compare_basis == "D"){
    # Get rid of duplicates within each series (if existing)
    df1 <- df1_orig %>% distinct(Value_dw) %>% arrange(Value_dw)
    df2 <- df2_orig %>% distinct(Value_dw) %>% arrange(Value_dw)
    # Make median of each data set 
    val1_c <- median(df1$Value_dw, na.rm = TRUE)
    val2_c <- median(df2$Value_dw, na.rm = TRUE)
    dw1 <- median(df1_orig$DRYWT, na.rm = TRUE)
    dw2 <- median(df2_orig$DRYWT, na.rm = TRUE)
    # Comparison
    pairwise_median_value <- median(c(val1_c, val2_c), na.rm = TRUE)
    pairwise_dev <- val1_c - val2_c
    pairwise_dev_perc <- pairwise_dev/pairwise_median_value*100
    dupl_1 <- nrow(df1) < nrow(df1_orig)   # TRUE if duplicates within dataset 1
    dupl_2 <- nrow(df2) < nrow(df2_orig)   # TRUE if duplicates within dataset 2
    diff_n <- nrow(df1) - nrow(df2)
    pairwise_basis <- "D"
  # If both data sets have at least one wet-weight value, we compare Value_ww
  } else if (compare_basis == "W"){
    # Get rid of duplicatges within each series (if existing)
    df1 <- df1_orig %>% distinct(Value_ww) %>% arrange(Value_ww)
    df2 <- df2_orig %>% distinct(Value_ww) %>% arrange(Value_ww)
    # Make median of each data set 
    val1_c <- median(df1$Value_ww, na.rm = TRUE)
    val2_c <- median(df2$Value_ww, na.rm = TRUE)
    dw1 <- median(df1_orig$DRYWT, na.rm = TRUE)
    dw2 <- median(df2_orig$DRYWT, na.rm = TRUE)
    # Comparison
    pairwise_median_value <- median(c(val1_c, val2_c), na.rm = TRUE)
    pairwise_dev <- val1_c - val2_c
    pairwise_dev_perc <- pairwise_dev/pairwise_median_value*100
    dupl_1 <- nrow(df1) < nrow(df1_orig)
    dupl_2 <- nrow(df2) < nrow(df2_orig)
    diff_n <- nrow(df1) - nrow(df2)
    pairwise_basis <- "W"
  # Else we give up (fat-weight basis not implemented)
  } else {
    if (nrow(df1_orig) > 0){
      pairwise_median_value <- median(df1_orig$Value_orig)
      pairwise_basis <- df1_orig$BASIS[1]
    } else {
      pairwise_median_value <- median(df2_orig$Value_orig)
      pairwise_basis <- df1_orig$BASIS[1]    }
    val1_c <- NA
    val2_c <- NA
    dw1 <- median(df1_orig$DRYWT, na.rm = TRUE)
    dw2 <- median(df2_orig$DRYWT, na.rm = TRUE)
    pairwise_dev <- NA
    pairwise_dev_perc <- NA
    dupl_1 <- NA
    dupl_2 <- NA
    diff_n <- NA
  }

  # Assessment
  if (assess_basis %in% "D"){
    df1 <- df1_orig %>% distinct(Value_dw) %>% arrange(Value_dw)
    df2 <- df2_orig %>% distinct(Value_dw) %>% arrange(Value_dw)
    val1_as <- median(df1$Value_dw, na.rm = TRUE)
    val2_as <- median(df2$Value_dw, na.rm = TRUE)
  } else if (assess_basis %in% "W"){
    df1 <- df1_orig %>% distinct(Value_ww) %>% arrange(Value_ww)
    df2 <- df2_orig %>% distinct(Value_ww) %>% arrange(Value_ww)
    val1_as <- median(df1$Value_ww, na.rm = TRUE)
    val2_as <- median(df2$Value_ww, na.rm = TRUE)
  } else {
    val1_as <- NA
    val2_as <- NA
  }
  assess_median_value <- median(c(val1_as, val2_as), na.rm = TRUE)
  
  tibble(
    id1 = id1,
    id2 = id2,
    species = species,
    matrix = matrix,
    param = param,
    metal = metal,
    year = year,
    t1 = t1,
    t2 = t2,
    days_diff = days_diff,
    val1 = val1,
    val2 = val2,
    ran1 = val1_range,
    ran2 = val2_range,
    bas1 = bas1,
    bas2 = bas2,
    u_loq1 = under_loq1,
    u_loq2 = under_loq2,
    n1 = n1,
    n2 = n2,
    val1_c = val1_c,
    val2_c = val2_c,
    dw1 = dw1,
    dw2 = dw2,
    common_basis = pairwise_basis,
    val1_as = val1_as,
    val2_as = val2_as,
    assess_basis = assess_basis,
    overall_median = pairwise_median_value,
    dev_min = min(abs(pairwise_dev)),
    dev_max = max(abs(pairwise_dev)),
    dev_mean = mean(abs(pairwise_dev)),
    devperc_min = min(abs(pairwise_dev_perc)) %>% round(1),
    devperc_max = max(abs(pairwise_dev_perc)) %>% round(1),
    devperc_mean = mean(abs(pairwise_dev_perc)) %>% round(1),
    dupl_1 = dupl_1,
    dupl_2 = dupl_2,
    diff_n = diff_n
  )
  
}

if (FALSE){
  # debugonce(check_pair_matrix_param)
  check_pair_matrix_param_year(dataset = dat_a, year = 2000)
  check_pair_matrix_param_year(dataset = dat_a, year = 2001)
  check_pair_matrix_param_year(dataset = dat_a, param = "PB", year = 2000)
}

#
# - function for checking duplicates for common years, for a given matrix and parameter  
#

check_pair_matrix_param <- function(id1 = "59.684_23.683", id2 = "59.68333333_23.63333333", var_id = "LatLong", 
                                    dataset, 
                                    species = "Clupea harengus",
                                    matrix = "LI", param = "CD"){

  year1 <- dataset[dataset[[var_id]] == id1,] %>%
    filter(WoRMS_scientificName == species & MATRX == matrix & PARAM %in% param) %>%
    count(Year, name = "n1")
  year2 <- dataset[dataset[[var_id]] == id2,] %>%
    filter(WoRMS_scientificName == species & MATRX == matrix & PARAM %in% param) %>%
    count(Year, name = "n2")
  
  df_join <- full_join(year1, year2, by = "Year")
  
  year_common <- df_join %>% filter(n1 > 0 & n2 > 0) %>% pull(Year)
  year_all <- df_join %>% pull(Year) %>% sort()
  
  df_check <- map_dfr(
    year_all, 
    ~check_pair_matrix_param_year(
      id1 = id1, id2 = id2, var_id = var_id, dataset = dataset, 
      species = species, matrix = matrix, param = param, year = .)
  )
  
  df_check %>%
    mutate(
      Year_common = length(year_common),
      Year_all = length(year_all),
      Year_1_only = df_join %>% filter(n1 > 0 & is.na(n2)) %>% nrow(),
      Year_2_only = df_join %>% filter(is.na(n1) & n2 > 0) %>% nrow()
    )

}

if (FALSE){
  debugonce(check_pair_matrix_param)
  #check_pair_matrix_param(dataset = dat_a)
  check_pair_matrix_param(dataset = dat_a, param = "PB")
  # check_pair_matrix_param_year(dataset = dat_a, param = "PB")
}


#
# - function for checking duplicate for common years and parameters, for a given matrix  
#

check_pair_matrix <- function(id1 = "59.684_23.683", id2 = "59.68333333_23.63333333", var_id = "LatLong", 
                              dataset,
                              species = "Clupea harengus",
                              matrix = "LI"){
  
  param1 <- dataset[dataset[[var_id]] == id1,] %>%
    filter(WoRMS_scientificName == species & MATRX == matrix) %>%
    count(PARAM, name = "n1")
  param2 <- dataset[dataset[[var_id]] == id2,] %>%
    filter(WoRMS_scientificName == species & MATRX == matrix) %>%
    count(PARAM, name = "n2")
  
  df_join <- full_join(param1, param2, by = "PARAM")
  
  param_common <- df_join %>% filter(n1 > 0 & n2 > 0) %>% pull(PARAM)
  param_all <- df_join %>% pull(PARAM) %>% sort()
  
  df_check <- map_dfr(
    param_all, 
    ~check_pair_matrix_param(
      id1 = id1, id2 = id2, var_id = var_id, dataset = dataset, 
      species = species, matrix = matrix, param = .)
  )
  
  df_check %>%
    mutate(
      PARAM_common = length(param_common),
      PARAM_all = length(param_all),
      PARAM_1_only = df_join %>% filter(n1 > 0 & is.na(n2)) %>% nrow(),
      PARAM_2_only = df_join %>% filter(is.na(n1) & n2 > 0) %>% nrow()
    )
  
  
}


if (FALSE){
  check_pair_matrix(dataset = dat_a)
  debugonce(check_pair_matrix)
  check_pair_matrix(dataset = dat_a, matrix = "LI")
  
}


#
# - function for checking duplicate, for all matrices (and parameters)    
#
  
check_pair_species <- function(id1 = "59.684_23.683", id2 = "59.68333333_23.63333333", var_id = "LatLong", 
                               species = "Clupea harengus",
                               dataset){
  
  matrx1 <- dataset[dataset[[var_id]] == id1,] %>%
    filter(WoRMS_scientificName == species) %>%
    count(MATRX, name = "n1")
  matrx2 <- dataset[dataset[[var_id]] == id2,] %>%
    filter(WoRMS_scientificName == species) %>%
    count(MATRX, name = "n2")
  
  df_join <- full_join(matrx1, matrx2, by = "MATRX")
  
  matrx_common <- df_join %>% filter(n1 > 0 & n2 > 0) %>% mutate() %>% pull(MATRX)
  matrx_all <- df_join %>% pull(MATRX) %>% sort()
  
  df_check <- map_dfr(
    matrx_all, 
    ~check_pair_matrix(
      id1 = id1, id2 = id2, var_id = var_id, dataset = dataset, 
      species = species, matrix =  .)
  )
}
  # df_check %>% mutate(MATRX_common = 2)
 # df_check
  
  # df_check %>%
  #   mutate(
  #     MATRX_common = length(matrx_common),
  #     MATRX_all = length(matrx_all),
  #     MATRX_1_only = df_join %>% filter(n1 > 0 & is.na(n2)) %>% nrow(),
  #     MATRX_2_only = df_join %>% filter(is.na(n1) & n2 > 0) %>% nrow()
  #   )
  #}


#
# - function for checking duplicate, for all matrices (and parameters)    
#

check_pair <- function(id1 = "59.684_23.683", id2 = "59.68333333_23.63333333", var_id = "LatLong", 
                               dataset){
  
  spec1 <- dataset[dataset[[var_id]] == id1,] %>%
    count(WoRMS_scientificName, name = "n1")
  spec2 <- dataset[dataset[[var_id]] == id2,] %>%
    count(WoRMS_scientificName, name = "n2")
  
  df_join <- full_join(spec1, spec2, by = "WoRMS_scientificName")
  
  spec_common <- df_join %>% filter(n1 > 0 & n2 > 0) %>% mutate() %>% pull(WoRMS_scientificName)
  spec_all <- df_join %>% pull(WoRMS_scientificName) %>% sort()
  
  df_check <- map_dfr(
    spec_all, 
    ~check_pair_species(
      id1 = id1, id2 = id2, var_id = var_id, dataset = dataset, 
      species = .)
  )
}
if (FALSE){
  # debugonce(check_pair)
  # debugonce(check_pair_matrix)
  check_pair(dataset = dat_a)
}


#
# As check_pair, but we assume that the data we test have a variable 'Station_dataset' formed by
#   adding "_EMODnet" or "_ICES" to the end of 'id'  
# Example, 'dat_test' for the Danish data (see DK4 and DK6)
#
# Note: dataset 1 (val1_c etc.) will always be EMODnet, dataet 2 = ICES 
#
check_pair_by_stationname <- function(station_name, data, logfile = NULL){
  
  stat1 <- paste0(station_name, "_EMODnet")  
  stat2 <- paste0(station_name, "_ICES") 
  
  # debugonce(check_pair_matrix)
  check_one_pair <- check_pair(
    id1 = stat1, id2 = stat2, 
    var_id = "Station_dataset", 
    dataset = data)
  
  result <- check_one_pair %>%
    mutate(id = station_name, .before = everything())
  
  if (!is.null(logfile)){
    logfile_conn <- file(logfile, "a")  # open for appending
    cat(station_name, "\n", file = logfile_conn)
    close(logfile_conn)
  }
  
  result
  
}

if (FALSE){
  # debugonce(check_one_pair_by_stationname)
  test <- check_one_pair_by_stationname("DMU Hornbk", dat_test)
}

#
# Plotting 'comparison values' (val1_c and val2_c)
#

plot_pair_by_stationname <- function(st, param_sel = NULL, checkdata = df_check_all){
  
  df_comparison <- checkdata %>%
    select(id, species, matrix, param, year, common_basis, val1_c, val2_c) %>%
    tidyr::pivot_longer(c(val1_c, val2_c), 
                        names_to = "Station", values_to = "Value") %>%
    filter(!is.na(Value)) 
  
  df_plot <- df_comparison %>%
    filter(id == st)  
  
  if (!is.null(param_sel)){
    df_plot <- df_plot %>% filter(param %in% param_sel)
  }
  
  ggplot(df_plot, aes(year, Value)) +
    geom_point(aes(color = Station, size = Station)) +
    scale_size_manual(values = c(2,1)) +    # make 'val1' dots bigger, so they are visible
    geom_text(aes(label = common_basis), vjust = -0.8) +
    # abs(title = param) +
    facet_wrap(vars(param), scales = "free_y") +
    labs(title = paste(st, 
                       paste(unique(df_plot$species), collapse = ", "),
                       paste(unique(df_plot$matrix), collapse = ", ")
    ))
  
}


#
# Get data for regression, using 'checkdata' output 
#   from either 'plot_pair_by_stationname' or 'plot_pair'
#
get_data_for_regression_checkdata <- function(id, id_column, species, param, data_dupl){
  
  param_selected = param
  species_selected = species
  
  sel <- data_dupl[[id_column]] %in% id
  df <- data_dupl[sel,] %>%
    filter(species %in% species_selected & param %in% param_selected)
  if (nrow(df) > 0){
    result <- df %>% mutate(
      # Assessment value  
      val_pick = case_when(
        !is.na(val2_as) ~ val2_as, 
        !is.na(val1_as) & val1_as == 0 ~ as.numeric(NA), 
        !is.na(val1_as) & val1_as > 0 ~ val1_as 
      ),
      # LOQ value  
      val_pick_flag = case_when(
        !is.na(val2_as) ~ ifelse(u_loq2/n2 > 0.5, "<", as.character(NA)), 
        !is.na(val1_as) & val1_as == 0 ~ as.character(NA), 
        !is.na(val1_as) & val1_as > 0 ~ ifelse(u_loq1/n1 > 0.5, "<", as.character(NA)) 
      ),
      # n value  
      val_pick_n = case_when(
        !is.na(val2_as) ~ n2, 
        !is.na(val1_as) & val1_as == 0 ~ as.integer(NA), 
        !is.na(val1_as) & val1_as > 0 ~ n1 
      )
    )
  } else {
    result <- NULL
  }
  result
}


#
# Get data for regression, using raw data
#
get_data_for_regression_rawdata <- function(id, id_column, species, param, data){
  
  param_selected = param
  sel <- data[[id_column]] %in% id
  
  data[sel,] %>%
    filter(WoRMS_scientificName == species & MATRX == matrix & PARAM == param & Year == year) %>%
    mutate(
      # Calculate dry and wet-weight basis (but DRYWT often lacking in one of the pair)
      Value_dw = case_when(
        BASIS %in% "W" & DRYWT_MUNIT %in% "%" ~ Value/(DRYWT*0.01),
        BASIS %in% "L" & DRYWT_MUNIT %in% "%" ~ (Value*FATWT_species*0.01)/(DRYWT*0.01),
        BASIS %in% "D" ~ Value),
      Value_ww = case_when(
        BASIS %in% "W" ~ Value,
        BASIS %in% "L" & FATWT_MUNIT %in% "%" ~ Value*(FATWT_species*0.01),
        BASIS %in% "D" & DRYWT_MUNIT %in% "%" ~ Value*(DRYWT*0.01))
    )
  
}

  