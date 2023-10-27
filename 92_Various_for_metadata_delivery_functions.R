
make_pie <- function(data, xvar, yvar, cols = c("1"="green", "2"="blue", "3"="red")){
  data$x <- factor(data[[xvar]])
  data$y <- data[[yvar]]
  ggplot(data, aes(x=1, y, fill = x)) + 
    geom_bar(stat="identity", color = "black", size = 1) + coord_polar(theta="y") +
    scale_fill_manual(values = cols) +
    theme_void() + theme(legend.position="none") + theme_transparent()
}



plot_contaminant <- function(param){
  
  pies <- dat_region_status %>%
    filter(PARAM == param) %>%
    split(~Region) %>%
    purrr::map(make_pie, xvar = "Status", yvar = "n")
  length(pies)
  
  # Make data set for pie center coordinates
  df_pies <- c(
    28.92, 56.97,     # Baltic
    8, 36,      # Mediterranean
    0.1, 66.1) %>%    # `North-East Atlantic Ocean` 
    matrix(ncol = 2, byrow = TRUE) %>%
    as.data.frame() %>%
    set_names(c("x", "y"))
  
  # Add pies as a variable
  df_pies$pie <- pies[c(1,3,4)]   # skip number 2, Black Sea
  
  # Sets size of pies
  df_pies$width = 12  
  df_pies$height = 12   
  
  # Plot points + pies for HG
  gg <- dat_status_trend %>%
    filter(PARAM == param) %>%
    mutate(Concentrations = case_when(
      Status == 1 ~ "Low", 
      Status == 2 ~ "Moderate", 
      Status == 3 ~ "High")) %>% 
    ggplot(aes(Longitude, Latitude, color = Concentrations)) +
    annotation_map(very_simple_map, fill = "navajowhite2") +
    geom_point() +
    annotate("text", x = -Inf, y = Inf, label = "Mercury", hjust = -0.3, vjust = 1.3, size = 7) +  # 
    geom_subview(data=df_pies, aes(x=x, y=y, subview=pie, width=width, height=height)) +
    easy_remove_axes() +
    theme(panel.background = element_rect(fill = "azure"))
  gg
  
}