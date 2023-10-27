
make_pie <- function(data, xvar, yvar, cols = c("1"="lightgreen", "2"="orange", "3"="red2")){
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
  df_pies$Region <- c("Baltic", "Mediterranean", "North-East Atlantic Ocean")
  
  # Add pies as a variable
  # df_pies$pie <- pies[c(1,3,4)]   # skip number 2, Black Sea
  df_pies$pie <- pies
  
  # Sets size of pies
  df_pies$width = 12  
  df_pies$height = 12   
  
  # Add trend metaanalysis result 
  df_pies <- df_pies %>%
    left_join(
      dat_region_trend %>% filter(PARAM == param),
      by = join_by(Region)
    ) %>%
    mutate(
      `Overall trend` = ifelse(is.na(`Overall trend`), "Trend unknown", `Overall trend`)
    )
  
  # Plot points + pies for HG
  gg <- dat_status_trend %>%
    filter(PARAM == param) %>%
    mutate(
      Concentrations = case_when(
        Status == 1 ~ "Low", 
        Status == 2 ~ "Moderate", 
        Status == 3 ~ "High"),
      Concentrations = factor(Concentrations, levels = c("Low", "Moderate", "High"))
    ) %>%
    ggplot(aes(Longitude, Latitude)) +
    annotation_map(very_simple_map, fill = "navajowhite2") +
    geom_point(aes(color = Concentrations)) +
    scale_color_manual(values = c(High="red2", Moderate="orange", Low="lightgreen")) +
    annotate("text", x = -Inf, y = Inf, label = param, hjust = -0.3, vjust = 1.3, size = 5) +  # 
    geom_subview(data=df_pies, aes(x=x, y=y, subview=pie, width=width, height=height)) +
    # The next line doesn't print anything visible, but needs to be there, otherwise the map "shrinks"
    geom_text(data = df_pies, aes(x = x, y = y-3.5, label = "")) +
    # Trend text:
    geom_label(data = df_pies, aes(x = x-3, y = y-3.5, label = `Overall trend`), size = 3) +
    coord_fixed(ratio = 1.4) +
    easy_remove_axes() +
    easy_remove_legend() +
    theme(panel.background = element_rect(fill = "azure"))
    # plot.margin = margin(0,32,0,0))
  gg
  
}