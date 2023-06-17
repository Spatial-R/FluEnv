
rt_predict_heatmap <- function(data,mid_point,region_var,rt_var){
  
dat.heatmap <- data
names(dat.heatmap)[which(names(dat.heatmap) == region_var)] <- "region"
names(dat.heatmap)[which(names(dat.heatmap) == rt_var)] <- "Rt"

fig_tem <- ggplot(data = dat.heatmap, aes(y = region, x = week)) + 
  geom_tile(aes(fill = Rt)) + 
  scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                       midpoint = mid_point,
                       name = "Standarized transmission intensity", 
                       guide = guide_colorbar(title.position = "right", 
                                              title.theme = element_text(angle = -90), 
                                              title.hjust = 0.5, 
                                              barwidth = 1, 
                                              barheight = 29)) + 
  scale_x_date(expand = c(0, 0), date_breaks = "3 month", 
               date_labels = "%Y-%m") + 
  scale_y_discrete(expand = c(0, 0)) + 
  labs(y = "",x ="Time (Weeks)"
       # , title = "Weekly Number of Influenza Cases per Sentinel for 47 Prefectures in Japan from 2012-09-02 to 2017-08-27"
  ) + 
  theme_bw() + 
  theme(legend.position = "right", 
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, -9), 
        axis.text.x = element_text(angle = 90, vjust = .5), 
        plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
        plot.title = element_text(vjust = 5, hjust = .5))
return(fig_tem)
}

normalize_target <- function(data){
  dat.tem <- data
  res_fin <- (dat.tem - min(dat.tem,na.rm = T))/(max(dat.tem,na.rm = T) - min(dat.tem,na.rm = T))
  return(res_fin)
}



rt_predict_heatmap_oneyear <- function(data,mid_point,region_var,rt_var){
  
  dat.heatmap <- data
  names(dat.heatmap)[which(names(dat.heatmap) == region_var)] <- "region"
  names(dat.heatmap)[which(names(dat.heatmap) == rt_var)] <- "Rt"
  
  fig_tem <- ggplot(data = dat.heatmap, aes(y = region, x = week)) + 
    geom_tile(aes(fill = Rt)) + 
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = mid_point,
                         name = "Standarized transmission intensity", 
                         guide = guide_colorbar(title.position = "right", 
                                                title.theme = element_text(angle = -90), 
                                                title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 29)) + 
    scale_y_discrete(expand = c(0, 0)) + 
    scale_x_continuous(breaks = seq(1,53,2),expand = c(0, 0)) +
    labs(x = "Time (Weeks)",y ="") + 
    theme_bw() + 
    theme(legend.position = "right", 
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -9), 
          axis.text.x = element_text(angle = 0, vjust = .5), 
          plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
          plot.title = element_text(vjust = 5, hjust = .5))
  return(fig_tem)
}


week_change <- function(data){
  week_dat <- data
  if(all(week_dat > 26)){
    median_week <- mean(week_dat)
  } else if (all(week_dat < 26)) {
    median_week <- mean(week_dat,na.rm = T)
  } else {
    min.week <- min(week_dat);  max.week <- max(week_dat); 
    if((max.week - min.week) > 25){
      tem_week <- ifelse(week_dat >26, 52- week_dat,week_dat)
      median_week <- mean(tem_week,na.rm = T)
    }
    median_week <- mean(week_dat,na.rm = T)
  }
  return(median_week)
}
