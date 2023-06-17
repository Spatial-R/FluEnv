
rm(list = ls())

library(dplyr)
library(ggplot2)
library(reshape)
library(epical)
library(EpiEstim)
library(ggdendro)
library(cowplot)
library(zoo)
library(mgcv)
library(splines)
library(lubridate)
library(ggsignif)
library(brms)
library(stringi)
library(humidity)
library(ggsci)
library(cowplot)
library(splines)

sun_ind <- TRUE
load("Data/Model/Basic_New.RData")
source("Codes/Functions/Epidemic_Peak_Predication_Functions.R")
cols <- pal_npg("nrc")(9);#scales::show_col(cols)
source("Codes/Functions/functions.R")

####################################################################################################
######################################  Enviromenta factors  #######################################
####################################################################################################

if(isTRUE(sun_ind)){
  
    load("Results/RH_Model/H3N2_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
    load("Process_Data/Environmental_Data.RData")
    dat.env.all <- mutate(dat.env.all,week = floor_date(date,"week"))
    dat.env.week <- aggregate(dat.env.all[,c(3:7)],
                              by = list(province = dat.env.all$province, 
                                        week = dat.env.all$week),mean,na.rm =T)
    dat.env.week <- mutate(dat.env.week,sunshine = sunshine*0.1, rain = ifelse(rain > 0,log(rain*0.1),-6))
    dat.tem <- dat.env.week
    dat.tem$idd <- rep(H3N2_RH_Rain_Sun_Region_ID_Work$data[1,"idd"],nrow(dat.tem))
    dat.tem$season <- rep(H3N2_RH_Rain_Sun_Region_ID_Work$data[1,"season"],nrow(dat.tem))
    dat.tem$work <- rep(H3N2_RH_Rain_Sun_Region_ID_Work$data[1,"work"],nrow(dat.tem))
    dat.tem$region_id <- rep(H3N2_RH_Rain_Sun_Region_ID_Work$data[1,"idd"],nrow(dat.tem))
    
    H3N2_RH_Rt <- exp(predict(H3N2_RH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/H3N2_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
    H3N2_AH_Rt <- exp(predict(H3N2_AH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/RH_Model/H1N1_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
    H1N1_RH_Rt <- exp(predict(H1N1_RH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/H1N1_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
    H1N1_AH_Rt <- exp(predict(H1N1_AH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/RH_Model/B_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
    B_RH_Rt <- exp(predict(B_RH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/B_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
    B_AH_Rt <- exp(predict(B_AH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/RH_Model/Total_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
    Total_RH_Rt <- exp(predict(Total_RH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/Total_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
    Total_AH_Rt <- exp(predict(Total_AH_Rain_Sun_Region_ID_Work, newdata = dat.tem))
    
} else {

    load("Results/RH_Model/H3N2_Adjust_RH_Rain_Region_ID_Work.RData")
    load("Process_Data/Environmental_Data.RData")
    dat.env.all <- mutate(dat.env.all,week = floor_date(date,"week"))
    dat.env.week <- aggregate(dat.env.all[,c(3:7)],
                              by = list(province = dat.env.all$province, 
                                            week = dat.env.all$week),mean,na.rm =T)
    dat.env.week <- mutate(dat.env.week,sunshine = sunshine*0.1, rain = ifelse(rain > 0,log(rain*0.1),-6))
    dat.tem <- dat.env.week
    dat.tem$idd <- rep(H3N2_RH_Rain_Region_ID_Work$data[1,"idd"],nrow(dat.tem))
    dat.tem$season <- rep(H3N2_RH_Rain_Region_ID_Work$data[1,"season"],nrow(dat.tem))
    dat.tem$work <- rep(H3N2_RH_Rain_Region_ID_Work$data[1,"work"],nrow(dat.tem))
    dat.tem$region_id <- rep(H3N2_RH_Rain_Region_ID_Work$data[1,"idd"],nrow(dat.tem))
    
    H3N2_RH_Rt <- exp(predict(H3N2_RH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/H3N2_Adjust_AH_Rain_Region_ID_Work.RData")
    H3N2_AH_Rt <- exp(predict(H3N2_AH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/RH_Model/H1N1_Adjust_RH_Rain_Region_ID_Work.RData")
    H1N1_RH_Rt <- exp(predict(H1N1_RH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/H1N1_Adjust_AH_Rain_Region_ID_Work.RData")
    H1N1_AH_Rt <- exp(predict(H1N1_AH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/RH_Model/B_Adjust_RH_Rain_Region_ID_Work.RData")
    B_RH_Rt <- exp(predict(B_RH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/B_Adjust_AH_Rain_Region_ID_Work.RData")
    B_AH_Rt <- exp(predict(B_AH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/RH_Model/Total_Adjust_RH_Rain_Region_ID_Work.RData")
    Total_RH_Rt <- exp(predict(Total_RH_Rain_Region_ID_Work, newdata = dat.tem))
    load("Results/AH_Model/Total_Adjust_AH_Rain_Region_ID_Work.RData")
    Total_AH_Rt <- exp(predict(Total_AH_Rain_Region_ID_Work, newdata = dat.tem))
}

dat.env.week_1 <- mutate(dat.env.week,
                        Total_AH = Total_AH_Rt[,1],Total_RH = Total_RH_Rt[,1],
                        H3N2_AH = H3N2_AH_Rt[,1],H3N2_RH = H3N2_RH_Rt[,1],
                        H1N1_AH = H1N1_AH_Rt[,1],H1N1_RH = H1N1_RH_Rt[,1],
                        B_AH = B_AH_Rt[,1],B_RH = B_RH_Rt[,1])

if(isTRUE(sun_ind)){
  write.csv(dat.env.week_1,file = "Results/Rt_Predication/China_Weekly_Rt_Sun.csv",row.names = F)
} else {
  write.csv(dat.env.week_1,file = "Results/Rt_Predication/China_Weekly_Rt.csv",row.names = F)
}


########################################################################################################
############################################### Heatmap  ###############################################
########################################################################################################

load("Process_Data/Basic_All.RData")

pro_china$name <- ifelse(as.character(pro_china$name) == "Shangxi","Shanxi",as.character(pro_china$name))
dat.env.week_1$province <- factor(dat.env.week_1$province,levels = as.character(rev(pro_china$name)))

dat.env.week_1 <- mutate(dat.env.week_1, week_num = week(dat.env.week_1$week))
dat.env.week_2 <- aggregate(dat.env.week_1[,c(8:15)],by = list(province = dat.env.week_1$province,
                                                              week_num  = dat.env.week_1$week_num),mean)

value.limits <- c(0.8,1.3)

(china_prdict_Total_RH <- ggplot(data = dat.env.week_2, aes(y = province, x = week_num)) + 
   #geom_tile(aes(fill = normalize_target(Total_RH))) + xlab("Time(weeks)") + ylab("")+
    geom_tile(aes(fill = Total_RH)) + xlab("Time(weeks)") + ylab("")+
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 1.1,#limits = value.limits,
                         name = expression(paste(R[t]^{"'"},sep="")), 
                         guide = guide_colorbar(title.position = "top", 
                                                #title.theme = element_text(angle = -90), 
                                                #title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 16)) +
    scale_x_continuous(breaks = seq(1,53,4),expand = c(0, 0)) + theme_bw(base_family = "serif",base_size = 11) +
   theme(legend.position = "right", 
         legend.margin = margin(0, 0, 0, 0),
         legend.box.margin = margin(0, 0, 0, -6), 
         axis.text.x = element_text(angle = 0, vjust = .5), 
         plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
         plot.title = element_text(vjust = 5, hjust = .5)))

(china_prdict_H3N2_RH <- ggplot(data = dat.env.week_1, aes(y = province, x = week_num)) + 
   # geom_tile(aes(fill = normalize_target(H3N2_RH))) + xlab("Time(weeks)") + ylab("")+
    geom_tile(aes(fill = H3N2_RH)) + xlab("Time(weeks)") + ylab("")+
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 1.1,#limits = value.limits,
                         name = expression(paste(R[t]^{"'"},sep="")), 
                         guide = guide_colorbar(title.position = "top", 
                                                #title.theme = element_text(angle = -90), 
                                                #title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 16)) +
    scale_x_continuous(breaks = seq(1,53,4),expand = c(0, 0)) + theme_bw(base_family = "serif",base_size = 11) +
    theme(legend.position = "right", 
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -6), 
          axis.text.x = element_text(angle = 0, vjust = .5), 
          plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
          plot.title = element_text(vjust = 5, hjust = .5)))


(china_prdict_H1N1_RH <- ggplot(data = dat.env.week_1, aes(y = province, x = week_num)) + 
    #geom_tile(aes(fill = normalize_target(H1N1_RH))) + xlab("Time(weeks)") + ylab("")+
    geom_tile(aes(fill = H1N1_RH)) + xlab("Time(weeks)") + ylab("")+
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 1.05,#limits = value.limits,
                         name = expression(paste(R[t]^{"'"},sep="")), 
                         guide = guide_colorbar(title.position = "top", 
                                                #title.theme = element_text(angle = -90), 
                                                #title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 16)) +
    scale_x_continuous(breaks = seq(1,53,4),expand = c(0, 0)) + theme_bw(base_family = "serif",base_size = 11) +
    theme(legend.position = "right", 
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -6), 
          axis.text.x = element_text(angle = 0, vjust = .5), 
          plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
          plot.title = element_text(vjust = 5, hjust = .5)))

(china_prdict_B_RH <- ggplot(data = dat.env.week_1, aes(y = province, x = week_num)) + 
    #geom_tile(aes(fill = normalize_target(B_RH))) + xlab("Time(weeks)") + ylab("")+
    geom_tile(aes(fill = B_RH)) + xlab("Time(weeks)") + ylab("")+
    scale_fill_gradient2(low = "white", high = "red", na.value = "transparent", 
                         midpoint = 1.1,#limits = value.limits,
                         name = expression(paste(R[t]^{"'"},sep="")), 
                         guide = guide_colorbar(title.position = "top", 
                                                #title.theme = element_text(angle = -90), 
                                                #title.hjust = 0.5, 
                                                barwidth = 1, 
                                                barheight = 16)) +
    scale_x_continuous(breaks = seq(1,53,4),expand = c(0, 0)) + theme_bw(base_family = "serif",base_size = 11) +
    theme(legend.position = "right", 
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(0, 0, 0, -6), 
          axis.text.x = element_text(angle = 0, vjust = .5), 
          plot.margin = unit(c(1.5, .5, 0, .5), "line"), 
          plot.title = element_text(vjust = 5, hjust = .5)))

(china_prdict <- plot_grid(china_prdict_Total_RH,china_prdict_H3N2_RH,china_prdict_H1N1_RH,china_prdict_B_RH,
          labels = c("Influenza","A(H3N2)","A(H1N1)pdm09","B-lineage")))

ggsave(china_prdict,filename = paste("Figures/","Rt_Prediction/China_MTRH.tiff",sep = ""),
       dpi = 300,units = "cm",width = 28,height = 20)


########################################################################################################
############################################### Three provinces  #######################################
########################################################################################################

province_dat <- filter(dat.env.week_1,province %in% c("Guangdong","Zhejiang","Beijing"))
province_dat$province <- factor(province_dat$province,levels = rev(c("Guangdong","Zhejiang","Beijing")))
province_plot <- melt(province_dat[,c(1,2,9,11,13,15)],id = c("province","week"))
province_plot$variable <- factor(province_plot$variable,levels =  c("Total_RH","H3N2_RH","H1N1_RH","B_RH"),
                                 labels = c("Total","A(H3N2)","A(H1N1)pdm09","B-lineage"))

plot_list <- lapply(c("Total","A(H3N2)","A(H1N1)pdm09","B-lineage"),function(data){
  dat_tem <- filter(province_plot,variable == data)
  ggplot(data = dat_tem,aes(x = week, y = value, color = factor(province))) + 
    #geom_point(size = 0.7) +
    geom_line() + xlab("") + ylab(expression(paste("Adjusted reproduction number (",R[t]^{"'"},")",sep=""))) +
    scale_color_manual(values = rev(cols[c(1,3,9)]),
                       guide = guide_legend(title = "")) +
    theme_classic(base_family = "serif",base_size = 11)+
    scale_x_date(date_breaks = "1 year",date_labels = "%Y") +
    theme(legend.position = c(0.1,0.9),plot.title = element_text(hjust = 0.5,face = "bold")) +
    scale_y_continuous(limits = c(0.96,1.3)) +
    ggtitle(data)
})

legend_province <- get_legend(
  plot_list[[1]] +
    guides(color = guide_legend(nrow = 1,title = "")) +
    theme(legend.position = "bottom",legend.box.margin = margin(-20, 0, 0, -6))
)

three_cities <- plot_grid(plot_list[[1]]+ theme(legend.position="none"),
                         plot_list[[2]]+ theme(legend.position="none"),
                         plot_list[[3]]+ theme(legend.position="none"),
                         plot_list[[4]]+theme(legend.position="none"),labels = c("A","B","C","D"))
(three_cities_new <- plot_grid(three_cities, legend_province, ncol = 1, rel_heights = c(1, 0.03)))

ggsave(three_cities_new,filename = paste("Figures/","Rt_Prediction/Three_cities.tiff",sep = ""),
       dpi = 300,units = "cm",width = 20,height = 16)


############################################################################################################
############################################### Correlation  ###############################################
############################################################################################################

Total_rho <- signif(cor.test(dat.env.week_1$Total_AH,dat.env.week_1$Total_RH)$estimate,3)
H3_rho <- signif(cor.test(dat_rt_sun$H3N2_RH,dat_rt$H3N2_RH)$estimate,3)
H1_rho <- signif(cor.test(dat_rt_sun$H1N1_RH,dat_rt$H1N1_RH)$estimate,3)
B_rho <- signif(cor.test(dat_rt_sun$B_RH,dat_rt$B_RH)$estimate,3)

Total_fit <- lm(Total_AH~Total_RH,data = dat.env.week_1);summary(Total_fit)
H3N2_fit <- lm(H3N2_AH~H3N2_RH,data = dat.env.week_1);summary(H3N2_fit)
H1N1_fit <- lm(H1N1_AH~H1N1_RH,data = dat.env.week_1);summary(H1N1_fit)
B_fit <- lm(B_AH~B_RH,data = dat.env.week_1);summary(B_fit)


fig_total <- ggplot(Total_fit$model,aes_string(x = names(Total_fit$model)[2],y = names(Total_fit$model)[1])) +
  geom_point(size =0.2) +
  stat_smooth(method = "lm",col = "red") +
  theme_classic(base_size = 12,base_family = "serif")+
  #annotate("text",x = 1.02,y = 1.13, label = paste("R2=",signif(summary(Total_fit)$adj.r.squared,2))) +
  annotate("text",x = 1.02,y = 1.13, label = paste0("rho ==",signif(Total_rho,2)),parse = TRUE) +
  xlab(expression(paste("Predicted ",R[t]^{"'"}," using RH Model",sep = ""))) +
  ylab(expression(paste("Predicted ",R[t]^{"'"}," using AH Model",sep = "")))

fig_h3n2 <- ggplot(H3N2_fit$model,aes_string(x = names(H3N2_fit$model)[2],y = names(H3N2_fit$model)[1])) +
  geom_point(size =0.2) +
  stat_smooth(method = "lm",col = "red") +
  theme_classic(base_size = 12,base_family = "serif")+
  #annotate("text",x = 1.06,y = 1.25, label = paste("R2=",signif(summary(H3N2_fit)$adj.r.squared,2))) +
  annotate("text",x = 1.06,y = 1.25, label = paste0("rho ==",signif(H3_rho,2)),parse = TRUE) +
  xlab(expression(paste("Predicted ",R[t]^{"'"}," using RH Model",sep = ""))) +
  ylab(expression(paste("Predicted ",R[t]^{"'"}," using AH Model",sep = "")))

fig_h1n1 <- ggplot(H1N1_fit$model,aes_string(x = names(H1N1_fit$model)[2],y = names(H1N1_fit$model)[1])) +
  geom_point(size =0.2) +
  stat_smooth(method = "lm",col = "red") +
  theme_classic(base_size = 12,base_family = "serif")+
  #annotate("text",x = 0.9,y = 1.115, label = paste("R2=",signif(summary(H1N1_fit)$adj.r.squared,2))) +
  annotate("text",x = 0.9,y = 1.115, label = paste0("rho ==",signif(H1_rho,2)),parse = TRUE)  +
  xlab(expression(paste("Predicted ",R[t]^{"'"}," using RH Model",sep = ""))) +
  ylab(expression(paste("Predicted ",R[t]^{"'"}," using AH Model",sep = "")))


fig_b <- ggplot(B_fit$model,aes_string(x = names(B_fit$model)[2],y = names(B_fit$model)[1])) +
  geom_point(size =0.2) +
  stat_smooth(method = "lm",col = "red") +
  theme_classic(base_size = 12,base_family = "serif")+
  #annotate("text",x = 1.04,y = 1.2, label = paste("R2=",signif(summary(B_fit)$adj.r.squared,2))) +
  annotate("text",x = 1.04,y = 1.2, label = paste0("rho ==",signif(B_rho,2)),parse = TRUE)  +
  xlab(expression(paste("Predicted ",R[t]^{"'"}," using RH Model",sep = ""))) +
  ylab(expression(paste("Predicted ",R[t]^{"'"}," using AH Model",sep = "")))


(fig_corre <- plot_grid(fig_total,fig_h3n2,fig_h1n1,fig_b,labels = c("Influenza","A(H3N2)","A(H1N1)pdm09","B-lineage")))

ggsave(fig_corre,filename = paste("Figures/","Rt_Prediction/Correlation.tiff",sep = ""),
       dpi = 300,units = "cm",width = 20,height = 16)


###############################################################################################################
#################################  predicted the incidence and transmisision rate #########################################
###############################################################################################################


source("Codes/Functions/duaggplot2.R")

plot_list <- lapply(c("Beijing","Zhejiang","Guangdong"),function(region.id){
  print(region.id)
  env.tem <- filter(province_plot,province == region.id)
  case.tem <- filter(model.dat,province == region.id)
  case.tem[,c(15,18:20)] <- case.tem[,c(15,18:20)]/case.tem$pop/10000
  case.tem.plot <- melt(case.tem[,c(1,15,18:20)],id = "week")
  case.tem.plot$variable <- factor(case.tem.plot$variable,levels = c("Total","H3_A","H1_A","B_A"),
                                   labels = c("Total","A(H3N2)","A(H1N1)pdm09","B-lineage"))
  plot.tem <- merge(env.tem,case.tem.plot,by = c("week","variable"))
  plot.tem <- filter(plot.tem,week > as.Date("2011-01-30"))
  flu_type_list <-lapply(c("A(H3N2)","A(H1N1)pdm09","B-lineage"),function(flu_id){
    cols_use <- cols[c(1,3,9)][which( c("A(H3N2)","A(H1N1)pdm09","B-lineage") == flu_id)]
    plot.tem.1 <- filter(plot.tem,variable == flu_id)
    p1 <- ggplot(data = plot.tem.1,aes(x = week,y = value.x,color = factor(variable))) + 
      geom_line(size = 0.4) +
      theme_bw(base_family = "serif",base_size = 10) + xlab("")+
      #xlab("Time (weeks)") +
      ylab("") + 
      scale_color_manual(values = cols_use,guide = "none") +
      geom_text(x = as.Date("2014-01-01"),y = 1.29,label = flu_id, fontface = "bold",family = "serif")+
      scale_y_continuous(limits = c(0.9,1.3))
    p2 <- ggplot(data = plot.tem.1,aes(x = week,y = value.y))+
      # geom_line(color = "black") +
      geom_bar(stat = "identity",fill = "black",width =3.2)+
      scale_y_continuous(limits = c(0,0.03),breaks = c(0,0.01,0.02,0.03))+
      theme_bw(base_size = 10,base_family = "serif")
  fig_mer <- dual_axis_graph(p1,p2,label_1 = (""),label_2 = c(""))
    return(fig_mer)
  })
  plot_grid(plotlist = flu_type_list,ncol = 3)
})

plot_all <- plot_grid(plotlist =plot_list,ncol = 1,labels = c("Beijing","Zhejiang","Guangdong"),label_size = 12)

ggsave(plot_all,filename = paste0("Figures/","Rt_Prediction/Inc_Rt.tiff"),
       dpi = 300,units = "cm",width = 20,height = 16)

save.image("Process_Data/Rt_Predication.RData")


##############################################################################
#######################  test for the variation explained by the climatic ###################
##############################################################################

scale_cases <- function(data){
  (data -min(data))/(max(data)-min(data))
}

case_tem_list <- lapply(unique(model.dat$province),function(region.id){
  print(region.id)
  env.tem <- filter(dat.env.week_1,province == region.id)
  case.tem <- filter(model.dat,province == region.id)
  target_col <- match(c(paste0(c("B_","H3_","H1_"),"A"),"Total"),names(case.tem))
  case.tem[,c(target_col)] <- case.tem[,c(target_col)]/case.tem$pop/10000
  #case.tem[,c(target_col)] <- apply(case.tem[,c(target_col)],2,scale_cases)
  plot.tem <- merge(env.tem,case.tem,by = c("week","province"))
  return(plot.tem)
})
case_model <- bind_rows(case_tem_list)

summary(lm(H1_A~H1N1_RH+factor(province),data = case_model))
summary(lm(H1_A~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3)+factor(province),data = case_model))

summary(lm(H3_A~H3N2_RH+factor(province),data = case_model))
summary(lm(H3_A~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3)+factor(province),data = case_model))

summary(lm(B_A~B_RH+factor(province),data = case_model))
summary(lm(B_A~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3)+factor(province),data = case_model))


target_cols <- c(cols[c(9,8,3,4)])


case_r2_list <- lapply(unique(model.dat$province),function(region.id){
  print(region.id)
  env.tem <- filter(dat.env.week_1,province == region.id)
  case.tem <- filter(model.dat,province == region.id)
  target_col <- match(c(paste0(c("B_","H3_","H1_"),"A"),"Total"),names(case.tem))
  case.tem[,c(target_col)] <- case.tem[,c(target_col)]/case.tem$pop/10000
  plot.tem <- merge(env.tem,case.tem,by = c("week","province"))
  h1_rt_rh_var <- summary(lm(H1_A~H1N1_RH,data = plot.tem))
  h1_rt_ah_var <- summary(lm(H1_A~H1N1_AH,data = plot.tem))
  h1_env_var <- summary(lm(H1_A~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3),data = plot.tem))
  
  h3_rt_rh_var <- summary(lm(H3_A~H3N2_RH,data = plot.tem))
  h3_rt_ah_var <- summary(lm(H3_A~H3N2_AH,data = plot.tem))
  h3_env_var <-summary(lm(H3_A~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3),data = plot.tem))
  
  b_rt_rh_var <- summary(lm(B_A~B_RH,data = plot.tem))
  b_rt_ah_var <- summary(lm(B_A~B_AH,data = plot.tem))
  b_env_var <-summary(lm(B_A~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3),data = plot.tem))
  
  total_rt_rh_var <- summary(lm(Total~Total_RH,data = plot.tem))
  total_rt_ah_var <- summary(lm(Total~Total_AH,data = plot.tem))
  total_env_var <- summary(lm(Total~ns(MT,3)+ns(RH,3)+ns(rain,3)+ns(sunshine,3),data = plot.tem))
  
  return(data.frame(h1_rt_rh = h1_rt_rh_var$adj.r.squared,h1_rt_ah = h1_rt_ah_var$adj.r.squared,
                    h3_rt_rh = h3_rt_rh_var$adj.r.squared,h3_rt_ah = h3_rt_ah_var$adj.r.squared,
                    b_rt_rh = b_rt_rh_var$adj.r.squared,b_rt_ah = b_rt_ah_var$adj.r.squared,
                    total_rt_rh = total_rt_rh_var$adj.r.squared,total_rt_ah = total_rt_ah_var$adj.r.squared,
                    h1_env = h1_env_var$adj.r.squared,h3_env = h3_env_var$adj.r.squared,
                    b_env = b_env_var$adj.r.squared,total_env = total_env_var$adj.r.squared,
                    province = region.id
                    ))
})
case_r2 <- bind_rows(case_r2_list)

case_r2_plot <- melt(case_r2[,c(1:8,13)],id = "province")
case_r2_plot$type <- ifelse(stri_detect_fixed(case_r2_plot$variable,"rt_rh"),"RH model","AH model")
case_r2_plot$variable <- unlist(lapply(case_r2_plot$variable,
                                       function(data)stri_split_fixed(data,"_")[[1]][1]))
case_r2_plot <- mutate(case_r2_plot,
                       variable = factor(variable,levels = c("h1","h3","b","total"),
                                         labels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","Influenza")))
case_r2_plot$value <- ifelse(case_r2_plot$value < 0,0,case_r2_plot$value)
anno_h1 <- wilcox.test(case_r2$h1_rt_rh,case_r2$h1_rt_ah)$p.value
anno_h3 <- wilcox.test(case_r2$h3_rt_rh,case_r2$h3_rt_ah)$p.value
anno_b <- wilcox.test(case_r2$b_rt_rh,case_r2$b_rt_ah)$p.value
anno_total <- wilcox.test(case_r2$total_rt_rh,case_r2$total_rt_ah)$p.value


Pt_size <- 2.8

(fig_sqr <- ggplot(data = case_r2_plot,aes(y = value,x = factor(variable))) + 
  geom_boxplot(aes(fill = factor(type)),position = "dodge") + 
  theme_classic(base_family = "serif",base_size = 10) +
  ylab("Adjusted R-square")+
  xlab("") +
  geom_signif(
      annotation = formatC(anno_h1, digits = 2),
      y_position = 0.4, xmin = 0.8, xmax = 1.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    geom_signif(
      annotation = formatC(anno_h3, digits = 2),
      y_position = 0.4, xmin = 1.8, xmax = 2.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    geom_signif(
      annotation = paste0(formatC(anno_b, digits = 2)),
      y_position = 0.4, xmin = 2.8, xmax = 3.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    geom_signif(
      annotation = formatC(anno_total, digits = 2),
      y_position = 0.6, xmin = 3.8, xmax = 4.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
  scale_fill_manual(values = cols[c(3,5)],guide = guide_legend(title = "")) +
  scale_y_continuous(limits = c(0,0.6))+
  theme(legend.position = c(0.2,0.9)))
  

ggsave(filename = paste0("Figures/","Rt_Prediction/R_square.tiff"),
       fig_sqr, width = 14,height = 10,units = "cm",dpi = 300)
       


north_name <- c("Beijing","Gansu","Tianjin","Hebei","Henan","Xinjiang",
                "Heilongjiang","Liaoning","Jilin","Inner Mongolia","Qinghai",
                "Ningxia","Shaanxi","Shandong","Shangxi")
case_r2_region <- mutate(case_r2,region = ifelse(province %in% north_name,"Northern provinces","Southern provinces"))
case_r2_region_plot <- melt(case_r2_region[,c(1,3,5,7,14)],id = "region")
case_r2_region_plot <- mutate(case_r2_region_plot,
                       variable = factor(variable,levels = c("h1_rt_rh","h3_rt_rh","b_rt_rh","total_rt_rh"),
                                         labels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","Influenza")))


(anno_reg_h1 <- wilcox.test(case_r2_region$h1_rt_rh~case_r2_region$region)$p.value)
(anno_reg_h3 <- wilcox.test(case_r2_region$h3_rt_rh~case_r2_region$region)$p.value)
(anno_reg_b <- wilcox.test(case_r2_region$b_rt_rh~case_r2_region$region)$p.value)
(anno_reg_total <- wilcox.test(case_r2_region$total_rt_rh~case_r2_region$region)$p.value)


Pt_size <- 2.8

(fig_reg_sqr <- ggplot(data = case_r2_region_plot,aes(y = value,x = factor(variable))) + 
    geom_boxplot(aes(fill = factor(region)),position = "dodge") + 
    theme_classic(base_family = "serif",base_size = 10) +
    ylab("Adjusted R-square")+
    xlab("") +
    geom_signif(
      annotation = formatC(anno_reg_h1, digits = 2),
      y_position = 0.4, xmin = 0.8, xmax = 1.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    geom_signif(
      annotation = formatC(anno_reg_h3, digits = 2),
      y_position = 0.4, xmin = 1.8, xmax = 2.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    geom_signif(
      annotation = paste0(formatC(anno_reg_b, digits = 2)),
      y_position = 0.4, xmin = 2.8, xmax = 3.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    geom_signif(
      annotation = formatC(anno_reg_total, digits = 2),
      y_position = 0.6, xmin = 3.8, xmax = 4.2,
      tip_length = c(0.1, 0.1),textsize = Pt_size
    )+
    scale_fill_manual(values = cols[c(3,5)],guide = guide_legend(title = "")) +
    scale_y_continuous(limits = c(0,0.6))+
    theme(legend.position = c(0.3,0.9)))


ggsave(filename = paste0("Figures/","Rt_Prediction/R_square_region.tiff"),
       fig_reg_sqr, width = 14,height = 10,units = "cm",dpi = 300)
