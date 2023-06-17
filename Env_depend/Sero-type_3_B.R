library(ggplot2)
library(dplyr)
library(ggsci)
library(lme4)
library(splines)
library(ggeffects)
library(cowplot)
library(brms)

cols <- pal_npg("nrc")(9);scales::show_col(cols)
y_limits <- c(0.8,1.2)  
normolized <- function(data){
  dat_res <- (data - min(data))/(max(data) - min(data))
  return(dat_res)
}

load("Data/Model/Basic_New.RData")
load("Results/RH_Model/B_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
load("Results/RH_Model/H1N1_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
load("Results/RH_Model/H3N2_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
load("Results/AH_Model/H3N2_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
load("Results/SH_Model/H3N2_Adjust_SH_Rain_Sun_Region_ID_Work.RData")
load("Results/AH_Model/H1N1_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
load("Results/SH_Model/H1N1_Adjust_SH_Rain_Sun_Region_ID_Work.RData")
load("Results/AH_Model/B_Adjust_AH_Rain_Sun_Region_ID_Work.RData")
load("Results/SH_Model/B_Adjust_SH_Rain_Sun_Region_ID_Work.RData")
load("Results/RH_Model/Total_Adjust_RH_Rain_Sun_Region_ID_Work.RData")
load("Results/AH_Model/Total_Adjust_AH_Rain_Sun_Region_ID_Work.RData")

total_dat_rh <- marginal_smooths(Total_RH_Rain_Sun_Region_ID_Work)
total_dat_ah <- marginal_smooths(Total_AH_Rain_Sun_Region_ID_Work)
h1n1_dat_rh <- marginal_smooths(H1N1_RH_Rain_Sun_Region_ID_Work)
h1n1_dat_ah <- marginal_smooths(H1N1_AH_Rain_Sun_Region_ID_Work)
h3n2_dat_rh <- marginal_smooths(H3N2_RH_Rain_Sun_Region_ID_Work)
h3n2_dat_ah <- marginal_smooths(H3N2_AH_Rain_Sun_Region_ID_Work)
b_dat_rh <- marginal_smooths(B_RH_Rain_Sun_Region_ID_Work)
b_dat_ah <- marginal_smooths(B_AH_Rain_Sun_Region_ID_Work)
b_dat_sh <- marginal_smooths(B_SH_Rain_Sun_Region_ID_Work)
h1n1_dat_sh <- marginal_smooths(H1N1_SH_Rain_Sun_Region_ID_Work)
h3n2_dat_sh <- marginal_smooths(H3N2_SH_Rain_Sun_Region_ID_Work)


h3n2_rh <- h3n2_dat_rh$`mu: s(RH)`
h1n1_rh <- h1n1_dat_rh$`mu: s(RH)`
b_rh <- b_dat_rh$`mu: s(RH)`
total_rh <- total_dat_rh$`mu: s(RH)`
  
h3n2_mt <- h3n2_dat_rh$`mu: s(MT)`
h1n1_mt <- h1n1_dat_rh$`mu: s(MT)`
b_mt <- b_dat_rh$`mu: s(MT)`
total_mt <- total_dat_rh$`mu: s(MT)`

h3n2_rain <- h3n2_dat_rh$`mu: s(rain)`
h1n1_rain <- h1n1_dat_rh$`mu: s(rain)`
b_rain <- b_dat_rh$`mu: s(rain)`
total_rain <- total_dat_rh$`mu: s(rain)`

h3n2_sun <- h3n2_dat_rh$`mu: s(sunshine)` 
h1n1_sun <- h1n1_dat_rh$`mu: s(sunshine)`
b_sun <- b_dat_rh$`mu: s(sunshine)`
total_sun <- total_dat_rh$`mu: s(sunshine)`

h3n2_sh <- h3n2_dat_sh$`mu: s(SH)` 
h1n1_sh <- h1n1_dat_sh$`mu: s(SH)`
b_sh <- b_dat_sh$`mu: s(SH)`

h3n2_ah <- h3n2_dat_ah$`mu: s(AH)`
h1n1_ah <- h1n1_dat_ah$`mu: s(AH)`
b_ah <- b_dat_ah$`mu: s(AH)`
total_ah <- total_dat_ah$`mu: s(AH)`


h3n2_ah_sun <- h3n2_dat_ah$`mu: s(sunshine)`
h1n1_ah_sun <- h1n1_dat_ah$`mu: s(sunshine)`
b_ah_sun <- b_dat_ah$`mu: s(sunshine)`
total_ah_sun <- total_dat_ah$`mu: s(sunshine)`
h3n2_ah_rain <- h3n2_dat_ah$`mu: s(rain)`
h1n1_ah_rain<- h1n1_dat_ah$`mu: s(rain)`
b_ah_rain <- b_dat_ah$`mu: s(rain)`
total_ah_rain <- total_dat_ah$`mu: s(rain)`


##################################################################################################

mt_data = b_mt;rh_data = b_rh;rain_data = b_rain;sun_data = b_sun;ah_data = b_ah

combn_plot <- function(mt_data,rh_data,rain_data,sun_data,ah_data,y_ranges = y_limits){
  
  if(!is.null(rh_data)){
    figrh <- ggplot(data = rh_data,aes(x= RH, y= exp(estimate__))) + geom_line(color = "red") +
      geom_ribbon(data=rh_data,aes(ymin=exp(lower__),ymax=exp(upper__)),alpha=0.3,fill = "blue") +
      #geom_hline(yintercept = 1,linetype =2) +
      theme_classic(base_size = 12,base_family = "serif") +
      xlab("Relative humidity (%)") + 
      ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'",")",sep = ""))) +
      scale_y_continuous(limits = y_ranges)
  } else {
    figrh <- NULL
  } 
  
  if(!is.null(mt_data)){
    figmt <- ggplot(data = mt_data,aes(x= MT,y=exp(estimate__))) + geom_line(color = "red") +
      geom_ribbon(data=mt_data,aes(ymin=exp(lower__),ymax=exp(upper__)),alpha=0.3,fill = "blue") +
      #geom_hline(yintercept = 1,linetype =2) +
      xlab(expression(paste("Mean temperature (",degree,"C)",sep = ""))) + 
      ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
      theme_classic(base_size = 12,base_family = "serif") +
      scale_y_continuous(limits = y_ranges)
  } else {
    figmt <- NULL
  }
  
  if(!is.null(sun_data)){
    figss <- ggplot(data = sun_data,aes(x= sunshine, y=exp(estimate__))) + geom_line(color = "red") +
      geom_ribbon(data=sun_data,aes(ymin=exp(lower__),ymax=exp(upper__)),alpha=0.3,fill = "blue") +
      # geom_hline(yintercept = 1,linetype =2) +
      xlab("Sunshine hours (hours/day)") + 
      ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
      theme_classic(base_size = 12,base_family = "serif") +
      scale_y_continuous(limits = y_ranges)
  }else {
    figss<- NULL
  }
  
  if(!is.null(rain_data)){
    figrn <- ggplot(data = rain_data,aes(x= rain, y=exp(estimate__))) + geom_line(color = "red") +
      geom_ribbon(data=rain_data,aes(ymin=exp(lower__),ymax=exp(upper__)),alpha=0.3,fill = "blue") +
      #geom_hline(yintercept = 1,linetype =2) +
      xlab("Precipitation(mm)") + 
      ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
      #scale_x_continuous(trans = "log")
      scale_x_continuous(breaks = c(-6,-3,0,3),labels = c(expression(paste(e^-6,"",sep = "")),
                                                          expression(paste(e^-3,"",sep = "")),
                                                          expression(paste(e^0,"",sep = "")),
                                                          expression(paste(e^3,"",sep = ""))))+
      theme_classic(base_size = 12,base_family = "serif") +
      scale_y_continuous(limits = y_ranges)
  } else {
    figrn <- NULL
  }
  
  if(!is.null(ah_data)){
    figah <- ggplot(data = ah_data,aes(x= AH, y=exp(estimate__))) + geom_line(color = "red") +
      geom_ribbon(data=ah_data,aes(ymin=exp(lower__),ymax=exp(upper__)),alpha=0.3,fill = "blue") +
      #geom_hline(yintercept = 1,linetype =2) +
      xlab(expression(paste("Absolute humidity (g/",m^3,")",sep="") )) + 
      ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
      theme_classic(base_size = 12,base_family = "serif") +
      scale_y_continuous(limits = y_ranges)
  } else {
    figah <- NULL
  }
  return(list(figmt,figrh,figss,figrn,figah))
}

h3n2_plot_list <- combn_plot(mt_data = h3n2_mt,rh_data = h3n2_rh,rain_data = h3n2_rain,
                             sun_data = h3n2_sun,ah_data = h3n2_ah)
h3n2_plot <- plot_grid(plotlist = h3n2_plot_list,ncol=1)

h1n1_plot_list <- combn_plot(mt_data = h1n1_mt,rh_data = h1n1_rh,rain_data = h1n1_rain,
                             sun_data = h1n1_sun,ah_data = h1n1_ah)
h1n1_plot <- plot_grid(plotlist = h1n1_plot_list,ncol=1)

b_plot_list <- combn_plot(mt_data = b_mt,rh_data = b_rh,rain_data = b_rain,
                          sun_data = b_sun,ah_data = b_ah)
b_plot <- plot_grid(plotlist = b_plot_list,ncol=1)

plot_all <- plot_grid(h3n2_plot,h1n1_plot,b_plot,ncol = 3)

ggsave(plot_all,filename = paste("Figures/","EnvDep/All.tiff",sep=""),
       width = 20,height = 30,units = "cm",dpi = 300)


total_plot_list <- combn_plot(mt_data = total_mt,rh_data = total_rh,rain_data = total_rain,
                              sun_data = total_sun,ah_data = total_ah,y_ranges = c(0.8,1.1))
total_plot <- plot_grid(plotlist = total_plot_list)

ggsave(total_plot,filename = paste("Figures/","EnvDep/Total.tiff",sep=""),
       width = 16,height = 14,units = "cm",dpi = 300)


target_cols <- cols[c(9,8,3)]


h3n2_mt$type <- "A(H3N2)";h1n1_mt$type <- "A(H1N1)pdm09";b_mt$type <- "B-lineage"
mt_dat <- rbind(h3n2_mt,h1n1_mt,b_mt)
sero_mt <- ggplot(data = mt_dat,aes(x= MT, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(data=mt_dat,aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  xlab(expression(paste("Mean temperature (",degree,"C)",sep = ""))) + 
  scale_y_continuous(limits = y_limits) +
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))

h3n2_rh$type <- "A(H3N2)";h1n1_rh$type <- "A(H1N1)pdm09";b_rh$type <- "B-lineage"
rh_dat <- rbind(h3n2_rh,h1n1_rh,b_rh)
sero_rh <- ggplot(data = rh_dat,aes(x= RH, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  xlab("Relative humidity (%)") + 
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  scale_y_continuous(limits = y_limits) +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))


h3n2_rain$type <- "A(H3N2)";h1n1_rain$type <- "A(H1N1)pdm09";b_rain$type <- "B-lineage"
rain_dat <- rbind(h3n2_rain,h1n1_rain,b_rain)
sero_rain <- ggplot(data = rain_dat,aes(x= rain, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Precipitation (mm)") + 
  scale_y_continuous(limits = y_limits) +
  scale_x_continuous(breaks = c(-6,-3,0,3),labels = c(expression(paste(e^-6,"",sep = "")),
                                                      expression(paste(e^-3,"",sep = "")),
                                                      expression(paste(e^0,"",sep = "")),
                                                      expression(paste(e^3,"",sep = ""))))+
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))


h3n2_sun$type <- "A(H3N2)";h1n1_sun$type <- "A(H1N1)pdm09";b_sun$type <- "B-lineage"
sun_dat <- rbind(h3n2_sun,h1n1_sun,b_sun)
sun_dat <- filter(sun_dat,sunshine < 12.2)
sero_sun <- ggplot(data = sun_dat,aes(x= sunshine, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Sunshine hours (hours/day)") + 
  scale_y_continuous(limits = y_limits)+
  scale_x_continuous(breaks = c(0,5,10))+
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))

h3n2_ah$type <- "A(H3N2)";h1n1_ah$type <- "A(H1N1)pdm09";b_ah$type <- "B-lineage"
ah_dat <- rbind(h3n2_ah,h1n1_ah,b_ah)
sero_ah <- ggplot(data = ah_dat,aes(x= AH, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab(expression(paste("Absolute humidity (g/",m^3,")",sep="") )) + 
  scale_y_continuous(limits = y_limits)+
  scale_color_manual(values = target_cols,guide = guide_legend(title = "",keyheight = 2)) +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.5,0.6))


fig_merge <- plot_grid(sero_mt,sero_rh,sero_sun,sero_rain,
          sero_ah+theme(legend.position = "none"),get_legend(sero_ah))

ggsave(fig_merge,filename = paste("Figures/","EnvDep/Sub_Types.tiff",sep=""),
       width = 20,height = 17,units = "cm",dpi = 300)


#################################  Merge the total #######################################


target_cols <- c(cols[c(9,8,3,4)])


h3n2_mt$type <- "A(H3N2)";h1n1_mt$type <- "A(H1N1)pdm09";b_mt$type <- "B-lineage";total_mt$type <- "All subtypes"
mt_dat <- rbind(h3n2_mt,h1n1_mt,b_mt,total_mt)
mt_dat <- mutate(mt_dat,type = factor(type,levels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","All subtypes")))

sero_mt <- ggplot(data = mt_dat,aes(x= MT, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(data=mt_dat,aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  xlab(expression(paste("Mean temperature (",degree,"C)",sep = ""))) + 
  scale_y_continuous(limits = y_limits) +
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))

h3n2_rh$type <- "A(H3N2)";h1n1_rh$type <- "A(H1N1)pdm09";b_rh$type <- "B-lineage";total_rh$type <- "All subtypes"
rh_dat <- rbind(h3n2_rh,h1n1_rh,b_rh,total_rh)
sero_rh <- ggplot(data = rh_dat,aes(x= RH, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  xlab("Relative humidity (%)") + 
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  scale_y_continuous(limits = y_limits) +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))


h3n2_rain$type <- "A(H3N2)";h1n1_rain$type <- "A(H1N1)pdm09";b_rain$type <- "B-lineage";total_rain$type <- "All subtypes"
rain_dat <- rbind(h3n2_rain,h1n1_rain,b_rain,total_rain)
rain_dat <- mutate(rain_dat,type = factor(type,levels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","All subtypes")))

sero_rain <- ggplot(data = rain_dat,aes(x= rain, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Precipitation (mm)") + 
  scale_y_continuous(limits = y_limits) +
  scale_x_continuous(breaks = c(-6,-3,0,3),labels = c(expression(paste(e^-6,"",sep = "")),
                                                      expression(paste(e^-3,"",sep = "")),
                                                      expression(paste(e^0,"",sep = "")),
                                                      expression(paste(e^3,"",sep = ""))))+
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))


h3n2_sun$type <- "A(H3N2)";h1n1_sun$type <- "A(H1N1)pdm09";b_sun$type <- "B-lineage"; total_sun$type <- "All subtypes"
sun_dat <- rbind(h3n2_sun,h1n1_sun,b_sun,total_sun)
sun_dat <- filter(sun_dat,sunshine < 12.2)
sun_dat <- mutate(sun_dat,type = factor(type,levels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","All subtypes")))

sero_sun <- ggplot(data = sun_dat,aes(x= sunshine, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Sunshine hours (hours/day)") + 
  scale_y_continuous(limits = y_limits)+
  scale_x_continuous(breaks = c(0,5,10))+
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))

h3n2_ah$type <- "A(H3N2)";h1n1_ah$type <- "A(H1N1)pdm09";b_ah$type <- "B-lineage";total_ah$type <- "All subtypes"
ah_dat <- rbind(h3n2_ah,h1n1_ah,b_ah,total_ah)
ah_dat <- mutate(ah_dat,type = factor(type,levels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","All subtypes")))
sero_ah <- ggplot(data = ah_dat,aes(x= AH, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab(expression(paste("Absolute humidity (g/",m^3,")",sep="") )) + 
  scale_y_continuous(limits = y_limits)+
  scale_color_manual(values = target_cols,guide = guide_legend(title = "",keyheight = 2)) +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.5,0.6))


fig_merge_rh <- plot_grid(sero_mt,sero_rh,sero_sun,sero_rain,ncol = 2)
fig_merge_ah <- plot_grid(sero_ah+theme(legend.position = "none"),get_legend(sero_ah),ncol = 1)
fig_merge <- plot_grid(fig_merge_rh,fig_merge_ah,rel_widths = c(0.6,0.3),labels = LETTERS[1:2])


ggsave(fig_merge,filename = paste("Figures/","EnvDep/Sero_type_flu.tiff",sep=""),
       width = 20,height = 17,units = "cm",dpi = 300)


ggsave(sero_mt+scale_y_continuous(limits = c(0.85,1.15)),
       filename = paste("Figures/","EnvDep/Sero_mt.tiff",sep=""),
       width = 10,height = 9,units = "cm",dpi = 300)


h3n2_ah_rain$type <- "A(H3N2)";h1n1_ah_rain$type <- "A(H1N1)pdm09";b_ah_rain$type <- "B-lineage"
total_ah_rain$type <- "All subtypes"
rain_ah_dat <- rbind(h3n2_rain,h1n1_rain,b_rain,total_ah_rain)
rain_ah_dat <- mutate(rain_ah_dat,type = factor(type,levels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","All subtypes")))

sero_ah_rain <- ggplot(data = rain_ah_dat,aes(x= rain, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Precipitation (mm)") + 
  scale_y_continuous(limits = y_limits) +
  scale_x_continuous(breaks = c(-6,-3,0,3),labels = c(expression(paste(e^-6,"",sep = "")),
                                                      expression(paste(e^-3,"",sep = "")),
                                                      expression(paste(e^0,"",sep = "")),
                                                      expression(paste(e^3,"",sep = ""))))+
  scale_color_manual(values = target_cols,guide = "none") +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.7,0.9))


h3n2_ah_sun$type <- "A(H3N2)";h1n1_ah_sun$type <- "A(H1N1)pdm09";b_ah_sun$type <- "B-lineage"
total_ah_sun$type <- "All subtypes"
sun_ah_dat <- rbind(h3n2_sun,h1n1_sun,b_sun,total_ah_sun)
sun_ah_dat <- filter(sun_ah_dat,sunshine < 12.2)
sun_ah_dat <- mutate(sun_ah_dat,type = factor(type,levels = c("A(H1N1)pdm09","A(H3N2)","B-lineage","All subtypes")))

sero_ah_sun <- ggplot(data = sun_ah_dat,aes(x= sunshine, y = exp(estimate__),group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 1) +
  geom_ribbon(aes(ymin=exp(lower__),ymax=exp(upper__),fill = factor(type)),alpha=0.3)+
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Adjusted reproduction number ( ",R[t]^"'"," )",sep = "")))+
  theme_classic(base_size = 12,base_family = "serif") +
  xlab("Sunshine hours (hours/day)") + 
  scale_y_continuous(limits = y_limits)+
  scale_x_continuous(breaks = c(0,5,10))+
  scale_color_manual(values = target_cols,guide = guide_legend(title = "",keyheight = 1.3)) +
  scale_fill_manual(values = target_cols,guide = "none") +
  theme(legend.position = c(0.35,0.85))

fig_merge_ah <- plot_grid(sero_ah_rain,sero_ah_sun,ncol = 2)
ggsave(fig_merge_ah,filename = paste("Figures/","EnvDep/Sun_rain_AH.tiff",sep=""),
       width = 15,height = 10,units = "cm",dpi = 300)



#################################################################################################
#################################################################################################
#################################################################################################

h3n2_mt_dat <- data.frame(x = h3n2_mt$MT,y = h3n2_mt$estimate__,type = "A(H3N2)")
h1n1_mt_dat <- data.frame(x = h1n1_mt$MT,y = h1n1_mt$estimate__,type = "A(H1N1)pdm09")
b_mt_dat <- data.frame(x = b_mt$MT,y = b_mt$estimate__,type = "B-lineage")

mt_dat <- rbind(h3n2_mt_dat,h1n1_mt_dat,b_mt_dat)

sero_mt <- ggplot(data = mt_dat,aes(x= x, y=y,group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 2) +
  xlab(expression(paste("Mean temperature (",degree,"C)",sep = ""))) + 
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  theme(legend.position = c(0.7,0.9))
#ggsave(sero_mt,filename = paste("Figures/","Serotype/MT.tiff",sep=""),
 #      width = 16,height = 16,units = "cm",dpi = 300)


h3n2_rh_dat <- data.frame(x = h3n2_rh$RH,y = h3n2_rh$estimate__,type = "A(H3N2)")
h1n1_rh_dat <- data.frame(x = h1n1_rh$RH,y = h1n1_rh$estimate__,type = "A(H1N1)pdm09")
b_rh_dat <- data.frame(x = b_rh$RH,y = b_rh$estimate__,type = "B-lineage")

rh_dat <- rbind(h3n2_rh_dat,h1n1_rh_dat,b_rh_dat)

sero_rh <- ggplot(data = rh_dat,aes(x= x, y=y,group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 2) +
  xlab("Relative humidity (%)")  + 
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  theme(legend.position = c(0.7,0.9))
#ggsave(sero_rh,filename = paste("Figures/","Serotype/RH.tiff",sep=""),
#       width = 16,height = 16,units = "cm",dpi = 300)


h3n2_rain_dat <- data.frame(x = h3n2_rain$rain,y = h3n2_rain$estimate__,type = "A(H3N2)")
h1n1_rain_dat <- data.frame(x = h1n1_rain$rain,y = h1n1_rain$estimate__,type = "A(H1N1)pdm09")
b_rain_dat <- data.frame(x = b_rain$rain,y = b_rain$estimate__,type = "B-lineage")

rain_dat <- rbind(h3n2_rain_dat,h1n1_rain_dat,b_rain_dat)

sero_rain <- ggplot(data = rain_dat,aes(x= x, y=y,group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 2) +
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  xlab("Precipitation(mm)") + 
  #scale_x_continuous(trans = "log")
  scale_x_continuous(breaks = c(-6,-3,0,3),labels = c(expression(paste(e^-6,"",sep = "")),
                                                      expression(paste(e^-3,"",sep = "")),
                                                      expression(paste(e^0,"",sep = "")),
                                                      expression(paste(e^3,"",sep = ""))))+
  theme_classic(base_size = 12,base_family = "serif")+
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  theme(legend.position = c(0.7,0.85))
#ggsave(sero_rain,filename = paste("Figures/","Serotype/Rain.tiff",sep=""),
#       width = 16,height = 16,units = "cm",dpi = 300)

h3n2_sun_dat <- data.frame(x = h3n2_sun$sunshine,y = h3n2_sun$estimate__,type = "A(H3N2)")
h1n1_sun_dat <- data.frame(x = h1n1_sun$sunshine,y = h1n1_sun$estimate__,type = "A(H1N1)pdm09")
b_sun_dat <- data.frame(x = b_sun$sunshine,y = b_sun$estimate__,type = "B-lineage")

sun_dat <- rbind(h3n2_sun_dat,h1n1_sun_dat,b_sun_dat)

sero_sun <- ggplot(data = sun_dat,aes(x= x, y=y,group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 2) +
  xlab("Sunshine hours(/hours/day)") + 
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  theme(legend.position = c(0.9,0.9))
#ggsave(sero_sun,filename = paste("Figures/","Serotype/Sun.tiff",sep=""),
 #      width = 16,height = 16,units = "cm",dpi = 300)




h3n2_ah_dat <- data.frame(x = h3n2_ah$AH,y = h3n2_ah$estimate__,type = "A(H3N2)")
h1n1_ah_dat <- data.frame(x = h1n1_ah$AH,y = h1n1_ah$estimate__,type = "A(H1N1)pdm09")
b_ah_dat <- data.frame(x = b_ah$AH,y = b_ah$estimate__,type = "B-lineage")

ah_dat <- rbind(h3n2_ah_dat,h1n1_ah_dat,b_ah_dat)

sero_ah <- ggplot(data = ah_dat,aes(x= x, y=y,group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 2) +
  xlab(expression(paste("Absolute humidity (g/",m^3,")",sep="") )) + 
  #ylab(expression(paste("Normolized adj  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  theme(legend.position = c(0.7,0.9))
#ggsave(sero_ah,filename = paste("Figures/","Serotype/AH.tiff",sep=""),
 #      width = 16,height = 16,units = "cm",dpi = 300)



h3n2_sh_dat <- data.frame(x = h3n2_sh$SH,y = h3n2_sh$estimate__,type = "A(H3N2)")
h1n1_sh_dat <- data.frame(x = h1n1_sh$SH,y = h1n1_sh$estimate__,type = "A(H1N1)pdm09")
b_sh_dat <- data.frame(x = b_sh$SH,y = b_sh$estimate__,type = "B-lineage")

sh_dat <- rbind(h3n2_sh_dat,h1n1_sh_dat,b_sh_dat)

sero_sh <- ggplot(data = sh_dat,aes(x= x, y=y,group = factor(type))) + 
  geom_line(aes(color = factor(type)),size = 2) +
  xlab(expression(paste("Specific humidity (g/",kg,")",sep="") )) + 
  #ylab(expression(paste("Normolized adjusted reproduction number ( ",R[t]^"'"," )",sep = ""))) +
  ylab(expression(paste("Normolized  ",R[t]^"'","",sep = ""))) +
  theme_classic(base_size = 12,base_family = "serif") +
  #scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = guide_legend(title = "")) +
  scale_color_manual(values = c(cols[c(8,9,3)],"black"),guide = "none") +
  theme(legend.position = c(0.7,0.9))
#ggsave(sero_ah,filename = paste("Figures/","Serotype/SH.tiff",sep=""),
 #      width = 16,height = 16,units = "cm",dpi = 300)

fig_all <- plot_grid(sero_mt,sero_rh,sero_rain,sero_sun,sero_ah,sero_sh)

ggsave(fig_all,filename = paste("Figures/","EnvDep/Sero_type_Scaled.tiff",sep=""),
       width = 24,height = 16,units = "cm",dpi = 300)

save.image("Process_Data/Serotype_Nonlinear.RData")
