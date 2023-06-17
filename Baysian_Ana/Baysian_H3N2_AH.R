###########################################################################################################
###########################       Environmental factors on flu   ##########################################
###########################   Email: zhangbing4502431@outlook.com #########################################
###########################################################################################################

library(dplyr)
library(ggplot2)
library(reshape)
library(epical)
library(EpiEstim)
library(ggdendro)
library(cowplot)
library(mgcv)
library(splines)
library(lubridate)
library(pomp)
library(lme4)
#library(EFTR)
library(stringi)
library(humidity)
library(ggsci)
library(brms)
cols <- pal_npg("nrc")(9);#scales::show_col(cols)
load("Code-Data/Basic_All.RData")
source("Codes/functions.R")
influenza_type <- "H3_A"
DicreteSIDistr <- Dicret_SI(mu = 3.48,sigma=1.88,day = 10)


daily.case <- Weekly_To_Daily(data = model.dat,date.v = "week",case.v = influenza_type,
                              pop.v = "pop",region.v = "province")

######################################################################################################
####################### we need to change the value for Rc for the case with zero ####################
######################################################################################################

change.rc <- T
change.value <- 0.5

daily.rr.list <- lapply(unique(daily.case$region), function(region.id){
  print(region.id)
  dat.tem <- filter(daily.case,region == region.id)
  rr.tem <- Rc(data = dat.tem,date.id = "date",case.v = "case",n.num = 1)[,c(3,12,13)]
  rr.tem$region <- region.id
  names(rr.tem)[1] <- "rc"
  rr.tem <- na.omit(rr.tem)
  if(isTRUE(change.rc)){
    for (i in c(2:(nrow(rr.tem) -1))){
      #print(i)
      if(rr.tem[i,1] > 5 & rr.tem[i,3] < 3){
        rr.tem[i,1] <- change.value
      }
    }
    rr.tem[which(rr.tem$cases == 0),1] <- change.value
  }
  return(rr.tem)
})
daily.rr <- bind_rows(daily.rr.list)
daily.rr$rc <- ifelse(daily.rr$rc > 4,4,daily.rr$rc)

#cum.case <- Peak_Series(daily.case,region.v = "region",interval.v = 5,season.adj = T)
#result.tem <- merge(daily.rr,cum.case,by = c("region","date"),all.y = T)

names(dat.env.all)[2] <- "region"
result.tem <- merge(daily.rr,dat.env.all,by=c("date","region"),all.x = T)
analysis.tem <- result.tem[!is.na(result.tem$rc) & !is.na(result.tem$MT),]

#############################################################################################
###################################### Rt' approach #########################################
#############################################################################################

vaccnation <- read.csv("Data/vaccination-true.csv",header = T,fileEncoding = "gb18030")
vaccnation <- mutate(vaccnation,date = as.Date(date))[,-4]

epidemic.period <- epidemic_period(data = analysis.tem,season.period = c("09","08"),
                                   region.v = "region",interval.v = 10,method ="varying")
epidemic.period <- mutate(epidemic.period,id = paste(region,season,sep=""))

plot_period_list <- lapply(unique(analysis.tem$region),function(region_id){
  dat.tem <- filter(analysis.tem,region == region_id)
  epidemi_dat <- filter(epidemic.period,region == region_id)
  epidemic_period_vis(data = dat.tem,epidemic_date = epidemi_dat$date) +
    ggtitle(region_id)
})
#cowplot::plot_grid(plotlist = plot_period_list)


epidemic.period_list <- lapply(unique(epidemic.period$id),function(id_tem){
  dat.tem <- filter(epidemic.period, id == id_tem)
  date_max <- dat.tem[which.max(dat.tem$cases),"date"]
  date_range <- seq.Date(date_max - 6*7,date_max+6*7,"day")
  dat.tem.1 <- filter(dat.tem,date %in% date_range)
  dat.tem.1$cum <- c(0,cumsum(dat.tem.1$case)[-nrow(dat.tem.1)])
  return(dat.tem.1)
})
epidemic.period.1 <- bind_rows(epidemic.period_list)


plot_period_list <- lapply(unique(analysis.tem$region),function(region_id){
  dat.tem <- filter(analysis.tem,region == region_id)
  epidemi_dat <- filter(epidemic.period.1,region == region_id)
  epidemic_period_vis(data = dat.tem,epidemic_date = epidemi_dat$date) +
    ggtitle(region_id)
})
#cowplot::plot_grid(plotlist = plot_period_list)

analysis.tem.1 <- env_relate_rt(epidemic.period,intercept = T,holiday = F,spline_set = F,Rsquare = F)
#analysis.tem.1 <- env_relate_rt_glm(data = epidemic.period,type = "scs")
analysis.tem.1 <- filter(analysis.tem.1,adjust.rc < 4 & RH > 29)    ##### remove the reproduction above 4


########################################################################################################
################################## Rt' vsiualiztion ####################################################
########################################################################################################

influenza_type <- "H3N2_Adjust"

####################################################################################################
################## Relationship between Environemntal factors and flu ##############################
####################################################################################################

library(splines);library(ggeffects)
analysis.tem.3 <- analysis.tem.1[analysis.tem.1$adjust.rc < 2 & !is.na(analysis.tem.1$sunshine),]
analysis.tem.3 <- mutate(analysis.tem.3,id = paste(region,season,winter,sep=""))
analysis.tem.3 <- filter(analysis.tem.3,MT > -20 & MT < 32)
id_dat <- data.frame(table(analysis.tem.3$id))
id_remian <- as.character(filter(id_dat,Freq > 100)$Var1)
analysis.tem.3 <- filter(analysis.tem.3, id %in% id_remian)
analysis.tem.3 <- merge(analysis.tem.3,vaccnation,by = "date",all.x = T)
analysis.tem.3$work <- as.factor(analysis.tem.3$work)
analysis.tem.3 <- mutate(analysis.tem.3,idd = paste(season,winter,sep="-1"))
#analysis.tem.3 <- mutate(analysis.tem.3,AH = 1000*AH(WVP2(RH,SVP(C2K(MT)))))
write.csv(analysis.tem.3,file = paste("Code-Data/","Rt_Adjusted/",influenza_type,".csv",sep = ""),row.names = F)


for (i in c(unique(analysis.tem.3$season))){
  

(prior <- get_prior(bf(log(adjust.rc)~s(AH)+s(rain)+(1|region)+(1|idd)+ work),data = analysis.tem.3))
prior$prior[1] <- "student_t(10, 0, 1)"   
prior$prior[5] <- "student_t(10, 0, 1)"   
prior$prior[6] <- "student_t(10, 0, 1)"   
prior$prior[11] <- "student_t(10, 0, 1)"   
prior$prior[14] <- "student_t(10, 0, 1)"   

H3N2_AH_Rain_Region_ID_Work <- brm(bf(log(adjust.rc)~s(AH)+s(rain)+(1|region)+(1|idd)+work),
                                   data = analysis.tem.3,
                                   prior = prior,
                                   chains = 6,cores = 6,control = list(adapt_delta = 0.99),iter = 2000)
#save.image("mtrhrainsunAHine_idd.RData")
save(H3N2_AH_Rain_Region_ID_Work,file = paste("Code-Data/Year_Specific/",influenza_type,"AH_Rain_Region_ID_Work_Region.RData",sep = "_"))
sss
#########################  SUN, AH, Rain

(prior <- get_prior(bf(log(adjust.rc)~s(AH)+s(rain)+s(sunshine)+(1|region)+(1|idd)+work),
                    data = analysis.tem.3))
prior$prior[1] <- "student_t(10, 0, 1)"   
prior$prior[6] <- "student_t(10, 0, 1)"   
prior$prior[7] <- "student_t(10, 0, 1)"   
prior$prior[12] <- "student_t(10, 0, 1)"   
prior$prior[16] <- "student_t(10, 0, 1)"   

H3N2_AH_Rain_Sun_Region_ID_Work <- brm(bf(log(adjust.rc)~s(AH)+s(rain)+s(sunshine)+(1|region)+(1|idd)+work),
                                       data = analysis.tem.3,
                                       prior = prior,
                                       chains = 6,cores = 6,control = list(adapt_delta = 0.99),iter = 2000)
#save.image("mtrhrainsunshine_idd.RData")
save(H3N2_AH_Rain_Sun_Region_ID_Work,file = paste("Code-Data/Year_Specific/",influenza_type,"AH_Rain_Sun_Region_ID_Work_Region.RData",sep = "_"))

}