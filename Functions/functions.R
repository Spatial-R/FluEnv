
glm_rt_plot <- function(data,var.v){
  names(data)[which(names(data) == var.v)] <- "target"
  model.target <- glm(log(adjust.rc)~ns(target,3),data = data)
  predict <- predict(model.target, newdata = data, type = "terms",se.fit=T)
  res <- data.frame(predict$fit);res.se <- data.frame(predict$se.fit)
  estimation.reg <- data.frame(target=data$target,est=res$ns.target..3.,se=res.se$ns.target..3.)
  estimation.reg <- arrange(estimation.reg,target) %>% 
    mutate(lower = est - 1.96*se, upper = est + 1.96*se)
  estimation.reg[,c(2,4,5)] <- apply(estimation.reg[,c(2,4,5)],2,exp) 
  
  ggplot(data = estimation.reg,aes( x= target, y=est)) + geom_line(color = "red") +
    geom_ribbon(data=estimation.reg,aes(ymin=lower,ymax=upper),alpha=0.3,fill = "blue") +
    geom_hline(yintercept = 1,linetype =2) 
}

glm_rt_estimation <- function(data,var.v){
  names(data)[which(names(data) == var.v)] <- "target"
  model.target <- glm(log(adjust.rc)~ns(target,3),data = data)
  predict <- predict(model.target, newdata = data, type = "terms",se.fit=T)
  res <- data.frame(predict$fit);res.se <- data.frame(predict$se.fit)
  estimation.reg <- data.frame(target=data$target,est=res$ns.target..3.,se=res.se$ns.target..3.)
  estimation.reg <- arrange(estimation.reg,target) %>% 
    mutate(lower = est - 1.96*se, upper = est + 1.96*se)
  estimation.reg[,c(2,4,5)] <- apply(estimation.reg[,c(2,4,5)],2,exp) 
  return(estimation.reg)
}


cums_seres <- function(data,region.v,season.adj = TRUE){
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4)))
    dat.target.final <- data.frame()
    for (i in seq(from_year,to_year)){
      date.range <- seq.Date(as.Date(paste(i-1,"-",season.period[1],"-01",sep = "")),
                             as.Date(paste(i,"-",season.period[2],"-01",sep = "")),by = "day")
      dat.target <- filter(dat.region,date %in% date.range)
      target.date <- dat.target[which.max(dat.target$case),"date"]
      target.range <-  seq.Date(target.date - interval.v*7,target.date + interval.v*7,by = "day")
      dat.target <-    filter(dat.target,date %in% target.range)
      dat.target$cum <- c(0,cumsum(dat.target$case)[-nrow(dat.target)]);
      dat.target$num <- c(1:dim(dat.target)[1])
      dat.target$season <- i
      dat.target.final <- rbind(dat.target.final,dat.target)
    }
    return(dat.target.final)
  })
  date_result <- bind_rows(date_list)
  return(date_result)
}


Peak_Series <- function(data,region.v,interval.v,season.adj = TRUE,season.period){
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4)))
    dat.target.final <- data.frame()
    for (i in seq(from_year,to_year)){
      date.range <- seq.Date(as.Date(paste(i-1,"-",season.period[1],"-01",sep = "")),
                             as.Date(paste(i,"-",season.period[2],"-01",sep = "")),by = "day")
      dat.date <- filter(dat.region,date %in% date.range)
      target.date <- dat.date[which.max(dat.date$case),"date"]
      target.range <-  seq.Date(target.date - interval.v*7,target.date + interval.v*7,by = "day")
      if(!isTRUE(season.adj)){
        dat.date$cum <- c(0,cumsum(dat.date$case)[-dim(dat.region)]);
        dat.target <- filter(dat.date,date %in% target.range)
        dat.target$num <- c(1:nrow(dat.target))
      } else{
        dat.target <- filter(dat.region,date %in% target.range)
        dat.target$cum <- c(0,cumsum(dat.target$case)[-nrow(dat.target)]);
        dat.target$num <- c(1:nrow(dat.target))
      }
      dat.target$season <- i
      dat.target.final <- rbind(dat.target.final,dat.target)
    }
    return(dat.target.final)
  })
  date_result <- bind_rows(date_list)
  return(date_result)
}


epidemic_period_region <- function(data,region.v,interval.v,season.period){
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4)))
    dat.target.final <- data.frame()
    for (i in seq(from_year,to_year)){
      date.range <- seq.Date(as.Date(paste(i-1,"-",season.period[1],"-01",sep = "")),
                             as.Date(paste(i,"-",season.period[2],"-01",sep = "")),by = "day")
      dat.date <- filter(dat.region,date %in% date.range)
      target.date <- dat.date[which.max(dat.date$case),"date"]
      
      if(!region.id %in% north_china_cluster){
        dat.target <- dat.date
        dat.target$cum <- c(0,cumsum(dat.region$case)[-nrow(dat.region)]);
        dat.target$num <- c(1:dim(dat.target)[1])
        dat.target$season <- i
      } else {
        target.date <- dat.date[which.max(dat.date$case),"date"]
        target.range <-  seq.Date(target.date - interval.v*7,target.date + interval.v*7,by = "day")
        dat.target <-    filter(dat.date,date %in% target.range)
        dat.target$cum <- c(0,cumsum(dat.region$case)[-nrow(dat.region)]);
        dat.target$num <- c(1:dim(dat.target)[1])
        dat.target$season <- i
      }
    }
  })
}



extract.conituous.serie <- function(data,type="min",length,frist = T){
  
  data <- sort(data); 
  ori.date <- data.frame(date = data,id = rep(0,length(data)))
  date <- seq(data[1],data[length(data)])
  tar.date <- data.frame(date)
  date.tem <- merge(ori.date,tar.date,by = "date",all.y = T)
  date.tem$id <- ifelse(is.na(date.tem$id),1,date.tem$id)
  date.tem <- mutate(date.tem, cumgroup = cumsum(id));
  date.tem %>%
    group_by(cumgroup) %>%
    mutate(grouplength = n()) %>%
    ungroup() %>%
    filter(grouplength > 1) %>%
    data.frame() -> res.1 

  if(type == "min"){
    res.1.tem <- filter(res.1,grouplength > length)
      if(nrow(res.1.tem) == 0){
        res.2 <- NULL
      } else {
        if(isTRUE(frist)){
          res.2 <- filter(res.1.tem,grouplength == res.1.tem[1,"grouplength"])$date   
        } else {
          res.2 <- filter(res.1.tem,grouplength == res.1.tem[nrow(res.1.tem),"grouplength"])$date
        }
      }
    } else if (type == "max"){
      res.1.tem <- filter(res.1,grouplength < length)
        if(nrow(res.1.tem) == 0){
          res.2 <- NULL
    } else {
      if(isTRUE(frist)){
        res.2 <- filter(res.1.tem,grouplength == res.1.tem[1,"grouplength"])$date   
      } else {
        res.2 <- filter(res.1.tem,grouplength == res.1.tem[nrow(res.1.tem),"grouplength"])$date
      }
    }
  } else {
    res.2 <- res.1;
    res.1.tem <- res.1
  }
  return(list( pure = res.2, all = res.1.tem$date))
}



period_update <- function(year_v,month_v){
  year_v <- as.numeric(year_v);
  month_v <- as.numeric(month_v)
  if(month_v[1] <= 0){
    month_v[1] <- month_v[1] + 12
    year_v[1] <- year_v[1] -1
  }
  if(month_v[2] <= 0) {
    month_v[2] <- month_v[2] + 12
    year_v[2] <- year_v[2] -1
  }
  if(month_v[1] > 12){
    month_v[1] <- month_v[1] - 12
    year_v[1] <- year_v[1] + 1
  }
  if(month_v[2] > 12) {
    month_v[2] <- month_v[2] - 12
    year_v[2] <- year_v[2] + 1
  }
  return(list(year = year_v, month = month_v))
}

case_update <- function(origin_case,select_period){
  year.range <- select_period$year;
  month.range <- select_period$month
  date.range <- seq.Date(as.Date(paste(year.range[1],"-",month.range[1],"-01",sep = "")),
                         as.Date(paste(year.range[2],"-",month.range[2],"-01",sep = "")),by = "day")
  dat.fin <- filter(origin_case,date %in% date.range)
  return(dat.fin)
}


period_reshape <- function(target_data,full_data,month_period, year_period){
  
  dat.date <- target_data;dat.region <- full_data
  month_period_0 <- month_period;
  year_period_0 <- year_period;
  
  if(which.max(dat.date$cases) < (nrow(dat.date)/4)) {
    month_period_0[1] <- as.numeric(month_period_0[1]) - 2
  } 
  
  if(which.max(dat.date$cases) > (nrow(dat.date) - nrow(dat.date)/4)) {
    month_period_0[2] <- as.numeric(month_period_0[2]) + 2
  } 
  
  month_period_1 <- period_update(year_v = year_period_0,month_v = month_period_0)$month
  year_period_1  <- period_update(year_v = year_period_0,month_v = month_period_0)$year
  
  dat.date <- case_update(origin_case = dat.region,
                          select_period = list(year = year_period_1,month=month_period_1))
  
  ######################## second part
  
  incre_pos_0 <- which(diff(dat.date$cases) > 0); decrs_pos_0 <- which(diff(dat.date$cases) < 0)  ####
  if(length(incre_pos_0) == 0 ){
    incres_length <- NULL
  } else {
  incres_length <- extract.conituous.serie(incre_pos_0,type = "min",length = 20,frist = T)$pure[1]
  }
  
  if(is.null(incres_length)){
    first_point_incr <- FALSE
  } else if (incres_length < 10) {
    first_point_incr <- TRUE
  } else {
    first_point_incr <- FALSE 
  }
  
  
  if(length(decrs_pos_0) == 0 ){
    desc_length <- NULL
  } else {
    desc_length <- extract.conituous.serie(decrs_pos_0,type = "min",length = 20,frist = F)$pure[1]
  }
  
  if(is.null(desc_length)){
    end_point_decr <- FALSE
  } else if (desc_length > (nrow(dat.date) - 10) ) {
    end_point_decr <- TRUE
  } else {
    end_point_decr <- FALSE
  }
  
  if(first_point_incr){
    month_period_1[1] <- as.character(as.numeric(month_period_1[1]) - 1)
  }
  if(end_point_decr){
    month_period_1[2] <- as.character(as.numeric(month_period_1[2]) + 1)
  }
  month_period_2 <- period_update(year_v = year_period_1,month_v = month_period_1)$month
  year_period_2  <- period_update(year_v = year_period_1,month_v = month_period_1)$year
  
  dat.tem <- case_update(origin_case = dat.region,
                          select_period = list(year = year_period_2,month=month_period_2))
  return(list(data = dat.tem,month_period = month_period_2, year_period = year_period_2))
}



epidemic_period_peak <- function(data,region.v,season.period,interval.v){
  
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4))) 
    
    dat.target.final <- data.frame()
    
    for (year_id in seq(from_year,to_year)){
      
      print(paste(region.id,": year == ",year_id,sep=""))
      
      month_period_0 <- season.period
      year_period_0 <- c(year_id - 1, year_id)

      dat.date_0 <- case_update(origin_case = dat.region,
                                select_period = list(year = year_period_0,month=month_period_0))
      
      ######################### changed the period according to the peak time and consecutive 
      
      dat_result_0  <- period_reshape(target_data = dat.date_0,full_data = dat.region,
                                      month_period = month_period_0,year_period = year_period_0)
      month_period_1 <- dat_result_0$month_period; year_period_1 <- dat_result_0$year_period
      dat.date_1 <- dat_result_0$data
      
      begin.date <- as.Date(paste(year_period_1[1],"-",month_period_1[1],"-01",sep = ""))
      center.date <- begin.date + floor(nrow(dat.date_1)/2)
      
      #############################################################################################
      #############   remove the former epidemic with zero ########################################
      #############################################################################################
      
      low.date.range <- seq.Date(begin.date,center.date,by = "day")
      dat.low <- filter(dat.date_1,date %in% low.date.range)
      dat.high <- filter(dat.date_1,!date %in% low.date.range)
      
      incre_des <- which(diff(dat.low$cases) < 0)
      
      if(length(incre_des) == 0){
        dat.low <- NULL
      } else if(length(incre_des) > 14){
        continu.pos <- extract.conituous.serie(incre_des,type = "min",length = 20)$pure 
        if(!length(continu.pos) == 0){
          if(continu.pos[1] == 1){
            dat.low.bef <- filter(dat.region, 
                                  date %in% seq.Date((begin.date - 30),begin.date,by = "day"))
            pos.bef <- max(which(dat.low.bef$cases == min(dat.low.bef$cases)))
            dat.low.bef <- dat.low.bef[-c(1:pos.bef),]
            dat.low <- rbind(dat.low.bef,dat.low) 
          }
        }
      } else {
        dat.low <- dat.low
      }
      
      dat.date.1 <- rbind(dat.low,dat.high)
      zero_position <- which(dat.date.1$cases == 0)
      
      if(length(zero_position) == 0){
        dat.date.1 <- dat.date.1
      } else {
        zero_oimt <- extract.conituous.serie(zero_position,type = "min",length = 20,frist = F)$all
        if(length(zero_oimt) == 0){
          dat.date.1 <- dat.date.1
        } else {
          pos_omit_max <- max(zero_oimt); pos_omit_min <- min(zero_oimt)
          peak_pos <- which.max(dat.date.1$cases)
          if(peak_pos > pos_omit_max) dat.date.1 <- dat.date.1[-(1:pos_omit_max),]
          if(peak_pos < pos_omit_min) dat.date.1 <- dat.date.1[(1:pos_omit_min),]
        } #else {
          #pos_omit <- zero_oimt[which.max(diff(zero_oimt))]:zero_oimt[which.max(diff(zero_oimt)) + 1]
          #dat.date.1 <- dat.date.1[pos_omit,]
        #}
      }
      
      month_period_1 <- month(range(dat.date.1$date));year_period_1 <- year(range(dat.date.1$date))
      
      #############################################################################################
      #############   update the time-series again for the epidemic ###############################
      #############################################################################################
      
      if(dat.date.1[1,"cases"] == 0){
        dat.date_2 <- dat.date.1;
        year_period_2 <- year_period_1; month_period_2 <- month_period_1
      } else {
      dat_result_1  <- period_reshape(target_data = dat.date.1,full_data = dat.region,
                                      month_period = month_period_1,year_period = year_period_1)
      month_period_2 <- dat_result_1$month_period; year_period_2 <- dat_result_1$year_period
      dat.date_2 <- dat_result_1$data
      }
      
      begin.date <- as.Date(paste(year_period_2[1],"-",month_period_2[1],"-01",sep = ""))
      center.date <- begin.date + floor(nrow(dat.date_2)/2)
      
      ############################################################################################
      #####################################  step two: select the maximie value for ##############
      ############################################################################################
      
      ### step 1: for the start and end period to search for the epidemic period
      
      incres_series <- which(diff(dat.date_2$cases) > 0)   ### select for the value with maxmium value
      
      if(length(incres_series) == 0){                      ### no increase series, so the final result is null
        dat.final <- NULL
      } else {
        start.peak <- extract.conituous.serie(incres_series,type = "min",length = 20,frist =T)$pure   ## the start points
        if(is.null(start.peak)){                           ### no consecutive period above 3 weeks for the increases
          frist.peak.start.date <- dat.date_2[1,"date"]    ### the frist.peak.start.date is the first data points
          start.peak <- 1
        } else {
          frist.peak.start.date <- dat.date_2[start.peak[1],"date"] #### the first continus
        }
        
        ##### for the end points   
        end_peak_continue <- which(diff(dat.date_2$cases) < 0)
        
        if(length(end_peak_continue) == 0){
          end.peak <- NULL
        } else {
          end.peak <- extract.conituous.serie(end_peak_continue,type = "min",length = 20,frist =F)$pure
        }
        
        ############# select the time period which include the start and end points
        
        if(is.null(end.peak)){
          end.peak.date <- as.Date(paste(year_period_2[2],"-",month_period_2[2],"-01",sep = ""))
          dat.res.1 <- dat.date_2[c(start.peak[1]:nrow(dat.date_2)),] 
        } else {
          end.peak.date <- dat.date_2[rev(end.peak)[1],"date"]
          dat.res.1 <- dat.date_2[c(start.peak[1]:rev(end.peak)[1]),] 
        }
        
        
        ############ step 2: cut the time period into two parts
        
        #### situation 1: for the start date and center_month
        if(dat.res.1[1,"date"] <  center.date){
          peak.date.range <- seq.Date(dat.res.1[1,"date"],center.date,by = "day")        
        } else{
          peak.date.range <- dat.date_2$date
        }
        dat.peak.1 <- filter(dat.date_2,date %in% peak.date.range)
        
        #### situation 2: the end date and center month
        
        if(!(dat.peak.1[nrow(dat.peak.1),"date"] < center.date)){
          dat.peak.2 <- filter(dat.date_2,date >  center.date)
        } else {
          dat.peak.2 <- dat.date_2
        }
        
        peak.1.value <- max(dat.peak.1$cases);peak.2.value <- max(dat.peak.2$cases); 
        peak.1.date <- dat.peak.1[which.max(dat.peak.1$cases),"date"]
        peak.2.date <- dat.peak.2[which.max(dat.peak.2$cases),"date"]
        dat.peak.between <- filter(dat.date_2,date < peak.2.date & date > peak.1.date)
        diff.peak.between <- which(diff(dat.peak.between$cases) > 0)
        
        if(length(diff.peak.between) == 0){
          second.peak.up <- NULL
        } else {
          second.peak.up <- extract.conituous.serie(which(diff(dat.peak.between$cases) > 0),
                                                    type = "min",length = 20,frist =T)$pure
        }
        
        if(is.null(second.peak.up)){
          second.peak.start.date <- dat.date_2[nrow(dat.date_2),"date"]
          second.peak.dat <- NULL
        }else{  ### the peak time to the season end time ####
          second.peak.start.date <- dat.peak.between[second.peak.up[1],"date"]
          second.peak.gap <- as.numeric(diff.Date(c(peak.2.date,dat.date_2[nrow(dat.date_2),"date"]))) 
          if(second.peak.gap < 30 | second.peak.start.date > end.peak.date){
            second.peak.date.period <- seq.Date(second.peak.start.date,peak.2.date + 30,by = "day")
            second.peak.dat <- filter(dat.region,date %in% second.peak.date.period)
          } else {
            second.peak.dat <- filter(dat.date_2,date %in% seq.Date(second.peak.start.date,end.peak.date,by ="day"))
          }
        }
        
        frist.peak.down <- filter(dat.date_2,date %in% seq.Date(peak.1.date,second.peak.start.date,by ="day"))
        continus_peak_down_c <- which(diff(frist.peak.down$cases) < 0)
        
        if(nrow(frist.peak.down) > 20 & !length(continus_peak_down_c) == 0){
          frist.peak.down.con <- extract.conituous.serie(which(diff(frist.peak.down$cases) < 0),
                                                         type = "min",length = 20,frist =F)$pure
        } else {
          frist.peak.down.con <- NULL
        }
        
        if(is.null(frist.peak.down.con)){
          frist.peak.down.date <- frist.peak.down[nrow(frist.peak.down),"date"]
        } else {
          frist.peak.down.date <- frist.peak.down[rev(frist.peak.down.con)[1],"date"]
        }
        
        if(frist.peak.start.date < frist.peak.down.date){
           frist.peak.dat <- filter(dat.date_2,date %in% seq.Date(frist.peak.start.date,frist.peak.down.date,by = "day"))
           frist_second <- abs(as.numeric(diff(c(frist.peak.down.date,second.peak.dat[1,"date"]))))
        } else {
          frist.peak.dat <- NULL; frist_second <- NULL
        }
        
        if(length(frist_second) == 0){
          dat.final <- rbind(frist.peak.dat,second.peak.dat)
          if(!is.null(dat.final)){
          dat.final$season <- year_id
          dat.final$cum <- c(0,cumsum(dat.final$case)[-nrow(dat.final)]);
          dat.final <- winter_change(dat.final)
          dat.final <- dat.final[-c(1:7,(ncol(dat.final)-7):ncol(dat.final)),]
          
          max_date <- dat.final[which.max(dat.final$cases),"date"]
          min_dist <- min(c(as.numeric(abs(max_date - range(dat.final$date)[1])),
                            as.numeric(abs(max_date - range(dat.final$date)[2]))))
          if(min_dist < (0.1*nrow(dat.final))){
            dat.final <- NULL 
          }
          }
          # } else if (frist_second < 14){
          #   dat.final <- rbind(frist.peak.dat,second.peak.dat)
          #   dat.final$season <- year_id
          #   dat.final <- winter_change(dat.final)
          #   dat.final$cum <- c(0,cumsum(dat.final$case)[-nrow(dat.final)]);
          #   dat.final <- dat.final[-c(1:7,(ncol(dat.final)-7):ncol(dat.final)),]
        } else{
          
          if(!(is.null(frist.peak.dat)| nrow(frist.peak.dat) == 0)){
             frist.peak.dat$season <- year_id
             frist.peak.dat <- winter_change(frist.peak.dat)
             frist.peak.dat$cum <- c(0,cumsum(frist.peak.dat$case)[-nrow(frist.peak.dat)]);
             
             max_date <- frist.peak.dat[which.max(frist.peak.dat$cases),"date"]
             
             min_dist <- min(c(as.numeric(abs(max_date - range(frist.peak.dat$date)[1])),
                               as.numeric(abs(max_date - range(frist.peak.dat$date)[2]))))
             
             if(min_dist < (0.1*nrow(frist.peak.dat))){
               frist.peak.dat <- NULL 
             } else {
               frist.peak.dat <- frist.peak.dat[-c(1:7,(ncol(frist.peak.dat)-7):ncol(frist.peak.dat)),]
             }
          }
          
          if(!(is.null(second.peak.dat) | nrow(second.peak.dat) == 0)){
            second.peak.dat <- winter_change(second.peak.dat)
            second.peak.dat$season <- year_id
            second.peak.dat$cum <- c(0,cumsum(second.peak.dat$case)[-nrow(second.peak.dat)]);
            second.peak.dat <- second.peak.dat[-c(1:7,(ncol(second.peak.dat)-7):ncol(second.peak.dat)),]
            
            max_date <- second.peak.dat[which.max(second.peak.dat$cases),"date"]
            min_dist <- min(c(as.numeric(abs(max_date - range(second.peak.dat$date)[1])),
                              as.numeric(abs(max_date - range(second.peak.dat$date)[2]))))
            if(min_dist < (0.1*nrow(second.peak.dat))){
              second.peak.dat <- NULL 
            }
          }
          
          dat.final <- rbind(frist.peak.dat,second.peak.dat)
      }
      }
      #if(max(dat.final$cases) < 20) dat.final <- NULL
      dat.target.final <- rbind(dat.target.final,dat.final)
    }
    dat.target.final$season <- as.character(dat.target.final$season)
    return(dat.target.final)
  })
  result.final <- bind_rows(date_list)
  result.final <- filter(result.final,cases > 10)
  return(result.final)
}


winter_change <- function(data){
  peak_time <- data[which.max(data$cases),"date"]
  time_round <- as.numeric(julian(peak_time,
                                  origin = as.Date(paste(substr(peak_time,1,4),"-01-01",sep = ""))))/365.22
  ### from 4 to 9
  if(time_round > 0.24 & time_round < 0.67){
    data$winter <- 0
  } else {
    data$winter <- 1
  }
  return(data)
}




epidemic_period_peak_PNAS <- function(data,region.v,interval.v){
  
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    
    dat.region <- filter(dat.tem, region == region.id)
    
    max.ah <- dat.region[which.max(dat.region$cases),"date"]
    season.mid <- substr(max.ah - 60,6,7)
    center.date <- paste(stri_pad_left(season.mid,2,pad="0"),"-01",sep = "")
    season.start <- substr(max.ah - 180,6,7)
    season.end <- substr(max.ah + 180,6,7)
    
    if(season.start == season.end){
      season.start <- substr(max.ah - 210,6,7)
    }
    season.period <- c(season.start,season.end)
      
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4))) 
    
    dat.target.final <- data.frame()
    
    for (year_id in seq(from_year,to_year)){
      dat.final <- NULL
      tryCatch(expr = {
      print(paste(region.id,": year == ",year_id,sep=""))
      
      month_period_0 <- season.period
      year_period_0 <- c(year_id - 1, year_id)
      
      dat.date_0 <- case_update(origin_case = dat.region,
                                select_period = list(year = year_period_0,month=month_period_0))
      
      ######################### changed the period according to the peak time and consecutive 
      
      dat_result_0  <- period_reshape(target_data = dat.date_0,full_data = dat.region,
                                      month_period = month_period_0,year_period = year_period_0)
      month_period_1 <- dat_result_0$month_period; year_period_1 <- dat_result_0$year_period
      dat.date_1 <- dat_result_0$data
      
      begin.date <- as.Date(paste(year_period_1[1],"-",month_period_1[1],"-01",sep = ""))
      center.date <- begin.date + floor(nrow(dat.date_1)/2)
      
      #############################################################################################
      #############   remove the former epidemic with zero ########################################
      #############################################################################################
      
      low.date.range <- seq.Date(begin.date,center.date,by = "day")
      dat.low <- filter(dat.date_1,date %in% low.date.range)
      dat.high <- filter(dat.date_1,!date %in% low.date.range)
      
      incre_des <- which(diff(dat.low$cases) < 0)
      
      if(length(incre_des) == 0){
        dat.low <- NULL
      } else if(length(incre_des) > 14){
        continu.pos <- extract.conituous.serie(incre_des,type = "min",length = 20)$pure 
        if(!length(continu.pos) == 0){
          if(continu.pos[1] == 1){
            dat.low.bef <- filter(dat.region, 
                                  date %in% seq.Date((begin.date - 30),begin.date,by = "day"))
            pos.bef <- max(which(dat.low.bef$cases == min(dat.low.bef$cases)))
            dat.low.bef <- dat.low.bef[-c(1:pos.bef),]
            dat.low <- rbind(dat.low.bef,dat.low) 
          }
        }
      } else {
        dat.low <- dat.low
      }
      
      dat.date.1 <- rbind(dat.low,dat.high)
      zero_position <- which(dat.date.1$cases == 0)
      
      if(length(zero_position) == 0){
        dat.date.1 <- dat.date.1
      } else {
        zero_oimt <- extract.conituous.serie(zero_position,type = "min",length = 20,frist = F)$all
        if(length(zero_oimt) == 0){
          dat.date.1 <- dat.date.1
        } else if(all(diff(zero_oimt) == 1)){
          pos_omit_max <- max(zero_oimt); pos_omit_min <- min(zero_oimt)
          peak_pos <- which.max(dat.date.1$cases)
          if(peak_pos > pos_omit_max) dat.date.1 <- dat.date.1[-(1:pos_omit_max),]
          if(peak_pos < pos_omit_min) dat.date.1 <- dat.date.1[(1:pos_omit_min),]
        } else {
          pos_omit <- zero_oimt[which.max(diff(zero_oimt))]:zero_oimt[which.max(diff(zero_oimt)) + 1]
          dat.date.1 <- dat.date.1[pos_omit,]
        }
      }
      
      month_period_1 <- month(range(dat.date.1$date));year_period_1 <- year(range(dat.date.1$date))
      
      #############################################################################################
      #############   update the time-series again for the epidemic ###############################
      #############################################################################################
      
      if(dat.date.1[1,"cases"] == 0){
        dat.date_2 <- dat.date.1;
        year_period_2 <- year_period_1; month_period_2 <- month_period_1
      } else {
        dat_result_1  <- period_reshape(target_data = dat.date.1,full_data = dat.region,
                                        month_period = month_period_1,year_period = year_period_1)
        month_period_2 <- dat_result_1$month_period; year_period_2 <- dat_result_1$year_period
        dat.date_2 <- dat_result_1$data
      }
      begin.date <- as.Date(paste(year_period_2[1],"-",month_period_2[1],"-01",sep = ""))
      center.date <- begin.date + floor(nrow(dat.date_2)/2)
      
      ############################################################################################
      #####################################  step two: select the maximie value for ##############
      ############################################################################################
      
      ### step 1: for the start and end period to search for the epidemic period
      
      incres_series <- which(diff(dat.date_2$cases) > 0)   ### select for the value with maxmium value
      
      if(length(incres_series) == 0){                      ### no increase series, so the final result is null
        dat.final <- NULL
      } else {
        start.peak <- extract.conituous.serie(incres_series,type = "min",length = 20,frist =T)$pure   ## the start points
        if(is.null(start.peak)){                           ### no consecutive period above 3 weeks for the increases
          frist.peak.start.date <- dat.date_2[1,"date"]    ### the frist.peak.start.date is the first data points
          start.peak <- 1
        } else {
          frist.peak.start.date <- dat.date_2[start.peak[1],"date"] #### the first continus
        }
        
        ##### for the end points   
        end_peak_continue <- which(diff(dat.date_2$cases) < 0)
        
        if(length(end_peak_continue) == 0){
          end.peak <- NULL
        } else {
          end.peak <- extract.conituous.serie(end_peak_continue,type = "min",length = 20,frist =F)$pure
        }
        
        ############# select the time period which include the start and end points
        
        if(is.null(end.peak)){
          end.peak.date <- as.Date(paste(year_period_2[2],"-",month_period_2[2],"-01",sep = ""))
          dat.res.1 <- dat.date_2[c(start.peak[1]:nrow(dat.date_2)),] 
        } else {
          end.peak.date <- dat.date_2[rev(end.peak)[1],"date"]
          dat.res.1 <- dat.date_2[c(start.peak[1]:rev(end.peak)[1]),] 
        }
        
        
        ############ step 2: cut the time period into two parts
        
        #### situation 1: for the start date and center_month
        if(dat.res.1[1,"date"] <  center.date){
          peak.date.range <- seq.Date(dat.res.1[1,"date"],center.date,by = "day")        
        } else{
          peak.date.range <- dat.date_2$date
        }
        dat.peak.1 <- filter(dat.date_2,date %in% peak.date.range)
        
        #### situation 2: the end date and center month
        
        if(!(dat.peak.1[nrow(dat.peak.1),"date"] < center.date)){
          dat.peak.2 <- filter(dat.date_2,date >  center.date)
        } else {
          dat.peak.2 <- dat.date_2
        }
        
        peak.1.value <- max(dat.peak.1$cases);peak.2.value <- max(dat.peak.2$cases); 
        peak.1.date <- dat.peak.1[which.max(dat.peak.1$cases),"date"]
        peak.2.date <- dat.peak.2[which.max(dat.peak.2$cases),"date"]
        dat.peak.between <- filter(dat.date_2,date < peak.2.date & date > peak.1.date)
        diff.peak.between <- which(diff(dat.peak.between$cases) > 0)
        
        if(length(diff.peak.between) == 0){
          second.peak.up <- NULL
        } else {
          second.peak.up <- extract.conituous.serie(which(diff(dat.peak.between$cases) > 0),
                                                    type = "min",length = 20,frist =T)$pure
        }
        
        if(is.null(second.peak.up)){
          second.peak.start.date <- dat.date_2[nrow(dat.date_2),"date"]
          second.peak.dat <- NULL
        }else{  ### the peak time to the season end time ####
          second.peak.start.date <- dat.peak.between[second.peak.up[1],"date"]
          second.peak.gap <- as.numeric(diff.Date(c(peak.2.date,dat.date_2[nrow(dat.date_2),"date"]))) 
          if(second.peak.gap < 30 | second.peak.start.date > end.peak.date){
            second.peak.date.period <- seq.Date(second.peak.start.date,peak.2.date + 30,by = "day")
            second.peak.dat <- filter(dat.region,date %in% second.peak.date.period)
          } else {
            second.peak.dat <- filter(dat.date_2,date %in% seq.Date(second.peak.start.date,end.peak.date,by ="day"))
          }
        }
        
        frist.peak.down <- filter(dat.date_2,date %in% seq.Date(peak.1.date,second.peak.start.date,by ="day"))
        continus_peak_down_c <- which(diff(frist.peak.down$cases) < 0)
        
        if(nrow(frist.peak.down) > 20 & !length(continus_peak_down_c) == 0){
          frist.peak.down.con <- extract.conituous.serie(which(diff(frist.peak.down$cases) < 0),
                                                         type = "min",length = 20,frist =F)$pure
        } else {
          frist.peak.down.con <- NULL
        }
        
        if(is.null(frist.peak.down.con)){
          frist.peak.down.date <- frist.peak.down[nrow(frist.peak.down),"date"]
        } else {
          frist.peak.down.date <- frist.peak.down[rev(frist.peak.down.con)[1],"date"]
        }
        
        frist.peak.dat <- filter(dat.date_2,
                                 date %in% seq.Date(frist.peak.start.date,frist.peak.down.date,by = "day"))
        
        frist_second <- abs(as.numeric(diff(c(frist.peak.down.date,second.peak.dat[1,"date"]))))
        
        if(length(frist_second) == 0){
          dat.final <- rbind(frist.peak.dat,second.peak.dat)
          dat.final$season <- year_id
          dat.final$cum <- c(0,cumsum(dat.final$case)[-nrow(dat.final)]);
        } else if (frist_second < 20){
          dat.final <- rbind(frist.peak.dat,second.peak.dat)
          dat.final$season <- year_id
          dat.final$cum <- c(0,cumsum(dat.final$case)[-nrow(dat.final)]);
        } else{
          #frist.peak.dat$season <- paste(year_id,"-a",sep = "")
          #frist.peak.dat$cum <- c(0,cumsum(frist.peak.dat$case)[-nrow(frist.peak.dat)]);
          #second.peak.dat$season <- paste(year_id,"-b",sep = "")
          #second.peak.dat$cum <- c(0,cumsum(second.peak.dat$case)[-nrow(second.peak.dat)]);
          #dat.final <- rbind(frist.peak.dat,second.peak.dat)
          max.before <- max(frist.peak.dat$cases);max.after <- max(second.peak.dat$cases);
          if(max.before > max.after){
            dat.final <- frist.peak.dat
          } else {
            dat.final <- second.peak.dat
          }
          dat.final$season <- year_id
          dat.final$cum <- c(0,cumsum(dat.final$case)[-nrow(dat.final)]);
          
          #dat.final
        }}},    error = function(e){ 
          dat.final <- NULL
        })

      #if(max(dat.final$cases) < 20) dat.final <- NULL
      dat.target.final <- rbind(dat.target.final,dat.final)
    }
    dat.target.final$season <- as.character(dat.target.final$season)
    return(dat.target.final)
  })
  result.final <- bind_rows(date_list)
  return(result.final)
}





epidemic_continus_period <-function(data,region.v,quantilevalue,season.period){
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4)))
    dat.target.final <- data.frame()
    for (i in seq(from_year,to_year)){
      date.range <- seq.Date(as.Date(paste(i-1,"-",season.period[1],"-01",sep = "")),
                             as.Date(paste(i,"-",season.period[2],"-01",sep = "")),by = "day")
      
      dat.date <- filter(dat.region,date %in% date.range)
      
      start.peak <- extract.conituous.serie(which(diff(dat.date$cases) > 0),type = "min",length = 21)$pure
      
      if(is.null(start.peak)){
       peak.start.date <- dat.date[1,"date"]
       start.peak <- 1
      } else {
       peak.start.date <- dat.date[start.peak[1],"date"]
      }
      
      end.peak <- extract.conituous.serie(which(diff(dat.date$cases) < 0),type = "min",length = 21,frist =F)$pure
      
      if(is.null(end.peak)){
        end.peak.date <- as.Date(paste(i,"-",season.period[2],"-01",sep = ""))
        period.range <- dat.date[(start.peak[1]:nrow(dat.date)),"date"]
      } else {
        end.peak.date <- dat.date[rev(end.peak)[1],"date"]
        period.range <- dat.date[(start.peak[1]:rev(end.peak)[1]),"date"]
      }
      dat.target <- mutate(dat.date,period =ifelse(date %in% period.range,1,0)) 
      #cutvalue <- as.numeric(quantile(dat.res.1$cases,at = quantilevalue))
      #period.start <- which(dat.res.1$cases > cutvalue)[1]
      #period.end <- which(dat.res.1[-c(1:period.start),"cases"] < cutvalue)[1]
      
      #dat.target <- dat.res.1[period.start:period.end,];
      #dat.target$season <- i
      dat.target.final <- rbind(dat.target.final,dat.target)
    }
    return(dat.target.final)
  })
  date_result <- bind_rows(date_list)
  return(date_result)
}




epidemic_continus_period_date <-function(data,region.v,quantilevalue,season.period){
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.v)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(as.numeric(substr(dat.region$date,1,4))) + 1
    to_year <-   max(as.numeric(substr(dat.region$date,1,4)))
    dat.target.final <- data.frame()
    for (i in seq(from_year,to_year)){
      date.range <- seq.Date(as.Date(paste(i-1,"-",season.period[1],"-01",sep = "")),
                             as.Date(paste(i,"-",season.period[2],"-01",sep = "")),by = "day")
      
      dat.date <- filter(dat.region,date %in% date.range)
      
      start.peak <- extract.conituous.serie(which(diff(dat.date$cases) > 0),type = "min",length = 21)
      
      if(is.null(start.peak)){
        peak.start.date <- dat.date[1,"date"]
        start.peak <- 1
      } else {
        peak.start.date <- dat.date[start.peak[1],"date"]
      }
      
      end.peak <- extract.conituous.serie(which(diff(dat.date$cases) < 0),type = "min",length = 21,frist =F)
      
      if(is.null(end.peak)){
        end.peak.date <- as.Date(paste(i,"-",season.period[2],"-01",sep = ""))
        period.range <- dat.date[(start.peak[1]:nrow(dat.date)),"date"]
      } else {
        end.peak.date <- dat.date[rev(end.peak)[1],"date"]
        period.range <- dat.date[(start.peak[1]:rev(end.peak)[1]),"date"]
      }
      dat.target <- mutate(dat.date,period =ifelse(date %in% period.range,1,0)) 
      dat.target <- data.frame(start = period.range[1],end = rev(period.range)[1],
                               season = i,region = region.id)
      dat.target.final <- rbind(dat.target.final,dat.target)
    }
    return(dat.target.final)
  })
  date_result <- bind_rows(date_list)
  return(date_result)
}


epidemic_period <- function(data,method = "varying",region.v = region.v,
                            season.period = c("09","08"),interval.v = 7){
  dat.tem <- data
  if( method ==  "fixed"){
    dat.1 <- Peak_Series(data = dat.tem,region.v = region.v,
                                    season.period = season.period,season.adj = T,
                                    interval.v = interval.v)
  } else {
    dat.1 <- epidemic_period_peak(data = dat.tem,region.v = region.v,
                                   season.period = season.period,
                                   interval.v = interval.v)
  }
  return(dat.1)
} 



epidemic_period_interaction <- function(data,region.v = region.v,
                            season.period = c("09","08"),period.range = 365){
  dat.tem <- data
  dat.1 <- epidemic_period_peak(data = dat.tem,region.v = region.v,
                                  season.period = season.period,
                                  interval.v = interval.v)
  start_date <- range(dat.1$date)[1]+period.range
  dat.2 <- filter(dat.1,date > start_date)

  region_list <- lapply(unique(dat.1$region),function(region_id){
    dat.region <- filter(dat.1,region == region_id)
    date_list <- lapply(dat.region$date,function(date_id){
      date_range <- seq.Date(date_id - period.range,date_id,by = "day")
      dat_sum <- filter(dat.tem,region == region_id & date %in% date_range)
      h1_cases <- sum(dat_sum$H1)
      hb_cases <- sum(dat_sum$HB)
      data.frame(date = date_id,h1 = h1_cases,hb = hb_cases)
    })
   date_res <- bind_rows(date_list)
   dat.region.1 <- merge(dat.region,date_res,by = "date")
   return((dat.region.1))
})
  region_res <- bind_rows(region_list)
  return(region_res)
} 



epidemic_period_PNAS <- function(data,method = "fixed",region.v = region.v,
                            interval.v = 7){
  dat.tem <- data
  if( method ==  "fixed"){
    dat.1 <- Peak_Series(data = dat.tem,region.v = region.v,
                         season.period = season.period,season.adj = T,
                         interval.v = interval.v)
  } else {
    dat.1 <- epidemic_period_peak_PNAS(data = dat.tem,region.v = region.v,
                                  #
                                  interval.v = interval.v)
  }
  return(dat.1)
} 


adjust_rt_rr <- function(data,intercept = T,holiday = T){
  
  dat.1 <- data
  
  result_list <- lapply(unique(dat.1$region), function(region.id){
    print(region.id)
    dat.region <- filter(dat.1, region == region.id)[-1,]
    dat.region$time <- dat.region$num
    dat.region <- mutate(dat.region, id = paste(season,winter,sep="-"))
    adjust_list <- lapply(unique(dat.region$id),function(year.id){
      print(year.id)
      dat.region.season <- filter(dat.region,id == year.id) 
       dat.region.target <- dat.region.season
      dat.region.target <- filter(dat.region.target,!cases == 0)
      dat.region.target$time <-   julian.Date(dat.region.target$date,origin = as.Date("2010-09-01"))/365.25+2010
    
      if(nrow(dat.region.target) < 21){
        return(NULL)
      } else {
      ################################# model performance #############################################
      
      if(isTRUE(holiday)){
        dat.region.target <- merge(dat.region.target,vaccnation,by.x ="date",by.y ="date",all.x=T)
          model.int <- summary(lm(log(rc)~cum + factor(work),data = dat.region.target))
          model.reint <- summary(lm(log(rc)~cum + factor(work)-1,data = dat.region.target))
      } else {
          model.int <- summary(lm(log(rc)~cum ,data = dat.region.target))
          model.reint <- summary(lm(log(rc)~cum -1,data = dat.region.target))
      }
      data.frame(region = region.id, 
                 year = year.id,
                 r2_int = model.int$r.squared,
                 r2_reint = model.reint$r.squared)  
      }})
    adjust_rc_resaon <- bind_rows(adjust_list)
    return(adjust_rc_resaon) 
  })
  adjust_rc_resaon_final <- bind_rows(result_list)
  return(adjust_rc_resaon_final)
}


env_relate_rt_glm <- function(data,type = ""){
  
  dat.1 <- data
  
  result_list <- lapply(unique(dat.1$region), function(region.id){
    print(region.id)
    dat.region <- filter(dat.1, region == region.id)
    dat.region.target <- merge(dat.region,vaccnation,by.x ="date",by.y ="date",all.x=T)
    sea_length <- length(unique(dat.region.target$season))
    
    if(type == "env"){
        ################################# model performance #############################################
      
        if(sea_length > 1){
          model.all <- glm(log(rc)~ns(MT,3)+ns(RH,3) + factor(season) + 
                           factor(work) + cum,data = dat.region.target)
        } else {
          model.all <- glm(log(rc)~ns(MT,3)+ns(RH,3) + 
                             factor(work) + cum,data = dat.region.target)
        }
        predict <- predict(model.all, newdata = dat.region.target, type = "terms",se.fit=T)
        res <- data.frame(predict$fit);res.se <- data.frame(predict$se.fit)
        estimation.reg <- data.frame(est.rh=res$ns.RH..3.,se.rh=res.se$ns.RH..3.,
                                     est.mt=res$ns.MT..3.,se.mt=res.se$ns.MT..3.)
        estimation.reg <- cbind(estimation.reg,dat.region.target)
        estimation.reg <- arrange(estimation.reg,RH) %>% 
          mutate(lower.rh = est.rh - 1.96*se.rh, upper.rh = est.rh + 1.96*se.rh,
                 lower.mt = est.mt - 1.96*se.mt, upper.mt = est.mt + 1.96*se.mt)
        estimation.reg[,c(2,5,7:10)] <- apply(estimation.reg[,c(2,5,7:10)],2,exp) 
    } else if (type == "scs") {
      
      date_num <- as.numeric(julian(dat.region.target$date,origin = as.Date("2009-01-01"))/365)
      dat.region.target$date_num_update <- date_num #- floor(date_num)
      
      if(sea_length > 1){
        model.all <- glm(log(rc)~ sin(2*pi*date_num_update) + cos(2*pi*date_num_update) +
                           sin(pi*date_num_update) + cos(pi*date_num_update)+
                         factor(work) + cum + factor(season),data = dat.region.target)
      } else {
        model.all <- glm(log(rc)~sin(2*pi*date_num_update) + cos(2*pi*date_num_update) +
                         sin(pi*date_num_update) + cos(pi*date_num_update)+
                         factor(work) + cum,data = dat.region.target)
      }
      predict <- predict(model.all, newdata = dat.region.target, type = "terms",se.fit=T)
      res <- data.frame(predict$fit);res.se <- data.frame(predict$se.fit)
      names(res)[c(1:4)] <- c("sin_o","cos_o","sin_1o","cos_1o")
      estimation.reg <- data.frame(adjust.rc=exp(res$sin_o+res$cos_o+res$sin_1o+res$cos_1o))
      estimation.reg <- cbind(estimation.reg,dat.region.target)
    } else {
      if(sea_length > 1){
        model.all <- glm(log(rc)~factor(work) + cum + factor(season),data = dat.region.target)
      } else {
        model.all <- glm(log(rc)~factor(work) + cum, data = dat.region.target)
      }
      estimation.reg <- data.frame(adjust.rc=exp(model.all$residuals))
      estimation.reg <- cbind(estimation.reg,dat.region.target)
    } 
      return(estimation.reg)
      })  
  adjust_rc_resaon_final <- bind_rows(result_list)
  return(adjust_rc_resaon_final)
}



env_relate_rt_antigenic <- function(data,holiday=T,intercept = T){
  
  dat.1 <- data

  result_list <- lapply(unique(dat.1$region), function(region.id){
    print(region.id)
    dat.region <- filter(dat.1, region == region.id)
    adjust_list <- lapply(unique(dat.region$season),function(year.id){
      print(year.id)
      dat.region.season <- filter(dat.region,season == year.id & !is.na(MT)) 
      dat.region.target <- dat.region.season
      dat.region.target <- filter(dat.region.target,!cases == 0)
      dat.region.target$time <-   julian.Date(dat.region.target$date,origin = as.Date("2010-09-01"))/365.25+2010
      dat.region.target <- merge(dat.region.target,vaccnation,by.x ="date",by.y ="date",all.x=T)
      
      if(region.id %in% south_china_cluster){
        dat.anti <- filter(antigenci,area == "South")
      } else {
        dat.anti <- filter(antigenci,area == "South")
      }
      anti_more <- plyr::summarize(dat.anti,
                                   time = seq(from = 2011,to = 2017,by = 1/(365)),
                                   dist0 = predict(smooth.spline(x = month,y = dist0),x = time)$y,
                                   dist1 = predict(smooth.spline(x = month,y = dist1),x = time)$y)
      start_num <- which.min(abs(dat.region.target[1,"time"] - anti_more$time))
      anti_target <- anti_more[(start_num:(start_num + nrow(dat.region.target)-1)),]
      dat.region.target <- cbind(dat.region.target,anti_target[,-1])
      
      if(nrow(dat.region.target) < 21){
        return(NULL)
      } else {
        ################################# model performance #############################################
        
        if(isTRUE(holiday)){

          model.lm <- summary(lm(log(rc)~cum + factor(work) + dist1 + dist0,data = dat.region.target))
        } else {
          model.lm <- summary(lm(log(rc)~cum + dist1 + dist0, data = dat.region.target))
        }
        
        if(isTRUE(intercept)){
          dat.region.target$adjust.rc <- exp(model.lm$residuals + model.lm$coefficients[1])        
        } else {
          dat.region.target$adjust.rc <- exp(model.lm$residuals)        
        }
          return(dat.region.target[-1,])
      }
    })
    adjust_rc_resaon <- bind_rows(adjust_list)
    return(adjust_rc_resaon) 
  })
  adjust_rc_resaon_final <- bind_rows(result_list)
  return(adjust_rc_resaon_final)
}


lm_eqn = function(m) {
  # Displays regression line equation and R^2 value on plot
  # Usage:
  # p + annotate("text", x=25, y=300, label=lm_eqn(lm(y ~ x, df)), parse=TRUE)
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
}



epidemic_period_vis <- function(data=lh.result,epidemic_date = epidemic.period$date){
  dat.tem <- data
  dat.tem$period <- ifelse(dat.tem$date %in% epidemic_date,1,0)
  
  ggplot(data = dat.tem, aes(x = date, y = cases)) +
    geom_point(aes(color = factor(period))) +
    scale_colour_manual(values = c("black","blue"),guide = "none") +
    theme_bw(base_family = "serif",base_size = 12) 
}


# 
# peak3<-c(100,125,145,199)
# 
# lpeak3<-length(peak3)
# r.square<-matrix(0,nrow=14,ncol =lpeak3)
# for(r in 1:(lpeak3)){
#   for(h in 2:15){
#     ny1<-mydata$count.smooth[(peak3[r] - h):(peak3[r] + h)]
#     y1<-mydata$count.smooth[(peak3[r] - h):(peak3[r] + h)]
#     s<-sum(y1)
#     lambda<-2200
#     alpha<-0.6
#     x1<-seq(1,length(y1),1)
#     for(ii in 2:length(y1)){
#       y1[ii]<-y1[ii]+y1[ii-1]
#     }
#     
#     
#     logistic<-function(t, lambda, alpha,s){
#       s/((1+lambda*exp(-alpha*t)))
#     }
#     
#     R.square<-function(para, x, y)
#     {
#       lambda<-exp(para[1])
#       alpha<-exp(para[2])
#       s<-exp(para[3])
#       y.hat<-logistic(x, lambda,alpha,s)
#       y.hat1<-y.hat
#       for(i in 2:(length(y.hat1))){
#         y.hat1[i]<-y.hat[i]-y.hat[(i-1)]
#       }
#       for(i in 1:length(y.hat1)){
#         cc<-rep(0,length(y.hat1))
#         dd<-rep(0,length(y.hat1))
#         for(i in 1:length(y.hat1)){
#           cc[i]<-((y.hat1[i]-ny1[i])*(y.hat1[i]-ny1[i]))
#           dd[i]<-((ny1[i]-mean(ny1))*(ny1[i]-mean(ny1)))
#         }
#         aaa<-1-sum(cc)/sum(dd)
#       }
#       aaa
#     }
#     
#     ini.par<-log(c(lambda, alpha,s))
#     fit<-optim(ini.par, R.square, x=x1, y=y1, method='BFGS', control=list(fnscale=-1, reltol=1e-5))
#     y.fit<-logistic(x1, exp(fit$par[1]), exp(fit$par[2]), exp(fit$par[3])) 
#     y.fit1<-y.fit
#     for(i in 2:(length(y.fit1))){
#       y.fit1[i]<-y.fit[i]-y.fit[(i-1)]
#     }
#     r.square[(h-1),r]<-fit$value
#     plot(x1,y.fit1,type="l",col="black")
#     lines(x1,ny1,type="l",col="red")  
#   }
# }

glm_model_estimation <- function(data,vars = c("AH","MT"),degree=c(3,3)){
  
  if(!length(vars) == length(degree)){
    stop("The number of variables are different from the degree settings")
  }
  
  var_col <- which(names(data) %in% vars);target_col <- which(names(data) == "adjust.rc")
  dat.tem <- data[,c(target_col,var_col)]

  formu_1 <- unlist(lapply(1:length(vars),function(id){
    paste("ns(",vars[id],",",degree[id],")",sep = "")
  }))
  formu_2 <- as.formula(paste("log(adjust.rc)~",paste(formu_1,collapse = "+"),sep = ""))
  model.tem <- glm(formu_2,data = dat.tem)
  predict <- predict(model.tem, newdata = dat.tem, type = "terms",se.fit=T)
  res <- data.frame(predict$fit);res.se <- data.frame(predict$se.fit)
  
  res_list <- lapply(1:length(vars), function(id){
    est.tem <- res[,id];se.tem <- res.se[,id];
    res.tem <- data.frame(vars_v = dat.tem[,id+1],est = est.tem,
               low = est.tem - 1.96*se.tem,high =  est.tem + 1.96*se.tem)
    names(res.tem) <- paste(vars[id],c("","_est","_low","_high"),sep = "")
    res.tem[,2:4] <- apply(res.tem[,2:4],2,exp)
    return(res.tem)
  })
  res_final <- bind_cols(res_list)
  res_final <- apply(res_final,2,round,3)
  return(data.frame(res_final))
}



env_relate_rt <- function(data,intercept = T,holiday = T,Rsquare = T,spline_set){
  
  dat.1 <- data
  
  #dat.1 <- Peak_Series(data = dat.tem,region.v = "region",interval.v = 15,season.adj = T,
  #                     season.period = c("09","08"))
  #
  result_list <- lapply(unique(dat.1$region), function(region.id){
    print(region.id)
    dat.region <- filter(dat.1, region == region.id)[-1,]
    dat.region$time <- dat.region$num
    #  p.dat <- pdata.frame(dat.region, index=c("season","num"), drop.index=TRUE, row.names=TRUE)
    #  model.plm <- plm(log(rc)~cum + cos(2*pi*time/365) + sin(2*pi*time/365) +
    #                     cos(4*pi*time/365) + sin(4*pi*time/365),data = p.dat,model = "within",effect = "individual")
    #  coef <- summary(model.plm)$coefficients
    dat.region$id <- paste(dat.region$season,dat.region$winter,sep = "")
    adjust_list <- lapply(unique(dat.region$id),function(year.id){
      print(year.id)
      dat.region.season <- filter(dat.region,id == year.id) 
      dat.region.target <- dat.region.season
      dat.region.target <- filter(dat.region.target,!cases == 0)
      dat.region.target$time <-   julian.Date(dat.region.target$date,origin = as.Date("2010-09-01"))/365.25+2010
      exbspline_basis <- periodic.bspline.basis(dat.region.target$time,nbasis=3,degree=2,period=1)
      colnames(exbspline_basis)<- paste("xi",1:3,sep="")
      dat.region.target <- cbind(dat.region.target,exbspline_basis)
      
      if(nrow(dat.region.target) < 21){
        return(NULL)
      } else {

        if(isTRUE(holiday)){
          dat.region.target <- merge(dat.region.target,vaccnation,by.x ="date",by.y ="date",all.x=T)
          if(isTRUE(spline_set)){
            model.lm <- summary(lm(log(rc)~cum + factor(work)+xi1+xi2+xi3,data = dat.region.target))
          } else {
            model.lm <- summary(lm(log(rc)~cum + factor(work),data = dat.region.target))
          }
        } else {
          if(isTRUE(spline_set)){
            model.lm <- summary(lm(log(rc)~cum+xi1+xi2+xi3,data = dat.region.target))
          } else {
            model.lm <- summary(lm(log(rc)~cum,data = dat.region.target))
          }
        }
        
        if(isTRUE(intercept)){
          dat.region.target$adjust.rc <- exp(model.lm$residuals + model.lm$coefficients[1])        
        } else {
          dat.region.target$adjust.rc <- exp(model.lm$residuals)        
        }
        if(isTRUE(spline_set)){
          adjust.tem <- model.lm$coefficients[4]*dat.region.target$xi1 + model.lm$coefficients[5]*dat.region.target$xi2 +
            model.lm$coefficients[6]*dat.region.target$xi3# + model.lm$coefficients[7]*dat.region.target$xi4
          dat.region.target$adjust.rc <- exp(adjust.tem)
        } 
        
        if(isTRUE(Rsquare) & (model.lm$r.squared > 0.3)){
          return(dat.region.target[-1,])
        } else if (!isTRUE(Rsquare)){
          return(dat.region.target[-1,])
        } else {
          return(NULL)
        }
      }
    })
    adjust_rc_resaon <- bind_rows(adjust_list)
    return(adjust_rc_resaon) 
  })
  adjust_rc_resaon_final <- bind_rows(result_list)
  adjust_rc_resaon_final$id <- NULL
  return(adjust_rc_resaon_final)
}


###' @title Convert the weekly case into daily basis using the splines method
###' @description Convert the weekly surveillance dataset into daily basis. In the most of public surveillance system, the
###' repoted case are all aggregated into weekly basis, however, it is difficult to estimate the daily reproduction number
###' directly from the weekly dataset. Whereas, we can interpreate the weekly case into daily basis using the spline method.
###' @param data data.frame includes the columns for the case, date, region, and population.
###' @param date.v character, column name for the date
###' @param case.v character, column name for the case
###' @param region.v character, column name for the region
###' @param pop.v character,column name for the population
###' @param continue logistcal, choose the whole observtion period to smooth
###' @examples
###' \dontrun{
###' Weekly_To_Daily(
###' data = influenza_fs,
###' date.v = "date",
###' case.v = "count",
###' region.v ="region",
###' pop.v ="population",
###' continue = TRUE)
###' }
###' @return a dataframe.
###' @author Bing Zhang (Spatial-R), \url{https://github.com/Spatial-R/EFTR}


Weekly_To_Daily <- function(data,    date.v="date",
                            case.v ="count",
                            region.v = "region",
                            pop.v ="population",
                            continue = TRUE){
  dat.tem <- data
  var.pos <- unlist(lapply(c(date.v,case.v,region.v,pop.v),function(data)which(names(dat.tem) == data)))
  dat.tem <- dat.tem[,var.pos];names(dat.tem) <- c("date","count","region","pop")
  dat.tem$year <- as.numeric(substr(dat.tem$date,1,4))
  
  if(isTRUE(continue)){
    
    daily.incidence.list <- lapply(unique(dat.tem$region),function(region.id){
      dat.region <- filter(dat.tem,region == region.id)
      dat.region <- mutate(dat.region,cum = cumsum(count),incid = cum/pop,number = as.numeric(date))
      date.range <- seq.Date(dat.region[1,"date"],dat.region[nrow(dat.region),"date"],by="day")
      incid.tem <- plyr::summarize(dat.region,
                                   week = seq(from = min(number),
                                              to = max(number),by = 1),
                                   incid = predict(smooth.spline(x = number,y = incid),x = week)$y,
                                   popul = predict(smooth.spline(x = number,y = pop),x = week)$y)
      result.inc <- diff(incid.tem[,2])
      result.inc <- ifelse(result.inc <0,0,result.inc)
      result.fin <- data.frame(date=date.range[-1],case = floor(result.inc * incid.tem[-1,3]),
                               population = incid.tem[-1,3])
      result.fin$region <- region.id
      return(result.fin)
    })
  } else {
    
    daily.incidence.list <- lapply(unique(dat.tem$region),function(region.id){
      
      dat.region <- filter(dat.tem,region == region.id)
      year.list <- lapply(unique(dat.region$year),function(year.id){
        dat.year <- filter(dat.region,year == year.id)
        dat.year <- mutate(dat.year,cum = cumsum(count),incid = cum/pop,number = as.numeric(date))
        
        min.number <- as.numeric(dat.year[1,1]-as.Date(paste(year.id,"-01-01",sep="")))
        max.number <- as.numeric(as.Date(paste(year.id,"-12-31",sep="")) - dat.year[nrow(dat.year),1])
        
        incid.tem <- plyr::summarize(dat.year,
                                     week = seq(from = min(number)-min.number-1,
                                                to = max(number)+max.number,by = 1),
                                     incid = predict(smooth.spline(x = number,y = incid),x = week)$y,
                                     popul = predict(smooth.spline(x = number,y = pop),x = week)$y)
        result.inc <- diff(incid.tem[,2])
        result.inc <- ifelse(result.inc <0,0,result.inc)
        result.date <- seq.Date(as.Date(paste(year.id,"-01-01",sep="")),
                                as.Date(paste(year.id,"-12-31",sep="")),by="day")
        result.fin <- data.frame(date=result.date,case = floor(result.inc * incid.tem[-1,3]),
                                 population = incid.tem[-1,3])
        
        return(result.fin)
      })
      dat.daily <- bind_rows(year.list)
      dat.daily$region <- region.id
      return(dat.daily)
      print(region.id)
    })}
  daily.incidence <- bind_rows(daily.incidence.list)
  return(daily.incidence)
}



###' @title Estimate the daily instantaneous reproduction number ($R_{t}$)
###' @description
###' Transmissibility of infectious disease can be measured by the instantaneous
###' reproduction number (Rt), the average number of secondary cases that each infected
###' individual would infect if the conditions remained as they were at time t.
###' we utilized the likelihood-based approach proposed in
###' "Estimating Individual and Household Reproduction Numbers in an Emerging Epidemic"
###' which modeled the transmission with a Poisson process.
###' @param data daily case dataset including at least two columns for the date and count
###' @param case.v character, column name for the case
###' @param pop.v character,column name for the population
###' @param n.num measure windows to analysis, seeing EpiEstim package
###' @param sid serial interval distribution
###' @examples
###' \dontrun{
###' Weekly_To_Daily(
###' data = influenza_fs,
###' date.v = "date",
###' case.v = "count",
###' region.v ="region",
###' pop.v ="population",
###' continue = TRUE)
###' }
###' @return a dataframe.
###' @author Bing Zhang (Spatial-R), \url{https://github.com/Spatial-R/EFTR}

Rc <- function(data,date.id,case.v,n.num,DicreteSIDistr){
  dat.1 <- data
  col.number <- which(names(dat.1) %in% c(date.id,case.v))
  dat.1 <- dat.1[,c(col.number)]
  names(dat.1) <- c("date","cases")
  #n.first <- ifelse(time.interval == "week",1, which(dat.1[,"date"] == start.date))
  #dat.1 <- dat.1[-c(1:(n.first - 10)),]
  n.total <- nrow(dat.1)
  result <- EstimateR(dat.1[,"cases"],T.Start = 2:(n.total - n.num),
                      T.End = (2 + n.num):n.total, method = "NonParametricSI",
                      SI.Distr=DicreteSIDistr, plot=F)
  Rc.R.1 <- result$R;
  Rc.R.1[,"date"] <- dat.1[c(2:(n.total - n.num)),"date"]
  Rc.R.1[,"cases"] <- dat.1[c(2:(n.total - n.num)),"cases"]
  return(Rc.R.1)
}


Dicret_SI <- function(mu,sigma,day){
  MeanFluSI <- mu;sigma = sigma
  DicreteSIDistr <- vector()
  for(i in 0:day)
  {
    DicreteSIDistr[i+1] <- DiscrSI(i, MeanFluSI, SdFluSI)
  }
  return(DicreteSIDistr)
}



epidemic_calculate <- function(data){
  dat_tem <- data
  dat_tem <- mutate(dat_tem,id = paste(region,season,winter,sep = "-"))
  epidemic_list <- lapply(unique(dat_tem$id),function(id_tem){
    dat_tem_1 <- filter(dat_tem, id == id_tem)

    epidemic_dur <- range(dat_tem_1$date)
    dat_res <- data.frame(region = dat_tem_1[1,"region"],season = dat_tem_1[1,"season"],
               winter = dat_tem_1[1,"winter"], id = id_tem,
               start_date = epidemic_dur[1],end_date = epidemic_dur[2],
               peak_date = dat_tem_1[which.max(dat_tem_1$cases),"date"])
    dat_res1 <- mutate(dat_res,start_day = as.numeric(peak_date - start_date)/7,
                           end_day = as.numeric(end_date - peak_date)/7,
                           peak_week_1 = 52 + as.numeric(peak_date - as.Date(paste(season,"-01-01",sep="")))/7)
    
    
    if(unique(dat_res1$winter) == 1){
      dat_res1  <- filter(dat_res1,peak_week_1 < (52+0.67*365/7))
    }
    return(dat_res1)
  })
  epidemic_dat <- bind_rows(epidemic_list)

  # epidemic_dat <- mutate(epidemic_dat,start_day = as.numeric(peak_date - start_date)/7,
  #                        end_day = as.numeric(end_date - peak_date)/7,
  #                        peak_week = week(peak_date),
  #                        peak_week_1 = ifelse(winter == 1 & peak_week < 26, (52 + peak_week),peak_week),
  #                        peak_week_1 = ifelse(winter == 2, (52 + peak_week_1),peak_week_1))
  
  return(epidemic_dat)
} 



env_relate_rt_coef <- function(data,intercept = T,holiday = T,Rsquare = T,spline_set){
  
  dat.1 <- data
  
  result_list <- lapply(unique(dat.1$region), function(region.id){
    print(region.id)
    dat.region <- filter(dat.1, region == region.id)[-1,]
    dat.region$time <- dat.region$num
    #  p.dat <- pdata.frame(dat.region, index=c("season","num"), drop.index=TRUE, row.names=TRUE)
    #  model.plm <- plm(log(rc)~cum + cos(2*pi*time/365) + sin(2*pi*time/365) +
    #                     cos(4*pi*time/365) + sin(4*pi*time/365),data = p.dat,model = "within",effect = "individual")
    #  coef <- summary(model.plm)$coefficients
    dat.region <- mutate(dat.region, id = paste(season,winter,sep="-"))
    adjust_list <- lapply(unique(dat.region$id),function(year.id){
      print(year.id)
      dat.region.season <- filter(dat.region,id == year.id) 
      #  target.date <- dat.region.season[which.max(dat.region.season$case),"date"]
      #  start.date <- which(dat.region.season$cases > 10)[1]
      #  interval.v <- 10
      #  target.range <-  seq.Date(target.date - interval.v*7,target.date + interval.v*7,by = "day")
      #  if(region.id %in% north_china){
      #    dat.region.season <- dat.region.season[dat.region.season$date %in% target.range,]
      #  }
      # zero.date <- which(dat.region.season$rc == change.value)
      #  if(length(zero.date) < 4 | all(diff(zero.date) == 1)){
      dat.region.target <- dat.region.season
      #  } else {
      #    discou.point <- which.max(diff(which(dat.region.season$rc == change.value)))
      #    date.range <- zero.date[c(discou.point,discou.point+1)]
      #    dat.region.target <- dat.region.season[date.range[1]:date.range[2],]        
      #  }
      #dat.region.target$cum <- cumsum(dat.region.target$cases)
      dat.region.target <- filter(dat.region.target,!cases == 0)
      
      dat.region.target$time <-   julian.Date(dat.region.target$date,origin = as.Date("2010-09-01"))/365.25+2010
      exbspline_basis <- periodic.bspline.basis(dat.region.target$time,nbasis=3,degree=2,period=1)
      colnames(exbspline_basis)<- paste("xi",1:3,sep="")
      dat.region.target <- cbind(dat.region.target,exbspline_basis)
      
      if(nrow(dat.region.target) < 21){
        return(NULL)
      } else {
        ################################# model performance #############################################
        
        if(isTRUE(holiday)){
          dat.region.target <- merge(dat.region.target,vaccnation,by.x ="date",by.y ="date",all.x=T)
          if(isTRUE(spline_set)){
            model.lm <- summary(lm(log(rc)~cum + factor(work)+xi1+xi2+xi3,data = dat.region.target))
          } else {
            model.lm <- summary(lm(log(rc)~cum + factor(work),data = dat.region.target))
          }
        } else {
          if(isTRUE(spline_set)){
            model.lm <- summary(lm(log(rc)~cum+xi1+xi2+xi3,data = dat.region.target))
          } else {
            model.lm <- summary(lm(log(rc)~cum,data = dat.region.target))
          }
        }
        
        inter_cept <- data.frame(model.lm$coefficients)[,c(1:2)]
        inter_cept$region <- region.id; inter_cept$year <- year.id 
        inter_cept$r2 <-  model.lm$r.squared 
        inter_cept$var <- row.names(inter_cept)  
        
        #model <- lm(log(rc)~cum+bs(time,df = 3),data = dat.region.target)
        #logrr <- predict(model,data = dat.region.target) - summary(model)$coefficients[1] - 
        # summary(model)$coefficients[2]*dat.region.target$cum
        #bs.dat  <- bs(dat.region.target$time,df = 3)
        #logrr <- bs.dat[,1]*model.aa$coefficients[3] +  bs.dat[,2]*model.aa$coefficients[4] + bs.dat[,3]*model.aa$coefficients[5]
        #model.lm <- lm(log(rc)~cum,data = dat.region.target)
        #dat.region.season$adjust.rc <- exp(model.lm$residuals + summary(model.lm)$coefficients[1])
        
        return(inter_cept)
      }
      #    if(model.lm$r.squared > 0.1){
      #     return(dat.region.target[-1,]) 
      #    } else {
      #      return(NULL)
      #    } 
      #  } else {
      #    return(dat.region.target[-1,])
      #  }
    })
    adjust_rc_resaon <- bind_rows(adjust_list)
    return(adjust_rc_resaon) 
  })
  adjust_rc_resaon_final <- bind_rows(result_list)
  return(adjust_rc_resaon_final)
}
