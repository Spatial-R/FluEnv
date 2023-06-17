data("Flu2009")

############################# define the DicreteSIDistr for Influenza

MeanFluSI <- 2.6
SdFluSI <- 1.9
DicreteSIDistr <- vector()
for(i in 0:8)
{
  DicreteSIDistr[i+1] <- DiscrSI(i, MeanFluSI, SdFluSI)
}

Rc <- function(data,date.id,case.id,n.num){
  ## data     :  the analysis dataset for
  ## date.id  :  which column has the information about date
  ## case.id  :  which column has the information about the case number
  ## n.num    :  time span
  ## start.date : 
  ## time.interval:  daily or weekly
  dat.1 <- data
  col.number <- which(names(dat.1) %in% c(date.id,case.id))
  dat.1 <- dat.1[,c(col.number)]
  names(dat.1) <- c("date","cases")
  #n.first <- ifelse(time.interval == "week",1, which(dat.1[,"date"] == start.date))
  #dat.1 <- dat.1[-c(1:(n.first - 10)),]
  n.total <- dim(dat.1)[1]
  result <- EstimateR(dat.1[,"cases"],T.Start = 1:(n.total - n.num), 
                      T.End = (1 + n.num):n.total, method = "NonParametricSI",
                      SI.Distr=DicreteSIDistr, plot=F)
  Rc.R.1 <- result$R; 
 # Rc.R.1 <- filter(Rc.R,T.Start %in% c(1:(n.total - n.num)))
  Rc.R.1[,"date"] <- dat.1[c(1:(n.total - n.num)),"date"]
  Rc.R.1[,"cases"] <- dat.1[c(1:(n.total - n.num)),"cases"]
  return(Rc.R.1)
}

WT.E <- function(data,date.id,case.id,n.num){
  ## data     :  the analysis dataset for
  ## date.id  :  which column has the information about date
  ## case.id  :  which column has the information about the case number
  ## n.num    :  time span
  ## start.date : 
  ## time.interval:  daily or weekly
  dat.1 <- data
  col.number <- which(names(dat.1) %in% c(date.id,case.id))
  dat.1 <- dat.1[,c(col.number)]
  names(dat.1) <- c("date","cases")
  #n.first <- ifelse(time.interval == "week",1, which(dat.1[,"date"] == start.date))
  #dat.1 <- dat.1[-c(1:(n.first - 10)),]
  n.total <- dim(dat.1)[1]
  result <- WT(dat.1[,"cases"],T.Start = 1:(n.total - n.num), 
                      T.End = (1 + n.num):n.total, method = "NonParametricSI",
                      SI.Distr=DicreteSIDistr, plot=F)
  Rc.R <- result$R; 
  Rc.R.1 <- filter(Rc.R,T.Start %in% c(seq(2,dim(Rc.R)[1],n.num)))
  Rc.R.1[,"date"] <- dat.1[seq(2,dim(Rc.R)[1],n.num),"date"]
  Rc.R.1[,"cases"] <- dat.1[seq(2,dim(Rc.R)[1],n.num),"cases"]
  return(Rc.R.1)
}



daily.incidence.caculate <- function(data,date.id="date",
                                     case.id ="count",
                                     city.id = "region",
                                     pop.id ="pop",year_sep = TRUE){
  dat.tem <- data
  var.pos <- unlist(lapply(c(date.id,case.id,city.id,pop.id),function(data)which(names(dat.tem) == data)))
  dat.tem <- dat.tem[,var.pos];names(dat.tem) <- c("date","count","region","pop")
  dat.tem$year <- as.numeric(substr(dat.tem$date,1,4))
  
  if(isTRUE(year_sep)){
    
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
      
      ##################### smooth line for the sum incidence up to week t
      
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




WeeklyRc.caculate <- function(data,case.id="case",region.id ="region",date.id="date"){
  dat.tem <- data
  var.pos <- unlist(lapply(c(date.id,case.id,region.id),function(data)which(names(dat.tem) == data)))
  dat.tem <- dat.tem[,var.pos]; names(dat.tem) <- c("date","case","region")
  #dat.tem <- arrange(dat.tem,date)
  rc.list <- lapply(unique(dat.tem$region),function(region.id){
    dat.tem.1 <- filter(dat.tem,region == region.id)
    rc.result <- Rc(data = dat.tem.1,date.id = "date",case.id = "case",n.num = 1)[,c(3,6,10,12)]
    names(rc.result) <- c("rc","rcl","rch","date")
    rc.result$week <- floor_date(rc.result$date,"week")
    rc.week <- data.frame(summarise(group_by(rc.result,week),
                                    rc = mean(rc,na.rm = T),
                                    rcl = mean(rcl,na.rm = T), 
                                    rch = mean(rch,na.rm = T)))
    rc.week$region <- region.id
    return(rc.week)
  })
  week.rc <- bind_rows(rc.list)
  return(week.rc)
}

analysis_range <- function(data,region.var,interval.v){
  dat.tem <- data
  names(dat.tem)[which(names(dat.tem) == region.var)] <- "region"
  
  date_list <- lapply(unique(dat.tem$region), function(region.id){
    dat.region <- filter(dat.tem, region == region.id)
    from_year <- min(year(dat.region$date)) + 1
    to_year <-   max(year(dat.region$date))
    dat.target.final <- data.frame()
    for (i in seq(from_year,to_year)){
      date.range <- seq.Date(as.Date(paste(i-1,"-10-01",sep = "")),
                             as.Date(paste(i,"-05-01",sep = "")),by = "day")
      dat.date <- filter(dat.region,date %in% date.range)
      target.date <- dat.date[which.max(dat.date$case),"date"]
      dat.result <- data.frame(from = target.date - week.range*interval.v,
                               to = target.date + week.range*interval.v,
                               region = region.id)
      dat.target.final <- rbind(dat.target.final,dat.result)
    }
    return(dat.target.final)
  })
  date_result <- bind_rows(date_list)
  return(date_result)
}

