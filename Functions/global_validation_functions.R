library(pomp)
#######################################################################################
################################# some structure for pomp #############################
#######################################################################################


scale_year <- function(data,case.id,date.id = "time"){
  dat.tem <- data.frame(data)
  dat.tem$target <- dat.tem[,which(names(dat.tem) == case.id)]
  dat.tem$date <- dat.tem[,which(names(dat.tem) == date.id)]
  dat.tem <- mutate(dat.tem,year = as.numeric(substr(date,1,4)))
  result.list <- lapply(unique(dat.tem$year), function(id){
    dat.tem.1 <- filter(dat.tem,year == id)
    dat.tem.1$target <- scale(dat.tem.1$target)
    return(dat.tem.1)
  })
  result.fin <- bind_rows(result.list)
  return(result.fin)
}

model_build_rt <- function(data){

rproc <- Csnippet("
                  double beta, seas, foi, gh, beta0,dw, births;
                  double rate[6], trans[6];
                  
                  // term-time seasonality
                  t = (t-floor(t))*365.25;
                  // wuyi: 5.1-5.3; shiyi: 2.1-2.7;summer:7.1-8.31;yuandan:1.1
                  if ( (t>=90&&t<=121) || (t>=274&&t<=304) || (t>=4&&t<40) || (t>=124 && t<212) || (t>= 314 && t<365))
                  seas = 1.0+amplitude*0.34/0.66;
                  else
                  seas = 1.0-amplitude;
                  
                  // transmission rate

                  beta = R0*(1 + beta1*cos(2*pi*time));
                  
                  // expected force of infection
                  gh   = pow(I+iota,alpha);
                  foi  = beta*gh/pop;    // add the effect of temp on suspect and
                  
                  // white noise (extrademographic stochasticity)
                  dw = rgammawn(sigmaSE,dt);
                  
                  rate[0] = foi*dw/dt;      // stochastic force of infection
                //  rate[0] = foi;
                  rate[1] = mu;             // natural S death
                  rate[2] = gamma;          // recovery
                  rate[3] = mu;             // natural I death
                  rate[4] = mu;             // natural R death
                  rate[5] = guta;           // natural R to S
                  
                  // Poisson births
                  births = rpois(birthrate*dt); 
                  
                  // transitions between classes
                  reulermultinom(2,S,&rate[0],dt,&trans[0]);
                  reulermultinom(2,I,&rate[2],dt,&trans[2]);
                  reulermultinom(2,R,&rate[4],dt,&trans[4]);
                  
                  S += births   + trans[5] - trans[0] - trans[1];
                  I += trans[0] - trans[2] - trans[3];
                  R += trans[2] - trans[4] - trans[5];
                  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                  C += trans[0];           // true incidence
                  ")

initlz <- Csnippet("
                   double m = pop/(S_0+I_0+R_0);
                   S = nearbyint(m*S_0);
                   I = nearbyint(m*I_0);
                   R = nearbyint(m*R_0);
                   W = 0;
                   C = 0;
                   ")

dmeas <- Csnippet("
                  double m = rho * C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                  }
                  ")

rmeas <- Csnippet("
                  double m = rho*C;
                  double v = m*(1.0-rho+psi*psi*m);
                  double tol = 1.0e-18;
                  cases = rnorm(m,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")

toEst <- Csnippet("
                  Talpha = log(alpha);
                  Tiota = log(80-iota);
                  Trho = logit(rho);
                  Tamplitude = logit(amplitude);
                  TsigmaSE = log(sigmaSE);
                  Tpsi = log(psi);
                  to_log_barycentric (&TS_0, &S_0, 3);
                  ")

fromEst <- Csnippet("
                    Talpha = exp(alpha);
                    Tiota = 80-exp(iota);
                    Trho = expit(rho);
                    Tamplitude = expit(amplitude);
                    TsigmaSE = exp(sigmaSE);
                    Tpsi = exp(psi);
                    from_log_barycentric (&TS_0, &S_0, 3);
                    ")

col.case <- which(names(data) %in% c("week","flu"));
col.covar <-  which(!names(data) %in% c("flu","date","country"))
cases.bj <- data[,c(col.case)];names(cases.bj)[1] <- "cases"
covar.bj <- data[,c(col.covar)]
covar.bj$pop <- 5800000; covar.bj$birthrate <- covar.bj$pop*0.008

cases.bj %>%
  pomp(t0 = cases.bj$week[1],
       time = "week",
       rprocess = euler.sim(step.fun=rproc,delta.t=1/365.25),
       initializer = initlz,
       dmeasure = dmeas,
       rmeasure = rmeas,
       covar = covar.bj,
       toEstimationScale = toEst,
       fromEstimationScale = fromEst ,
       tcovar = "week",
       zeronames = c("C","W"),
       statenames = c("S","I","R","C","W"),
       paramnames = c("R0","mu","gamma","alpha","iota",
                      "rho","sigmaSE","psi","amplitude","guta",
                      "S_0","I_0","R_0")
  ) -> origin.model
return(origin.model)
}





model_build_pnas <- function(data){
  
  rproc <- Csnippet("
                    double beta, seas, foi, gh, beta0,dw, births;
                    double rate[6], trans[6];
                    
                    // term-time seasonality
                    t = (t-floor(t))*365.25;
                    // wuyi: 5.1-5.3; shiyi: 2.1-2.7;summer:7.1-8.31;yuandan:1.1
                    if ( (t>=90&&t<=121) || (t>=274&&t<=304) || (t>=4&&t<40) || (t>=124 && t<212) || (t>= 314 && t<365))
                    seas = 1.0+amplitude*0.34/0.66;
                    else
                    seas = 1.0-amplitude;
                    
                    // transmission rate
                    
                    beta = Rt*seas/dt;
                    
                    // expected force of infection
                    gh   = pow(I+iota,alpha);
                    foi  = beta*gh/pop;    // add the effect of temp on suspect and
                    
                    // white noise (extrademographic stochasticity)
                    dw = rgammawn(sigmaSE,dt);
                    
                    rate[0] = foi*dw/dt;      // stochastic force of infection
                    //  rate[0] = foi;
                    rate[1] = mu;             // natural S death
                    rate[2] = gamma;          // recovery
                    rate[3] = mu;             // natural I death
                    rate[4] = mu;             // natural R death
                    rate[5] = guta;           // natural R to S
                    
                    // Poisson births
                    births = rpois(birthrate*dt); 
                    
                    // transitions between classes
                    reulermultinom(2,S,&rate[0],dt,&trans[0]);
                    reulermultinom(2,I,&rate[2],dt,&trans[2]);
                    reulermultinom(2,R,&rate[4],dt,&trans[4]);
                    
                    S += births   + trans[5] - trans[0] - trans[1];
                    I += trans[0] - trans[2] - trans[3];
                    R += trans[2] - trans[4] - trans[5];
                    W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                    C += trans[0];           // true incidence
                    ")
  
  initlz <- Csnippet("
                     double m = pop/(S_0+I_0+R_0);
                     S = nearbyint(m*S_0);
                     I = nearbyint(m*I_0);
                     R = nearbyint(m*R_0);
                     W = 0;
                     C = 0;
                     ")
  
  dmeas <- Csnippet("
                    double m = rho * C;
                    double v = m*(1.0-rho+psi*psi*m);
                    double tol = 1.0e-18;
                    if (cases > 0.0) {
                    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
                    } else {
                    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
                    }
                    ")
  
  rmeas <- Csnippet("
                    double m = rho*C;
                    double v = m*(1.0-rho+psi*psi*m);
                    double tol = 1.0e-18;
                    cases = rnorm(m,sqrt(v)+tol);
                    if (cases > 0.0) {
                    cases = nearbyint(cases);
                    } else {
                    cases = 0.0;
                    }
                    ")
 
col.case <- match(c("week","flu"),names(data))
col.covar <-  which(!names(data) %in% c("flu","date","country"))
cases.bj <- data[,c(col.case)];names(cases.bj)[2] <- "cases"
covar.bj <- data[,c(col.covar)]
#covar.bj$pop <- 5800000; covar.bj$birthrate <- covar.bj$pop*0.008

cases.bj %>%
  pomp(t0 = cases.bj$week[1],
       time = "week",
       rprocess = euler(step.fun=rproc,delta.t=1/365.25),
       rinit = initlz,
       dmeasure = dmeas,
       rmeasure = rmeas,
       covar =  covariate_table(covar.bj,times = "week"),
       accumvars = c("C","W"),
       statenames = c("S","I","R","C","W"),
       paramnames = c("R0","mu","gamma","alpha","iota",
                      "rho","sigmaSE","psi","amplitude","guta",
                      "S_0","I_0","R_0")
  ) -> origin.model

  return(origin.model)
}



season_measure <- function(data,case_var,date_var){
  lh.result <- data
  names(lh.result)[which(names(lh.result) == case_var)] <- "cases"
  names(lh.result)[which(names(lh.result) == date_var)] <- "time"
  lh.result$target <- scale(ifelse(lh.result$cases == 0, 0, log(lh.result$cases)))
  lh.result$time.1 <- floor((lh.result$time - floor(lh.result$time))*53)
  lm.result <- summary(lm(target~cos(2*pi*time.1/52.17) + sin(2*pi*time.1/52.17)+ cos(4*pi*time.1/52.17)  + sin(4*pi*time.1/52.17), data = lh.result))
  lm.cof <- data.frame(lm.result$coefficients)
  seas <-data.frame(ann.1 = sqrt(lm.cof[2,1]^2 + lm.cof[3,1]^2),
                    ann.2 = sqrt(lm.cof[4,1]^2 + lm.cof[5,1]^2),
                    peak.t = atan2(lm.cof[3,1], lm.cof[2,1]))
  return(seas)
}



season_peak <- function(data,case_var,date_var){
  dat_tem <- data
  names(dat_tem)[which(names(dat_tem) == case_var)] <- "cases"
  names(dat_tem)[which(names(dat_tem) == date_var)] <- "time"
  #lh.result$year <- as.numeric(substr(lh.result$date,1,4))
  dat_tem$year <- as.numeric(floor(dat_tem$time))
  week_res <- unlist(lapply(unique(dat_tem$year),function(year_id){
    dat.tem <- filter(dat_tem,year == year_id)
    if(nrow(dat.tem) > 40){
      week_num <- week(dat.tem[which.max(dat.tem$cases),"date"])
    } else {
      week_num <- NULL
    }
    return(week_num)
  }))
  if(!is.null(week_res)){
    seas <- data.frame(peak = week_res) 
  } else {
    seas <- data.frame(peak = NA)
  }
  
  return(seas)
}




season_peak_gisaid <- function(data,case_var,date_var,sample_var){
  dat_tem <- data
  names(dat_tem)[which(names(dat_tem) == case_var)] <- "cases"
  names(dat_tem)[which(names(dat_tem) == date_var)] <- "time"
  names(dat_tem)[which(names(dat_tem) == sample_var)] <- "sample"
  #lh.result$year <- as.numeric(substr(lh.result$date,1,4))
  dat_tem$year <- as.numeric(floor(dat_tem$time))
  week_res <- unlist(lapply(unique(dat_tem$year),function(year_id){
    #dat.tem <- filter(dat_tem, (Year == year_id & Week %in% c(1:39)) | (Year == (year_id-1) & Week %in% c(41:53)))
    dat.tem <- filter(dat_tem,year == year_id)
    dat.tem <- filter(dat.tem,sample > 20 & sample > ALL_INF)
    if(nrow(dat.tem) > 20){
      week_num <- week(dat.tem[which.max(dat.tem$cases),"date"])
    } else {
      week_num <- NULL
    }
    #print(paste(year_id,"||",week_num,sep = ""))
    return(week_num)
  }))
  if(!is.null(week_res)){
    seas <- data.frame(peak = week_res) 
  } else {
    seas <- data.frame(peak = NA)
  }
  return(seas)
}



peak_diff <- function(dat_1,dat_2){
res_fin <- c()
  for (i in c(1:nrow(dat_1))){
    for (j in c(1:nrow(dat_2))){
      #print(paste(i,":",j,sep = ""))
     def <-  abs(dat_1[i,1] - dat_2[j,1]) > 26
      if(def){
        dat_1_dis <- abs(52 - dat_1[i,1])
        dat_2_dis <- abs(52 - dat_2[j,1]) 
        if(dat_1_dis > dat_2_dis){
          res_tem <- dat_2_dis + dat_1[i,1]
        } else {
          res_tem <- dat_1_dis + dat_2[j,1]
        }
      } else {
        res_tem <- abs(dat_1[i,1] - dat_2[j,1])
      }
     res_fin <- c(res_fin,res_tem)
    }
  }
  return(min(res_fin))
}
