dat_hist <- read.csv("historical_North_Sea_intpp_annual_1850_2014.csv",row.names = 1)

other<-read.csv("ersem_obs.csv")

mod_out_pp_past <- cbind(other[5:35,2],dat_hist[paste(1984:2014),])
dat45 <- rbind(dat_hist[paste(1984:2014),],read.csv("ssp370_North_Sea_intpp_annual_2015_2100.csv",row.names = 1))* 3600*24*365 * 12.0107
dat85 <- rbind(dat_hist[paste(1984:2014),],read.csv("ssp585_North_Sea_intpp_annual_2015_2100.csv",row.names = 1))* 3600*24*365 * 12.0107

# ciso<- read.csv("PP/CISO.csv",row.names = 1)
# 
# plot(c(1984,2017),c(min(mod_out_pp_past,ciso,na.rm=T),max(mod_out_pp_past,na.rm=T)),type="n",ylab="PP",xlab="year",log="y")
# lines(1984:2017,mod_out_pp_past[,1])
# for(i in 1:ncol(dat45)){
#   lines(1984:2017,mod_out_pp_past[,i+1],col="red")
# }
# lines(1984:2017,ciso[paste(1984:2017),1],col="red")

library(rstan)
st_mod2 <- stan_model("ensemble_KF.stan")
last_year <- 2100
mod_run_45 <- log(cbind(dat45[paste(1984:last_year),]))
mod_run_me_45 <- colMeans(mod_run_45)

centre <- mean(log(mod_out_pp_past[5:30,1]))

dat_45 <- list(
  N=length(1984:last_year),
  M=ncol(dat45),
  M_runs = t(t(mod_run_45)-mod_run_me_45),
  N_os = 5,
  N_oe = 30,
  obs = log(mod_out_pp_past[5:30,1]) - centre,
  obs_sd = sqrt((var(log(mod_out_pp_past[5:30,1]))/0.86) * 0.14), ## R^2 = 0.86
  l_disc_p=2,
  s_ldisc_p=2,
  init_t_p=2,
  tr_shape=2,tr_rate=5/2
)

fit45 <- sampling(st_mod2,data=dat_45,chains=1,iter=4000,control=list(adapt_delta =0.99,max_treedepth=12))
ex.fit_45 <- extract(fit45)


dyn_kappa_raw_45 <- ex.fit_45$latent[,,1] 


mod_run_85 <- log(cbind(dat85[paste(1984:last_year),]))
mod_run_me_85 <- colMeans(mod_run_85)

dat_85 <- list(
  N=length(1984:last_year),
  M=ncol(dat85),
  M_runs = t(t(mod_run_85)-mod_run_me_85),
  N_os = 5,
  N_oe = 30,
  obs = log(mod_out_pp_past[5:30,1]) - centre,
  obs_sd = sqrt((var(log(mod_out_pp_past[5:30,1]))/0.86) * 0.14), ## R^2 = 0.86
  l_disc_p=2,
  s_ldisc_p=2,
  init_t_p=2,
  tr_shape=2,tr_rate=5/2
)


fit_85 <- sampling(st_mod2,data=dat_85,chains=1,iter=4000,control=list(adapt_delta =0.99,max_treedepth=12))
ex.fit_85 <- extract(fit_85)

dyn_kappa_raw_85 <- ex.fit_85$latent[,,1]

