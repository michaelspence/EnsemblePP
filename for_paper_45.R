#setwd("C:/Users/MS23/OneDrive - CEFAS/BX033/North Sea")
dat45 <- read.csv("PP_RCP45X.csv",row.names = 1)

other<-read.csv("ersem_obs.csv")

mod_out_pp_past <- cbind(other[5:38,2],dat45[paste(1984:2017),])

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
last_year <- 2050
mod_run <- log(cbind(dat45[paste(1984:last_year),]))
mod_run_me <- colMeans(mod_run)

centre <- mean(log(mod_out_pp_past[5:30,1]))

dat <- list(
  N=length(1984:last_year),
  M=ncol(dat45),
  M_runs = t(t(mod_run)-mod_run_me),
  N_os = 5,
  N_oe = 30,
  obs = log(mod_out_pp_past[5:30,1]) - centre,
  obs_sd = sqrt((var(log(mod_out_pp_past[5:30,1]))/0.86) * 0.14), ## R^2 = 0.86
  l_disc_p=2,
  s_ldisc_p=2,
  init_t_p=2,
  tr_shape=2,tr_rate=5/2
)

fit45 <- sampling(st_mod2,data=dat,chains=1,iter=2000,control=list(adapt_delta =0.99,max_treedepth=12))
ex.fit_45 <- extract(fit45)

#save(ex.fit_45,file="PP/fit_45.Rdata")

dyn_kappa_raw_45 <- ex.fit_45$latent[,,1] 
#save(dyn_kappa_raw_45, file="NorthSeaModel/bayes_kappa_45.Rdata")

