##---- libraries ----

library(tidyverse)
# library(qcc)
library(foreach)
library(doParallel)
library(reshape2)

##---- functions ----

# source utilities.R, qcc.R and rules.R from fork of qcc

source("utilities.R")
source("qcc.R")
source("rules.R")
source("ewma.R")
source("cusum.R")
source("describe.R")


### simulation functions

#' Create jump change in set of process data sequences
#'
#' @param x -- matrix of random sequences; each column is a process data sequence 
#' @param sigma -- process standard deviation
#' @param jump_size -- size of jump in terms of sigma
#' @param jump_time -- integer, time in sequence in which jump occurs
#'
#' @return -- a matrix, x + described jump
#' @export
#'
#' @examples
jump_change = function(x,sigma,jump_size,jump_time) {
  n_seq = ncol(x)
  len_seq = nrow(x)
  jmp = c(rep(0,jump_time-1),
          rep(sigma*jump_size,len_seq-jump_time+1))
  
  ret = sweep(x,1,jmp,"+")
  
  return(ret)
}


#' Add trend to set of process data sequences
#'
#' @param x -- matrix of random sequences; each column is a process data sequence
#' @param sigma -- process standard deviation
#' @param trend_slope -- slope of trend in terms of sigma
#' @param trend_start_time -- integer, time period when trend starts
#'
#' @return -- a matrix, x + described trend
#' @export
#'
#' @examples
trend_change = function(x,sigma,trend_slope,trend_start_time) {
  n_seq = ncol(x)
  len_seq = nrow(x)
  trnd = c(rep(0,trend_start_time-1),
           (1:(len_seq-trend_start_time+1))*sigma*trend_slope)
  
  ret = sweep(x,1,trnd,"+")
  
  return(ret)
}


all_rules_violations_comb = function(dat,n_est,newdat,
                                nsigmas = 3,
                                nMR = 3.668,
                                lambda= 0.2,
                                h=5,
                                k=0.25) {
  
  # return min index or NA
  min_viol = function(viol,n1) if(length(viol)>0) min(viol) else n1

  
  n = length(newdat)
  n1 = n+1
  
  # first the moving range
  MR = abs(diff(dat))
  avgMR = mean(MR)
  MR_lim = nMR * avgMR
  
  new_MR = abs(diff(newdat))
  vMR = which(new_MR > MR_lim) + 1 # returns numeric(0) or violation indices + 1
  MR1 = min_viol(vMR,n1)
  
  # create qcc object
  qcc_obj = qcc(data=dat[1:n_est],type="xbar.one",newdata = newdat,
                     rules = NULL,
            plot=F) #, nsigmas = nsigmas)
  # find first violations of alll rules 
  v1 = min_viol(qccRulesViolatingWER1(qcc_obj),n1)
  v2 = min_viol(qccRulesViolatingWER2(qcc_obj),n1)
  v3 = min_viol(qccRulesViolatingWER3(qcc_obj),n1)
  v4 = min_viol(qccRulesViolatingWER4(qcc_obj),n1)

  # ewma
  ewma_obj = ewma(data=dat[1:n_est],newdata = newdat, lambda = lambda,
                  plot = F) # ,nsigmas=nsigmas)
  e1 = min_viol(ewma_obj$violations,n1)
 
  # cusum
  cusum_obj = cusum(data=dat[1:n_est],newdata = newdat, 
                    decision.interval = h, se.shift = k,
                    plot=F)
  c1 = min_viol(c(cusum_obj$violations$upper,cusum_obj$violations$lower),n1)

  
  out = data.frame(r1=v1,r2=v2,r3=v3,r4=v4,e1=e1,c1=c1,mr1 = MR1) %>%
    mutate(R1 = r1,
           R1.R2 =  min(r1,r2,na.rm = T),
           R1.MR =   min(r1,mr1,na.rm = T),
           R1.R2.MR =   min(r1,r2,mr1, na.rm = T),
           R1.R3.MR =   min(r1,r3,mr1,na.rm = T),
           R1.R4.MR =   min(r1,r4,mr1,na.rm = T),
           R1.R2.R3.MR =   min(r1,r2,r3,mr1,na.rm = T),
           R1.R2.R4.MR =   min(r1,r2,r4,mr1,na.rm = T),
           R1.R3.R4.MR =   min(r1,r3,r4,mr1,na.rm = T),
           R1.R2.R3.R4.MR =   min(r1,r2,r3,r4,mr1,na.rm = T),
           R1.E1 =   min(r1,e1,na.rm = T),
           R1.MR.E1 =   min(r1,mr1,e1,na.rm = T),
           R1.C1 =   min(r1,mr1,na.rm = T),
           R1.MR.C1 =   min(r1,mr1,c1,na.rm = T)) %>%
          unlist()
  
  
  return(out)  
}



vec_rules_violations = function(x, dat,n_est,newdat,
                                nsigmas = 3, nMR = 3.668,
                                lambda= 0.2,
                                h=5,
                                k=0.25) {
  
  dat = dat[,x]
  newdat = newdat[,x]
  
  out=all_rules_violations_comb(dat,n_est=n_est,newdat,
                                nsigmas = nsigmas,
                                nMR = nMR,
                                lambda= lambda,
                                h=h,
                                k=k) 
  return(out)
}



sim_jump = function(base_dat, est_dat, est_size, step_size,
                    nsigmas = 3, nMR = 3.668,
                    lambda= 0.2,
                    h=5,
                    k=0.25) {
  
  step_dat = base_dat + step_size
  
  out = sapply(1:ncol(step_dat),
               FUN = vec_rules_violations,
               dat = est_dat,
               n_est = est_size,
               newdat = step_dat,
               nsigmas = nsigmas,
               nMR = nMR,
               lambda= lambda,
               h=h,
               k=k)
  
  out = t(out)
  rnames = dimnames(out)[[2]]
  
  nr = nrow(out)
  out = cbind(rep(step_size,nr),
              rep(est_size,nr),
              rep(nsigmas,nr),
              rep(nMR,nr),
              rep(lambda,nr),
              rep(h,nr),
              rep(k,nr),
              out)
  
  dimnames(out)[[2]] = c("step_size","est_size","n_sigmas",
                         "k_MR", "lambda","h", "k",   rnames)
  
  return(out)
}




##---- parameters ----

rand_seed = 19671031
set.seed(rand_seed)


seq_len = 1e3
n_seq = 1e4
nsims = seq_len * n_seq
sigma = 1
mu = 0

# sample sizes for estimating chart parameters
est_n = c(15,30, 60, 120)

# lambda values for EWMA
lambda_vals = c(0.1,0.2,0.3)

# h and k values for CUSUM
h_vals = 4:5
k_vals = c(0,0.25,0.5)

# step change
jmp_size = c(seq(0,to=1.5, by = 0.25),2,3)

# linear trend
slope_size = c(0.1,0.2,0.5,1)




##---- generate base  data ----

# chart parameter estimation data
xs = matrix(rnorm(n_seq * max(est_n), mu,sigma), ncol=n_seq)


# base data for charts
xx = matrix(rnorm(nsims, mu,sigma), ncol=n_seq)



## ---- set up cluster ----


n_cores = detectCores()
use_cores = 1
if(n_cores > 2) use_cores = n_cores-1
if(n_cores > 3) use_cores = n_cores-2

registerDoParallel(cores=use_cores)



##---- do step change sims ----



### change to CC code
 step_sims =
   foreach(es = est_n,js = jmp_size, 
           ls = lambda_vals,
           hs = h_vals,
           ks = k_vals,
          .combine=rbind,
          .multicombine = T,
          .inorder = F,
          .packages = c("qcc")
          ) %dopar%   sim_jump(base_dat = xx, est_dat =  xs,
                               est_size =  es, step_size = js,
                               nsigmas = 3, nMR = 3.668,
                               lambda = ls,h = hs,k = ks)

system.time({
step_sims = NULL

for(es in est_n) {
  for(js in jmp_size) {
    for(ls in lambda_vals) {
      for(hs in h_vals) {
        for(ks in k_vals) {
      tmp = sim_jump(base_dat = xx, est_dat =  xs,
                     est_size =  es, step_size = js,
                     nsigmas = 3, nMR = 3.668,
                     lambda = ls,h = hs,k = ks)
      step_sims = rbind(step_sims,tmp)
      gc()
        }}}}
  print(es)
  }
})

write.csv(step_sims,"Jump_sims_EC2.csv")



##---- deregister implicit cluster ----------

 stopImplicitCluster()

##---- summarize ----


step_sims_tbl = step_sims %>% as_tibble() 

step_sims_narrow = step_sims_tbl %>%
  melt(id.vars=1:2, measure.vars=3:12,variable.name="Rules",value.name = "RL")

step_sims_rules = step_sims_narrow %>%
  group_by(step_size, est_size, Rules) %>%
  summarize(prob_max = mean(RL == 1001),
            mean = mean(RL),
            median = median(RL))


write.csv(step_sims_sum,"Jump_sims_sum.csv")

