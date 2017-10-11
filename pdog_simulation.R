
source('pdog_utility.R')

package_load(c("reshape2", "dplyr", "LaplacesDemon", "runjags", "mcmcplots",
  "coda", "igraph", "foreach", "doParallel"))

psim <- readRDS("model17_4_2017-10-06.RDS")

psim <- as.mcmc.list(psim) %>% as.matrix(., chains = TRUE)

# remove p_prob stuff

psim <- psim[,-grep("p_prob", colnames(psim))]

# get intercepts

ints <- psim[,grep("0", colnames(psim))]

# betas
betas <- psim[,-grep("0", colnames(psim))]

# standard deviations
sds <- betas[,grep("sd", colnames(betas))]

betas <- betas[,-grep("sd", colnames(betas))]
# remove insignificant terms


# will equal 2 if significant
signs <- apply(betas, 2, HDIofMCMC) %>% sign() %>% apply(., 2, sum) %>% abs

to_med <- which(signs == 0) %>% as.numeric()

my_meds <- as.numeric(apply(betas[,to_med], 2, median)) %>% 
  rep(., each = nrow(betas)) %>% 
  matrix(., ncol = length(to_med), nrow = nrow(betas))

betas[1:nrow(betas),to_med] <- my_meds

extn <- cbind(ints[,colnames(ints)=="e0"], betas[,grep("e_", colnames(betas))])
colnames(extn) <- 1:ncol(extn)

coln <- cbind(ints[,colnames(ints)=="c0"], betas[,grep("c_", colnames(betas))])
colnames(extn) <- 1:ncol(coln)

surv <- cbind(ints[,colnames(ints)=="d0"], betas[,grep("d_", colnames(betas))])
colnames(surv) <- 1:ncol(surv)

# make covariates for projections
cc <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "nearest_pd", "pck"))) %>% tail(., 384)

cc2 <- cc[rep(row.names(cc), 10),]

to_add_age <- rep(1:10, each = 384)
cc2$frag.age <- cc2$frag.age + to_add_age
# scale by the way the previous covariates were

X <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "nearest_pd", "pck")))

m_means <- apply(X, 2, mean)
m_sd <- apply(X, 2, sd)

# make covariates for historic simulations

hist_covs <- covs[,,c(1,1,2,3,4,6,10)]
hist_covs[,,1] <- 1


cc2 <- sweep(cc2, 2, m_means)
cc2 <- sweep(cc2, 2, m_sd, "/")
# add ones
cc2 <- cbind(1, cc2)

pcov <- array(0, dim = c(384, 10, 7))

# fill pcov
for(i in 1:7){
  pcov[,,i] <- cc2[,i]
}

# covaraites used for each transition
e_sp <- c(1,3,4,6)
d_sp <- c(1,7)
c_sp <- c(1,3,4,6)

# get fragment data from the last year

yr_5_frags <- data.frame(FRAG.ID = unique(pd$FRAG.ID),
                         tail(pd_temp, 384))

coords <- pd_temp %>% select(one_of(c("easting", "northing")))

# state in 2017
pstate = status[,6]


start <-Sys.time()
one_samp = sample(1:nrow(extn), 100, replace = FALSE)

one_go <- array(0, dim = c(384, 10, length(one_samp)))

do_sim <- function(eb= NULL, cb = NULL, sb = surv, pcov = NULL,
  e_sp = NULL, d_sp = NULL, c_sp = NULL, coords = NULL, m_means = NULL,
  sds = NULL, m_sd = NULL, pstate = NULL, my_samp = NULL, raw_frag = NULL){
  single_samp <- matrix(0, ncol = dim(pcov)[2], nrow = 384)
  single_samp[,1] <- make_tpm_once(eb = eb, cb = cb, sb = sb, 
    pcov = pcov, yr = 1,
    e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, 
    m_means = m_means, sds = sds, m_sd = m_sd,
    pstate = pstate, my_samp = my_samp, raw_frag = raw_frag)
  for(sim in 2:dim(pcov)[2]){
    single_samp[,sim] <- make_tpm_once(eb = eb, cb = cb, sb = sb, 
      pcov = pcov, yr = sim,
      e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, 
      m_means = m_means, sds = sds, m_sd = m_sd,
      pstate = single_samp[,sim-1], my_samp = my_samp, raw_frag = raw_frag)
    
  }
  return(single_samp)
}


cores <- detectCores()-2
cl <- makeCluster(cores)
registerDoParallel(cl)
start <- Sys.time()
proj_2007 <- foreach(i = 1:10002, 
  .packages = c('magrittr', 'dplyr', 'igraph', 'parallel', 'doParallel',
    'LaplacesDemon')) %dopar% {
   do_sim(eb = extn, cb = coln, sb = surv, 
    pcov = pcov,
    e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, 
    m_means = m_means, sds = sds, m_sd = m_sd,
    pstate = pstate, my_samp = i, raw_frag = yr_5_frags)
    }

saveRDS(proj_2007, "proj_2007.RDS")
stopCluster(cl)
end <- Sys.time()
cores <- detectCores()-2
cl <- makeCluster(cores)
registerDoParallel(cl)
proj_2002 <- foreach(i = 1:10002, 
  .packages = c('magrittr', 'dplyr', 'igraph', 'parallel', 'doParallel',
    'LaplacesDemon')) %dopar% {
      do_sim(eb = extn, cb = coln, sb = surv, 
        pcov = hist_covs[,-6,],
        e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, 
        m_means = m_means, sds = sds, m_sd = m_sd,
        pstate = status[,1], my_samp = i, raw_frag = yr_5_frags)
    }
saveRDS(proj_2002, "proj_2002.RDS")
stopCluster(cl)





proj_2002 <- foreach(i = 1:10)
end <- Sys.time()

for(iter in 1:){
  
  one_go[,1,iter] <- make_tpm_once(eb = extn, cb = coln, sb = surv, 
    pcov = pcov, yr = 1,
    e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, 
    m_means = m_means, sds = sds, m_sd = m_sd,
    pstate = pstate, my_samp = one_samp[iter], raw_frag = yr_5_frags)
  for(sim in 2:10){
    one_go[,sim,iter] <- make_tpm_once(eb = extn, cb = coln, sb = surv, 
      pcov = pcov, yr = sim,
      e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, 
      m_means = m_means, sds = sds, m_sd = m_sd,
      pstate = one_go[,sim-1,iter], my_samp = one_samp[iter], raw_frag = yr_5_frags)
  }
  
}

end <- Sys.time()

ack <- array(unlist(proj_2002), dim = c(nrow(proj_2002[[1]]), 
  ncol(proj_2002[[1]]), length(proj_2002)))

my_ans <- array(0, dim = c(3, 10, length(one_samp)))

ones <- twos <- tre <- matrix(0, ncol = 10, nrow = 10)

on <- tw <- tr <- matrix(0, nrow = 5, ncol = 5) # sim by year matrix
for(i in 1:5){
  on[i,] <- colSums(ack[,,i]==1)
  tw[i,]<-  colSums(ack[,,i]==2)
  tr[i,]<-  colSums(ack[,,i]==3)
}

for(i in 1:10){
  ones[i,] <- colSums(one_go[,,i]==1)
  twos[i,] <- colSums(one_go[,,i]==2)
  tre[i,] <- colSums(one_go[,,i]==3)
}

oo <- apply(ones, 2, quantile, probs = c(0.025,0.5,0.975))
tt <- apply(twos, 2, quantile, probs = c(0.025,0.5,0.975))
tr <- apply(tre, 2, quantile, probs = c(0.025,0.5,0.975))

windows(10, 5)
par(mfrow= c(1, 3))

plot(oo[2,], ylim = range(oo), type = 'l', ylab = "sites developed")
lines(oo[1,], lty = 2)
lines(oo[3,], lty = 2)
points(x = 10, y = sum(status[,7]==1), pch = 19)

plot(tt[2,], ylim = range(tt), type = 'l', ylab = "fragments no pdogs")
lines(tt[1,], lty = 2)
lines(tt[3,], lty = 2)
points(x = 10, y = sum(status[,7]==2), pch = 19)

plot(tr[2,], ylim = range(tr), type = 'l', ylab = "fragments w/ pdogs")
lines(tr[1,], lty = 2)
lines(tr[3,], lty = 2)
points(x = 10, y = sum(status[,7]==3), pch = 19)





table(status[,7])


for()
  
  ns <- apply(new_state, 2, median)

hm <- apply(new_state, 2, table)

# previous p dogs
p_pd_year <- pd[pd$year == uyear[i-1] & pd$pd.status == 1,]
# previous frag
p_frag_year <- pd[pd$year == uyear[i-1] & pd$frag.status == 1,]
# current p dogs
myyear <- pd[pd$year == uyear[i],]
# store distance info
nearest_dogs <- nearest_frag <- rep(0, nrow(myyear))
for(j in 1:nrow(myyear)) {
  # attach one sample to the previous year
  to_dist_pd <- rbind(myyear[j,], p_pd_year) %>% 
    filter(!duplicated(FRAG.ID))
  
  to_dist_frag <- rbind(myyear[j,], p_frag_year) %>% 
    filter(!duplicated(FRAG.ID))
  
  # calculate distances
  my_dist_pd <- to_dist_pd %>% select(one_of(c("easting", "northing"))) %>% 
    dist(, diag = TRUE, upper = TRUE) %>% as.matrix
  
  my_dist_frag <- to_dist_frag %>% select(one_of(c("easting", "northing"))) %>% 
    dist(, diag = TRUE, upper = TRUE) %>% as.matrix
  
  # sorting, then grab second value (first value is 0)
  nearest_dogs[j] <- sort(as.numeric(my_dist_pd[,1]))[2]
  nearest_frag[j] <- sort(as.numeric(my_dist_frag[,1]))[2]
  # if the 2nd value is 0 (should not be), grab 3rd
  if(nearest_dogs[j] == 0) {
    nearest_dogs[j] <- sort(my_dist_pd[,1])[3]
    
  }
  
  if(nearest_frag[j] == 0) {
    nearest_frag[j] <- sort(as.numeric(my_dist_frag[,1]))[3]
    
  }
  
}

}

# last state
sl <- status[,6]

pr_extn <- extn %*%
  
  for(i in 1:nrow(coln)){
    
  }
