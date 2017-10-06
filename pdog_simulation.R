
source('pdog_utility.R')

package_load(c("reshape2", "dplyr", "LaplacesDemon", "runjags", "mcmcplots",
  "coda", "igraph", "foreach", "doParallel"))

psim <- readRDS("model13_1.RDS")

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

# make covariates
cc <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "nearest_pd"))) %>% tail(., 384)

cc2 <- cc[rep(row.names(cc), 10),]

to_add_age <- rep(1:10, each = 384)
cc2$frag.age <- cc2$frag.age + to_add_age
# scale by the way the previous covariates were

X <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "nearest_pd")))

m_means <- apply(X, 2, mean)
m_sd <- apply(X, 2, sd)


cc2 <- sweep(cc2, 2, m_means)
cc2 <- sweep(cc2, 2, m_sd, "/")
# add ones
cc2 <- cbind(1, cc2)

pcov <- array(0, dim = c(384, 10, 6))

# fill pcov
for(i in 1:6){
  pcov[,,i] <- cc2[,i]
}

# to use for extinction
e_sp <- c(1,3,4,6)
d_sp <- c(1,5)
c_sp <- c(1,3,4,6)

coords <- pd_temp %>% select(one_of(c("easting", "northing")))

pstate = status[,6]

make_tpm_once <- function(eb = NULL, cb = NULL, sb=NULL, 
  pcov = NULL, yr = NULL, e_sp = NULL, 
  d_sp = NULL, c_sp = NULL,
  pstate = NULL, coords = NULL, m_means = NULL, m_sd = NULL, my_samp = NULL,
  sds = NULL){
  
  tpm <- array(0, dim = c(3, 3, 384))
  
  # change covariates a bit based on last state
  
  my_dogs <- which(pstate == 3)
  my_frags <- which(pstate %in% c(2:3))
  
  nearest_dogs <- nearest_frag <- rep(0, 384)
  for(j in 1:384) {
    # attach one sample to the previous year
    if(j %in% my_dogs){
      mds <- my_dogs[-which(my_dogs == j)]
    } else {
      mds <- my_dogs
    }
    
    if(j %in% my_frags){
      mfs <- my_frags[-which(my_frags == j)]
    } else {
      mfs <- my_frags
    }
    
    to_dist_pd <- rbind(coords[j,], coords[mds,])
    
    to_dist_frag <- rbind(coords[j,], coords[mfs,])
    
    # calculate distances
    my_dist_pd <- to_dist_pd %>% select(one_of(c("easting", "northing"))) %>% 
      dist(, diag = TRUE, upper = TRUE) %>% as.matrix
    
    my_dist_frag <- to_dist_frag %>% select(one_of(c("easting", "northing"))) %>% 
      dist(, diag = TRUE, upper = TRUE) %>% as.matrix
    
    # sorting, then grab second value (first value is 0)
    nearest_dogs[j] <- sort(as.numeric(my_dist_pd[,1]))[2]
    nearest_frag[j] <- sort(as.numeric(my_dist_frag[,1]))[2]
    # if the 2nd value is 0 (should not be), grab 3rd
    
  }
  
  nearest_dogs <- (nearest_dogs - m_means[5])/m_sd[5]
  
  pcov[,yr,6] <- nearest_dogs
  
  samp_sb <- sb[my_samp,]
  samp_cb <- cb[my_samp,]
  samp_eb <- eb[my_samp,]
  
  # fill dvlp
  tpm[1, 1, ] <- 1
  
  for(i in 1:dim(tpm)[3]){
    # from frag exists to
    tpm[1,2,i] <- 1 - plogis((samp_sb %*% pcov[i,yr,d_sp])) # dvlp
    tpm[2,2,i] <- plogis((samp_sb %*% pcov[i,yr,d_sp])) * 
      (1 - plogis((samp_cb %*% pcov[i,yr,c_sp]))) # frag
    tpm[3,2,i] <- plogis((samp_sb %*% pcov[i,yr,d_sp])) * 
      plogis((samp_cb %*% pcov[i,yr,c_sp]))# pdogs
    # from pdogs to
    tpm[1,3,i] <- 1 - plogis((samp_sb %*% pcov[i,yr,d_sp])) # developed
    tpm[2,3,i] <- plogis((samp_sb %*% pcov[i,yr,d_sp]))*
      plogis((samp_eb %*% pcov[i,yr,e_sp])) # frag
    tpm[3,3,i] <- plogis((samp_sb %*% pcov[i,yr,d_sp]))*
      (1 - plogis((samp_eb %*% pcov[i,yr,e_sp]))) #pdogs
  }
  
  # sample based off of the previous state
  new_state = rep(0, 384)
  for(i in 1:384){
    new_state[i] <- rcat(1, tpm[1:3, pstate[i],i ])
  }
  
  return(new_state)
  
  
}

one_samp = sample(1:nrow(eb), 100, replace = TRUE)

one_go <- array(0, dim = c(384, 10, 100))
for(iter in 1:100){
  
  one_go[,1,iter] <- make_tpm_once(eb = eb, cb = cb, sb = sb, pcov = pcov, yr = 1,
    e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, msd = msd,
    pstate = pstate, my_samp = one_samp[iter])
  for(sim in 2:10){
    one_go[,sim,iter] <- make_tpm_once(eb = eb, cb = cb, sb = sb, pcov = pcov, yr = sim,
      e_sp = e_sp,d_sp = d_sp, c_sp = c_sp, coords = coords, msd = msd,
      pstate = one_go[,sim-1,iter], one_samp[iter])
  }
  
}

my_ans <- array(0, dim = c(3, 10, 100))

ones <- twos <- tre <- matrix(0, ncol = 10, nrow = 28)

for(i in 1:28){
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
