model{
  # set up probabilities for transition matrix
  for(site in 1:nsite){
    for(year in 1:5){
      logit(pr_surv[site, year]) <- d_time[year] + 
        inprod(d_beta, dcovs[site, year])
      logit(pr_coln[site, year]) <- c_time[year] + 
        inprod(c_beta, ccovs[site, year,])
      logit(pr_extn[site, year]) <- e_time[year] + 
        inprod(e_beta, ecovs[site, year,])
    }
  }
  # fill transition matrix
  for(site in 1:nsite){
    for(year in 1:5){
      # from developed to
      tpm[1,1, site, year] <- 1 # developed
      tpm[2,1, site, year] <- 0 # fragmented habitat
      tpm[3,1, site, year] <- 0 # pdogs
      # from frag exists to
      tpm[1,2, site, year] <- 1 - pr_surv[site, year] # developed
      tpm[2,2, site, year] <- pr_surv[site,year] * (1 - pr_coln[site,year]) # frag
      tpm[3,2, site, year] <- pr_surv[site, year] * pr_coln[site, year] # pdogs
      # from pdogs to
      tpm[1,3, site, year] <- 1 - pr_surv[site,year] # developed
      tpm[2,3, site, year] <- pr_surv[site,year]*pr_extn[site,year] # frag
      tpm[3,3, site, year] <- pr_surv[site,year]*(1 - pr_extn[site,year]) #pdogs
    }
  }
  
  
  for(si in 1:nsite){
    p2002[si] ~ dbern(pr_occu[si])
    logit(pr_occu[si]) <- o0 + inprod(o_beta, covs02[si,])
    for(yr in 2:nyear){
      status[si,yr] ~ dcat(tpm[1:3, status[si,yr-1], si, yr-1])
    }
  }
  
  # priors
  
  d0 ~ dlogis(0, 1)
  c0 ~ dlogis(0, 1)
  e0 ~ dlogis(0, 1)
  o0 ~ dlogis(0, 1)
  
  # ranef time
  for(i in 1:5){
    d_time[i] ~ dnorm(d0, d_tau)
    c_time[i] ~ dnorm(c0, c_tau)
    e_time[i] ~ dnorm(e0, e_tau)
  }
  
  d_sd ~ dgamma(1, 1)
  c_sd ~ dgamma(1, 1)
  e_sd ~ dgamma(1, 1)
  
  d_tau <- 1 / (d_sd^2)
  c_tau <- 1 / (c_sd^2)
  e_tau <- 1 / (e_sd^2)
  
  for(co in 1:ncov_c){
    c_beta[co] ~ dt(0, 2.5, 1)
    
  }
  
  for(eo in 1:ncov_e){
    e_beta[eo] ~ dt(0, 2.5, 1)
  }
  
  for(db in 1:ncov_d){
  d_beta[db] ~  dt(0, 2.5, 1)
  }
  
  for(k in 1:4){
    
    o_beta[k] ~  dt(0, 2.5, 1)
  }
  
  for( site in 1:(nsite)){
    p_prob[site, 1] <- (pr_occu[site] ^ p2002[site]) * 
      ((1 - pr_occu[site])^(1 - p2002[site]))
    for( yr in 2:6){
      p_prob[site, yr] <- tpm[status[site,yr], status[site,yr-1], site, yr-1]
    }
  }
  
  
  
}