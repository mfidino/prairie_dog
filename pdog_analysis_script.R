######################
# pdog analysis script
######################

source("pdog_utility.R")

library(reshape2)
library(dplyr)
library(LaplacesDemon)
library(runjags)
library(mcmcplots)
library(coda)

all_dat <- readRDS("cleaned_pdog.RDS")

# write each element in all_dat as an object
for(i in 1:length(all_dat)){
  
  assign(names(all_dat)[i], all_dat[[i]])
  
}

# covariate order
corder <- c("frag.age", "easting", "northing",
  "area", "time", "nearest_pd", "nearest_frag","aw_pd", "aw_frag", "had_pd")

# model
model_covs_ce <- list(list(c(2,3), ranef = TRUE) ,
  list(8, ranef = TRUE),
  list(9, ranef = TRUE),
  list(6, ranef = TRUE),
  list(7, ranef = TRUE),
  list(c(2,3), ranef = FALSE),
  list(8, ranef = F),
  list(9, ranef = F),
  list(6, ranef = F),
  list(7, ranef = F),
  list(ranef = TRUE),
  list(ranef = F),
  list(c(2,3,6), ranef = TRUE),
  list(c(2,3,8), ranef = TRUE),
  list(c(2,3,9), ranef = TRUE),
  list(c(2,3,6), ranef = TRUE),
  list(c(2,3,7), ranef = TRUE),
  list(c(2,3,8), ranef = F),
  list(c(2,3,9), ranef = F),
  list(c(2,3,6), ranef = F),
  list(c(2,3,7), ranef = F),
  list(c(2,3,1), ranef = TRUE),
  list(1, ranef = TRUE),
  list(c(2,3,1), ranef = F),
  list(1, ranef = F)
  )






model_covs_surv <- list(list(4), list(10), list(1))

models <- c("pdog_model_covariate_model_ranef.R",
            "pdog_model_ranef.R",
            "pdog_model_covariate_model.R",
            "pdog_model_intercept.R")

my_waic <- matrix(0, ncol = length(model_covs_surv), 
  nrow = length(model_covs_ce))

for(j in 1:length(model_covs_surv)){
for(i in 1:length(model_covs_ce)){
  # if ranef with covariates
 if(length(model_covs_ce[[i]])> 1 & model_covs_ce[[i]]$ranef){
   
   to_model <- list(status = status, ccovs = covs[,,model_covs_ce[[i]][[1]]], 
     p2002 = status[,1]-2, ecovs = covs[,,model_covs_ce[[i]][[1]]],
     covs02 = occ_covs, nsite = 384, 
     nyear = 6, ncov_c = length(model_covs_ce[[i]][[1]])
       , ncov_d = 1, ncov_e = length(model_covs_ce[[i]][[1]]),
     dcovs = covs[,,model_covs_surv[[j]][[1]]])
   
   inits <- function(chain){
     gen_list <- function(chain = chain){
       list( 
         c0 = rnorm(1),
         d0 = rnorm(1),
         e0 = rnorm(1),
         o0 = rnorm(1),
         c_beta = rnorm(to_model$ncov_c),
         e_beta = rnorm(to_model$ncov_e),
         d_beta = rnorm(to_model$ncov_d),
         o_beta = rnorm(4),
         .RNG.name = switch(chain,
           "1" = "base::Wichmann-Hill",
           "2" = "base::Marsaglia-Multicarry",
           "3" = "base::Super-Duper",
           "4" = "base::Mersenne-Twister",
           "5" = "base::Wichmann-Hill",
           "6" = "base::Marsaglia-Multicarry",
           "7" = "base::Super-Duper",
           "8" = "base::Mersenne-Twister"),
         .RNG.seed = sample(1:1e+06, 1)
       )
     }
     return(switch(chain,           
       "1" = gen_list(chain),
       "2" = gen_list(chain),
       "3" = gen_list(chain),
       "4" = gen_list(chain),
       "5" = gen_list(chain),
       "6" = gen_list(chain),
       "7" = gen_list(chain),
       "8" = gen_list(chain)
     )
     )
   }
   
   if(to_model$ncov_c > 1){
   
   mout <- run.jags(model = "pdog_model_covariate_model_ranef.R",
     monitor = c("d0", "d_beta", "e0", "e_beta", "c0", "c_beta",
       "o0", "o_beta", "d_sd","e_sd", "c_sd", "p_prob"),
     data = to_model,
     inits = inits,
     n.chains = 6,
     adapt = 2000,
     burnin = 5000,
     sample = ceiling(10000/6),
     thin = 2,
     summarise = FALSE,
     plots = FALSE,
     method = "parallel")
   }
   if(to_model$ncov_c == 1){
     mout <- run.jags(model = "pdog_covariate_model_ranef_one.R",
       monitor = c("d0", "d_beta", "e0", "e_beta", "c0", "c_beta",
         "o0", "o_beta", "d_sd","e_sd", "c_sd", "p_prob"),
       data = to_model,
       inits = inits,
       n.chains = 6,
       adapt = 2000,
       burnin = 5000,
       sample = ceiling(10000/6),
       thin = 2,
       summarise = FALSE,
       plots = FALSE,
       method = "parallel")
     
   }
   
   my_waic[i,j] <- as.mcmc.list(mout) %>% as.matrix(.,chains = TRUE) %>% 
     calc_waic(.)
   
   saveRDS(mout, paste0("model",i,"_",j,".RDS"))
   
   
 }
  # only ranef model
  if(length(model_covs_ce[[i]])==1 & model_covs_ce[[i]]$ranef){
    to_model <- list(status = status, 
      p2002 = status[,1]-2,
      covs02 = occ_covs, nsite = 384, 
      nyear = 6)
    
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          c0 = rnorm(1),
          d0 = rnorm(1),
          e0 = rnorm(1),
          o0 = rnorm(1),
          c_sd = rgamma(1,1,1),
          d_sd = rgamma(1,1,1),
          e_sd = rgamma(1,1,1),
          o_beta = rnorm(4),
          .RNG.name = switch(chain,
            "1" = "base::Wichmann-Hill",
            "2" = "base::Marsaglia-Multicarry",
            "3" = "base::Super-Duper",
            "4" = "base::Mersenne-Twister",
            "5" = "base::Wichmann-Hill",
            "6" = "base::Marsaglia-Multicarry",
            "7" = "base::Super-Duper",
            "8" = "base::Mersenne-Twister"),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,           
        "1" = gen_list(chain),
        "2" = gen_list(chain),
        "3" = gen_list(chain),
        "4" = gen_list(chain),
        "5" = gen_list(chain),
        "6" = gen_list(chain),
        "7" = gen_list(chain),
        "8" = gen_list(chain)
      )
      )
    }
    
    mout <- run.jags(model = "pdog_model_ranef.R",
      monitor = c("d0", "e0",  "c0", 
        "o0", "o_beta", "d_sd","e_sd", "c_sd", "p_prob"),
      data = to_model,
      inits = inits,
      n.chains = 6,
      adapt = 2000,
      burnin = 5000,
      sample = ceiling(10000/6),
      thin = 2,
      summarise = FALSE,
      plots = FALSE,
      method = "parallel")
    
    my_waic[i,j] <- as.mcmc.list(mout) %>% as.matrix(.,chains = TRUE) %>% 
      calc_waic(.)
    
    saveRDS(mout, paste0("model",i,"_",j,".RDS"))
    
    
  }
  # if homog with covariates
  if(length(model_covs_ce[[i]])> 1 & !model_covs_ce[[i]]$ranef){
    
    to_model <- list(status = status, ccovs = covs[,,model_covs_ce[[i]][[1]]], 
      p2002 = status[,1]-2, ecovs = covs[,,model_covs_ce[[i]][[1]]],
      covs02 = occ_covs, nsite = 384, 
      nyear = 6, ncov_c = length(model_covs_ce[[i]][[1]])
      , ncov_d = 1, ncov_e = length(model_covs_ce[[i]][[1]]),
      dcovs = covs[,,model_covs_surv[[j]][[1]]])
    
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          c0 = rnorm(1),
          d0 = rnorm(1),
          e0 = rnorm(1),
          o0 = rnorm(1),
          c_beta = rnorm(to_model$ncov_c),
          e_beta = rnorm(to_model$ncov_e),
          d_beta = rnorm(to_model$ncov_d),
          o_beta = rnorm(4),
          .RNG.name = switch(chain,
            "1" = "base::Wichmann-Hill",
            "2" = "base::Marsaglia-Multicarry",
            "3" = "base::Super-Duper",
            "4" = "base::Mersenne-Twister",
            "5" = "base::Wichmann-Hill",
            "6" = "base::Marsaglia-Multicarry",
            "7" = "base::Super-Duper",
            "8" = "base::Mersenne-Twister"),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,           
        "1" = gen_list(chain),
        "2" = gen_list(chain),
        "3" = gen_list(chain),
        "4" = gen_list(chain),
        "5" = gen_list(chain),
        "6" = gen_list(chain),
        "7" = gen_list(chain),
        "8" = gen_list(chain)
      )
      )
    }
    
    if(to_model$ncov_c > 1){
      
      mout <- run.jags(model = "pdog_model_covariate_model.R",
        monitor = c("d0", "d_beta", "e0", "e_beta", "c0", "c_beta",
          "o0", "o_beta", "d_sd","e_sd", "c_sd", "p_prob"),
        data = to_model,
        inits = inits,
        n.chains = 6,
        adapt = 2000,
        burnin = 5000,
        sample = ceiling(10000/6),
        thin = 2,
        summarise = FALSE,
        plots = FALSE,
        method = "parallel")
    }
    if(to_model$ncov_c == 1){
      mout <- run.jags(model = "pdog_covariate_model_one.R",
        monitor = c("d0", "d_beta", "e0", "e_beta", "c0", "c_beta",
          "o0", "o_beta", "d_sd","e_sd", "c_sd", "p_prob"),
        data = to_model,
        inits = inits,
        n.chains = 6,
        adapt = 2000,
        burnin = 5000,
        sample = ceiling(10000/6),
        thin = 2,
        summarise = FALSE,
        plots = FALSE,
        method = "parallel")
      
    }
    
    my_waic[i,j] <- as.mcmc.list(mout) %>% as.matrix(.,chains = TRUE) %>% 
      calc_waic(.)
    
    saveRDS(mout, paste0("model",i,"_",j,".RDS"))
    
    
  }
  # only homog model
  if(length(model_covs_ce[[i]])==1 & !model_covs_ce[[i]]$ranef){
    to_model <- list(status = status, 
      p2002 = status[,1]-2,
      covs02 = occ_covs, nsite = 384, 
      nyear = 6)
    
    inits <- function(chain){
      gen_list <- function(chain = chain){
        list( 
          c0 = rnorm(1),
          d0 = rnorm(1),
          e0 = rnorm(1),
          o0 = rnorm(1),
          o_beta = rnorm(4),
          .RNG.name = switch(chain,
            "1" = "base::Wichmann-Hill",
            "2" = "base::Marsaglia-Multicarry",
            "3" = "base::Super-Duper",
            "4" = "base::Mersenne-Twister",
            "5" = "base::Wichmann-Hill",
            "6" = "base::Marsaglia-Multicarry",
            "7" = "base::Super-Duper",
            "8" = "base::Mersenne-Twister"),
          .RNG.seed = sample(1:1e+06, 1)
        )
      }
      return(switch(chain,           
        "1" = gen_list(chain),
        "2" = gen_list(chain),
        "3" = gen_list(chain),
        "4" = gen_list(chain),
        "5" = gen_list(chain),
        "6" = gen_list(chain),
        "7" = gen_list(chain),
        "8" = gen_list(chain)
      )
      )
    }
    
    mout <- run.jags(model = "pdog_model_intercept.R",
      monitor = c("d0", "e0",  "c0", 
        "o0", "o_beta", "p_prob"),
      data = to_model,
      inits = inits,
      n.chains = 6,
      adapt = 2000,
      burnin = 5000,
      sample = ceiling(10000/6),
      thin = 2,
      summarise = FALSE,
      plots = FALSE,
      method = "parallel")
    
    my_waic[i,j] <- as.mcmc.list(mout) %>% as.matrix(.,chains = TRUE) %>% 
      calc_waic(.)
    
    saveRDS(mout, paste0("model",i,"_",j,".RDS"))
    
    
  }
  
    
  }
}
  
write.csv(my_waic, "waic_out_10_2.csv")

mo <- readRDS("model1_1.RDS") %>% as.mcmc.list()

caterplot(mo, "c_beta")


plot(a[-c(1,7)] ~ b)



m3 <- read.csv("waic_out_model_mult.csv")
