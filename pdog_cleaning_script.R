#######################
# Cleaning script
#######################

# this starts formatting the data for actual analysis

source("pdog_utility.R")

package_load(c("reshape2", "dplyr", "LaplacesDemon", "runjags", "mcmcplots",
  "coda", "igraph", "foreach", "doParallel"))

pd <- read.csv("base_pdog_data.csv", header = TRUE, 
  stringsAsFactors = FALSE)

# calculate the three types of events that could occur.

uyear <- unique(pd$year)
usites <- unique(pd$FRAG.ID)

# going to store results for each year into a list, then 
# combine all of them together

each_year <- vector("list", length = length(uyear)-1)

for(i in 2:length(uyear)){
  # pastes together the fragments previous status to its current status
  frag.pro <- paste(pd$frag.status[pd$year==uyear[i-1]], 
    pd$frag.status[pd$year==uyear[i]], sep = "-" )
  # pastes together the pds previous status to its current status
  pd.pro <- paste(pd$pd.status[pd$year==uyear[i-1]], 
    pd$pd.status[pd$year==uyear[i]], sep = "-" )
  
  # previous p dogs
  p_pd_year <- pd[pd$year == uyear[i-1] & pd$pd.status == 1,]
  # previous frag
  p_frag_year <- pd[pd$year == uyear[i-1] & pd$frag.status == 1,]
  
  # current fragments
  myyear <- pd[pd$year == uyear[i],]
  
  # calculate connectivity stuff
pck <- calc_pck(p_frag_year)
pck <- data.frame(FRAG.ID = p_frag_year$FRAG.ID, pck = pck,
  stringsAsFactors = FALSE)
pck <- left_join(data.frame(FRAG.ID = usites,
  stringsAsFactors = FALSE), pck, by = "FRAG.ID")
disties <- calc_distances(p_frag_year, p_pd_year, myyear)

  each_year[[i-1]] <- pd[pd$year == uyear[i],]
  each_year[[i-1]]$frag.status <- frag.pro
  each_year[[i-1]]$pd.status <- pd.pro
  each_year[[i-1]]$nearest_pd <- disties$nearest_dogs
  each_year[[i-1]]$nearest_frag <- disties$nearest_frag
  each_year[[i-1]]$aw_pd <- disties$aw_dogs
  each_year[[i-1]]$aw_frag <- disties$aw_frag
  each_year[[i-1]]$pck <- pck$pck.1
  each_year[[i-1]]$pck_pd <- pck$pck.2
  
}

# bind all of the years together again
pd_temp <- plyr::rbind.fill(each_year)
pd_temp$status <- 0

# for fragment. 1-0 = 1 (developed) or 0-0
pd_temp$status[pd_temp$frag.status == "1-0"] <- 1
pd_temp$status[pd_temp$frag.status == "0-0"] <- 1

# for fragment. 1-1 = 2 (frag exists)
pd_temp$status[pd_temp$frag.status == "1-1"] <- 2

# for pd. 1-1 = 3
pd_temp$status[pd_temp$pd.status == "1-1"] <- 3


# for pd. 1-0 is good from the fragment stuff.

# for pd. 0-1 = 3
pd_temp$status[pd_temp$pd.status == "0-1"] <- 3


status <- matrix(pd_temp$status, ncol = 6, nrow = 384)

pd_2002 <- pd[pd$year == 2002, ] %>% 
  select(one_of(c("pd.status", "frag.age", "easting", "northing","area")))

# add 2002 data onto status
status <- cbind(pd_2002$pd.status + 2, status)


pd_temp <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "time", "nearest_pd", "nearest_frag","aw_pd", "aw_frag", "pck",
  "pck_pd","status")))

pd_temp$nearest_frag <- 1 / pd_temp$nearest_frag
pd_temp$nearest_pd <- 1 / pd_temp$nearest_pd

# will give warnings, can be ignored
had_pd <- apply(status, 1, function(x) min(which(x==3)))
had_pd[is.infinite(had_pd)] <- 0
had_pd[had_pd>0] <- 1

# if there is an NA value in pck it means that the site is gone. We
# will then give it an importance value of zero

pd_temp$pck[is.na(pd_temp$pck)] <- 0
pd_temp$pck_pd[is.na(pd_temp$pck_pd)] <- 0

# covariates
X <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "time", "nearest_pd", "nearest_frag", "aw_pd", "aw_frag", "pck",
  "pck_pd")))




X[,] <- scale(X) 

# get attributes

msd <- attributes(X)

X <- as.matrix(X)
covs <- array(NA, dim = c(384, 6, 12))
for(i in 1:12){
  if(i == 12){
    covs[,,i] <- had_pd
  } else{
  covs[,,i] <- X[,i]
}
}




Y <- pd_2002 %>% select(one_of(c("frag.age", "easting", "northing", "area")))
Y <- scale(Y) %>% as.matrix()


all_data <- list(covs = covs, occ_covs = Y, status = status, 
                 occ_status = pd_2002$pd.status, msd = msd)

saveRDS(all_data, "cleaned_pdog.RDS")

# everything should now be ready for analysis


