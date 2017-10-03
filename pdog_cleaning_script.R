#######################
# Cleaning script
#######################


source("pdog_utility.R")
pdat <- read.csv("pdog_data.csv", header = TRUE)


library(reshape2)
library(dplyr)
library(LaplacesDemon)
library(runjags)
library(mcmcplots)
library(coda)
library(igraph)

# get the status data
fragstatus <- melt(pdat, id = "FRAG.ID", 
  measure.vars = grep("^st\\.\\d\\d\\d\\d", colnames(pdat))  ) %>% 
  filter( nchar(as.character(FRAG.ID)) > 1 )

# convert the variable data
fragstatus$variable <- pull_year(fragstatus)

# convert column names
colnames(fragstatus)[2:3] <- c("year", "frag.status")

# do the same with pd
pdstatus <- melt(pdat, id = "FRAG.ID", 
  measure.vars = grep("^pd\\.\\d\\d\\d\\d", colnames(pdat))  )%>% 
  filter( nchar(as.character(FRAG.ID)) > 1 )

# remove variable because it is already in fragstatus
pdstatus <- pdstatus %>% select(-one_of("variable"))

# change column names 
colnames(pdstatus)[2] <- "pd.status"


# variables we want to keep
to_keep <- c("FRAG.ID", "frag.age", "easting", "northing", 
  "frag.area.2002")

red_pdat <- pdat %>% select(one_of(to_keep)) %>% 
   filter(complete.cases(.))

# make column names shorter
colnames(red_pdat) <- c("FRAG.ID", "frag.age", "easting", "northing", "area")

# combine all of these data. First, get pdstatus and fragstatus together.
# pretty simple because they are ordered identically

fragstatus$pd.status <- pdstatus$pd.status

# join them together

pd <- left_join( fragstatus, red_pdat, by = "FRAG.ID")

# make frag.age change with each year
pd$frag.age <- pd$frag.age + (pd$year - min(pd$year))

# time since study started
pd$time <- pd$year - (min(pd$year) - 1)

# calculate the three types of events that could occur.

uyear <- unique(pd$year)

each_year <- vector("list", length = length(uyear)-1)

for(i in 2:length(uyear)){
  frag.pro <- paste(pd$frag.status[pd$year==uyear[i-1]], 
    pd$frag.status[pd$year==uyear[i]], sep = "-" )
  pd.pro <- paste(pd$pd.status[pd$year==uyear[i-1]], 
    pd$pd.status[pd$year==uyear[i]], sep = "-" )
  
  # previous p dogs
  p_pd_year <- pd[pd$year == uyear[i-1] & pd$pd.status == 1,]
  # previous frag
  p_frag_year <- pd[pd$year == uyear[i-1] & pd$frag.status == 1,]
  
  frag_areas <- tcrossprod(p_frag_year$area)
  # current p dogs
  myyear <- pd[pd$year == uyear[i],]
  # store distance info
  nearest_dogs <- nearest_frag <- aw_dogs <- aw_frag <- rep(0, nrow(myyear))
  
  # calculate connectivity stuff
  
  ppd <- p_frag_year %>% select(one_of(c("easting", "northing"))) %>% 
    dist(, diag = TRUE, upper = TRUE) %>% as.matrix
  ppd[ppd>2000] <- 0
  colnames(ppd) <- p_frag_year$FRAG.ID
  
  net2 <- graph_from_adjacency_matrix(ppd, weighted = TRUE)

  sq.m <- 347000000
  
 
  hm <- distances(net2, V(net2), 
                  weights = E(net2)$weight, to = V(net2))
  diag(hm) <- 1

  
  hm3 <- exp(-0.00006*hm) * frag_areas

  
  PC <- (sum(hm3)) / (sq.m^2)
  
  ack3 <- lose1(net2, PC, p_frag_year$area, 0.00006)
  

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
    aw_dogs[j] <- (prod(to_dist_pd$area[c(1,order(my_dist_pd[,1])[2])])^0.7)/
      (nearest_dogs[j]^1.7)
    nearest_frag[j] <- sort(as.numeric(my_dist_frag[,1]))[2]
    aw_frag[j] <- (prod(to_dist_frag$area[c(1,order(my_dist_frag[,1])[2])])^0.7)/
      (nearest_frag[j]^1.7)
    # if the 2nd value is 0 (should not be), grab 3rd
    if(nearest_dogs[j] == 0) {
      nearest_dogs[j] <- sort(my_dist_pd[,1])[3]
      
    }
    
    if(nearest_frag[j] == 0) {
      nearest_frag[j] <- sort(as.numeric(my_dist_frag[,1]))[3]
      
    }
    
  }
    
  
  each_year[[i-1]] <- pd[pd$year == uyear[i],]
  each_year[[i-1]]$frag.status <- frag.pro
  each_year[[i-1]]$pd.status <- pd.pro
  each_year[[i-1]]$nearest_pd <- nearest_dogs
  each_year[[i-1]]$nearest_frag <- nearest_frag
  each_year[[i-1]]$aw_pd <- aw_dogs
  each_year[[i-1]]$aw_frag <- aw_frag
  
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
  "area", "time", "nearest_pd", "nearest_frag","aw_pd", "aw_frag", "status")))

pd_temp$nearest_frag <- 1 / pd_temp$nearest_frag
pd_temp$nearest_pd <- 1 / pd_temp$nearest_pd

had_pd <- apply(status, 1, function(x) min(which(x==3)))
had_pd[is.infinite(had_pd)] <- 0
had_pd[had_pd>0] <- 1

# covariates
X <- pd_temp %>% select(one_of(c("frag.age", "easting", "northing",
  "area", "time", "nearest_pd", "nearest_frag", "aw_pd", "aw_frag")))


X[,] <- scale(X) 

# get attributes

msd <- attributes(X)

X <- as.matrix(X)
covs <- array(NA, dim = c(384, 6, 10))
for(i in 1:10){
  if(i == 10){
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




