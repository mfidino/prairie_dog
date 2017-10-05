# First data munge

# This script just does some very initial formatting of the prairie dog data

source("pdog_utility.R")
pdat <- read.csv("pdog_data.csv", header = TRUE)

package_load(c("reshape2", "dplyr", "LaplacesDemon", "runjags", "mcmcplots",
  "coda", "igraph"))



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


write.csv(pd, "base_pdog_data.csv")