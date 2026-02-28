library(dplyr)
library(tidyr)
library(tibble)
library(adehabitatHR)

#### 1. Read in the data (change qp <-> lad as needed) ####
d <- read.csv("data\\placements\\qp_placements.csv")

# stimuli
stims <- unique(d$stimulus)

# descriptive statistics
mean(d$error_in_km)
max(d$error_in_km)
boxplot(d$error_in_km)

#### 2. Compute the overlap coefficient #### 

# initialize matrix
m <- matrix(data=NA,
            nrow=length(stims),
            ncol=length(stims),
            dimnames=list(stims,stims))

# loop over all combinations of stimuli
for (i in stims) {
  for (j in stims) {
    # get stimuli pair
    stimuli <- as.data.frame(d[d$stimulus %in% c(i,j),1])
    # get placements coordinates
    coords <- d[d$stimulus %in% c(i,j),2:3]
    # make spatial dataframe
    placements_sp <- SpatialPointsDataFrame(coords,stimuli)
    # compute BA values and store in matrix
    if (i!=j) {
      m[i,j] <- round(kerneloverlap(placements_sp, method='BA')[1,2],2)
    } else {
      m[i,j] <- 1
    }
  }
}

#### 3. Make distances symmetric #### 

ba2hd <- function(x) round(sqrt(1-x),2)

m_sym <- apply(m,1:2,ba2hd)

#### 4. Write results to file ####

#final dataset fed to Gabmap
write.table(m_sym
            ,"distances\\lad_distances.txt"
            ,sep="\t"
            ,quote = F
            ,row.names = T)
