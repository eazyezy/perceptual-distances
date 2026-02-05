library(dplyr)
library(tidyr)
library(tibble)
library(adehabitatHR)

#### 1. Read in the data (change qp <-> lad as needed) ####
d <- read.csv("data\\placements\\qp_placements.csv") %>%
  arrange(stimulus) #%>%
  mutate(stimulus = factor(stimulus))

#### 2. Compute the overlap coefficient #### 

# initialize matrix
m <- matrix(NA,17,17,dimnames=list(unique(d$stimulus),unique(d$stimulus)))

# loop over all combinations of stimuli
for (i in unique(d$stimulus)) {
  for (j in unique(d$stimulus)) {
    # get stimuli pair
    stimuli <- d %>%
      filter(stimulus %in% c(i,j)) %>%
      select(stimulus)
    # get placements coordinates
    coords <- d %>%
      filter(stimulus %in% c(i,j)) %>%
      select(placement_x,placement_y)
    
    # make spatial dataframe
    placements_sp <- SpatialPointsDataFrame(stimuli,coords)
    
    if (i!=j) {
      # compute BA values and store in matrix
      m[i,j] <- round(kerneloverlap(placements_sp, method='BA')[1,2],3) #%>% print()
    } else {
      m[i,j] <- 1
    }
      
  }
}

#### 3. Make distances symmetric #### 
symm <- function(x) round(sqrt(1-x^2),2)

m_sym <- apply(m,1:2,symm)

#### 4. Write results to file ####

#final dataset fed to Gabmap
write.csv(m_sym,"distances\\demo\\qp_distances.csv")
