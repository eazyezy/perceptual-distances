library(dplyr)
library(tidyr)
library(tibble)
library(adehabitatHR)

#### 1. Read in the data ####
d <- read.csv("placements\\qp_placements.csv") # change qp <-> lad as needed

#### 2. Get all combinations of stimuli ####
combos <- as.data.frame(t(combn(unique(d$stimulus),2))) %>%
  rename(A = V1, B = V2)

#### 3. Compute the overlap coefficient over all pairs of stimuli #### 
for (i in 1:nrow(combos)) {
  # get stimuli pair
  stimuli <- d %>%
    filter(stimulus %in% c(combos[i,1:2])) %>%
    select(stimulus)
  # get placements coordinates
  coords <- d %>%
    filter(stimulus %in% c(combos[i,1:2])) %>%
    select(placement_x,placement_y)
  
  # make spatial dataframe
  placements_sp <- SpatialPointsDataFrame(coords, stimuli)
  
  # compute BA values and store in combos
  combos[i,3] <- round(kerneloverlap(placements_sp, method='BA')[1,2],3)
}

#### 4. Derive the distance matrix #### 
dist <- combos %>%
  pivot_wider(names_from = B, values_from = V3,values_fill = NA) %>%
  column_to_rownames("A") %>%
  mutate(Vicobarone=NA,.before = 1) %>%
  add_row("Romagnese" = NA)
rownames(dist) <- colnames(dist)

bhatt <- dist
bhatt[lower.tri(bhatt)] = t(bhatt)[lower.tri(bhatt)]
diag(bhatt) <- 1

#### 5. Make distances symmetric #### 
symm <- function(x) sqrt(1-x^2)

bhatt_symm <- bhatt %>%
  mutate_all(symm) 

# final dataset fed to Gabmap
write.csv(bhatt_symm,"distances\\qp_distances.csv")