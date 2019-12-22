
# helpful for data operations
library(tidyverse)
library(here)
# some custom functions
source(here::here("R","Functions.R"))

# Read in the data (Daniel Schneider gave these to me, haha)
TR <- readRDS(here::here("Data","TR_v06_collapsed.rds"))

# take a peek: we have many strata
head(TR)

# These are the things we can subset on.
unique(TR$time) 
unique(TR$sex) 
unique(TR$edu) 
unique(TR$age)

# define a subset
TRsub <- TR %>% filter(sex == "f",
              edu == "terciary",
              time == 1996)

# starting proportions in each state
init <- TRsub[1,c("s1_prop","s2_prop")] %>% unlist()
names(init) <- c("H","D")

# Make the submatrices of U, the transient matrix
HH <- pi2u(pivec = TRsub[,"m11"], from = "H", to = "H")
HU <- pi2u(pivec = TRsub[,"m12"], from = "H", to = "U")
UH <- pi2u(pivec = TRsub[,"m21"], from = "U", to = "H")
UU <- pi2u(pivec = TRsub[,"m22"], from = "U", to = "U")

# we need to bind these like this:
# |-------|
# | HH UH |
# | HU UU |
# |-------|

U <- u2U(HH = HH, # healthy to healthy
         HU = HU, # healthy to unhealthy
         UH = UH, # unhealthy to healthy
         UU = UU) # unhealthy to unhealthy

# so far this is all matrix architecture, now we have 
# the transient matrix, which we can transform to the 
# fundamental matrix

N <- U2N(U)
N
# This thing has what we want: conditional expected time spent in each
# state and age (conditional on survival and starting state!)
# There are matrix tricks to select what we want. But this is where
# we instead use tidyverse to grab what we need:
library(reshape2)
N %>% 
  melt(varnames = c("to","from"),
       value.name = "time") %>% 
  mutate(to = as.character(to),
         from = as.character(from)) %>% 
  separate(col = "to", 
           sep = "::",
           into = c("to","age2"),
           convert = TRUE) %>% 
  separate(col = "from", 
           sep = "::",
           into = c("from","age1"),
           convert = TRUE) %>% 
  filter(age1 == 50,
         from != "D",
         to != "D") %>% 
  group_by(from,to) %>% 
  summarize(Ex_cond = sum(time)) %>% 
  mutate(init = ifelse(from == "H", init[1],init[2]),
         Ex = Ex_cond * init) %>% 
  group_by(to) %>% 
  summarize(Ex = sum(Ex))
  


