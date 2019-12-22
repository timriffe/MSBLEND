#try four different starting age blendng approaches

# Get prev, transfers, etc, and try to remake the Brouard graph. Is it possible with
# data just estimated from chronological data, or do we need a separate estimation procedure?


# helpful for data operations
library(tidyverse)
library(here)
library(reshape2)
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

m11 <- TRsub[,"m11"]
m12 <- TRsub[,"m12"]
m13 <- TRsub[,"m14"]

m21 <- TRsub[,"m21"]
m22 <- TRsub[,"m22"]
m23 <- TRsub[,"m24"]

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

Nlong <- U %>% 
  U2N() %>% 
  melt(varnames = c("to","from"),
       value.name = "time") %>% 
  mutate(to = as.character(to),
         from = as.character(from)) %>% 
  separate(col = "to", 
           sep = "::",
           into = c("state_to","age_to"),
           convert = TRUE) %>% 
  separate(col = "from", 
           sep = "::",
           into = c("state_from","age_from"),
           convert = TRUE) %>% 
  rename(age=age_to) %>% 
  filter(age > age_from,
         state_to != "D")

# calculate total time spend in each age and state,
# should have as many versions as origin ages.
time_to <- Nlong %>% 
  group_by(state_from, age, age_from) %>% 
  summarize(time = sum(time)) %>% 
  arrange(state_from, age_from, age)

head(time_to)
head(Nlong)

# tempting to do lifetable averaging of lx to get Lx
# but will resist for the time being. Still to double
# check age alignments.

TRsub %>% 
  rename(H_H = m11,
         H_U = m12,
         H_D = m14,
         U_U = m22,
         U_H = m21,
         U_D = m24
         ) %>% 
  pivot_longer(cols = H_H:U_D,
               names_to = c("state_from", "state_to"),
               names_sep = "_",
               values_to = "prob") %>% 
  # ensure exact composition, possibly redundant
  group_by(state_from, age) %>% 
  mutate(prob = prob / sum(prob))

# TODO:
# 1) merge probabilities
# 2) get event counts relative to age 48 starting pop
# 3) figure out how to get backwards probabilities. 
# is it just as simple as event/time_to ?




