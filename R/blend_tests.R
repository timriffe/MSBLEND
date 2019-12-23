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

# starting proportions in each state,
# based on observed prevalence around age 50
init <- TRsub[1,c("s1_prop","s2_prop")] %>% unlist()
names(init) <- c("H","U")

# Make the submatrices of U, the transient matrix
U <- sub2U(TRsub)

# so far this is all matrix architecture, now we have 
# the transient matrix, which we can transform to the 
# fundamental matrix and then post-process to get what
# is needed from it

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

Nlong %>% head()


# calculate total time spent in each age and state,
# should have as many versions as origin ages.
time_to <- Nlong %>% 
  group_by(state_from, age, age_from) %>% 
  summarize(time = sum(time)) %>% 
  arrange(state_from, age_from, age)

# this can only happen with blending already in place.
# Hmmmm. Could do a blend with an assumed age 50 prev, turn
# it into an optimization problem.
Nlong %>% 
  group_by(state_to, age, age_from) %>% 
  summarize(time = sum(time)) %>% 
  arrange(state_to, age_from, age)


# cnvergence within 5 time steps 
# (see chronological convergence plots later in script)
time_to %>% 
  filter(age_from == 48) %>% 
  group_by(state_from) %>% 
  mutate(time = time / time[age == 60]) %>% 
  ggplot(aes(x = age, y = time, color = state_from)) + 
  geom_line() 
  dim(time_to)

# tempting to do lifetable averaging of lx to get Lx
# but will resist for the time being. Still to double
# check age alignments.

tidysub <- TRsub %>% 
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

head(tidysub)
Nlong %>% 
  filter(age_from == 48)


# TODO:
# 1) merge probabilities
# 2) get event counts relative to age 48 starting pop
# 3) figure out how to get backwards probabilities. 
# is it just as simple as event/time_to ?


# ------------------------------------ #
# chronologocal prevalence convergence #
# ------------------------------------ #
pop <- matrix(0,ncol=33,nrow=64)
pop[1:32,1] <- 1
for (i in 1:32){
  pop[,i+1] <-  U %*% pop[,i]
}
PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
PH[upper.tri(PH)] <- NA
PU[upper.tri(PU)] <- NA

pi1 <- PH / (PH + PU)

pop <- matrix(0,ncol=33,nrow=64)
pop[33:64,1] <- 1
for (i in 1:32){
  pop[,i+1] <-  U %*% pop[,i]
}
PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
PH[upper.tri(PH)] <- NA
PU[upper.tri(PU)] <- NA

pi2 <- PH / (PH + PU)

plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1))
for (i in 1:32){
  lines(a2[i:32],pi1[row(pi1) == (col(pi1) -1 + i)], col = "#FF000080")
  lines(a2[i:32],pi2[row(pi2) == (col(pi2) -1 + i)], col = "#0000FF80")
}


# Same thing: closer initial prev is to age x prev the faster it converges.
# and with .8/.2 ?
pop <- matrix(0,ncol=33,nrow=64)
pop[1:32,1] <- .8
pop[33:64,1] <- .2
for (i in 1:32){
  pop[,i+1] <-  U %*% pop[,i]
}
PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
PH[upper.tri(PH)] <- NA
PU[upper.tri(PU)] <- NA

pi1_d <- PH / (PH + PU)

pop <- matrix(0,ncol=33,nrow=64)
pop[1:32,1]  <- .2
pop[33:64,1] <- .8
for (i in 1:32){
  pop[,i+1] <-  U %*% pop[,i]
}
PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
PH[upper.tri(PH)] <- NA
PU[upper.tri(PU)] <- NA

pi2_d <- PH / (PH + PU)

a2 <- seq(50,112,by=2)
plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1))
for (i in 1:32){
  lines(a2[i:32],pi1_d[row(pi1) == (col(pi1) -1 + i)], col = "#FF000080")
  lines(a2[i:32],pi2_d[row(pi2) == (col(pi2) -1 + i)], col = "#0000FF80")
}

# --------------------------------------------------- #
# optimized starting prevalence                       #
# --------------------------------------------------- #

# pretty good, but not optimal. Really, should just optimize on first jump.
init_pop <- function(init, U){
  pop <- matrix(0,ncol=33,nrow=64)
  pop[1,1]  <- init
  pop[33,1] <- 1 - init
  for (i in 1:32){
    pop[,i+1] <-  U %*% pop[,i]
  }
  PH <- pop[1:32,1:32]
  PU <- pop[33:64,1:32]
  
  ph <- diag(PH)
  pu <- diag(PU)
  
  list(ph=ph,pu=pu)
}

init_min <- function(par =.9, U){
  p <- init_pop(init = par, U = U)
  ph <- p$ph
  pu <- p$pu
  # really will only affect first part of curve.
  sum(abs(diff(ph / (ph + pu))))
}

phoptim   <- optimize(interval = c(.7,1), f=init_min, U = U)$min
initoptim <- c(phoptim, 1 - phoptim)

p <- init_pop(init = phoptim, U = U)
ph <- p$ph
pu <- p$pu

plot(a2, ph / (ph + pu), ylim = c(0,1))

# --------------------------------------------------- #
# mix first probs many times?                         #
# --------------------------------------------------- #

u1 <- matrix(c(TRsub[1,"m11"],TRsub[1,"m12"],TRsub[1,"m21"],TRsub[1,"m22"]),2)
init_it <- c(.5,.5)
for (i in 1:30){
  init_it <- u1 %*% init_it
  init_it <- init_it / sum(init_it)
}
init_it

p <- init_pop(init = init_it[1], U = U)
ph <- p$ph
pu <- p$pu

plot(a2, ph / (ph + pu), ylim = c(0,1))
# this works better



