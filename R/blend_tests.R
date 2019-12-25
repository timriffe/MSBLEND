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
a2 <- seq(50,112,by=2)
plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1))
for (i in 1:32){
  lines(a2[i:32],pi1[row(pi1) == (col(pi1) -1 + i)], col = "#FF000080")
  lines(a2[i:32],pi2[row(pi2) == (col(pi2) -1 + i)], col = "#0000FF80")
}

# try different starting prevalence.

phoptim   <- optimize(interval = c(.7,1), f=init_min, U = U)$min
initoptim <- c(phoptim, 1 - phoptim)

p  <- init_pop(init = phoptim, U = U)
ph <- p$ph
pu <- p$pu

lines(a2, ph / (ph + pu))
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
# again using the back probability-derived starting prop
rp <- init_pop(init = 0.9158685 , U = U)
rph <- rp$ph
rpu <- rp$pu
plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1))
for (i in 1:32){
  lines(a2[i:32],pi1[row(pi1) == (col(pi1) -1 + i)], col = "#FF000080")
  lines(a2[i:32],pi2[row(pi2) == (col(pi2) -1 + i)], col = "#0000FF80")
}
lines(a2, ph / (ph + pu), lty = 1,lwd=1,col="magenta")
lines(a2, rph / (rph + rpu), lty = 2,lwd=1)


plot(a2,ph / (ph + pu) - rph / (rph + rpu))
# this works better

# --------------------------------------------------- #
# get 

head(Nlong)
time_to <- Nlong %>% 
  filter(age_from == 48) %>% 
  group_by(state_from, age) %>% 
  mutate(time= time * ifelse(state_from == "H", init_it[1], init_it[2])) %>% 
  ungroup() %>% 
  group_by(state_to,age) %>% 
  summarize(time = sum(time)) %>% 
  arrange(state_to, age)

# take peek at occupancy curves
time_to %>% 
  ggplot(aes(x = age, y = time, color = state_to)) +
  geom_line()

tidysub <-
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
  mutate(prob = prob / sum(prob)) %>% 
  filter(state_to != "D") %>% 
  rename(year = time,
         state_in = state_to) %>% 
  glimpse()


# try calc transfers, presently these are wrong:
# transfers need deno of time in state_from, not state_in.
# these need to be old transfers not new ones. Somewhere we
# have an alignment problem.
  time_to %>% 
  rename(state_in = state_to) %>% 
  left_join(tidysub) %>% 
  mutate(transfers = prob * time) %>% 
  select(state_in, state_from, age, time, prob, transfers) %>% 
  glimpse()

# need new notation. (in)_(from)
dat <- TRsub

ID <- 
  TRsub %>% 
  IDLT() %>% 
  mutate(tr_hhx = hhx * Hx,
         tr_hux = hux * Hx,
         tr_uux = uux * Ux,
         tr_uhx = uhx * Ux,
         # reverse probabilities: 
         # (note interstate transfer directions)
         r_hh = tr_hhx / lead(Hx),
         r_hu = tr_uhx / lead(Hx),
         r_uu = tr_uux / lead(Ux),
         r_uh = tr_hux / lead(Ux))

# check
# ID$tr_hhx[-32] + ID$tr_uhx[-32] - ID$Hx[-1]
# 
# # now figure out rU
# 
# ID$tr_uhx


rHH <- rpi2u(rpivec=ID$r_hh[-32],"H","H")
rUH <- rpi2u(rpivec=ID$r_uh[-32],"U","H")
rUU <- rpi2u(rpivec=ID$r_uu[-32],"U","U")
rHU <- rpi2u(rpivec=ID$r_hu[-32],"H","U")
# does this need to change?
rU <- u2U(rHH,rHU,rUH,rUU)

0.9158685 # rev prev 1
0.9130396 # chrono prev 1

# ------------------------------------ #
# reverse prevalence convergence       #
# ------------------------------------ #
rU[,64]
pop <- matrix(0,ncol=33,nrow=64)
pop[1:32,1] <- 1
for (i in 1:32){
  pop[,i+1] <- rU %*% pop[,i]
}

PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
rpi1 <- PH / (PH + PU)

pop <- matrix(0,ncol=33,nrow=64)
pop[33:64,1] <- 1
for (i in 1:32){
  pop[,i+1] <-  rU %*% pop[,i]
}
PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
rpi2 <- PH / (PH + PU)

a2 <- seq(50,112,by=2)
plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1),
     main = "not yet identical to chrono convergence")

for (i in 1:32){
  lines(a2[1:(33-i)],rev(rpi1[row(rpi1) + col(rpi1) == (34 - i)]), col = "#FF000080")
  lines(a2[1:(33-i)],rev(rpi2[row(rpi2) + col(rpi2) == (34 - i)]), col = "#0000FF80")
}
stationary_init <- ph / (ph + pu)
lines(a2, stationary_init, lty = 2,lwd=2)
check <- c(rep(0,31),1,rep(0,32))

rU%*% (rU %*% check)


# Test shows convergence toward a different prevalence.
# so left with the question of how to define the backward 
# probabilities properly
pop <- matrix(0,ncol=33,nrow=64)
pop[,1] <- c(stationary_init,1-stationary_init)
for (i in 1:32){
  pop[,i+1] <- rU %*% pop[,i]
}

PH <- pop[1:32,1:32]
PU <- pop[33:64,1:32]
rpitest <- PH / (PH + PU)

plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1),
     main = "not yet identical to chrono convergence")

for (i in 1:32){
  lines(a2[1:(33-i)],rev(rpitest[row(rpitest) + col(rpitest) == (34 - i)]), col = "#FF000080")
}
lines(a2, stationary_init, lty = 2,lwd=2)

rpitest

stationary_init
