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

N <- U %>% 
  U2N()
N[,c("H::48","U::48")] %*% diag(c(.9158685, 1-.9158685)) %>% rowSums()

0.0000336526
# ------------------------------------ #
# chronologocal prevalence convergence #
# ------------------------------------ #
pi1 <- proj_prev(initH = rep(1, 32),U)
pi0 <- proj_prev(initH = rep(0, 32),U)
a2  <- seq(50, 112 ,by = 2)
plot(NULL, type = 'n', xlim = c(50, 110), ylim = c(0, 1))
for (i in 1:32){
  lines(a2[i:32], pi1[row(pi1) == (col(pi1) - 1 + i)], col = "#FF000080")
  lines(a2[i:32], pi0[row(pi0) == (col(pi0) - 1 + i)], col = "#0000FF80")
}

# --------------------------------------------------- #
# alternative starting prevs                          #
# --------------------------------------------------- #
phoptim   <- optimize(interval = c(.7,1), f=init_min, U = U)$min
pi_optim  <- proj_prev(c(phoptim,rep(0,31)), U) %>% diag()

init_it   <- init_constant(TRsub)
pi_it     <- proj_prev(c(init_it[1],rep(0,31)), U) %>% diag()
# again using the back probability-derived starting prop
rpi       <- proj_prev(c(0.9158685,rep(0,31)), U) %>% diag()

plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1))
for (i in 1:32){
  lines(a2[i:32],pi1[row(pi1) == (col(pi1) -1 + i)], col = "#FF000080")
  lines(a2[i:32],pi0[row(pi0) == (col(pi0) -1 + i)], col = "#0000FF80")
}
lines(a2, pi_it, lty = 1,lwd=1,col="magenta")
lines(a2, rpi, lty = 2,lwd=1)
lines(a2, pi_optim)
# this works better

# --------------------------------------------------- #

time_to <- 
  U %>% 
  U2N() %>% 
  melt(varnames = c("to","from"),
       value.name = "time") %>% 
  mutate(to = as.character(to),
         from = as.character(from)) %>% 
  separate(col = "to", 
           sep = "::",
           into = c("state_to","age"),
           convert = TRUE) %>% 
  separate(col = "from", 
           sep = "::",
           into = c("state_from","age_from"),
           convert = TRUE) %>% 
  filter(age >= age_from,
         age_from == 48,
         state_from != "D") %>%
  mutate(state_to = case_when(
    state_to == "UUD" & state_from == "U" ~ "U",
    state_to == "HHD" & state_from == "H" ~ "H",
    TRUE ~ state_to),
    age = ifelse(is.infinite(age), 110, age)
  ) %>% 
  filter(state_to %in% c("U","H")) %>% 
  group_by(state_from, age) %>% 
  mutate(time_it = time * ifelse(state_from == "H", init_it[1], init_it[2]),
         time_sprop = time *  ifelse(state_from == "H", init[1], init[2]),
         time_r = time * ifelse(state_from == "H", .9158685, 1-.9158685),
         time_test = time * ifelse(state_from == "H", .8, .2)) %>% 
  ungroup() %>% 
  group_by(state_to,age) %>% 
  summarize(time_it = sum(time_it),
            time_sprop = sum(time_sprop),
            time_r = sum(time_r),
            time_test = sum(time_test)) %>% 
  arrange(state_to, age)

# take peek at occupancy curves
time_to %>% 
  ggplot(aes(x = age, y = time_sprop, color = state_to)) +
  geom_line() +
  geom_line(aes(y = time_it), linetype = "dashed") +
  geom_line(aes(y = time_r), linetype = "dotted",size=1.5) +
  geom_line(aes(y = time_test), linetype = "dotted",size=1.5)
  

# tidy up the input probabilities
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
  filter(state_from != "D",
         state_to != "D") %>% 
  # if we end up doing this with all subsets then adjust here to keep year, sex, edu
  select(age, state_from, state_to, prob) %>% 
  ungroup() %>% 
  mutate(age = as.integer(age - 2))

# compare expectancies:
time_to %>% 
  group_by(state_to) %>% 
  summarize(time_r = sum(time_r),
            time_sprop = sum(time_sprop),
            time_it = sum(time_it),
            time_test = sum(time_test))

# calculate transfers, and use these to get back probabilities
back_prob <-  
  time_to %>% 
  rename(state_from = state_to) %>% 
  left_join(tidysub) %>% 
  mutate(tr = prob * time_r,
         state_to = ifelse(is.na(state_to), state_from, state_to)) %>% 
  select(state_to, state_from, age, time_r, prob, tr) %>% 
  arrange(age, state_to)  %>% 
  rename(time = time_r) %>% 
  pivot_wider(names_from = c(state_from, state_to), id_cols = age,
              values_from = c(prob,time,tr)) %>% 
  select(age, Hx = time_H_H, Ux = time_U_U, tr_H_H, tr_U_H, tr_H_U, tr_U_U) %>% 
  mutate(r_hh = tr_H_H / lead(Hx),
         r_hu = tr_U_H / lead(Hx),
         r_uu = tr_U_U / lead(Ux),
         r_uh = tr_H_U / lead(Ux)) 
  
tail(back_prob)
tail(ID)
# need new notation. (in)_(from)
# build out increment-decrement LT

ID <- 
  TRsub %>% 
  IDLT(init = c(.9158685, 1-.9158685)) %>% 
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

# build back prob projection matrices
rHH <- rpi2u(rpivec=ID$r_hh[-32],"H","H")
rUH <- rpi2u(rpivec=ID$r_uh[-32],"U","H")
rUU <- rpi2u(rpivec=ID$r_uu[-32],"U","U")
rHU <- rpi2u(rpivec=ID$r_hu[-32],"H","U")
# does this need to change?
rU <- u2U(rHH,rHU,rUH,rUU)


rHH <- rpi2u(rpivec=back_prob$r_hh[-32],"H","H")
rUH <- rpi2u(rpivec=back_prob$r_uh[-32],"U","H")
rUU <- rpi2u(rpivec=back_prob$r_uu[-32],"U","U")
rHU <- rpi2u(rpivec=back_prob$r_hu[-32],"H","U")
# does this need to change?
rU2 <- u2U(rHH,rHU,rUH,rUU)

0.9158685 # rev prev 1
0.9130396 # chrono prev 1

# ------------------------------------ #
# reverse prevalence convergence       #
# ------------------------------------ #
pop <- proj_pop(initH = rep(1, 32),rU)

rpi1 <- proj_prev(initH = rep(1, 32),rU, TRUE)
rpi0 <- proj_prev(initH = rep(0, 32),rU, TRUE)

rpi1_2 <- proj_prev(initH = rep(1, 32),rU2, TRUE)
rpi0_2 <- proj_prev(initH = rep(0, 32),rU2, TRUE)

a2 <- seq(50,112,by=2)
plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1),
     main = "not yet identical to chrono convergence")

for (i in 1:32){
  lines(a2[1:(33-i)],rev(rpi1[row(rpi1) + col(rpi1) == (34 - i)]), col = "#FF000080")
  lines(a2[1:(33-i)],rev(rpi0[row(rpi0) + col(rpi0) == (34 - i)]), col = "#0000FF80")
}
lines(a2, rpi,lwd=2)


