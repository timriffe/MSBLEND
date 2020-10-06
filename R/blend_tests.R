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
                       edu == "all_edu",
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
# these proportions need to come from somewhere.
N[,c("H::48","U::48")] %*% diag(c(0.9190056, 1-0.9190056)) %>% sum()
N[,c("H::48","U::48")] %*% diag(init) %>% sum()
N[,c("H::48","U::48")] %*% diag(init_constant(TRsub)) %>% sum() 

N[,c("H::48","U::48")] %*% diag(c(1,0)) %>% sum() 
N[,c("H::48","U::48")] %*% diag(c(0,1)) %>% sum() 

x <- seq(48,110,by=2)
hsx <- N[1:32,"H::48"]+N[33:64,"H::48"]
usx <- N[1:32,"U::48"]+N[33:64,"U::48"]
hin <- N[1:32,"H::48"] / (N[1:32,"H::48"]+N[33:64,"H::48"])
uin <- N[1:32,"U::48"] / (N[1:32,"U::48"]+N[33:64,"U::48"])
plot(x, hin, type= 'l')
lines(x, uin)
polygon(c(x,rev(x)),c(hin,rev(uin)),col = "#9E87BB")

plot(x,hsx/2,type = 'l')
lines(x,usx/2)
sum(hsx) - sum(usx)



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
rpi       <- proj_prev(c(0.9190056,rep(0,31)), U) %>% diag()

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
         time_r = time * ifelse(state_from == "H", 0.9190056, 1-0.9190056),
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
         U_D = m24) %>% 
  pivot_longer(cols = H_H:U_D,
               names_to = c("state_from", "state_to"),
               names_sep = "_",
               values_to = "prob") %>% 
  # ensure exact composition, possibly redundant
  group_by(state_from, age) %>% 
  mutate(prob = prob / sum(prob)) %>% 
  filter(state_from != "D") %>% 
  # previosuly also removes state_to = "D", but testing here.
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
  pivot_wider(names_from = c(state_from, state_to), 
              id_cols = age,
              values_from = c(prob,time,tr)) %>% 
  select(age, 
         Hx = time_H_H, 
         Ux = time_U_U, 
         tr_H_H, 
         tr_U_H, 
         tr_H_D, 
         tr_U_D, 
         tr_H_U, 
         tr_U_U) %>% 
  mutate(r_hh = tr_H_H / lead(Hx),
         r_hu = tr_U_H / lead(Hx),
         r_dh = tr_H_D / (tr_H_D + tr_U_D),   # test
         r_uu = tr_U_U / lead(Ux),
         r_uh = tr_H_U / lead(Ux),
         r_du = tr_U_D / (tr_H_D + tr_U_D) )  # test
# test
  

# need new notation. (in)_(from)
# build out increment-decrement LT

ID <- 
  TRsub %>% 
  IDLT(init = c(0.9190056, 1-0.9190056)) %>% 
  mutate(CDx = 2 - (Hx + Ux),
         Dx = diff(c(0,CDx)), # test
         tr_hhx = hhx * Hx,
         tr_hux = hux * Hx,
         tr_hdx = hdx * Hx,
         tr_uux = uux * Ux,
         tr_uhx = uhx * Ux,
         tr_udx = udx * Ux,
         # reverse probabilities: 
         # (note interstate transfer directions)
         r_hh = tr_hhx / lead(Hx) ,
         r_hu = tr_uhx / lead(Hx),
         r_dh = tr_hdx / lead(Dx), # test
         r_uu = tr_uux / lead(Ux),
         r_uh = tr_hux / lead(Ux),
         r_du = tr_udx / lead(Dx),
         # experiments:
         r_hh = tr_hhx / lead(Hx) ,
         r_hu = tr_uhx / lead(Hx),
         r_dh = tr_hdx / lead(Dx), # test
         r_uu = tr_uux / lead(Ux),
         r_uh = tr_hux / lead(Ux),
         r_du = tr_udx / lead(Dx),
         )


# check
# ID$tr_hhx[-32] + ID$tr_uhx[-32] - ID$Hx[-1]
# 
# # now figure out rU
# 
# ID$tr_uhx

# build back prob projection matrices
rHH <- rpi2u(rpivec = ID$r_hh[-32], "H", "H")
rUH <- rpi2u(rpivec = ID$r_uh[-32], "U", "H")
rDH <- rpi2u(rpivec = ID$r_dh[-32], "D", "H")
rUU <- rpi2u(rpivec = ID$r_uu[-32], "U", "U")
rHU <- rpi2u(rpivec = ID$r_hu[-32], "H", "U")
rDU <- rpi2u(rpivec = ID$r_du[-32], "D", "U")

# does this need to change?
rU <- u2U(rHH,rHU,rUH,rUU)

# maybe should be -1 ?
rHH <- rpi2u(rpivec=back_prob$r_hh[-32],"H","H")
rUH <- rpi2u(rpivec=back_prob$r_uh[-32],"U","H")
rDH <- rpi2u(rpivec=back_prob$r_dh[-32],"D","H")
rUU <- rpi2u(rpivec=back_prob$r_uu[-32],"U","U")
rHU <- rpi2u(rpivec=back_prob$r_hu[-32],"H","U")
rDU <- rpi2u(rpivec=back_prob$r_du[-32],"D","U")

# does this need to change?
# rU2 <- u2U(rHH,rHU,rUH,rUU)
# dim(rHH)

# Include resucisitation origins:
rU2 <- rbind(
  cbind(rHH, rUH, rDH),
  cbind(rHU, rUU, rDU),
  matrix(0,nrow=32,ncol=32*3))
rU2[is.na(rU2)] <- 0

# rDH <- rpi2u(rpivec=back_prob$r_dh[-32],"D","H")
# rDU <- rpi2u(rpivec=back_prob$r_du[-32],"D","U")
# 
# rU3 <- rbind(cbind(rHH, rUH, rDH),
#       cbind(rHU, rUU, rDU))

0.9190056 # rev prev 1
0.91586850 # chrono prev 1

# just starts with 1s.


# ------------------------------------ #
# reverse prevalence convergence       #
# ------------------------------------ #
# pop    <- proj_pop(initH = rep(1, 32),rU)
# 
# rpi1   <- proj_prev(initH = rep(1, 32),rU, TRUE)
# rpi0   <- proj_prev(initH = rep(0, 32),rU, TRUE)

rpiH <- proj_prev_back(rU2,"H")
rpiU <- proj_prev_back(rU2,"U")
rpiD <- proj_prev_back(rU2,"D")
# these need to be different. How many deaths do we add in?
# stationary deaths per age? Interesting catch 22. Deaths per age
# requires an initial mixture! The initial mixture requires deaths
# per age! So looks like an iterative solution would be nice? i.e.
# these back probabilities are just a first approximation! 

# 
a2     <- seq(50,112,by=2)
plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1),
     main = "back-convergence")

for (i in 1:32){
  lines(a2[1:(33-i)],rev(rpiH[row(rpiH) + col(rpiH) == (34 - i)]), col = "#FF000080")
  lines(a2[1:(33-i)],rev(rpiU[row(rpiU) + col(rpiU) == (34 - i)]), col = "#0000FF80")
  # lines(a2[1:(33-i)],rev(rpiD[row(rpiD) + col(rpiD) == (34 - i)]), col = "#00000080")
}


lines(a2, rpi,lwd=2)


rpiH 

back_prob_int <- function(TRsub, forward_init = c(.9,.1)){
  
  # produce increment-decrement lifetable
  # derive transfers and back probs
  ID <- 
    TRsub %>% 
    IDLT(init = forward_init) %>% 
    mutate(CDx = 2 - (Hx + Ux),
           Dx = diff(c(0,CDx)), # test
           tr_hhx = hhx * Hx,
           tr_hux = hux * Hx,
           tr_hdx = hdx * Hx,
           tr_uux = uux * Ux,
           tr_uhx = uhx * Ux,
           tr_udx = udx * Ux,
           # reverse probabilities: 
           # (note interstate transfer directions)
           r_hh = tr_hhx / lead(Hx) ,
           r_hu = tr_uhx / lead(Hx),
           r_dh = tr_hdx / lead(Dx), # test
           r_uu = tr_uux / lead(Ux),
           r_uh = tr_hux / lead(Ux),
           r_du = tr_udx / lead(Dx),
           # experiments:
           r_hh = tr_hhx / lead(Hx) ,
           r_hu = tr_uhx / lead(Hx),
           r_dh = tr_hdx / lead(Dx), # test
           r_uu = tr_uux / lead(Ux),
           r_uh = tr_hux / lead(Ux),
           r_du = tr_udx / lead(Dx),
    )
  
  # create projection blocks
  rHH <- rpi2u(rpivec = ID$r_hh[-32], "H", "H")
  rUH <- rpi2u(rpivec = ID$r_uh[-32], "U", "H")
  rDH <- rpi2u(rpivec = ID$r_dh[-32], "D", "H")
  rUU <- rpi2u(rpivec = ID$r_uu[-32], "U", "U")
  rHU <- rpi2u(rpivec = ID$r_hu[-32], "H", "U")
  rDU <- rpi2u(rpivec = ID$r_du[-32], "D", "U")
  
  # compose reverse projection matrix
  rU2 <- rbind(
    cbind(rHH, rUH, rDH),
    cbind(rHU, rUU, rDU),
    matrix(0,nrow=32,ncol=32*3))
  rU2[is.na(rU2)] <- 0
  
  # back-project prevalence (H) from H, U, D
  rpiH <- proj_prev_back(rU2,"H")
  rpiU <- proj_prev_back(rU2,"U")
  rpiD <- proj_prev_back(rU2,"D")
  
  # They all agree, but just in case, take the converged-upon
  # prevalence healthy as the mean of those from each starting state.
  (rpiH[1,ncol(rpiH)-1] + 
    rpiU[1,ncol(rpiU)-1] + 
    rpiD[1,ncol(rpiD)-1]) / 3
}


# rpi1_2 <- proj_prev(initH = rep(1, 32),rU2, TRUE)
# rpi0_2 <- proj_prev(initH = rep(0, 32),rU2, TRUE)
# 
# a2     <- seq(50,112,by=2)
# plot(NULL, type = 'n',xlim = c(50,110),ylim=c(0,1),
#      main = "not yet identical to chrono convergence")
# 
# for (i in 1:32){
#   lines(a2[1:(33-i)],rev(rpi1[row(rpi1) + col(rpi1) == (34 - i)]), col = "#FF000080")
#   lines(a2[1:(33-i)],rev(rpi0[row(rpi0) + col(rpi0) == (34 - i)]), col = "#0000FF80")
# }
# lines(a2, rpi,lwd=2)


# it seems these back-probability calculations always give back more
# optimistic starting conditions than what they start with, asymptotically
# leading to everyone in good health at the start...
N <- 200
init_i <- c(.1,.9)

init_mat<- matrix(NA,nrow=N+1,ncol=2)
init_mat[1, ] <- init_i
for (i in 1:N){
  init_mat[i+1, ] <- init_back_prob(TRsub, forward_init = init_mat[i, ])
}
init_mat
plot(init_mat[,1])

