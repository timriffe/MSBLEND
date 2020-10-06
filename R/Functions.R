# Function's we'll use

# Make a single U submatrix from a pi (transfer probs) vector
pi2u <- function(pivec, 
                 from ="H",
                 to = "H",
                 start_age = 50,
                 interval = 2) {
  out           <- cbind(rbind(0, diag(pivec)), 0)
  n             <- length(pivec)
  # the final subtraction of the interval is particular to
  # the way these probabilities were estimated and labelled.
  # to technically our first one goes from 48 to 50, not from 50 to 52.
  ages          <- ((0:n) * interval) + start_age - interval
  from_names    <- c(paste(from,ages[-length(ages)],sep="::"),"D::Inf")
  # to_names      <-c(paste(to,ages[-1],sep="::"),"D::Inf")
  # TR: double checking alignment of age
  to_names      <- c(paste(to,ages[-length(ages)],sep="::"),paste0(from,to,"D::Inf"))
  dimnames(out) <- list(to_names, from_names)
  out
}

# same thing, but assuming back probabilities are supplied
rpi2u <- function(rpivec, 
                  from ="H",
                  to = "H",
                  start_age = 50,
                  interval = 2) {
  out           <- cbind(0,rbind(diag(rpivec),0))
  n             <- length(rpivec)
  # the final subtraction of the interval is particular to
  # the way these probabilities were estimated and labelled.
  # to technically our first one goes from 48 to 50, not from 50 to 52.
  ages          <- ((0:n) * interval) + start_age - interval
  from_names    <- c(paste(from,ages,sep="::"))
  # to_names      <-c(paste(to,ages[-1],sep="::"),"D::Inf")
  # TR: double checking alignment of age
  to_names      <- c(paste(to,ages,sep="::"))
  dimnames(out) <- list(to_names, from_names)
  out
}
# Compose u blocks into U
u2U <- function(HH, HU, UH, UU){
  rbind(
    cbind(HH, UH),
    cbind(HU, UU))
}

# convert transient dynamics into outcomes: the fundamental matrix, N
U2N <- function(U, interval = 2) {
  I   <- diag(nrow(U))
  Nsx <- solve(I - U) * interval
  dimnames(Nsx) <- dimnames(U)
  Nsx
}

# create U from data subset
sub2U <- function(X, start_age = 50, interval = 2){
  HH <- pi2u(pivec = X[,"m11"], 
             from = "H", 
             to = "H", 
             start_age = start_age, 
             interval = interval)
  HU <- pi2u(pivec = X[,"m12"], 
             from = "H", 
             to = "U", 
             start_age = start_age, 
             interval = interval)
  UH <- pi2u(pivec = X[,"m21"], 
             from = "U", 
             to = "H", 
             start_age = start_age, 
             interval = interval)
  UU <- pi2u(pivec = X[,"m22"], 
             from = "U", 
             to = "U", 
             start_age = start_age, 
             interval = interval)
  
  U <- u2U(HH = HH, # healthy to healthy
           HU = HU, # healthy to unhealthy
           UH = UH, # unhealthy to healthy
           UU = UU) # unhealthy to unhealthy
}

# --------------------------------------------------- #
# optimized starting prevalence                       #
# --------------------------------------------------- #

# pretty good, but not optimal. Really, should just optimize on first jump.

# project forward
proj_pop <- function(initH = rep(1,32), U){
  N <- length(initH)
  init <- c(initH,1-initH)
  
  
  pop <- matrix(0, ncol= N+1,nrow = N*2)
  pop[, 1]  <- init
  for (i in 1:N){
    pop[,i+1] <-  U %*% pop[,i]
  }
  pop
}

# get back projection matrix of prevalences, in diagonals,
# direction depends on whether forward or backward.
proj_prev <- function(initH, U, r = FALSE){
  N   <- length(initH)
  pop <- proj_pop(initH,U)
  PH  <- pop[1:N,1:N]
  PU  <- pop[(N+1):(2*N),1:N]
  
  if (!r){
  PH[upper.tri(PH)] <- NA
  PU[upper.tri(PU)] <- NA
  }
  pi  <- PH / (PH + PU)
  
  pi
}

# do a proj and pick out just the prevalence
init_pop <- function(initH, U){
  N   <- length(initH)
  pop <- proj_pop(initH, U = U)
  PH  <- pop[1:N, 1:N]
  PU  <- pop[(N+1):(2*N), 1:N]
  
  ph  <- diag(PH)
  pu  <- diag(PU)
  
  list(ph = ph, pu = pu)
}

init_min <- function(par =.9, U){
  p  <- init_pop(init = par, U = U)
  ph <- p$ph
  pu <- p$pu
  # really will only affect first part of curve.
  sum(abs(diff(ph / (ph + pu))))
}

# based on assumption of frozen probs before first age group
init_constant <- function(TRsub){
  u <- matrix(c(TRsub[1,"m11"],TRsub[1,"m12"],TRsub[1,"m21"],TRsub[1,"m22"]),2)
  v <- eigen(u)$vectors[,1]
  v / sum(v)
}

# increment decrement lifetable, two states only.
# dat shoul dhave columns m11, m12, m22, m21, age
IDLT <- function(dat, init, interval = 2){
  n  <- nrow(dat)
  Hx <- rep(0, n+1)
  Ux <- rep(0, n+1)
  
  hhx <- dat %>% pull(m11)
  hux <- dat %>% pull(m12)
  uux <- dat %>% pull(m22)
  uhx <- dat %>% pull(m21)
  
  hdx <- dat %>% pull(m14)
  udx <- dat %>% pull(m24)
  # if not given then assume constant.
  if (missing(init)){
    u1   <- matrix(c(hhx[1],hux[1],uhx[1],uux[1]),2)
    v1   <- eigen(u1)$vectors[,1]
    init <- v1 / sum(v1)
  }
  #cat(init)
  Hx[1] <- init[1] * interval
  Ux[1] <- init[2] * interval
  
  for (i in 1:n){
    Hx[i+1] <- Hx[i] * hhx[i] + Ux[i] * uhx[i]
    Ux[i+1] <- Ux[i] * uux[i] + Hx[i] * hux[i]
  }
  ages <- c(min(dat$age)-interval,dat$age)
  data.frame(age = ages,
             Hx=Hx,
             Ux=Ux,
             hhx=c(hhx, 0),
             hux=c(hux, 0),
             hdx=c(hdx, 0),
             uux=c(uux, 0),
             uhx=c(uhx, 0),
             udx=c(udx, 0))
}
proj_pop_back <- function(rU, orig = "H"){
  N     <- nrow(rU)
  n     <- N / 3
  ones  <- rep(1, n)
  zeros <- ones * 0
  if (orig == "H") init <- c(ones, zeros, zeros)
  if (orig == "U") init <- c(zeros, ones, zeros)
  if (orig == "D") init = c(zeros, zeros, ones)
  
  pop <- matrix(0, 
                ncol = n + 1, 
                nrow = N)
  pop[, 1]  <- init
  for (i in 1:n){
    pop[, i + 1] <-  rU %*% pop[, i]
  }
  pop
}

proj_prev_back <- function(rU, orig = "H"){
  pop   <- proj_pop_back(rU, orig = orig)
  N     <- nrow(rU)
  n     <- N / 3
  # now get proportion healthy
  PH    <- pop[1:n,1:n]
  PU    <- pop[(n + 1):(2 * n), 1:n]
  
  pi    <- PH / (PH + PU)
  
  pi
}

logit <- function(x){log(x / (1 - x))}
expit <-  function(x){exp(x) / (1 + exp(x))}
# back-and-forth projected starting prevalence:
# this drop TRsub to a lower age 
# using some constrained linear assumptions
down_proj <- function(TRsub, 
                      start_age = 40, 
                      max_age_fit = 65, 
                      interval = 2){
  # extract pieces, easier that way
  m14     <- TRsub$m14
  m11     <- TRsub$m11
  m12     <- TRsub$m12
  
  m24     <- TRsub$m24
  m22     <- TRsub$m22
  m21     <- TRsub$m21
  
  age     <- TRsub$age
  
  fit_i   <- age <= max_age_fit
  age_fit <- age[fit_i]
  # 1) healthy mortality down project
  mod14    <- lm(log(m14[fit_i]) ~ age_fit)
  age_new <- seq(start_age, min(age) - interval, by = interval)
  m14new  <- exp(predict(mod14, 
                         newdata = data.frame(age_fit = age_new)))
  
  # 2) down extrap relation between conditional bla bla
  
  # plot(age, m12 / m11, log = 'y')
  # plot(age, m12 / (m12 + m11), log = 'y')
  m12frac  <- m12[fit_i] / (m12[fit_i] + m11[fit_i])
  
  mod12    <- lm(log(m12frac) ~ age_fit)
  m12_frac_new   <- exp(predict(mod12, 
                                newdata = data.frame(age_fit = age_new)))
  
  m12new <- m12_frac_new * (1 - m14new)
  m11new <- (1 - m14new - m12new)
  
  # 3) unhealthy mortality down project
  mod24   <- lm(log(m24[fit_i]) ~ age_fit)
  age_new <- seq(start_age, min(age) - interval, by = interval)
  m24new  <- exp(predict(mod24, 
                         newdata = data.frame(age_fit = age_new)))
  
  # 2) down extrap relation between conditional bla bla
  
  # plot(age, m21, log = 'y')
  # plot(age, m21 / (m22 + m21))
  # plot(age, logit(m22 / (m22 + m21)))
  m22fraclogit <-  logit(m22[fit_i] / (m22[fit_i] + m21[fit_i])) 
  
  mod22    <- lm(m22fraclogit ~ age_fit)
  m22frac_new   <- expit(predict(mod22, 
                                 newdata = data.frame(age_fit = age_new)))
  
  m22new <- m22frac_new * (1 - m24new)
  m21new <- (1 - m24new - m22new)
  
  # append
  data.frame(
    age = c(age_new, age),
    m11 = c(m11new, m11),
    m12 = c(m12new, m12),
    m14 = c(m14new, m14),
    m22 = c(m22new, m22),
    m21 = c(m21new, m21),
    m24 = c(m24new, m24))
}

# Get a converged-upon starting prevalence based on down-projected
# transtion probabilities.
init_back_proj <- function(
  TRsub, 
  start_age = 40, 
  max_age_fit = 65, 
  interval = 2){
  TRsub_proj <- 
    TRsub %>% 
    down_proj()
  Uproj <- sub2U(TRsub_proj, 
                 start_age = start_age, 
                 interval = interval)
  Nproj <- U2N(Uproj)
  
  # an arbitrary, but close initial mixture
  init_it_const <- init_constant(TRsub = TRsub_proj)
  
  Stocks <- 
    Nproj[, c("H::38","U::38")] %*% 
    diag(init_it_const) %>% 
    rowSums() 
  
  H <- Stocks["H::48"] / (Stocks["H::48"] + Stocks["U::48"])
  c(H, 1 - H)
}

# starting prevalence from back-probabilities:
init_back_prob <- function(TRsub, forward_init = c(.9,.1)){
  
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
  out <- (rpiH[1,ncol(rpiH)-1] + 
      rpiU[1,ncol(rpiU)-1] + 
      rpiD[1,ncol(rpiD)-1]) / 3
  c(out, 1-out)
}




# ----------------------------------------------------
# turnover experiments:

# turnover functions
plot_turnover <- function(TRsub, scale = 1){
  x     <- rescale_turnover(TRsub[,c("m11","m12","m14")],scale)[,"m12"]
  y     <- rescale_turnover(TRsub[,c("m22","m21","m24")],scale)[,"m21"]
  age   <- TRsub$age
  age10 <- age %% 10 == 0
  plot(
    x, 
    y,
    asp = 1,
    xlab = "onset",
    ylab = "recovery",
    las = 1,
    type = 'l')
  points(x[age10],y[age10],pch=16)
  text(x[age10],y[age10],labels = age[age10],pos=4)
  abline(a = 0, b = 1)
}

rescale_turnover <- function(mat,scale=.5){
  to_d             <- mat[,3]
  
  # if scale < 1 then we scale DOWN moving.
  if (scale <= 1){
    # conditional prob of switching states
    cond_move        <- mat[,2] / rowSums(mat[,1:2])
    cond_move_scaled <- cond_move * scale
    out              <- cbind((1-to_d) * (1-cond_move_scaled),
                              (1-to_d) * (cond_move_scaled),
                              to_d)
  }
  # if scale > 1 then we scale DOWN staying.
  if (scale > 1){
    cond_stay        <- mat[,1] / rowSums(mat[,1:2])
    cond_stay_scaled <- cond_stay * (1 / scale)
    out              <- cbind((1-to_d) * (cond_stay_scaled),
                              (1-to_d) * (1 - cond_stay_scaled),
                              to_d)
  }
  colnames(out) <- colnames(mat)
  out
}

# harmonic_mean_vec <- function(x,y){
#   2 / (1/x + 1/y)
# }
# geometric_mean_vec <- function(x,y){
#   sqrt(x*y)
# }
# 
# rescale_turnover <- function(mat,
#                              scale=.5,
#                              mean_type = "harmonic"){
#     to_d             <- mat[,3]
#   
#   # if scale < 1 then we scale DOWN moving.
#  
#     # conditional prob of switching states
#     cond_move        <- mat[,2] / rowSums(mat[,1:2])
#     cond_stay        <- mat[,1] / rowSums(mat[,1:2])
#     
#     
#     if (scale <= 1){
#       cond_move_scaled <- logit(cond_move) * scale
#     }
#     cond_move_scaled <- cond_move * scale
#     cond_stay_scaled <- cond_stay * (1 / scale)
# 
#     
#     
#     # cond_move_scaled[cond_move_scaled > 1] <- 1
#     # cond_move_scaled[cond_move_scaled < 0] <- 0
#     # 
#     # cond_stay_scaled[cond_stay_scaled > 1] <- 1
#     # cond_stay_scaled[cond_stay_scaled < 0] <- 0
#     # if (mean_type == "harmonic"){
#     #   cond_move_hm     <- harmonic_mean_vec(logit(cond_move_scaled), 
#     #                                         logit(1 - cond_stay_scaled))
#     # }
#     # if (mean_type == "geometric"){
#     #   cond_move_hm     <- geometric_mean_vec(logit(cond_move_scaled), 
#     #                                          logit(1 - cond_stay_scaled))
#     # }
#     cond_stay_hm     <- 1 - cond_move_hm 
#     out              <- cbind((1-to_d) * cond_stay_hm,
#                               (1-to_d) * cond_move_hm,
#                               to_d)
#     
#   colnames(out) <- colnames(mat)
#   out
# }
# 

rescale_turnover_U <- function(mat, scale){
  mat[,c("m11","m12","m14")] <- 
    rescale_turnover(mat[,c("m11","m12","m14")],scale)
  mat[,c("m22","m21","m24")] <- 
    rescale_turnover(mat[,c("m22","m21","m24")],scale)
  sub2U(mat)
}
