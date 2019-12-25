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
  to_names      <- c(paste(to,ages[-length(ages)],sep="::"),"D::Inf")
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
sub2U <- function(X){
  HH <- pi2u(pivec = TRsub[,"m11"], from = "H", to = "H")
  HU <- pi2u(pivec = TRsub[,"m12"], from = "H", to = "U")
  UH <- pi2u(pivec = TRsub[,"m21"], from = "U", to = "H")
  UU <- pi2u(pivec = TRsub[,"m22"], from = "U", to = "U")
  
  U <- u2U(HH = HH, # healthy to healthy
           HU = HU, # healthy to unhealthy
           UH = UH, # unhealthy to healthy
           UU = UU) # unhealthy to unhealthy
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


# increment decrement lifetable, two states only.
IDLT <- function(dat, init, interval = 2){
  n <- nrow(dat)
  Hx <- rep(0, n+1)
  Ux <- rep(0, n+1)
  
  hhx <- dat %>% pull(m11)
  hux <- dat %>% pull(m12)
  uux <- dat %>% pull(m22)
  uhx <- dat %>% pull(m21)
  
  if (missing(init)){
    u1 <- matrix(c(hhx[1],hux[1],uhx[1],uux[1]),2)
    init <- c(.5,.5)
    for (i in 1:30){
      init <- u1 %*% init
      init <- init / sum(init)
    }
  }
  #cat(init)
  Hx[1] <- init[1] * interval
  Ux[1] <- init[2] * interval
  
  for (i in 1:n){
    Hx[i+1] <- Hx[i] * hhx[i] + Ux[i] * uhx[i]
    Ux[i+1] <- Ux[i] * uux[i] + Hx[i] * hux[i]
  }
  
  data.frame(age= c(48,dat$age),Hx=Hx,Ux=Ux,hhx=c(hhx,0),hux=c(hux,0),uux=c(uux,0),uhx=c(uhx,0))
}
