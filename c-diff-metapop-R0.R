library(tidyverse)

# model without vaccination or risk
R_0_cdiff_metapop_no_vacc_no_risk <- function(S0, parameters) {
  with(as.list(c(S0, parameters)),
       {
         BC <- matrix(
           c(betaC1*S1, 0,
             0, betaC2*S2),
           nrow = 2, byrow = TRUE
         )
         
         BD <- matrix(
           c(betaD1*S1, 0,
             0, betaD2*S2),
           nrow = 2, byrow = TRUE
         )
         
         d <- (delta+d2)*N2/N1
         TC1 <- d1+xi1+phi1+d
         TC2 <- d2+xi2+phi2+delta
         TD1 <- d1+p1*eps1+d
         TD2 <- d2+p2*eps2+delta
         
         TC <- matrix(
           c(TC1, -delta,
             -d, TC2),
           nrow = 2, byrow = TRUE
         )
         
         TD <- matrix(
           c(TD1, -delta,
             -d, TD2),
           nrow = 2, byrow = TRUE
         )
         
         PHI <- diag(c(phi1, phi2), 2, 2)
         FVinv <- (BC + BD %*% solve(TD) %*% PHI) %*% solve(TC)
         spectral_radius(FVinv)
       })
}

get_S0 <- function(parameters) {
  with(as.list(c(parameters)),
       {
         d <- (delta+d2)*N2/N1
         TS1 <- d1+theta1+d+alpha1
         TS2 <- d2+theta2+delta+alpha2
         S <- matrix(
           c(TS1, -delta,
             -d, TS2),
           nrow = 2, byrow = TRUE
         )
         b <- matrix(
           c(alpha1*N1, alpha2*N2),
           nrow = 2, byrow = TRUE
         )
         S0 <- solve(S) %*% b
         S0 <- S0[c(1,2)]
         names(S0) <- c("S1", "S2")
         S0
       })
}

spectral_radius <- function(X) {
  eigs <- eigen(X, only.values = TRUE)$values
  eigs[which.max(abs(eigs))]
}

parameters <- c(N1 = 100000, N2 = 2000, d1 = 1/(78.5*365), d2 = 0.0068,
                
                alpha1 = 0.01, theta1 = 0.033, xi1 = 0.0165, phi1 = 0.2, 
                p1 = 0.8, eps1 = 0.1, betaC1 = 0.007/50, betaD1 = 0.007/50,
                
                alpha2 = 0.5, theta2 = 0.033, xi2 = 0.0165, phi2 = 0.2, 
                p2 = 0.8, eps2 = 0.1, betaC2 = 0.007, betaD2 = 0.007,
                
                delta = 0.135)

S0 <- get_S0(parameters)

R_0_cdiff_metapop_no_vacc_no_risk(S0, parameters)

# model without vaccination
get_R0 <- function(S0, parameters) {
  with(as.list(c(S0, parameters)),
       {
         BC <- matrix(
           c(betaC1L*S1L, betaC1L*S1L, 0, 0,
             betaC1H*S1H, betaC1H*S1H, 0, 0,
             0, 0, betaC2L*S2L, betaC2L*S2L,
             0, 0, betaC2H*S2H, betaC2H*S2H),
           nrow = 4, byrow = TRUE
         )
         
         BD <- matrix(
           c(betaD1L*S1L, betaD1L*S1L, 0, 0,
             betaD1H*S1H, betaD1H*S1H, 0, 0,
             0, 0, betaD2L*S2L, betaD2L*S2L,
             0, 0, betaD2H*S2H, betaD2H*S2H),
           nrow = 4, byrow = TRUE
         )
         
         TC1L <- d1+xi1L+phi1L+(delta+d2)*N2/N1
         TC1H <- d1+xi1H+phi1H+(delta+d2)*N2/N1+rhoC
         TC2L <- d2+xi2L+phi2L+delta+rC
         TC2H <- d2+xi2H+phi2H+delta
         
         TD1L <- d1+p1L*eps1L+(delta+d2)*N2/N1
         TD1H <- d1+p1H*eps1H+(delta+d2)*N2/N1+rhoD
         TD2L <- d2+p2L*eps2L+delta+rD
         TD2H <- d2+p2H*eps2H+delta
         
         TC <- matrix(
           c(TC1L, -rhoC, -delta, 0,
             0, TC1H, 0, -delta,
             -(delta+d2)*N2/N1, 0, TC2L, 0,
             0, -(delta+d2)*N2/N1, -rC, TC2H),
           nrow = 4, byrow = TRUE
         )
         
         TD <- matrix(
           c(TD1L, -rhoD, -delta, 0,
             0, TD1H, 0, -delta,
             -(delta+d2)*N2/N1, 0, TD2L, 0,
             0, -(delta+d2)*N2/N1, -rD, TD2H),
           nrow = 4, byrow = TRUE
         )
         
         PHI <- diag(c(phi1L, phi1H, phi2L, phi2H), 4, 4)
         FVinv <- (BC + BD %*% solve(TD) %*% PHI) %*% solve(TC)
         spectral_radius(FVinv)
       })
}

get_S0 <- function(parameters) {
  with(as.list(c(parameters)),
       {
         d <- (delta+d2)*N2/N1
         TS1L <- d1+theta1L+d+alpha1L
         TR1H <- d1+alpha1H+d+rhoR
         TS1H <- d1+theta1H+d+rhoS
         TS2L <- d2+theta2L+delta+alpha2L+rS
         TR2H <- d2+alpha2H+delta+rR
         TS2H <- d2+theta2H+delta
         S <- matrix(
           c(-TS1L, -alpha1L, rhoS-alpha1L, delta, 0, 0,
             0, -TR1H, theta1H, 0, delta, 0,
             0, alpha1H, -TS1H, 0, 0, delta,
             d, 0, 0, -TS2L, -alpha2L, -alpha2L,
             0, d, 0, -rR, -TR2H, theta2H-rR,
             0, 0, d, rS, alpha2H, -TS2H),
           nrow = 6, byrow = TRUE
         )
         b <- matrix(
           c(-alpha1L*N1, 0, 0, -alpha2L*N2, -rR*N2, 0),
           nrow = 6, byrow = TRUE
         )
         S0 <- solve(S) %*% b
         S0 <- S0[c(1,3,4,6)]
         names(S0) <- c("S1L", "S1H", "S2L", "S2H")
         S0
       })
}

parameters <- c(N1 = 100000, N2 = 2000, d1 = 1/(78.5*365), d2 = 0.0068,
                
                alpha1L = 0, theta1L = 0.033, xi1L = 0.0165, phi1L = 0.2, 
                p1L = 0.8, eps1L = 0.1, betaC1L = 0.007/50, betaD1L = 0.007/50,
                
                alpha1H = 0, theta1H = 0.033, xi1H = 0.0165, phi1H = 0.2, 
                p1H = 0.8, eps1H = 0.1, betaC1H = 0.007/50, betaD1H = 0.007/50,
                
                alpha2L = 0.4, theta2L = 0.033, xi2L = 0.0165, phi2L = 0.2, 
                p2L = 0.8, eps2L = 0.1, betaC2L = 0.007, betaD2L = 0.007,
                
                alpha2H = 0, theta2H = 0.033, xi2H = 0.0165, phi2H = 0.2, 
                p2H = 0.8, eps2H = 0.1, betaC2H = 0.007, betaD2H = 0.007,
                
                rR = 0.25, rS = 0.25, rC = 0.25, rD = 0.25,
                rhoR = 0.001, rhoS = 0.001, rhoC = 0.001, rhoD = 0.001,
                delta = 0.135)

get_new_param_vec <- function(x, y, p) {
  y[p] <- as.numeric(x)
  y
}

vary_par <- "alpha2L"

x <- tibble(a2L = seq(0.01, 0.1, by = 0.01)) %>% 
  mutate(
    parms = a2L %>% map(get_new_param_vec, parameters, vary_par),
    S0 = parms %>% map(get_S0),
    R0 = map2(parms, S0, get_R0)
  )

A <- matrix(
  c(unname(x$S0[[1]]), 
    unname(x$S0[[2]]), 
    unname(x$S0[[3]]), 
    unname(x$S0[[4]])
  ), byrow = TRUE, nrow = 4)

b <- x$R0

