library(tidyverse)
library(deSolve)

# model specification
parameters <- c(
  aR = 0.75, aS = 0.09, aCn = 0.06, aCp = 0.09, aD = 0.01, aV = 0.0, aCv = 0.0,
  kR = 0.33, kS = 0.15, kCn = 0.15, kCp = 0.15, kD = 0.068, kV = 0.15, kCv = 0.15,
  alpha = 0.5, theta = 0.033, phi = 0.06, p = 0.8, eps = 0.10,
  beta = 0.000001, f = 0.6,
  c0 = 3, c1 = 0.1, c2 = 15, c3 = 0.5, c4 = 0.5, c5 = 15, M = 0.3
)
state_0 <- c(R = 112, S = 14, Cn = 9, Cp = 13, D = 2, V = 0, Cv = 0)
adjoint_T <- c(lambdaR = 0, lambdaS = 0, lambdaCn = as.numeric(parameters["c0"]), lambdaCp = 0, lambdaD = 0, lambdaV = 0, lambdaCv = 0)
control_0 <- c(v = 0)
t_int <- c(0,30)

state <- function(t, y, p, u) { # u should be return of approxfun
  with(as.list(c(y, p)),
       {
         betaC <- beta
         betaD <- beta
         u <- u(t)
         deltaN <- kR*R + kS*S + kCn*Cn + kCp*Cp + kD*D + kV*V + kCv*Cv
         lambda <- betaC*(Cn+Cp+Cv) + betaD*D
         dR <- aR*deltaN + theta*S - (alpha+kR)*R
         dS <- aS*deltaN + alpha*R + p*eps*D - (theta+lambda+u+kS)*S
         dCn <- aCn*deltaN + (1-f)*lambda*S - (phi+u+kCn)*Cn
         dCp <- aCp*deltaN + f*lambda*S - (u+kCp)*Cp
         dD <- aD*deltaN + phi*Cn - (p*eps+kD)*D
         dV <- aV*deltaN + u*S - (lambda+kV)*V
         dCv <- aCv*deltaN + lambda*V + u*(Cn+Cp) - kCv*Cv
         list(c(dR, dS, dCn, dCp, dD, dV, dCv))
       })
}

adjoint <- function(t, y, p, u, S, Cn, Cp, D, V, Cv) { # u, etc. should be return of approxfun
  with(as.list(c(y, p)),
       {
         betaC <- beta
         betaD <- beta
         u <- u(t)
         S <- S(t)
         Cn <- Cn(t)
         Cp <- Cp(t)
         D <- D(t)
         V <- V(t)
         Cv <- Cv(t)
         lambda <- betaC*(Cn+Cp+Cv) + betaD*D
         dlambdaR <- - lambdaR*(aR*kR-alpha-kR) - lambdaS*(aS*kR+alpha) - kR*(aCn*lambdaCn+aCp*lambdaCp+aD*lambdaD+aV*lambdaV+aCv*lambdaCv)
         dlambdaS <- - c1*u - c5*lambda - lambdaR*(aR*kS+theta) - lambdaS*(aS*kS-(theta+lambda+u+kS)) - lambdaCn*(aCn*kS+(1-f)*lambda) - lambdaCp*(aCp*kS+f*lambda) - lambdaD*aD*kS - lambdaV*(aV*kS+u) - lambdaCv*aCv*kS 
         dlambdaCn <- - c1*u - c4 - c5*betaC*(S+V) - lambdaR*aR*kCn - lambdaS*(aS*kCn-betaC*S) - lambdaCn*(aCn*kCn+(1-f)*betaC*S-(phi+kCn+u)) - lambdaCp*(aCp*kCn+f*betaC*S) - lambdaD*(aD*kCn+phi) - lambdaV*(aV*kCn-betaC*V) - lambdaCv*(aCv*kCn+u+betaC*V)
         dlambdaCp <- - c1*u - c5*betaC*(S+V) - lambdaR*aR*kCp - lambdaS*(aS*kCp-betaC*S) - lambdaCn*(aCn*kCp+(1-f)*betaC*S) - lambdaCp*(aCp*kCp+f*betaC*S-kCp-u) - lambdaD*(aD*kCp) - lambdaV*(aV*kCp-betaC*V) - lambdaCv*(aCv*kCp+u+betaC*V)
         dlambdaD <- - c3 - c5*betaD*(S+V) - lambdaR*aR*kD - lambdaS*(aS*kD+p*eps-betaD*S) - lambdaCn*(aCn*kD+(1-f)*betaD*S) - lambdaCp*(aCp*kD+f*betaD*S) - lambdaD*(aD*kD-p*eps-kD) - lambdaV*(aV*kD-betaD*V) - lambdaCv*(aCv*kD+betaD*V)
         dlambdaV <- -c5*lambda - lambdaR*aR*kV - lambdaS*aS*kV - lambdaCn*aCn*kV - lambdaCp*aCp*kV - lambdaD*aD*kV - lambdaV*(aV*kV-kV-lambda) - lambdaCv*(aCv*kV+lambda)
         dlambdaCv <- - c5*betaC*(S+V) - lambdaR*(aR*kCv) - lambdaS*(aS*kCv-betaC*S) - lambdaCn*(aCn*kCv+(1-f)*betaC*S) - lambdaCp*(aCp*kCv+f*betaC*S) - lambdaD*(aD*kCv) - lambdaV*(aV*kCv-betaC*V) - lambdaCv*(aCv*kCv-kCv+betaC*V)
         list(c(dlambdaR, dlambdaS, dlambdaCn, dlambdaCp, dlambdaD, dlambdaV, dlambdaCv))
       })
}

fbsweep <- function(x_0, lambda_T, u_0, t_int, p, x, lambda, delta = 0.001, t_length = 1000) {
  times <- seq(t_int[1], t_int[2], length = t_length + 1)
  df <- tibble(
    system = c("state", "adjoint", "control"),
    data = list(
      as_tibble(t(x_0)),
      as_tibble(t(lambda_T)),
      as_tibble(t(u_0))
    )
  ) %>% 
    mutate(
      data = pluck(., "data") %>% 
        map(map_df,~c(.x,rep(0, t_length))),
      data_mat = map(data, as.matrix)
    )
  test <- -1
  while (test < 0) {
    df <- df %>%
      mutate(old_data_mat = data_mat)
    v_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 3)), rule = 2)
    df$data_mat[[1]] <- ode(
      times = times, y = x_0, func = x, parms = p, u = v_fun
    ) %>% .[,2:8]
    S_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 1)[,"S"]), rule = 2)
    Cn_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 1)[,"Cn"]), rule = 2)
    Cp_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 1)[,"Cp"]), rule = 2)
    D_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 1)[,"D"]), rule = 2)
    V_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 1)[,"V"]), rule = 2)
    Cv_fun <- approxfun(times, as.vector(pluck(df, "data_mat", 1)[,"Cv"]), rule = 2)
    df$data_mat[[2]] <- ode(
      times = rev(times), y = lambda_T, func = lambda, parms = p, u = v_fun, S = S_fun, Cn = Cn_fun, Cp = Cp_fun, D = D_fun, V = V_fun, Cv = Cv_fun
    ) %>% .[nrow(.):1,2:8]
    temp <- (df$data_mat[[1]][,"S"]*df$data_mat[[2]][,"lambdaS"] + df$data_mat[[1]][,"Cn"]*df$data_mat[[2]][,"lambdaCn"] + df$data_mat[[1]][,"Cp"]*df$data_mat[[2]][,"lambdaCp"] - df$data_mat[[1]][,"S"]*df$data_mat[[2]][,"lambdaV"] - (df$data_mat[[1]][,"Cn"]+df$data_mat[[1]][,"Cp"])*df$data_mat[[2]][,"lambdaCv"] - p["c1"]*(df$data_mat[[1]][,"S"]+df$data_mat[[1]][,"Cn"]+df$data_mat[[1]][,"Cp"])) / (2*p["c2"])
    df$data_mat[[3]][,"v"] <- (df$data_mat[[3]][,"v"] + pmin(p["M"], pmax(0, temp))) / 2
    df <- df %>%
      mutate(
        convergence_vals = map2_dbl(data_mat, old_data_mat, ~ delta*norm(.x) - norm(.y - .x))
      )
    test <- min(df$convergence_vals)
  }
  df %>%
    mutate(data = map(data_mat, as_tibble)) %>%
    filter(system != "adjoint") %>%
    select(data) %>%
    map_df(bind_cols) %>%
    mutate(
      t = times,
      sumC = Cn + Cp + Cv,
      sumV = V + Cv
    ) %>%
    select(t, R, S, Cn, sumC, D, sumV, v) %>% 
    gather(func, val, -t, factor_key = TRUE)
}

get_new_par_val_vec <- function(x) {
  for (i in seq_along(x)) {
    names(x)[i] <- as.character(x[i])
  }
  x
}
get_new_param_vec <- function(x, y, p) {
  y[p] <- as.numeric(x)
  y
}

vary_par <- "beta"
par_vals <- c(0.000001, 0.0007, 0.007)
par_vals <- get_new_par_val_vec(par_vals)

out_df <- par_vals %>% 
  map(get_new_param_vec, parameters, vary_par) %>%
  map(~fbsweep(state_0,adjoint_T,control_0,t_int,.,state, adjoint)) %>% 
  bind_rows(.id = "par_vals") %>% 
  mutate(par_vals = as.factor(par_vals))

out_df %>% 
  ggplot(aes(t, val, col = par_vals)) + 
  geom_line() +
  facet_wrap(~func, nrow = 4, scales = "free_y") +
  scale_color_discrete(name=vary_par)
