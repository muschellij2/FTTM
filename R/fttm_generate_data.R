
beta_s <- function(s) {
  1 * cos (pi * s)
}

#' Generate Data for FTTM Option 1
#'
#' @param n number of participants/sample size
#' @param m Number of functional time points
#' @param mu_c censoring times mean
#' @param seed seed for sampling/data generation
#'
#' @returns A list of observed and true time observations, including censoring
#' times
#' @export
fttm_generate_data = function(
    n = 500,
    m = 101,
    mu_c = 20,   #
    seed = 123) {
  withr::with_seed(seed = seed, {
    Z1 <- rbinom(n, 1, 0.5)
    Z2 <- rnorm(n, 0, 1)
  })
  X <- cbind(rep(1, n), Z1, Z2)
  ###############################
  S <- seq(0, 1, length = m)
  rankX <- 10
  Phi <- cbind(1 / sqrt(m), stats::poly(S, degree = rankX - 1))
  lambda <- rankX:1
  Xi <- sapply(lambda, function(l)
    scale(rnorm(n, sd = 2 * sqrt(l)), scale = FALSE))
  Xs <- Xi %*% t(Phi)

  beta <- c(log(0.2), -0.5, 0.4)

  betaS <- beta_s(S)
  lfterm2 <- as.numeric((1 / m) * Xs %*% betaS) #xs=zs

  ratesurv <- function(i) {
    exp(sum(beta * X[i, ]) + lfterm2[i])
  }
  ratelat <- unlist(lapply(1:n, ratesurv)) #
  T_tr <- c() ##survival times true
  C <- c() #censoring time
  T_obs <- c() #observed survival time
  for (i in 1:n)
  {
    temp <- stats::rexp(1, rate = ratelat[i])
    T_tr[i] <- temp #true survival time
    C[i] <- stats::rexp(1, rate = (1 / mu_c))
    T_obs[i] = min(T_tr[i], C[i])
  }
  delta <- as.numeric(T_tr < C)

  list(
    n = n,
    data_scalar = X,
    data_functional = Xs,
    functional_domain = S,
    T_obs = T_obs,
    T_true = T_tr,
    C = C,
    delta = delta
  )
}
