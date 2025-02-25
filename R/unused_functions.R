
ChooseAICval <- function(Nrvec,
                         p1,
                         T_obs, delta, Xmat, Xs, S, tau,
                         verbose = TRUE)
{
  N0 <- as.numeric(Nrvec[1])
  N1 <- as.numeric(Nrvec[2])
  r <- as.numeric(Nrvec[3])
  initval <- choose.init(T_obs, delta, Xmat, Xs, S, N1 = N1)#N0 not needed here
  if (verbose) {
    print(initval)
  }
  betainit <- initval$betainit
  thetainit <- initval$thetainit
  gammainit <- sort((rnorm((N0 + 1), 0, 1)))#c(-0.5,0.1,0.3,0.5) #has to be increasing #H(t)=beta0+logt
  #need to give better values by understanding H#
  etainit <- c()
  etainit[1] <- gammainit[1]
  for (k in 2:length(gammainit))
  {
    etainit[k] <- log(gammainit[k] - gammainit[k - 1]) #these cannot be all same.. chose baseline hazard differently
  }
  psi <- c(betainit, etainit, thetainit)
  out <- optim(
    psi,
    loglikfun2,
    T_obs = T_obs,
    delta = delta,
    Xmat = Xmat,
    Xs = Xs,
    S = S,
    tau = tau,
    r = r,
    N0 = N0,
    N1 = N1,
    method = "BFGS",
    control = list(trace = 1, maxit = 500)
  ) #inc maxit
  psiest <- out$par
  AICval <- ICfunc(psiest, N0, N1, r,
                   p1, T_obs, delta, Xmat, Xs, S, tau)

  #BICval
  result <- list(AICval = AICval, psiest = psiest)
  #print(Nrvec)
  return(result)
}


beta_s_estfun <- function(s, theta, N1)
{
  basismat <- NULL
  for (k in 0:N1) {
    basismat <- cbind(basismat, stats::dbeta((s), shape1 = (k + 1), shape2 = N1 - k +
                                               1))
  }
  #basis should add up to 1
  basismat <- (1 / (N1 + 1)) * basismat
  val <- as.numeric(basismat %*% theta)
  val
}

