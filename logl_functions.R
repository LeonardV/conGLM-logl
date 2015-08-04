########
#gaussian log-likelihood distribution
logl.gauss.identity <- function(theta, y, X) {
  if (theta[1] <= 0) return(0.01)
  ll = sum(suppressWarnings(dnorm(y, mean=X%*%theta[-1], sd=sqrt(theta[1]), log=TRUE)))
    return(-ll)
}

logl.gauss.identity.grad <- function(theta, y, X) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  e <- y - X%*%beta
  n <- nrow(X)
  g <- numeric(length(theta))
  g[1] <- (-n/(2*sigma2)) + (t(e)%*%e)/(2*sigma2*sigma2) 
  g[-1] <- (t(X) %*% e)/sigma2                           
    return(-g)
}

##########################################################
#logit log-likelihood function #
logl.binomial.logit <- function(theta, y, X) {
  ll <- sum(y*log(exp(X%*%theta)/(1 + exp(X%*%theta))) + 
                 (1-y)*log(1/(1 + exp(X%*%theta))))
  return(-ll)
}

logl.binomial.logit.grad <- function(theta, y, X) {
  gr <- theta*0
  eXb <- exp(X%*%theta) 
  prob1 <- eXb/(1+eXb)
  for (k in 1:ncol(X)) { 
    gr[k] <- sum(X[,k] * (y - prob1))
  }
  return(-gr)
}

######## Econometrics I, Green ########
#poisson log-likelihood distribution
logl.poisson.log <- function(theta, y, X) {
  Xb <- X%*%theta
  lambda <- exp(Xb)
  ll <- sum(-lambda + y*log(lambda) - log(factorial(y)))
  return(-ll)
}

logl.poisson.log.grad <- function(theta, y, X) {
  Xb <- X%*%theta
  lambda <- exp(Xb)
  #gr[1] <- sum(-n + sum(y) / lambda)
  gr <- t(X)%*%(y-lambda)
  return(-gr)
}

#gamma log-likelihood distribution (link=identiy)
logl.gamma.identiy <- function(theta, y, X, shape) {
  # log likelihood for Gamma(alpha,scale=beta)
  #ll <- sum(dgamma(y, shape=theta[1], scale=theta[2], log=TRUE))
  Xb <- X%*%theta
  #rate = shape / mean:
  rate <- shape / Xb
  # sum of negative log likelihoods:
  ll <- sum(dgamma(y, rate=rate, shape=shape, log=TRUE))
  #cat("ll =", ll, "...theta =", theta, "\n")
    return(-ll)
}


logl.gamma.idenity.grad <- function(theta, y, X) {
  Xb <- X%*%theta
  
  gr <- 
  return(-gr)
}


#gamma log-likelihood distribution (link=log)
logl.gamma.log <- function(theta, y, X, logshape) {
  # log likelihood for Gamma(alpha,scale=beta)
  #ll <- sum(dgamma(y, shape=theta[1], scale=theta[2], log=TRUE))
  Xb <- X%*%theta
  #rate = shape / mean:
  rate <- exp(shape) / exp(Xb)
  ll <- sum(dgamma(y, rate=rate, shape=exp(logshape), log=TRUE))
  return(-ll)
}

logl.gamma.log.grad <- function(theta, y, X) {
  Xb <- X%*%theta
  
  gr <- 
    return(-gr)
}

#gamma log-likelihood distribution (link=inverse)                               #FIXME
logl.gamma.inverse <- function(theta, y, X, shape) {
  # log likelihood for Gamma(alpha,scale=beta)
  #ll <- sum(dgamma(y, shape=theta[1], scale=theta[2], log=TRUE))
  Xb <- X%*%theta
  #rate = shape / mean:
  rate <- shape / Xb
  # sum of negative log likelihoods:
  #ll <- sum(1/dgamma(y, rate=rate, shape=shape, log=TRUE))
  ll <- sum(1/(scale*gamma(shape)) * (scale / (y-Xb))^(shape+1) * exp(-scale/(y-Xb) ))
  return(-ll)
}

logl.gamma.inverse.grad <- function(theta, y, X) {
  Xb <- X%*%theta
  
  gr <- 
    return(-gr)
}


