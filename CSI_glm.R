glmCSI <- function(object, constraints=NULL, verbose=FALSE)  {
  
  if(!is.null(constraints)) {
    stop("no constraints specified")
  }
  if(!class(object)[1]=="glm") {
    stop("it only works for class is glm()")
  }
  
  fit.glm <- object
    cnames <- variable.names(fit.glm)
  
  X <- model.matrix(fit.glm)
  y <- model.response(fit.glm$model)
  n <- nrow(X)
  p <- ncol(X)
  
  family <- fit.glm$family[1]
  link <- fit.glm$family[2]
  
  stopifnot(family %in% c("gaussian", "binomial", "poisson"))
  
  if(family=="gaussian" && link=="identity") { #FIX constraint Amat x2 is now x1
    logl <- logl.gauss.identity
    logl.gr <- logl.gauss.identity.grad
    startv <- c(1,rep(0L, length(coef(fit.glm))))
  }
  else if(family=="binomial" && link=="logit") {
    logl <- logl.binomial.logit
    logl.gr <- logl.binomial.logit.grad
    startv <- rep(0L, length(coef(fit.glm)))
  }
  else if(family=="poisson" && link=="log") {
    logl <- logl.poisson.log
    logl.gr <- logl.poisson.log.grad
    startv <- rep(0L, length(coef(fit.glm)))
  }
  else {
    stop("The following family/link combinations are supported (for now): 
          \"family=Gaussian and link=identity\"  
          \"family=Binomial and link=logit\"  
          \"family=Poisson  and link=log\"")
  }
  
  
  lavpartable <- lav_partable_from_lm(fit.glm, est=TRUE, label=TRUE)
  constraints0 <- gsub("[<>]","==", constraints) 
  
  CON1 <- lavaan:::lav_constraints_parse(constraints=constraints, partable=lavpartable, theta=)
  CON0 <- lavaan:::lav_constraints_parse(constraints=constraints0, partable=lavpartable)
  heq0 <- CON0$ceq.function
  hin1 <- CON1$cin.function
  
  
  mle0 <- nlminb.constr(start=startv, objective=logl, X=X, y=y,
                        gradient=logl.gr, hessian=NULL, 
                        ceq=heq0, ceq.jac=NULL,
                        cin=NULL, cin.jac=NULL, control.outer=list(verbose=verbose))
    
  heq0 <- if(nrow(CON1$ceq.JAC) > 0) { CON1$ceq.function } else { NULL }
  
  mle1 <- nlminb.constr(start=startv, objective=logl, X=X, y=y,
                        gradient=logl.gr, hessian=FALSE, 
                        ceq=heq0, ceq.jac=NULL,
                        cin=hin1, cin.jac=NULL, control.outer=list(verbose=verbose))
  
  Amat <- rbind(CON1$ceq.JAC, CON1$cin.JAC)
  #Amat <- rbind(c(0,1,0))                                                       #FIXME!!!
  #Wald test statistic
  if(family=="gaussian") {
    ##QLR test
    s2 <- summary(fit.glm)$dispersion
    QLR <- -2/s2*( logl(mle1$par, y, X) - logl(mle0$par, y, X) )
    
    ##Wald test statistic
    #residuals under the ic. constrained model
    res1 <- y - X%*%mle1$par[-1]
    #information matrix under the ic. model
    I <- solve(1/(n-p) * as.numeric(t(res1)%*%res1) * solve(t(X)%*%X))
    Wald <- 1/s2*(t(Amat%*%mle1$par[-1]) %*% solve(Amat%*%solve(I)%*%t(Amat)) 
                  %*% (Amat%*%mle1$par[-1]))
    
    ##Score test statistic
    #residuals under the eq. constrained model
    res0 <- y - X%*%mle0$par[-1]
    I <- solve(1/(n-p) * as.numeric(t(res0)%*%res0) * solve(t(X)%*%X))
    s20 <- sum(res0^2)/(n-p)
    d0 <- t(1/s20 * t(X)%*%(y-X%*%mle0$par[-1]))
    d1 <- t(1/s2 * t(X)%*%(y-X%*%mle1$par[-1]))
    G0 <- d0
    G1 <- d1
    Score <- 1/s2*(G1-G0)%*%solve(I)%*%t(G1-G0)      
  } 
  else if(family=="binomial") {
      ##QLR test statistic
      QLR <- -2*( logl(mle1$par, y, X) - logl(mle0$par, y, X) )
      
      ##Wald test statistic
      pr = 1/(1+exp(-X%*%mle1$par))
      V = array(0,dim=c(dim(X)[1],dim(X)[1]))
      diag(V) = pr*(1-pr)
      I = t(X)%*%V%*%X    
      Wald <- t(Amat%*%mle1$par) %*% solve(Amat%*%solve(I)%*%t(Amat)) %*% (Amat%*%mle1$par)
    
      ##Score test statistic
      # The variance of a Bernouilli distribution is given by p(1-p)
      pr = 1/(1+exp(-X%*%mle0$par))
      V = array(0,dim=c(dim(X)[1],dim(X)[1]))
      diag(V) = pr*(1-pr)
      I = t(X)%*%V%*%X    
      G0 <- (y%*%X) - t(t(X)%*%(exp(X%*%mle0$par)/(1 + exp(X%*%mle0$par))))
      G1 <- (y%*%X) - t(t(X)%*%(exp(X%*%mle1$par)/(1 + exp(X%*%mle1$par))))
      Score <- 1*(G1-G0)%*%solve(I)%*%t(G1-G0)
  }  
  else if(family=="poisson") {
    ##QLR test statistic
    QLR <- -2*( logl(mle1$par, y, X) - logl(mle0$par, y, X) )
    
    ##Wald test statistic
    #information matrix under the ic. model
    I <- matrix(0L, nrow=dim(X)[2], ncol=dim(X)[2])
    for(i in 1:nrow(X)){
      bZ <- t(mle1$par)%*%(X[i,])
      Z <- X[i,]
      ZZ <- Z%*%t(Z)
      I <- I + ZZ*as.numeric(exp(bZ))
    }
    Wald <- t(Amat%*%mle1$par) %*% solve(Amat%*%solve(I)%*%t(Amat)) %*% (Amat%*%mle1$par)
    
    #information matrix under the eq. model
    I <- matrix(0L, nrow=dim(X)[2], ncol=dim(X)[2])
    for(i in 1:nrow(X)){
      bZ <- t(mle0$par)%*%(X[i,])
      Z <- X[i,]
      ZZ <- Z%*%t(Z)
      I <- I + ZZ*as.numeric(exp(bZ))
    }
    G0 <- t(t(X)%*%(y-exp(X%*%mle0$par)))
    G1 <- t(t(X)%*%(y-exp(X%*%mle1$par)))
    Score <- 1*(G1-G0)%*%solve(I)%*%t(G1-G0)
  }
  
  cat(" ...QLR:", formatC(QLR, digits=4), " ...Wald:", formatC(Wald, digits=4), 
      " ...Score:", formatC(Score, digits=4), "\n")

  meq <- nrow(CON1$ceq.JAC)
  df <- 0:(nrow(Amat) - meq) 
  wt_bar <- wt(cov=Amat%*%solve(I)%*%t(Amat), meq=meq)
  p_QLR   <- sapply(QLR,  function(x) 1-pchibar(x, df1a=df, wt=rev(wt_bar)))  
  p_Wald  <- sapply(Wald, function(x) 1-pchibar(x, df1a=df, wt=rev(wt_bar)))
  p_Score <- sapply(Score,  function(x) 1-pchibar(x, df1a=df, wt=rev(wt_bar)))
  
  out <- data.frame(QLR=QLR, Wald=Wald, Score=Score, pvalue_QLR=p_QLR, pvalue_Wald=p_Wald, 
                    pvalue_Score=p_Score)
  
  out
}  

