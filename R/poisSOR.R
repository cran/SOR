
pois.SOR <- function(y.formula,
                     w1.formula,
                     w2.formula,
                     y0,
                     hfunc, 
                     id,
                     support,
                     pi1.pi0.ratio,  
                     data = parent.frame(),
                     init.beta=NULL,
                     est.var = T,
                     CORSTR="independence",
                     weights=weights){
  
  # Process the arguments
  DAT.ods <- data 
  ny <- length(support)  
  DAT.ods <- data
  id <- DAT.ods$id
  
  nobs <- length(id)
  nsubj <- length(unique(id))
  obspersubj <- as.numeric(summary(split(id, id))[,1])
  maxtime <- max(obspersubj)
  pi1.pi0.ratio <- DAT.ods$pi.ratio
  phi <- 1
  
  aux <- fitz(DAT.ods, w1.formula, w2.formula, y.formula, hfunc, weights=weights)
  
  mod.z <- aux$mod.z
  W1 <- aux$W1
  W2 <- aux$W2
  gamma.hat <- mod.z$coefficients
  
  terms.y <- terms(y.formula, data=DAT.ods)
  y.matrix <- model.matrix(y.formula, data=DAT.ods)
  
  Y <- model.frame(y.formula, data=DAT.ods)[,attr(terms.y, "response")]
  Z <- mod.z$y
  X <- model.matrix(y.formula,data=DAT.ods)
  
  # Calculate linear predictor for \lambda_S
  # Argument y is a vector, in almost all cases the same as support
  
  
  lambda.p.y0 <- lambda.p.y(y0, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  lpy <- lambda.p.y(support, W1, W2, hfunc, obspersubj, pi1.pi0.ratio, gamma.hat, DAT.ods)
  
  rho.scaled.y0   <- as.vector(1-lambda.p.y0+pi1.pi0.ratio*lambda.p.y0)
  rho.y <- rho.scaled.y(lpy, pi1.pi0.ratio)
  
  ## Calculate F_i (y) for y=1 and y=0
  Fy0    <- Fy(lambda.p.y0, pi1.pi0.ratio)
  F.y <- Fy(lpy, pi1.pi0.ratio)
  
  
  # The odds funcion, includes the poisson stuff as reference distribution
  odds.s.y <- function(y, eta){
    ny <- length(y)
    rhomat <- log( apply(rho.y ,2,"/", rho.scaled.y0))
    ymat <- apply(matrix(rep(y, nobs), nrow=nobs, byrow=T), 2, "-", y0)
    etaymat <- apply(ymat, 2, "*", eta)
    referencedist <- -lfactorial(matrix(rep(y, nobs), nrow=nobs, byrow=T)) + lfactorial(y0)
    #exp(( etaymat/phi + referencedist/phi + rhomat ))
    exp(weights*( etaymat/phi + referencedist/phi + rhomat ))
  }
  
  
  # Inverse link is int(y*odds)/int(odds)
  InvLink <- function(eta){
    oddsout <- odds.s.y(support, eta)
    as.vector((oddsout%*%support)/(oddsout%*%rep(1, ny)))
  }
  
  LinkFun <- log
  
  # Variance function is int((y^2)*odds)/int(odds) - (int(y*odds)/int(odds))^2
  VarFun <- function(eta){
    oddsout <- odds.s.y(support, eta)
    as.vector((oddsout%*%(support^2))/(oddsout%*%rep(1, ny))) - as.vector((oddsout%*%support)/(oddsout%*%rep(1, ny)))^2
  }
  
  # InvLinkDeriv is the same as Variance function
  InvLinkDeriv <- function(eta){
    VarFun(eta)
  }
  
  
  y.minus.lam <- function(phi, lambda.p){
    rhomat <- log(rho.y)
    yfacmat <- apply(matrix(rep(support, nobs), nrow=nobs),2, lfactorial)
    yloglam <- apply(matrix(rep(support, nobs), nrow=nobs), 2, "*", log(lambda.p))
    ymat <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", lambda.p) - yfacmat
    toint <- -ymat*(exp(ymat/phi + rhomat))
    as.vector(apply(toint, 1, sum))
  }
  
  y.minus.lam.2 <- function(phi, lambda.p){
    rhomat <- log(rho.y)
    yfacmat <- apply(matrix(rep(support, nobs), nrow=nobs),2, lfactorial)
    yloglam <- apply(matrix(rep(support, nobs), nrow=nobs), 2, "*", log(lambda.p))
    ymat <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", lambda.p) - yfacmat
    toint <- (ymat^2)*(exp(ymat/phi + rhomat))
    as.vector(apply(toint, 1, sum))
  }
  
  
  normalize <- function(phi, lambda.p){
    rhomat <- log(rho.y)
    yfacmat <- apply(matrix(rep(support, nobs), nrow=nobs),2, lfactorial)
    yloglam <- apply(matrix(rep(support, nobs), nrow=nobs), 2, "*", log(lambda.p))
    ymat <- apply(matrix(rep(support, nobs), nrow=nobs, byrow=T), 2, "-", lambda.p) - yfacmat
    toint <- (exp(ymat/phi + rhomat))
    as.vector(apply(toint, 1, sum))
  }
  
  
  FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)
  
  if(est.var){
    f.sig <- function(phi, lambda.p){
      #sum(-Y*log(lambda.p)+lambda.p + lfactorial(Y)) - sum(y.minus.lam(phi, lambda.p)/normalize(phi, lambda.p))
      sum((-Y*log(lambda.p)+lambda.p + lfactorial(Y))*weights) - sum(weights*y.minus.lam(phi, lambda.p)/normalize(phi, lambda.p))
    }
    
    df.sig <- function(phi, lambda.p){
      #-(1/(phi^2))*(sum(y.minus.lam.2(phi, lambda.p)/normalize(phi, lambda.p)) + sum((y.minus.lam(phi, lambda.p)^2)/(normalize(phi, lambda.p)^2)))
      -(1/(phi^2))*(sum(weights*y.minus.lam.2(phi, lambda.p)/normalize(phi, lambda.p)) + sum(weights*(y.minus.lam(phi, lambda.p)^2)/(normalize(phi, lambda.p)^2)))
    }
    sigstop <- F
    phi.new <- 0
    
    while(!sigstop){
      odds.s.y <- function(y, eta){
        ny <- length(y)
        rhomat <- log( apply(rho.y ,2,"/", rho.scaled.y0))
        ymat <- apply(matrix(rep(y, nobs), nrow=nobs, byrow=T), 2, "-", y0)
        etaymat <- apply(ymat, 2, "*", eta)
        referencedist <- -lfactorial(matrix(rep(y, nobs), nrow=nobs, byrow=T)) + lfactorial(y0)
        #exp( etaymat/phi + referencedist/phi + rhomat )
        exp(weights* (etaymat/phi + referencedist/phi + rhomat ))
      }
      mod.y <- geemR(y.formula, family=FunList, id = id, data=DAT.ods, corstr = CORSTR, Mv = 1, init.beta=init.beta, tol=1e-6, scale.fix=T, init.phi=phi, weights=weights)    
      init.beta <- mod.y$beta
      lambda.p <- exp(X%*%init.beta)
      
      phi.new <- phi - f.sig(phi, lambda.p)/df.sig(phi, lambda.p)
      
      if(abs(phi - phi)/phi < 0.001){
        sigstop <- T
      }
      phi <- phi.new
      
    }
  }
  #print(phi)
  
  
  
  phi <- 1
  # Run through geemR to calculate coefficients
  mod.y <- geemR(y.formula, family=FunList, id = id, data=DAT.ods, corstr = CORSTR, Mv = 1, init.beta = init.beta, scale.fix=T, init.phi=phi, weights=weights)
  
  se.coefs <- getvar(mod.z, mod.y, odds.s.y, VarFun, Y,"d", hfunc, support, F.y, Fy0, y0, W1, W2, obspersubj)#, weights=weights)
  
  return(list("coefs" = mod.y$beta, "se.coefs" = se.coefs))
}