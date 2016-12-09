sor <- function(y.formula,
                w1.formula,
                w2.formula = ~1,
                id,
                waves=NULL,
                family = "binomial",
                y0 = 0,
                hfunc = identity, 
                support = c(0,1),
                pi1.pi0.ratio = 1,  
                data = parent.frame(),
                init.beta=NULL,
                init.sig.2 = 1,
                weights=NULL,
                est.var = TRUE,
                CORSTR="independence"){
  call <- match.call()
  
  fam <- as.character(charmatch(family, c("gaussian", "normal", "poisson", "binomial")))
  if(fam=="2") fam <- "1"
  
  if(typeof(data) == "environment"){
    id = id
    waves <- waves
    weights <- weights
    DAT.ods <- data.frame(model.frame(y.formula), model.frame(w1.formula), model.frame(w2.formula))
  }else{
    DAT.ods <- data 
    nn <- dim(DAT.ods)[1]
    if(length(call$id) == 1){
      subj.col <- which(colnames(data) == call$id)
      if(length(subj.col) > 0){
        id <- data[,subj.col]
      }else{
        id <- eval(call$id, envir=parent.frame())
      }
    }else if(is.null(call$id)){
      id <- 1:nn
    }
    if(length(call$weights) == 1){
      weights.col <- which(colnames(data) == call$weights)
      if(length(weights.col) > 0){
        weights <- data[,weights.col]
      }else{
        weights <- eval(call$weights, envir=parent.frame())
      }
    }else if(is.null(call$weights)){
      weights <- rep.int(1,nn)
    }
    if(length(call$waves) == 1){
      waves.col <- which(colnames(data) == call$waves)
      if(length(waves.col) > 0){
        waves <- data[,waves.col]
      }else{
        waves <- eval(call$waves, envir=parent.frame())
      }
    }else if(is.null(call$waves)){
      waves <- NULL
    }
  }
  ny <- length(support)  
  DAT.ods$id       <- id
  DAT.ods$waves <- waves
  DAT.ods$weights <- weights
  DAT.ods$offset.z <- log(pi1.pi0.ratio) 
  DAT.ods$pi.ratio <- pi1.pi0.ratio
  
  if(!is.numeric(DAT.ods$waves) & !is.null(DAT.ods$waves)) stop("waves must be either an integer vector or NULL")
  

  # W is diagonal matrix of weights, sqrtW = sqrt(W)
  # included is diagonal matrix with 1 if weight > 0, 0 otherwise
  # includedvec is logical vector with T if weight > 0, F otherwise
  # Note that we need to assign weight 0 to rows with NAs
  # in order to preserve the correlation structure
  na.inds <- NULL
  
  if(any(is.na(DAT.ods))){
    na.inds <- which(is.na(DAT.ods), arr.ind=T)
  }
  
  #SORT THE DAT.odsA ACCORDING TO WAVES
  if(!is.null(waves)){
    DAT.ods <- DAT.ods[order(id, waves),]
  }else{
    DAT.ods <- DAT.ods[order(id),]
  }
  
  
  # Figure out the correlation structure
  cor.vec <- c("independence", "ar1", "exchangeable", "m-dependent", "unstructured", "fixed", "userdefined")
  cor.match <- charmatch(CORSTR, cor.vec)
  
  if(is.na(cor.match)){stop("Unsupported correlation structure")}
  
  if(!is.null(DAT.ods$waves)){
    wavespl <- split(DAT.ods$waves, DAT.ods$id)
    idspl <- split(DAT.ods$id, DAT.ods$id)
    
    maxwave <- rep(0, length(wavespl))
    incomp <- rep(0, length(wavespl))
    
    for(i in 1:length(wavespl)){
      maxwave[i] <- max(wavespl[[i]]) - min(wavespl[[i]]) + 1
      if(maxwave[i] != length(wavespl[[i]])){
        incomp[i] <- 1
      }
    }
    
    #If there are gaps and correlation isn't exchangeable or independent
    #then we'll add some dummy rows
    if(!is.element(cor.match, c(1,3)) & (sum(incomp) > 0) ){
      DAT.ods <- dummyrows(y.formula, DAT.ods, incomp, maxwave, wavespl, idspl)
      id <- DAT.ods$id
      waves <- DAT.ods$waves
      weights <- DAT.ods$weights
    }
  }
  
  if(!is.null(na.inds)){
    weights[unique(na.inds[,1])] <- 0
    for(i in unique(na.inds)[,2]){
      if(is.factor(DAT.ods[,i])){
        DAT.ods[na.inds[,1], i] <- levels(DAT.ods[,i])[1]
      }else{
        DAT.ods[na.inds[,1], i] <- median(DAT.ods[,i], na.rm=T)
      }
    }
  }
  
  
  includedvec <- weights>0
  
  
  inclsplit <- split(includedvec, id)
  
  dropid <- NULL
  allobs <- T
  if(any(!includedvec)){
    allobs <- F
    for(i in 1:length(unique(id))){
      if(all(!inclsplit[[i]])){
        dropid <- c(dropid, i)
      }
    }
  }
  
  
  if(length(dropid)>0){
    dropind <- which(is.element(id, dropid))
    DAT.ods <- DAT.ods[-dropind,]
    includedvec <- includedvec[-dropind]
    weights <- weights[-dropind]
    
    id <- id[-dropind]
  }
  nn <- dim(DAT.ods)[1]
  K <- length(unique(id))
  
  
  
  
  yframe <- model.frame(y.formula, DAT.ods)
  resp <- model.response(yframe)
  w1frame <- model.frame(w1.formula, DAT.ods)
  ref <- model.response(w1frame)
  w2frame <- model.frame(w2.formula, DAT.ods)
  
  if(!all( is.element(ref, c(0,1)))){stop("Response in w1.formula must be binary")}
  if(!all(ref==model.response(w2frame)) && !is.null(model.response(w2frame))){stop("There should be no response for w2.formula")}
  
  if(max(resp) > max(support) | min(resp) < min(support)) warning("Some values of response are outside support")
  
  switch(fam,
          #### NORMAL ####
          "1"= {
            results <- norm.SOR(y.formula, w1.formula, w2.formula, y0, hfunc, id, support, pi1.pi0.ratio, DAT.ods, init.beta, init.sig.2, est.var, CORSTR, weights=DAT.ods$weights)                
            },
          #### POISSON ####
          "3" = {
            if(!all( is.wholenumber(resp))){stop("Response in y.formula must be positive integers if family is Poisson")}
            if(init.sig.2 != 1){warning("Initial variance specification ignored")}
            results <- pois.SOR(y.formula, w1.formula, w2.formula, y0, hfunc, id, support, pi1.pi0.ratio, DAT.ods, init.beta,  est.var=est.var, CORSTR=CORSTR, weights=DAT.ods$weights)
           },
          #### BINOMIAL ####
           "4" = {
             if(!all( is.element(resp, c(0,1)))){stop("Response in y.formula must be binary if family is binomial")}
             if(y0 != 0){warning("Setting y0 to 0")}
             results <- binom.SOR(y.formula, w1.formula, w2.formula, DAT.ods, init.beta,  CORSTR)
           }
         )
  class(results) <- "sor"
  results$coefnames <- colnames(model.matrix(y.formula, DAT.ods))
  results$call <- call
  results$family <- c("normal","normal", "Poisson", "binomial")[as.numeric(fam)]
  return(results)
}

print.sor <- function(x, ...){
  coefdf <- data.frame(x$coefs)
  rownames(coefdf) <- x$coefnames
  colnames(coefdf) <- ""
  print(x$call)
  cat("\n", "Coefficients:", "\n")
  print(t(coefdf))
  cat("\n Reference Distribution: ", x$family, "\n")
}

summary.sor <- function(object, ...)  {
  Coefs <- matrix(0,nrow=length(object$coefs),ncol=4)
  Coefs[,1] <- c(object$coefs)
  Coefs[,2] <- object$se.coefs
  Coefs[,3] <- Coefs[,1]/Coefs[,2]
  Coefs[,4] <- round(2*pnorm(abs(Coefs[,3]), lower.tail=F), digits=8)
  colnames(Coefs) <- c("Estimates","Std. Err.", "wald", "p")
  
  summ <- list(coefs = Coefs[,1], se.coefs = Coefs[,2],  wald.test = Coefs[,3], p = Coefs[,4],
                coefnames = object$coefnames, family=object$family, call=object$call)
  
  class(summ) <- 'summary.sor'
  return(summ)
}

print.summary.sor <- function(x, ...){
  Coefs <- matrix(0,nrow=length(x$coefnames),ncol=4)
  rownames(Coefs) <- c(x$coefnames)
  colnames(Coefs) <- c("Estimates","Std. Err.", "wald", "p")
  Coefs[,1] <- x$coefs
  Coefs[,2] <- x$se.coefs
  Coefs[,3] <- x$wald.test
  Coefs[,4] <- x$p
  
  print( x$call)
  cat("\n")
  print(Coefs)
  cat("\n Reference Distribution: ", x$family, "\n")
}

