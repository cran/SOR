\name{sor}
\alias{sor}
\alias{SOR}
\alias{SOR-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sequentially Offsetted Regression
}
\description{
Fits model for data which was sampled based on a variable associated with the outcome.  This function works for binary, count, and continuous responses.
}
\usage{
	    sor(y.formula,
              	w1.formula,
                w2.formula = ~1,
                id,
                waves = NULL,
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
                CORSTR="independence")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.formula}{
Regression formula for response
}
  \item{w1.formula}{
Formula for Z, not interacted with hfunc(Y).  Of form Z~terms
}
  \item{w2.formula}{
Formula for Z, interacted with hfunc(Y).  Of form ~terms
}
  \item{id}{
a vector identifying the clusters. By default, data are assumed to be sorted such that observations in a cluster are in consecutive rows and higher numbered rows in a cluster are assumed to be later.  If NULL, then each observation is assigned its own cluster.
}
  \item{waves}{
an integer vector identifying components of a cluster. For example, this could be a time ordering. If integers are skipped within a cluster, then dummy rows with weight 0 are added in an attempt to preserve the correlation structure (except if \code{corstr = "exchangeable"} or \code{"independent"}). This can be skipped by setting \code{nodummy=TRUE}.
}
  \item{family}{
Character string representing reference distribution for the response.  Can be one of "normal", "poisson", or "binomial".
}
  \item{y0}{
Representative value of response.  Ignored if \code{family="binomial"}.
}
  \item{hfunc}{
Function h, used with Y.  Set to identity if \code{family="binomial"}.
}

  \item{support}{
Values on which to evaluate the integrals.  The lowest value should be less than the minimum response and the highest should be higher than the maximum response.  If response is binary, support should be c(0,1).  If response is count data, support should be an integer vector, for instance 0:50.  If response is continuous, support should be a vector of points on which to integrate.
}
  \item{pi1.pi0.ratio}{
The referral ratio
}
  \item{data}{
Data frame or environment with all the data
}
  \item{init.beta}{
Initial values for parameters in \code{y.formula}.  Convergence may depend heavily on the initial values used.  If \code{family="binomial"}, the default is recommended.
}
  \item{init.sig.2}{
Initial value for sigma^2.  Only for \code{family="normal"}.
}
  \item{weights}{
  A vector of weights for each observation.  If an observation has weight 0, it is excluded from the calculations of any parameters.  Observations with a \code{NA} anywhere (even in variables not included in the model) will be assigned a weight of 0.  This should normally be used to preserve the correlation structure.
}
  \item{est.var}{
Logical.  Should the variance be estimated.  Only for \code{family="normal"}.
}
  \item{CORSTR}{
Correlation structure
}
}
\value{
Returns a list with values from the fit.
%% ...
}
\references{
This package relies heavily on code from geeM:

McDaniel, L. S., Henderson, N. C., & Rathouz, P. J. (2013). Fast pure R implementation of GEE: application of the matrix package. The R journal, 5(1), 181.
}
\author{
Lee S. McDaniel, Jonathan S. Schildcrout
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
generatedata <- function(beta,alpha,X,ntime,nsubj, betat, betat1)  {
  
  mean.vec <- exp(crossprod(t(X), beta))
  y <- matrix(0,nrow=nsubj, ncol=ntime)
  y[,1] <- rpois(nsubj ,lambda = mean.vec)
  old.mean <- mean.vec
  new.mean <- old.mean*exp(betat + betat1*X[,2])
  for (t in 1:(ntime-1))  {
    lambda.t <- new.mean - alpha*sqrt(old.mean*new.mean)
    theta.t <- alpha*sqrt(new.mean/old.mean)
    I <- rpois(nsubj, lambda = lambda.t)
    W <- rbinom(nsubj, y[,t], theta.t)
    
    y[,t+1] = W + I
    old.mean <- new.mean
    new.mean <- old.mean*exp(betat + betat1*X[,2])
  }
  longform <- c(t(y))
  time <- rep(1:ntime,times=nsubj)
  subject <- rep(c(1:nsubj),each=ntime)
  
  simdata <- data.frame(count = longform, time = time, subject=subject)
  return(simdata)
}
logit <- function(p) log(p)-log(1-p)
expit <- function(x) exp(x)/(1+exp(x))
set.seed(1)

npop <- 10000
beta0 <- -1.4
beta1 <- 0.4
alpha <- 0.9
gam0 <- -3.15
gam1 <- 6.3
nsubj <- 200
ntime <- 8
betat <- -0.1; betat1 <- 0.1
thresh <- 1

x0 <- rep(1, npop)
x1 <- rbinom(npop, 1, 0.5)

Xmat <- cbind(x0, x1)
timevec <- 0:(ntime-1)

testdat <- generatedata(c(beta0, beta1), alpha, Xmat, ntime, npop, betat = betat, betat1 = betat1)
Y <- matrix(testdat$count, nrow=npop, ncol=ntime, byrow=TRUE)
lambdap <- expit(gam0 + gam1*as.numeric(Y[,1]>=thresh))
Z <- rbinom(npop, 1, lambdap)

casesamp <- rep(0, npop)
casesamp[Z==1] <- rbinom(sum(Z), 1, nsubj/(2*sum(Z)))
controlsamp <- rep(0, npop)
controlsamp[Z==0] <- rbinom(sum(1-Z), 1, nsubj/(2*sum(1-Z)))


case <- which(casesamp==1)
control <- which(controlsamp==1)
id <- sort(c(case, control))
nsubj <- length(control) + length(case)
Ysamp <- NULL
lamsamp <- NULL
zsamp <- NULL
x1samp <- NULL
idsamp <- NULL
time <- NULL
obspersubj <- sample(3:ntime, size=nsubj, replace=TRUE)
for(i in 1:nsubj){
  Ysamp <- c(Ysamp, Y[id[i],1:obspersubj[i]])
  zsamp <- c(zsamp, rep(as.numeric(Z[id[i]]), obspersubj[i]))
  x1samp <- c(x1samp, rep(x1[id[i]], obspersubj[i]))
  time <- c(time, 0:(obspersubj[i]-1))
  idsamp <- c(idsamp, rep(i, obspersubj[i]))
}
p1p0 <- sum((1-Z))/sum(Z)

timemax <- pmax(time-2, 0)
y0 <- 1
betas <- c(beta0, beta1, betat, betat1)
init <-  runif(4, betas -0.1, betas + 0.1)

y.formula <- y~x1+time + x1:time
w1 <- z ~ x1+ as.factor(time) +  x1:time + x1:timemax 
w2 <- ~x1 + time + timemax + x1:time + x1:timemax 

DAT.ods <- data.frame("x1"= x1samp, "time" = time,
                      "timemax" = timemax, "z" = zsamp, "y" = Ysamp, "id" = idsamp)

sor(y.formula, w1, w2, id, family="poisson",y0=1,
    support=0:25, pi1.pi0.ratio=p1p0, data=DAT.ods, init.beta=init, CORSTR="ar1")



}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

