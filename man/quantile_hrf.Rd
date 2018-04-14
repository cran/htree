\name{quantile_hrf}
\alias{quantile_hrf}
\title{Quantiles for historical random forest}
\description{
  Quantiles for historical random forest
}
\usage{
quantile_hrf(object, x=NULL,time=NULL,id=NULL,
	ntrees=NULL,prob=seq(.1,.9,length=10))
}
\arguments{
  \item{object}{The return list from running \code{hrf}.}
  \item{x}{a data frame or matrix containing new data.  If \code{NULL} then training data in \code{object} is used. }
   \item{time}{observation times.}
   \item{id}{subject identifiers}
   \item{ntrees}{number of trees to use in prediction.}
  \item{prob}{Probabilities at which quantiles are to be computed}

}
\details{
Meinhausen (2006) showed how a random forest model could be used to estimate conditional quantiles. 
\code{quantile_hrf} implements this technique for a random forest with historical regression trees. 
}

\value{
Returns a matrix of quantiles, with columns corresponding values of \code{prob}.
 }
\references{
 Meinhausen (2006) 'Quantile regression forests' Journal of Machine Learning Research. 
}

\seealso{\code{\link{hrf}}}

\examples{

\dontrun{




# -------------------- ##
# Boston Housing data 
# -------------------- ##
library(mlbench)
library(randomForest)

data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))


## omitting arguments time/id assumes rows are iid 
h=hrf(x=dat,yindx="medv",ntrees=500)

## get deciles on training data 
quant=quantile_hrf(h,prob=seq(.1,.9,length=10))


# ------------------------------------------------ ##
# Simulated data from model with random intercept 
# ------------------------------------------------ ##

p=5;sigma_e=.5;sigma_a=.5;v=rep(1,p);n=500;pnoise=2
random_intercept=rnorm(n,sd=sigma_a)
random_intercept=as.numeric(matrix(random_intercept,nrow=p,ncol=n,byrow=TRUE))
x=random_intercept+rnorm(n*p,sd=sigma_e)
id=sort(rep(1:n,p))
time<-rep(1:p,n)
znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise)
xx=cbind(time,x,znoise)

# fit historical random forest
hb=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100,mtry=4,nsamp=5)

# get deciles on training data 
quant=quantile_hrf(hb,prob=seq(.1,.9,length=10))
}
}
\keyword{quantile}% at least one, from doc/KEYWORDS
\keyword{regression}
\keyword{longitudinal}

