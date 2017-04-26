\name{partdep_hrf}
\alias{partdep_hrf}
\alias{partdep_htb}
\title{Partial dependence}
\description{
  Partial dependence for historical tree ensembles.
}
\usage{
partdep_htb(object,xindx,xlim=NULL,ngrid=25,subsample=.1,
	ntrees=NULL,plot.it=TRUE,cat.plot=FALSE)
partdep_hrf(object,xindx,xlim=NULL,ngrid=25,subsample=.1,
	plot.it=TRUE,cat.plot=FALSE)

}
\arguments{
  \item{object}{an object of class \code{htree}.}
  \item{xindx}{Column index of training data, \code{object$x}, for which partial dependence is computed.}
  \item{xlim}{Optional range for partial dependence.}
  \item{ngrid}{Number of values in grid for partial dependence.}
   \item{subsample}{Fraction of training data sampled for calculation.}
   \item{plot.it}{If \code{TRUE} then a plot is produced.}
   \item{ntrees}{Number of trees to in predictions (only for \code{htb}).}
   \item{cat.plot}{If \code{TRUE} then a plot useful for categorical values is shown, defaults to \code{FALSE}.}

}


\value{
Returns \code{x} and \code{y}, the grid and partial dependence values. If \code{htb} run includes cross-validation then 
approximate standard errors 
are also output (based on leave-out-m (subjects) jackknifing). If \code{hrf} estimation is run with \code{se=TRUE} then approximate bootstrap 
standard errors are returned (resampling subjects with replacement). 
 }
\references{
 J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

}

\seealso{\code{\link{hrf}}, \code{\link{htb}}}

\examples{

\dontrun{
# -- random intercept example --- #

p=5;sigma_e=.5;sigma_a=.5;v=rep(1,p);n=500;pnoise=2
random_intercept=rnorm(n,sd=sigma_a)
random_intercept=as.numeric(matrix(random_intercept,nrow=p,ncol=n,byrow=TRUE))
x=random_intercept+rnorm(n*p,sd=sigma_e)
z_cat=sample(c(0,.3),length(x),replace=TRUE)
x=x+z_cat
id=sort(rep(1:n,p))
time<-rep(1:p,n)
znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise)
xx=cbind(time,x,as.numeric(z_cat>0),znoise)

# fit historical random forest
hb=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100,mtry=4,nsamp=5,se=TRUE,B=50)

# partial dependence of second predictor (the historical values of response) 
pd=partdep_hrf(hb,xindx=2,ngrid=25,subsample=.1)

# partial dependence of categorical predictor 
pd=partdep_hrf(hb,xindx=3,ngrid=25,subsample=.1,cat.plot=TRUE)
}

}
\keyword{regression}
