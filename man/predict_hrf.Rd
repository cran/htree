\name{predict_hrf}
\alias{predict_htb}
\alias{predict_hrf}
\title{Prediction}
\description{
  Prediction functions for \code{hrf} and \code{htb}.
}
\usage{

predict_htb(object,x=NULL,yindx=NULL,time=NULL,id=NULL,
	ntrees=NULL,type="response",se=FALSE)
predict_hrf(object,x=NULL,yindx=NULL,time=NULL,id=NULL,
	all.trees=FALSE,se=FALSE)
}
\arguments{
  \item{object}{an object of class \code{htree}.}
  \item{x}{a data frame or matrix containing new data.  If \code{NULL} then training data in \code{object} is used. }
  \item{yindx}{column number corresponding response.}
   \item{time}{observation times.}
   \item{id}{subject identifiers}
   \item{ntrees}{number of trees to use in prediction.}
   \item{type}{If \code{response} then prediction on same scale as response.}
   \item{all.trees}{If \code{TRUE} then prediction for each tree is returned. Only for \code{hrf}.}
  \item{se}{If \code{TRUE} then standard errors of predictions are returned.}

}

\value{
Returns predictions.
 }


\seealso{\code{\link{hrf}},\code{\link{htb}}}

\examples{

\dontrun{
# -- random intercept example --- #

p=5;sigma_e=.5;sigma_a=.5;v=rep(1,p);n=500;pnoise=2
random_intercept=rnorm(n,sd=sigma_a)
random_intercept=as.numeric(matrix(random_intercept,nrow=p,ncol=n,byrow=TRUE))
x=random_intercept+rnorm(n*p,sd=sigma_e)
id=sort(rep(1:n,p))
time<-rep(1:p,n)
znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise)
xx=cbind(time,x,znoise)

# fit historical random forest
hb=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100,mtry=4,nsamp=5,se=TRUE,B=50)

# get predictions with standard errors
pred=predict_hrf(hb,se=TRUE)
}
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{regression}
