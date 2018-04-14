\name{predict_hrf}
\alias{predict_htb}
\alias{predict_hrf}
\title{Prediction}
\description{
  Prediction functions for \code{hrf} and \code{htb}.
}
\usage{

predict_htb(object,x=NULL,time=NULL,id=NULL,
	ntrees=NULL,type="response",se=FALSE)
predict_hrf(object,x=NULL,time=NULL,id=NULL,
	all.trees=FALSE,se=FALSE)
}
\arguments{
  \item{object}{an object of class \code{htree}.}
  \item{x}{a data frame or matrix containing new data.  If \code{NULL} then training data in \code{object} is used. }
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

# ------------------------------- #
# Simulated data example          #
# ------------------------------- #

p=5;sigma_e=.5;sigma_a=.5;v=rep(1,p);n=500;pnoise=2
random_intercept=rnorm(n,sd=sigma_a)
random_intercept=as.numeric(matrix(random_intercept,nrow=p,ncol=n,byrow=TRUE))
x=random_intercept+rnorm(n*p,sd=sigma_e)
id=sort(rep(1:n,p))
time<-rep(1:p,n)
znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise)
xx=cbind(time,x,znoise)

# fit historical random forest
hb=hrf(x=xx,time=time,id=id,yindx=2,se=TRUE,B=50)

# get predictions with standard errors
pred=predict_hrf(hb,se=TRUE)


# ------------------------------------------------------------------ #
# Comparison of SE-estimates with actual standard errors for 'hrf'
# ------------------------------------------------------------------ #  


## -- evaluation points
n=200
datp=data.frame(y=rep(0,n),w=seq(-2,2,length=n),z=rep(0,n))

## -- estimate model on 50 simulated data sets
pred=NULL
pred_se=NULL
nsim=50
for(k in 1:nsim){
	
	if(is.element(k,seq(1,nsim,by=10)))
		cat(paste("simulation: ",k," of ",nsim," \n",sep=""))
	# -- simulation model -- #
	dat=data.frame(y=(4*datp$w+rnorm(n)),x=datp$w,z=rnorm(n))
	# ---------------------- #
	h=hrf(x=dat,yindx="y",ntrees=500,B=50,R=10,se=T)
	mm=predict_hrf(object=h,x=datp,se=TRUE)
	pred=cbind(pred,mm[,1])
	pred_se=cbind(pred_se,mm[,2])
}


# --- Actual Standard errors at datp 
pred_se_true=apply(pred,1,sd)

# --- Mean of estimated standard errors
pred_se_est=apply(pred_se,1,mean)
pred_se_lower=apply(pred_se,1,quantile,prob=.1)
pred_se_upper=apply(pred_se,1,quantile,prob=.9) 


# -- Plot estimated SE and true SE (+smooth)
z=c(pred_se_true,pred_se_est,pred_se_lower,pred_se_upper)
ylim=c(min(z),max(z))
plot(datp$w,pred_se_est,ylim=ylim,col="blue",xlab="w",
	ylab="Standard error",type="l",main=" SE-true (red) SE-est (blue)")
points(datp$w,pred_se_lower,col="blue",type="l",lty=2)
points(datp$w,pred_se_upper,col="blue",type="l",lty=2)

points(datp$w,pred_se_true,col="red",type="l")

 



# ------------------------------------------------------------------ #
# Comparison of SE-estimates with actual standard errors for 'htb'
# ------------------------------------------------------------------ #  


## -- evaluation points
n=200
datp=data.frame(y=rep(0,n),w=seq(-2,2,length=n),z=rep(0,n))

## -- estimate model on 50 simulated data sets
pred=NULL
pred_se=NULL
nsim=50
for(k in 1:nsim){
	
	if(is.element(k,seq(1,nsim,by=10)))
		cat(paste("simulation: ",k," of ",nsim," \n",sep=""))
	# -- simulation model -- #
	dat=data.frame(y=(4*datp$w+rnorm(n)),x=datp$w,z=rnorm(n))
	# ---------------------- #
	h=htb(x=dat,yindx="y",ntrees=200,cv.fold=10)
	mm=predict_htb(object=h,x=datp,se=TRUE)
	pred=cbind(pred,mm[,1])
	pred_se=cbind(pred_se,mm[,2])
}


# --- Actual Standard errors at datp 
pred_se_true=apply(pred,1,sd)

# --- Mean of estimated standard errors
pred_se_est=apply(pred_se,1,mean)
pred_se_lower=apply(pred_se,1,quantile,prob=.1)
pred_se_upper=apply(pred_se,1,quantile,prob=.9) 


# -- Plot estimated SE and true SE (+smooth)
z=c(pred_se_true,pred_se_est,pred_se_lower,pred_se_upper)
ylim=c(min(z),max(z))
plot(datp$w,pred_se_est,ylim=ylim,col="blue",xlab="w",
	ylab="Standard error",type="l",main=" SE-true (red) SE-est (blue)")
points(datp$w,pred_se_lower,col="blue",type="l",lty=2)
points(datp$w,pred_se_upper,col="blue",type="l",lty=2)

points(datp$w,pred_se_true,col="red",type="l")




}
}
\keyword{classif}% at least one, from doc/KEYWORDS
\keyword{regression}
