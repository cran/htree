\name{htb}
\alias{htb}
\title{Tree boosting for longitudinal data}
\description{Fits a boosted ensemble of historical regression trees to longitudinal data.}
\usage{


htb(x,
    time,
    id,
    yindx,
    ntrees = 100,
    nsplit=2,
    lambda=.1,
    family="gaussian",
    cv.fold=0,
    cv.rep=NULL,
    nsamp=5,
    rsplit=FALSE,
    historical=TRUE,
    keep_data=TRUE,
    vh=NULL,
    vc=NULL,
    order.cat=FALSE)


}
\arguments{
\item{x}{a data frame containing response and predictors}
\item{time}{vector of observation times}
\item{id}{subject identifier}
\item{yindx}{Column number in \code{x} corresponding response variable}
\item{ntrees}{Number of trees in ensemble}
\item{nsplit}{Number of splits in each regression tree.}
\item{lambda}{Shrinkage parameter.}
\item{family}{Either "gaussian" (default) or "bernoulli".}
\item{cv.fold}{Number of cross-validation folds, if \code{cv<=1} no cross-validation is run. }
\item{cv.rep}{Number of times to repeat the cross-validation. If not given set to \code{cv.fold}.}
\item{historical}{If \code{TRUE} then historical splitting is done, else standard splitting.}
\item{nsamp}{Number of sampled \code{delta} values, see below}
\item{keep_data}{If \code{TRUE} training data is returned in fit object.}
\item{rsplit}{ If \code{TRUE}, then randomized historical splitting is done, defaults to \code{FALSE}.}
\item{vh}{Optional vector of indexes giving the historical predictors.}
\item{vc}{Optional vector of indexes giving strictly concurrent effects.}
\item{order.cat}{If \code{TRUE} then categorical predictors are converted into ordinal predictors according to response mean in categories, else these predictors are converted according to order of appearance in data (default).  }
}

\details{

The \code{htb} function estimates a gradient boosted model to longitudinal data. In particular, it considers data of the form: 
\eqn{(y_{ij},t_{ij},x_{ij})} for \eqn{i=1,..,n} and \eqn{j=1,..,n_i}. Here \eqn{y_{ij}} is the response for the \eqn{i}-th subject 
at the \eqn{j}-th observation time \eqn{t_{ij}}. The predictors at time \eqn{t_{ij}} are \eqn{x_{ij}}. The number of observations can vary across 
subjects, and the sampling in time can be irregular (vary in intensity both within and across subjects).


The novelty of \code{htb} is that it can construct a predictive model for the response \eqn{y_{ij}} that makes use of \eqn{(t_{ij},x_{ij})} (the observations concurrent with \eqn{y_{ij}}) but also all preceeding observations of the \eqn{i}-th subject upto time \eqn{t_{ij}}. This is 
accomplished using a modification of the standard tree growing algorithm. The modification is to allow nodes to be split based on 
all previous observations of the corresponding subject observation in the node. This is done by creating (invertible) summaries of the observation 
history. In particular, for the  \eqn{k}-th predictor the summary is 
\deqn{n_{ij}(\delta,\tau,k)=\sum_{h: t_{ij}-\delta\leq t_{ih}<t_{ij}} I(x_{ihk}< \tau)}. 
The knowledge of \eqn{n_{ij}(\delta,\tau,k)} for all \eqn{\delta} and \eqn{\tau} is equivalent to knowing 
all values of predictor \eqn{k} upto time \eqn{t_{ij}}.


A node in a regression tree can be split using the \eqn{n_{ij}()} summaries. In particular, a potential historical split can be detemined as 
\deqn{\mbox{argmin}\sum_{(ij)\in Node} (y_{ij}-\mu_L I(n_{ij}(\delta,\tau,k)<c)-\mu_R I(n_{ij}(\delta,\tau,k)\geq c))^2}

where the minimization is across \eqn{\delta, \tau, \mu, c, k}. For a fixed value of \eqn{delta}, the minimization is solved using a single pass through 
the vector of historic values of predictor \eqn{k}. The \code{htb} function samples \code{nsamp} values of \code{delta}.The best split among all possible historical splits as well as concurrent splits (ie those 
determined by splitting on measurements made concurrently with the observations of the node) is chosen as the partition. Proceeding recursively, the 
historical regression tree is built. The gradient boosting model is fit as described in Friedman (2001), using historical regression trees instead of 
standard regression trees. 

When \code{cv.fold>1} then cross-validation is performed. In subsequent summaries of the model, say the partial dependence plots from \code{partdep_htb}, 
these cross-validation model fits are used to estimate the standard error. This is done using the delete-d jackknife estimator (where deletion refers to 
subjects, instead of rows of the training data). Each cross-validation model fit is performed by removing a random sample of \code{1/cv.fold} of the  subjects. The number of cross-validation runs is determined by \code{cv.rep} which defaults to the value of \code{cv.fold}. 



 
}

\value{ Returns a \code{list} with elements: \code{trees} giving the tree ensemble, \code{cv} (if \code{cv.fold>1}) giving cross-validation model estimates, \code{cv_error} cross-validated error, as well as arguments of the function call.}

\references{

J.H. Friedman (2001). \dQuote{Greedy Function Approximation: A Gradient Boosting Machine,} \emph{Annals of Statistics} 29(5):1189-1232.

J.H. Friedman (2002). \dQuote{Stochastic Gradient Boosting,} \emph{Computational Statistics and Data Analysis} 38(4):367-378.
}

\author{Joe Sexton \email{joesexton0@gmail.com}}
\seealso{ \code{\link{predict_htb}}, \code{\link{partdep_htb}},
          \code{\link{varimp_htb}}. }
\examples{

\dontrun{
data(mscm) 
mscm=na.omit(mscm)

# logistic regression
ff=htb(x=mscm[,-1],id=mscm$id,time=mscm$day,yindx=3,ntrees=100,
	lambda=.1,nsplit=2,nsamp=5,vh=c(2,3),vc=c(1,4:14),cv.fold=10,family="bernoulli")

# cross-validated negative log-likelihood
plot(1:ff$ntrees,ff$cv_error,type="l",col="blue",ylab="",
	xlab="iterations",main="Cross-validated error")



# -- random intercept example --- #

p=5;sigma_e=.5;sigma_a=.5;v=rep(1,p);n=500;pnoise=2
random_intercept=rnorm(n,sd=sigma_a^.5)
random_intercept=as.numeric(matrix(random_intercept,nrow=p,ncol=n,byrow=TRUE))
x=random_intercept+rnorm(n*p,sd=sigma_e^.5)
id=sort(rep(1:n,p))
time<-rep(1:p,n)
znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise)
xx=cbind(time,x,znoise)

# fit historical random forest
ff=htb(x=xx,time=time,id=id,yindx=2,ntrees=100,nsamp=5,lambda=.1,cv.fold=10)

# plot cv-error 
plot(1:ff$ntrees,ff$cv_error,type="l",col="blue",xlab="iterations",
	ylab="",main="Cross-validated mse")


# plot partial dependence of response on its past 
pd=partdep_htb(ff,xindx=2,ngrid=25,subsample=.1)

# variable importance z-scores (barplot)
barplot(varimp_htb(ff,nperm=20),main="Importance z-scores")
}

}
\keyword{nonparametric}
\keyword{tree}
\keyword{longitudinal}


