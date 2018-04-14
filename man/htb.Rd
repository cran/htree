\name{htb}
\alias{htb}
\title{Tree boosting for longitudinal data}
\description{Fits a boosted ensemble of historical regression trees to longitudinal data.}
\usage{


htb(x,
    time=NULL,
    id=NULL,
    yindx,
    ntrees = 100,
    method="freqw",
    nsplit=1,
    lambda=.1,
    family="gaussian",
    cv.fold=0,
    cv.rep=NULL,
    nsamp=5,
    historical=TRUE,
    keep_data=TRUE,
    vh=NULL,
    vc=NULL,
    delta=NULL)


}
\arguments{
\item{x}{a data frame containing response and predictors}
\item{time}{vector of observation times}
\item{id}{subject identifier, if \code{NULL} observations are assumed independent}
\item{yindx}{Column number in \code{x} corresponding response variable}
\item{ntrees}{Number of trees in ensemble}
\item{method}{Historical summary method, can be \code{freq}, \code{frac}, \code{mean} or \code{sum}}
\item{nsplit}{Number of splits in each regression tree.}
\item{lambda}{Shrinkage parameter.}
\item{family}{Either "gaussian" (default) or "bernoulli".}
\item{cv.fold}{Number of cross-validation folds, if \code{cv<=1} no cross-validation is run. }
\item{cv.rep}{Number of times to repeat the cross-validation. If not given set to \code{cv.fold}.}
\item{historical}{If \code{TRUE} then historical splitting is done, else standard splitting.}
\item{nsamp}{Number of sampled \code{delta} values, see below}
\item{keep_data}{If \code{TRUE} training data is returned in fit object.}
\item{vh}{Optional vector of indexes giving the historical predictors.}
\item{vc}{Optional vector of indexes giving strictly concurrent effects.}
\item{delta}{A vector of history lags to be used (see below), defaults to \code{NULL} in which case all possible observed lags are available for splitting.}
}

\details{

The \code{htb} function fits a random forest model to longitudinal data. Data is assumed to be of form: 
\eqn{z_{ij}=(y_{ij},t_{ij},x_{ij})} for \eqn{i=1,..,n} and \eqn{j=1,..,n_i}, with \eqn{y_{ij}} being the response for the \eqn{i}-th subject 
at the \eqn{j}-th observation time \eqn{t_{ij}}. The predictors at time \eqn{t_{ij}} are \eqn{x_{ij}}. The number of observations can vary across 
subjects, and the sampling in time can be irregular.


\code{htb} estimates a model for the response \eqn{y_{ij}} using both \eqn{(t_{ij},x_{ij})} (the observations concurrent with \eqn{y_{ij}}) and all preceeding observations of the \eqn{i}-th subject up to (but not including) time \eqn{t_{ij}}. The model is fit using \code{historical} regression (alt. classification) trees. Here a predictor is one 
of two types, either \code{concurrent} or \code{historical}. The concurrent predictors for \eqn{y_{ij}} are the elements of the vector (\eqn{(t_{ij},x_{ij})}), while a 
historic predictor is the set of all preceeding values (preceeding time \eqn{t_{ij}}) of a given element of \eqn{(y_{ij},t_{ij},x_{ij})}, for subject \eqn{i}. In 
a historic regression tree, node splitting on a \code{concurrent} predictor follows the approach in standard regression (classification) trees. For \code{historical} predictors 
the splitting is modified since, associated with each observed response \eqn{y_{ij}}, the set of observations of a \code{historical} predictor will vary according to 
\eqn{i} and \eqn{j}. For these, the splitting is done by first transforming the preceeding values of a predictor using a \code{summary function}. This summary is 
invertible, in the sense that knowledge of it is equivalent to knowledge of the covariate history. Letting \eqn{\bar{z}_{ijk}} denote the set of historical values of the 
\eqn{k}-th element of \eqn{z_{ij}}, the summary function is denoted \eqn{s(\eta;\bar{z}_{ijk})} where \eqn{\eta} is the argument vector of the summary 
function. Node splitting based on a historical predictor is done by solving \deqn{\mbox{argmin}\sum_{(ij)\in Node} (y_{ij}-\mu_L I(s(\eta;\bar{z}_{ijk})<c)-\mu_R I(s(\eta;\bar{z}_{ijk})\geq c))^2,} where the minimization is over the vector \eqn{(k,\mu,c,\eta)}. Each node of \code{historical} regression tree is split according to the best split among all possible splits based on both concurrent and historical predictors. 




Different \code{summary functions} are available in \code{hrf}, specified by the argument \code{method}, where \code{method="freq"} and corresponds the summary 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} I(z_{ihk}<\eta_2)};  \code{method="frac"}:  
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} I(z_{ihk}<\eta_2)/n_{ij}(\eta)};  \code{method="mean0"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}} z_{ihk}/n_{ij}(\eta)};  \code{method="freqw"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} I(z_{ihk}<\eta_3)};  \code{method="fracw"}: 
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} I(z_{ihk}<\eta_3)/n_{ij}(\eta)};  \code{method="meanw0"}:  
\deqn{s(\eta;\bar{z}_{ijk})=\sum_{h: t_{ij}-\eta_1\leq t_{ih}<t_{ij}-\eta_2} z_{ihk}/n_{ij}(\eta)}. Here \eqn{n_{ij}(\eta)} denotes the number of observations of subject 
\eqn{i} in the 
time window \eqn{[t_{ij}-\eta_1,t_{ij}-\eta_2)}. In the case  \eqn{n_{ij}(\eta)=0}, the summary function is set to zero, i.e \eqn{s(\eta;\bar{z}_{ijk})=0}.   The default is \code{method="freqw"}. The possible values of \eqn{\eta_1} in the summary function can be set by the argument \code{delta}. If not supplied, the set of possible values of 
\eqn{\eta_1} is 
determined by the difference in time between successive observations in the data. When a split is attempted on a historical predictor, a sample of this set is 
taken from which the best split is selected. The size of this set equals that of the \code{nsamp} argument. 




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


# ----------------------------------------------------------------------------------------- ##
# Mother's stress on child illness:
# 	Investigate whether mother's stress is (Granger) causal on child illness 
#	'hrf' model is fit using previous observations of mother's stress to predict 
#	child's illness at given time point, but not mother's stress at that time point
#
#	Predictor variables are classified into "historical" and "concurrent"
#
#	A predictor is "historical" if its prior realizations can be used to predict 
#	the outcome.   
#
#	A predictor is "concurrent" if its realization at the same timepoint as the outcome
#	can be used to predict the outcome at that timepoint
#
#	A predictor can be both "concurrent" and "historical", the status of the predictors 
#	can be set by the 'vh' and 'vc' arguments of 'hrf'. 
#	(if not set these are automatically determined) 
#	 
# ------------------------------------------------------------------------------------------- ##


data(mscm) 
mscm=as.data.frame(na.omit(mscm))

## first 2 weeks
mscm=mscm[mscm$day<=14,] 

# -- set concurrent and historical predictors 
historical_predictors=match(c("illness","stress"),names(mscm))
concurrent_predictors=which(names(mscm)!="stress")


# logistic regression
ff=htb(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",ntrees=100,lambda=.5,method="freqw",
	vh=historical_predictors,vc=concurrent_predictors,cv.fold=10,family="bernoulli")

# cross-validated negative log-likelihood
plot(1:ff$ntrees,ff$cv_error,type="l",col="blue",ylab="",
	xlab="iterations",main="Cross-validated deviance")

vi=varimp_htb(ff)
vi

par(mfrow=c(2,2))
for(k in 1:4)
	partdep_htb(ff,xindx=rownames(vi)[k])

par(mfrow=c(1,1))

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

# variable importance table
varimp_htb(ff,nperm=20)




# ------------------------------------------ #
# Non-longitudinal data examples
# ------------------------------------------ #


# ------------------------------------------ #
# Boston Housing data 
# ------------------------------------------ #

library(htree)
library(mlbench)
data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))


# omitting arguments 'time' and 'id' assumes rows are iid 
h=htb(x=dat,yindx="medv",ntrees=500,cv.fold=10)

# -- plot cross-validated Mean-squared error --- #
plot(1:length(h$cv_error),h$cv_error,type="l",xlab="Boosting iterations",
	ylab="MSE",main="Cross-validated Mean-squared error")

# -- iteration with lowest error
iter_hat=order(h$cv_error)[1]

# -- variable importance table
vi=varimp_htb(h,nperm=20)
vi 

# -- partial dependence of top 4 predictors with +/-2 S.E.
#	(standard errors are based on the 
#	delete-d jack-knife estimator, using the cross-validation runs)
par(mfrow=c(2,2))
for(k in 1:4)
	partdep_htb(h,xindx=rownames(vi)[k])


par(mfrow=c(1,1))

}

}
\keyword{nonparametric}
\keyword{tree}
\keyword{longitudinal}


