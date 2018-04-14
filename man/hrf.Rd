\name{hrf}
\alias{hrf}
\title{Random forest for longitudinal data}
\description{Fits a random forest of historical regression trees to longitudinal data.}
\usage{



hrf(x,
    time=NULL,
    id=NULL,
    yindx,
    ntrees = 100,
    method="freqw",
    mtry=NULL,
    se=FALSE,
    B=100,
    R=10,
    nsplit=NULL,
    nsamp=5,
    historical=TRUE,
    keep_data=TRUE,
    vh=NULL,
    vc=NULL,
    delta=NULL,
    control=list(nmax=10,nodesize=1))
}
\arguments{
\item{x}{a data frame containing response and predictors}
\item{time}{vector of observation times}
\item{id}{subject identifier, if \code{NULL} observations are assumed independent}
\item{yindx}{Name of response variable, alt. its column number in \code{x}}
\item{ntrees}{Number of trees in ensemble}
\item{method}{Historical summary method, can be \code{freq}, \code{frac}, \code{mean0}, \code{freqw}, \code{fracw} and \code{mean0w}}
\item{mtry}{Number of predictors sampled at each split}
\item{se}{If \code{TRUE} then bootstrap standard errors are computed. Total number of trees for fit for bootstrapping is \code{B*R}. }
\item{B}{Only used if \code{se==TRUE}, number of bootstrap samples, defaults to \code{100}.}
\item{R}{Only used if \code{se==TRUE}, forest size for each bootstrap sample, defaults to \code{R=10}.}
\item{nsplit}{Number of splits in regression trees. Defaults to \code{NULL}, in which case tree size determined by requiring terminal nodes to contain atleast 5 observations.}
\item{nsamp}{Number of sampled \code{delta} values, see below}
\item{historical}{If \code{TRUE} then historical splitting is done, else only standard (ie concurrent predictor) splitting.}
\item{keep_data}{If \code{TRUE} training data is returned in fit object.}
\item{vh}{Optional vector of indexes giving the historical predictors. Indexes correspond column numbers of \code{x}.}
\item{vc}{Optional vector of indexes giving concurrent predictors.}
\item{delta}{A vector of history lags to be used (see below), defaults to \code{NULL} in which case all possible observed lags are available for splitting.}
\item{control}{A list of control parameters.}

}

\details{

The \code{hrf} function fits a random forest model to longitudinal data. Data is assumed to be of form: 
\eqn{z_{ij}=(y_{ij},t_{ij},x_{ij})} for \eqn{i=1,..,n} and \eqn{j=1,..,n_i}, with \eqn{y_{ij}} being the response for the \eqn{i}-th subject 
at the \eqn{j}-th observation time \eqn{t_{ij}}. The vector of predictors at time \eqn{t_{ij}} are \eqn{x_{ij}}. The number of observations can vary across 
subjects, and the sampling in time can be irregular.


\code{hrf} estimates a model for the response \eqn{y_{ij}} using both \eqn{(t_{ij},x_{ij})} (the observations concurrent with \eqn{y_{ij}}) and all preceeding observations of the \eqn{i}-th subject up to (but not including) time \eqn{t_{ij}}. The model is fit using \code{historical} regression (alt. classification) trees. Here a predictor is one 
of two types, either \code{concurrent} or \code{historical}. The concurrent predictors for \eqn{y_{ij}} are the elements of the vector (\eqn{(t_{ij},x_{ij})}), while a 
historic predictor is the set of all preceeding values (preceeding time \eqn{t_{ij}}) of a given element of \eqn{(y_{ij},t_{ij},x_{ij})}, for subject \eqn{i}. In 
a historic regression tree, node splitting on a \code{concurrent} predictor follows the approach in standard regression (classification) trees. For \code{historical} predictors 
the splitting is modified since, associated with each observed response \eqn{y_{ij}}, the number (and observation times) of observations of a \code{historical} predictor will vary according to 
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



Setting \code{se=TRUE} performs standard error estimation. The number of bootstrap samples (sampling subjects with replacement) is determined by 
\code{B}. For each bootstrap sample a random forest with \code{R} trees is built, which defaults to \code{R=10}. The bias induced by 
using   smaller bootstrap ensemble sizes is corrected for in the estimate. Using \code{se=TRUE} will influence summaries from 
the fitted model, such as providing approximate confidence intervals for partial dependence plots (when running \code{partdep_hrf}), and give 
standard errors for predictions when \code{predict_hrf} is used.  






 
}

\value{ Returns a \code{list} with elements: \code{trees} giving the forest, \code{error} the OOB mse, \code{pred} the OOB predictions, \code{boot} 
giving bootstrapped models, as well as arguments of the function call.}  

\references{
L. Breiman (2001). \dQuote{Random Forests,} \emph{Machine Learning} 45(1):5-32.

Zhang and Singer (2010) \dQuote{Recursive Partioning and Applications} \emph{Springer}.

Sexton and Laake (2009) \dQuote{Standard errors for bagged and random forest estimators,} \emph{Computational Statistics and Data Analysis}.
}

\author{Joe Sexton \email{joesexton0@gmail.com}}

\seealso{ \code{\link{quantile_hrf}}, \code{\link{predict_hrf}}, \code{\link{partdep_hrf}},
          \code{\link{varimp_hrf}}. }

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

# -- first two weeks
mscm=mscm[mscm$day<=14,]

# -- set concurrent and historical predictors 
historical_predictors=match(c("stress","illness"),names(mscm))
concurrent_predictors=which(names(mscm)!="stress")

# -- fit historical random forest
#	(NOTE: response is 0/1 so a regression tree is equivalent 
##	 to a classification tree with Gini-index splitting) 
ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",
	vh=historical_predictors,vc=concurrent_predictors)

# -- variable importance table
vi=varimp_hrf(ff)
vi


# -- fit historical random forest with 'se=TRUE'
ff=hrf(x=mscm,id=mscm$id,time=mscm$day,yindx="illness",
	vh=historical_predictors,vc=concurrent_predictors,se=T,B=50)

# -- partial dependence with top 4 predictors
par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_hrf(ff,xindx=rownames(vi)[k])

par(mfrow=c(1,1))




# -------------------------------------------- ##
# Simulated data from random intercept model
# -------------------------------------------- ##

## -- simulated data -------------------------------------------------- #
p=5;sigma_e=.5;sigma_a=.5;v=rep(1,p);n=500;pnoise=2
random_intercept=rnorm(n,sd=sigma_a^.5)
random_intercept=as.numeric(matrix(random_intercept,nrow=p,ncol=n,byrow=TRUE))
x=random_intercept+rnorm(n*p,sd=sigma_e^.5)
id=sort(rep(1:n,p))
time<-rep(1:p,n)
znoise=matrix(rnorm(n*p*pnoise),ncol=pnoise)
xx=as.data.frame(cbind(time,x,znoise))
xx$rfact=as.factor(sample(letters,nrow(xx),replace=TRUE))
## --------------------------------------------------------------------- #

# fit historical random forest
ff=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100)

# plot oob-error 
plot(1:ff$ntrees,ff$error,type="l",col="blue")


# fit forest with standard error selection
ff=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100,se=TRUE)

# plot partial dependence of response on its past 
pd=partdep_hrf(ff,xindx=2)

# variable importance table
varimp_hrf(ff)




# -------------------------------------------- ##
# Examples with non-longitudinal data
# -------------------------------------------- ##



# -------------------- ##
# Boston Housing data 
# -------------------- ##
library(mlbench)
library(randomForest)

data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))


## omitting arguments time/id assumes rows are iid 
h=hrf(x=dat,yindx="medv",ntrees=500)

## randomForest comparison
r=randomForest(medv~.,data=dat)

## plot oob-error for both
plot(1:length(r$mse),r$mse,type="l",ylim=c(min(r$mse,h$error),max(r$mse,h$error)),
	main="BostonHousing",xlab="forest size",ylab="out-of-bag mse")
points(1:length(h$error),h$error,type="l",col="blue")

## -- variable importance table
vi=varimp_hrf(h)
vi

## -- partial dependence plots with approximate 95% C.I 
h=hrf(x=dat,yindx="medv",ntrees=500,se=TRUE)

par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_hrf(h,xindx=rownames(vi)[k])

par(mfrow=c(1,1))

# ------------------------------------- ##
# Ionosphere data (classification)  
# ---------------------------------- ##
data(Ionosphere)
dat=as.data.frame(na.omit(Ionosphere))


h=hrf(x=dat,yindx="Class",ntrees=500)
r=randomForest(Class~.,data=dat)

plot(1:length(r$err[,1]),r$err[,1],type="l",ylim=c(min(r$err[,1],h$error),max(r$err[,1],h$error)),
	ylab="mis-classification rate",xlab="forest size",main="Out-of-bag mis-classification")
points(1:length(h$error),h$error,type="l",col="blue")

# -- variable importance table
vi=varimp_hrf(h)


par(mfrow=c(2,2))
for(k in 1:4)
	pd=partdep_hrf(h,xindx=rownames(vi)[k])


par(mfrow=c(1,1))




}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{longitudinal}


