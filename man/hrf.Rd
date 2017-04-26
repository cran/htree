\name{hrf}
\alias{hrf}
\title{Random forest for longitudinal data}
\description{Fits a random forest of historical regression trees to longitudinal data.}
\usage{
hrf(x,
    time,
    id,
    yindx,
    ntrees = 100,
    mtry=NULL,
    se=FALSE,
    B=100,
    R=10,
    nsplit=NULL,
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
\item{mtry}{Number of predictors sampled at each split}
\item{se}{If \code{TRUE} then bootstrap standard errors are computed. Total number of trees for fit for bootstrapping is \code{B*R}. }
\item{B}{Only used if \code{se==TRUE}, number of bootstrap samples, defaults to \code{100}.}
\item{R}{Only used if \code{se==TRUE}, forest size for each bootstrap sample, defaults to \code{R=10}.}
\item{nsplit}{Number of splits in regression trees. Defaults to \code{NULL}, in which case tree size determined by requiring terminal nodes to contain atleast 5 observations.}
\item{nsamp}{Number of sampled \code{delta} values, see below}
\item{rsplit}{ If \code{TRUE}, then randomized historical splitting is done, defaults to \code{FALSE}.}
\item{historical}{If \code{TRUE} then historical splitting is done, else only standard (ie concurrent predictor) splitting.}
\item{keep_data}{If \code{TRUE} training data is returned in fit object.}
\item{vh}{Optional vector of indexes giving the historical predictors. Indexes correspond column numbers of \code{x}.}
\item{vc}{Optional vector of indexes giving concurrent predictors.}
\item{order.cat}{If \code{TRUE} then categorical predictors are converted into ordinal predictors according to response mean in categories, else these predictors are converted according to order of appearance in data (default).  }
}

\details{

The \code{hrf} function fits a random forest model to longitudinal data. In particular, it considers data of the form: 
\eqn{(y_{ij},t_{ij},x_{ij})} for \eqn{i=1,..,n} and \eqn{j=1,..,n_i}. Here \eqn{y_{ij}} is the response for the \eqn{i}-th subject 
at the \eqn{j}-th observation time \eqn{t_{ij}}. The predictors at time \eqn{t_{ij}} are \eqn{x_{ij}}. The number of observations can vary across 
subjects, and the sampling in time can be irregular (vary in intensity both within and across subjects).


The novelty of \code{hrf} is that it can construct a predictive model for the response \eqn{y_{ij}} that makes use of \eqn{(t_{ij},x_{ij})} (the observations concurrent with \eqn{y_{ij}}) but also all preceeding observations of the \eqn{i}-th subject upto time \eqn{t_{ij}}. This is 
accomplished using a modification of the standard tree growing algorithm. The modification is to allow nodes to be split based on 
all previous observations of the corresponding subject observation in the node. This is done by creating (invertible) summaries of the observation 
history. In particular, for the  \eqn{k}-th predictor the summary is 
\deqn{n_{ij}(\delta,\tau,k)=\sum_{h: t_{ij}-\delta\leq t_{ih}<t_{ij}} I(x_{ihk}< \tau)}. 
The knowledge of \eqn{n_{ij}(\delta,\tau,k)} for all \eqn{\delta} and \eqn{\tau} is equivalent to knowing all values of predictor \eqn{k} upto time \eqn{t_{ij}}.


A node in a regression tree can be split using the \eqn{n_{ij}()} summaries. In particular, a potential historical split can be detemined as 
\deqn{\mbox{argmin}\sum_{(ij)\in Node} (y_{ij}-\mu_L I(n_{ij}(\delta,\tau,k)<c)-\mu_R I(n_{ij}(\delta,\tau,k)\geq c))^2}

where the minimization is across \eqn{\delta, \tau, \mu, c, k}. For a fixed value of \eqn{delta}, the minimization is solved using a single pass through 
the vector of historic values of predictor \eqn{k}. The \code{hrf} function samples \code{nsamp} values of \code{delta}. The best split among all possible historical splits as well as concurrent splits (ie those 
determined by splitting on measurements made concurrently with the observations of the node) is chosen as the partition. Proceeding recursively, the 
historical regression tree is built. The random forest model is fit using historical regression trees rather than standard regression trees. 


Setting \code{se=TRUE} performs standard error estimation. The number of bootstrap samples (sampling subjects with replacement) is determined by 
\code{B}. For each bootstrap sample a random forest with \code{R} trees is built, which defaults to \code{R=10}. The bias induced by 
using   smaller bootstrap ensemble sizes is corrected for as described in Sexton and Laake (2009). Using \code{se=TRUE} will influence summaries from 
the fitted model, such as providing approximate confidence intervals for partial dependence plots (when running \code{partdep_hrf}), and give 
standard errors for predictions (if desired) when \code{predict_hrf} is used.  




 
}

\value{ Returns a \code{list} with elements: \code{trees} giving the forest, \code{error} the OOB mse, \code{pred} the OOB predictions, \code{boot} 
giving bootstrapped models, as well as arguments of the function call.}  

\references{
L. Breiman (2001). \dQuote{Random Forests,} \emph{Machine Learning} 45(1):5-32.

Sexton and Laake (2009) \dQuote{Standard errors for bagged and random forest estimators,} \emph{Computational Statistics and Data Analysis}.
}

\author{Joe Sexton \email{joesexton0@gmail.com}}

\seealso{ \code{\link{quantile_hrf}}, \code{\link{predict_hrf}}, \code{\link{partdep_hrf}},
          \code{\link{varimp_hrf}}. }

\examples{

\dontrun{
data(mscm) 
mscm=na.omit(mscm)


ff=hrf(x=mscm[,-1],id=mscm$id,time=mscm$day,yindx=3,ntrees=100,nsamp=5,vh=c(2,3),vc=c(1,4:14))
plot(1:length(ff$error),ff$error,type="l",col="blue")


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
ff=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100,mtry=4,nsamp=5)

# plot oob-error 
plot(1:ff$ntrees,ff$error,type="l",col="blue")


# fit forest, also do noisy-bootstrap for standard error estimates
ff=hrf(x=xx,time=time,id=id,yindx=2,ntrees=100,mtry=4,nsamp=5,se=TRUE)

# plot partial dependence of response on its past 
pd=partdep_hrf(ff,xindx=2,ngrid=25,subsample=.1)

# variable importance z-scores (barplot)
barplot(varimp_hrf(ff,nperm=20),main="Importance z-scores")

}
}
\keyword{nonparametric}
\keyword{tree}
\keyword{longitudinal}


