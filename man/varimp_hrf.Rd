\name{varimp_hrf}
\alias{varimp_htb}
\alias{varimp_hrf}
\title{Variable importance}
\description{Z-score variable importance for \code{hrf} and \code{htb}}
\usage{
varimp_hrf(object,nperm=10)
varimp_htb(object,nperm=10,ntrees)

}
\arguments{
\item{object}{Return list from \code{hrf} or \code{htb}}
\item{nperm}{Vector of observation times}
\item{ntrees}{Number of trees. Only for \code{varimp_htb}.}

}

\details{
To measure the importance of a predictor, \code{varimp_hrf} and \code{varimp_htb} compare the prediction errors of the estimated model with the 
prediction errors obtained after integrating the predictor out of the model. If \eqn{F} denotes the estimated model, the model obtained by 
integrating out predictor \code{k} is \eqn{F_k(x)=\int F(x) dP(x_k)}. Here \eqn{P(x_k)} is the marginal distribution of \eqn{x_k}. In 
practice, the integration is done by averaging over multiple predictions from \eqn{F},  where each has been obtained using a random permutation of the observed values of \eqn{x_k}. The number of 
permutations is determined by \code{nperm}. Letting \eqn{L(y,y_{hat})}) be the loss of predicting \eqn{y} with \eqn{y_{hat}}, one obtains 
the vector \eqn{w_i=L(y_i,F_k(x_i))-L(y_i,F(x_i))} for \eqn{i=1,..,n}. The corresponding z-score is 
\eqn{z=mean(w_i)/se(w_i)}, which is an approximate paired test for the equality of the prediction errors. Larger z-score values indicate that the prediction error increases if \eqn{x_k} is marginalized out, and thus that \eqn{x_k} is useful. On the other hand, large negative values of the z-score indicate that the integrated model is more accurate.  For longitudinal data the \code{w_i} are computed by averaging across all observations from the \code{i}-th subject. For \code{htb} the prediction error is calculated based on the cross-validation model estimates, for \code{hrf} out-of-bag predictions 
are used. 
}
\value{z-scores for the predictors.}

\references{
L. Breiman (2001). \dQuote{Random Forests,} \emph{Machine Learning} 45(1):5-32.
}

\seealso{ \code{\link{hrf}}, \code{\link{htb}}}

\examples{
\dontrun{
data(mscm) 
mscm=na.omit(mscm)

# -- random forest model (predicting illness, with stress and illness as historical predictors)
ff=hrf(x=as.matrix(mscm),id=mscm$id,time=mscm$day,yindx=4,vh=c(3,4),vc=c(1,2,5:14))
vi=varimp_hrf(ff,nperm=20)
vi=vi[vi>0]
vi=vi[order(vi,decreasing=TRUE)]
barplot(vi,main="Importance z-scores") 

}

}

