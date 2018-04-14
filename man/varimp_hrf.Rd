\name{varimp_hrf}
\alias{varimp_htb}
\alias{varimp_hrf}
\title{Variable importance}
\description{Z-score variable importance for \code{hrf} and \code{htb}}
\usage{
varimp_hrf(object,nperm=20)
varimp_htb(object,nperm=20,ntrees)

}
\arguments{
\item{object}{Return list from \code{hrf} or \code{htb}}
\item{nperm}{Number of permutations.}
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
\value{A \code{data.frame} with  \code{Relative change}: relative change in OOB error due to variable permutation;  \code{Mean change}: mean change in OOB error due to variable permutation;  \code{SE}: standard error of \code{Mean change}; \code{Z-value}:  \code{Mean change}/\code{SE};\code{P-value}: Approximate p-value of \code{Z-value}.   }

\references{
L. Breiman (2001). \dQuote{Random Forests,} \emph{Machine Learning} 45(1):5-32.
}

\seealso{ \code{\link{hrf}}, \code{\link{htb}}}

\examples{
\dontrun{


# --------------------------------------------------------------------------------------------- ##
# Boston Housing data 
#	Comparison of Z-score variable importance with coefficient Z-scores from linear model
# --------------------------------------------------------------------------------------------- ##

# Boston Housing data 
library(mlbench)
data(BostonHousing)
dat=as.data.frame(na.omit(BostonHousing))
dat$chas=as.numeric(dat$chas)

# -- random forest 
h=hrf(x=dat,yindx="medv",ntrees=500)


# -- tree boosting
hb=htb(x=dat,yindx="medv",ntrees=500,cv.fold=10)


# -- Comparison of variable importance Z-scores and Z-scores from linear model 
vi=varimp_hrf(h)
vb=varimp_htb(hb)
dvi=data.frame(var=rownames(vi),Z_hrf=vi$Z)
dvb=data.frame(var=rownames(vb),Z_htb=vb$Z)

dlm=summary(lm(medv~.,dat))$coeffi
dlm=data.frame(var=rownames(dlm),Z_lm=round(abs(dlm[,3]),3))
dlm=merge(dlm[-1,],dvi,by="var",all.x=TRUE)

# -- Z-scores of hrf and lm for predictor variables 
merge(dlm,dvb,by="var",all.x=TRUE)



}

}

