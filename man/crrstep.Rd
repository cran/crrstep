\name{crrstep}
\alias{crrstep}
\title{Stepwise regression for competing risks regression}
\description{
Performs forward and backward stepwise regression for the Fine & Gray regression model in competing risks. Procedure uses AIC, BIC and BICcr as selection criteria. BICcr has a penalty of k = log(n*), where n* is the number of Type I events.
}
\usage{
crrstep(formula, scope.min = ~1, etype, ..., subset, 
data, direction = c("backward", "forward"), 
criterion = c("AIC", "BICcr", "BIC"), crr.object = FALSE, 
trace = TRUE, steps = 100)
}
\arguments{
  \item{formula}{
formula object where LHS is failure time and RHS is linear predictors; intercept `1' should always be included.
}
  \item{scope.min}{
formula object denoting final model for backward selection and starting model for forward selection.
}
  \item{etype}{
integer variable that denotes type of failure for each person.
}
  \item{\dots}{
variables passed to `crr' function; two key variables are \emph{failcode} and \emph{cencode}; see below in \code{Description}.
}
  \item{subset}{
subset of data to be used for model selection.
}
  \item{data}{
data-frame containing all the variables. Only complete cases are used in the analysis, i.e. rows of dataframe with missing values in any of the predictors are deleted.
}
  \item{direction}{
forward or backward direction for model selection.
}
  \item{criterion}{
selection criterion; default is AIC.  BIC uses log(n) as penalty, where `n' is total sample size, and BICcr uses log(n*) as the penalty where n* is the number of primary events.
}
  \item{crr.object}{
logical variable indicating whether to return final `crr' object.
}
  \item{trace}{
logical indicating whether to display stepwise model selection process.
}
  \item{steps}{
maximum number of steps in stepwise selection.
}
}
\details{
Based on the existing code of stepAIC in the MASS package. Variables passed to 'crr' function include two key variables: \emph{failcode} and \emph{cencode}. \emph{failcode} is an integer value that denotes primary failure, and \emph{cencode} is an integer denoting censoring event.
}
\value{
  \item{variables}{Variables in the final model}
  \item{coefficients}{The estimated coefficients of the variables}
  \item{std.errors}{Standard errors of the estimated coefficients}
  \item{log.lik}{The partial log-likelihood of the model}
}
\references{
Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for
the subdistribution of a competing risk. Journal of the American Statistical
Association.

Volinsky, C. T. and Raftery, A. E. (2000). Bayesian information criterion
for censored survival models. Biometrics.

}
\author{
Ravi Varadhan & Deborah Kuk.
}

\seealso{
\code{\link{crr}}
}
\examples{
library(cmprsk)

set.seed(123)
n <- 500
ftime <- rexp(n)
fstatus <- sample(0:2,n,replace=TRUE)
cov1 <- matrix(runif(5*n),nrow=n)
x61 <- as.factor(sample(3, size=n, rep=TRUE))
x71 <- as.factor(sample(5, size=n, rep = TRUE))
cov1 <- cbind(cov1, x61, x71)
dimnames(cov1)[[2]] <- c('x1','x2','x3','x4','x5', 'x6', 'x7')
formula1 <- ftime ~ 1 + x1 + x2 + x3 + x4 + x5 + as.factor(x6) + as.factor(x7) 

crrstep(formula1, , fstatus, data = as.data.frame(cov1), direction = "backward", criterion = "BIC")

ans <- crrstep(formula1, , fstatus, data = as.data.frame(cov1), direction = "forward", criterion = "AIC")

ans$coef

}
\keyword{stepwise}
\keyword{competing risks}
