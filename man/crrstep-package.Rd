\name{crrstep-package}
\alias{crrstep-package}
\docType{package}
\title{
Stepwise regression procedure for the Fine & Gray regression model in competing risks 
}
\description{
Performs forward and backward stepwise regression for the Fine & Gray regression model in competing risks. Procedure uses AIC, BIC and BICcr as selection criteria. BICcr has a penalty of k = log(n*), where n* is the number of Type I events. }
\details{
\tabular{ll}{
Package: \tab crrstep\cr
Type: \tab Package\cr
Version: \tab 2014-01.08\cr
Date: \tab 2014-01.08\cr
License: \tab GPL (version 2)\cr
LazyLoad: \tab yes\cr
}
The package contains a singe function \code{crrstep}, which implements backward and forward stepwise regression for the Fine & Gray regression model. The Fine & Gray model (Fine & Gray, 1999) estimates the hazard that corresponds to the cumulative incidence function of a certain event type. Selection criteria that are can be used are: AIC, BIC and BICcr. BICcr is a selection criteria based on the work by Volinksy and Raftery in which the penalty is k = log(n*), where n* is the total number of Type I events.
}
\author{
Ravi Varadhan & Deborah Kuk.

Maintainers: Ravi Varadhan <rvaradhan@jhmi.edu>
}
\references{
Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for
the subdistribution of a competing risk. \emph{Journal of the American Statistical
Association}.


Volinsky, C. T. and Raftery, A. E. (2000). Bayesian information criterion
for censored survival models. \emph{Biometrics}.


Kuk, D. and Varadhan, R. (2013). Model selection in competing risks regression. \emph{Statistics in Medicine}. 
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
%\keyword{package}
\examples{
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
crrstep(formula1, , fstatus, data = as.data.frame(cov1), direction = "backward", criterion = "AIC")

ans2 <- crrstep(formula1, , fstatus, data = as.data.frame(cov1), direction = "forward", 
		failcode=2, criterion = "AIC")
ans2
}
