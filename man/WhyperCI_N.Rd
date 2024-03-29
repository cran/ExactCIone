% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WhyperCI_N.R
\name{WhyperCI_N}
\alias{WhyperCI_N}
\title{An Admissible Exact Confidence Interval for N, the Number of Balls in an Urn.}
\usage{
WhyperCI_N(x, n, M, conf.level, details = FALSE)
}
\arguments{
\item{x}{integer representing the number of white balls in the drawn balls.}

\item{n}{integer representing the number of balls we draw in the urn without
replacement, i.e., the sample size.}

\item{M}{the number of white balls in the urn.}

\item{conf.level}{the confidence level of confidence interval.}

\item{details}{TRUE/FALSE, can be abbreviate. If choose FALSE, the confidence
interval at the observed X will be returned. If choose TRUE, the confidence
intervals for all sample points and  the infimum coverage probability will
be returned. Default is FALSE.}
}
\value{
a list which contains i) the confidence interval for N and ii) the
infimum coverage probability of the intervals.
}
\description{
An admissible exact confidence interval for the number of balls
in an urn, which is the population number of a hypergeometric distribution. This
function can be used to calculate the interval constructed method proposed
by Wang (2015).
}
\details{
Suppose X~Hyper(M,N,n). When M and n are known, Wang (2015) construct an
admissible confidence interval for N by uniformly shrinking the initial 1-alpha
Clopper-Pearson type interval from  0 to min(M,n). This interval is admissible so
that any proper sub-interval of it cannot assure the confidence coefficient. This
means the interval cannot be shortened anymore.
}
\examples{
WhyperCI_N(10,50,800,0.95,details=TRUE)
WhyperCI_N(50,50,800,0.95)
}
\references{
Wang, W. (2015). Exact Optimal Confidence Intervals for
Hypergeometric Parameters. "Journal of the American Statistical
Association" 110 (512): 1491-1499.
}
