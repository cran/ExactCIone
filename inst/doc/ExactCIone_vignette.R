## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup--------------------------------------------------------------------
library(ExactCIone)

## -----------------------------------------------------------------------------
# Compute the 95% confidence interval when x=2, n=5.
WbinoCI(x=2,n=5,conf.level=0.95)

## -----------------------------------------------------------------------------
WbinoCI(x=2,n=5,conf.level=0.95,details=TRUE)

## -----------------------------------------------------------------------------
WbinoCI_lower(x=2,n=5,conf.level=0.95)
WbinoCI_lower(x=2,n=5,conf.level=0.95,details=TRUE)
WbinoCI_upper(x=2,n=5,conf.level=0.95)
WbinoCI_upper(x=2,n=5,conf.level=0.95,details=TRUE)

## -----------------------------------------------------------------------------
# The admissible CI for poisson mean when the observed sample is x=3. 
WpoisCI(x=3,conf.level = 0.95)
#We show the intervals from 0 to the sample of interest when "details=TRUE".
WpoisCI(x=3,conf.level = 0.95,details = TRUE)

## -----------------------------------------------------------------------------
WpoisCI_lower(x=3,conf.level = 0.95)
WpoisCI_lower(x=3,conf.level = 0.95,details = TRUE)
WpoisCI_upper(x=3,conf.level = 0.95)
WpoisCI_upper(x=3,conf.level = 0.95,details = TRUE)

## -----------------------------------------------------------------------------
# For hyper(M,N,n), construct 95% CI for N on the observed sample x when n,M are known.
WhyperCI_N(x=5,n=10,M=800,conf.level = 0.95)
# It shows CIs for all the sample point When "details=TRUE".
WhyperCI_N(x=5,n=10,M=800,conf.level = 0.95,details=TRUE)

## -----------------------------------------------------------------------------
WhyperCI_N_lower(x=0,n=10,M=800,conf.level = 0.95)
WhyperCI_N_lower(x=0,n=10,M=800,conf.level = 0.95,details=TRUE)
WhyperCI_N_upper(x=0,n=10,M=800,conf.level = 0.95)
WhyperCI_N_upper(x=0,n=10,M=800,conf.level = 0.95,details=TRUE)

## -----------------------------------------------------------------------------
# For Hyper(M,N,n), construct the CI for M on the observed sample x when n, N are known. 
# Also output CI for p=M/N.
WhyperCI_M(x=0,n=10,N=2000,conf.level = 0.95)
WhyperCI_M(x=0,n=10,N=2000,conf.level = 0.95,details = TRUE)

## -----------------------------------------------------------------------------
WhyperCI_M_lower(X=0,n=10,N=2000,conf.level = 0.95)
WhyperCI_M_lower(X=0,n=10,N=2000,conf.level = 0.95,details = TRUE)
WhyperCI_M_upper(X=0,n=10,N=2000,conf.level = 0.95)
WhyperCI_M_upper(X=0,n=10,N=2000,conf.level = 0.95,details = TRUE)


