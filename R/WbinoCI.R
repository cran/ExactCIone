lcp_bino<-function(x,n,alpha,prec2=1e-12){ # 1-alpha/2 lower CP interval
  if(length(x)>1){
    re<-x
    # the lower limit of 1-alpha CP interval using F
    ind_x<-which(x<0.5)
    re[ind_x]<-0
    re[-ind_x]<-(1+(n-x[-ind_x]+1)/(x[-ind_x]*stats::qf(alpha/2,2*x[-ind_x],2*(n-x[-ind_x]+1))))^(-1)-prec2

  }else{
    if(x<0.5) {re<-0}else{
      re<-(1+(n-x+1)/(x*stats::qf(alpha/2,2*x,2*(n-x+1))))^(-1)-prec2
    }
  }
  re
}

# the upper limit of 1-alpha CP interval
ucp_bino<- function(x,n,alpha,prec2=1e-12){ 1-lcp_bino(n-x,n,alpha,prec2)}
# indicator function of interval [a,b]
ind<- function(x,a,b){(x>=a)*(x<=b)}

# the coverage probability function for the interval [lcin,ucin].
cpcig_bino <- function(p,lcin,ucin,n){
  xx <- c(0:n)
  indp <- xx
  uu <- 1
  while (uu<length(xx)+0.5) {
    indp[uu]<-ind(p,lcin[uu],ucin[uu])*stats::dbinom(xx[uu],n,p)
    uu<-uu+1}
  re<-sum(indp)
  re
}

# The icp of
Icp<-function(CIM,prec2=1e-12,n){
  lcpw<-CIM[,1]
  ucpw<-CIM[,2]
  pp1<-which(lcpw>prec2)
  cp1<-min(sapply(lcpw[pp1]-prec2,cpcig_bino,lcin=lcpw,ucin=ucpw,n=n))
  cp1
}

# The main function
WbinoCITwo_core<-function(x,n,alpha){
  prec1<-1e-7
  prec2<-1e-12

  xx<-c(0:n)

  # The clopper-pearson interval
  lcpn<-lcp_bino(xx,n,alpha)   #the lower limit of CP interval
  ucpn<-ucp_bino(xx,n,alpha)   #the upper limit of CP interval

  # Initialize the lower and upper limit of Wang exact interval.
  lcpw<-lcpn                      #the lower limit of CP refinement, should be between lcpn and lcp1n   (2)
  ucpw<-ucpn                      #the upper limit of CP refinement, should be between ucpn and ucp1n
  lcp1n<-lcp_bino(xx,n,2*alpha)  #the lower limit of  one-sided CP interval >=lcpn
  ucp1n<-ucp_bino(xx,n,2*alpha)  #the upper limit of one-sided CP interval <=ucpn

  if(n/2==floor(n/2)){ # If n is odds
    ss<-n/2+1
    a0<-lcpn[ss]
    a1<-lcp1n[ss]  # a1>a0. the final lower limit is between a1 and a0.

    while (a1-a0>prec1){ # The refined interval of the longest interval
      lcpw0<-lcpw
      ucpw0<-ucpw

      aa<-(a0+a1)/2
      lcpw[ss]<-aa
      ucpw[ss]<-1-aa

      ind_ss<-which(lcpw[(ss+1):(n+1)]<aa)
      if(length(ind_ss)>0){
        lcpw[(ss+1):(n+1)][ind_ss]<-aa
        ucpw<-1-lcpw[(n+1):1]
      }

      L0<-aa
      pp<-c(L0-prec2,(ucpw+prec2)[which(ucpw<=L0 & ucpw>=a0)])   #(4)
      cp=0*pp
      kkk=1
      while (kkk<n+1+0.5){
        cp=cp+ind(pp,lcpw[kkk],ucpw[kkk])*stats::dbinom(kkk-1,n,pp)
        kkk=kkk+1
      }

      #cp<-sapply(pp,cpcig_bino,lcin=lcpw,ucin=ucpw,n=n)
      minic<-min(cp)
      if (minic>=1-alpha){a0<-aa}else{
        a1<-aa
        lcpw[(ss+1):(n+1)][ind_ss]<-lcpw0[(ss+1):(n+1)][ind_ss]
        ucpw<-1-lcpw[(n+1):1]
      }
    }

    lcpw[ss]<-a0
    ucpw[ss]<-1-lcpw[ss]

    ss<-n/2} else {ss<-floor(n/2)+1}

  maxi<-c(min(x+1,n-x+1):max(x+1,n-x+1))

  while (ss>min(x,n-x)){# construct the lower and upper limit from X=n/2 to 0

    a0<-lcpw[ss]
    a1<-lcp1n[ss]  # a1>a0. the final lower limit is between a1 and a0.
    b0<-ucp1n[ss]
    b1<-ucpw[ss]   # b1>b0. the final upper limit is between b1 and b0.

    while (a1-a0>prec1){ # solve the lower limit
      lcpw0<-lcpw
      ucpw0<-ucpw

      aa<-(a0+a1)/2
      lcpw[ss]<-aa
      ucpw[n+2-ss]<-1-aa
      ind_ss1<-which(lcpw[(ss+1):(n+2-ss-1)]<aa)
      if(length(ind_ss1)>0){
        lcpw[(ss+1):(n+2-ss-1)][ind_ss1]<-aa
        ucpw<-1-lcpw[(n+1):1]}


      L0<-aa
      pp<-c(L0-prec2,(ucpw+prec2)[which(ucpw<=L0 & ucpw>=a0)])  #(4)
      cp=0*pp
      kkk=1
      while (kkk<n+1+0.5){
        cp=cp+ind(pp,lcpw[kkk],ucpw[kkk])*stats::dbinom(kkk-1,n,pp)
        kkk=kkk+1
      }

      #cp<-sapply(pp,cpcig_bino,lcin=lcpw,ucin=ucpw,n=n)
      minic<-min(cp)
      if (minic>1-alpha){
         a0<-aa
      }else{
        a1<-aa
        lcpw[(ss+1):(n+2-ss-1)][ind_ss1]<-lcpw0[(ss+1):(n+2-ss-1)][ind_ss1]
        ucpw<-1-lcpw[(n+1):1]
       }
    }
    lcpw[ss]<-a0
    ucpw[n+2-ss]<-1-lcpw[ss]
    while (b1-b0>prec1){
      lcpw0<-lcpw
      ucpw0<-ucpw
      bb<-(b0+b1)/2
      ucpw[ss]<-bb
      lcpw[n+2-ss]<-1-bb
      ind_ss2<-which(lcpw[(n+2-ss+1):(n+1)]<1-bb)
      if(length(ind_ss2)>0){
        lcpw[(n+2-ss+1):(n+1)][ind_ss2]<-1-bb
        ucpw<-1-lcpw[(n+1):1]}

      u0<-bb;
      pp<-c(u0+prec2,(lcpw-prec2)[which(lcpw>=u0 & lcpw<=b1)])   #(4)
      cp=0*pp
      kkk=1
      while (kkk<n+1+0.5){
        cp=cp+ind(pp,lcpw[kkk],ucpw[kkk])*stats::dbinom(kkk-1,n,pp)
        kkk=kkk+1
      }

      #cp<-sapply(pp,cpcig_bino,lcin=lcpw,ucin=ucpw,n=n)
      minic<-min(cp)
      if (minic>=1-alpha){
        b1<-bb
      }else {
        b0<-bb
        lcpw[min((n+2-ss+1),n+1):(n+1)][ind_ss2]<-lcpw0[min((n+2-ss+1),n+1):(n+1)][ind_ss2]
        ucpw<-1-lcpw[(n+1):1]
      }
    }


    ucpw[ss]<-b1
    lcpw[n+2-ss]<-1-ucpw[ss]

    ss<-ss-1
  }
  re<-cbind(xx[maxi],lcpw[maxi],ucpw[maxi])
  colnames(re)<-c("x","lower","upper")
  re
}

#' An Admissible Exact Confidence Interval for the Bnomial Proportion
#' @description An admissible exact confidence interval of level 1-alpha is
#'   constructed for the binomial proportion p. This function can be used to
#'   calculate the interval constructed method proposed by Wang (2014).
#' @details Suppose X~bino(n,p), the sample space of X is \{0,1,...,n\}. Wang
#'   (2014) proposed an admissible interval which is obtained by uniformly
#'   shrinking the initial 1-alpha Clopper-Pearson interval from the middle to
#'   both sides of the sample space iteratively. This interval is admissible so
#'   that any proper sub-interval of it cannot assure the confidence coefficient.
#'   This means the interval cannot be shortened anymore.
#' @references Clopper, C. J. and Pearson, E. S. (1934). The use of confidence
#'   or fiducial limits in the case of the binomial. "Biometrika" 26: 404-413.
#' @references Wang, W. (2014). An iterative construction of confidence
#'   intervals for a proportion. "Statistica Sinica" 24: 1389-1410.
#' @param x the number of success or the observed data.
#' @param n the sample size.
#' @param conf.level Confidence level. The default is 0.95.
#' @param details TRUE/FALSE, can be abbreviated. To choose whether to compute
#'   the confidence interval for the whole sample points and output the infimum
#'   coverage probability. The default is FALSE.
#' @return A list which contains the confidence interval (CI) of the sample
#'   point and the confidence intervals (CIM) for all the points and the icp.
#' @export
#' @examples WbinoCI(x=2,n=5,conf.level=0.95,details=TRUE)
#' @examples WbinoCI(x=2,n=5,conf.level=0.95)
WbinoCI<-function(x,n,conf.level=0.95,details=FALSE){
  alpha<-1-conf.level
  if(details==TRUE){
    SSpace_CIM<-WbinoCITwo_core(0,n,alpha)
    icp<-Icp(SSpace_CIM[,-1],prec2=1e-12,n)
    # sample<-x
    re<-list(CI=cbind(x,t(SSpace_CIM[x+1,-1])),CIM=SSpace_CIM,icp=icp)
    return(re)

  }else{
    CI_x<-WbinoCITwo_core(x,n,alpha)
    CI<-cbind(t(CI_x[1,1]),t(CI_x[1,-1]))
    if(x>n/2) CI<-cbind(t(CI_x[dim(CI_x)[1],1]),t(CI_x[dim(CI_x)[1],-1]))
    re<-list(CI=CI)
    return(re)
  }
}



#' An Admissible Exact Lower Interval for the Binomial Proportion
#' @description The 1-alpha Clopper-Pearson lower interval for the binomial
#'   proportion p.
#' @param x the number of success or the observed data.
#' @param n the sample size.
#' @param conf.level Confidence level. The default is 0.95.
#' @param details TRUE/FALSE, can be abbreviated. To choose whether to compute
#'   the confidence interval for the whole sample points. The default is FALSE.
#' @return A list which contains the confidence interval (CI) of the sample
#'   point and the confidence intervals (CIM) for all the points.
#' @references Clopper, C. J. and Pearson, E. S. (1934). The use of confidence
#'   or fiducial limits in the case of the binomial. "Biometrika" 26: 404-413.
#' @export
#' @examples WbinoCI_lower(x=2,n=5,conf.level=0.95,details=TRUE)
#' @examples WbinoCI_lower(x=2,n=5,conf.level=0.95)
WbinoCI_lower<-function(x,n,conf.level=0.95,details=FALSE){
  alpha<-1-conf.level
  if(details==T){
    sample<-c(0:n)
    CIM<-cbind(lcp_bino(0:n,n,2*alpha,1e-12),1)
    colnames(CIM)<-c("lower","upper")
    #icp<-Icp(CIM,1e-12,n)
    re<-list(CI=cbind(sample,CIM))
  }else{
    CI<-cbind(lcp_bino(x,n,2*alpha,1e-12),1)
    colnames(CI)<-c("lower","upper")
    sample<-x
    re<-list(CI=cbind(sample,CI))
  }
  re
}

#' An Admissible Exact Upper Interval for the Binomial Proportion
#' @description The 1-alpha Clopper-Pearson upper interval for the binomial
#'   proportion p.
#' @param x the number of success or the observed data.
#' @param n the sample size.
#' @param conf.level Confidence level. The default is 0.95.
#' @param details TRUE/FALSE, can be abbreviated. To choose whether to compute
#'   the confidence interval for the whole sample points. The default is FALSE.
#' @return A list which contains the confidence interval (CI) of the sample
#'   point and the confidence intervals (CIM) for all the points.
#' @references Clopper, C. J. and Pearson, E. S. (1934). The use of confidence
#'   or fiducial limits in the case of the binomial. "Biometrika" 26: 404-413.
#' @export
#' @examples WbinoCI_upper(x=2,n=5,conf.level=0.95,details=TRUE)
#' @examples WbinoCI_upper(x=2,n=5,conf.level=0.95)
WbinoCI_upper<-function(x,n,conf.level=0.95,details=FALSE){
  alpha<-1-conf.level
  if(details==T){
    sample<-c(0:n)
    CIM<-cbind(0,ucp_bino(0:n,n,2*alpha,1e-12))
    colnames(CIM)<-c("lower","upper")
    #icp<-Icp(CIM,1e-12,n)
    re<-list(CI=cbind(sample,CIM))
  }else{
    CI<-cbind(0,ucp_bino(x,n,2*alpha,1e-12))
    colnames(CI)<-c("lower","upper")
    sample<-x
    re<-list(CI=cbind(sample,CI))
  }
  re
}
