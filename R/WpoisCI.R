# the lower limit of 1-alpha/2 one-sided CI based on sufficient statistic
lcp<- function(x,alpha,prec2=1e-10){
  if(length(x)>1){
    aa<-x
    ind_x<-which(x<0.5)
    aa[ind_x]<-0
    aa[-ind_x]<-stats::qchisq(alpha/2,2*x[-ind_x])/2-prec2/2
  }else{
    if(x<0.5){aa<-0}else{
      aa<-stats::qchisq(alpha/2,2*x)/2-prec2/2
    }
  }
  aa
}

# the upper limit of 1-alpha/2 one-sided CI based on sufficient statistic
ucp<- function(x,alpha,prec2=1e-10){
  aa<-stats::qchisq(1-alpha/2,2*(x+1))/2+prec2/2
  aa
}

# indicator function of interval [a,b]
ind<- function(x,a,b){(x>=a)*(x<=b)}

# the coverage probability function for the interval [lcin,ucin] with nondecreasing lcin and ucin.
cpcig<- function(lambda,lcin,ucin,sx){
  xx<-c(0:sx)
  kk<-1:length(lambda)
  for (i in kk){
    kkk<-xx+1
    indd<-ind(lambda[i],lcin[kkk],ucin[kkk])
    kkkk<-kkk[which(indd>0)]-1
    mink<-min(kkkk)
    maxk<-max(kkkk)
    kk[i]<-stats::ppois(maxk,lambda[i])-stats::ppois(mink-1,lambda[i]) # use the monotonicity property on confidence limits
  }
  kk}
# this is for an increasing lciw
cpcig1<- function(lambda,lcin,ucin,ss,sx){
  xx<-c(0:sx)
  kk<-1:length(lambda)
  for (i in kk){
    kkk<-xx[1:ss]+1
    indd<-ind(lambda[i],lcin[kkk],ucin[kkk])
    kkkk<-kkk[which(indd>0)]-1
    mink<-min(kkkk)
    maxk<-max(kkkk)
    kk[i]<-stats::ppois(maxk,lambda[i])-stats::ppois(mink-1,lambda[i]) # use the monotonicity property on confidence limits
  }
  kk}

# The main function for poisson mean
WpoisCI_core<-function(x,alpha){

  prec<-1e-8
  prec2<-1e-10

  x0<-x
  u_x<-ucp(x0,alpha/2)
  l_x<-lcp(x0,alpha/2)
  while (l_x<=u_x){
    l_x<-lcp(x0+1,alpha/2)
    x0<-x0+1}
  sx<-x0+1

  xx<-c(0:sx)
  lciw2<-lcp(xx,alpha)     #  1-alpha/2 lower interval, smaller
  lciw1<-lcp(xx,(2*alpha)) #  1-alpha lower interval, larger
  lciw<-lciw2              #  the new 1-alpha interval's lower limit \in [lciw2,lciw1]

  uciw2<-ucp(xx,alpha)     # 1-alpha/2 upper interval, larger
  uciw1<-ucp(xx,(2*alpha)) # 1-alpha upper interval, smaller
  uciw<-uciw2              # the new 1-alpha interval's upper limit \in [uciw1,uciw2]

  # result<-cbind(xx,lciw2,lciw,lciw1,uciw1,uciw,uciw2)


  #### compute uciw at x=0 (where ss=1)
  ss<-1
  d1<-uciw[ss]
  d0<-uciw1[ss]  # d1>d0. the final upper limit is between d1 and d0.
  while (d1-d0>prec) {
    da<-(d0+d1)/2
    uciw[ss]<-da
    U0<-da
    rr<-c(U0+prec2,(lciw-prec2)[which(lciw>=U0 & lciw<=d1)])                                 #(3)
    ff<-cpcig(rr,lciw,uciw,sx)
    if (min(ff)>=1-alpha) {d1<-da} else {d0<-da}
  }
  uciw[ss]<-d1
  lciw[ss]<-0
  ##### compute lciw and uciw for the observed x(<=sx), where ss=x+1(<=sx+1)
  ss<-ss+1
  #sss<-x+1
  sss<-x+1
  while (ss<sss+0.5){  # construct lower and upper limits for the sample <=x.

    a0<-lciw[ss]
    a1<-lciw1[ss]  # a1>a0. the final lower limit is between a1 and a0.
    #la0=log(a0);la1=log(a1);
    while (a1-a0>prec){ # solve the lower limit
      la<-(a0+a1)/2
      lciw[ss]<-la
      L0<-la
      rr<-c(L0-prec2,(uciw+prec2)[which(uciw<=L0 & uciw>=a0)])                                 #(3)
      ff<-cpcig1(rr,lciw,uciw,ss,sx)
      minic<-min(ff)
      if (minic>=1-alpha){a0=la}else{a1=la}
    }
    lciw[ss]<-a0

    b0<-max(uciw1[ss],uciw[ss-1]) # for an increasing uciw
    b1<-uciw[ss]                  # b1>b0. the final upper limit is between b1 and b0.

    while (b1-b0>prec){           # solve the upper limit
      lb<-(b0+b1)/2
      uciw[ss]<-lb
      U0<-lb
      rr<-c(U0+prec2,(lciw-prec2)[which(lciw>=U0 & lciw<=b1)])                                #(3)
      ff<-cpcig(rr,lciw,uciw,sx)
      minic<-min(ff)
      if (minic>=1-alpha){b1=lb}else{b0=lb}
    }
    uciw[ss]<-b1
    ss<-ss+1
  }
  # result=cbind(xx,lciw2,lciw,lciw1,uciw1,uciw,uciw2)
  # result=cbind(xx,lciw2,lciw,uciw,uciw2)
  # options(digits=10)
  rr1<-c((uciw+prec2)[1:(x+1)],(lciw-prec2)[2:(x+1)])                                #(3)
  ff<-cpcig(rr1,lciw,uciw,sx)
  icp<-min(ff)
  result<-cbind(xx,lciw,uciw)
  colnames(result)<-c("x","lower","upper")
  re<-list(result=result,icp=icp)
  re
}


#' An Admissible Exact Confidence Interval for the Poisson Mean
#' @description An admissible exact confidence interval for the Poisson mean.
#'   This function can be used to calculate the interval constructed method
#'   proposed by Wang (2014).
#' @details Suppose X~poi(lambda), the sample space of X is \{0,1,...\}. Wang
#'   (2014) proposed an admissible interval which is obtained by uniformly
#'   shrinking the initial 1-alpha Clopper-Pearson interval from 0 to the sample
#'   point of interest. This interval is admissible so that any proper
#'   sub-interval of it cannot assure the confidence coefficient. This means the
#'   interval cannot be shortened anymore.
#' @param x the sample or the observed point.
#' @param conf.level confidence level. The default is 0.95.
#' @param details TRUE/FALSE, can be abbreviated. To choose whether to compute
#'   the confidence intervals for all the sample points. Default is FALSE.
#' @references Wang, W. (2014). An iterative construction of confidence
#'   intervals for a proportion. "Statistica Sinica" 24: 1389-1410.
#' @return a list which contain the confidence interval and the ICP.
#' @export
#' @examples WpoisCI(1)
#' @examples WpoisCI(3,details = TRUE)
WpoisCI<-function(x,conf.level=0.95,details=FALSE){
  alpha<-1-conf.level
  if(details==T){
    # sample<-x
    CIM_icp<-WpoisCI_core(x,alpha)
    CIM<-CIM_icp$result
    CI<-CIM[x+1,-1]
    result<-list(CI=cbind(x,t(CI)),CIM=CIM[c(1:(x+1)),],icp=CIM_icp$icp)
  }else{
    # sample<-x
    CI<-WpoisCI_core(x,alpha)$result[x+1,-1]
    result<-list(CI=cbind(x,t(CI)))
  }
  result
}



#' An Admissible Exact One-sided Lower Confidence Interval for Poisson Mean
#' @description The 1-alpha Clopper-Pearson type lower interval for the Poisson
#'   mean.
#' @param x the sample or the observed point.
#' @param conf.level confidence level. The default is 0.95.
#' @param details TRUE/FALSE, can be abbreviated. To choose whether to compute
#'   the confidence intervals for all the sample points. Default is FALSE.
#' @references  Garwood, F. (1936). Fiducial Limits for the Poisson
#'   Distribution. "Biometrika" 28:	437-442.
#' @return a list which contain the one-sided lower confidence interval.
#' @export
#' @examples WpoisCI_lower(1)
#' @examples WpoisCI_lower(3,details = TRUE)
WpoisCI_lower<-function(x,conf.level=0.95,details=FALSE){
  alpha<-1-conf.level
  if(details==T){
    sample<-0:x
    CIM<-cbind(lcp(sample,2*alpha,prec2=1e-10),Inf)
    colnames(CIM)<-c("lower","upper")
    CI<-CIM[x+1,]
    result<-list(CI=cbind(x,t(CI)),CIM=cbind(sample,CIM))

  }else{
    sample<-x
    CI<-cbind(lcp(sample,2*alpha,prec2=1e-10),Inf)
    result<-list(CI=cbind(sample,CI))
  }
  result
}


#' An Admissible Exact One-sided Upper Confidence Interval for Poisson Mean
#' @description The 1-alpha Clopper-Pearson type upper interval for the Poisson
#'   mean.
#' @param x the sample or the observed point.
#' @param conf.level confidence level. The default is 0.95.
#' @param details TRUE/FALSE, can be abbreviated. To choose whether to compute
#'   the confidence intervals for all the sample points. Default is FALSE.
#' @references  Garwood, F. (1936). Fiducial Limits for the Poisson
#'   Distribution. "Biometrika" 28:	437-442.
#' @return a list which contain the one-sided upper confidence interval.
#' @export
#' @examples WpoisCI_upper(1)
#' @examples WpoisCI_upper(3,details = TRUE)
WpoisCI_upper<-function(x,conf.level=0.95,details=FALSE){
  alpha<-1-conf.level
  if(details==T){
    sample<-0:x
    CIM<-cbind(0,ucp(sample,2*alpha,prec2=1e-10))
    colnames(CIM)<-c("lower","upper")
    CI<-CIM[x+1,]
    result<-list(CI=cbind(x,t(CI)),CIM=cbind(sample,CIM))

  }else{
    sample<-x
    CI<-cbind(0,ucp(sample,2*alpha,prec2=1e-10))
    result<-list(CI=cbind(sample,CI))
  }
  result
}
