# indicator fucntion of interval [a,b]
ind<- function(x,a,b){(x>=a)*(x<=b)}

ffnl<-function(N,n,M){1-stats::phyper(0,M,N-M,n)}

ffnu<-function(N,n,M){stats::phyper(0,M,N-M+1,n)}

ffnl2<-function(N,n,M){1-stats::phyper(0,M,N-M,n)}

ffnu2<-function(N,n,M){stats::phyper(0,M,N-M+1,n)}


# the lower limit of 1-2alpha CP interval
lcit<- function(x,n,M,alpha,Nl){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]>min(n,M)-0.5){
      kk[i]<-max(n,M)} else {
        aa<-max(M,n):Nl
        bb<-aa*0+1
        bb[2:length(aa)]<-1-stats::phyper(x[i],M,aa[2:length(aa)]-M-1,n)
        dd<-cbind(aa,bb)
        dd<-dd[which(dd[,2]>=1-alpha),]
        if (length(dd)==2){
          kk[i]<-dd[1]} else {kk[i]<-max(dd[,1])}
      }
  }
  kk
}

# The upper limit of 1-2alpha CP interval
ucit<- function(x,n,M,alpha,Nu){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]<0.5){kk[i]=Inf} else {
      aa<-max(M,n):(Nu+1)
      bb<-stats::phyper(x[i]-1,M,aa[1:length(aa)]-M+1,n)
      dd<-cbind(aa,bb)
      dd<-dd[which(dd[,2]>=1-alpha),]
      if (length(dd)==2){kk[i]<-dd[1]} else {kk[i]<-min(dd[,1])}
    }
  }
  kk
}


# The upper limit of 1-2alpha CP interval, same as ucit
ucitnew<- function(x,n,M,alpha,Nu2,step){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]<0.5){kk[i]=Inf} else {
      aa<-seq(max(M,n),(Nu2+step),by=step)
      bb<-stats::phyper(x[i]-1,M,aa[1:length(aa)]-M+1,n)
      dd<-cbind(aa,bb)
      dd<-dd[which(dd[,2]>=1-alpha),]
      if (length(dd)==2){kk[i]<-dd[1]} else {
        kk[i]<-min(dd[,1])
      }
      aa1<-max(max(M,n),(kk[i]-step-1)):kk[i]
      bb1<-stats::phyper(x[i]-1,M,aa1[1:length(aa1)]-M+1,n)
      dd1<-cbind(aa1,bb1)
      dd1<-dd1[which(dd1[,2]>=1-alpha),]
      if (length(dd1)==2){kk[i]<-dd1[1]} else {kk[i]<-min(dd1[,1])}
    }
  }
  kk
}

# The lower limit of 1-alpha CP interval
lcit2<- function(x,n,M,alpha,Nl2){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]>min(n,M)-0.5){kk[i]<-max(n,M)} else {
      aa<-max(M,n):Nl2
      bb<-aa*0+1
      bb[2:length(aa)]<-1-stats::phyper(x[i],M,aa[2:length(aa)]-M-1,n)
      dd<-cbind(aa,bb)
      dd<-dd[which(dd[,2]>=1-alpha/2),]
      if (length(dd)==2){kk[i]<-dd[1]} else {kk[i]<-max(dd[,1])}
    }
  }
  kk
}

# the upper limit of 1-alpha CP interval
ucit2<- function(x,n,M,alpha,Nu2){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]<0.5){kk[i]=Inf} else {
      aa<-max(M,n):(Nu2+1)
      bb<-stats::phyper(x[i]-1,M,aa[1:length(aa)]-M+1,n)
      dd<-cbind(aa,bb)
      dd<-dd[which(dd[,2]>=1-alpha/2),]
      if (length(dd)==2){kk[i]<-dd[1]} else {kk[i]<-min(dd[,1])}
    }
  }
  kk
}


# the upper limit of 1-alpha CP interval, same as ucit2
ucit2new<- function(x,n,M,alpha,Nu2,step){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]<0.5){kk[i]=Inf} else {
      aa<-seq(max(M,n),(Nu2+step),by=step)
      bb<-stats::phyper(x[i]-1,M,aa[1:length(aa)]-M+1,n)
      dd<-cbind(aa,bb)
      dd<-dd[which(dd[,2]>=1-alpha/2),]
      if (length(dd)==2){kk[i]<-dd[1]} else {kk[i]<-min(dd[,1])}
      aa1<-max(max(M,n),(kk[i]-step-1)):kk[i]
      bb1<-stats::phyper(x[i]-1,M,aa1[1:length(aa1)]-M+1,n)
      dd1<-cbind(aa1,bb1)
      dd1<-dd1[which(dd1[,2]>=1-alpha/2),]
      if (length(dd1)==2){kk[i]<-dd1[1]} else {kk[i]<-min(dd1[,1])}
    }
  }
  kk
}

# The coverage probability
cpcit<- function(N,lcin,ucin,n,M){
  kk<-1:length(N)
  for (i in kk){
    indp<-0:min(n,M)
    uu<-0
    while (uu<min(n,M)+0.5) {
      indp[uu+1]<-ind(N[i],lcin[uu+1],ucin[uu+1])*stats::dhyper(uu,M,N[i]-M,n)
      uu<-uu+1}
    kk[i]<-sum(indp)
  }
  kk
}

# the minimum of cpcit
cpcits<- function(lcin,ucin,n,M){
  kk<-length(lcin)
  N11<-1:(kk-1)
  for (ii in 1:(kk-1)){N11[ii]<-max(lcin[ii]-1,max(n,M))}
  N12<-ucin[2:kk]+1
  N1<-c(N11,N12)
  ccp<-N1
  for (i in 1:length(N1)){ #%
    indp<-0:(kk-1)
    uu<-0
    while (uu<kk-1+0.5) {
      indp[uu+1]<-ind(N1[i],lcin[uu+1],ucin[uu+1])*stats::dhyper(uu,M,N1[i]-M,n)
      uu<-uu+1
    }
    ccp[i]<-sum(indp)
  } #%
  min(ccp)
}

# The main function
WhyperCITwo_N_core<-function(x,n,M,alpha){

  # General settings
  prec2<-1 # do not change it
  step1<-100
  step<-1000  # step should be large when M and n are large, e.g., M=n=250, step=1000
  minx<-min(M,n)
  maxx<-max(M,n)


  # An upper bound for the lower limit for N
  ii<-maxx
  ffnL<-ffnl(ii,n,M)
  while(ffnL>1-alpha){ffnL<-ffnl(ii,n,M);ii<-ii+step1}
  Nl<-ii            # L(X)<=Nl

  # An upper bound for the upper limit for N
  ii<-maxx
  ffnU<-ffnu(ii,n,M)
  while(ffnU<1-alpha){ffnU<-ffnu(ii,n,M);ii<-ii+step1}
  Nu<-ii            # U(X)<=NU for X>0

  # An upper bound for the lower limit for N
  ii<-maxx
  ffnL2<-ffnl2(ii,n,M)
  while(ffnL2>1-alpha/2){ffnL2<-ffnl2(ii,n,M);ii<-ii+step1}
  Nl2<-ii            # L(X)<=Nl2

  # An upper bound for the upper limit for N
  ii<-maxx
  ffnU2<-ffnu2(ii,n,M)
  while(ffnU2<1-alpha/2){ffnU2<-ffnu2(ii,n,M);ii<-ii+step1}
  Nu2<-ii            # U(X)<=NU2 for X>0


  # The sample space
  xx<-0:min(n,M)

  lcin1<-lcit2(xx,n,M,alpha,Nl2)   # this is the lower one-sided 1-alpha/2 interval
  lcin2<-lcit(xx,n,M,alpha,Nl)    # this is the lower one-sided 1-alpha interval
  ucin1<-ucit2new(xx,n,M,alpha,Nu2,step)   # this is the upper one-sided 1-alpha/2 interval
  ucin2<-c(ucin1[2:(min(n,M)+1)],ucit(min(n,M),n,M,alpha,Nu))            # <=ucpn

  # so [lcin2,ucin2] is a two-sided 1-alpha interval
  lciw<-lcin1  #  lcin1<=lciw<=lcin2
  uciw<-ucin1  #  ucin2<=uciw<=ucin1

  #  Need to compute lciw[1] (ss=0) since uciw[1]=Inf.
  ss<-1
  a0<-lciw[ss]
  a1<-lcin2[ss]  # a1>a0
  while (a1-a0>1.5){  #$
    aa<-ceiling((a0+a1)/2)
    lciw[ss]<-aa
    N1<-c(max(aa-prec2,maxx),(uciw+prec2)[which(uciw<=aa & uciw>=lcin1[ss]-prec2)])   #(4)
    cpw<-min(cpcit(N1,lciw,uciw,n,M))
    if (cpw>=1-alpha){a0<-aa} else {a1<-aa}
  } #$

  lciw[ss]<-a1
  N1<-c(max(a1-prec2,maxx),(uciw+prec2)[which(uciw<=a1 & uciw>=lcin1[ss]-prec2)])   #(4)
  cpw<-min(cpcit(N1,lciw,uciw,n,M))
  if (cpw<1-alpha){lciw[ss]<-a0}

  # Compute lciw[ii] and uciw[ii] for min(n,M)+1>=ii>=2 when ii goes large
  # ii=ss+1. ss=1 iff ii=2
  ss<-ss+1
  while (ss<x+1.5){ ##
    #while (ss<42.5){ ##
    b0<-ucin2[ss]
    b1<-uciw[ss]
    while (b1-b0>1.5){  #$
      bb<-ceiling((b0+b1)/2)
      uciw[ss]<-bb
      N1<-c(bb+prec2,(lciw-prec2)[which(lciw>=bb & lciw<=b1)])   #(4)
      cpw<-min(cpcit(N1,lciw,uciw,n,M))
      if (cpw>=1-alpha){b1<-bb} else {b0<-bb}
    } #$

    uciw[ss]<-b0            # this may be improved. Here b1-b0<=1
    N1<-(lciw-prec2)[which(lciw>=b0 & lciw<=b1)]
    N1<-c(b0+prec2,N1[which(N1>=maxx)])   #(4)
    cpw<-min(cpcit(N1,lciw,uciw,n,M))
    if (cpw<1-alpha){uciw[ss]<-b1}

    a0<-lciw[ss]
    a1<-lcin2[ss]
    while (a1-a0>1.5){  #$
      aa<-ceiling((a0+a1)/2)
      lciw[ss]<-aa
      N1<-c(max(aa-prec2,maxx),(uciw+prec2)[which(uciw<=aa & uciw>=lcin1[ss]-prec2)])   #(4)
      cpw<-min(cpcit(N1,lciw,uciw,n,M))
      if (cpw>=1-alpha){a0<-aa} else {a1<-aa}
    } #$

    lciw[ss]<-a1
    N1<-c(max(a1-prec2,maxx),(uciw+prec2)[which(uciw<=a1 & uciw>=lcin1[ss]-prec2)])   #(4)
    cpw<-min(cpcit(N1,lciw,uciw,n,M))
    if (cpw<1-alpha){lciw[ss]<-a0}
    ss<-ss+1
  }  ##

  re<-cbind(xx[1:(x+1)],lciw[1:(x+1)],uciw[1:(x+1)])
  colnames(re)<-c("x","lower","upper")
  re
}


#' An Admissible Exact Confidence Interval for N, the Number of Balls in an Urn.
#'
#' @description An admissible exact confidence interval for the number of balls
#'   in an urn, which is the population number of a hypergeometric distribution. This
#'   function can be used to calculate the interval constructed method proposed
#'   by Wang (2015).
#' @details Suppose X~Hyper(M,N,n). When M and n are known, Wang (2015) construct an
#' admissible confidence interval for N by uniformly shrinking the initial 1-alpha
#' Clopper-Pearson type interval from  0 to min(M,n). This interval is admissible so
#' that any proper sub-interval of it cannot assure the confidence coefficient. This
#' means the interval cannot be shortened anymore.
#' @param x integer representing the number of white balls in the drawn balls.
#' @param n integer representing the number of balls we draw in the urn without
#'   replacement, i.e., the sample size.
#' @param M the number of white balls in the urn.
#' @param conf.level the confidence level of confidence interval.
#' @param details TRUE/FALSE, can be abbreviate. If choose FALSE, the confidence
#'   interval at the observed X will be returned. If choose TRUE, the confidence
#'   intervals for all sample points and  the infimum coverage probability will
#'   be returned. Default is FALSE.
#' @references Wang, W. (2015). Exact Optimal Confidence Intervals for
#'   Hypergeometric Parameters. "Journal of the American Statistical
#'   Association" 110 (512): 1491-1499.
#' @return a list which contains i) the confidence interval for N and ii) the
#'   infimum coverage probability of the intervals.
#' @export
#' @examples WhyperCI_N(10,50,800,0.95,details=TRUE)
#' @examples WhyperCI_N(50,50,800,0.95)
WhyperCI_N<-function(x,n,M,conf.level,details=FALSE){
  alpha<-1-conf.level
  if(details==TRUE){
    SSpace_CIM<-WhyperCITwo_N_core(min(n,M),n,M,alpha)
    icp<-cpcits(SSpace_CIM[,2],SSpace_CIM[,3],n,M)
    #sample<-x
    re<-list(CI=cbind(x,t(SSpace_CIM[x+1,-1])),CIM=SSpace_CIM,icp=icp)
    return(re)

  }else{
    CI_x<-WhyperCITwo_N_core(x,n,M,alpha)
    CI<-cbind(t(CI_x[x+1,1]),t(CI_x[x+1,-1]))
    re<-list(CI=CI)
    return(re)
  }

}

#' An Admissible Exact One-sided Lower Interval for the Population Number of
#' Hypergeometric Distribution
#' @description The 1-alpha Clopper-Pearson type lower interval for the
#'   population number of hypergeometric distribution.
#' @references Konijn, H. S. (1973). Statistical Theory of Sample Survey Design
#'   and Analysis, Amsterdam: North-Holland.
#' @param x integer representing the number of white balls we observed when
#'   drawn without replacement from an urn which contains both black and white
#'   balls.
#' @param n the number we drawn.
#' @param M the number of the white balls.
#' @param conf.level the confidence level of confidence interval.
#' @param details TRUE/FALSE, can be abbreviate. Default is FALSE. If choose
#'   TRUE, the confidence intervals for the whole sample space will be returned.
#' @return a list which contains the confidence interval.
#' @export
#' @examples WhyperCI_N_lower(0,50,800,0.95,details=TRUE)
#' @examples WhyperCI_N_lower(0,50,800,0.95)
WhyperCI_N_lower<-function(x,n,M,conf.level,details=FALSE){
  alpha<-1-conf.level

  # An upper bound for the lower limit for N
  ii<-max(n,M)
  ffnL<-ffnl(ii,n,M)
  while(ffnL>1-alpha){ffnL<-ffnl(ii,n,M);ii<-ii+100}
  Nl<-ii            # L(X)<=Nl

  if(details==TRUE){
    # The sample space
    xx<-0:min(n,M)
    SSpace_CIM<-cbind(xx,lcit(xx,n,M,alpha,Nl),Inf)
    colnames(SSpace_CIM)<-c("sample","lower","upper")
    # icp<-cpcits(SSpace_CIM[,2],SSpace_CIM[,3],n,M)
    sample<-x
    # re<-list(CI=cbind(sample,t(SSpace_CIM[x+1,-1])),CIM=SSpace_CIM,icp=icp)
    re<-list(CI=cbind(sample,t(SSpace_CIM[x+1,-1])),CIM=SSpace_CIM)
    return(re)

  }else{
    CI<-c(lcit(x,n,M,alpha,Nl),Inf)
    re<-list(CI=cbind(x,t(CI)))
    return(re)
  }
}


#' An Admissible Exact One-sided Upper Interval for the Population Number of Hypergeometric Distribution
#' @description The 1-alpha Clopper-Pearson type upper interval for the
#'   population number of hypergeometric distribution.
#' @references Konijn, H. S. (1973). Statistical Theory of Sample Survey Design
#'   and Analysis, Amsterdam: North-Holland.
#' @param x integer representing the number of white balls we observed when drawn without replacement from an urn which contains both black and white balls.
#' @param n the number we drawn.
#' @param M the number of the white balls.
#' @param conf.level the confidence level of confidence interval.
#' @param details TRUE/FALSE, can be abbreviate. Default is FALSE. If choose TRUE, the confidence intervals for the whole sample space will be returned.
#' @return a list which contains the confidence interval.
#' @export
#' @examples WhyperCI_N_upper(0,50,800,0.95,details=TRUE)
#' @examples WhyperCI_N_upper(0,50,800,0.95)
WhyperCI_N_upper<-function(x,n,M,conf.level,details=FALSE){
  alpha<-1-conf.level

  # An upper bound for the upper limit for N
  ii<-max(n,M)
  ffnU2<-ffnu2(ii,n,M)
  while(ffnU2<1-alpha/2){ffnU2<-ffnu2(ii,n,M);ii<-ii+100}
  Nu2<-ii            # U(X)<=NU2 for X>0

  if(details==TRUE){
    # The sample space
    xx<-0:min(n,M)
    SSpace_CIM<-cbind(xx,0,ucit2new(xx,n,M,2*alpha,Nu2,step=1000))
    colnames(SSpace_CIM)<-c("sample","lower","upper")
    # icp<-cpcits(SSpace_CIM[,2],SSpace_CIM[,3],n,M)
    sample<-x
    # re<-list(CI=cbind(sample,t(SSpace_CIM[x+1,-1])),CIM=SSpace_CIM,icp=icp)
    re<-list(CI=cbind(sample,t(SSpace_CIM[x+1,-1])),CIM=SSpace_CIM)
    return(re)

  }else{
    CI<-c(0,ucit2new(x,n,M,2*alpha,Nu2,step=1000))
    re<-list(CI=cbind(x,t(CI)))
    return(re)
  }
}
