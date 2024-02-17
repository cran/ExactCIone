# indicator fucntion of interval [a,b]
ind<- function(x,a,b){(x>=a)*(x<=b)}

# the lower limit of 1-alpha CP interval
lci<- function(x,n,N,alpha){
  kk<-1:length(x)
  for (i in kk){
    if (x[i]<0.5){kk[i]<-0} else {
      aa<-0:N
      bb<-aa+1
      bb[2:(N+1)]<-stats::phyper(x[i]-1,aa[2:(N+1)]-1,N-aa[2:(N+1)]+1,n)
      dd<-cbind(aa,bb)
      dd<-dd[which(dd[,2]>=1-alpha),]
      if (length(dd)==2){kk[i]<-dd[1]} else {kk[i]<-max(dd[,1])}
    }
  }
  kk
}

# upper limit of 1-alpha CP interval
uci<- function(x,n,N,alpha){ N-lci(n-x,n,N,alpha)}

# the coverage probability function for the combined interval.
cpci<- function(M,n,N,lcl,ucl){
  kk<-1:length(M)
  for (i in kk){
    xx<- 0:n
    indp<-xx
    uu<-0
    while (uu<n+0.5) {
      indp[uu+1]<-ind(M[i],lcl[uu+1],ucl[uu+1])*stats::dhyper(uu,M[i],N-M[i],n)
      uu<-uu+1}
    kk[i]<-sum(indp)
  }
  kk
}

# The main function
WhyperCI_M_core<-function(N,n,X,alpha){
  # General settings
  prec2=1
  xx=0:n
  lcin1<-lci(xx,n,N,alpha/2);ucin1<-uci(xx,n,N,alpha/2) #two-sided 1-alpha interval
  lcin2<-lci(xx,n,N,alpha);  ucin2<-uci(xx,n,N,alpha)   #two-sided 1-2alpha interval
  lciw<-lcin1; uciw<-ucin1  #  lcin1<=lciw<=lcin2; ucin2<=uciw<=ucin1

  if(n/2==floor(n/2)){#**  # for an even number n
    xvalue<-n/2+1  # start from the center
    a0<-lcin1[xvalue]
    a1<-lcin2[xvalue]  #a1>a0. the final lower limit is between a1 and a0.
    aa<-a0:a1
    ii<-1
    ii1<-ii
    while (ii <length(aa)+0.5){ #ii
      lciw[xvalue]<-aa[ii]
      uciw[xvalue]<-N-aa[ii]
      L0<-aa[ii]
      M<-c(L0-prec2,(uciw+prec2)[which(uciw<=min(L0,N-1) & uciw>=a0)])   #(4)
      bb<-min(cpci(M,n,N,lciw,uciw))
      if (bb>=1-alpha){ii1<-ii;ii<-ii+1} else {ii<-length(aa)+1}
    }                           #ii
    lciw[xvalue]<-aa[ii1]
    uciw[xvalue]<-N-lciw[xvalue]

    xvalue<-xvalue-1  # xvalue >=1
    while (xvalue>min(X,n-X)+0.5){  ##
      a0<-lcin1[xvalue]
      a1<-lcin2[xvalue] # a1>a0. the final lower limit is between a1 and a0.
      aa<-a0:a1   # solve the lower limit
      ii<-1
      ii1<-ii
      while (ii<length(aa)+0.5){ #i
        lciw[xvalue]<-aa[ii]
        uciw[n+2-xvalue]<-N-aa[ii]
        L0<-aa[ii]
        M<-c(max(0,L0-prec2),(uciw+prec2)[which(uciw<=min(L0,N-1) & uciw>=a0)])   #(4)
        ppu<-min(cpci(M,n,N,lciw,uciw))
        if (ppu>=1-alpha){ii1<-ii; ii<-ii+1} else{ii<-length(aa)+1}
      }  #i
      lciw[xvalue]<-aa[ii1]
      uciw[n+2-xvalue]<-N-lciw[xvalue]

      b0<-ucin2[xvalue]
      b1<-ucin1[xvalue] # b1>b0. the final lower limit is between b1 and b0.
      bb<-b1:b0   # solve the upper limit
      ii<-1
      ii0<-ii
      while (ii<length(bb)+0.5){
        uciw[xvalue]<-bb[ii]
        lciw[n+2-xvalue]<-N-bb[ii]
        u0<-bb[ii]
        M<-c(min(u0+prec2,N),(lciw-prec2)[which(lciw>=max(u0,1) & lciw<=b1)])   #(4)
        ppu<-min(cpci(M,n,N,lciw,uciw))
        if (ppu>=1-alpha){ii0<-ii; ii<-ii+1} else{ii<-length(bb)+1}
      }
      uciw[xvalue]<-bb[ii0]
      lciw[n+2-xvalue]<-N-bb[ii0]
      xvalue<-xvalue-1
    }
  }else{
    # for an odd number n
    xvalue<-floor(n/2)+1  # xvalue >=1
    while (xvalue>min(X,n-X)+0.5){  ##
      a0<-lcin1[xvalue]
      a1<-lcin2[xvalue] # a1>a0. the final lower limit is between a1 and a0.
      aa<-a0:a1   # solve the lower limit
      ii<-1
      ii0<-ii
      while (ii<length(aa)+0.5){ #i
        lciw[xvalue]<-aa[ii]
        uciw[n+2-xvalue]<-N-aa[ii]
        L0<-aa[ii]
        M<-c(max(0,L0-prec2),(uciw+prec2)[which(uciw<=min(L0,N-1) & uciw>=a0)])   #(4)
        ppu<-min(cpci(M,n,N,lciw,uciw))
        if (ppu>=1-alpha){ii1<-ii; ii<-ii+1} else{ii<-length(aa)+1}
      }  #i
      lciw[xvalue]<-aa[ii1]
      uciw[n+2-xvalue]<-N-lciw[xvalue]

      b0<-ucin2[xvalue]
      b1<-ucin1[xvalue] # b1>b0. the final lower limit is between b1 and b0.
      bb<-b1:b0   # solve the upper limit
      ii<-1
      ii0<-ii
      while (ii<length(bb)+0.5){
        uciw[xvalue]<-bb[ii]
        lciw[n+2-xvalue]<-N-bb[ii]
        u0<-bb[ii]
        M<-c(min(u0+prec2,N),(lciw-prec2)[which(lciw>=max(u0,1) & lciw<=b1)])   #(4)
        ppu<-min(cpci(M,n,N,lciw,uciw))
        if (ppu>=1-alpha){ii0<-ii; ii<-ii+1} else{ii<-length(bb)+1}
      }
      uciw[xvalue]<-bb[ii0]
      lciw[n+2-xvalue]<-N-bb[ii0]
      xvalue<-xvalue-1
    }
  }#***

  indx<-c(min(X+1,n-X+1):max(X+1,n-X+1))
  result<-cbind(xx[indx],lciw[indx],uciw[indx])
  colnames(result)<-c("x","lower","upper")
  result
}


#' An Admissible Exact Confidence Interval for M, the Number of White Balls in
#' an Urn
#'
#' @description The confidence interval for the number of white balls in an urn
#'   that contains M white balls and  N-M black balls when sampling without
#'   replacement. This function can be used to calculate the interval
#'   constructed method proposed by Wang (2015).
#' @details Suppose X~Hyper(M,N,n). When N and n are known, Wang (2015) construct an
#' admissible confidence interval for N by uniformly shrinking the initial 1-alpha
#' Clopper-Pearson type interval from  the mid-point of the sample space to 0. This interval
#' is admissible so that any proper sub-interval of it cannot assure the confidence coefficient. This
#' means the interval cannot be shortened anymore.
#'
#' @references Wang, W. (2015). Exact Optimal Confidence Intervals for
#'   Hypergeometric Parameters. "Journal of the American Statistical
#'   Association" 110 (512): 1491-1499.
#' @param N integer representing the number of all balls in an urn, i.e., the
#'   population size.
#' @param n integer representing the number of balls we draw in the urn without
#'   replacement, i.e., the sample size.
#' @param x integer representing the number of white balls in the drawn balls.
#' @param conf.level the confidence level of confidence interval.
#' @param details TRUE/FALSE, can be abbreviate. If choose FALSE, the confidence
#'   interval at the observed X will be returned. If choose TRUE, the confidence
#'   intervals for all sample points and  the infimum coverage probability will
#'   be returned. Default is FALSE.
#' @return a list which contains i) the confidence interval for M, ii)the
#'   confidence interval for p=M/N (this interval is equal to the previous
#'   interval divided by N) and iii) the infimum coverage probability of the two
#'   intervals.
#' @export
#' @examples WhyperCI_M(0,50,2000,0.95,details = TRUE)
#' @examples WhyperCI_M(0,50,2000,0.95)
WhyperCI_M<-function(x,n,N,conf.level,details=FALSE){
  X<-x
  alpha<-1-conf.level
  if(details==TRUE){
    SSpace_CIM<-WhyperCI_M_core(N,n,0,alpha)
    CIM_p<-cbind(SSpace_CIM[,1]/n,SSpace_CIM[,c(2,3)]/N)
    colnames(CIM_p)<-c("p","lower_p","upper_p")
    M<-SSpace_CIM[,2]
    M<-M[which(M>0)]-1
    icp<-min(cpci(M,n,N,SSpace_CIM[,2],SSpace_CIM[,3]))           #the ICP
    re<-list(CI=cbind(X,t(SSpace_CIM[X+1,-1])),CIM=cbind(SSpace_CIM),CIM_p=CIM_p,icp=icp)
  }else{
    CI_x<-WhyperCI_M_core(N,n,X,alpha)
    CI<-cbind(t(CI_x[1,1]),t(CI_x[1,-1]))
    if(X>n/2) CI<-cbind(t(CI_x[dim(CI_x)[1],1]),t(CI_x[dim(CI_x)[1],-1]))
    CI_p<-CI/N
    p<-X/n
    re<-list(CI=CI,CI_p=cbind(p,t(CI_p[,-1])))
  }
  return(re)

}

#' An Admissible Exact One-sided Lower Interval for the Number of White Balls in Hypergeometric Distribution
#' @description The 1-alpha Clopper-Pearson type lower interval for the
#'   number of white balls in an urn.
#' @param N integer representing the number of the whole balls in an urn.
#' @param n the number we drawn.
#' @param X integer representing the number of white balls we observed when drawn without replacement from an urn which contains both black and white balls.
#' @param conf.level the confidence level of confidence interval.
#' @param details TRUE/FALSE, can be abbreviate. Default is FALSE. If choose TRUE, the confidence intervals for the whole sample space and  the icp will be returned.
#' @references Konijn, H. S. (1973). Statistical Theory of Sample Survey Design
#'   and Analysis, Amsterdam: North-Holland.
#' @return a list which contains the confidence interval.
#' @export
#' @examples WhyperCI_M_lower(0,50,2000,0.95,details = TRUE)
#' @examples WhyperCI_M_lower(0,50,2000,0.95)
WhyperCI_M_lower<-function(X,n,N,conf.level,details=FALSE){
  alpha<-1-conf.level
  if(details==TRUE){
    xx=0:n
    lcin2<-lci(xx,n,N,alpha)
    result<-cbind(xx,lcin2,N)
    colnames(result)<-c("sample","lower","upper")
    re<-list(CI=cbind(X,lcin2[X+1],N),CIM=result)
    return(re)
  }else{
    re<-list(CI=cbind(X,lci(X,n,N,alpha),N))
    return(re)
  }
}

#' An Admissible Exact One-sided Upper Interval for the Number of White Balls in Hypergeometric Distribution
#' @description The 1-alpha Clopper-Pearson type upper interval for the
#'   number of white balls in an urn.
#' @param N integer representing the number of the whole balls in an urn.
#' @param n the number we drawn.
#' @param X integer representing the number of white balls we observed when drawn without replacement from an urn which contains both black and white balls.
#' @param conf.level the confidence level of confidence interval.
#' @param details TRUE/FALSE, can be abbreviate. Default is FALSE. If choose TRUE, the confidence intervals for the whole sample space and  the icp will be returned.
#' @references Konijn, H. S. (1973). Statistical Theory of Sample Survey Design
#'   and Analysis, Amsterdam: North-Holland.
#' @return a list which contains the confidence interval.
#' @export
#' @examples WhyperCI_M_upper(0,50,2000,0.95,details = TRUE)
#' @examples WhyperCI_M_upper(0,50,2000,0.95)
WhyperCI_M_upper<-function(X,n,N,conf.level,details=FALSE){
  alpha<-1-conf.level
  if(details==TRUE){
    xx<-0:n
    ucin2<-uci(xx,n,N,alpha)
    result<-cbind(xx,0,ucin2)
    colnames(result)<-c("sample","lower","upper")
    re<-list(CI=cbind(X,0,ucin2[X+1]),CIM=result)
    return(re)
  }else{
    re<-list(CI=cbind(X,0,uci(X,n,N,alpha)))
    return(re)
  }
}


