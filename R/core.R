#2020.01.27


library(pROC)


mROC_class_template<-list(p=NA,FPs=NA,TPs=NA)

class(mROC_class_template)<-"mROC"


#' @export
plot.mROC<-function(mROC_obj,...)
{
  #sf<-stepfun(mROC_obj$FPs,c(0,mROC_obj$TPs))
  plot(mROC_obj$FPs,mROC_obj$TPs,xlim=c(0,1),ylim=c(0,1),type='l',xlab="Falst Positive",ylab="True Positive") #TODO: let possible xlim and ylim from ... override the default
}



#' Takes in a vector of probabilities and returns mROC values (Tp,TP, in an object of class mROC)
#' @param p A numberic vector of probabilities.
#' @param ordered Optional, if the vector p is ordered from small to large (if not the function will do it; TRUE is to facilitate fast computations).
#' @return This function returns an object of class mROC.
#' @examples
#' @export
mROC<-function(p, ordered=F)
{
  if(min(p)<0 || max(p)>1) {stop("Error: invalid probability vector."); return(-1); }
  if(!ordered) p<-p[order(p)]

  sumP1<-sum(p)
  sumP0<-sum(1-p)

  n<-length(p)

  tp<-0
  fp<-0

  roc<-cbind(fp=rep(NA,length(p)),tp=rep(NA,length(p)))

  for(i in n:1)
  {
    fp<-fp+(1-p[i])/sumP0
    tp<-tp+p[i]/sumP1
    roc[n-i+1,]<-c(fp,tp)
  }

  out<-mROC_class_template
  out$p=p
  out$FPs<-c(0,roc[,1])
  out$TPs<-c(0,roc[,2])
  
  return(out)
}








#' Takes in a mROC object and calculates the area under the curve
#' @param mROC_obj An object of class mROC
#' @return Returns the aurea under the mROC curve
#' @export
mAUC<-function(mROC_obj)
{
  l<-length(mROC_obj$FPs)
  x<-mROC_obj$FPs[-1]-mROC_obj$FPs[-l]
  a<-sum(x*mROC_obj$TPs[-1])
  b<-sum(x*mROC_obj$TPs[-l])
  return((a+b)/2)
}





#Calculates the asolute surface between the empirical and expected ROCs
calc_mROC_stats_old<-function(p,y, ordered=F, fast=T)
{
  if(!ordered)
  {
    o<-order(p)
    p<-p[o]
    y<-y[o]
  }

  if(fast)
  {
    tmp<-Ccalc_mROC_stats(p,y)
    return(c(A=tmp[1],B=tmp[2]))
  }

  n0<-length(which(y==0))
  n1<-length(which(y==1))

  n<-n0+n1

  tpO<-0
  fpO<-0

  sumP1<-sum(p)
  sumP0<-sum(1-p)
  tpE<-0
  fpE<-0

  j<-n

  delta<-0

  for(i in (n0+n1):1)
  {
    if(y[i]==1)
    {
      tpO<-tpO+1/n1
    }
    else
    {
      fpO<-fpO+1/n0

      while(fpO>fpE && j>0)
      {
        fpE<-fpE+(1-p[j])/sumP0
        tpE<-tpE+p[j]/sumP1
        j<-j-1
        delta<-delta+abs(tpO-tpE)
        #cat(paste(fpO,",",fpE,";"))
      }
      #cat("\n")
    }
  }
  return(delta/n)
}








#Calculates the asolute surface between the empirical and expected ROCs
calc_mROC_stats<-function(p,y, ordered=F, fast=T)
{
  if(!ordered)
  {
    o<-order(p)
    p<-p[o]
    y<-y[o]
  }
  
  if(fast)
  {
    tmp<-Ccalc_mROC_stats(p,y)
    return(c(A=tmp[1],B=tmp[2]))
  }
  
  n0<-length(which(y==0))
  n1<-length(which(y==1))
  n<-n0+n1
  sumP1<-sum(p)
  sumP0<-sum(1-p)
  
  xo<-0
  xe<-0
  yo<-0
  ye<-0
  io<-n
  ie<-n
  
  B<-0
  
  step<-0
  
  #plot(c(0,1),c(0,1))
  
  while(io>0 && ie>0)
  {
    if(xo<xe) #xo is behind and has to make a jump
    {
      if(y[io]==1)
      {
        step<-0
        yo<-yo+1/n1
      }
      else
      {
        step<-1/n0
        B<-B+abs(yo-ye)*min(step,xe-xo)
      }
      xo<-xo+step
      io<-io-1
    }
    else #now xe is behind
    {
      step<-(1-p[ie])/sumP0
      B<-B+abs(yo-ye)*min(step,xo-xe)
      xe<-xe+step
      ye<-ye+p[ie]/sumP1
      ie<-ie-1
    }
    #lines(xe,ye,col='red',type='o')
    #lines(xo,yo,type='o')
  }
  
  return(list(A=abs(mean(p-y)),B=B))
}




aux<-environment()




#' Statistical inference for comparing empirical and expectec ROCs. If CI=TRUE then also returns pointwise CIs
#' @param p vector of probabilities
#' @param y vector of binary response values
#' @param n_sim number of Monte Carlo simulations to calculate p-value
#' @param CI whether confidence interval should be alculated for each point of mROC
#' @return Returns the aurea under the mROC curve
#' @export

mROC_inference<-function(p,y,n_sim=100000,CI=FALSE,aux=FALSE,fast=TRUE,conditional=FALSE)
{
  out<-list()

  n<-length(p)

  if(aux)
  {
    aux$AB<<-matrix(NA, nrow =n_sim, ncol=2)
  }

  o<-order(p)
  p<-p[o]
  y<-y[o]
  
  if(conditional)
  {
    p1<-t.test(p-y)$p.value
    
    n1<-sum(y)
    a<-p
    b<-p
    
    for(i in 1:n)
    {
      a[i]<-p[i]*dpoisbinom(n1-1,p[-i])
      b[i]<-(1-p[i])*dpoisbinom(n1,p[-i])
    }
    
    pik<-a/(a+b)
    
    pikt=UPMEpiktildefrompik(pik)
  
    w=pikt/(1-pikt)
    q=UPMEqfromw(w,n1)
    
    if(fast)
    {
      tmp=Csimulate_null_ds_conditional_crazy(p,n1,n_sim)
      #tmp=Csimulate_null_ds_conditional(q,n_sim)
      stats<-calc_mROC_stats(p,y)
      out$stats<-stats
      
      out$null_stats<-c(distance.se=sqrt(var(tmp[,1]/length(tmp[,1]))),surface.mu=mean(tmp[,2]),surface.se=sqrt(var(tmp[,2]/length(tmp[,2]))))
      
      cdf2<-ecdf(x = tmp[,2])
      
      p2<-1-cdf2(stats[2])
      d<- -2*(log(p1)+log(p2))
      
      out$stat<-c(value=d,df=4)
      
      p3<-1-pchisq(q = d, df = 4)
      
      out$pval<-p3
      out$pvals<-c(distance=p1,surface=p2)
      
    }
    else
    {
      #TODO
    }
  }
  else #If conditional
  {
    if(fast)
    {
      n1<-sum(y)
      
      tmp<-Csimulate_null_mROC_stats_unconditional(p,n_sim)
      
      stats<-calc_mROC_stats(p,y)
      out$stats<-stats
  
      out$null_stats<-c(A.mu=mean(tmp[,1]),A.se=sqrt(var(tmp[,1]/length(tmp[,1]))),B.mu=mean(tmp[,2]),B.se=sqrt(var(tmp[,2]/length(tmp[,2]))))
  
      cdf1<-ecdf(x = tmp[,1])
      cdf2<-ecdf(x = tmp[,2])
  
      p1<-1-cdf1(stats[1])
      p2<-1-cdf2(stats[2])
      d<- -2*(log(p1)+log(p2))
  
      p1s<-cdf1(tmp[,1])
      p2s<-cdf2(tmp[,2])
      ds<--2*(log(p1s)+log(p2s))
      var_ds<-var(ds)
      e_ds<-mean(ds)
      c<-var_ds/2/e_ds
      k<-2*e_ds^2/var_ds
  
      out$stat<-c(value=d/c,df=k)
  
      p3<-1-pchisq(q = d/c, df = k)
  
      out$pval<-p3
      out$pvals<-c(A=p1,B=p2)
  
      if(aux)
      {
        p3s<-1-pchisq(q = ds/c, df = k)
        aux$ds<<-tmp
        aux$pvals<<-cbind(p1s,p2s,p3s)
      }
    }
    else #if (fast)
    {
      tmp<-matrix(NA,n_sim,2)
      sns<-matrix(NA,nrow = n+1, ncol = n_sim)
      for(i in 1:n_sim)
      {
        y<-rbinom(n,size = 1,prob=p)
        if(CI)
        {
          res<-roc(y,p,quiet = T)
          sns[,i]<-coords(res,(0:n)/n,input="specificity",transpose=FALSE)[,"sensitivity",]
        }
        tmp[i,]<-calc_mROC_stats(p,y,ordered = TRUE)
  
        if(aux)
        {
          aux$AB[i,]<<-c(mean(y)-mean(m),tmp[i])
        }
      }
    }
  } #else conditional
  if(CI)
  {
    tmp2<-apply(sns,MARGIN = 1,FUN = ecdf)
    out$low<-unlist(lapply(tmp2,quantile,0.025))
    out$high<-unlist(lapply(tmp2,quantile,0.975))
  }


  return(out)
}







#Main eRoc analysis: draws the ROC and eROC. inference=0: no inference, inference=1: p-value, inference=2: p-value and 95%CI
mROC_analysis<-function(p,y,inference=0, n_sim, fast=TRUE)
{
  if(inference==2 && fast) stop("Confidence intervals are currently only available when fast=FALSE")
  
  roc_data<-pROC::roc(y,p)
  plot(roc_data)

  out<-list()
  out$roc_data<-roc_data

  message("AUC is ",roc_data$auc)
  
  res<-mROC(p)
  out$mROC_data<-res
  
  message("mAUC is ",calc_mAUC(res))

  lines(1-res[,'fp'],res[,'tp'],col='red')

  if(inference)
  {
    inf<-mROC_inference(p=p, y=y, CI=(inference==2), n_sim = n_sim,  fast=fast)
    if(inference==2)
    {
      n<-length(p)
      lines((0:n)/n,inf$low,type='l',col="gray")
      lines((0:n)/n,inf$high,type='l',col="gray")
    }
    out$inference<-inf
    message("Test statistic is ",inf$stat,"\n")
    message("p-value is ",inf$pval,"\n")
  }

  return(out)
}











