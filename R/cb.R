

#the S3 class!
Cb_output.template<-list(Cb=NA,e_b=NA,e_max_b1b2=NA,Gini=NA,AUCi=NA,p=NA,q=NA)
class(Cb_output.template)<-"Cb_output"

#' @export
print.Cb_output<-function(x, ...)
{
  cat("Cb=",x$Cb,"\ne_b=",x$e_b,"\ne_max_b1b2=",x$e_max_b1b2,"\nGini=",x$Gini,"\nAUCi=",x$AUCi,"\nData length:",length(x$q))
  return(invisible(x$Cb))
}

#' @export
plot.Cb_output<-function(x, ...)
{
  plot(x$p,x$q,type='l',xlab="p",ylab="q")
}


#' @export
lines.Cb_output<-function(x, ...)
{
  lines(x$p,x$q,type='l',xlab="p",ylab="q")
}


e_max_b1_b2<-function(B, dumb=FALSE, ordered=FALSE)
{
  out<-0
  
  n<-length(B)
  
  if(dumb)
  {
    N<-n*n
    
    for(i in 1:(n-1))
      for(j in (i+1):n)
      {
        out<-out+2*max(B[i],B[j])
      }
    
    out<-(out+sum(B))
    
    return(out/N)
  }
  else
  {
    if(!ordered)  B<-B[order(B,decreasing = TRUE)]
    return((2*sum(cumsum(B))-sum(B))/(n^2))
  }
}



#' Simple, parametric calculation of Cb.
#' @param B A numberic vector.
#' @return This function returns Cb as 1-E(B1)/E(max(B2,B3)), where B1, B2, and B3 are random draws from the empirical distribution of B.
#' @examples
#' b <- runif(1000)
#' Cb.simple(b)
#' @export
Cb.simple<-function(B)
{
  a<-mean(B)
  
  if(a<0)
  {
    B<--B
    a<--a
  }
  
  o<-order(B,decreasing = T)
  B<-B[o]
  n<-length(B)
  
  b<-e_max_b1_b2(B,ordered = T)
  
  out<-Cb_output.template
  
  out$Cb<-1-a/b
  out$e_b<-a
  out$e_max_b1b2<-b
  out$Gini<-(b-a)/a
  out$AUCi<-b/a/2
  out$p<-(1:n)/n
  out$q<-cumsum(B)
  return(out)
}







#' Cb calculations for a logistic regression model.
#' @param reg_object An object of class 'glm' that contains the resuls of the logit model.
#' @param tx_var A string containing the name of the treatment indicator variable.
#' @param semi_parametric Optional (default=FALSE). If TRUE, the semi-parametric estimator for Cb will be returned.
#' @return This function returns an object of class Cb_output, which includes Cb as a member.
#' @examples
#' data("rct_data")
#' #Creating a binary variable indicating whether an exacerbation happened during the first 6 months.
#' #Because everyone is followed for at least 6 months, there is no censoring.
#' rct_data[,'b_exac']<-rct_data[,'tte']<0.5
#' rct_data[which(is.na(rct_data[,'b_exac'])),'b_exac']<-FALSE
#'
#' reg.logostic<-glm(formula = b_exac ~ 
#'                   tx + sgrq + prev_hosp + prev_ster + fev1, 
#'                   data = rct_data, 
#'                   family = binomial(link="logit"))
#' res.logistic<-Cb.logistic(reg.logostic,tx_var = "tx", semi_parametric = T)
#' print(res.logistic)
#' @export
Cb.logistic<-function(reg_object,tx_var,semi_parametric=FALSE)
{
  if(!inherits(reg_object,"glm")) stop("reg_object should be an object of class glm.")
  if(is.null(tx_var)) stop("Treatment variable label (tx_var) is not speficied.")
  
  tx_values<-unique(reg_object$model[,tx_var])
  if(length(tx_values)!=2) stop("Treatment variable must have two and only two values.")
  
  out<-Cb_output.template
  
  data<-reg_object$model
  
  n<-dim(data)[1]
  
  newdata0<-data
  newdata0[,tx_var]<-tx_values[1]
  y0<-predict.glm(reg_object,newdata = newdata0, type = "response")
  
  newdata1<-data
  newdata1[,tx_var]<-tx_values[2]
  y1<-predict.glm(reg_object,newdata = newdata1, type="response")
  
  if(semi_parametric)
  {
    outcomes<-reg_object$y
    
    B<-y0-y1
    
    if(mean(B)<0)
    {
      B<- -B
      data[,tx_var]<- sum(tx_values)-data[,tx_var]
    }
    
    o<-order(B,runif(n),decreasing=TRUE)
    
    id0<-which(data[,tx_var]==tx_values[1])
    data[,'y0__']<-0
    data[,'t0__']<-0
    data[id0,'t0__']<-1
    data[id0,'y0__']<-outcomes[id0]
    
    id1<-which(data[,tx_var]==tx_values[2])
    data[,'y1__']<-0
    data[,'t1__']<-0
    data[id1,'t1__']<-1
    data[id1,'y1__']<-outcomes[id1]
    
    data<-data[o,]
    B<-B[o]
    
    tmp0<-cumsum(data[,'y0__'])*(1:n)/cumsum(data[,'t0__'])
    tmp1<-cumsum(data[,'y1__'])*(1:n)/cumsum(data[,'t1__'])
    
    n_nan<-sum(is.nan(tmp0))
    if(n_nan>0)
    {
      tmp0[1:n_nan]<-(1:n_nan)*tmp0[n_nan+1]/n_nan
    }
    else
    {
      n_nan<-sum(is.nan(tmp1))
      if(n_nan>0)
      {
        tmp1[1:n_nan]<-(1:n_nan)*tmp1[n_nan+1]/n_nan
      }
    }
    
    data[,'q_sp__']<-tmp0-tmp1
    
    out$e_b<-sum(data[,'y0__'])/sum(data[,'t0__'])-sum(data[,'y1__'])/sum(data[,'t1__'])
    out$e_max_b1b2<-(2*sum(data[,'q_sp__']))/n/n-out$e_b/n
    out$Cb<-1-out$e_b/out$e_max_b1b2
    out$p<-(1:n)/n
    out$q<-data[,'q_sp__']
    
    return(out)
  }
  else
  {
    return(Cb.simple(y0-y1))
  }
}







#' Cb calculations for a Poisson (or negative binomial) regression model.
#' @param reg_object An object of class 'glm' that contains the resuls of the Poisson/NegBin model.
#' @param tx_var A string containing the name of the treatment indicator variable.
#' @param semi_parametric Optional (default=FALSE). If TRUE, the semi-parametric estimator for Cb will be returned.
#' @param time Optional (default=1). The value of time at which Cb is calculated.
#' @return This function returns an object of class Cb_output, which includes Cb as a member.
#' @examples
#' data("rct_data")
#' reg<-glm(formula = n_exac ~ tx + sgrq + prev_hosp + prev_ster + fev1, 
#'          data = rct_data, 
#'          family = poisson(link="log"), 
#'          offset=ln_time)
#' res.Poisson<-Cb.Poisson(reg,tx_var = "tx", semi_parametric = T)
#' res.Poisson
#' @export
Cb.poisson<-function(reg_object,tx_var,semi_parametric=FALSE,time=1)
{
  if(!inherits(reg_object,"glm")) stop("reg_object should be an object of class glm.")
  if(is.null(tx_var)) stop("Treatment variable label (tx_var) is not speficied.")
  
  out<-Cb_output.template
  
  data<-reg_object$data
  
  n<-dim(data)[1]
  
  newdata0<-data
  newdata0[,tx_var]<-0
  y0<-predict.glm(reg_object,newdata = newdata0, type = "link")-reg_object$offset+reg_object$family$linkfun(time)
  y0<-reg_object$family$linkinv(y0)
  newdata1<-data
  newdata1[,tx_var]<-1
  y1<-predict.glm(reg_object,newdata = newdata1, type="link")-reg_object$offset+reg_object$family$linkfun(time)
  y1<-reg_object$family$linkinv(y1)
  
  if(semi_parametric)
  {
    outcome_var<-as.character(reg_object$call$formula[[2]])
    
    B<-y0-y1
    
    if(mean(B)<0)
    {
      B<- -B
      data[,tx_var]<- 1-data[,tx_var]
    }
    
    o<-order(B,runif(n),decreasing=TRUE)
    
    id0<-which(data[,tx_var]==0)
    data[,'y0__']<-0
    data[,'t0__']<-0
    data[id0,'y0__']<-data[id0,outcome_var]
    data[id0,'t0__']<-exp(data[id0,offset_var])
    
    id1<-which(data[,tx_var]==1)
    data[,'y1__']<-0
    data[,'t1__']<-0
    data[id1,'y1__']<-data[id1,outcome_var]
    data[id1,'t1__']<-exp(data[id1,offset_var])
    
    data<-data[o,]
    B<-B[o]
    
    tmp0<-cumsum(data[,'y0__'])*(1:n)/cumsum(data[,'t0__'])
    tmp1<-cumsum(data[,'y1__'])*(1:n)/cumsum(data[,'t1__'])
    
    n_nan<-sum(is.nan(tmp0))
    if(n_nan>0)
    {
      tmp0[1:n_nan]<-(1:n_nan)*tmp0[n_nan+1]/n_nan
    }
    else
    {
      n_nan<-sum(is.nan(tmp1))
      if(n_nan>0)
      {
        tmp1[1:n_nan]<-(1:n_nan)*tmp1[n_nan+1]/n_nan
      }
    }
    
    data[,'q_sp__']<-tmp0-tmp1
    
    out$e_b<-(sum(data[,'y0__'])*n/sum(data[,'t0__'])-sum(data[,'y1__'])*n/sum(data[,'t1__']))/n
    out$e_max_b1b2<-(2*sum(data[,'q_sp__']))/n/n-out$e_b/n
    out$Cb<-1-out$e_b/out$e_max_b1b2
    out$p<-(1:n)/n
    out$q<-data[,'q_sp__']
    
    return(out)
    
  }
  else
  {
    return(Cb.simple(y0-y1))
  }
}










#' Cb calculations for a Cox proportional hazard model.
#' @param reg_object An object of class 'coxph' that contains the model. IMPORTANT: the coxph function call for fitting the model must have model=TRUE as an input argument such that the underlying data becomes part of the returned object (this is the case by default for glm)
#' @param tx_var A string containing the name of the treatment indicator variable.
#' @param semi_parametric Optional (default=FALSE). If TRUE, the semi-parametric estimator for Cb will be returned.
#' @param time Optional (default=1). The value of time at which Cb is calculated.
#' @return This function returns an object of class Cb_output, which includes Cb as a member.
#' @examples
#' library(survival)
#' data("rct_data")
#' # Create an event indicator and update the tte (time-to-event) variable to be 
#' # equal to follow-up time for censored individuals.
#' event<-(!is.na(rct_data[,'tte']))*1
#' ids<-which(event==0)
#' rct_data[ids,'tte']<-rct_data[ids,'time']
#' rct_data['event']<-event
#'
#' reg.coxph<-coxph(Surv(time=tte,event=event) ~ 
#' tx + tx:female + tx:age + sgrq + prev_hosp + prev_ster + fev1, 
#' data=rct_data, model=TRUE)
#' res.coxph<-Cb.coxph(reg.coxph,tx_var = "tx",semi_parametric = T)
#' @export
Cb.cox<-function(reg_object,tx_var,semi_parametric=FALSE, time=1)
{
  if(!inherits(reg_object,"coxph")) stop("reg_object should be an object of class coxph.")
  if(is.null(tx_var)) stop("Treatment variable label (tx_var) is not speficied.")
  if(is.null(reg_object$model)) stop("No model data available in the regression object. Run coxph with model=TRUE argument.")
  
  out<-Cb_output.template
  
  data<-reg_object$model
  time_var<-as.character(reg_object$formula[[2]][[2]])
  
  n<-dim(data)[1]
  
  newdata0<-data
  newdata0[,tx_var]<-0
  newdata0[,time_var]<-time
  y0<-predict(reg_object,newdata = newdata0, type = "expected")
  
  newdata1<-newdata0
  newdata1[,tx_var]<-1
  y1<-predict(reg_object,newdata = newdata1, type="expected")
  
  if(semi_parametric)
  {
    x<-as.matrix(reg_object$y)
    times<-x[,1]
    events<-x[,2]
    
    B<-y0-y1
    
    if(mean(B)<0)
    {
      B<- -B
      data[,tx_var]<- 1-data[,tx_var]
    }
    
    o<-order(B,runif(n),decreasing=TRUE)
    
    id0<-which(data[,tx_var]==0)
    data[,'y0__']<-0
    data[,'t0__']<-0
    data[id0,'y0__']<-events[id0]
    data[id0,'t0__']<-times[id0]
    
    id1<-which(data[,tx_var]==1)
    data[,'y1__']<-0
    data[,'t1__']<-0
    data[id1,'y1__']<-events[id1]
    data[id1,'t1__']<-times[id1]
    
    data<-data[o,]
    B<-B[o]
    
    tmp0<-cumsum(data[,'y0__'])*(1:n)/cumsum(data[,'t0__'])
    tmp1<-cumsum(data[,'y1__'])*(1:n)/cumsum(data[,'t1__'])
    
    n_nan<-sum(is.nan(tmp0))
    if(n_nan>0)
    {
      tmp0[1:n_nan]<-(1:n_nan)*tmp0[n_nan+1]/n_nan
    }
    else
    {
      n_nan<-sum(is.nan(tmp1))
      if(n_nan>0)
      {
        tmp1[1:n_nan]<-(1:n_nan)*tmp1[n_nan+1]/n_nan
      }
    }
    
    data[,'q_sp__']<-tmp0-tmp1
    
    out$e_b<-(sum(data[,'y0__'])*n/sum(data[,'t0__'])-sum(data[,'y1__'])*n/sum(data[,'t1__']))/n
    out$e_max_b1b2<-(2*sum(data[,'q_sp__']))/n/n-out$e_b/n
    out$Cb<-1-out$e_b/out$e_max_b1b2
    out$p<-(1:n)/n
    out$q<-data[,'q_sp__']
    
    return(out)
  }
  else
  {
    return(Cb.simple(y0-y1))
  }
}