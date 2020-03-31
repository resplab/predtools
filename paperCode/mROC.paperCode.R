#Code for the main paper
#Authors: Mohsen Sadatsafavi
#Last update: 2020.02.17




library(mROC)
library(pROC)
library(sqldf)


project_folder<-"M:/Projects/2018/Project.GhostROC/Code/"
setwd(project_folder)


settings<-list(
  n_sim_inference_fast=100000,
  n_sim_inference_slow=10000
  )



aux<-list() #For global objects that functions might return;




recalibrate<-function(p,y)
{
  x<-log(p/(1-p))
  reg<-glm(y~offset(x),family=binomial(link="logit"))
  x<-x+coefficients(reg)[1]
  return(1/(1+exp(-x)))
}




calibration_plot<-function(p,y,n_bins=10,type_prob=TRUE)
{
  n<-length(p)
  o<-order(p)
  p<-p[o]
  y<-y[o]
  
  out_x<-rep(NA,n_bins)
  out_y<-rep(NA,n_bins)
  out_y_l<-rep(NA,n_bins)
  out_y_h<-rep(NA,n_bins)
  
  for(i in 1:n_bins)
  {
    i_l<-floor((i-1)*n/n_bins)+1
    i_h<-floor(i*n/n_bins)
    mp<-mean(p[i_l:i_h])
    my<-mean(y[i_l:i_h])
    sey<-sqrt(var(y[i_l:i_h])/(i_h-i_l+1))
    out_x[i]<-mp
    out_y[i]<-my
    out_y_l[i]<-my-1.96*sey
    out_y_h[i]<-my+1.96*sey
  }
  
  if(type_prob)
  {
    max_x<-1
    max_y<-1
    xlab<-"Predicted risk"
    ylab<-"Observed risk"
  }
  else
  {
    max_x<-max(out_x)
    max_y<-max(out_y_h)
    xlab<-"Predicted rate"
    ylab<-"Observed rate"
  }
  
  plot(c(0,max_x),c(0,max_y),type='l',xlim=c(0,max_x),ylim=c(0,max_y),xlab=xlab,ylab=ylab)
  
  for(i in 1:n_bins)
  {
    lines(mp,my,type="p")
  }
  
  for(i in 1:n_bins)
  {
    lines(out_x[i],out_y[i],type="p")
    lines(c(out_x[i],out_x[i]),c(out_y_l[i],out_y_h[i]),type="l")
  }
  
  return(list(out_x=out_x,out_y=cbind(out_y,out_y_l,out_y_h)))
}














############Section 1 code: stylized example######################


stylized_example<-function()
{
  t<-seq(0,1,length.out = 100)
  y<-rep(NA,length(t))

  #ROC
  plot(t,2*sqrt(t)-t,type='l',xlab = "False Positive", ylab="True Positive", lwd=2)
  lines(c(0,1),c(0,1),type="l",col="grey")

  #If pi* = sqrt(pi):
  for(i in 1:length(y))
  {
    y[i]<-uniroot(function(x) 3*x^2-2*x^3-(1-t[i]),interval = c(0,1))$root
    y[i]<-1-y[i]^3
  }
  lines(t,y,type='l',col="red", lwd=2)

  #If pi* = pi^2:
  for(i in 1:length(y))
  {
    y[i]<-uniroot(function(x) (3*sqrt(x)-(x)^(3/2))/2-(1-t[i]),interval = c(0,1))$root
    y[i]<-1-y[i]^(3/2)
  }
  lines(t,y,type='l',col="blue", lwd=2)

  plot(t,t,type="l",xlab="Predicted risk",ylab="Actual risk",lwd=2)
  lines(t,t^2,col="red",lwd=2)
  lines(t,sqrt(t),col="blue",lwd=2)
}











#################################Section 2: stylized simulation#################

stylized_simulation<-function(n=10000)
{
  
  x<-rnorm(n,0,1)
  b0<-0
  b1<-1


  #Internal model
  p<-1/(1+exp(-(b0+b1*x)))
  y<-rbinom(n,size=1,p=p)
  r<-pROC::roc(y,p)
  plot(1-r$specificities,r$sensitivities,xlim=c(0,1),ylim=c(0,1),type='l',xlab="False Positive",ylab="True Positive", main="Development and validation sample exchangable")
  message("Development c=",r$auc)
  
  
  #Validaiton - external exchangable
  xx<-rnorm(n,0,1)
  pp<-1/(1+exp(-(b0+b1*xx)))
  yy<-rbinom(n,size=1,p=pp)
  zz<-mROC(pp)
  rr<-pROC::roc(yy,pp)
  message("Scenario 1 c=",rr$auc)
  plot(1-rr$specificities,rr$sensitivities,xlim=c(0,1),ylim=c(0,1),type='l',xlab="False Positive",ylab="True Positive", main="Development and validation sample exchangable")
  lines(c(0,1),c(0,1),type='l',col='gray')
  lines(1-r$specificities,r$sensitivities,col="blue")
  lines(zz$FPs,zz$TPs,type='l',col="red")
  calibration_plot(pp,yy)
  
  

  #Validation - Underdispersed predictor, correct model
  xx<-rnorm(n,0,0.5)
  pp<-1/(1+exp(-(b0+b1*xx)))
  yy<-rbinom(n,size=1,p=pp)
  zz<-mROC(pp)
  rr<-pROC::roc(yy,pp)
  message("Scenario 2 c=",rr$auc)
  plot(1-rr$specificities,rr$sensitivities,xlim=c(0,1),ylim=c(0,1),type='l',xlab="False Positive",ylab="True Positive", main="Underdispersed predictor")
  lines(c(0,1),c(0,1),type='l',col='gray')
  lines(1-r$specificities,r$sensitivities,col="blue")
  lines(zz$FPs,zz$TPs,type='l',col="red")
  calibration_plot(pp,yy)



  #Validation - incorrect (optimistic) model
  xx<-rnorm(n,0,1)
  pp<-1/(1+exp(-(b0+b1*xx)))
  pp_real<-1/(1+exp(-(b0+b1/2*xx)))
  yy<-rbinom(n,size=1,p=pp_real)
  zz<-mROC(pp)
  rr<-pROC::roc(yy,pp)
  message("Scenario 3 c=",rr$auc)
  plot(1-rr$specificities,rr$sensitivities,xlim=c(0,1),ylim=c(0,1),type='l',xlab="False Positive",ylab="True Positive", main="Optimistic model")
  lines(c(0,1),c(0,1),type='l',col='gray')
  lines(1-r$specificities,r$sensitivities,col="blue")
  lines(zz$FPs,zz$TPs,type='l',col="red")
  calibration_plot(pp,yy)




  #Validation: Underdispersed predictor + optimistic model
  xx<-rnorm(n,0,0.5)
  pp<-1/(1+exp(-(b0+b1*xx)))
  pp_real<-1/(1+exp(-(b0+b1/2*xx)))
  yy<-rbinom(n,size=1,p=pp_real)
  zz<-mROC(pp)
  rr<-pROC::roc(yy,pp)
  message("Scenario 4 c=",rr$auc)
  plot(1-rr$specificities,rr$sensitivities,xlim=c(0,1),ylim=c(0,1),type='l',xlab="False Positive",ylab="True Positive", main="Underdispersed predictor & optimistic model")
  lines(c(0,1),c(0,1),type='l',col='gray')
  lines(1-r$specificities,r$sensitivities,col="blue")
  lines(zz$FPs,zz$TPs,type='l',col="red")
  calibration_plot(pp,yy)
}






















#################################Section 3: Detailed simulation######################



#X_dist:mean and SD of the distirbution of the simple predictor. If NULL, then directly samples pi from standard uniform. 
detailed_sim<-function(sample_sizes=c(100,250,1000), X_dist=c(0,1), b0s=2*c(-0.25,-0.125,0,0.125,0.25),b1s=c(0.5,0.75,1,1.5,2),n_sim=1000, draw_plots="", GRuse=FALSE)
{
  pi<-runif(100)
  y<-rbinom(100,1,pi)
  template<-unlist(mROC_inference(y=y,p=pi,CI=FALSE, n_sim = 100,fast=TRUE))
  
  out<-as.data.frame(matrix(NA, nrow =n_sim*length(sample_sizes)*length(b0s)*length(b1s),ncol = 4+length(template)+1))
  colnames(out)<-c("i_sim","sample_size", "b0", "b1", names(template),"pval.LRT")
  
  if(draw_plots!="")
  {
    par(mar=c(1,1,1,1))
    par(mfrow=c(length(b0s),length(b1s)))
  }
  
  index<-1
  for(i in 1:n_sim)
  {
    cat(".")
    for(j in 1:length(sample_sizes))
    {
      ss<-sample_sizes[j]
      if(is.null(X_dist)) pi<-runif(ss) else  pi<-1/(1+exp(-rnorm(ss,X_dist[1],X_dist[2])))
      for(k0 in 1:length(b0s))  
      {
        b0<-b0s[k0]
        for(k1 in 1:length(b1s))
        {
          #message(paste(ss,b0,b1,sep=","))
          b1<-b1s[k1]
          pi_star<-1/(1+exp(-(b0+b1*log(pi/(1-pi)))))
        
          repeat
          {
            y=rbinom(ss,size = 1,prob = pi)
            if(sum(y)>=3 && sum(1-y)>=3) 
              break 
            else
            {
                message("OUCH - bad sample, too few outcomes, repeat!")
            }
          }
          
          if(draw_plots!="")
            if(i==1 && j==length(sample_sizes))
            {
              if(draw_plots=="pp")
              {
                o<-order(pi_star)
                pi_<-pi[o]
                pi_star_<-pi_star[o]
                plot(c(0,pi_star_,1),c(0,pi_,1),xlab=(paste(ss,b0,b1,sep=",")), ann=FALSE, xaxt='n', yaxt='n', type='l', xlim=c(0,1), ylim=c(0,1), col="blue", lwd=2)
                lines(c(0,1),c(0,1),col="gray",type='l')
                text(0.50,0.075, sprintf("E(\U03C0*)=%s",format(mean(pi_star),digits =  3)))
                #text(0.75,0.5, sprintf("E(pi)=%s",format(mean(pi),digits = 3)))
                title(sprintf(paste0("\U03B2","0=%s,\U03B2","1=%s"),format(b0,3),format(b1,3)))
              }
              if(draw_plots=="roc")
              {
                tmp<-pROC::roc(y,pi_star)
                plot(1-tmp$specificities,tmp$sensitivities, type='l', ann=FALSE, xaxt='n', yaxt='n')
                AUC=tmp$auc*1
                tmp<-mROC(pi_star)
                mAUC=mAUC(tmp)
                B=calc_mROC_stats(y,pi_star)['B']
                #text(0.5,0.5, sprintf("AUC:%s",format(AUC, digits = 2)),cex = 1) AUC is constant! 
                text(0.5,0.3, sprintf("mAUC:%s",format(mAUC, digits = 2)),cex=1)
                text(0.5,0.1, sprintf("B:%s",format(B, digits = 2)),cex=1)
                lines(tmp$FPs,tmp$TPs,col="red")
                lines(c(0,1),c(0,1),col="grey")
                #title(sprintf("b0=%s,b1=%s",format(b0,3),format(b1,3)))
              }
            }
            
          #message(paste(ss,b0,b1,sep=","))
          
          tmp<-unlist(mROC_inference(y=y,p=pi_star,CI=FALSE, n_sim = settings$n_sim_inference_fast, fast=TRUE))
          out[index,1]<-i
          out[index,2]<-ss
          out[index,3]<-b0
          out[index,4]<-b1
          out[index,5:(5+length(template)-1)]<-tmp
          logit.pi_star<-log(pi_star/(1-pi_star))
          f.0<-glm(y~ -1+offset(logit.pi_star),family="binomial")
          #f.a<-glm(y~offset(logit.pi_star),family="binomial")
          f.ab<-glm(y~logit.pi_star,family="binomial")
          #message(coefficients(f.ab))
          p.val.ab<-1-pchisq(f.0$deviance-f.ab$deviance,2)
          out[index,"pval.LRT"]<-p.val.ab
          index<-index+1
        }
      }
    }
    
    if((i%%10)==0 && GRuse)
    {
      GRpush(out,overWrite = T)
      cat("Pushed at i=",i," \n")
    }
  }
  
  aux$out<<-out
  
  return(out)
}





#X_dist:mean and SD of the distirbution of the simple predictor. If NULL, then directly samples pi from standard uniform. 
detailed_sim_power<-function(sample_sizes=c(100,250,1000), X_dist=c(0,1), b0s=c(0,0.25,0.5), b1s=c(0.5,0.75,1,1.5,2), b2s=NULL, n_sim=1000, draw_plots="", GRuse=FALSE)
{
  pi<-runif(100)
  y<-rbinom(100,1,pi)
  template<-unlist(mROC_inference(y=y,p=pi,CI=FALSE, n_sim = 100,fast=TRUE))
  
  out<-as.data.frame(matrix(NA, nrow =n_sim*length(sample_sizes)*length(b0s)*length(b1s),ncol = 5+length(template)+1))
  colnames(out)<-c("i_sim","sample_size", "b0", "b1", "b2", names(template), "pval.LRT")
  
  if(draw_plots!="")
  {
    par(mar=c(1,1,1,1))
    if(is.null(b2s)) 
      par(mfrow=c(length(b0s),length(b1s)))
    else
      par(mfrow=c(length(b1s),length(b2s)))
  }

  index<-1
  for(i in 1:n_sim)
  {
    cat(".")
    for(j in 1:length(sample_sizes))
    {
      ss<-sample_sizes[j]
      if(is.null(X_dist)) pi<-runif(ss) else  pi<-1/(1+exp(-rnorm(ss,X_dist[1],X_dist[2])))
      for(k0 in 1:length(b0s))  
      {
        b0<-b0s[k0]
        for(k1 in 1:length(b1s))
        {
          b1<-b1s[k1]
          if(is.null(b2s) || b2s<0) b2s<- -1/b1   #If b2s is null, it is set as reciprocal of b1
          for(k2 in 1:length(b2s))
          {
            b2<-abs(b2s[k2])
            
            x<-log(pi/(1-pi))
            xx<-b0+sign(x)*b1*(abs(x)^b2)
            pi_star<-1/(1+exp(-xx))
            
            repeat
            {
              y=rbinom(ss,size = 1,prob = pi)
              if(sum(y)>=3 && sum(1-y)>=3) 
                break 
              else
              {
                message("OUCH - bad sample, too few outcomes, repeat!")
              }
            }
            
            if(draw_plots!="")
              if(i==1 && j==length(sample_sizes))
              {
                aux$y <<- y
                aux$pi_star <<- pi_star
                aux$pi <<- pi
                
                if(draw_plots=="pp")
                {
                  q<-(0:100)/100
                  q_x<-log(q/(1-q))
                  
                  p_x<-(abs(q_x-b0)/b1)^(1/b2)
                  p_x[which(q_x-b0<0)]<- -abs(p_x[which(q_x-b0<0)])
                  p<-1/(1+exp(-p_x))
                  
                  plot(q,p ,xlab=(paste(ss,b0,b1,b2,sep=",")), ann=FALSE, xaxt='n', yaxt='n', type='l', xlim=c(0,1), ylim=c(0,1), col="blue", lwd=2)
                  lines(c(0,1),c(0,1),col="gray",type='l')
                  text(0.50,0.075, sprintf("E(\U03C0*)=%s",format(mean(pi_star),digits =  3)))
                  #text(0.75,0.5, sprintf("E(pi)=%s",format(mean(pi),digits = 3)))
                  title(sprintf(paste0("a=%s,b=%s"),format(b0,3),format(b1,3)))
                }
                if(draw_plots=="roc")
                {
                  plot(c(0,1),c(0,1),col="grey")
                  #for(i in 1:100)
                  {
                  #  yy<-rbinom(length(pi_star),1,pi_star)
                  #  tmp<-roc(yy,pi_star)
                  #  lines(1-tmp$specificities,tmp$sensitivities,col="grey")
                  }
                  tmp<-pROC::roc(y,pi_star)
                  lines(1-tmp$specificities,tmp$sensitivities, type='l', ann=FALSE, xaxt='n', yaxt='n')
                  AUC=tmp$auc*1
                  tmp<-mROC(pi_star)
                  mAUC=mAUC(tmp)
                  B=calc_mROC_stats(y,pi_star)['B']
                  #text(0.5,0.5, sprintf("AUC:%s",format(AUC, digits = 2)),cex = 1) AUC is constant! 
                  text(0.5,0.3, sprintf("mAUC:%s",format(mAUC, digits = 2)),cex=1)
                  text(0.5,0.1, sprintf("B:%s",format(B, digits = 2)),cex=1)
                  sf<-stepfun(tmp$FPs,c(0,tmp$TPs))
                  lines(sf,col="red")
                  #lines(tmp$FPs,tmp$TPs,col="red")
                  #title(sprintf("b0=%s,b1=%s",format(b0,3),format(b1,3)))
                }
                if(draw_plots=="logit")
                {
                  plot(log(pi_star/(1-pi_star)),log(pi/(1-pi)),xlim=c(-3,3),ylim=c(-3,3))
                  lines(c(-3,3),c(-3,3))
                  text(0.50,0.075, sprintf("E(\U03C0*)=%s",format(mean(pi_star),digits =  3)))
                  #text(0.75,0.5, sprintf("E(pi)=%s",format(mean(pi),digits = 3)))
                  title(sprintf("a=%s,b=%s",format(b0,3),format(b1,3)))
                }
              }
            
            #message(paste(ss,b0,b1,sep=","))
            
            tmp<-unlist(mROC_inference(y=y,p=pi_star,CI=FALSE, n_sim = settings$n_sim_inference_fast, fast=TRUE))
            out[index,1]<-i
            out[index,2]<-ss
            out[index,3]<-b0
            out[index,4]<-b1
            out[index,5]<-b2
            out[index,6:(6+length(template)-1)]<-tmp
            
            logit.pi_star<-log(pi_star/(1-pi_star))
            f.0<-glm(y~ -1+offset(logit.pi_star),family="binomial")
            #f.a<-glm(y~offset(logit.pi_star),family="binomial")
            f.ab<-glm(y~logit.pi_star,family="binomial")
            #message(coefficients(f.ab))
            p.val.ab<-1-pchisq(f.0$deviance-f.ab$deviance,2)
            out[index,"pval.LRT"]<-p.val.ab
            index<-index+1
          }
        }
      }
    }
    
    if((i%%10)==0 && GRuse)
    {
      GRpush(out,overWrite = T)
      cat("Pushed at i=",i," \n")
    }
  }
  
  aux$out<<-out
  
  return(out)
}







detailed_sim_quadratic<-function(sample_sizes=c(100,250,1000), X_dist=c(0,1), b0s=c(-1,-0.5,0,0.5,1),b1s=c(-1.41,0.7,0,0.7,1.41),n_sim=1000, draw_plots="", GRuse=FALSE)
{
  b1<-1 
  pi<-runif(100)
  y<-rbinom(100,1,pi)
  template<-unlist(mROC_inference(y=y,p=pi,CI=FALSE, n_sim = 100,fast=TRUE))
  
  out<-as.data.frame(matrix(NA, nrow =n_sim*length(sample_sizes)*length(b0s)*length(b1s),ncol = 4+length(template)+1))
  colnames(out)<-c("i_sim","sample_size", "b0", "b1", names(template),"pval.LRT")
  
  index<-1
  for(i in 1:n_sim)
  {
    cat(".")
    for(j in 1:length(sample_sizes))
    {
      ss<-sample_sizes[j]
      if(is.null(X_dist)) pi<-runif(ss) else  pi<-1/(1+exp(-rnorm(ss,X_dist[1],X_dist[2])))
      for(k0 in 1:length(b0s))  
      {
        b0<-b0s[k0]
        for(k1 in 1:length(b1s))
        {
          #message(paste(ss,b0,b1,sep=","))
          b2<-b1s[k1]
          logit.pi<-log(pi/(1-pi))
          logit.pi_star<-b0+b1*logit.pi+b2*logit.pi*logit.pi
          pi_star<-1/(1+exp(-logit.pi_star))
          
          repeat
          {
            y=rbinom(ss,size = 1,prob = pi)
            if(sum(y)>=3 && sum(1-y)>=3) 
              break 
            else
            {
              message("OUCH - bad sample, too few outcomes, repeat!")
            }
          }
          
          if(draw_plots!="")
            if(i==1 && j==length(sample_sizes))
            {
              #par(mar=c(1,1,1,1))
              #par(mfrow=c(length(b0s),length(b1s)))
              
              aux$y <<- y
              aux$pi_star <<- pi_star
              aux$pi <<- pi
              
              if(draw_plots=="pp")
              {
                q<-(0:100)/100
                q_x<-log(q/(1-q))
                
                if(b2==0)
                {
                  root<-(q_x-b0)/b1
                  q_p<-1/(1+exp(-root))
                }
                else
                {
                  delta<- b1^2-4*b2*(b0-q_x)
                  good<-which(delta>0)
                  q_p<-q
                  q_p[-good]<-NA
                  root1 <- (-b1-sqrt(delta[good]))/2/b2
                  root2 <- (-b1+sqrt(delta[good]))/2/b2
                  q_p[good] <- (dnorm(root1)*1/(1+exp(-root1))+dnorm(root2)*1/(1+exp(-root2)))/(dnorm(root1)+dnorm(root2))
                }
                
                plot(q,q_p,xlab=(paste(ss,b0,b1,sep=",")), ann=FALSE, xaxt='n', yaxt='n', type='l', xlim=c(0,1), ylim=c(0,1), col="blue", lwd=2)
                lines(c(0,1),c(0,1),col="gray",type='l')
                text(0.50,0.075, sprintf("E(\U03C0*)=%s",format(mean(pi_star),digits =  3)))
                #text(0.75,0.5, sprintf("E(pi)=%s",format(mean(pi),digits = 3)))
                title(sprintf(paste0("\U03B2","0=%s,\U03B2","1=%s"),format(b0,3),format(b2,3)))
              }
              if(draw_plots=="roc")
              {
                plot(c(0,1),c(0,1),col="grey")
                for(i in 1:100)
                {
                  yy<-rbinom(length(pi_star),1,pi_star)
                  tmp<-roc(yy,pi_star)
                  lines(1-tmp$specificities,tmp$sensitivities,col="grey")
                }
                tmp<-pROC::roc(y,pi_star)
                lines(1-tmp$specificities,tmp$sensitivities, type='l', ann=FALSE, xaxt='n', yaxt='n')
                AUC=tmp$auc*1
                tmp<-mROC(pi_star)
                mAUC=mAUC(tmp)
                B=calc_mROC_stats(y,pi_star)['B']
                #text(0.5,0.5, sprintf("AUC:%s",format(AUC, digits = 2)),cex = 1) AUC is constant! 
                text(0.5,0.3, sprintf("mAUC:%s",format(mAUC, digits = 2)),cex=1)
                text(0.5,0.1, sprintf("B:%s",format(B, digits = 2)),cex=1)
                sf<-stepfun(tmp$FPs,c(0,tmp$TPs))
                lines(sf,col="red")
                #title(sprintf("b0=%s,b1=%s",format(b0,3),format(b1,3)))
              }
              if(draw_plots=="logit")
              {
                plot(logit.pi_star,logit.pi,xlim=c(-3,3),ylim=c(-3,3))
                lines(c(-3,3),c(-3,3))
                text(0.50,0.075, sprintf("E(\U03C0*)=%s",format(mean(pi_star),digits =  3)))
                #text(0.75,0.5, sprintf("E(pi)=%s",format(mean(pi),digits = 3)))
                title(sprintf(paste0("\U03B2","0=%s,\U03B2","1=%s"),format(b0,3),format(b1,3)))
              }
              if(draw_plots=="logitE")
              {
                q_x<-(-300:300)/100
                if(b2==0)
                {
                  root<-(q_x-b0)/b1
                  q_y<-root
                }
                else
                {
                  delta<- b1^2-4*b2*(b0-q_x)
                  good<-which(delta>0)
                  q_y<-q_x
                  q_y[-good]<-NA
                  root1 <- (-b1-sqrt(delta[good]))/2/b2
                  root2 <- (-b1+sqrt(delta[good]))/2/b2
                  q_y[good] <- (dnorm(root1)*root1+dnorm(root2)*root2)/(dnorm(root1)+dnorm(root2))
                }
                
                plot(q_x,q_y,xlab=(paste(ss,b0,b1,sep=",")), ann=FALSE, xaxt='n', yaxt='n', type='l', xlim=c(-3,3), ylim=c(-3,3), col="green", lwd=2)
                lines(c(-3,3),c(-3,3),col="gray",type='l')
                text(0.50,0.075, sprintf("E(\U03C0*)=%s",format(mean(pi_star),digits =  3)))
                #text(0.75,0.5, sprintf("E(pi)=%s",format(mean(pi),digits = 3)))
                title(sprintf(paste0("\U03B2","0=%s,\U03B2","1=%s"),format(b0,3),format(b2,3)))
              }
            }
          
          #message(paste(ss,b0,b1,sep=","))
          
          tmp<-unlist(mROC_inference(y=y,p=pi_star,CI=FALSE, n_sim = settings$n_sim_inference_fast, fast=TRUE))
          out[index,1]<-i
          out[index,2]<-ss
          out[index,3]<-b0
          out[index,4]<-b1
          out[index,5:(5+length(template)-1)]<-tmp
          logit.pi_star<-log(pi_star/(1-pi_star))
          f.0<-glm(y~ -1+offset(logit.pi_star),family="binomial")
          #f.a<-glm(y~offset(logit.pi_star),family="binomial")
          f.ab<-glm(y~logit.pi_star,family="binomial")
          message(coefficients(f.ab))
          p.val.ab<-1-pchisq(f.0$deviance-f.ab$deviance,2)
          out[index,"pval.LRT"]<-p.val.ab
          index<-index+1
        }
      }
    }
    
    if((i%%10)==0 && GRuse)
    {
      GRpush(out,overWrite = T)
      cat("Pushed at i=",i," \n")
    }
  }
  
  aux$out<<-out
  
  return(out)
}









process_detailed_sim_results<-function(detailed=F,dec_points=3)
{
  internal_formatter<-function(data)
  {
    beautify<-function(value)
    {
      return(format(round(value,dec_points),digits = dec_points,nsmall=dec_points))
    }
    
    out<-"D:"
    if(detailed)
    {
      out<-paste0(out, beautify(mean(data[,'stats.distance'])))
      out<-paste0(out,"(",beautify(sd(data[,'stats.distance'])),")")
    }
    out<-paste0(out,"",beautify(sum(data[,'pvals.distance']<0.05)/dim(data)[1]),"")
    out<-paste0(out,";A:")
                
    if(detailed)
    {
      out<-paste0(out, beautify(mean(data[,'stats.surface'])))
      out<-paste0(out,"(",beautify(sd(data[,'stats.surface'])),")")
    }
    out<-paste0(out,"",beautify(sum(data[,'pvals.surface']<0.05)/dim(data)[1]),"")
    out<-paste0(out,";C:")
    out<-paste0(out,beautify(sum(data[,'pval']<0.05)/dim(data)[1]),"")
    
    return(out)
  }
  
  out<-GRcomp:::GRformatOutput(aux$out, rGroupVars = c("b0") , cGroupVars = c("b1"),func = internal_formatter)
  write.table(out,"clipboard")
  return(out)
}






process_detailed_sim_results_graph<-function(data=aux$out, detailed=F, dec_points=3, dim1="b0", dim2="b1")
{
  n_d1<-dim(sqldf(paste0("SELECT DISTINCT ",dim1," FROM data")))[1]
  n_d2<-dim(sqldf(paste0("SELECT DISTINCT ",dim2," FROM data")))[1]
  
  sample_sizes<-sqldf("SELECT DISTINCT sample_size FROM data")
  
  par(mfrow=c(n_d1,n_d2))
  par(mar=0.25*c(1,1,1,1))
  
  out<-GRcomp:::GRformatOutput(data, rGroupVars = c(dim1) , cGroupVars = c(dim2),func = internal_formatter)
  write.table(out,"clipboard")
  return(out)
}



internal_formatter<-function(data)
{
  beautify<-function(value)
  {
    return(format(round(value,dec_points),digits = dec_points,nsmall=dec_points))
  }
  
  #browser()
  y<-sqldf("SELECT sample_size, AVG([pvals.A]<0.05), AVG([pvals.B]<0.05), AVG([pval]<0.05), AVG([pval.LRT]<0.05) FROM data GROUP BY sample_size")
  #y[,1]<-0
  #y<-rbind(y,y/2,y)
  sample_sizes<-y[,1]
  z<-as.vector(t(cbind(y[,-1],0)))
  bp<-barplot(z,xaxt='n', yaxt='n', space=0, ylim=c(-0.25,1.5),col=c("white","pink","orange","purple","black"))
  text(x=0.4+c(0:15)*1,y=z+0.25,ifelse(z==0,"",format(z,3,3)),cex=1, srt=90)
  text(x=c(1.5,5.5,10),y=-0.1,paste(t(sample_sizes)))
  return("")
}





















##########################################Section 5: new case study#####################################################

library("haven")


new_case_study<-function(covars=c("gender","age10","oxygen","hosp1yr","sgrq10","fev1","nowsmk","LABA","LAMA"),only_severe=FALSE, do_recalibrate=FALSE)
{
  results<-list()

  trials_data<<-read_sas(data_file = "M:\\DATA\\2018\\MACRO+STATCOPE+OPTIMAL\\exacevents.sas7bdat")
  trials_data[which(trials_data[,'trial']=="OPTIMAL"),'hosp1yr']<<-1

  #eclipse_data<-readRDS("validatedECLIPSE.RData")

  dev_data<-trials_data[which(trials_data[,'trial']=="MACRO"),]
  dev_data<-as.data.frame(dev_data)
  dev_data<-dev_data[which(dev_data[,'period']==1),]
  results$dev_n_initial<-dim(dev_data)[1]
  missing_data<-which(is.na(rowSums(dev_data[,covars])))
  message("N removed due to missing data from development set:",length(missing_data))
  dev_data<-dev_data[-missing_data,]
  short_fu<-which(dev_data[,'event']==0 & dev_data[,'tte0']<0.5)
  message("N removed due to censoring from development set:",length(short_fu),"(",length(short_fu)/dim(dev_data)[1],")")
  results$dev_n_shortfu<-length(short_fu)
  dev_data<-dev_data[-short_fu,]
  dev_data['event_bin']<-(dev_data['event']>only_severe*1)*1
  results$dev_n_missing<-length(missing_data)
  results$dev_n<-dim(dev_data)[1]
  results$dev_n_event<-sum(dev_data['event_bin'])

  val_data<-trials_data[which(trials_data[,'trial']=="STATCOPE"),]
  val_data<-as.data.frame(val_data)
  val_data<-val_data[which(val_data[,'period']==1),]
  results$val_n_initial<-dim(val_data)[1]
  missing_data<-which(is.na(rowSums(val_data[,covars])))
  message("N removed due to missing data from validation set:",length(missing_data))
  val_data<-val_data[-missing_data,]
  short_fu<-which(val_data[,'event']==0 & val_data[,'tte0']<0.5)
  message("N removed due to censoring from validation set:",length(short_fu),"(",length(short_fu)/dim(val_data)[1],")")
  results$val_n_shortfu<-length(short_fu)
  val_data<-val_data[-short_fu,]
  val_data['event_bin']<-(val_data['event']>only_severe*1)*1
  results$val_n_missing<-length(missing_data)
  results$val_n<-dim(val_data)[1]
  results$val_n_event<-sum(val_data['event_bin'])
  
  #val_data<-eclipse_data; val_data[,'hosp1yr']<-val_data[,'obsExac_yr1'];val_data[,'fev1']<-val_data[,'FEV1']; val_data[,'bmi10']<-val_data[,'BMI10']
  #missing_data<-which(is.na(rowSums(val_data[,covars])))
  #val_data<-val_data[-missing_data,]
  #short_fu<-which(val_data[,'event']==0 & val_data[,'tte0']<0.5)
  #val_data<-val_data[-short_fu,]
  
  formula<-as.formula(paste0("event_bin~",paste(covars,collapse="+"),"+randomized_azithromycin"))
  
  reg<-glm(data=dev_data,formula = formula, family=binomial(link="logit"))
  
  results$model<-reg$coefficients
  
  message(paste(round(results$model,3),names(results$model),sep="*",collapse="+"))
  
  write.table(apply(round(summary(reg)$coefficients[,c(1,2)],3),1,paste0,collapse=" / "),"clipboard")
  message("Reg table copied to clipboard")
  
  dev_preds<-predict(reg,type="response")
  dev_roc<-roc(as.vector(dev_data[,'event_bin']),dev_preds)
  message("Development C is ",dev_roc$auc)
  calibration_plot(dev_preds,dev_data[,'event_bin'])
  
  val_preds<-predict.glm(reg,newdata = val_data,type="response")
  ################IF you want to rcalibrate
  if(do_recalibrate) {message("Reclaibration done!"); val_preds<-recalibrate(val_preds,val_data[,'event_bin'])}
  val_roc<-roc(as.vector(val_data[,'event_bin']),val_preds)
  message("Validation C is ",val_roc$auc)
  calibration_plot(val_preds,val_data[,'event_bin'])
  
  plot(1-dev_roc$specificities,dev_roc$sensitivities,col='blue',xlim=c(0,1),ylim=c(0,1),lwd=2,type='l',xlab="False Positive",ylab="True Positive")
  lines(1-val_roc$specificities,val_roc$sensitivities,lwd=2,type='l')
  mres<-mROC(val_preds)
  lines(mres$FPs,mres$TPs,col="red",lwd=2,type='l')
  lines(c(0,1),c(0,1),col="gray",type='l')
  
  results$Pexac_p_dev<-mean(dev_preds)
  results$Oexac_p_dev<-mean(dev_data[,'event_bin'])
  results$Pexac_p_val<-mean(val_preds)
  results$Oexac_p_val<-mean(val_data[,'event_bin'])
  results$Dn=results$Oexac_p_val-results$Pexac_p_val
  results$ttest_p<-t.test(val_preds,val_data[,'event_bin'])$p.value
  
  x<-mROC_analysis(y=as.vector(val_data[,'event_bin']), p = val_preds, n_sim = settings$n_sim_inference_fast, inference = 1)

  return(c(results,x))
}







