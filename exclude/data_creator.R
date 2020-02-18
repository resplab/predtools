set.seed(3344)

create_data<-function(n,c_mus=c(0,1),c_sds=c(1,2),b_ps=c(0.25,0.5),covar_betas=c(b0=-1,b1=0.1,b2=0.2,b3=0.3,b4=0.4),covar_names=c("age","severity","sex","comorbidity"))
{
  covars<-matrix(NA, nrow=n, ncol=length(c_mus)+length(b_ps))
  
  #colnames(covars)<-c(paste0("x",1:length(c_mus)),paste0("c",1:length(b_ps)))
  colnames(covars)<-covar_names
  
  for(i in 1:length(c_mus))
  {
    covars[,i]<-rnorm(n,c_mus[i],c_sds[i])
  }
  
  for(i in 1:length(b_ps))
  {
    covars[,i+length(c_mus)]<-rbinom(n,size=1,prob=b_ps[i])
  }
  
  lin_pred<-cbind(1,covars)%*%covar_betas
  
  y<-rbinom(n,1,1/(1+exp(-lin_pred)))
  
  return(as.data.frame(cbind(covars,y=y)))
}



covar_betas=c(b0=-1,b1=0.5,b2=0.5,b3=1,b4=1)

dev_data<-create_data(500, covar_betas=covar_betas)
save(dev_data,file=paste0(getwd(),"/data/dev_data.RData"))

covar_betas[c('b2','b3')]=0.3; 
val_data<-create_data(400,covar_betas=covar_betas)
save(val_data,file=paste0(getwd(),"/data/val_data.RData"))


covars<-colnames(dev_data)
covars<-covars[-length(covars)]
formula<-paste0("y~",paste0(covars,collapse="+"))

dev_reg<-glm(formula,data=dev_data,family=binomial(link="logit"))
pred<-predict.glm(dev_reg, type='response')
plot(roc(dev_data[,'y'],pred))

pred2<-predict.glm(dev_reg,newdata = val_data, type='response')
lines(roc(val_data[,'y'],pred2),col='blue')


mres<-mROC(pred2)
lines(mres,col='red')

mROC_inference(val_data[,'y'],pred2)
