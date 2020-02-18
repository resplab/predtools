
############################OLD code (Remove?)#################################

OLD_stylized_sim<-function(sample_sizes=c(100,1000,10000),n_sim=1000)
{
  pi<-runif(100)
  y<-rbinom(100,1,pi)
  template<-unlist(mROC_inference(p=pi,y=y,CI=FALSE, n_sim = 100,fast=TRUE))
  
  out<-as.data.frame(matrix(NA, nrow = n_sim*length(template)*length(sample_sizes)*3,ncol=5))
  colnames(out)<-c("id","sample_size","scenario","Parameter","Value")
  
  index<-1
  id<-1
  for(i in 1:n_sim)
  {
    cat(".")
    for(j in 1:length(sample_sizes))
    {
      ss<-sample_sizes[j]
      
      pi=runif(ss)
      y=rbinom(ss,size = 1,prob = pi)
      tmp<-unlist(mROC_inference(p=pi,y=y,CI=FALSE, n_sim = 10000,fast=TRUE))
      out[index:(index+length(template)-1),1]<-id
      out[index:(index+length(template)-1),2]<-ss
      out[index:(index+length(template)-1),3]<-"Correct"
      out[index:(index+length(template)-1),4]<-names(tmp)
      out[index:(index+length(template)-1),5]<-tmp
      index<-index+length(template)
      
      pi<-sqrt(pi)
      tmp<-unlist(mROC_inference(p=pi,y=y,CI=FALSE, n_sim = 10000,fast=TRUE))
      out[index:(index+length(template)-1),1]<-id
      out[index:(index+length(template)-1),2]<-ss
      out[index:(index+length(template)-1),3]<-"^0.5"
      out[index:(index+length(template)-1),4]<-names(tmp)
      out[index:(index+length(template)-1),5]<-tmp
      index<-index+length(template)
      
      pi<-pi^4
      tmp<-unlist(mROC_inference(p=pi,y=y,CI=FALSE, n_sim = 10000,fast=TRUE))
      out[index:(index+length(template)-1),1]<-id
      out[index:(index+length(template)-1),2]<-ss
      out[index:(index+length(template)-1),3]<-"^2"
      out[index:(index+length(template)-1),4]<-names(tmp)
      out[index:(index+length(template)-1),5]<-tmp
      index<-index+length(template)
    }
  }
  
  aux$out<<-out
  
  return(list(out=out))
}






#Relies on aux$out
process_results<-function(dec_points=3)
{
  
  my_summary<-function(values,mu=T,se=T,p_reject=T)
  {
    out<-""
    if(mu==T) out<-paste0(out, beautify(mean(values)))
    if(se==T) {if(out!="") out=paste0(out,","); out<-paste0(out, beautify(sd(values)))}
    
    return(out)
  }
  
  beautify<-function(value)
  {
    return(format(round(value,dec_points),digits = dec_points,nsmall=dec_points))
  }
  
  data<-as.data.frame(aux$out)
  rows<-sqldf("SELECT DISTINCT sample_size FROM data")[,1]
  cols<-sqldf("SELECT DISTINCT scenario FROM data")[,1]
  
  out<-structure(paste(rep(1:(length(rows)*length(cols)))),.Dim=c(3,3))
  
  for(i in 1:length(rows))
  {
    this_row<-rows[[i]]
    for(j in 1:length(cols))
    {
      
      this_col<<-cols[[j]]
      
      ssql<-sprintf("SELECT value FROM data WHERE sample_size=%d AND scenario='%s'",this_row,this_col)
      q<-sqldf(paste0(ssql," AND parameter='stats.distance'"))[,1]
      out[i,j]<-paste0(out[i,j],"|zdelta|:",my_summary(q) , "")
      q<-sqldf(paste0(ssql," AND parameter='pvals.distance'"))[,1]
      out[i,j]<-paste0(out[i,j],beautify(sum(q<0.05)/length(q)))
      
      ssql<-sprintf("SELECT value FROM data WHERE sample_size=%d AND scenario='%s'",this_row,this_col)
      q<-sqldf(paste0(ssql," AND parameter='stats.surface'"))[,1]
      out[i,j]<-paste0(out[i,j],"S:",my_summary(q) , "")
      q<-sqldf(paste0(ssql," AND parameter='pvals.surface'"))[,1]
      out[i,j]<-paste0(out[i,j],beautify(sum(q<0.05)/length(q)))
      
      ssql<-sprintf("SELECT value FROM data WHERE sample_size=%d AND scenario='%s'",this_row,this_col)
      q<-sqldf(paste0(ssql," AND parameter='pval'"))[,1]
      out[i,j]<-paste0(out[i,j],"C:",beautify(sum(q<0.05)/length(q)))
    }
  }
  
  write.table(out,"clipboard", row.names = FALSE,  col.names = FALSE)
  
  return(out)
  
  message("Results are copied into clipboard.")
}











################################################Section 4: ECLIPSE##################################



analyze_eclipse<-function()
{
  eclipse<<-readRDS("validatedECLIPSE.RData")
  
  results<-list()
  results$Pexac_rate<-mean(eclipse$predicted_exac_rate)
  results$Oexac_rate<-mean(eclipse$obsExac_yr2)
  results$Pexac_prob<-mean(1-ppois(q=0,lambda=eclipse$predicted_exac_rate))
  results$Oexac_prob<-mean(eclipse$obsExac_yr2>=1)
  results$Dn=results$Pexac_prob-results$Oexac_prob
  results$ttest_p<-t.test(1-ppois(q=0,lambda=eclipse$predicted_exac_rate),eclipse$obsExac_yr2>=1)$p.value
  results$Pexac_rate_S<-mean(eclipse$predicted_severe_exac_rate)
  results$Oexac_rate_S<-mean(eclipse$obsExac_yr2)
  results$Pexac_prob_S<-mean(1-ppois(q=0,lambda=eclipse$predicted_severe_exac_rate))
  results$Oexac_prob_S<-mean(eclipse$obsExac_severe_yr2>=1)
  results$Dn_S=results$Pexac_prob_S-results$Oexac_prob_S
  results$ttest_p_S<-t.test(1-ppois(q=0,lambda=eclipse$predicted_severe_exac_rate),eclipse$obsExac_severe_yr2>=1)$p.value
  
  
  x<-mROC_analysis(p = 1-ppois(q=0,lambda=eclipse$predicted_exac_rate), y=eclipse$obsExac_yr2>=1, n_sim = settings$n_sim_inference_fast, inference = 1)
  #mROC_analysis(p = 1-ppois(q=0,lambda=eclipse$predicted_exac_rate), y=eclipse$obsExac_yr2>=1, n_sim = settings$n_sim_inference_slow, inference = 2, fast=FALSE)
  plot.new()
  calibration_plot(1-ppois(q=0,lambda=eclipse$predicted_exac_rate),eclipse$obsExac_yr2>=1,type_prob=TRUE)
  
  mROC_analysis(p = 1-ppois(q=0,lambda=eclipse$predicted_severe_exac_rate), y=eclipse$obsExac_severe_yr2>=1, n_sim = settings$n_sim_inference_fast, inference = 1)
  plot.new()
  calibration_plot(1-ppois(q=0,lambda=eclipse$predicted_severe_exac_rate),eclipse$obsExac_severe_yr2>=1,type_prob=TRUE)
}





































###LASSO!!!
#model_data<-model.matrix(formula,data=dev_data)
#reg_lasso<-cv.glmnet(model_data,y=dev_data[,'event_bin'][[1]],family="binomial")
#dev_preds<-predict(reg_lasso,newx=model_data, s="lambda.min",type="response")
#dev_roc<-roc((dev_data[,'event_bin'])[[1]],val_preds)
#calibration_plot(dev_preds,dev_data[,'event_bin'][[1]])

#tmp<-model.matrix(formula,data=val_data)
#val_preds<-predict(reg_lasso,newx=tmp, s="lambda.min",type="response")
#val_roc<-roc((val_data[,'event_bin'])[[1]],val_preds)
#calibration_plot(val_preds,val_data[,'event_bin'][[1]])

#plot(dev_roc,col='blue')
#lines(val_roc)
#mres<-mROC(val_preds)
#lines(1-mres[,1],mres[,2],col="red")









