#include <Rcpp.h>
using namespace Rcpp;




enum errors
{
  ERR_INCORRECT_SETTING_VARIABLE=-1,
  ERR_INCORRECT_VECTOR_SIZE=-2,
  ERR_INCORRECT_INPUT_VAR=-3,
  ERR_EVENT_STACK_FULL=-4,
  ERR_MEMORY_ALLOCATION_FAILED=-5
} errors;
/*** R
errors<-c(
  ERR_INCORRECT_SETTING_VARIABLE=-1,
  ERR_INCORRECT_VECTOR_SIZE=-2,
  ERR_INCORRECT_INPUT_VAR=-3,
  ERR_EVENT_STACK_FULL=-4,
  ERR_MEMORY_ALLOCATION_FAILED=-5
)
*/




#define AS_VECTOR_DOUBLE(src) std::vector<double>(&src[0],&src[0]+sizeof(src)/sizeof(double))
#define AS_VECTOR_DOUBLE_SIZE(src,size) std::vector<double>(&src[0],&src[0]+size)

#define AS_VECTOR_BOOL_SIZE(src,size) std::vector<bool>(&src[0],&src[0]+size)

#define AS_MATRIX_DOUBLE(src)  array_to_Rmatrix(std::vector<double>(&src[0][0],&src[0][0]+sizeof(src)/sizeof(double)),sizeof(src[0])/sizeof(double))
#define AS_MATRIX_DOUBLE_SIZE(src,size)  array_to_Rmatrix(std::vector<double>(&src[0][0],&src[0][0]+size*sizeof(src[0])/sizeof(double)),sizeof(src[0])/sizeof(double))

#define AS_MATRIX_INT(src)  array_to_Rmatrix(std::vector<int>(&src[0][0],&src[0][0]+sizeof(src)/sizeof(int)),sizeof(src[0])/sizeof(int))
#define AS_MATRIX_INT_SIZE(src,size)  array_to_Rmatrix(std::vector<int>(&src[0][0],&src[0][0]+size*sizeof(src[0])/sizeof(int)),sizeof(src[0])/sizeof(int))

#define AS_VECTOR_INT(src) std::vector<int>(&src[0],&src[0]+sizeof(src)/sizeof(int))
#define AS_VECTOR_INT_SIZE(src,size) std::vector<int>(&src[0],&src[0]+size)

#define READ_R_VECTOR(src,dest) std::copy(src.begin(),src.end(),&dest[0])
#define READ_R_MATRIX(src,dest) {if(src.size()==sizeof(dest)/sizeof(dest[0][0])) {std::copy(src.begin(),src.end(),&dest[0][0]);} else return(ERR_INCORRECT_VECTOR_SIZE);}







double rand_unif()
{
  return((double)rand()/RAND_MAX);
}


typedef struct mROC_stats
{
  double A;
  double B;
} mROC_stats;


int int_buffer[1000000];
double double_buffer[1000000];
bool bool_buffer[1000000];
mROC_stats mROC_stats_buffer[1000000];


//Calculates the A (distance) and B (ROC equality) statistics for a vector of predicted probabilities (PI) and observed outcomes (Y)//
mROC_stats calc_mROC_stats(int n, double *M, int *Y)
{
  mROC_stats out;
  out.A=0;
  out.B=0;

  int n0=0,n1=0;
  double sumP0=0,sumP1=0;

  for(int i=0;i<n;i++)
  {
    if(Y[i]==0) ++n0; else ++n1;
    sumP1+=M[i];
    sumP0+=1-M[i];
    out.A+=Y[i]-M[i];
  }
  out.A=fabs(out.A/n);


  double xo=0, xe=0;
  double yo=0,ye=0;

  double step=0;
  double B=0;

  int io=n-1,ie=n-1;
  
  while(io>=0 && ie>=0)
  {
    if(xo<xe) //xo is behind and has to make a jump
    {
      if(Y[io]==1)
      {
        step=0;
        yo=yo+(double)1/n1;
      }
      else
      {
        step=(double)1/n0;
        //B=B+fabs(yo-ye)*std::min(step,xe-xo);
        B=B+fabs(yo-ye)*std::min(step,xe-xo);
      }
      xo=xo+step;
      io=io-1;
    }
    else //now xe is behind
    {
      step=(1-M[ie])/sumP0;
      B=B+fabs(yo-ye)*std::min(step,xo-xe);
      xe=xe+step;
      ye=ye+M[ie]/sumP1;
      ie=ie-1;
    }
  }
  
  out.B=B;
  return(out);
}






// [[Rcpp::export]] 
std::vector<double> Ccalc_mROC_stats(NumericVector M, NumericVector Y)
{
  int n=M.size();
  double m[n];
  int y[n];
  READ_R_VECTOR(M,m);
  READ_R_VECTOR(Y,y);
  mROC_stats out=calc_mROC_stats(n,m,y);
  double tmp[2];
  tmp[0]=out.A;
  tmp[1]=out.B;
  return(AS_VECTOR_DOUBLE_SIZE(tmp,2));
}










//Conditionaly samples from vector pof probabilities p, such that m individuals are selected.
int *conditional_sample(int n, double *p, int size)
{
  double odds[n];
  double sum_odds=0;

  for(int i=0;i<n;i++)
  {
    odds[i]=p[i]/(1-p[i]);
    sum_odds=sum_odds+odds[i];
    int_buffer[i]=0;
  }

  for(int i=0;i<size;i++)
  {
    double r=runif(1)[0]*sum_odds;
    int j=-1;
    double rsum=0;
    while(rsum<r)
    {
      ++j;
      if(int_buffer[j]==0)  rsum=rsum+odds[j];
    }
    int_buffer[j]=1;
    sum_odds=sum_odds-odds[j];
  }

  return(int_buffer);
}






//CRAZY: conditionally samples through try and error!!!
int *conditional_sample_crazy(int n, double *p, int size)
{
  int current_size=-1;
  int *y=int_buffer;

  while(current_size!=size)
  {
    current_size=0;
    for(int j=0;j<n;j++)
    {
      if(rand_unif()<p[j]) {y[j]=1; ++current_size;} else y[j]=0;
    }
  }

  return(y);
}














mROC_stats *simulate_null_mROC_stats_unconditional(int n, double *M, int n_sim)
{
  int y[n];
  for(int i=0;i<n_sim;i++)
  {
    for(int j=0;j<n;j++)
    {
      if(rand_unif()<M[j]) y[j]=1; else y[j]=0;
    }
    mROC_stats_buffer[i]=calc_mROC_stats(n,M,y);
  }
  return(mROC_stats_buffer);
}





// [[Rcpp::export]]
NumericMatrix Csimulate_null_mROC_stats_unconditional(NumericVector M, int n_sim)
{
  int n=M.size();
  double m[n];
  std::copy(M.begin(),M.end(),m);

  mROC_stats *tmp=simulate_null_mROC_stats_unconditional(n,m,n_sim);

  NumericMatrix out(n_sim,2);
  for(int i=0;i<n_sim;i++)
  {
    out(i,0)=tmp[i].A;
    out(i,1)=tmp[i].B;
  }
  return(out);
}







/*


int *UPMEsfromq(double *q, int n_row, int n_col)
{
  //printf("UPMEsfromq: You want me to sample %d from a vector of size %d\n",n_col,n_row);
  int *out=int_buffer;
  
  int col_index=n_col-1;
  
  for (int k=0;k<n_row;k++)
  {
    out[k]=0;
    //printf("%d,%d,%f\n",k,col_index,q[col_index*(n_row)+k]);
    if (col_index != -1)
    {
      if (rand_unif() < q[col_index*(n_row)+k]) 
      {
        out[k] = 1;
        col_index = col_index - 1;
      }
    }
  }
  
  int sum=0;
  for(int i=0;i<n_row;i++)  sum=sum+out[i];
  
  //printf("Final: %d,%d\n",col_index,sum);
  
  return(out);
}





// [[Rcpp::export]]
NumericVector CUPMEsfromq(NumericMatrix q)
//For test purposes.
{
  int n_row=q.nrow();
  int n_col=q.ncol();
  
  NumericVector outV(n_row);
  
  for (int k=0;k<n_row;k++)
  {
    outV(k)=0;
    if (n_col != -1)
    {
      //printf("%d,%d,%f\n",k,n_col,q(n_row,n_col-1));
      if (rand_unif() < q(k,n_col-1))
      {
        outV(k) = 1;
        n_col = n_col - 1;
      }
    }
  }
  
    return(outV);
}





mROC_stats *simulate_null_ds_conditional(double *q, int n_row, int n_col, int n_sim)
{
  //printf("simulate_null_ds_conditional: You want me to sample %d from a vector of size %d\n",n_col,n_row);
  int *y;
  for(int i=0;i<n_sim;i++)
  {
    y=UPMEsfromq(q, n_row, n_col);
    mROC_stats_buffer[i]=calc_mROC_stats(n_row,q,y);
  }
  return(mROC_stats_buffer);
}




// [[Rcpp::export]]
NumericMatrix Csimulate_null_ds_conditional(NumericMatrix q, int n_sim)
{
//Number of rows in q is the size of vector; number of oclumns is desdired sample size
  int n_row=q.nrow();
  int n_col=q.ncol();
  
  //printf("Csimulate_null_ds_conditional: You want me to sample %d from a vector of size %d\n",n_col,n_row);
  
  double _q[n_row*n_col];
  std::copy(q.begin(),q.end(),_q);
  
  //return(NumericMatrix(n_row,n_col,_q));
  
  mROC_stats *tmp=simulate_null_ds_conditional(_q, n_row, n_col, n_sim);
  
  NumericMatrix out(n_sim,2);
  for(int i=0;i<n_sim;i++)
  {
    
    out(i,0)=tmp[i].A;
    out(i,1)=tmp[i].B;
  }
  return(out);
}






mROC_stats *simulate_null_ds_conditional_crazy(double *p, int vlength, int size, int n_sim)
{
  //printf("simulate_null_ds_conditional: You want me to sample %d from a vector of size %d\n",n_col,n_row);
  int *y;
  for(int i=0;i<n_sim;i++)
  {
    y=conditional_sample_crazy(vlength, p, size);
    mROC_stats_buffer[i]=calc_mROC_stats(vlength,p,y);
  }
  return(mROC_stats_buffer);
}




// [[Rcpp::export]]
NumericMatrix Csimulate_null_ds_conditional_crazy(NumericVector p, int size, int n_sim)
{
  double _p[p.length()];
  
  std::copy(p.begin(),p.end(),_p);
  
  mROC_stats *tmp=simulate_null_ds_conditional_crazy(_p, p.length(), size, n_sim);
  
  NumericMatrix out(n_sim,2);
  for(int i=0;i<n_sim;i++)
  {
    
    out(i,0)=tmp[i].A;
    out(i,1)=tmp[i].B;
  }
  return(out);
}






// [[Rcpp::export]]
std::vector<int> Cconditional_sample(NumericMatrix q)
{
  int n_row=q.nrow();
  int n_col=q.ncol();
  
  double _q[n_row*n_col];
  std::copy(q.begin(),q.end(),_q);
  
  int *out=int_buffer;
  out=UPMEsfromq(_q, n_row, n_col);
  
  return(AS_VECTOR_INT_SIZE(out,n_row));
}




// [[Rcpp::export]]
std::vector<int> Cconditional_sample_crazy(NumericVector M, int size)
{
  int n=M.size();
  double m[n];
  std::copy(M.begin(),M.end(),m);
  int *tmp=conditional_sample_crazy(n,m,size);
  
  return(AS_VECTOR_INT_SIZE(tmp,n));
}


*/

