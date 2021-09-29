#include <Rcpp.h>
#include <R.h>

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
  return((double)unif_rand());
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
  double *m = new double[n];
  int *y = new int[n];
  READ_R_VECTOR(M,m);
  READ_R_VECTOR(Y,y);
  mROC_stats out=calc_mROC_stats(n,m,y);
  double tmp[2];
  tmp[0]=out.A;
  tmp[1]=out.B;
  
  delete[] m;
  delete[] y;
  
  return(AS_VECTOR_DOUBLE_SIZE(tmp,2));
}










//Conditionaly samples from vector pof probabilities p, such that m individuals are selected.
int *conditional_sample(int n, double *p, int size)
{
  double *odds = new double[n];
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

  delete[] odds;
  
  return(int_buffer);
}
















mROC_stats *simulate_null_mROC_stats_unconditional(int n, double *M, int n_sim)
{
  int *y = new int[n];
  for(int i=0;i<n_sim;i++)
  {
    for(int j=0;j<n;j++)
    {
      if(rand_unif()<M[j]) y[j]=1; else y[j]=0;
    }
    mROC_stats_buffer[i]=calc_mROC_stats(n,M,y);
  }
  
  delete[] y;
  
  return(mROC_stats_buffer);
}





// [[Rcpp::export]]
NumericMatrix Csimulate_null_mROC_stats_unconditional(NumericVector M, int n_sim)
{
  int n=M.size();
  double *m = new double[n];
  std::copy(M.begin(),M.end(),m);

  mROC_stats *tmp=simulate_null_mROC_stats_unconditional(n,m,n_sim);

  NumericMatrix out(n_sim,2);
  for(int i=0;i<n_sim;i++)
  {
    out(i,0)=tmp[i].A;
    out(i,1)=tmp[i].B;
  }
  
  delete[] m;
  
  return(out);
}






