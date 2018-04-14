#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "uthash.h"
#include "utarray.h"
//#include "utarray.h"

#define EPSILON 0.000001
#define MAXCAT 50

// structures for fast historical splitting
struct nid{
  int id;
  double y;
  double time;
  UT_hash_handle hh;
};





// structures for keeping track of node affiliation
struct rnumber{
  int row_number;
  UT_hash_handle hh;
}; 

struct node{
  int node_number;
  int nrows; // number of rows stored in hash
  UT_array *rows_array;
  struct rnumber *rows;
  UT_hash_handle hh;
};

struct nodelist{ // used to store multiple node structures
  int list_number;
  struct node *nodes;
  UT_hash_handle hh;

};

struct node *nodes=NULL;
struct node *oob=NULL; // for out-of-bag 
struct node *id_rows=NULL; // row numbers associated with each id


struct nodelist *tnode_row=NULL; // Stores rows corresponding each terminal node for each tree, used for quantile estimation



// ---------------------------------------------------------



struct split_result{
  // structure for packing results of a split 
  int success; // 0 if unsuccessful, else 1
  int varindx;
  double nL;
  double nR;		
  double predL;
  double predR;
  double point;
  double error;
  double error_reduction;
  double tau;
  double delta;
  double delta0; // used in window summary (ie prev obstime within [delta0,delta] current obstime)
  //
  double sumL;
  double sumR;
  double sumsqL;
  double sumsqR;
  int missing; // used with "missing together" approach    
               // also used for indicating whether split is done
               // on number of observations or fraction of obs above/below
               // some threshold (missing=1, =0 respectively)

  // CLASSIFICATION
  double sumL_cat[MAXCAT];
  double sumR_cat[MAXCAT];
  
  
};

int find_tnode(int tree_indx,int row_number);


void print_split_result(struct split_result *sr);
void create_history_summary_uthash(int nodeid,int vindx,struct split_result *sr,double *x);
void permute(int *p,int n);
void create_history_summary_meansum(int nodeid,int vindx,struct split_result *sr,double *x);
void create_meansum_summary(int nodeid,int vindx,double *delta,int ndelta,int *nobs);
void create_history_summary_tilde(int nodeid,int vindx,struct split_result *sr,double *x,double *nx);
void create_history_summary_window(int nodeid,int vindx,struct split_result *sr,double *x);
void create_history_summary_meansum0(int nodeid,int vindx,struct split_result *sr,double *x);
void create_history_summary_meansum0_window(int nodeid,int vindx,struct split_result *sr,double *x);
void create_history_summary_tilde_window(int nodeid,int vindx,struct split_result *sr,double *x,double *nx);
void delta_window(int nodeid,double *delta,double *delta0);
void assign_method(int *method); 
int majority_vote(double *x);


struct data_options{

  // .. testing data
  double **split_matrix;

  //
  int nsamp_old;

  // delta
  int set_delta;
  double *delta_unique;
  int window_delta; // regulates how (delta0,delta)-candidates are generated
  int window_summary; // 1 if window summary, else 0 
  int n_max;  // grid length for tilde splitting
  // total loop count 
  int ncc;
  // .. testing data finished.... 


  int ncat; // number of categories if classification
  
  double **x;
  double *response;
  int n;
  int p;
  int yindx;
  double *time;
  int *id;

  double *lower;
  double *upper;

  double lower_time;
  double upper_time;

  // scratch space
  double *daux1,*daux2,*daux3,*daux4,*daux5,*daux6;
  int *iaux1,*iaux2,*iaux3,*iaux4,*iaux5,*iaux6,*iaux7;  //iaux7 added, used for classification response

  // Split info for fast splitting
  //  int *nstat;
  //  double *ssqL,*ssqR,*sumL,*sumR,*nL,*nR;
 

  
  // predictions
  double *predictions;
  double *predictions_response;
  int *oob;
  double **predictions_gini;



  double mean_y;
  // split info and parameters 
  int min_nodesize;
  int nsamp;
  int nsamp_window;  // for window method .. 
  int nsamp_var;
  double lambda;
  double best_sumsq;
  double best_frac;
  double best_cut;
  double best_delta;
  int best_n_right;
  int best_n_left;
  int vtype;
  int valid_split;
  int subsample;
  double sample_fraction;
  int random_split;
  int method;
  int method_family;
  
  int *node;
  int nboost;
  int nsplit;
  int rf;  // if rf=1, then random forest, else boosting
  int mtry;
  int time_split; // determines whether splitting is done on time ... if =1, then yes, otherwise no 
  double delta_multiple; // multiply delta   
  // oob matrix
  //int **oob;
  int max_obs_history;
  
  // Tree matrix
  double **tree_matrix;
  int row_counter;
  int tree_counter;
  int *tree_start;
  int *tree_end;

  int *train;

  double **varimp; // varimp[k][j] is variable importance of k-th variable on j-th tree. 

  // information regarding historical and concurrent variable status
  int *splitvar_history;
  int *splitvar_concurrent;	
  int nvar_history;
  int nvar_concurrent;
  double **n_row_number;
  double **n_row_accumulator;
  int **n_changed;       // keeps track of which n_i(tau) change (for given tau)
  int *counter_changed; // keeps track of number of n_i(tau) changed (given tau)

  struct split_result **sr; 


  // mean/sum historical summaries
  double **xmat;
  double **nxmat;
  int meansummary;
  int *column_double;
  double *y_meansum;
  double LARGE_NUMBER; // Used for "Missing Together" splitting 
  
};

struct data_options daop;

void sampleWOR(int n,int size,int *res)
{
  int i,k,j,nElements;
  int *sampVec=(int*) malloc(n*sizeof(int));


  double nA,h;

  for(i=1;i<=n;i++)
    sampVec[i-1]=i-1;

  nElements=n;
  nA=(double) n;
 
  for(i=0;i<size;i++)
    {

      h=nA*unif_rand(); 
#ifdef DEBUG
      Rprintf("sampleWOR-function: unif=%lf \n",h);
#endif /*DEBUG*/

      k=(int) h;
      res[i]=sampVec[k];

      for(j=(k+1);j<nElements;j++)
	{
	  sampVec[j-1]=sampVec[j];
	}
      nElements=nElements-1;
      nA=nA-1;
    }

  free(sampVec);

}



void set_variable_status(int *vhistory,int *nhistory,int *vconcurrent,int *nconcurrent)
{
  int k;

  for(k=0;k<nhistory[0];k++)
    daop.splitvar_history[k]=vhistory[k];
  daop.nvar_history=nhistory[0];


  for(k=0;k<nconcurrent[0];k++)
    daop.splitvar_concurrent[k]=vconcurrent[k];
  daop.nvar_concurrent=nconcurrent[0];

}

void set_mtry(int *mtry)
{
  daop.mtry=mtry[0];
}


void set_minnodesize(int *n)
{
  daop.min_nodesize=n[0];
  //Rprintf("daop.min_nodesize=%d \n",daop.min_nodesize);
}


void set_method_family(int *n)
{
  daop.method_family=n[0];
  //Rprintf("daop.min_nodesize=%d \n",daop.min_nodesize);
}




void initialize_delta(int *ndelta,double *delta_u)
{
  int k;
  
  daop.set_delta=ndelta[0];
  daop.delta_unique=(double*) malloc(sizeof(double)*ndelta[0]); //freed
  for(k=0;k<daop.set_delta;k++)
    daop.delta_unique[k]=delta_u[k]+EPSILON;

 

  if(daop.set_delta<daop.nsamp){
    daop.nsamp=daop.set_delta;
    daop.nsamp_window=daop.nsamp*(daop.nsamp+1)/2;
  }
    

}

void read_data(double *x,int *n,int *p,double *time,int *id,int *yindx,double *lambda,int *nsplit,int *nboost,int *rf,int *nsamp,int *time_split,int *window_summary)
{
  int i,j,ncol,k;
  struct data_options *d;
  double *v;
  int max_history_length;
  int nsampw;
  d=&daop;


  daop.method=1;
  daop.method_family=1;
  daop.meansummary=1;
  daop.window_delta=1;
  daop.window_summary=window_summary[0];
  
  daop.ncc=0;

  daop.ncat=-10; // not classification, until set. 
  
  /* // start: testing data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  /* d->split_matrix=(double**) malloc(sizeof(double*)*(1000)); */
  /* v=(double*) malloc(sizeof(double)*1000*13); */
  /* d->split_matrix[0]=v; */
  /* for(i=1;i<1000;i++) */
  /*   d->split_matrix[i]=d->split_matrix[i-1]+13; */
  /* // end: testing data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */



  
  // .. get max number of subject level observations
  i=(-99);
  j=id[0];
  ncol=1;  
  for(k=1;k<n[0];k++)
    {
      if(id[k]==j){
	ncol++;
      }else{
	if(ncol>i)
	  i=ncol;
	ncol=1;
	j=id[k];
      }
   
    }
  
  if(ncol>i)  // last id 
    i=ncol;   // 'i' now has max number of observations of any 'id'

  daop.max_obs_history;


  nsampw=nsamp[0];
  if(window_summary[0]==1)
    nsampw=nsampw*(nsampw+1)/2;

  // n_max
  d->n_max=10; 
  
  // matrix of split_result structures
  if(i<50)
    i=50;   // needs to be changed ... for tilde-spliting
  
  d->sr=(struct split_result**) malloc(sizeof(struct split_result*)*nsampw);
  for(k=0;k<nsampw;k++)
    d->sr[k]=(struct split_result*) malloc(sizeof(struct split_result)*(i+2)); // way to many... //fr


  d->n_changed=(int**) malloc(sizeof(int*)*nsampw);       // keeps track of which n_i(tau) change (for given tau)
  for(k=0;k<nsampw;k++)
    d->n_changed[k]=(int*) malloc(sizeof(int)*i*n[0]);  //fr
  
  d->counter_changed=(int*) malloc(sizeof(int)*nsampw); // keeps track of number of n_i(tau) changed (given tau) //fr

  d->n_row_number=(double**) malloc(sizeof(double*)*nsampw); //fr
  for(k=0;k<nsampw;k++)
    d->n_row_number[k]=(double*) malloc(sizeof(double)*i*n[0]);

  d->n_row_accumulator=(double**) malloc(sizeof(double*)*nsampw); //fr
  for(k=0;k<nsampw;k++)
    d->n_row_accumulator[k]=(double*) malloc(sizeof(double)*i*n[0]);


  
  max_history_length=i;
  // -------------------------------------------------------
  ncol=15+MAXCAT;  // incremented for window-summaries (prior was 14)

  d->splitvar_history=(int*) malloc(sizeof(int)*p[0]); //fr
  d->splitvar_concurrent=(int*) malloc(sizeof(int)*p[0]);//fr
  for(k=0;k<p[0];k++)
    {
      d->splitvar_history[k]=k;
      d->splitvar_concurrent[k]=k;
    }
  
  d->nvar_concurrent=p[0];
  d->nvar_history=p[0];
  
  
  d->delta_multiple=2; 
  d->time_split=time_split[0];
  d->mtry=p[0];
  d->lambda=lambda[0]; // shrinkage parameter 
  d->nsplit=nsplit[0];
  d->nboost=nboost[0];
  // tree matrix
  j=(d->nboost*(3*(d->nsplit)));
  d->tree_matrix=(double**) malloc(sizeof(double*)*(d->nboost*(3*(d->nsplit))));
  v=(double*) malloc(sizeof(double)*j*ncol);
  d->tree_matrix[0]=v;
  for(i=1;i<j;i++)
    d->tree_matrix[i]=d->tree_matrix[i-1]+ncol;  //fr

  for(i=0;i<j;i++)
    for(k=0;k<ncol;k++)
      d->tree_matrix[i][k]=0;


  

  
  // Rprintf(" read in data \n");
  d->varimp=(double**) malloc(sizeof(double*)*p[0]);
  v=(double*) malloc(sizeof(double)*d->nboost*p[0]);
  d->varimp[0]=v;
  for(i=1;i<p[0];i++)
    d->varimp[i]=d->varimp[i-1]+d->nboost; //fr
  
  d->row_counter=-1; // this is incremented as rows are added
  d->tree_counter=0;

  d->tree_start=(int*) malloc(sizeof(int)*d->nboost); //fr
  d->tree_end=(int*) malloc(sizeof(int)*d->nboost); //fr

  d->rf=rf[0];
  if(d->rf==1)
    {
      d->subsample=1;
      d->sample_fraction=.5;
    }else{
    d->subsample=0;
    d->sample_fraction=1;
  }
  
  d->nsamp=nsamp[0];
  d->nsamp_old=d->nsamp;  // nsamp can be changed 
  d->nsamp_var=d->nsamp*d->nsamp*d->nsamp;
  d->min_nodesize=5;
  
  d->p=p[0];
  d->n=n[0];
  d->x=(double**) malloc(sizeof(double*)*d->p);  //fr
  for(i=0;i<d->p;i++)
    {
    d->x[i]=(double*) malloc(sizeof(double)*d->n);
    for(j=0;j<d->n;j++)
      d->x[i][j]=x[(d->n)*i+j];
    }

  // response 
  d->response=(double*) malloc(sizeof(double)*d->n);
  for(i=0;i<d->n;i++)
    d->response[i]=d->x[yindx[0]][i];  //fr


  daop.mean_y=0;
  for(i=0;i<d->n;i++)
    daop.mean_y+=d->response[i];

  daop.mean_y=daop.mean_y/((double) d->n);
  
  d->train=(int*) malloc(sizeof(int)*d->n);//fr
     for(j=0;j<d->n;j++)
    d->train[j]=1;

   
     d->time=(double*) malloc(sizeof(double)*d->n);//fr
  d->yindx=yindx[0];
  for(j=0;j<d->n;j++)
    d->time[j]=time[j];

  d->id=(int*) malloc(sizeof(int)*d->n);//fr
  for(j=0;j<d->n;j++)
    d->id[j]=id[j];

  d->predictions=(double*) malloc(sizeof(double)*d->n);//fr
  for(j=0;j<d->n;j++)
    d->predictions[j]=0;


  d->predictions_gini=(double**) malloc(sizeof(double)*d->n); //fr
    for(j=0;j<d->n;j++)
      d->predictions_gini[j]=(double*) malloc(sizeof(double)*MAXCAT);
    
    for(j=0;j<d->n;j++){
      for(i=0;i<MAXCAT;i++)
	d->predictions_gini[j][i]=0;
    }

    d->predictions_response=(double*) malloc(sizeof(double)*d->n);
    for(j=0;j<d->n;j++)  //fr
    d->predictions_response[j]=0;

// oob counter vector 
    d->oob=(int*) malloc(sizeof(int)*d->n);//fr
  for(j=0;j<d->n;j++)
    d->oob[j]=0;
  
  // allocate and get range of predictor variables 
  d->lower=(double*) malloc(sizeof(double)*d->p);//fr
  d->upper=(double*) malloc(sizeof(double)*d->p);//fr

  for(j=0;j<d->p;j++)
    {
      d->lower[j]=d->x[j][0];
      d->upper[j]=d->x[j][0];
      
      for(i=0;i<d->n;i++)
	{
	  if(d->x[j][i]>d->upper[j])
	    d->upper[j]=d->x[j][i];

	  if(d->x[j][i]<d->lower[j])
	    d->lower[j]=d->x[j][i];	  
	}
    }
  
  // range for time
  d->lower_time=d->time[0];
  d->upper_time=d->time[0];
  for(i=0;i<d->n;i++)
    {
      if(d->lower_time>d->time[i])
	d->lower_time=d->time[i];

     if(d->upper_time<d->time[i])
	d->upper_time=d->time[i];
    }
  

  // node vector
  d->node=(int*) malloc(sizeof(int)*d->n);//fr
  for(i=0;i<d->n;i++)
    d->node[i]=0;

  
  // scratch space  fr-all
  d->daux1=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux2=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux3=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux4=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux5=(double*) malloc(sizeof(double)*d->n*max_history_length);
  d->daux6=(double*) malloc(sizeof(double)*d->n*max_history_length);

  d->iaux1=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux2=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux3=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux4=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux5=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux6=(int*) malloc(sizeof(int)*d->n*max_history_length);
  d->iaux7=(int*) malloc(sizeof(int)*d->n*max_history_length);


  // mean/sum historical summary
  // WORKING HERE 
  d->xmat=(double**) malloc(sizeof(double*)*(2*nsampw));
  d->nxmat=(double**) malloc(sizeof(double*)*(2*nsampw));
  for(i=0;i<(nsampw*2);i++){
      d->xmat[i]=(double*) malloc(sizeof(double)*d->n);
      d->nxmat[i]=(double*) malloc(sizeof(double)*d->n);

  }
  d->column_double=(int*) malloc(sizeof(int)*nsampw);  
  d->y_meansum=(double*) malloc(sizeof(double)*d->n);

 

  d->LARGE_NUMBER=999999;

  // DEBUG
  //Rprintf("Exiting 'read_data'\n");
}

void set_LARGE_NUMBER(double *x)
{
  daop.LARGE_NUMBER=x[0];
}
  

void set_n_max(int *n_max)
{
  int nm;

  nm=n_max[0];
  if(n_max[0]>50)
    {
      nm=50;
      Rprintf("n_max must be <=50, setting n_max=50.\n");
    }

  if(n_max[0]<=0)
    {
      nm=10;
    }

  daop.n_max=nm;

}


void get_tree(int *tree_number,double *res,int *n)
{
  int i,j,k;
  i=0;
  for(k=0;k<=daop.row_counter;k++)
    {

      if(daop.tree_matrix[k][0]==tree_number[0])
	{
	  for(j=0;j<=14;j++)
	    {
	      res[i]=daop.tree_matrix[k][j];
	      i++;
	    }
	}

    }
  n[0]=i;
}


void get_tree_all(int *tree_number,double *res,int *n)
{
  int i,j,k;
  i=0;
  for(k=0;k<=daop.tree_end[tree_number[0]];k++)
    {
	  for(j=0;j<=14;j++)
	    {
	      res[i]=daop.tree_matrix[k][j];
	      i++;
	    }
    }
  n[0]=i;
}

void get_tree_all_gini(int *tree_number,int *ncat,double *res,int *n)
{
  int i,j,k;
  i=0;
  for(k=0;k<=daop.tree_end[tree_number[0]];k++)
    {
      for(j=0;j<=(14+ncat[0]);j++)
	    {
	      res[i]=daop.tree_matrix[k][j];
	      i++;
	    }
    }
  n[0]=i;
}

void set_ncat(int *ncat)
{

  daop.ncat=ncat[0];

}





double test_error()
{
  double error,h;
  int i,ntest;
  ntest=0;
  error=0;

  if(daop.rf!=1)
    {
  for(i=0;i<daop.n;i++)
    {
      if(daop.train[i]==0)
	{
	  error+=(daop.response[i]*daop.response[i]);
	  ntest++;
	}
    }
   }else{
    
  for(i=0;i<daop.n;i++)
    {
      if(daop.oob[i]>0)
	{
	  h=daop.response[i]-daop.predictions[i]/((double) daop.oob[i]);
	  error+=(h*h);
	  ntest++;
	}
    }
  }
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}




double test_error_gini()
{
  double error,h;
  int i,ntest,yint,hp;
  ntest=0;
  error=0;

  if(daop.rf!=1)
    {
  for(i=0;i<daop.n;i++)
    {
      if(daop.train[i]==0)
	{
	  error+=(daop.response[i]*daop.response[i]);
	  ntest++;
	}
    }
   }else{
    
  for(i=0;i<daop.n;i++)
    {
      if(daop.oob[i]>0)
	{
	  hp=majority_vote(daop.predictions_gini[i]);
	  h=1;
	  yint=(int) daop.response[i];
	  if(yint==hp)
	    h=0;
	  // h=daop.response[i]-daop.predictions[i]/((double) daop.oob[i]);
	  error+=h; // (h*h);
	  ntest++;
	}
    }
  }
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}


int majority_vote(double *x)
{
  int class,k;
  double mv;
  
  class=0;
  mv=0;
  for(k=0;k<daop.ncat;k++)
    {
      if(mv<x[k])
	{
	  mv=x[k];
	  class=k;
	}

    }
  return(class);

}

void get_data_options_info(int *n,int *p,int *b)
{
  n[0]=daop.n;
  p[0]=daop.p;
  b[0]=daop.nboost;
}

void subsample()
{
  int j;
 GetRNGstate();

 for(j=0;j<daop.n;j++)
   {
     if(unif_rand()<=daop.sample_fraction)
       {
	 daop.train[j]=1;
       }else{
       daop.train[j]=0;
     }
   }

 PutRNGstate();
}








void set_train(int *train)
{
  int i;
  for(i=0;i<daop.n;i++)
    daop.train[i]=train[i];
}


void free_daop()
{
 int i,j;
 struct data_options *d;

 daop.nsamp=daop.nsamp_old;
 if(daop.window_summary==1)
   daop.nsamp=daop.nsamp*(daop.nsamp+1)/2;

 d=&daop;

 if(daop.set_delta>0)
   free(daop.delta_unique);
 
 for(i=0;i<d->nsamp;i++)
   free(d->n_row_number[i]);
 free(d->n_row_number);

 for(i=0;i<d->nsamp;i++)
   free(d->n_row_accumulator[i]);
 free(d->n_row_accumulator);

 
  for(i=0;i<d->nsamp;i++)
   free(d->sr[i]);
 free(d->sr);

  for(j=0;j<daop.nsamp;j++)
    free(d->n_changed[j]);
  free(d->n_changed);

  free(d->counter_changed); 
 free(d->oob);

 free(d->predictions);
 free(d->predictions_response);
 free(d->splitvar_history);
 free(d->splitvar_concurrent);

 for(i=0;i<d->n;i++)
   free(d->predictions_gini[i]);
 free(d->predictions_gini);
 
 free(d->tree_matrix[0]);
 free(d->tree_matrix);

 free(d->train); 
   free(d->tree_start);
  free(d->tree_end);
   

  
  for(i=0;i<d->p;i++)
      free(d->x[i]);
  free(d->x);

  // varimp
  free(d->varimp[0]);
  free(d->varimp);
  
  // response 
  free(d->response);

  free(d->time);
  free(d->id);

  // allocate and get range of predictor variables 
  free(d->lower);
  free(d->upper);

  // node vector
  free(d->node);

 
  // scratch space
  free(d->daux1); 
  free(d->daux2); 
  free(d->daux3); 
  free(d->daux4); 
  free(d->daux5); 
  free(d->daux6);

 
  free(d->iaux1);
  free(d->iaux2);
  free(d->iaux3);
  free(d->iaux4);
  free(d->iaux5);
  free(d->iaux6);
  free(d->iaux7);


  // mean/sum
  free(d->y_meansum);
  for(i=0;i<(d->nsamp*2);i++)
    {
      free(d->xmat[i]);
      free(d->nxmat[i]);
    }
  free(d->xmat);
  free(d->nxmat);
  free(d->column_double);

}


void set_lambda(double *lambda)
{
  daop.lambda=lambda[0];

}


void set_data(double *x,int *n,int *p,double *time,int *id,int *yindx)
{
  // Delete old data, set new data (for prediction)
  //void read_data(double *x,int *n,int *p,double *time,int *id,int *yindx)

  struct data_options *d;
  int i,j,k;

  d=&daop;

  for(i=0;i<d->p;i++)
      free(d->x[i]);
  free(d->x);

  // response 
  free(d->response);

  free(d->time);
  free(d->id);
  free(d->node);
  
  d->p=p[0];
  d->n=n[0];
  d->x=(double**) malloc(sizeof(double*)*d->p);
  for(i=0;i<d->p;i++)
    {
    d->x[i]=(double*) malloc(sizeof(double)*d->n);
    for(j=0;j<d->n;j++)
      d->x[i][j]=x[(d->n)*i+j];
    }

  // response 
  d->response=(double*) malloc(sizeof(double)*d->n);
  for(i=0;i<d->n;i++)
    d->response[i]=d->x[yindx[0]][i];

  
  d->time=(double*) malloc(sizeof(double)*d->n);
  d->yindx=yindx[0];
  for(j=0;j<d->n;j++)
    d->time[j]=time[j];

  d->id=(int*) malloc(sizeof(int)*d->n);
  for(j=0;j<d->n;j++)
    d->id[j]=id[j];

d->node=(int*) malloc(sizeof(int)*d->n);
 
  
}



double tsummary_row_fast(int rownumber,double cut,double delta,int vindx)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta) // used to be j-i , what???? 
	    {
	      //nrec++;

	      if(d->x[vindx][i-j]<cut)
		{
		  nrec_condition++;
		}

	    }else{
	    break;
	  }
	}


      ff=nrec_condition;
      return(ff);
  
}







double tsummary_row_fast_window(int rownumber,double cut,double delta,double delta0,int vindx)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  double dd;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  
	  dd=(d->time[i]-d->time[i-j]);
	  if(dd<=delta&&dd>=delta0) // used to be j-i , what???? 
	    {
	      //nrec++;

	      if(d->x[vindx][i-j]<cut)
		{
		  nrec_condition++;
		}

	    }else{
	    if(dd>delta)
	      break;
	  }
	}


      ff=nrec_condition;
      return(ff);
  
}




double tsummary_row_permuted(int rownumber,double cut,double delta,int vindx,int *permuted_row)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta) // used to be j-i , what???? 
	    {
	      nrec++;

	      if(d->x[vindx][permuted_row[i-j]]>=cut)
		{
		  nrec_condition++;
		}

	    }else{
	    break;
	  }
	}


      ff=nrec_condition/(nrec+.00001);
      return(ff);
  
}





double tsummary_row_wrapper(int k,int row_number,int split_var);



double find_tnode_predict_fast(int tree_indx,int row_number)
{
  // For subject in row_number, find terminal node for subject

  int k,tnode,go,j,i,split_var;
  double ff;
  int counter;
  double pred;

  tnode=-99;
  //for(k=daop.tree_start[tree_indx];k<=daop.tree_end[tree_indx];k++)
  k=daop.tree_start[tree_indx];
  //Rprintf(" Row corresponding tree start %d \n",k);
  go=1;
  counter=0;
   while(go==1)
    {

      //Rprintf(" split-var =%d \n",((int) daop.tree_matrix[k][4]));
      if(daop.tree_matrix[k][4]<(-1))
	{
	  // terminal node
	  pred=(daop.tree_matrix[k][3]);
	  break;
	}else{
	split_var=((int) daop.tree_matrix[k][4]);

	if(daop.tree_matrix[k][5]<(-1)) // not a summary variable
	  {
	    if(daop.x[split_var][row_number]<daop.tree_matrix[k][6])
	      {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }

	  }else{
	  //ff=tsummary_row_fast(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
	  ff=tsummary_row_wrapper(k,row_number,split_var);
	  
	  if(ff<daop.tree_matrix[k][5])
	    {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }
	}

      }
    }

   return(pred);
}






void find_tnode_predict_gini(int tree_indx,int row_number,double *pred)
{
  // For subject in row_number, find terminal node for subject

  int k,tnode,go,j,i,split_var,m;
  double ff;
  int counter;
  //double pred;

  tnode=-99;
  //for(k=daop.tree_start[tree_indx];k<=daop.tree_end[tree_indx];k++)
  k=daop.tree_start[tree_indx];
  //Rprintf(" Row corresponding tree start %d \n",k);
  go=1;
  counter=0;
   while(go==1)
    {

      //Rprintf(" split-var =%d \n",((int) daop.tree_matrix[k][4]));
      if(daop.tree_matrix[k][4]<(-1))
	{
	  // terminal node
	  for(m=0;m<daop.ncat;m++)
	    pred[m]=(daop.tree_matrix[k][15+m]);

	  
	  break;
	}else{
	split_var=((int) daop.tree_matrix[k][4]);

	if(daop.tree_matrix[k][5]<(-1)) // not a summary variable
	  {
	    if(daop.x[split_var][row_number]<daop.tree_matrix[k][6])
	      {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }

	  }else{
	  //ff=tsummary_row_fast(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
	  ff=tsummary_row_wrapper(k,row_number,split_var);
	  
	  if(ff<daop.tree_matrix[k][5])
	    {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }
	}

      }
    }

   //return(pred);
}





double tsummary_row_tilde(int rownumber,double cut,double delta,int vindx,int nobs)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta) // used to be j-i , what???? 
	    {
	      nrec++;

	      if(d->x[vindx][i-j]<cut)
		{
		  nrec_condition++;
		}

	    }else{

	    break;
	  }
	}

      // nobs=1 if splitting is on number of observations 
	    if(nobs==0)
	      {
		if(nrec>EPSILON)
		  {
		    ff=(nrec_condition)/(nrec);
		  }else{
		  ff=1; // changed from 0 to 1, 2/3/2018
		}
	      }else{
	      // splitting is on number of observations
	      ff=nrec;
	    }
     
      return(ff);
  
}




double tsummary_row_tilde_window(int rownumber,double cut,double delta,double delta0,int vindx,int nobs)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff,dd;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  
	  dd=(d->time[i]-d->time[i-j]);
	  if((dd<=delta)&&(dd>=delta0)) // used to be j-i , what???? 
	    {
	      nrec++;

	      if(d->x[vindx][i-j]<cut)
		{
		  nrec_condition++;
		}

	    }else{
	    if(dd>delta)
	      break;
	  }
	}

      // nobs=1 if splitting is on number of observations 
	    if(nobs==0)
	      {
		if(nrec>EPSILON)
		  {
		    ff=(nrec_condition)/(nrec);
		  }else{
		  ff=1; // changed from 0 to 1, 2/3/2018
		}
	      }else{
	      // splitting is on number of observations
	      ff=nrec;
	    }
     
      return(ff);
  
}







double tsummary_row_meansum(int rownumber,double cut,double delta,int vindx,int nobs)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta) // used to be j-i , what???? 
	    {
	      nrec++;
	      nrec_condition+=d->x[vindx][i-j];
		

	    }else{


	  
	    break;
	  }
	}

		if(nrec>EPSILON)
		  {
		    if(daop.meansummary==1){
		    ff=(nrec_condition)/(nrec);
		    }else{
		      ff=nrec_condition;
		    }
		  }else{
		  if(nobs==-1)
		    {
		      // Missing is imputed by large negative number
		      ff=-daop.LARGE_NUMBER;
		    }else{
		    // Missing imputed by large positive number 
		    ff=daop.LARGE_NUMBER;
		  }
		}
     
      return(ff);
  
}









double tsummary_row_meansum0(int rownumber,double cut,double delta,int vindx,int nobs)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta) // used to be j-i , what???? 
	    {
	      nrec++;
	      nrec_condition+=d->x[vindx][i-j];
		

	    }else{


	  
	    break;
	  }
	}


      if(nobs==0)
	{
		if(nrec>EPSILON)
		  {
		 
		    if(daop.meansummary==1){
		    ff=(nrec_condition)/(nrec);
		    }else{
		      ff=nrec_condition;
		    }
		  }else{
		  ff=0;
		}
	}else{
	if(nobs==1)
	  ff=nrec;
      }
      return(ff);
  
}





double tsummary_row_meansum0_window(int rownumber,double cut,double delta,double delta0,int vindx,int nobs)
{
  // ff=tsummary_row(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);
  // Computes tsummary statistics for a single row (rownumber) 

 
  int i,j,cid;
  double nrec,nrec_condition,ff;
  double dd;
  
  struct data_options *d;
  d=&daop;

  i=rownumber;


      nrec=0;
      nrec_condition=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  
	  dd=(d->time[i]-d->time[i-j]);
	  if((dd<=delta)&&(dd>=delta0)) // within window?? 
	    {
	      nrec++;
	      nrec_condition+=d->x[vindx][i-j];
		

	    }else{

	    if(dd>delta)
	      break;
	  }
	}

      if(nobs==0)
	{
		if(nrec>EPSILON)
		  {
		    if(daop.meansummary==1){
		    ff=(nrec_condition)/(nrec);
		    }else{
		      ff=nrec_condition;
		    }
		  }else{
		  ff=0;
		}
	}else{
	if(nobs==1)
	  ff=nrec;
      }
      return(ff);
  
}



double tsummary_row_wrapper(int k,int row_number,int split_var)
{

  double p;
  int nobs;
  
  nobs=(int) daop.tree_matrix[k][13];
  
  if(daop.method==1||daop.method==2)
    p=tsummary_row_fast(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);

  if(daop.method==3)
    p=tsummary_row_meansum(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var,nobs);
 

  if(daop.method==4)
    p=tsummary_row_tilde(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var,nobs);

    if(daop.method==8)
    p=tsummary_row_tilde_window(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],daop.tree_matrix[k][14],split_var,nobs);

    // double tsummary_row_tilde_window(int rownumber,double cut,double delta,double delta0,int vindx,int nobs)
    
      //ff=tsummary_row_fast(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);

  if(daop.method==5)
    p=tsummary_row_meansum0(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var,nobs);
 

   if(daop.method==6)
     p=tsummary_row_meansum0_window(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],daop.tree_matrix[k][14],split_var,nobs);
  

  if(daop.method==7)
    p=tsummary_row_fast_window(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],daop.tree_matrix[k][14],split_var);
  
  
    return(p);
    }




void predict_trees_fast(int *ntrees,double *pred)
{

  int i,k;

  for(i=0;i<daop.n;i++)
    pred[i]=0;


  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
  for(i=0;i<daop.n;i++)
    {
      pred[i]+=find_tnode_predict_fast(k,i);
    }

    }

}



void predict_trees_gini(int *ntrees,int *ncat,double *pred)
{

  int i,k,m;
  double pr[MAXCAT];

  for(i=0;i<(daop.n*ncat[0]);i++)
    pred[i]=0;


  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
  for(i=0;i<daop.n;i++)
    {
      find_tnode_predict_gini(k,i,pr);
      for(m=0;m<ncat[0];m++)
	pred[i+m*daop.n]+=pr[m];
    }

    }

}





void predict_trees_all(int *ntrees,double *pred)
{

  int i,k;

 
  // find terminal node
  for(k=0;k<ntrees[0];k++)
    {
  for(i=0;i<daop.n;i++)
    {
      pred[i+(k*daop.n)]=find_tnode_predict_fast(k,i);
    }

    }

}








void read_predict(double *x,int *n,int *p,double *time,int *id,int *yindx,double *trees,int *nrow_trees,int *ncol_trees,int *nboost,int *time_split,int *method)
{
  // Read in data needed for predictions 

 
  
  int i,j,ncol,ctree;
  struct data_options *d;
  double *v;
  int ht;
  d=&daop;
  d->LARGE_NUMBER=999999;

  assign_method(method);
  
  ncol=ncol_trees[0];
  d->time_split=time_split[0];
  d->nboost=nboost[0]; // number of trees 
  d->tree_start=(int*) malloc(sizeof(int)*d->nboost);
  d->tree_end=(int*) malloc(sizeof(int)*d->nboost);

  
  // tree matrix
  j=nrow_trees[0];
  d->tree_matrix=(double**) malloc(sizeof(double*)*(j));
  v=(double*) malloc(sizeof(double)*j*ncol);
  d->tree_matrix[0]=v;
  for(i=1;i<j;i++)
    d->tree_matrix[i]=d->tree_matrix[i-1]+ncol;

  d->tree_start[0]=0;
  ctree=0;
  for(i=0;i<nrow_trees[0];i++)
    {

      ht=((int) trees[i*ncol]);
      if(ht!=ctree)
	{
	  d->tree_end[ctree]=i-1;
	  ctree++;
	  d->tree_start[ctree]=i;
	}
      
      for(j=0;j<ncol;j++)
	d->tree_matrix[i][j]=trees[i*ncol+j];
    }

  d->tree_end[d->nboost-1]=nrow_trees[0]-1;
  
  d->row_counter=-1; // this is incremented as rows are added
  d->tree_counter=0;
 
  d->p=p[0];
  d->n=n[0];
  d->x=(double**) malloc(sizeof(double*)*d->p);
  for(i=0;i<d->p;i++)
    {
    d->x[i]=(double*) malloc(sizeof(double)*d->n);
    for(j=0;j<d->n;j++)
      d->x[i][j]=x[(d->n)*i+j];
    }

  d->varimp=(double**) malloc(sizeof(double*)*p[0]);
  v=(double*) malloc(sizeof(double)*d->nboost*p[0]);
  d->varimp[0]=v;
  for(i=1;i<p[0];i++)
    d->varimp[i]=d->varimp[i-1]+d->nboost;
 
  d->time=(double*) malloc(sizeof(double)*d->n);
  d->yindx=yindx[0];
  for(j=0;j<d->n;j++)
    d->time[j]=time[j];

  d->id=(int*) malloc(sizeof(int)*d->n);
  for(j=0;j<d->n;j++)
    d->id[j]=id[j];

  d->predictions=(double*) malloc(sizeof(double)*d->n);
  for(j=0;j<d->n;j++)
    d->predictions[j]=0;

  // node vector
  d->node=(int*) malloc(sizeof(int)*d->n);
  for(i=0;i<d->n;i++)
    d->node[i]=0;
  
  // scratch space
  d->daux1=(double*) malloc(sizeof(double)*d->n);
  d->daux2=(double*) malloc(sizeof(double)*d->n);
  d->daux3=(double*) malloc(sizeof(double)*d->n);
  
  d->iaux1=(int*) malloc(sizeof(int)*d->n);
  d->iaux2=(int*) malloc(sizeof(int)*d->n);
  d->iaux3=(int*) malloc(sizeof(int)*d->n);

}


void free_predict()
{

  int i;
  
  free(daop.daux1);
  free(daop.daux2);
  free(daop.daux3);
  
  free(daop.iaux1);
  free(daop.iaux2);
  free(daop.iaux3);

  free(daop.varimp[0]);
  free(daop.varimp);
  
  free(daop.tree_matrix[0]);
  free(daop.tree_matrix);

  free(daop.node);
  free(daop.id);
  free(daop.predictions);
  free(daop.time);
  free(daop.tree_start);
  free(daop.tree_end);

  
  for(i=0;i<daop.p;i++)
      free(daop.x[i]);

  free(daop.x);


}


void add_bag(int treeindx,int r_number)
{
  struct node *s;
  struct rnumber *r;

 
  HASH_FIND_INT(oob,&treeindx,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=treeindx;
      s->rows=NULL;
      s->nrows=0;
      HASH_ADD_INT(oob,node_number,s);
      //utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&r_number,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=r_number;
      s->nrows++;
      //utarray_push_back(s->rows_array,&r_number); 
    }

 
}





void add_id_row(int id_indx,int r_number)
{
  struct node *s;
  struct rnumber *r;

 
  HASH_FIND_INT(id_rows,&id_indx,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=id_indx;
      s->rows=NULL;
      s->nrows=0;
      HASH_ADD_INT(id_rows,node_number,s);
      //utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&r_number,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=r_number;
      s->nrows++;
      //utarray_push_back(s->rows_array,&r_number); 
    }

 
}



void get_oob(int *oob_vec)
{

  //
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int k;

  for(k=0;k<(daop.nboost*daop.n);k++)
    oob_vec[k]=0;

  k=0;
  HASH_ITER(hh,oob,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      oob_vec[k*daop.n+r->row_number]=1;
    }
	k++;
  }

}




void subsample_oob(int treeindx)
{
  int j;
  // Checked
 GetRNGstate();

 for(j=0;j<daop.n;j++)
   {
     if(unif_rand()<=daop.sample_fraction)
       {	
	 daop.train[j]=1;
	 
       }else{
       daop.train[j]=0;
        add_bag(treeindx,j);
     }
   }


 PutRNGstate();
}


void subsample_id(int treeindx)
{
 
  int j;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  GetRNGstate();

  // for(j=0;j<daop.n;j++)
  HASH_ITER(hh,id_rows,s,stemp)
    {
      if(unif_rand()<=daop.sample_fraction)
	{
	  HASH_ITER(hh,s->rows,r,rtemp){	 
	    daop.train[r->row_number]=1;	  
	  }
	}else{
	HASH_ITER(hh,s->rows,r,rtemp){	 
	  daop.train[r->row_number]=0;
	  add_bag(treeindx,r->row_number);
	}
      
       
      }
    }


  PutRNGstate();


}



void delete_oob()
{

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,oob,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      HASH_DEL(s->rows,r);
      free(r);
    }
    HASH_DEL(oob,s);
    free(s);
  }

}


void delete_id_rows()
{

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,id_rows,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      HASH_DEL(s->rows,r);
      free(r);
    }
    HASH_DEL(id_rows,s);
    free(s);
  }

}



void set_oob(int *oob_vec)
{
  int i,b;
  delete_oob();

  daop.subsample=0; // turn it off 
  for(b=0;b<daop.nboost;b++)
    {
      for(i=0;i<daop.n;i++)
	{
	  if(oob_vec[b*daop.n+i]==1)
	    add_bag(b,i);
	}

    }
  

}


//void add_id_row(int id_indx,int r_number)


void set_id_rows(int *id_vec)
{
  int i;
  delete_id_rows();

      for(i=0;i<daop.n;i++)
	{
	  add_id_row(id_vec[i],i);
	}

}



void varimp_predictor(int j,double *pred)
{

  // Randomly permute predictor 'j' and compute oob predictions for 

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  double *x_orig;
  int *perm_indx;
  double *oobv;
  int i,k,m;
  
  x_orig=daop.daux1;
  perm_indx=daop.iaux1;
  oobv=daop.daux2;

  for(i=0;i<daop.n;i++)
    perm_indx[i]=i;
  
    permute(perm_indx,daop.n);  
 
    //  Rprintf(" vi: 1 \n");
  for(i=0;i<daop.n;i++){
    pred[i]=0;
    oobv[i]=0;
    x_orig[i]=daop.x[j][i];
  }

  //Rprintf(" vi: 2 \n");
  // permute predictor variable ..
  for(i=0;i<daop.n;i++)
    daop.x[j][i]=x_orig[perm_indx[i]];
    
  //Rprintf(" vi: 3 \n");

  // loop through trees and form oob predictions with variable permuted 

  for(k=0;k<daop.nboost;k++)
    {
      // HASH find ...

      HASH_FIND_INT(oob,&k,s);
      HASH_ITER(hh,s->rows,r,rtemp){
	m=r->row_number;
	pred[m]=pred[m]+find_tnode_predict_fast(k,m);
	oobv[m]=oobv[m]+1;
      }

    }
  //Rprintf(" vi: 4 \n");

  for(m=0;m<daop.n;m++)
    {
      if(oobv[m]>.1)
	pred[m]=pred[m]/oobv[m];
    }

  
  //Rprintf(" vi: 5 \n");
  
  // restore permuted predictor
    for(i=0;i<daop.n;i++) 
      daop.x[j][i]=x_orig[i];
    
}


void varimp(int *nperm,double *res)
{
  int j,k,i;
  double *pred;
  double a;

   GetRNGstate();
  pred=(double*) malloc(sizeof(double)*daop.n);

  //  Rprintf(" nperm=%d n=%d p=%d \n",nperm[0],daop.n,daop.p);
  
  //Rprintf(" 1 \n");
  for(i=0;i<(daop.n*daop.p);i++)
    res[i]=0;
  //Rprintf(" 2 \n");
  for(j=0;j<daop.p;j++)
    {
      //      Rprintf(" 3 \n");
      for(k=0;k<nperm[0];k++){
	varimp_predictor(j,pred);
	//Rprintf(" 3 b \n");
	for(i=0;i<daop.n;i++)
	  {
	    res[daop.n*j+i]+=pred[i];
	  }
      }
    }

  //Rprintf(" 4 \n");
  a=(double) nperm[0];
  for(i=0;i<(daop.n*daop.p);i++)
    res[i]=res[i]/a;

  free(pred);
   PutRNGstate();
}






void varimp_predictor_boost(int j,double *pred)
{

  // Randomly permute predictor 'j' and compute oob predictions for 

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  double *x_orig;
  int *perm_indx;
  double *oobv;
  int i,k,m;
  
  x_orig=daop.daux1;
  perm_indx=daop.iaux1;
  oobv=daop.daux2;

  for(i=0;i<daop.n;i++)
    perm_indx[i]=i;
  
    permute(perm_indx,daop.n);  
 
    //  Rprintf(" vi: 1 \n");
  for(i=0;i<daop.n;i++){
    pred[i]=0;
    oobv[i]=0;
    x_orig[i]=daop.x[j][i];
  }

  //Rprintf(" vi: 2 \n");
  // permute predictor variable ..
  for(i=0;i<daop.n;i++)
    daop.x[j][i]=x_orig[perm_indx[i]];
    
  //Rprintf(" vi: 3 \n");

  // loop through trees and form permuted oob predictions 

  for(k=0;k<daop.nboost;k++)
    {
      // HASH find ...

      for(i=0;i<daop.n;i++)
	{
	  pred[i]=pred[i]+find_tnode_predict_fast(k,i);
	}

    }
  //Rprintf(" vi: 4 \n");


  //Rprintf(" vi: 5 \n");
  
  // restore permuted predictor
    for(i=0;i<daop.n;i++) 
      daop.x[j][i]=x_orig[i];
    
}


void varimp_boost(int *nperm,double *res)
{
  int j,k,i;
  double *pred;
  double a;

   GetRNGstate();
  pred=(double*) malloc(sizeof(double)*daop.n);

  //Rprintf(" nperm=%d n=%d p=%d \n",nperm[0],daop.n,daop.p);
  
  //Rprintf(" 1 \n");
  for(i=0;i<(daop.n*daop.p);i++)
    res[i]=0;
  //Rprintf(" 2 \n");
  for(j=0;j<daop.p;j++)
    {
      //Rprintf(" 3 \n");
      for(k=0;k<nperm[0];k++){
	varimp_predictor_boost(j,pred);
	//Rprintf(" 3 b \n");
	for(i=0;i<daop.n;i++)
	  {
	    res[daop.n*j+i]+=pred[i];
	  }
      }
    }

  //Rprintf(" 4 \n");
  a=(double) nperm[0];
  for(i=0;i<(daop.n*daop.p);i++)
    res[i]=res[i]/a;

  free(pred);
   PutRNGstate();
}



void add_row(int n_number,int r_number)
{
  struct node *s;
  struct rnumber *r;
  HASH_FIND_INT(nodes,&n_number,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=n_number;
      HASH_ADD_INT(nodes,node_number,s);
      s->rows=NULL;
      s->nrows=0;
      utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&r_number,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=r_number;
      s->nrows++;
      utarray_push_back(s->rows_array,&r_number); 
    }

}

void delete_node(int n)
{
  struct node *s;
  struct rnumber *r,*rtemp;

  HASH_FIND_INT(nodes,&n,s);
  if(s==NULL){}else{
    HASH_ITER(hh,s->rows,r,rtemp)
      {
	HASH_DEL(s->rows,r);
	free(r);
      }
    utarray_free(s->rows_array);
    
    HASH_DEL(nodes,s);
    free(s);
  }
}

void delete_all_nodes()
{
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,nodes,s,stemp){
      HASH_ITER(hh,s->rows,r,rtemp)
	{
	  HASH_DEL(s->rows,r);
	  free(r);
	}
      utarray_free(s->rows_array);
      HASH_DEL(nodes,s);
      free(s);
    }

}

void delete_node_wrapper(int *n)
{
  delete_node(n[0]);
}

void add_row_wrapper(int *n,int *r)
{
  add_row(n[0],r[0]);
}

void print_node_members(int *n)
{
  int n_number;
  struct node *s;
  struct rnumber *r,*rtemp;
  n_number=n[0];
  HASH_FIND_INT(nodes,&n_number,s);
  if(s==NULL)
    {
      Rprintf(" Node %d does not exist...\n",n[0]);
    }else{
    HASH_ITER(hh,s->rows,r,rtemp){
      Rprintf(" i=%d \n",r->row_number);
    }
  }

}


int split_attempt_sample(int nodeid,int vindx,double *delta,double *tau,struct split_result *sr);
void tsummary_uthash_intermediate(int nodeid,double *nrec,double *y,int *nnode,double cut,double delta,int vindx,int *nmin,int *nmax);
void sample_tau(int nodeid,int vindx,double delta_max,double *tau);
void sample_delta(int nodeid,double *delta);




void split_historical_sample(int nodeid,struct split_result *best_split)
{

  // Sampling based splitting of historical predictors 
  struct split_result sr;
  double *delta,*tau;
  int j,success,init,k;
  double delta_max;
  delta=daop.daux5;
  tau=daop.daux6;
  
  init=0;
  best_split->success=0; // until proven otherwise....

  sample_delta(nodeid,&(delta[0]));
  delta_max=delta[0];
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      k=daop.splitvar_history[j];
      sample_tau(nodeid,k,delta_max,tau); // sample tau 
      //    success=split_attempt_history(nodeid,k,&(delta[0]),&sr);
      success=split_attempt_sample(nodeid,k,delta,tau,&sr);

      if(success==1)
	{
	  if(init==0||(sr.error<(best_split->error)))
	    {
	      
	      best_split[0]=sr;
	      best_split->varindx=k;
	      best_split->success=1;
	      init=1;
	    }
       
	}
    } 
}

//var highscore=0;


int split_attempt_sample(int nodeid,int vindx,double *delta,double *tau,struct split_result *sr)
{
  int k,i;
  int nmin,nmax;
  int nnode;
  double nsplit;
  double nL,nR,sumL,sumR,sumsqL,sumsqR,error,error_ns;
  double a1,a2,a3;
  int init;
  double nmin_d,nmax_d,y;
  init=0;

  
  for(k=0;k<daop.nsamp;k++)
    {
      tsummary_uthash_intermediate(nodeid,daop.daux1,daop.daux2,&nnode,tau[k],delta[k],vindx,&nmin,&nmax);
      nmin_d=(double) nmin;
      nmax_d=(double) nmax;
      //   void tsummary_uthash_intermediate(int nodeid,double *nrec,double *nrec_condition,int *nnode,double cut,double delta,int vindx,int *nmin,int *nmax)
      // sample n in (nmin,nmax)

      a1=(nmin_d)+(nmax_d-nmin_d)*unif_rand();
      nsplit=a1;

      // compute split at n, update
      sumL=0;
      sumR=0;
      sumsqL=0;
      sumsqR=0;
      nL=0;
      nR=0;
      // nrec<n=>Left; nrec>=n=>Right
      for(i=0;i<nnode;i++)
	{
	  if((daop.daux1[i]-nsplit)<(-EPSILON))
	    {
	      nL++;
	      y=daop.daux2[i];
	      sumL+=y;
	      sumsqL+=(y*y);
	      
	    }else{

	      nR++;
	      y=daop.daux2[i];
	      sumR+=y;
	      sumsqR+=(y*y);
	  }

	}

      if(nL>=daop.min_nodesize&&nR>=daop.min_nodesize)
	{
      // computer error
	  a2=(nL+nR);
	  a1=(sumL+sumR)/a2;
	  error_ns=((sumsqL+sumsqR)/a2-a1*a1);
	  a2=sumL/nL;
	  a3=nL/(nR+nL);
	  a1=(sumsqL/nL-(a2*a2))*a3;
	  error=a1;

	  a2=sumR/nR;
	  a3=1-a3;
	  a1=(sumsqR/nR-(a2*a2))*a3;
	  error+=a1;
	  if(init==0||(error<sr->error))
	    {

	      init=1;
	      sr->success=1;
	      sr->varindx=vindx;
	      sr->nL=nL;
	      sr->nR=nR;
	      sr->predL=sumL/nL;
	      sr->predR=sumR/nR;
	      sr->point=nsplit;
	      sr->error=error;
	      sr->error_reduction=error_ns-error;
	      sr->tau=tau[k];
	      sr->delta=delta[k];
	      sr->sumL=sumL;
	      sr->sumR=sumR;
	      sr->sumsqL=sumsqL;
	      sr->sumsqR=sumsqR;

	      

	    }
	}

      

    }

  return(init);
}







void tsummary_uthash(int nodeid,double *nrec,double *nrec_condition,double cut,double delta,int vindx,int *nmin,int *nmax)
{
  // Need use this with a sampling approach ... alternative to full-blown ... 
  
  // For each record, this function counts the number of previous records satisfying x[vindx]>cut
  // for all records whose time falls within  (time[record]-delta,time[record])
  //
  // The total number of records within time interval is written to 'nrec'
  // The number satisfying condition is written to 'nrec_condition'

  int i,j,cid;
  struct node *s;
  struct rnumber *r,*rtemp;
  int min_n,max_n;
  
  struct data_options *d;
  d=&daop;

  min_n=999999;
  max_n=0;

HASH_FIND_INT(nodes,&nodeid,s);

 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      
      nrec[i]=0;
      nrec_condition[i]=0;
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      nrec[i]++;

	      if(d->x[vindx][i-j]>=cut)
		{
		  nrec_condition[i]++;
		  if(nrec_condition[i]<min_n)
		    min_n=nrec_condition[i];
		  if(nrec_condition[i]>max_n)
		    max_n=nrec_condition[i];
		  
		  
		}

	    }else{
	    break;
	  }
	}


    }
 }
 nmin[0]=min_n;
 nmax[0]=max_n;
}





void tsummary_uthash_intermediate(int nodeid,double *nrec,double *y,int *nnode,double cut,double delta,int vindx,int *nmin,int *nmax)
{
  // Need use this with a sampling approach ... alternative to full-blown ... 

  // THIS function is similar to tsummary_uthash() except,
  // indexes are 0:number-in-node (not actual observation indexes)
  
  // For each record, this function counts the number of previous records satisfying x[vindx]>cut
  // for all records whose time falls within  (time[record]-delta,time[record])
  //
  // The total number of records within time interval is written to 'nrec'
  // The number satisfying condition is written to 'nrec_condition'

  int i,j,cid;
  struct node *s;
  struct rnumber *r,*rtemp;
  int min_n,max_n;
  int counter;
  struct data_options *d;
  d=&daop;

  counter=0;
  min_n=999999;
  max_n=0;

HASH_FIND_INT(nodes,&nodeid,s);

 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      
      nrec[counter]=0;
      // ... 
      y[counter]=daop.response[i];//daop.x[daop.yindx][i];
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	     

	      if(d->x[vindx][i-j]>=cut)
		{
		  nrec[counter]++;
		  
		}

	    }else{
	    break;
	  }
	}

      if(nrec[counter]<min_n)
	min_n=nrec[counter];
      if(nrec[counter]>max_n)
	max_n=nrec[counter];
		  

      
      counter++;
    }
 }
 nnode[0]=counter; 
 nmin[0]=min_n;
 nmax[0]=max_n;
}


void sample_tau(int nodeid,int vindx,double delta_max,double *tau)
{
  // Need use this with a sampling approach ... alternative to full-blown ... 
  
  // For each record, this function counts the number of previous records satisfying x[vindx]>cut
  // for all records whose time falls within  (time[record]-delta,time[record])
  //
  // The total number of records within time interval is written to 'nrec'
  // The number satisfying condition is written to 'nrec_condition'

  int i,j,cid,k;
  struct node *s;
  struct rnumber *r,*rtemp;
  double min_tau,max_tau;
  double a;
  
  struct data_options *d;
  d=&daop;

  min_tau=99999999;
  max_tau=-99999999;

HASH_FIND_INT(nodes,&nodeid,s);

 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      
 
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta_max)  // used to be j-i 
	    {

	      a=daop.x[vindx][i-j];
	      if(a<min_tau)
		min_tau=a;
	      if(a>max_tau)
		max_tau=a;

	    }else{
	    break;
	  }
	}


    }
 }


   GetRNGstate();
  // daop.lower_time is the smallest observation time point
  for(k=0;k<daop.nsamp;k++){
    tau[k]=min_tau+(max_tau-min_tau)*unif_rand()+.1;
    //Rprintf(" delta[%d]=%lf lower_time=%lf max_time=%lf \n",k,delta[k],daop.lower_time,max_time);
  }
  PutRNGstate();
 
  //  revsort(delta,daop.iaux1,daop.nsamp);

}





int sample_indx(int n)
{
  // sample an index from 0:(n-1)
  int i;
  double nd,u,uu;
  u=unif_rand();
  nd=((double) n);
  uu=nd*u;
  i=((int) uu);
  // Rprintf(" u=%lf nd=%lf uu=%lf i=%d \n",u,nd,uu,i);

  if(i==n)
    i=n-1;
  
  return(i);
}

void sample_indx_wrapper(int *n,int *res)
{

 GetRNGstate();

 res[0]=sample_indx(n[0]);
 
 PutRNGstate();
  
}


void permute(int *p,int n)
{

  // Take p and permute it
  int k,i,pa;
  
  for(k=0;k<n;k++)
    {
      i=sample_indx((n-k));
      pa=p[(n-k-1)];
      p[(n-k-1)]=p[i];
      p[i]=pa;
    }

}

void permute_wrapper(int *n,int *res)
{
  int *p;
  int k;
  //p=(int*) malloc(sizeof(int)*n[0]);
  for(k=0;k<n[0];k++)
    res[k]=k;
  
  GetRNGstate();
  permute(res,n[0]);
  PutRNGstate();

}

void update_nodevector(int nodeid,int nodeL,int nodeR,int vindx,double cut)
{
  //void tsummary(double *nrec,double *nrec_condition,double cut,double delta,int vindx)
  struct data_options *d;
  int i;
  double ff;
  struct node *s;
  struct rnumber *r,*rtemp;
  
  d=&daop;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){
    Rprintf(" Updating node list, node %d does not exists ...\n",nodeid);
  }else{

    // for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
	  ff=d->x[vindx][i];
	  if(ff<cut)
	    {
	      // left node 
	      add_row(nodeL,i);
	      //d->node[i]=nodeL;
	    }else{
	    // right node 
	    add_row(nodeR,i);
	    //d->node[i]=nodeR;
	  }
	
    }
  

    // Delete node
    delete_node(nodeid);
  }
}





double node_prediction_uthash(int nodeid)
{
  // return node mean
  int i;
  double a;
  double n;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  a=0;
  n=0;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    //for(i=0;i<daop.n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
	  if(daop.train[i]==1)
	    {
	  a+=daop.response[i];
	  n++;
	    }
       
    }
  }
  // Rprintf(" n=%lf \n",n);
  
  if(n<.5)
    n=1;
  a=a/n;

  return(a);
}




void node_prediction_gini(int nodeid,double *a)
{
  // return node mean
  int i,yint;
  //double a;
  double n;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  for(i=0;i<daop.ncat;i++)
    a[i]=0;
  
  n=0;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    //for(i=0;i<daop.n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
	  if(daop.train[i]==1)
	    {
	      yint=(int) daop.response[i];
	      a[yint]++; //=daop.response[i];
	      n++;
	    }
       
    }
  }
  // Rprintf(" n=%lf \n",n);

  
  /* if(n<.5) */
  /*   n=1; */
  /* for(i=0;i<daop.ncat;i++) */
  /*   a[i]=a[i]/n; */

}






void update_nodevector_time(int nodeid,int nodeL,int nodeR,struct split_result *sr)
{
  //void tsummary(double *nrec,double *nrec_condition,double cut,double delta,int vindx)
  struct data_options *d;
  int i;
  double ff;
  struct node *s;
  struct rnumber *r,*rtemp;
  double *x, *nx;
  double cut,cc;
  int vindx;
  
  x=daop.daux6;
  d=&daop;
  nx=daop.daux5;

  cut=sr->point;
  
  vindx=sr->varindx;
  if(daop.method==1||daop.method==2){
    create_history_summary_uthash(nodeid,vindx,sr,x);
  }

  if(daop.method==3){
    //  Rprintf(" create_history_summary_meansum\n");
    create_history_summary_meansum(nodeid,vindx,sr,x);

  }

  if(daop.method==4)
    create_history_summary_tilde(nodeid,vindx,sr,x,nx);

  if(daop.method==8)
    create_history_summary_tilde_window(nodeid,vindx,sr,x,nx);


  if(daop.method==7||daop.method==9)  // window-count and window-count gini 
    create_history_summary_window(nodeid,vindx,sr,x);

  if(daop.method==5)
    create_history_summary_meansum0(nodeid,vindx,sr,x);
 
  if(daop.method==6)
    create_history_summary_meansum0_window(nodeid,vindx,sr,x);
  


  /* Rprintf(" cut=%lf \n",cut); */
  /* cc=0; */

  //for(i=0;i<d->n;i++)
  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    HASH_ITER(hh,s->rows,r,rtemp)
  {

    i=r->row_number;
	  ff=x[i];

	  /* cc++; */
	  /* if(cc<20) */
	  /*   Rprintf("x[%d]=%lf \n",i,x[i]); */
	  
	  if((ff)<cut)
	    {
	      //d->node[i]=nodeL;
	       add_row(nodeL,i);
	    }else{
	    //d->node[i]=nodeR;
	     add_row(nodeR,i);
	  }
	
    }
  
  }


  
}







void update_nodevector_time_meansum(int nodeid,int nodeL,int nodeR,struct split_result *sr)
{
  //void tsummary(double *nrec,double *nrec_condition,double cut,double delta,int vindx)
  struct data_options *d;
  int i;
  double ff;
  struct node *s;
  struct rnumber *r,*rtemp;
  double *x;
  double cut;
  int vindx;
  
  x=daop.daux6;
  d=&daop;


  vindx=sr->varindx;
  create_history_summary_meansum(nodeid,vindx,sr,x);

  cut=sr->point;


  //for(i=0;i<d->n;i++)
  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    HASH_ITER(hh,s->rows,r,rtemp)
  {

    i=r->row_number;
	  ff=x[i];
	  if((ff)<cut)
	    {
	      //d->node[i]=nodeL;
	       add_row(nodeL,i);
	    }else{
	    //d->node[i]=nodeR;
	     add_row(nodeR,i);
	  }
	
    }
  
  }
  
}



void update_nodevector_time_wrapper(int nodeid,int nodeL,int nodeR,struct split_result *sr)
{


  if(daop.method==1||daop.method==2||daop.method>=4) // 4,5,6,7 
    {
      update_nodevector_time(nodeid,nodeL,nodeR,sr);

    }

    if(daop.method==3)
    {
      update_nodevector_time_meansum(nodeid,nodeL,nodeR,sr);
    }

 /* if(daop.method==4) */
 /*    { */
 /*      update_nodevector_time_tilde(nodeid,nodeL,nodeR,sr); */

 /*    } */

    
}

void update_residual_uthash(int tree_indx)
{
  int k,i,start,go,ni,ti,end;
  double node_predictions[10000]; // no more than 100 terminal nodes
  int tnodes[10000];// terminal nodes
  int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[k][0];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[k][2]);
	  //Rprintf(" node id =%d k=%d \n",ni,k);
	  node_predictions[ni]=daop.tree_matrix[k][3];
	  if(daop.tree_matrix[k][4]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	}else{
	go=0;
	break;
      }

    }

  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
      
	  i=r->row_number;
      if(daop.rf!=1)
	{
	  daop.response[i]=(daop.response[i]-node_predictions[ni]);
	  daop.predictions[i]+=(node_predictions[ni]);
	}else{

	if(daop.train[i]==0) // out-of-bag
	  {
	    daop.predictions[i]+=(node_predictions[ni]);
	    daop.oob[i]++;
	  }
      }

      }
      }
    }

 
}


void update_residual_gini(int tree_indx)
{
  int k,i,start,go,ni,ti,end,m;
  double node_predictions[10000][MAXCAT]; // no more than 100 terminal nodes
  int tnodes[10000];// terminal nodes
  int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[k][0];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[k][2]);
	  //Rprintf(" node id =%d k=%d \n",ni,k);
	  for(m=0;m<daop.ncat;m++)
	    node_predictions[ni][m]=daop.tree_matrix[k][15+m];
	  
	  if(daop.tree_matrix[k][4]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	}else{
	go=0;
	break;
      }

    }

  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
      
	  i=r->row_number;
      if(daop.rf!=1)
	{
	  //	  daop.response[i]=(daop.response[i]-node_predictions[ni]);
	  //daop.predictions[i]+=(node_predictions[ni]);
	}else{

	if(daop.train[i]==0) // out-of-bag
	  {
	    for(m=0;m<daop.ncat;m++)
	      daop.predictions_gini[i][m]+=(node_predictions[ni][m]);
	    daop.oob[i]++;
	  }
      }

      }
      }
    }

 
}


/* void variable_importance(int treeindx) */
/* { */
/*   struct node *s; */
/*   struct rnumber *r,*rtemp; */
/*   double h,h_perm,y_orig,y_hat,y_hat_perm; */
/*   int j,i; */
  
/*   HASH_FIND_INT(oob,&treeindx,s); */

/*   // void permute(int *p,int n) */

/*   // double find_tnode_predict_permuted(int tree_indx,int row_number,int permuted_col,int *permuted_row) */
/*   GetRNGstate(); */
/*   if(s==NULL){Rprintf(" Did not find requested oob \n");}else{ */

  

  
/*   for(j=0;j<daop.p;j++) */
/*     { */
/*       daop.varimp[j][treeindx]=0; */
/*       // need permute variable 'j' on the oob sample */
/*       for(i=0;i<daop.n;i++) */
/* 	daop.iaux1[i]=i; */
     
/*       permute(daop.iaux1,daop.n); */
/*       HASH_ITER(hh,s->rows,r,rtemp) */
/* 	{ */
/* 	  i=r->row_number; */
	  
/* 	  y_hat_perm=find_tnode_predict_permuted(treeindx,i,j,daop.iaux1); */
/* 	  y_hat=find_tnode_predict(treeindx,i); */

/* 	  y_orig=daop.x[daop.yindx][i]; */
/* 	  h_perm=y_hat_perm-y_orig; */
/* 	  h=y_hat-y_orig; */
	  
/* 	  //	  daop.varimp[j][treeindx]+=((h_perm*h_perm)-(h*h)); */
/* 	   daop.varimp[j][treeindx]+=(abs(h_perm)-abs(h)); */
	  
/* 	} */
/*       h=(double) s->nrows; */
/*       daop.varimp[j][treeindx]=(daop.varimp[j][treeindx]/h); */
      
      
/*     } */
/*   } */
/*   PutRNGstate(); */
/* } */



void set_train_oob(int treeindx)
{

  // Set daop.train
  int i,b;
  struct node *s;
  struct rnumber *r,*rtemp;
  
  for(i=0;i<daop.n;i++)
    daop.train[i]=1;


  HASH_FIND_INT(oob,&treeindx,s);
  if(s==NULL){Rprintf(" Did not find oob-list element ...\n");}else{

    HASH_ITER(hh,s->rows,r,rtemp)
      {
	daop.train[r->row_number]=0;
      }

  }

  
}

void node_summary_stat(int nodeid,double *sum,double *sumsq,double *n,double *nid)
{
  // Function computes sum of responses in node, sum-sq and total number.

  // nid: records number of historical observations per id.  
  
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  double suma,sumsqa,na,y;
  
HASH_FIND_INT(nodes,&nodeid,s);

 suma=0;
 sumsqa=0;
 na=0;
 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      y=daop.x[daop.yindx][i];
      suma+=y;
      sumsqa+=(y*y);
      na++;
      nid[daop.id[i]]++;
    }  

 }

 sum[0]=suma;
 sumsq[0]=sumsqa;
 n[0]=na;
 
}



void extract_history_uthash(int nodeid,int vindx,double *x,double *time_past,int *id_past,int *npast,double *y,double *time,int *ri_node,double delta,int *n)
{
  
  // This function extracts history of variable vindx for all observations falling in node with 'nodeid'
  // written to:  x,time_past,id_past 
  //
  // node repsonse, time and id: y,time_current,id

  // Function is run at each node splitting, once for each longitudinal predictor 
  // .... would be faster to extract all histories at once.....   
  
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  d=&daop;


HASH_FIND_INT(nodes,&nodeid,s);


 counter=0;
 counter_obs=0;
 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      if(daop.train[i]==1)
	{
 
      ri_node[counter_obs]=i;  // all row-indexes in node 
      counter_obs++;

      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      
	      //y[counter]=d->x[d->yindx][i];
	      y[counter]=d->response[i];
	      time[counter]=d->time[i];
	      x[counter]=d->x[vindx][i-j];
	      id_past[counter]=i; // should record row-number, not  cid;
	      time_past[counter]=d->time[i-j];
	      counter++;
	    }else{
	    break;
	  }
	}

	}
    }
 }

 npast[0]=counter;
 n[0]=counter_obs;
}






void extract_meansum_delta(int nodeid,int vindx,int dindx,double delta,int *n)
{
  
  // This function extracts history of variable vindx for all observations falling in node with 'nodeid'
  // written to:  x,time_past,id_past 
  //
  // node repsonse, time and id: y,time_current,id

  // Function is run at each node splitting, once for each longitudinal predictor 
  // .... would be faster to extract all histories at once.....   
  
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  d=&daop;


HASH_FIND_INT(nodes,&nodeid,s);


 counter=0;
 counter_obs=0;
 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      if(daop.train[i]==1)
	{
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];
      
      cid=d->id[i];

      daop.y_meansum[counter_obs]=d->response[i];
      daop.xmat[dindx][counter_obs]=0;
      daop.nxmat[dindx][counter_obs]=0;
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      
	      //y[counter]=d->x[d->yindx][i];
	      //y[counter_obs]=d->response[i];
	      daop.xmat[dindx][counter_obs]+=d->x[vindx][i-j];
	      daop.nxmat[dindx][counter_obs]+=1;
	    
	    }else{
	    break;
	  }
	    
	}
           counter_obs++;


	}
    }
 }

 
 n[0]=counter_obs;
}




void extract_meansum_window(int nodeid,int vindx,int dindx,double delta,double delta0,int *n)
{
  
  // This function extracts history of variable vindx for all observations falling in node with 'nodeid'
  // written to:  x,time_past,id_past 
  //
  // node repsonse, time and id: y,time_current,id

  // Function is run at each node splitting, once for each longitudinal predictor 
  // .... would be faster to extract all histories at once.....   
  
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  d=&daop;


HASH_FIND_INT(nodes,&nodeid,s);


 counter=0;
 counter_obs=0;
 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      if(daop.train[i]==1)
	{
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];
      
      cid=d->id[i];

      daop.y_meansum[counter_obs]=d->response[i];
      daop.xmat[dindx][counter_obs]=0;
      daop.nxmat[dindx][counter_obs]=0;
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if(((d->time[i]-d->time[i-j])<=delta)&&((d->time[i]-d->time[i-j])>=delta0))  // used to be j-i 
	    {
	      
	      //y[counter]=d->x[d->yindx][i];
	      //y[counter_obs]=d->response[i];
	      daop.xmat[dindx][counter_obs]+=d->x[vindx][i-j];
	      daop.nxmat[dindx][counter_obs]+=1;
	    
	    }else{
	    if((d->time[i]-d->time[i-j])>delta)
	      break;
	  }
	    
	}
           counter_obs++;


	}
    }
 }

 
 n[0]=counter_obs;
}




void check_meansum_summary(int *vi,double *delta,int *ndelta,double *xmat,double *nxmat,int *n)
{

  int nodeid;
  int nobs,j,i,k;

  nodeid=0;

  
  // initialize root node 
  j=0;
  for(i=0;i<daop.n;i++)
    add_row(j,i);  
  
  create_meansum_summary(nodeid,vi[0],delta,ndelta[0],&nobs);

  for(k=0;k<(ndelta[0]*2);k++)
    {
      for(j=0;j<nobs;j++)
	{
	  xmat[k*nobs+j]=daop.xmat[k][j];
	  nxmat[k*nobs+j]=daop.nxmat[k][j];

	}

    }

  Rprintf(" nobs=%d \n",nobs);
  n[0]=nobs;
  

}







void create_meansum_summary(int nodeid,int vindx,double *delta,int ndelta,int *nobs)
{

  // This creates a matrix of summary values for node with id=nodeid and variable 'vindx'. 
  //
  //

  int k,j;
  double *x,*y,*time_current,*time_past;
  int *ri_past,*ri_node;
  double max_delta;
  int n,npast,counter_obs,cid;
  double LARGE_NUMBER;

  LARGE_NUMBER=daop.LARGE_NUMBER;

  //Rprintf(" Entering create_meansum_summary \n");

  max_delta=0;
  for(k=0;k<ndelta;k++)
    {

      extract_meansum_delta(nodeid,vindx,k,delta[k],&n);
     
    }

   nobs[0]=n;


 

  //Rprintf(" create_mean_sum: counter_obs=%d \n",counter_obs);
  // Apply missing together approach, that is columns with missing are duplicated one where missing is replaced by LARGE_NUMER
  // another replaced by -LARGE_NUMBER

  
  for(k=0;k<ndelta;k++)
    {

      daop.column_double[k]=0;
      for(j=0;j<n;j++)
	{

	  if(daop.meansummary==1)
	    {
	      // Mean summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j]/(daop.nxmat[k][j]);
		daop.xmat[k+ndelta][j]=daop.xmat[k][j];
	      }else{
		daop.column_double[k]=1;  // column k has been duplicated. 
		daop.xmat[k][j]=LARGE_NUMBER;
		daop.xmat[k+ndelta][j]=-LARGE_NUMBER;

	      }
	    }else{

	    // Sum summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j];
		daop.xmat[k+ndelta][j]=daop.xmat[k][j];
	      }else{
		daop.column_double[k]=1;
		daop.xmat[k][j]=LARGE_NUMBER;
		daop.xmat[k+ndelta][j]=-LARGE_NUMBER;

	      }


	  }
	}

    }
  


}








void create_meansum0_summary(int nodeid,int vindx,double *delta,int ndelta,int *nobs)
{


  // Mean/sum summary with all missing set to 0 (ie not MISSING-TOGETHER) 
  // This creates a matrix of summary values for node with id=nodeid and variable 'vindx'. 
  //
  //

  int k,j;
  double *x,*y,*time_current,*time_past;
  int *ri_past,*ri_node;
  double max_delta;
  int n,npast,counter_obs,cid;

  //Rprintf(" Entering create_meansum_summary \n");

  max_delta=0;
  for(k=0;k<ndelta;k++)
    {

      extract_meansum_delta(nodeid,vindx,k,delta[k],&n);
     
    }

   nobs[0]=n;


 
  
  for(k=0;k<ndelta;k++)
    {

      daop.column_double[k]=0;
      for(j=0;j<n;j++)
	{

	  if(daop.meansummary==1)
	    {
	      // Mean summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j]/(daop.nxmat[k][j]);
	      }else{

		daop.xmat[k][j]=0;
	      }
	      
	    }else{

	    // Sum summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j];
	      }else{
		daop.xmat[k][j]=0;
	      }


	  }
	}

    }
  


}



void create_meansum0_window(int nodeid,int vindx,double *delta,double *delta0,int ndelta,int *nobs)
{


  // Mean/sum summary with all missing set to 0 (ie not MISSING-TOGETHER) 
  // This creates a matrix of summary values for node with id=nodeid and variable 'vindx'. 
  //
  //

  int k,j;
  double *x,*y,*time_current,*time_past;
  int *ri_past,*ri_node;
  double max_delta;
  int n,npast,counter_obs,cid;

  //Rprintf(" Entering create_meansum_summary \n");

  max_delta=0;
  for(k=0;k<ndelta;k++)
    {

      extract_meansum_window(nodeid,vindx,k,delta[k],delta0[k],&n);
     
    }

   nobs[0]=n;


 
  
  for(k=0;k<ndelta;k++)
    {

      daop.column_double[k]=0;
      for(j=0;j<n;j++)
	{

	  if(daop.meansummary==1)
	    {
	      // Mean summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j]/(daop.nxmat[k][j]);
	      }else{

		daop.xmat[k][j]=0;
	      }
	      
	    }else{

	    // Sum summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j];
	      }else{
		daop.xmat[k][j]=0;
	      }


	  }
	}

    }
  


}




void create_meansum_summary_old(int nodeid,int vindx,double *delta,int ndelta,int *nobs)
{

  // This creates a matrix of summary values for node with id=nodeid and variable 'vindx'. 
  //
  //

  int k,j;
  double *x,*y,*time_current,*time_past;
  int *ri_past,*ri_node;
  double max_delta;
  int n,npast,counter_obs,cid;
  double LARGE_NUMBER;

  LARGE_NUMBER=daop.LARGE_NUMBER;

  //Rprintf(" Entering create_meansum_summary \n");
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  // extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);

  max_delta=0;
  for(k=0;k<ndelta;k++)
    {
      if(delta[k]>max_delta)
	max_delta=delta[k];
    }

  Rprintf(" create_mean_summary: max_delta=%d \n",max_delta);
  
  //Rprintf(" create_meansum_summary: extract_history \n");
  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,max_delta,&n);
  nobs[0]=n;


 
  //Rprintf(" create_mean_sum: loop\n");
  for(k=0;k<ndelta;k++)
    {


      
      cid=ri_past[0];
      counter_obs=0;
      daop.xmat[k][0]=0;
      daop.nxmat[k][0]=0;
      daop.y_meansum[0]=y[0];

 

      for(j=0;j<npast;j++)
	{

	  if(cid!=ri_past[j])
	    {
	      cid=ri_past[j];
	      counter_obs++;
	      daop.xmat[k][counter_obs]=0;
	      daop.nxmat[k][counter_obs]=0;
	      daop.y_meansum[counter_obs]=y[j];
	    }
	  
	  if((time_current[j]-time_past[j])<=delta[k])
	    {
	      daop.xmat[k][counter_obs]+=x[j];
	      daop.nxmat[k][counter_obs]+=1;
	    }
	}

    }
  
  //Rprintf(" create_mean_sum: counter_obs=%d \n",counter_obs);
  // Apply missing together approach, that is columns with missing are duplicated one where missing is replaced by LARGE_NUMER
  // another replaced by -LARGE_NUMBER

  
  for(k=0;k<ndelta;k++)
    {

      daop.column_double[k]=0;
      for(j=0;j<counter_obs;j++)
	{

	  if(daop.meansummary==1)
	    {
	      // Mean summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j]/(daop.nxmat[k][j]);
		daop.xmat[k+ndelta][j]=daop.xmat[k][j];
	      }else{
		daop.column_double[k]=1;  // column k has been duplicated. 
		daop.xmat[k][j]=LARGE_NUMBER;
		daop.xmat[k+ndelta][j]=-LARGE_NUMBER;

	      }
	    }else{

	    // Sum summary 
	      if(daop.nxmat[k][j]>EPSILON){
		daop.xmat[k][j]=daop.xmat[k][j];
		daop.xmat[k+ndelta][j]=daop.xmat[k][j];
	      }else{
		daop.column_double[k]=1;
		daop.xmat[k][j]=LARGE_NUMBER;
		daop.xmat[k+ndelta][j]=-LARGE_NUMBER;

	      }


	  }
	}

    }
  


}





void create_history_summary_uthash(int nodeid,int vindx,struct split_result *sr,double *x)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //
  

  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta;
  d=&daop;


  tau=sr->tau;
  delta=sr->delta;

HASH_FIND_INT(nodes,&nodeid,s);


 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];

      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      if(d->x[vindx][i-j]<tau)
		x[i]++;
	    }else{
	    break;
	  }
	}


    }
 }


}






void create_history_summary_window(int nodeid,int vindx,struct split_result *sr,double *x)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //
  

  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta,delta0;
  d=&daop;


  tau=sr->tau;
  delta=sr->delta;
  delta0=sr->delta0;
  
HASH_FIND_INT(nodes,&nodeid,s);


 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];

      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta&&((d->time[i]-d->time[i-j])>=delta0))  // used to be j-i 
	    {
	      if(d->x[vindx][i-j]<tau)
		x[i]++;
	    }else{
	    if((d->time[i]-d->time[i-j])>delta)
	      break;
	  }
	}


    }
 }


}


void create_history_summary_meansum0_window(int nodeid,int vindx,struct split_result *sr,double *x)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //
  

  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta,delta0,dd;
  int n,nobs;
  d=&daop;

  nobs=sr->missing;
  tau=sr->tau;
  delta=sr->delta;
  delta0=sr->delta0;
  
HASH_FIND_INT(nodes,&nodeid,s);


 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      n=0;
 
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  
	  dd=(d->time[i]-d->time[i-j]);
	  if((dd<=delta)&&(dd>=delta0))  // used to be j-i 
	    {
	      n++;
	      x[i]=x[i]+d->x[vindx][i-j];
	    }else{

	    if((dd)>delta)
	      {
		break;
	      }
	  }
	}
      if(nobs==0)
	{
      if(n>0&&(d->meansummary==1))
	x[i]=x[i]/((double) n);
	}else{
	x[i]=((double) n);
      }

    }
 }


}
void create_history_summary_meansum0(int nodeid,int vindx,struct split_result *sr,double *x)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //
  

  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta,delta0;
  int n,nobs;
  d=&daop;

  nobs=sr->missing; // 0 if split on mean/sum, 1 if split on number of obs

  tau=sr->tau;
  delta=sr->delta;
  delta0=sr->delta0;
  
HASH_FIND_INT(nodes,&nodeid,s);


 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      n=0;
 
      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      n++;
	      x[i]=x[i]+d->x[vindx][i-j];
	    }else{
		break;
	  }
	}

      if(nobs==0)
	{
      if(n>0&&(d->meansummary==1))
	x[i]=x[i]/((double) n);
	}else{

	if(nobs==1)
	  x[i]=((double) n);
      }

    }
 }


}




void create_history_summary_tilde(int nodeid,int vindx,struct split_result *sr,double *x,double *nx)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //


  // FOR EACH row-index associated with the node:
  //       Find fraction of observations (historical) below sr->tau
  

  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta;
  d=&daop;


  tau=sr->tau;
  delta=sr->delta;

HASH_FIND_INT(nodes,&nodeid,s);


 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      nx[i]=0;
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];

      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      nx[i]++;
	      if(d->x[vindx][i-j]<(tau+EPSILON))
		x[i]++;
	    }else{
	    break;
	  }
	}

      if(sr->missing==0)  // indicates whether splitting is on fraction<=tau, or on total number of obs
	{
	  if(nx[i]>EPSILON){
	    x[i]=x[i]/nx[i];
	  }else{
	    x[i]=1.0; // ADDED 2/3/2018 (before missing values were set to 0, at odds with splitting)
	  }

      
      
	}else{
	x[i]=nx[i];

      }

    }
 }


}








void create_history_summary_tilde_window(int nodeid,int vindx,struct split_result *sr,double *x,double *nx)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //


  // FOR EACH row-index associated with the node:
  //       Find fraction of observations (historical) below sr->tau
  

  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta,delta0,dd,cc;
  d=&daop;

  cc=0;

  tau=sr->tau;
  delta=sr->delta;
  delta0=sr->delta0;

HASH_FIND_INT(nodes,&nodeid,s);

 
 //Rprintf(" create_tilde_window; delta=%lf delta0=%lf tau=%lf sr.missing=%d \n",delta,delta0,tau,sr->missing);
 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      nx[i]=0;
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];

      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  
	  dd=(d->time[i]-d->time[i-j]);
	  if((dd<=delta)&&(dd>=(delta0-EPSILON)))  // used to be j-i 
	    {
	      nx[i]++;
	      if(d->x[vindx][i-j]<(tau))
		x[i]++;
	    }else{
	    if(dd>delta)
	      break;
	  }
	}

      if(sr->missing==0)  // indicates whether splitting is on fraction<=tau, or on total number of obs
	{
	  if(nx[i]>EPSILON){
	    x[i]=x[i]/nx[i];
	  }else{
	    x[i]=1.0; // ADDED 2/3/2018 (before missing values were set to 0, at odds with splitting)
	  }

      
      
	}else{
	x[i]=nx[i];

      }
      /* cc++; */
      /* if(cc<20) */
      /* 	Rprintf("cc[%d]=%lf \n",i,x[i]); */
      
    }
 }


}








void create_history_summary_meansum(int nodeid,int vindx,struct split_result *sr,double *x)
{

  // For split 'sr' create history variable for all in node 'nodeid'
  //
  // The history is returned in 'x' such that for response in row 'r' is returned in x[r]
  //
  

  double *nx;
  int i,j,cid,counter;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter_obs;
  struct data_options *d;
  double tau,delta;
  d=&daop;

  nx=daop.daux1;

  
  delta=sr->delta;

HASH_FIND_INT(nodes,&nodeid,s);


 
 if(s==NULL){}else{ 

   //for(i=0;i<d->n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {

      i=r->row_number;
      x[i]=0;
      nx[i]=0;
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];

      cid=d->id[i];
      
      for(j=1;j<d->n;j++)
	{

	  if((i-j)<0)
	    break;

	  if(d->id[i-j]!=cid)
	    break;
	  

	  if((d->time[i]-d->time[i-j])<=delta)  // used to be j-i 
	    {
	      if(sr->missing==0)
		{ //sr->missing==0 means that number of observations is being split on
		  x[i]+=1;
		}else{
		x[i]+=d->x[vindx][i-j];
	      }
		nx[i]++;
	    }else{
	    break;
	  }
	}

	  if(d->meansummary==1)
	    {
	      if(nx[i]>EPSILON)
		{
		  x[i]=x[i]/(nx[i]);
		}
	    }

      if(sr->missing!=0)
	{
	  
	  if(nx[i]<EPSILON)
	    {
	      if(sr->missing>0){
		x[i]=d->LARGE_NUMBER;
	      }else{
		x[i]=-d->LARGE_NUMBER;
	      }

	    }

	}else{
	// split was on number of observations 
	x[i]=nx[i];
      }
    }
 }


}




void initialize_split_result(struct split_result *spr)
{
  int k;
  spr->success=0;
  spr->nL=0;
  spr->nR=0;
  spr->sumL=0;
  spr->sumR=0;
  spr->sumsqL=0;
  spr->sumsqR=0;
  spr->predL=0;
  spr->predR=0;
  spr->point=0;
  spr->error=0;
  spr->delta0=-99;
  spr->missing=0;
  spr->delta=-99;
  spr->tau=-99;
  if(daop.ncat>0)
    {
      // ncat>0 => classification
      for(k=0;k<daop.ncat;k++)
	{
	  spr->sumL_cat[k]=0;
	  spr->sumR_cat[k]=0;
	}

    }
}

void error_split_result(struct split_result *spr)
{

  // Compute error of split result (also predL and predR are formed) 
  spr->predL=spr->sumL/spr->nL;
  spr->predR=spr->sumR/spr->nR;
  spr->error=(spr->nL/(spr->nL+spr->nR))*(spr->sumsqL/spr->nL-spr->predL*spr->predL)+(spr->nR/(spr->nL+spr->nR))*(spr->sumsqR/spr->nR-spr->predR*spr->predR);



}



void error_split_gini(struct split_result *spr)
{

  // Compute error of split result (also predL and predR are formed)
  double eL,eR,p;
  int k;
  eL=0;
  eR=0;
  for(k=0;k<daop.ncat;k++)
    {
      p=spr->sumL_cat[k]/spr->nL;
      eL+=(p*(1-p));
      p=spr->sumR_cat[k]/spr->nR;
      eR+=(p*(1-p));
      

    }
  spr->predL=eL;
  spr->predR=eR;
  spr->error=(spr->nL/(spr->nL+spr->nR))*eL+(spr->nR/(spr->nL+spr->nR))*eR;
  //Rprintf(" nL=%lf nR=%lf eL=%lf eR=%lf \n",spr->nL,spr->nR,eL,eR);


}


void error_reduction_split_result(struct split_result *spr)
{
  // Compute error reduction, this is error prior to node splitting minus error post split 
  double error_prior;

  error_prior=(spr->sumL+spr->sumR)/(spr->nL+spr->nR);
  error_prior=error_prior*error_prior;
  error_prior=(spr->sumsqL+spr->sumsqR)/(spr->nL+spr->nR)-error_prior;
  
  spr->error_reduction=error_prior-spr->error;
      
}

void error_reduction_split_gini(struct split_result *spr)
{
  // Compute error reduction, this is error prior to node splitting minus error post split 
  double error_prior;
  int k;
  double p;

  error_prior=0;
  for(k=0;k<daop.ncat;k++)
    {
      p=(spr->sumL_cat[k]+spr->sumR_cat[k])/(spr->nL+spr->nR);
      error_prior+=(p*(1-p));

    }

  spr->error_reduction=error_prior-spr->error;
      
}


void initialize_sr(int k,int varindx,double *x,double *y,int *xi,int n) 
{
  //  initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);
  
  // initialize daop.sr[k][1...]
  // Return best split in spr_best
  // Function returs 1 if spr_best was identified, otherwise 0
  
  // Assumes 'x' in decreasing order, (represents n[id])
  // y: response, such that y[xi[k]] is the response associated with x[k]

  // initialize_split_result(struct split_result *spr)
  int init,n_max,j,nx,cnx,i;
  double s,s2,yaux,error;
  double epsilon=.000001;
  struct split_result *sp;
  
  // set to zero
  n_max=(int) x[0]; // highest possible value (these are counts, so should work) 
  for(j=0;j<=n_max;j++)
    {
      initialize_split_result(&(daop.sr[k][j]));
    }

  s=0;
  s2=0;
  for(j=0;j<n;j++)
    {
      s+=y[j];
      s2+=(y[j]*y[j]);
    }

  // sr[k][n] contains split info for split where n_ij(delta_k)<n => left, and n_ij>=n => right
    daop.sr[k][n_max].nL=(double) n;
    daop.sr[k][n_max].sumL=s;
    daop.sr[k][n_max].sumsqL=s2;
    daop.sr[k][n_max].varindx=varindx;
  
  // ...
    init=0;
    cnx=(int) x[0];
  for(j=0;j<n;j++)
    {

      nx=(int) x[j];
      if(cnx!=nx)
	{
	  init=1;
	  // new value ..
	  sp=&(daop.sr[k][cnx]);
	  error=sp->nL/(sp->nL+sp->nR)*(sp->sumsqL/sp->nL-sp->sumL*sp->sumL/(sp->nL*sp->nL));
	  error=error+sp->nR/(sp->nL+sp->nR)*(sp->sumsqR/sp->nR-sp->sumR*sp->sumR/(sp->nR*sp->nR));
	  sp->error=error;
	  for(i=nx;i<cnx;i++)
	    daop.sr[k][i]=daop.sr[k][cnx];
	  cnx=nx;
	}

     
      daop.sr[k][nx].nR++;
      yaux=y[xi[j]];
      daop.sr[k][nx].sumR+=yaux;
      daop.sr[k][nx].sumsqR+=(yaux*yaux);
      daop.sr[k][nx].sumsqL-=(yaux*yaux);
      daop.sr[k][nx].sumL-=yaux;
      daop.sr[k][nx].nL--;
      
       
      

    }
  
      
}

void initialize_sr_gini(int k,int varindx,double *x,double *y,int *xi,int n) 
{
  //  initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);
  
  // initialize daop.sr[k][1...]
  // Return best split in spr_best
  // Function returs 1 if spr_best was identified, otherwise 0
  
  // Assumes 'x' in decreasing order, (represents n[id])
  // y: response, such that y[xi[k]] is the response associated with x[k]

  // initialize_split_result(struct split_result *spr)
  int init,n_max,j,nx,cnx,i,yint;
  double s,s2,yaux,error;
  double epsilon=.000001;
  struct split_result *sp;
  
  // set to zero
  n_max=(int) x[0]; // highest possible value (these are counts, so should work) 
  for(j=0;j<=n_max;j++)
    {
      initialize_split_result(&(daop.sr[k][j]));
    }

 
  for(j=0;j<n;j++)
    {
      yint=(int) y[j];
      daop.sr[k][n_max].sumL_cat[yint]++; //+=y[j];
    
    }

  // sr[k][n] contains split info for split where n_ij(delta_k)<n => left, and n_ij>=n => right
    daop.sr[k][n_max].nL=(double) n;
    //daop.sr[k][n_max].sumL=s;
    
    daop.sr[k][n_max].varindx=varindx;
  
  // ...
    init=0;
    cnx=(int) x[0];
  for(j=0;j<n;j++)
    {

      nx=(int) x[j];
      if(cnx!=nx)
	{
	  init=1;
	  // new value ..
	  sp=&(daop.sr[k][cnx]);
	  error_split_gini(sp);
	  for(i=nx;i<cnx;i++)
	    daop.sr[k][i]=daop.sr[k][cnx];
	  cnx=nx;
	}

     
      daop.sr[k][nx].nR++;
      yint=(int) y[xi[j]];
      daop.sr[k][nx].sumR_cat[yint]++;//=yaux;
      //   daop.sr[k][nx].sumsqR+=(yaux*yaux);
      //daop.sr[k][nx].sumsqL-=(yaux*yaux);
      daop.sr[k][nx].sumL_cat[yint]--;
      daop.sr[k][nx].nL--;
      
       
      

    }
  
      
}





void initialize_sr_tilde(int varindx,double *y,int n,int n_max) 
{
  //  initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);
  
  // initialize daop.sr[k][1...]

  // sr[0][j] contains information on split if split had been made at j/n_max ( f_ij(c)<j/n_max=> left; f_ij(c)>=j/n_max=> right)  

  

  // y: response, such that y[xi[k]] is the response associated with x[k]

 
  int init,j,nx,cnx,i;
  double s,s2,yaux,error;
  double epsilon=.000001;
  struct split_result *sp;
  int k;

  k=0;  // left-over from 'sr_initialize' not used in this function. 
  
  // set to zero
  //  n_max=(int) x[0]; // highest possible value (these are counts, so should work) 
  for(j=0;j<=n_max;j++)
    {
      initialize_split_result(&(daop.sr[k][j]));
    }

  s=0;
  s2=0;
  for(j=0;j<n;j++)
    {
      s+=y[j];
      s2+=(y[j]*y[j]);
    }

  // Prior to looping through the sorted historical values (from largest to smallest), all observations have fraction value
  // equal to 1, so all sr contain the same information.

  
  
  for(j=0;j<=n_max;j++)
    {
    daop.sr[k][j].nR=(double) n;
    daop.sr[k][j].sumR=s;
    daop.sr[k][j].sumsqR=s2;
    daop.sr[k][j].varindx=varindx;
    }
 
}



void initialize_sr_tilde_k(int k,int varindx,double *y,int n,int n_max) 
{
  //  initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);
  
 
  int init,j,nx,cnx,i;
  double s,s2,yaux,error;
  double epsilon=.000001;
  struct split_result *sp;

  // set to zero
  //  n_max=(int) x[0]; // highest possible value (these are counts, so should work) 
  for(j=0;j<=n_max;j++)
    {
      initialize_split_result(&(daop.sr[k][j]));
    }

  s=0;
  s2=0;
  for(j=0;j<n;j++)
    {
      s+=y[j];
      s2+=(y[j]*y[j]);
    }

  
  
  for(j=0;j<=n_max;j++)
    {
    daop.sr[k][j].nR=(double) n;
    daop.sr[k][j].sumR=s;
    daop.sr[k][j].sumsqR=s2;
    daop.sr[k][j].varindx=varindx;
    }
 
}




void basic_split(double *x,double *y,int n,int *xi,int minnodesize,struct split_result *s)
{
  // basic splitting. given predictor x, response y, find best split point, return in structure split_result

  // Splitting is defined as:
  // x<best_point -> go left
  // x>=best_point -> go right 
  
  double sumL,sumR,nL,nR,sumsqL,sumsqR,predL,predR,split_error;
  int k;
  double min_error,best_nL,best_nR,best_predL,best_predR,best_point,best_ssqL,best_ssqR,best_sumL,best_sumR;
  double xc,yaux,mnsd;
  int init;
  struct split_result csplit,bsplit;

 initialize_split_result(&csplit);
  
  mnsd=(double) minnodesize;
  
  for(k=0;k<n;k++){
    xi[k]=k;
    csplit.nL++;
    csplit.sumL+=y[k];
    csplit.sumsqL+=(y[k]*y[k]);
  }
  
  revsort(x,xi,n);  // sorts descending order

  xc=x[0];
  csplit.point=xc;
  init=0;
  for(k=0;k<n;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{

	  if(csplit.nL>=mnsd&&csplit.nR>=mnsd)
	    {
	      // valid split
	      error_split_result(&csplit); // compute error at split point 	
	      if(init==0||(csplit.error<bsplit.error))
		{
		  bsplit=csplit;
		  init=1;
		}
	    }
	  
	}
      csplit.nR++;
      csplit.nL--;
      yaux=y[xi[k]];
      csplit.sumL-=yaux;
      csplit.sumR+=yaux;
      csplit.sumsqR+=(yaux*yaux);
      csplit.sumsqL-=(yaux*yaux);
      xc=x[k];
      csplit.point=xc;
      
      if(csplit.nL<mnsd)
	break;
    }


  if(init!=0)
    {
      s[0]=bsplit;
      s->success=1;
      
    }else{
    s->success=0;
    }
  
}





void basic_split_gini(double *x,double *y,int n,int *xi,int minnodesize,struct split_result *s)
{
  // basic splitting. given predictor x, response y, find best split point, return in structure split_result

  // Splitting is defined as:
  // x<best_point -> go left
  // x>=best_point -> go right 
  
  double sumL,sumR,nL,nR,sumsqL,sumsqR,predL,predR,split_error;
  int k;
  double min_error,best_nL,best_nR,best_predL,best_predR,best_point,best_ssqL,best_ssqR,best_sumL,best_sumR;
  double xc,yaux,mnsd;
  int init;
  int yint;
  struct split_result csplit,bsplit;

  initialize_split_result(&csplit);  // should work for classification
  
  mnsd=(double) minnodesize;
  
  for(k=0;k<n;k++){
    xi[k]=k;
    csplit.nL++;
    yint=(int) y[k];
    csplit.sumL_cat[yint]++; //+=y[k];
    //csplit.sumsqL+=(y[k]*y[k]);
  }
  
  revsort(x,xi,n);  // sorts descending order

  xc=x[0];
  csplit.point=xc;
  init=0;
  for(k=0;k<n;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{

	  if(csplit.nL>=mnsd&&csplit.nR>=mnsd)
	    {
	      // valid split
	      error_split_gini(&csplit); // compute error at split point
	      //Rprintf(" error_split=%lf x[%d]=%lf best-error=%lf \n",csplit.error,k,x[k],bsplit.error);
	      if(init==0||(csplit.error<bsplit.error))
		{
		  bsplit=csplit;
		  init=1;
		}
	    }
	  
	}
      csplit.nR++;
      csplit.nL--;
      yint=(int) y[xi[k]];
      csplit.sumL_cat[yint]--; //=yaux;
      csplit.sumR_cat[yint]++; //=yaux;
      //      csplit.sumsqR+=(yaux*yaux);
      //      csplit.sumsqL-=(yaux*yaux);
      xc=x[k];
      csplit.point=xc;
      
      if(csplit.nL<mnsd)
	break;
    }


  if(init!=0)
    {
      s[0]=bsplit;
      s->success=1;
      
    }else{
    s->success=0;
    }
  
}





int split_attempt_history(int nodeid,int vindx,double *delta,struct split_result *return_split)
{
  // Finds best split of node 'nodeid' by historical splitting by variable 'vindx', across 'ndelta' values of delta
  // Results returned to split_xx variables

  // NOTE: delta assumed to be in decreasing order
  
  int npast,n;
  double *x,*time_past,*y,*time_current;
  int *id_past,*id_current;
  double sum,sumsq,nd,dd;
  int *ni; // record number of historical obs per subject. 
  double *nid; // same as ni, but in as double 
  int *ri_node,*ri_past;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;
  int j,k,cni,init,i,n_i;
  double xc,error,yaux;
  double mnsd;
  double dn_i;
  int indx;
  int record_counter;
  int ncc;
  
  return_split->success=0;
  init=0;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
  //	if(0){ 
     if(npast>0)
      {
	for(k=0;k<npast;k++)
	  daop.iaux3[k]=k;


	revsort(x,daop.iaux3,npast);
 
	// Initialize count matrix for each delta 
	for(k=0;k<n;k++) // for each row in node, set counter to zero for each 'delta'
	  {
	    for(j=0;j<daop.nsamp;j++)
	      daop.n_row_number[j][ri_node[k]]=0;
	  }


	// ... record number of historical observations (for each delta) associated with each row-index in node
	for(k=0;k<npast;k++)
	  {
	    
	    indx=daop.iaux3[k];
	    cni=ri_past[indx];
	    dd=(time_current[indx]-time_past[indx]);
      
	    for(j=0;j<daop.nsamp;j++){
	      
	      if(delta[j]>=dd)
		daop.n_row_number[j][cni]++;
	      
	    }
	  }



  // for each value of delta attempt split 
  init=0;
  for(k=0;k<daop.nsamp;k++)
    {
      for(j=0;j<n;j++)
	{
	  daop.daux5[j]=daop.n_row_number[k][ri_node[j]];
	  daop.daux6[j]=daop.response[ri_node[j]];
	  daop.iaux5[j]=j;
	}


      basic_split(daop.daux5,daop.daux6,n,daop.iaux5,daop.min_nodesize,&sresult);

      sresult.tau=3*(x[0]+100); // some large value that all values are less than .... x[0] is max value in node.
      sresult.varindx=vindx;

      initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);

      
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      
	      best_split=sresult;
	      best_split.delta=(delta[0]+100);  // 
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
  

    }


  
  if(1){
  // loop through historical values
  xc=x[0];

  mnsd=(double) daop.min_nodesize;
  
  for(k=0;k<daop.nsamp;k++)
    daop.counter_changed[k]=0;

  
  record_counter=0;
  for(k=0;k<npast;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{
	  /* if(unif_rand()<.01) */
	  /*   { */
	  daop.ncc++;
	  // new split point, evaluate previous one 
	  for(j=0;j<daop.nsamp;j++)
	    {
	      if(daop.counter_changed[j]>0)
		{
		  for(i=0;i<daop.counter_changed[j];i++)
		    {
		      //daop.ncc++;
		      n_i=daop.n_changed[j][i];
		      sp=&daop.sr[j][n_i];
		      if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
			{
			  
			  error_split_result(sp);
			  sp->tau=xc;
			  sp->delta=delta[j];
			  sp->point=n_i;
			  if(init==0||(sp->error<best_split.error))
			    {
			      //Rprintf(" 3\n");
			      daop.sr[j][n_i].tau=xc;
			      daop.sr[j][n_i].error=sp->error;
			      best_split=daop.sr[j][n_i];
			      best_split.delta=delta[j];
			      best_split.point=(double) n_i;
			      init=1;
			  
			    }

			}
		      daop.counter_changed[j]=0;

		    }

		}

	  
	      //} // this is for if(unif_rand()<.01)... testing 
	    }
	}

      xc=x[k];
      dd=(time_current[daop.iaux3[k]]-time_past[daop.iaux3[k]]);
      for(j=0;j<daop.nsamp;j++)
	{

	  if(dd<=delta[j])
	    {

	      // daop.ncc++;
	      dn_i= daop.n_row_number[j][ri_past[daop.iaux3[k]]]; //
	      n_i=(int) dn_i;
	      yaux=y[daop.iaux3[k]];
	      daop.sr[j][n_i].nR--;
	      daop.sr[j][n_i].nL++;
	      daop.sr[j][n_i].sumL+=yaux;
	      daop.sr[j][n_i].sumR-=yaux;
	      daop.sr[j][n_i].sumsqL+=(yaux*yaux);
	      daop.sr[j][n_i].sumsqR-=(yaux*yaux);
	      
	      daop.n_row_number[j][ri_past[daop.iaux3[k]]]=((dn_i-1));

	      daop.n_changed[j][daop.counter_changed[j]]=n_i;
	      daop.counter_changed[j]++;
	    }
	}
 	
    }
  //Rprintf("Total loop count=%d\n",ncc);
  

      
  } // matches if(0){.. 
  return_split[0]=best_split;

 
      }
     return init;

}





int split_attempt_window(int nodeid,int vindx,double *delta,double *delta0,struct split_result *return_split)
{
  // Finds best split of node 'nodeid' by historical splitting by variable 'vindx', across 'ndelta' values of delta
  // Results returned to split_xx variables

  // NOTE: delta assumed to be in decreasing order
  
  int npast,n;
  double *x,*time_past,*y,*time_current;
  int *id_past,*id_current;
  double sum,sumsq,nd,dd;
  int *ni; // record number of historical obs per subject. 
  double *nid; // same as ni, but in as double 
  int *ri_node,*ri_past;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;
  int j,k,cni,init,i,n_i;
  double xc,error,yaux;
  double mnsd;
  double dn_i;
  int indx;
  int record_counter;
  int ncc;
  
  return_split->success=0;
  init=0;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
  //	if(0){ 
     if(npast>0)
      {
	for(k=0;k<npast;k++)
	  daop.iaux3[k]=k;


	revsort(x,daop.iaux3,npast);
 
	// Initialize count matrix for each delta 
	for(k=0;k<n;k++) // for each row in node, set counter to zero for each 'delta'
	  {
	    for(j=0;j<daop.nsamp_window;j++)
	      daop.n_row_number[j][ri_node[k]]=0;
	  }


	// ... record number of historical observations (for each delta) associated with each row-index in node
	for(k=0;k<npast;k++)
	  {
	    
	    indx=daop.iaux3[k];
	    cni=ri_past[indx];
	    dd=(time_current[indx]-time_past[indx]);
            
	    for(j=0;j<daop.nsamp_window;j++){
	      
	      if((delta[j]>=dd)&&(delta0[j]<=dd))  // within window
		daop.n_row_number[j][cni]++;
	      
	    }
	  }



  // for each value of delta attempt split 
  init=0;
  for(k=0;k<daop.nsamp_window;k++)
    {
      for(j=0;j<n;j++)
	{
	  daop.daux5[j]=daop.n_row_number[k][ri_node[j]];
	  daop.daux6[j]=daop.response[ri_node[j]];
	  daop.iaux5[j]=j;
	}

      //if(0) // TURNING OFF 
      basic_split(daop.daux5,daop.daux6,n,daop.iaux5,daop.min_nodesize,&sresult);

      sresult.tau=3*(x[0]+100); // some large value that all values are less than .... x[0] is max value in node.
      sresult.varindx=vindx;

       initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);

      
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      // delta updating here cant be right, or?? 
	      best_split=sresult;
	      best_split.delta=delta[k]; //(delta[0]+100);
	      best_split.delta0=delta0[k]; //EPSILON;
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
  

    }


  
  if(1){
  // loop through historical values
  xc=x[0];

  mnsd=(double) daop.min_nodesize;
  
  for(k=0;k<daop.nsamp_window;k++)
    daop.counter_changed[k]=0;

  
  record_counter=0;
  for(k=0;k<npast;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{
	  /* if(unif_rand()<.01) */
	  /*   { */
	  daop.ncc++;
	  // new split point, evaluate previous one 
	  for(j=0;j<daop.nsamp_window;j++)
	    {
	      if(daop.counter_changed[j]>0)
		{
		  for(i=0;i<daop.counter_changed[j];i++)
		    {
		      //daop.ncc++;
		      n_i=daop.n_changed[j][i];
		      sp=&daop.sr[j][n_i];
		      if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
			{
			  
			  error_split_result(sp);
			  sp->tau=xc;
			  sp->delta=delta[j];
			  sp->delta0=delta0[j];
			  sp->point=n_i;

			  if(init==0||(sp->error<best_split.error))
			    {
			      //Rprintf(" 3\n");
			      daop.sr[j][n_i].tau=xc;
			      daop.sr[j][n_i].error=sp->error;
			      best_split=daop.sr[j][n_i];
			      best_split.delta=delta[j];
			      best_split.delta0=delta0[j];
			      best_split.point=(double) n_i;
			      init=1;
			  
			    }

			}
		      daop.counter_changed[j]=0;

		    }

		}

	  
	      //} // this is for if(unif_rand()<.01)... testing 
	    }
	}

      xc=x[k];
      dd=(time_current[daop.iaux3[k]]-time_past[daop.iaux3[k]]);
      for(j=0;j<daop.nsamp_window;j++)
	{

	  if((dd<=delta[j])&&(dd>=delta0[j]))  // within window j 
	    {

	      // daop.ncc++;
	      dn_i= daop.n_row_number[j][ri_past[daop.iaux3[k]]]; //
	      n_i=(int) dn_i;
	      yaux=y[daop.iaux3[k]];
	      daop.sr[j][n_i].nR--;
	      daop.sr[j][n_i].nL++;
	      daop.sr[j][n_i].sumL+=yaux;
	      daop.sr[j][n_i].sumR-=yaux;
	      daop.sr[j][n_i].sumsqL+=(yaux*yaux);
	      daop.sr[j][n_i].sumsqR-=(yaux*yaux);
	      
	      daop.n_row_number[j][ri_past[daop.iaux3[k]]]=((dn_i-1));

	      daop.n_changed[j][daop.counter_changed[j]]=n_i;
	      daop.counter_changed[j]++;
	    }
	}
 	
    }
  //Rprintf("Total loop count=%d\n",ncc);
  

      
  } // matches if(0){.. 
  return_split[0]=best_split;

 
      }
     return init;

}







int split_attempt_window_gini(int nodeid,int vindx,double *delta,double *delta0,struct split_result *return_split)
{
  // Finds best split of node 'nodeid' by historical splitting by variable 'vindx', across 'ndelta' values of delta
  // Results returned to split_xx variables

  // NOTE: delta assumed to be in decreasing order
  
  int npast,n;
  double *x,*time_past,*y,*time_current;
  int *id_past,*id_current;
  double sum,sumsq,nd,dd;
  int *ni; // record number of historical obs per subject. 
  double *nid; // same as ni, but in as double 
  int *ri_node,*ri_past;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;
  int j,k,cni,init,i,n_i;
  double xc,error,yaux;
  double mnsd;
  double dn_i;
  int indx;
  int record_counter;
  int ncc;
  int yint;
  
  return_split->success=0;
  init=0;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  //  Rprintf(" split_attempt_window_gini: extract_history_uthash\n");
  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
  //	if(0){
  // Rprintf(" split_attempt_window_gini: done extract_history \n");
  
     if(npast>0)
      {
	for(k=0;k<npast;k++)
	  daop.iaux3[k]=k;


	revsort(x,daop.iaux3,npast);
 
	// Initialize count matrix for each delta 
	for(k=0;k<n;k++) // for each row in node, set counter to zero for each 'delta'
	  {
	    for(j=0;j<daop.nsamp_window;j++)
	      daop.n_row_number[j][ri_node[k]]=0;
	  }


	// ... record number of historical observations (for each delta) associated with each row-index in node
	for(k=0;k<npast;k++)
	  {
	    
	    indx=daop.iaux3[k];
	    cni=ri_past[indx];
	    dd=(time_current[indx]-time_past[indx]);
            
	    for(j=0;j<daop.nsamp_window;j++){
	      
	      if((delta[j]>=dd)&&(delta0[j]<=dd))  // within window
		daop.n_row_number[j][cni]++;
	      
	    }
	  }

	// Rprintf(" split_attempt_window_gini: n_row_number \n");
	


  // for each value of delta attempt split 
  init=0;
  for(k=0;k<daop.nsamp_window;k++)
    {
      for(j=0;j<n;j++)
	{
	  daop.daux5[j]=daop.n_row_number[k][ri_node[j]];
	  daop.daux6[j]=daop.response[ri_node[j]];
	  daop.iaux5[j]=j;
	}

       basic_split_gini(daop.daux5,daop.daux6,n,daop.iaux5,daop.min_nodesize,&sresult);
 
      sresult.tau=3*(x[0]+100); // some large value that all values are less than .... x[0] is max value in node.
      sresult.varindx=vindx;
      
      initialize_sr_gini(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);
      
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      // delta updating here cant be right, or?? 
	      best_split=sresult;
	      best_split.delta=delta[k]; //(delta[0]+100);
	      best_split.delta0=delta0[k]; //EPSILON;
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
  

    }

  //	} // matches if(0) ... turning off every thing but extract and sort. 
     // Print best split result ...
  /*  Rprintf(" ####   Print best-split result for delta=%lf ######################## \n",delta[k]);
      print_split_result(&best_split);
  */

  
  if(1){
  // loop through historical values
  xc=x[0];

  mnsd=(double) daop.min_nodesize;
  
  for(k=0;k<daop.nsamp_window;k++)
    daop.counter_changed[k]=0;

  //Rprintf("split_attempt_window_gini: loop; \n");
  
  record_counter=0;
  for(k=0;k<npast;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{
	  /* if(unif_rand()<.01) */
	  /*   { */
	  daop.ncc++;
	  // new split point, evaluate previous one 
	  for(j=0;j<daop.nsamp_window;j++)
	    {
	      if(daop.counter_changed[j]>0)
		{
		  for(i=0;i<daop.counter_changed[j];i++)
		    {
		      //daop.ncc++;
		      n_i=daop.n_changed[j][i];
		      sp=&daop.sr[j][n_i];
		      if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
			{
			  
			  error_split_gini(sp);
			  // ########### START: TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			  sp->tau=xc;
			  sp->delta=delta[j];
			  sp->delta0=delta0[j];
			  sp->point=n_i;
			  /* if(record_counter<1000){ */
			  /*   record_split(sp,record_counter); */
			  /*   record_counter++; */
			  /* } */
			  // ########### END: TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			  if(init==0||(sp->error<best_split.error))
			    {
			      //Rprintf(" 3\n");
			      daop.sr[j][n_i].tau=xc;
			      daop.sr[j][n_i].error=sp->error;
			      best_split=daop.sr[j][n_i];
			      best_split.delta=delta[j];
			      best_split.delta0=delta0[j];
			      best_split.point=(double) n_i;
			      init=1;
			  
			    }

			}
		      daop.counter_changed[j]=0;

		    }

		}

	  
	      //} // this is for if(unif_rand()<.01)... testing 
	    }
	}

      xc=x[k];
      dd=(time_current[daop.iaux3[k]]-time_past[daop.iaux3[k]]);
      for(j=0;j<daop.nsamp_window;j++)
	{

	  if((dd<=delta[j])&&(dd>=delta0[j]))  // within window j 
	    {

	      // daop.ncc++;
	      dn_i= daop.n_row_number[j][ri_past[daop.iaux3[k]]]; //
	      n_i=(int) dn_i;
	      yaux=y[daop.iaux3[k]];
	      yint=(int) yaux;
	      daop.sr[j][n_i].nR--;
	      daop.sr[j][n_i].nL++;
	      daop.sr[j][n_i].sumL_cat[yint]++;//=yaux;
	      daop.sr[j][n_i].sumR_cat[yint]--;// -=yaux;
	      //daop.sr[j][n_i].sumsqL+=(yaux*yaux);
	      //daop.sr[j][n_i].sumsqR-=(yaux*yaux);
	      
	      daop.n_row_number[j][ri_past[daop.iaux3[k]]]=((dn_i-1));

	      daop.n_changed[j][daop.counter_changed[j]]=n_i;
	      daop.counter_changed[j]++;
	    }
	}
 	
    }
  //Rprintf("Total loop count=%d\n",ncc);
  

      
  } // matches if(0){.. 
  return_split[0]=best_split;

 
      }
     return init;

}















int split_history_tilde_aux(int nodeid,int vindx,double *delta,struct split_result *return_split)
{

  // This splits on n(theta;delta)/n(
  // fills up daop.sr[j][k] giving split info for splitting on n_ih(theta_j)<k => left, and n_ih>=k => right
  // Finds best split of node 'nodeid' by historical splitting by variable 'vindx', across 'ndelta' values of delta
  // Results returned to split_xx variables

  // NOTE: delta assumed to be in decreasing order
  
  int npast,n;
  double *x,*time_past,*y,*time_current;
  int *id_past,*id_current;
  double sum,sumsq,nd,dd;
  int *ni; // record number of historical obs per subject. 
  double *nid; // same as ni, but in as double 
  int *ri_node,*ri_past;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;
  int j,k,cni,init,i,n_i;
  double xc,error,yaux;
  double mnsd;
  double dn_i;
  int indx;
  int record_counter;
  int *lpos; // used to keep track of where
  int n_max,n_changed_counter,curr_indx,prev_indx;
  double dn_max,cv,pv;
  int *n_changed_vec;
  int kk,ncc;


  return_split->success=0;
  init=0;


  // n_max sets the grid of fractions (0:n_max)/n_max
  n_max=5;  // this should be determined otherwise, but for now... 
  dn_max=(double) n_max;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  lpos=daop.iaux6;
  n_changed_vec=daop.iaux2;

  
  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
  if(npast>0)
    {
      for(k=0;k<npast;k++)
	daop.iaux3[k]=k;


      revsort(x,daop.iaux3,npast);
 
      // Initialize count matrix for each delta 
      for(k=0;k<n;k++) // for each row in node, set counter to zero for each 'delta'
	{
	  //for(j=0;j<daop.nsamp;j++)
	  j=0;
	  daop.n_row_number[j][ri_node[k]]=0;
	}


      // ... record number of historical observations (for each delta) associated with each row-index in node
      for(k=0;k<npast;k++)
	{
	    
	  indx=daop.iaux3[k];
	  cni=ri_past[indx];
	  dd=(time_current[indx]-time_past[indx]);
      
	  j=0;
	  if(delta[j]>=dd)
	    daop.n_row_number[j][cni]++;
	      
	    
	}



      // attempt split on number of observations 
      init=0;

      k=0;
      for(j=0;j<n;j++)
	{
	  daop.daux5[j]=daop.n_row_number[k][ri_node[j]];
	  daop.daux6[j]=daop.response[ri_node[j]];
	  daop.iaux5[j]=j;

	}

  
      basic_split(daop.daux5,daop.daux6,n,daop.iaux5,daop.min_nodesize,&sresult);
      
      sresult.tau=3*(x[0]+100); // some large value that all values are less than .... x[0] is max value in node.
      sresult.varindx=vindx;
      sresult.missing=1; // 1 for n_obs splitting, else 0 
      sresult.delta=delta[0];
	
 	
      initialize_sr_tilde(vindx,daop.daux6,n,n_max);


          k=0;
      for(j=0;j<n;j++)
	{
	  daop.daux5[ri_node[j]]=daop.n_row_number[k][ri_node[j]];
	}

      for(j=0;j<n;j++)
	lpos[ri_node[j]]=daop.daux5[ri_node[j]];

     

   
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      
	      best_split=sresult;
	      best_split.delta=delta[0]; //(delta[0]+100);
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
      

      

      // Splitting is left-node: f_ij(c;delta)< f  right-node  f_ij(c;delta)>= f

      // f_{ij}(c;delta) is fraction of observations for i-th subject in (t_ij-delta,t_ij) that are less-than-equal-to 'c'

      // Current implementation only considers fractions 1:n_max/n_max (n_max= max observed number of observations in interval)

      // 
  
      if(1){ // TESTING 
	// loop through historical values
	xc=x[0];

	mnsd=(double) daop.min_nodesize;
  

	
	n_changed_counter=0;
	ncc=0;
	for(k=0;k<npast;k++)
	  {

	    if(fabs(x[k]-xc)>EPSILON)
	      {

		// New potential split point => evaluate 
	  
		// new split point, evaluate previous one
		j=0; // left over from copy (here only single delta is considered)
		daop.ncc++;

		for(i=0;i<n_changed_counter;i++)
		  {
		    //daop.ncc++;
		    sp=&daop.sr[0][n_changed_vec[i]];
		    if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
		      {
			  
			error_split_result(sp);
			sp->tau=xc;
			sp->delta=delta[0];
			sp->point=((double) n_changed_vec[i])/((dn_max));
			sp->missing=0;
			if(init==0||(sp->error<best_split.error))
			  {
			  
			    daop.sr[0][n_changed_vec[i]].tau=xc;
			    daop.sr[0][n_changed_vec[i]].error=sp->error;
			    best_split=daop.sr[0][n_changed_vec[i]];
			    best_split.delta=delta[0];
			    best_split.point=((double) n_changed_vec[i])/((dn_max));
			    init=1;
			  
			  }

		      }
		  }
		xc=x[k];
		n_changed_counter=0;
	      }

	    // for each obs equal to x[k] perform change to sr,
	    // record which elements of sr have been changed
	    // when new x[k] encountered, evaluate all changed sr-elements
	    // also set n-changed counter to 0


 
	    // previous value
	    //Rprintf(" daop.iaux3[%d]=%d n_past=%d delta=%lf daux5=%lf \n",k,daop.iaux3[k],npast,delta[0],daop.daux5[daop.iaux3[k]]);

	    kk=ri_past[daop.iaux3[k]];
	    if(daop.daux5[kk]>EPSILON)
	      {
		//Rprintf("k=%d\n",k);
	    pv=(lpos[kk]/daop.daux5[kk]);
	    lpos[kk]=lpos[kk]-1;
	    cv=((lpos[kk])/daop.daux5[kk]);
	    //Rprintf(" DD: pv=%lf cv=%lf \n",pv,cv);
	    if(lpos[kk]<0)
	      {
		Rprintf("negative:  pv=%lf cv=%lf delta=%lf k=%d npast=%d\n",pv,cv,delta[0],k,npast);
		Rprintf(" lp=%lf daop.daux5=%lf daop.iaux3[k]=%d \n",lpos[kk],daop.daux5[kk],daop.iaux3[k]);
		break;

	      }
	    
	    prev_indx=(int) (pv*(dn_max));  // smallest grid value for which k was in right node
	    curr_indx=(int) (cv*dn_max); // smallest grid value for which k IS in right node
	    //Rprintf(" DD: prev_indx=%d curr_indx=%d delta=%lf npast=%d daux5=%lf\n",prev_indx,curr_indx,delta[0],npast,daop.daux5[daop.iaux3[k]]);

	    // From (but not including) curr_indx to prev_indx, k needs to be moved left ... 

	    yaux=y[daop.iaux3[k]];
	    for(j=(curr_indx+1);j<=prev_indx;j++)
	      {
		//daop.ncc++;
	    n_changed_vec[n_changed_counter]=j;
	    n_changed_counter++;
	     
	    daop.sr[0][j].nR--;
	    daop.sr[0][j].nL++;
	    daop.sr[0][j].sumL+=yaux;
	    daop.sr[0][j].sumR-=yaux;
	    daop.sr[0][j].sumsqL+=(yaux*yaux);
	    daop.sr[0][j].sumsqR-=(yaux*yaux);      

	  }
      

	      }
	  }

	//Rprintf(" Total loop count=%d \n",ncc);
  

      
      } // matches if(0){..

      //Rprintf("DEBUG: writing best_split to return_split (error=%lf)\n",best_split.error);
      //print_split_result(&best_split);
      
	    return_split[0]=best_split;

 
    }
	    return init;
}


 









int split_history_tildeNEW_aux(int nodeid,int vindx,double *delta,struct split_result *return_split)
{

  // RECODING of "split_history_tilde_aux"
  //
  // 

  
  // This splits on n(theta;delta)/n(
  // fills up daop.sr[j][k] giving split info for splitting on n_ih(theta_j)<k => left, and n_ih>=k => right
  // Finds best split of node 'nodeid' by historical splitting by variable 'vindx', across 'ndelta' values of delta
  // Results returned to split_xx variables

  // NOTE: delta assumed to be in decreasing order
  
  int npast,n;
  double *x,*time_past,*y,*time_current;
  int *id_past,*id_current;
  double sum,sumsq,nd,dd;
  int *ni; // record number of historical obs per subject. 
  double *nid; // same as ni, but in as double 
  int *ri_node,*ri_past;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;
  int j,k,cni,init,i,n_i;
  double xc,error,yaux;
  double mnsd;
  double dn_i;
  int indx;
  int record_counter;
  int *lpos; // used to keep track of where
  int n_max,n_changed_counter,curr_indx,prev_indx;
  double dn_max,cv,pv;
  int *n_changed_vec;
  int kk,ncc,current_rindx;
  double ff,ff_split,ff_prior;
  int counter;
  
  // Rprintf("Entering split_history_tilde_aux\n");
  return_split->success=0;
  init=0;


 


  
  // n_max sets the grid of fractions (0:n_max)/n_max
  n_max=daop.n_max;  // this should be determined otherwise, but for now... 
  dn_max=(double) n_max;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  lpos=daop.iaux6;
  n_changed_vec=daop.iaux2;


  // void extract_history_uthash(int nodeid,int vindx,double *x,double *time_past,int *id_past,int *npast,double *y,double *time,int *ri_node,double delta,int *n)
  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
  //   Rprintf("split_history_tilde_aux: finished  extract_history_uthash\n");
  if(npast>0)
    {
      for(k=0;k<npast;k++)
	daop.iaux3[k]=k;


      revsort(x,daop.iaux3,npast);
 
      // Initialize count matrix for each delta 
      for(k=0;k<n;k++) // for each row in node, set counter to zero for each 'delta'
	{
	  //for(j=0;j<daop.nsamp;j++)
	  j=0;
	  daop.n_row_number[j][ri_node[k]]=0;
	}


      // ... record number of historical observations (for each delta) associated with each row-index in node
      for(k=0;k<npast;k++)
	{
	    
	  indx=daop.iaux3[k]; 
	  cni=ri_past[indx];
	  dd=(time_current[indx]-time_past[indx]);
      
	  j=0;
	  if(delta[j]>=dd)
	    daop.n_row_number[j][cni]++;
	      
	    
	}



      // attempt split on number of observations 
      init=0;

      k=0;
      for(j=0;j<n;j++)
	{
	  daop.daux5[j]=daop.n_row_number[k][ri_node[j]];
	  daop.daux6[j]=daop.response[ri_node[j]];
	  daop.iaux5[j]=j;

	  //if(j<10)
	  //Rprintf(" %lf %lf  %d \n",daop.daux5[j],daop.daux6[j],daop.iaux5[j]);
	}

 
      basic_split(daop.daux5,daop.daux6,n,daop.iaux5,daop.min_nodesize,&sresult);
     
      sresult.tau=3*(x[0]+100);
      sresult.varindx=vindx;
      sresult.missing=1; // 1 for n_obs splitting, else 0 
      sresult.delta=delta[0];
	
      // Print result q
      //print_split_result(&sresult);
 	
      initialize_sr_tilde(vindx,daop.daux6,n,n_max);

      //  for(k=0;k<=n_max;k++)
      //	print_split_result(&daop.sr[0][k]);


      // The decrementor 
      k=0;
      for(j=0;j<n;j++)
	{
	  daop.daux5[ri_node[j]]=daop.n_row_number[k][ri_node[j]];
	}

  
      // 
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      
	      best_split=sresult;
	      best_split.delta=delta[0]; //(delta[0]+100);
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
  

  
 
  
      if(1){ // TESTING 
	// loop through historical values
	xc=x[0];

	mnsd=(double) daop.min_nodesize;
  

	
	n_changed_counter=0;
	ncc=0;

	counter=0;
	for(k=0;k<npast;k++)
	  {

	    current_rindx=ri_past[daop.iaux3[k]];
	    
	    daop.daux5[current_rindx]--;
	    if(daop.daux5[current_rindx]<(-EPSILON))
	      Rprintf(" Error, decrementor negative: daux[%d]=%lf \n",current_rindx,daop.daux5[current_rindx]);
	    
	    // fraction of obs less than or equal to x[k]
	    if(daop.n_row_number[0][current_rindx]>EPSILON)
	      {
		// For observations with no history, these are here defined to have ff=1 for all tau... 
		
		// Fraction of observations < x[k]
	    ff=daop.daux5[current_rindx]/daop.n_row_number[0][current_rindx];
	    ff_prior=ff+1/daop.n_row_number[0][current_rindx];

	    if(1){
	    for(j=1;j<=n_max;j++)
	      {
		ff_split=(((double) j)/((double) n_max));

		sp=&daop.sr[0][j];
		if((ff<(ff_split))&&(ff_prior>=(ff_split)))  // Added EPSILON .... 
		  {
		    
		    yaux=y[daop.iaux3[k]];
		    sp->nR--;
		    sp->nL++;
		    sp->sumL+=yaux;
		    sp->sumR-=yaux;
		    sp->sumsqL+=(yaux*yaux);
		    sp->sumsqR-=(yaux*yaux);
		    //		    sp->error=10;  // SOMETHING IS UP 
		    error_split_result(sp);
		    sp->tau=x[k];
		    sp->delta=delta[0];
		    sp->point=ff_split;
		    sp->missing=0;
		    
		  }
	      }
	    }
	      
	      }  
	    if(fabs(x[k]-xc)>EPSILON)
	      {


		for(j=1;j<=n_max;j++)
		  {
		    
		    sp=&daop.sr[0][j];

		    
		    if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
		      {
			if(init==0||(sp->error<best_split.error))
			  {
			  
			    best_split=daop.sr[0][j];
			    init=1;
			    //print_split_result(&best_split);
			  }
		    
		      }
		  }

		xc=x[k];
	      }
	  }


      
      } // matches if(1){..

      
	    return_split[0]=best_split;
	    //print_split_result(&best_split);
 
    }
	    return init;
}











int split_history_tilde_aux_window(int nodeid,int vindx,double *delta,double *delta0,struct split_result *return_split)
{

  // RECODING of "split_history_tilde_aux"
  //
  // 

  
  // This splits on n(theta;delta)/n(
  // fills up daop.sr[j][k] giving split info for splitting on n_ih(theta_j)<k => left, and n_ih>=k => right
  // Finds best split of node 'nodeid' by historical splitting by variable 'vindx', across 'ndelta' values of delta
  // Results returned to split_xx variables

  // NOTE: delta assumed to be in decreasing order
  
  int npast,n;
  double *x,*time_past,*y,*time_current;
  int *id_past,*id_current;
  double sum,sumsq,nd,dd;
  int *ni; // record number of historical obs per subject. 
  double *nid; // same as ni, but in as double 
  int *ri_node,*ri_past;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;
  int j,k,cni,init,i,n_i;
  double xc,error,yaux;
  double mnsd;
  double dn_i;
  int indx;
  int record_counter;
  int *lpos; // used to keep track of where
  int n_max,n_changed_counter,curr_indx,prev_indx;
  double dn_max,cv,pv;
  int *n_changed_vec;
  int kk,ncc,current_rindx;
  double ff,ff_split,ff_prior,ddw;
  int counter;
  
  // Rprintf("Entering split_history_tilde_aux\n");
  return_split->success=0;
  init=0;


 


  
  // n_max sets the grid of fractions (0:n_max)/n_max
  n_max=daop.n_max;  // this should be determined otherwise, but for now... 
  dn_max=(double) n_max;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  //  lpos=daop.iaux6;
  // n_changed_vec=daop.iaux2;


  // void extract_history_uthash(int nodeid,int vindx,double *x,double *time_past,int *id_past,int *npast,double *y,double *time,int *ri_node,double delta,int *n)
  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
  //   Rprintf("split_history_tilde_aux: finished  extract_history_uthash\n");
  if(npast>0)
    {
      for(k=0;k<npast;k++)
	daop.iaux3[k]=k;


      revsort(x,daop.iaux3,npast);
 
      // Initialize count matrix for each delta 
      for(k=0;k<n;k++) // for each row in node, set counter to zero for each 'delta'
	{
	  //for(j=0;j<daop.nsamp;j++)
	  j=0;
	  daop.n_row_number[j][ri_node[k]]=0;
	}


      // ... record number of historical observations (for each delta) associated with each row-index in node
      for(k=0;k<npast;k++)
	{
	    
	  indx=daop.iaux3[k]; 
	  cni=ri_past[indx];
	  dd=(time_current[indx]-time_past[indx]);
      
	  j=0;
	  if((delta[j]>=dd)&&(dd>=delta0[j]))
	    daop.n_row_number[j][cni]++;
	      
	    
	}



      // attempt split on number of observations 
      init=0;

      k=0;
      for(j=0;j<n;j++)
	{
	  daop.daux5[j]=daop.n_row_number[k][ri_node[j]];
	  daop.daux6[j]=daop.response[ri_node[j]];
	  daop.iaux5[j]=j;

	  //if(j<10)
	  //Rprintf(" %lf %lf  %d \n",daop.daux5[j],daop.daux6[j],daop.iaux5[j]);
	}

 
      basic_split(daop.daux5,daop.daux6,n,daop.iaux5,daop.min_nodesize,&sresult);
     
      sresult.tau=3*(x[0]+100);
      sresult.varindx=vindx;
      sresult.missing=1; // 1 for n_obs splitting, else 0 
      sresult.delta=delta[0];
      sresult.delta0=delta0[0];
	
      // Print result q
      //print_split_result(&sresult);
 	
      initialize_sr_tilde(vindx,daop.daux6,n,n_max);

      //  for(k=0;k<=n_max;k++)
      //	print_split_result(&daop.sr[0][k]);


      // The decrementor 
      k=0;
      for(j=0;j<n;j++)
	{
	  daop.daux5[ri_node[j]]=daop.n_row_number[k][ri_node[j]];
	}

  
      // 
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      
	      best_split=sresult;
	      best_split.delta=delta[0]; //(delta[0]+100);
	      best_split.delta0=delta0[0]; //(delta[0]+100);
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
  

  
 
  
      if(1){ // TESTING 
	// loop through historical values
	xc=x[0];

	mnsd=(double) daop.min_nodesize;
  

	
	n_changed_counter=0;
	ncc=0;

	counter=0;
	for(k=0;k<npast;k++)
	  {
	    // frac(x_h<x) 

	    // increment accumulator for row-index associated with k

	    // up-date sum/sumsq if fraction-transition takes place, otherwise not

	    // For ease of imp, do separately for each fraction.

	    // RECORD OF CHANGES

	    // loop to n_max changed start from  0 to 1

	    current_rindx=ri_past[daop.iaux3[k]];

	    ddw=(time_current[daop.iaux3[k]]-time_past[daop.iaux3[k]]);
	    if((ddw<=delta[0])&&(ddw>=delta0[0]))
	      {
	      daop.daux5[current_rindx]--;
	    
	    if(daop.daux5[current_rindx]<(-EPSILON))
	      Rprintf(" Error, decrementor negative: daux[%d]=%lf \n",current_rindx,daop.daux5[current_rindx]);
	    
	    // fraction of obs less than or equal to x[k]
	    if(daop.n_row_number[0][current_rindx]>EPSILON)
	      {
		// For observations with no history, these are here defined to have ff=1 for all tau... 
		
		// Fraction of observations < x[k]
	    ff=daop.daux5[current_rindx]/daop.n_row_number[0][current_rindx];
	    ff_prior=ff+1/daop.n_row_number[0][current_rindx];

	    if(1){
	    for(j=1;j<=n_max;j++)
	      {
		ff_split=(((double) j)/((double) n_max));

		sp=&daop.sr[0][j];
		if((ff<(ff_split))&&(ff_prior>=(ff_split)))  // Added EPSILON .... 
		  {
		    
		    yaux=y[daop.iaux3[k]];
		    sp->nR--;
		    sp->nL++;
		    sp->sumL+=yaux;
		    sp->sumR-=yaux;
		    sp->sumsqL+=(yaux*yaux);
		    sp->sumsqR-=(yaux*yaux);
		    //		    sp->error=10;  // SOMETHING IS UP 
		    error_split_result(sp);
		    sp->tau=x[k];
		    sp->delta=delta[0];
		    sp->delta0=delta0[0];
		    sp->point=ff_split;
		    sp->missing=0;
		    
		  }
	      }
	    }
	      
	      }  
	    if(fabs(x[k]-xc)>EPSILON)
	      {


		for(j=1;j<=n_max;j++)
		  {
		    
		    sp=&daop.sr[0][j];

		    
		    if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
		      {
			if(init==0||(sp->error<best_split.error))
			  {
			  
			    best_split=daop.sr[0][j];
			    init=1;
			    //print_split_result(&best_split);
			  }
		    
		      }
		  }

		xc=x[k];
	      }
	      }
	  }

      
      } // matches if(1){..

      
	    return_split[0]=best_split;
	    // print_split_result(&best_split);
 
    }
	    return init;
}











int split_history_tilde(int nodeid,int vindx,double *delta,int ndelta,struct split_result *return_split)
{
  int k,init,success;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;

  // Rprintf("Entering split_history_tilde\n");
  init=0;

  if(1){ // TESTING
   
    for(k=0;k<ndelta;k++)
    {
    success=split_history_tildeNEW_aux(nodeid,vindx,&(delta[k]),&sresult);
 

    
      if(success==1)
	{
	  //Rprintf("DEBUG: successful split error=%lf \n",sresult.error);

	  if(init==0)
	    {
	      
	      best_split=sresult;
	      init=1;
	      //Rprintf(" DEBUG: init was 0, now best_split.error=%lf \n",best_split.error);
	    }else{
	    
	    if(sresult.error<best_split.error)
	      {
		best_split=sresult;
	        //Rprintf(" DEBUG: init was 1, now best_split.error=%lf \n",best_split.error);

	      }
	  }


	}
      }

  }

  if(init==1)
    return_split[0]=best_split;
  
  return(init);
  
}





int split_history_tilde_window(int nodeid,int vindx,double *delta,double *delta0,struct split_result *return_split)
{
  int k,init,success;
  struct split_result sresult;
  struct split_result best_split;
  struct split_result *sp;

  // Rprintf("Entering split_history_tilde\n");
  init=0;

  if(1){ // TESTING 
  for(k=0;k<daop.nsamp_window;k++)
    {
      success=split_history_tilde_aux_window(nodeid,vindx,&(delta[k]),&(delta0[k]),&sresult);

      if(success==1)
	{
	  //Rprintf("DEBUG: successful split error=%lf \n",sresult.error);

	  if(init==0)
	    {
	      
	      best_split=sresult;
	      init=1;
	      //Rprintf(" DEBUG: init was 0, now best_split.error=%lf \n",best_split.error);
	    }else{
	    
	    if(sresult.error<best_split.error)
	      {
		best_split=sresult;
	        //Rprintf(" DEBUG: init was 1, now best_split.error=%lf \n",best_split.error);

	      }
	  }


	}
    }

  }

  if(init==1)
    return_split[0]=best_split;

  //print_split_result(&best_split);
  return(init);
  
}






void split_historical_tilde(int nodeid,struct split_result *best_split)
{
  // split node in using standardized counts

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  int success;


  // Rprintf("Entering: split_historical_tilde\n");
  init=0;
  best_split->success=0; // until proven otherwise....

  sample_delta(nodeid,&(delta[0]));
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      k=daop.splitvar_history[j];
     
      //      success=split_attempt_history(nodeid,k,&(delta[0]),&sr);
      success=split_history_tilde(nodeid,k,&(delta[0]),daop.nsamp,&sr);
      if(success==1)
	{
	  if(init==0||(sr.error<(best_split->error)))
	    {
	      //Rprintf("sr.error=%lf best_split->error=%lf \n",sr.error,best_split->error);
	      best_split[0]=sr;
	      best_split->varindx=k;
	      best_split->success=1;
	      init=1;
	    }
       
	}
    } 
}





void split_window_tilde(int nodeid,struct split_result *best_split)
{
  // split node in using standardized counts

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  double delta0[1000];
  int success;


  // Rprintf("Entering: split_historical_tilde\n");
  init=0;
  best_split->success=0; // until proven otherwise....
  delta_window(nodeid,&(delta[0]),&(delta0[0]));

 
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      k=daop.splitvar_history[j];
     
      success=split_history_tilde_window(nodeid,k,&(delta[0]),&(delta0[0]),&sr);
      if(success==1)
	{
	  if(init==0||(sr.error<(best_split->error)))
	    {
	      //Rprintf("sr.error=%lf best_split->error=%lf \n",sr.error,best_split->error);
	      best_split[0]=sr;
	      best_split->varindx=k;
	      best_split->success=1;
	      init=1;
	    }
       
	}
    } 
}



















void split_standard(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  int predictors[1000];

  
  x=daop.daux5;
  y=daop.daux4;
  xi=daop.iaux5;

  init=0;
  best_split->success=0; // until proven otherwise....

  if(daop.mtry<daop.p){
  counter=0;
  for(k=0;k<daop.p;k++)
    {
      if(k!=daop.yindx)
	{
	  predictors[counter]=k;
	  counter++;
	}
    }

  permute(predictors,(daop.p-1));
  }
  
  for(j=0;j<daop.mtry;j++)
    {
      if(daop.mtry==daop.p){
	k=j;
      }else{
	k=predictors[j];
      }
      
      if(k!=daop.yindx)
	{
	  // get data
	  HASH_FIND_INT(nodes,&nodeid,s);

	  counter=0;
	  if(s==NULL){}else{ 

	    HASH_ITER(hh,s->rows,r,rtemp)
	      {
		i=r->row_number;
		if(daop.train[i]==1){
		  //y[counter]=daop.x[daop.yindx][i];
		  y[counter]=daop.response[i];
		  x[counter]=daop.x[k][i];
		  counter++;
		}
	      }
	  }
	  
	  if(counter>0)
	    {
	      basic_split(x,y,counter,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=k;
		      init=1;
		    }
		  
		}
	      

	    }
	}
      
    }

}









void split_standard_gini(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  int predictors[1000];
  int yint;
  
  x=daop.daux5;
  y=daop.daux4;
  xi=daop.iaux5;

  init=0;
  best_split->success=0; // until proven otherwise....

  if(daop.mtry<daop.p){
  counter=0;
  for(k=0;k<daop.p;k++)
    {
      if(k!=daop.yindx)
	{
	  predictors[counter]=k;
	  counter++;
	}
    }

  permute(predictors,(daop.p-1));
  }
  
  for(j=0;j<daop.mtry;j++)
    {
      if(daop.mtry==daop.p){
	k=j;
      }else{
	k=predictors[j];
      }
      
      if(k!=daop.yindx)
	{
	  // get data
	  HASH_FIND_INT(nodes,&nodeid,s);

	  counter=0;
	  if(s==NULL){}else{ 

	    HASH_ITER(hh,s->rows,r,rtemp)
	      {
		i=r->row_number;
		if(daop.train[i]==1){
		  //y[counter]=daop.x[daop.yindx][i];
		  y[counter]=daop.response[i];
		  x[counter]=daop.x[k][i];
		  counter++;
		}
	      }
	  }
	  
	  if(counter>0)
	    {
	      basic_split_gini(x,y,counter,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=k;
		      init=1;
		    }
		  
		}
	      

	    }
	}
      
    }

}






// ... basic split test ...

void split_meansum(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  
  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  int success,vi,nobs;

  //Rprintf("Entering split_meansum \n");
  // WORKING ON THIS FUNCTION 10/11/2017
  init=0;
  best_split->success=0; // until proven otherwise....

  xi=daop.iaux1;
  sample_delta(nodeid,&(delta[0]));
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      vi=daop.splitvar_history[j];

      create_meansum_summary(nodeid,vi,delta,daop.nsamp,&nobs);
      // this will fill in daop.xmat[][]
      //Rprintf(" nobs=%d \n",nobs);
      for(k=0;k<daop.nsamp;k++)
	{
	  basic_split(daop.xmat[k],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=vi;
		      best_split->delta=delta[k];
		      best_split->missing=1;
		      init=1;
		    }
		  
		}
	      // WORKING HERE 
	      if(daop.column_double[k]==1)
		{
		  // there were missing values, split using "missing together"
		  basic_split(daop.xmat[k+daop.nsamp],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
		  if(sr.success==1)
		    {
		      if(init==0||((best_split->error)>sr.error))
			{
			  best_split[0]=sr;
			  best_split->varindx=vi;
			  best_split->delta=delta[k];
			  best_split->missing=-1;
			  init=1;
			}
		  
		    }    
	  
		}


	      if(0)
		{
	      // Attempt split on nxmat[k], i.e the number of observations between observation time and delta_k
	      
	      basic_split(daop.nxmat[k],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=vi;
		      best_split->delta=delta[k];
		      best_split->missing=0;
		      init=1;
		    }
		  
		}

		}

	}
    }
   

}





// ... basic split test ...

void split_meansum0(int nodeid,struct split_result *best_split)
{
  // split function for meansum where missing values set to zero (not MISSING-TOGETHER)

  
  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  int success,vi,nobs;


  init=0;
  best_split->success=0; // until proven otherwise....

  xi=daop.iaux1;
  sample_delta(nodeid,&(delta[0]));
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      vi=daop.splitvar_history[j];

      create_meansum0_summary(nodeid,vi,delta,daop.nsamp,&nobs);


      for(k=0;k<daop.nsamp;k++)
	{
	  basic_split(daop.xmat[k],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=vi;
		      best_split->delta=delta[k];
		      best_split->missing=0;  // ok?
		      init=1;
		    }
		  
		}



	      if(1){
		// turning this off .. better to add a column of 1's to data-matrix (split only once) 
	      // Attempt split on nxmat[k], i.e the number of observations between observation time and delta_k
		// Nope to the above. wont work without restructuring 

		
	      basic_split(daop.nxmat[k],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=vi;
		      best_split->delta=delta[k];
		      best_split->missing=1;  // CHANGED 0 to 1 (3/3/2018)
		      init=1;
		    }
		  
		}

	      }
	      

	}
    }
   

}





// ... basic split test ...

void split_meansum0_window(int nodeid,struct split_result *best_split)
{
  // split function for meansum where missing values set to zero (not MISSING-TOGETHER)

  
  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  double delta0[1000];
  int success,vi,nobs;


  init=0;
  best_split->success=0; // until proven otherwise....

  xi=daop.iaux1;
  //sample_delta(nodeid,&(delta[0]));
  delta_window(nodeid,&(delta[0]),&(delta0[0]));

  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      vi=daop.splitvar_history[j];

      //create_meansum0_summary(nodeid,vi,delta,daop.nsamp,&nobs);
      create_meansum0_window(nodeid,vi,delta,delta0,daop.nsamp_window,&nobs);


      for(k=0;k<daop.nsamp_window;k++)
	{
	  basic_split(daop.xmat[k],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=vi;
		      best_split->delta=delta[k];
		      best_split->delta0=delta0[k];
		      best_split->missing=0;  // ok?
		      init=1;
		    }
		  
		}




	      // Attempt split on nxmat[k], i.e the number of observations between observation time and delta_k
	      if(1){
		// TURNING THIS off.
		// (add column of 1's to data matrix instead)
	      basic_split(daop.nxmat[k],daop.y_meansum,nobs,xi,daop.min_nodesize,&sr);
	      
	      if(sr.success==1)
		{
		  if(init==0||((best_split->error)>sr.error))
		    {
		      best_split[0]=sr;
		      best_split->varindx=vi;
		      best_split->delta=delta[k];
		      best_split->delta0=delta0[k];
		      best_split->missing=1;  // changed 0 to 1 (3/3/2018)
		      init=1;
		    }
		  
		}
	      }
	      

	}
    }
   

}












void sample_delta(int nodeid,double *delta)
{

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int init,i,k;
  double max_time,a;


  if(daop.set_delta==0)
    {
 
  init=0;
  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    //for(i=0;i<daop.n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
	  if(daop.train[i]==1)
	    {
	      a=daop.time[i];
	      if(init==0||(max_time<a))
		{
		  max_time=a;
		  init=1;
		}
	    }
       
    }
  }

  
  GetRNGstate();
  // daop.lower_time is the smallest observation time point
  for(k=0;k<daop.nsamp;k++){
    delta[k]=daop.lower_time+(max_time-daop.lower_time)*unif_rand()+.1;
    //Rprintf(" delta[%d]=%lf lower_time=%lf max_time=%lf \n",k,delta[k],daop.lower_time,max_time);
  }
  PutRNGstate();

    }else{

    //
    sampleWOR(daop.set_delta,daop.nsamp,daop.iaux1);
    for(k=0;k<daop.nsamp;k++)
      delta[k]=daop.delta_unique[daop.iaux1[k]];
    

  }
  
  revsort(delta,daop.iaux1,daop.nsamp);
  
}





// ... basic split test ...

void split_historical(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  int success;

  
  init=0;
  best_split->success=0; // until proven otherwise....

  sample_delta(nodeid,&(delta[0]));
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      k=daop.splitvar_history[j];
     
      success=split_attempt_history(nodeid,k,&(delta[0]),&sr);

      if(success==1)
	{
	  if(init==0||(sr.error<(best_split->error)))
	    {
	      
	      best_split[0]=sr;
	      best_split->varindx=k;
	      best_split->success=1;
	      init=1;
	    }
       
	}
    } 
}


void delta_window(int nodeid,double *delta,double *delta0)
{

  double delta_new[1000];
  int counter;
  int k,j;

  // sample delta 
   sample_delta(nodeid,&(delta[0]));

  if(daop.window_delta==1)
    {
      
      // All possible according to delta, that is [delta[j],delta[k]] k=0,..,nsamp-1 and j=k,..,nsamp-1
      counter=0;
      for(k=0;k<daop.nsamp;k++)
	{
	  for(j=k;j<daop.nsamp;j++)
	    {
	      delta_new[counter]=delta[k]+EPSILON;
	      delta0[counter]=delta[j]-EPSILON;
	      counter++;
	    }

	}
      
      for(k=0;k<counter;k++)
	{
	  delta[k]=delta_new[k];
	}

      daop.nsamp_window=counter;

      //for(k=0;k<counter;k++)
      //Rprintf("delta[%d]=%lf delta0=%lf \n",k,delta[k],delta0[k]);
      
    }


  if(daop.window_delta==2)
    {

      for(k=0;k<daop.nsamp;k++)
	delta0[k]=3*EPSILON;

      daop.nsamp_window=daop.nsamp;


    }


  if(daop.window_delta==3)
    {




    }


}

void split_window(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  double delta0[1000];
  int success;

  
  init=0;
  best_split->success=0; // until proven otherwise....

  //sample_delta(nodeid,&(delta[0]));
  delta_window(nodeid,&(delta[0]),&(delta0[0]));
  // ... (delta,delta0) sampling functions with different methods
  //     1. sample delta, and single delta0
  //     2. exhaustive search
  //     3. sample delta, sample delta0
  
  // sample delta0 
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      k=daop.splitvar_history[j];
     
      success=split_attempt_window(nodeid,k,&(delta[0]),&(delta0[0]),&sr);

      if(success==1)
	{
	  if(init==0||(sr.error<(best_split->error)))
	    {
	      
	      best_split[0]=sr;
	      best_split->varindx=k;
	      best_split->success=1;
	      init=1;
	    }
       
	}
    } 
}


void split_window_gini(int nodeid,struct split_result *best_split)
{
  // split node in standard manner

  struct split_result sr;
  double *x,*y;
  int *xi;
  struct node *s;
  struct rnumber *r,*rtemp;
  int counter;
  int init,k,i,j;
  double delta[1000];
  double delta0[1000];
  int success;

  
  init=0;
  best_split->success=0; // until proven otherwise....

  //sample_delta(nodeid,&(delta[0]));
  delta_window(nodeid,&(delta[0]),&(delta0[0]));
  // ... (delta,delta0) sampling functions with different methods
  //     1. sample delta, and single delta0
  //     2. exhaustive search
  //     3. sample delta, sample delta0
  
  // sample delta0 
  
  init=0;
  for(j=0;j<daop.nvar_history;j++)
    {
      k=daop.splitvar_history[j];

      //  Rprintf(" split_attempt_window_gini \n");
      success=split_attempt_window_gini(nodeid,k,&(delta[0]),&(delta0[0]),&sr);

      if(success==1)
	{
	  if(init==0||(sr.error<(best_split->error)))
	    {
	      
	      best_split[0]=sr;
	      best_split->varindx=k;
	      best_split->success=1;
	      init=1;
	    }
       
	}
    } 
}






void update_treematrix_time(struct split_result *sr,int nodeid,int node_row,int node_counter)
{

  // update daop.tree_matrix according to split found (historical splitting) 

  // Information for parent node is filled as is information for offspring

  // sr: pointer to split_result containing split information of best split
  // nodeid: id of node just split
  // node_row: row in tree_matrix of node just split
  // node_counter: node_counter value just after split (incremented by two from whatever it was) 


  

  int m,ni,k;
  double a;
  struct data_options *d;
  double avec[MAXCAT];
  d=&daop;
  
  // update split info for nodeid (node just split)
  d->tree_matrix[node_row][4]=((double) sr->varindx); // split variable
  d->tree_matrix[node_row][5]=sr->point; // fraction where split is performed
  //if(d->method==1||d->method==2||d->method==4)
    d->tree_matrix[node_row][6]=sr->tau; //  variable cut

  d->tree_matrix[node_row][7]=sr->delta; // delta
  d->tree_matrix[node_row][8]=(d->row_counter+1); // row number corresponding left node 
  d->tree_matrix[node_row][9]=(d->row_counter+2); // row number corresponding right node
  d->tree_matrix[node_row][10]=sr->error_reduction;
  d->tree_matrix[node_row][11]=((double) (sr->nL+sr->nR));
  d->tree_matrix[node_row][12]=-99;
  d->tree_matrix[node_row][13]=(double) sr->missing;  // information for together-missing and for number of obs (tilde-method)
  d->tree_matrix[node_row][14]=sr->delta0;


  // left node prediction
  ni=node_counter-1;
  a=node_prediction_uthash(ni);
  a=d->lambda*a;
  d->row_counter++;

  // add to tree matrix 
  d->tree_matrix[d->row_counter][0]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[d->row_counter][1]=((double) nodeid); // parent node id
  d->tree_matrix[d->row_counter][2]=((double) (node_counter-1)); // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  for(m=4;m<=10;m++)
    d->tree_matrix[d->row_counter][m]=-99; // split info
  d->tree_matrix[d->row_counter][11]=((double) sr->nL);
  d->tree_matrix[d->row_counter][12]=-99;
  d->tree_matrix[d->row_counter][13]=-99;
 d->tree_matrix[d->row_counter][14]=-99;

 if(daop.ncat>0)
   {
     node_prediction_gini(ni,avec);
     for(k=0;k<daop.ncat;k++)
       d->tree_matrix[d->row_counter][15+k]=avec[k];

   }
 
 
  // right node prediction
  a=node_prediction_uthash((node_counter));
  a=d->lambda*a;
  d->row_counter++;
  // add to tree matrix 
  d->tree_matrix[d->row_counter][0]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[d->row_counter][1]=((double) nodeid); // parent node id
  d->tree_matrix[d->row_counter][2]=((double) (node_counter)); // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  for(m=4;m<=10;m++)
    d->tree_matrix[d->row_counter][m]=-99; // split info 
  //Rprintf(" inserted right node \n");
  d->tree_matrix[d->row_counter][11]=((double) sr->nR);
  d->tree_matrix[d->row_counter][12]=-99;
  d->tree_matrix[d->row_counter][13]=-99;
  d->tree_matrix[d->row_counter][14]=-99; // delta0 
 if(daop.ncat>0)
   {
     node_prediction_gini(node_counter,avec);
     for(k=0;k<daop.ncat;k++)
       d->tree_matrix[d->row_counter][15+k]=avec[k];

   }
 
  
}



void update_treematrix(struct split_result *sr,int nodeid,int node_row,int node_counter)
{


  int m,ni,k;
  double a;
  double avec[MAXCAT];
  struct data_options *d;
  d=&daop;
 
  // update split info for nodeid (node just split)
  d->tree_matrix[node_row][4]=((double) sr->varindx); // split variable
  d->tree_matrix[node_row][5]=-99; // fraction where split is performed 
  d->tree_matrix[node_row][6]=sr->point; //  variable cut
  d->tree_matrix[node_row][7]=-99; // delta
  d->tree_matrix[node_row][8]=(d->row_counter+1); // row number corresponding left node 
  d->tree_matrix[node_row][9]=(d->row_counter+2); // row number corresponding right node
  d->tree_matrix[node_row][10]=sr->error_reduction;
  d->tree_matrix[node_row][11]=sr->nL+sr->nR;
  d->tree_matrix[node_row][12]=0;
  d->tree_matrix[node_row][13]=-99;

						     
  // left node prediction
  ni=node_counter-1;
  a=node_prediction_uthash(ni);
  a=d->lambda*a;
  //Rprintf(" prediction left node=%lf \n",a);
  d->row_counter++;
  //Rprintf(" row counter=%d \n",d->row_counter);
  // add to tree matrix 
  d->tree_matrix[d->row_counter][0]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[d->row_counter][1]=((double) nodeid); // parent node id
  d->tree_matrix[d->row_counter][2]=((double) (node_counter-1)); // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  for(m=4;m<=10;m++)
    d->tree_matrix[d->row_counter][m]=-99; // split info
  d->tree_matrix[d->row_counter][11]=sr->nL; 
  d->tree_matrix[d->row_counter][12]=-99; 
  d->tree_matrix[d->row_counter][13]=-99;

   if(daop.ncat>0)
   {
     node_prediction_gini(ni,avec);
     for(k=0;k<daop.ncat;k++)
       d->tree_matrix[d->row_counter][15+k]=avec[k];

   }
	  

  // right node prediction
  a=node_prediction_uthash((node_counter));
  a=d->lambda*a;
  d->row_counter++;
  // add to tree matrix 
  d->tree_matrix[d->row_counter][0]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[d->row_counter][1]=((double) nodeid); // parent node id
  d->tree_matrix[d->row_counter][2]=((double) (node_counter)); // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  for(m=4;m<=10;m++)
    d->tree_matrix[d->row_counter][m]=-99; // split info 
  d->tree_matrix[d->row_counter][11]=sr->nR; // prediction
  d->tree_matrix[d->row_counter][12]=-99; // prediction
  d->tree_matrix[d->row_counter][13]=-99;

 if(daop.ncat>0)
   {
     node_prediction_gini(node_counter,avec);
     for(k=0;k<daop.ncat;k++)
       d->tree_matrix[d->row_counter][15+k]=avec[k];

   }
}









/* --------------- testing functions ------------------------------------------------------------ */

 

 void print_split_result(struct split_result *s)
 {
   Rprintf(" nL=%lf nR=%lf sumL=%lf sumR=%lf ssqL=%lf ssqR=%lf tau=%lf point=%lf delta=%lf \n",s->nL,s->nR,s->sumL,s->sumR,s->sumsqL,s->sumsqR,s->tau,s->point,s->delta);
      }





/* --------------- wrapper functions ------------------------------------------------------------ */

 
void tree_wrapper(int *ptr_nsplit);
int split_history_wrapper(int nodeid,struct split_result *sr);
int split_standard_wrapper(int nodeid,struct split_result *sr);
void update_residual_wrapper(int tree_counter);
double test_error_wrapper();
double node_prediction_wrapper(int nodeid);
double node_prediction_logistic_uthash(int nodeid);
void update_residual_logistic(int tree_indx);
void update_terminal_nodes_wrapper(int tree_indx);
double test_error_logistic();
void update_terminal_nodes(int tree_indx);




void assign_method(int *method)
{

  daop.method=method[0];                // type of regression 
  if(method[0]==5)
    {
      daop.method=3;
      daop.meansummary=0; // sum of historical 
    }
  if(method[0]==6)
    {
      daop.method=5;
      daop.meansummary=1;
    }
    if(method[0]==7)
    {
      daop.method=5;
      daop.meansummary=0;
    }
     if(method[0]==8)
    {
      daop.method=6;
      daop.meansummary=1;
    }
    if(method[0]==9)
    {
      daop.method=6;
      daop.meansummary=0;
    }
    if(method[0]==10)
      daop.method=7;

    if(method[0]==11)
      daop.method=8;

    if(method[0]==12)
      daop.method=9;
}


void initialize_prediction_response()
{
  int i;

  if(daop.method_family==2) // logistic 
    {

      for(i=0;i<daop.n;i++)
	{
	  daop.predictions[i]=0;
	  daop.predictions_response[i]=1/(1+exp(-daop.predictions[i]));
	  daop.response[i]=(daop.x[daop.yindx][i]-daop.predictions_response[i]);
	}

    }
  
}


void boost_wrapper(int *random_split,int *method,int *ptr_nboost,double *te,double *predictions,int *vimp,double *rvimp)
{
  int b,i,j;
  int B;

  
  GetRNGstate();
  B=ptr_nboost[0];

  // DEBUG
  //  Rprintf("Entering 'boost_wrapper'\n");

  // methods: daop.method
  //     1: count
  //     2: boosting-logistic
  //     3: meansum-summary (missing together) (mean vs sum determined by meansummary=1/0)
  //     4: tilde (standardized count)
  //     5: meansum0 (method=6,7)
  //     6: meansum0_window (method=8,9)
  //     7: count_window  (method=10)
  
  // set some parameters 

  daop.random_split=random_split[0];    // random history splitting
  assign_method(method);                // sets daop.method 


  // .. only active for logistic regression
  initialize_prediction_response();
      

  delete_oob();
   
  set_id_rows(daop.id);

  if(1)
    {
      //Rprintf(" Entering 0-loop \n");
  for(b=0;b<B;b++)
    {
      
      
      // Rprintf(" Tree b=%d \n",b);
     if(daop.subsample==1)
       {
	 //if(daop.id_sampling==1){
	   subsample_id(b);
	   /*}else{
	   subsample_oob(b);
	   }*/
	 
       }

     
      if(daop.rf==1&&daop.subsample==0)
	set_train_oob(b);

      //Rprintf(" Going into tree_wrapper \n");
       tree_wrapper(&(daop.nsplit));

      
      update_terminal_nodes_wrapper(daop.tree_counter);
      // does 1-step NR for non-L2 loss, ie if method_family==2, else nothing 
      
      update_residual_wrapper(daop.tree_counter);
      daop.tree_counter++;

      te[b]=test_error_wrapper();
      //train_error[b]=train_error_wrapper();
      
      delete_all_nodes();      
     
    }

  if(vimp[0]==1){
  for(j=0;j<daop.p;j++)
    {
      for(i=0;i<B;i++)
	rvimp[j*B+i]=daop.varimp[j][i];
    }

  // DEBUG
  // Rprintf("Entering 'boost_wrapper'\n");

  }
    }  // if(0)
  //Rprintf("Total loop count=%d \n",daop.ncc);
  for(i=0;i<daop.n;i++)
    predictions[i]=daop.predictions[i];

   PutRNGstate();
}





double test_error_wrapper()
{

  double terr;

  terr=0;
  if(daop.method_family==1){
  if(daop.method==1||((daop.method>2)&&(daop.method<9)))
    {
      // daop.method=1: regression standard
      terr=test_error();
    }
  }
  if(daop.method_family==2)
    {
      terr=test_error_logistic();
    }
  // for other methods .... 

  if(daop.method==9)
    terr=test_error_gini();

  
  return(terr);
}


void update_residual_wrapper(int tree_counter)
{

  if(daop.method_family==1){
  if(daop.method==1||((daop.method>2)&&(daop.method<9))) // regression 
    update_residual_uthash(tree_counter);
  }
  
   if(daop.method_family==2) // regression 
    update_residual_logistic(tree_counter);

   if(daop.method==9)
      update_residual_gini(tree_counter);
   
 

  
}


int split_standard_wrapper(int nodeid,struct split_result *sr)
{
  int success;

  success=0;

  if(daop.method<=8)
    {
      // regression (standard) and logistic
      split_standard(nodeid,sr);
      success=sr->success;
      if(success==1)
	error_reduction_split_result(sr);
    
    }


  if(daop.method==9)
    {
      
      split_standard_gini(nodeid,sr);
      success=sr->success;
      if(success==1)
	error_reduction_split_gini(sr);
      //Rprintf("split_gini_standard: %d\n",success);
    }

  if(daop.method==4)
    {
      // classification tree 

    }

  return(success);
  
}


int split_history_wrapper(int nodeid,struct split_result *sr)
{

  int success;


  success=0;
     if(daop.time_split==1)
	{

	  if(daop.method==1||daop.method==2)
	    {
	      //  regression or logistic regression  
	      if(daop.random_split==0){
		split_historical(nodeid,sr);
	      }else{
		split_historical_sample(nodeid,sr);
	      }
	  
	      success=sr->success;
	      if(success==1)
		error_reduction_split_result(sr);
	    }



	  if(daop.method==3)
	    {
	      // history summarization mean(delta) or sum(delta)
	      //	      Missing values using MISSING-TOGETHER method 
	      split_meansum(nodeid,sr);
	      success=sr->success;

	      if(success==1)
		error_reduction_split_result(sr);
	    }

	  if(daop.method==4)
	    {
	      //
	      //split_history_tilde();
	      split_historical_tilde(nodeid,sr);
	      success=sr->success;
	      //Rprintf("DEBUG: success=%d \n",sr->success);
	      
	      if(success==1){
		error_reduction_split_result(sr);
		//Rprintf("DEBUG: successful split\n");
		//print_split_result(sr);
	      }
		
	    }

	   if(daop.method==5)
	    {
	      // history summarization mean(delta) or sum(delta)
	      //	      Missing values set to 0 
	      split_meansum0(nodeid,sr);
	      success=sr->success;

	      if(success==1)
		error_reduction_split_result(sr);
	    }


	    if(daop.method==6)
	    {
	      // windowed-history summarization mean(delta) or sum(delta)
	      //	      Missing values set to 0 
	      split_meansum0_window(nodeid,sr);
	      success=sr->success;

	      if(success==1)
		error_reduction_split_result(sr);
	    }


	   if(daop.method==7)
	    {
	      // windowed-history splitting using window count
	      split_window(nodeid,sr);
	      success=sr->success;

	      if(success==1)
		error_reduction_split_result(sr);
	    }


	  if(daop.method==8)
	    {
	      // windowed-history splitting using window count
	      split_window_tilde(nodeid,sr);
	      success=sr->success;

	      if(success==1)
		error_reduction_split_result(sr);
	    }


	   if(daop.method==9)
	    {
	      // windowed-history splitting for classification (gini-index)
	      //Rprintf(" Start split_window_gini\n");
	      split_window_gini(nodeid,sr);
	      success=sr->success;

	      if(success==1)
		error_reduction_split_gini(sr);
	    }


	    
	}

     return(success);
}


void tree_wrapper(int *ptr_nsplit)
{

  int i,m,k,j;
  int nsplit,node_counter;
  int selector;
  int best_var,best_var_2;
  int success,success_2,nodeid;
  int *nodes_to_split;
  int *node_rownumber;
  int nodes_not_split;
  double a;
  struct data_options *d;
  int ni;
  struct split_result best_split;
  struct split_result best_split_time;
  struct split_result spr;
  double ssq_reduction,ssq_reduction_time;
  
  d=&daop;
  
  nsplit=ptr_nsplit[0];
  nodes_to_split=(int*) malloc(sizeof(int)*10*nsplit);
  node_rownumber=(int*) malloc(sizeof(int)*10*nsplit);
  
  d->row_counter++; // increment  
  d->tree_start[d->tree_counter]=d->row_counter;

  // initialize root node 
  j=0;
  for(i=0;i<d->n;i++)
    add_row(j,i);

  // Set up root-node info 
  nodes_to_split[0]=0;
  node_rownumber[0]=d->row_counter;
  // mean of response
  ni=0;
  a=node_prediction_uthash(ni); 
  a=d->lambda*a;       // shrink mean by shrinkage factor lambda.
  
  d->tree_matrix[d->row_counter][0]=d->tree_counter; // tree indicator
  d->tree_matrix[d->row_counter][1]=-99; // parent node id
  d->tree_matrix[d->row_counter][2]=0; // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  d->tree_matrix[d->row_counter][4]=-99; // .. this is the terminal node until we know more!!! 

 
  
  nodes_not_split=1;
  node_counter=0;
  nodeid=0;
  for(k=0;k<nsplit;k++)
    {
     
      if(k>node_counter)
	break;
      
      nodeid=nodes_to_split[k];

      
      success=split_history_wrapper(nodeid,&best_split_time);
      success_2=split_standard_wrapper(nodeid,&best_split);
        
	
      // Determine which split: best historical or best standard 
      selector=-1;
      if(success_2>0&&success>0&&(best_split.error_reduction>=best_split_time.error_reduction)) // choose simpler when the same.. 
        selector=2;
      if(success_2>0&&success==0)
        selector=2;
      if(success_2>0&&success>0&&(best_split.error_reduction<best_split_time.error_reduction))
        selector=1;
      if(success_2==0&&success>0)
        selector=1;

      if(selector>0)
	{
	  for(i=1;i<=2;i++){
	    node_counter++;
	    nodes_to_split[node_counter]=node_counter;
	    node_rownumber[node_counter]=(d->row_counter+i);
	  }

	}
 
      if(selector==1)
	{
	  update_nodevector_time_wrapper(nodeid,(node_counter-1),node_counter,&best_split_time);
	  update_treematrix_time(&best_split_time,nodeid,node_rownumber[nodeid],node_counter);
	}


        if(selector==2)
	{
	  update_nodevector(nodeid,(node_counter-1),node_counter,best_split.varindx,best_split.point);
	  update_treematrix(&best_split,nodeid,node_rownumber[nodeid],node_counter); 
	}
    }

  d->tree_end[d->tree_counter]=d->row_counter;
  free(nodes_to_split);
  free(node_rownumber);

}


double node_prediction_wrapper(int nodeid)
{
  // Called by: update_terminal_nodes(_wrapper)
  double prediction;
  prediction=0;
  if(daop.method_family==1){
  if(daop.method==1||daop.method==3||daop.method==4)
    {
      // method=1: regression
      prediction=node_prediction_uthash(nodeid);
    }
  }

  if(daop.method_family==2)
    {
      // logistic regression 
      prediction=node_prediction_logistic_uthash(nodeid);
    }

  return(prediction);
}






double node_prediction_logistic_uthash(int nodeid)
{
  // This is the node prediction at the terminal nodes for logistic tree boosting
  // Trees themselves are fit using least squares 
  
  // return node mean
  int i;
  double a,b;
  double n;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  
  a=0;
  b=0;
  n=0;

  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    //for(i=0;i<daop.n;i++)
    HASH_ITER(hh,s->rows,r,rtemp)
    {
      i=r->row_number;
	  if(daop.train[i]==1)
	    {
	      // THIS DONT LOOK RIGHT 
	  a+=daop.response[i];
	  b+=((daop.predictions_response[i])*(1-(daop.predictions_response[i])));
	  n++;
	    }
       
    }
  }

  if(n<.5)
    n=1;
  //  if(b<.00001)
    //b=1;




  // testing
  //a=a/(n*(daop.mean_y*(1-daop.mean_y)));
  //Rprintf(" a=%lf b=%lf mf=%d \n",a,b,daop.method_family);
  a=a/(b+.0001);
  return(a);
}






void update_treematrix_wrapper(struct split_result *sr,int nodeid,int node_row,int node_counter)
{


  int m,ni;
  double a;
  struct data_options *d;
  d=&daop;
 
  // update split info for nodeid (node just split)
  d->tree_matrix[node_row][4]=((double) sr->varindx); // split variable
  d->tree_matrix[node_row][5]=-99; // fraction where split is performed 
  d->tree_matrix[node_row][6]=sr->point; //  variable cut
  d->tree_matrix[node_row][7]=-99; // delta
  d->tree_matrix[node_row][8]=(d->row_counter+1); // row number corresponding left node 
  d->tree_matrix[node_row][9]=(d->row_counter+2); // row number corresponding right node
  d->tree_matrix[node_row][10]=sr->error_reduction;
  d->tree_matrix[node_row][11]=sr->nL+sr->nR;
  d->tree_matrix[node_row][12]=0;
  d->tree_matrix[node_row][13]=-99;

						     
  // left node prediction
  ni=node_counter-1;
  a=node_prediction_uthash(ni);
  a=d->lambda*a;
  //Rprintf(" prediction left node=%lf \n",a);
  d->row_counter++;
  //Rprintf(" row counter=%d \n",d->row_counter);
  // add to tree matrix 
  d->tree_matrix[d->row_counter][0]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[d->row_counter][1]=((double) nodeid); // parent node id
  d->tree_matrix[d->row_counter][2]=((double) (node_counter-1)); // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  for(m=4;m<=10;m++)
    d->tree_matrix[d->row_counter][m]=-99; // split info
  d->tree_matrix[d->row_counter][11]=sr->nL; 
  d->tree_matrix[d->row_counter][12]=-99; 
  d->tree_matrix[d->row_counter][13]=-99;
	  

  // right node prediction
  a=node_prediction_uthash((node_counter));
  a=d->lambda*a;
  d->row_counter++;
  // add to tree matrix 
  d->tree_matrix[d->row_counter][0]=((double) d->tree_counter); // tree indicator
  d->tree_matrix[d->row_counter][1]=((double) nodeid); // parent node id
  d->tree_matrix[d->row_counter][2]=((double) (node_counter)); // node id
  d->tree_matrix[d->row_counter][3]=a; // prediction
  for(m=4;m<=10;m++)
    d->tree_matrix[d->row_counter][m]=-99; // split info 
  d->tree_matrix[d->row_counter][11]=sr->nR; // prediction
  d->tree_matrix[d->row_counter][12]=-99; // prediction
  d->tree_matrix[d->row_counter][13]=-99;


}



void update_terminal_nodes_wrapper(int treeindx)
{

  // For loss functions other than squared error, a single newton step is performed to fit the terminal node fits.
  // This function performs the update.


  // if(daop.method==1) : do nothing. 

  // This function just decides whether to run update_terminal_nodes() function
  // (not needed for regression, but needed for logistic regression, fex)
  
  if(daop.method_family==2)
    {
      update_terminal_nodes(treeindx);
    }

  

}





void update_terminal_nodes(int tree_indx)
{
  int k,i,start,go,ni,ti,end;
  double node_predictions[10000]; // no more than 100 terminal nodes
  int tnodes[10000];// terminal nodes
  int ntnodes;
  double h;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[k][0];
      if(ti==tree_indx)
	{

	  if(daop.tree_matrix[k][4]<(-1))
	    {
	      // THIS is a terminal node ....
	      ni=((int) daop.tree_matrix[k][2]);
	      // node_predictions[ni]=daop.tree_matrix[k][3];
	      h=node_prediction_wrapper(ni);
	      daop.tree_matrix[k][3]=(daop.lambda*h);


	    }
	  
	}else{
	go=0;
	break;
      }

    }



 
}

  




double test_error_logistic()
{

  
  double error,h;
  int i,ntest;
  ntest=0;
  error=0;


  for(i=0;i<daop.n;i++)
    {
      if(daop.train[i]==0)
	{
	  h=daop.x[daop.yindx][i]*log(daop.predictions_response[i]);
	  h=h+(1-daop.x[daop.yindx][i])*log((1-daop.predictions_response[i]));
	  error+=h;
	  ntest++;
	}
    }
    
  if(ntest>0)
    error=error/((double) ntest);

  return(error);
}





void update_residual_logistic(int tree_indx)
{
  int k,i,start,go,ni,ti,end;
  double node_predictions[10000]; // no more than 100 terminal nodes
  int tnodes[10000];// terminal nodes
  int ntnodes;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  // Rprintf(" update_residual_logistic \n");
  ntnodes=0;
  start=daop.tree_start[tree_indx];
  end=daop.tree_end[tree_indx];
  //Rprintf(" start=%d end=%d \n",start,end);  
  for(k=start;k<=end;k++)
    {

      ti=(int) daop.tree_matrix[k][0];
      if(ti==tree_indx)
	{
	  ni=((int) daop.tree_matrix[k][2]);
	  //Rprintf(" node id =%d k=%d \n",ni,k);
	  node_predictions[ni]=daop.tree_matrix[k][3];
	  if(daop.tree_matrix[k][4]<(-1))
	    {
	      tnodes[ntnodes]=ni;
	      ntnodes++;
	    }
	  
	}else{
	go=0;
	break;
      }

    }
 
  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
    
 	  i=r->row_number;
	  daop.predictions[i]+=(node_predictions[ni]);
	  daop.predictions_response[i]=1/(1+exp(-daop.predictions[i]));
	  daop.response[i]=(daop.x[daop.yindx][i]-daop.predictions_response[i]);

	}
      }
    }

  //  Rprintf(" finished last loop \n");
  
}






// ----- quantile regression functions ---------------------------------



// tnode_row

void add_tnode_row(int treeindx,int nodeindx,int rowindx)
{
 


  struct nodelist *t;
  struct node *s;
  struct rnumber *r;

  HASH_FIND_INT(tnode_row,&treeindx,t);
  if(t==NULL)
    {
      t=(struct nodelist*) malloc(sizeof(struct nodelist));
      t->list_number=treeindx;
      t->nodes=NULL;
      HASH_ADD_INT(tnode_row,list_number,t);
      //utarray_new(s->rows_array,&ut_int_icd);
    }
  
 
  HASH_FIND_INT(t->nodes,&nodeindx,s);
  if(s==NULL)
    {
      s=(struct node*) malloc(sizeof(struct node));
      s->node_number=nodeindx;
      s->rows=NULL;
      s->nrows=0;
      HASH_ADD_INT(t->nodes,node_number,s);
      //utarray_new(s->rows_array,&ut_int_icd);
    }

  HASH_FIND_INT(s->rows,&rowindx,r);
  if(r==NULL)
    {
      r=(struct rnumber*) malloc(sizeof(struct rnumber));
      HASH_ADD_INT(s->rows,row_number,r);
      r->row_number=rowindx;
      
      //utarray_push_back(s->rows_array,&r_number); 
    }

 
}



void delete_tnode_row()
{
  struct nodelist *t,*ttemp;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;

  HASH_ITER(hh,tnode_row,t,ttemp){
  HASH_ITER(hh,t->nodes,s,stemp){
    HASH_ITER(hh,s->rows,r,rtemp){
      HASH_DEL(s->rows,r);
      free(r);
    }
   
    HASH_DEL(t->nodes,s);
    free(s);
  }
  HASH_DEL(tnode_row,t);
  free(t);
  }
}





int find_tnode(int tree_indx,int row_number)
{
  // For subject in row_number, find terminal node for subject, return terminal node number

  int k,tnode,go,j,i,split_var;
  double ff;
  int counter;
  int pred;

  tnode=-99;
  //for(k=daop.tree_start[tree_indx];k<=daop.tree_end[tree_indx];k++)
  k=daop.tree_start[tree_indx];
  //Rprintf(" Row corresponding tree start %d \n",k);
  go=1;
  counter=0;
   while(go==1)
    {

      //Rprintf(" split-var =%d \n",((int) daop.tree_matrix[k][4]));
      if(daop.tree_matrix[k][4]<(-1))
	{
	  // terminal node
	  pred=((int) daop.tree_matrix[k][2]);
	  break;
	}else{
	split_var=((int) daop.tree_matrix[k][4]);

	if(daop.tree_matrix[k][5]<(-1)) // not a summary variable
	  {
	    if(daop.x[split_var][row_number]<daop.tree_matrix[k][6])
	      {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }

	  }else{
	  ff=tsummary_row_fast(row_number,daop.tree_matrix[k][6],daop.tree_matrix[k][7],split_var);

	  if(ff<daop.tree_matrix[k][5])
	    {
		// Go Left
		k=((int) daop.tree_matrix[k][8]);
	      }else{
	      // Go Right
	      k=((int) daop.tree_matrix[k][9]);
	    }
	}

      }
    }

   return(pred);
}





void quantile_aux(double *prob,double *y,int *n,double *weights,double *quant)
{
  double cc,cc_new;
  int k;

  quant[0]=-99;
  cc=0;
  for(k=0;k<n[0];k++)
    {
      cc_new=cc+weights[k];
      if(cc<=prob[0]&&cc_new>=prob[0])
	{
	  quant[0]=y[k-1]+(y[k]-y[k-1])*(prob[0]-cc)/(cc_new-cc);
	  break;
	}
	    cc=cc_new;
    }
  
}



void quantile_R(double *prob,int *np,double *y,int *n,double *weights,int *m,double *quant)
{

  // find quantiles prob[1:np] for m observations, weights in weights
  int k,j,counter,i;
  double q;

  counter=0;
  for(k=0;k<m[0];k++)
    {

      for(j=0;j<np[0];j++)
	{

	  quantile_aux(&(prob[j]),y,n,&(weights[k*n[0]]),&q);
	  quant[counter]=q;
	  counter++;
	}


    }


}



void set_tnode_row()
{
  // Fill up tnode_row list, mapping rows to terminal nodes of each tree. 
  
  int i,k,tnode;
  delete_tnode_row();

  for(k=0;k<daop.nboost;k++)
    {
      for(i=0;i<daop.n;i++){
	tnode=find_tnode(k,i);
	add_tnode_row(k,tnode,i);
      }

    }

}

void get_rows(int *treeindx,int *nodeindx,int *res,int *n)
{
  int i,k,tnode;

  struct nodelist *t,*ttemp;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int counter;

  HASH_FIND_INT(tnode_row,treeindx,t);
  counter=0;
  if(t==NULL){}else{

    HASH_FIND_INT(t->nodes,nodeindx,s);
    if(s==NULL){}else{
      
	      HASH_ITER(hh,s->rows,r,rtemp){
		res[counter]=r->row_number;
		counter++;
	      }
    }

  }

  n[0]=counter;
}

void get_weights(int *n,int *res)
{
  int i,k,tnode;

  struct nodelist *t,*ttemp;
  struct node *s,*stemp;
  struct rnumber *r,*rtemp;


  for(i=0;i<(n[0]*daop.n);i++)
    res[i]=0;
  
  for(i=0;i<daop.n;i++)
    {

      for(k=0;k<daop.nboost;k++)
	{

	  tnode=find_tnode(k,i);
	  HASH_FIND_INT(tnode_row,&k,t);
	  if(t==NULL){}else{
	    HASH_FIND_INT(t->nodes,&tnode,s);
	    if(s==NULL){}else{
	      HASH_ITER(hh,s->rows,r,rtemp){
		res[i*n[0]+r->row_number]+=1;
	      }
	    }
	  }
	}

    }

}
