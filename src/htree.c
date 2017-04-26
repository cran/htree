#include <R.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "uthash.h"
#include "utarray.h"

#define EPSILON 0.000001

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
  //
  double sumL;
  double sumR;
  double sumsqL;
  double sumsqR;
  
  
};

int find_tnode(int tree_indx,int row_number);


void print_split_result(struct split_result *sr);
void create_history_summary_uthash(int nodeid,int vindx,struct split_result *sr,double *x);
void permute(int *p,int n);


struct data_options{

  // .. testing data
  double **split_matrix;

  // .. testing data finished.... 
  
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
  int *iaux1,*iaux2,*iaux3,*iaux4,*iaux5,*iaux6;

  // Split info for fast splitting
  //  int *nstat;
  //  double *ssqL,*ssqR,*sumL,*sumR,*nL,*nR;
 

  
  // predictions
  double *predictions;
  double *predictions_response;
  int *oob;


  double mean_y;
  // split info and parameters 
  int min_nodesize;
  int nsamp;
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
  int **n_changed;       // keeps track of which n_i(tau) change (for given tau)
  int *counter_changed; // keeps track of number of n_i(tau) changed (given tau)

  struct split_result **sr; 
};

struct data_options daop;


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



void read_data(double *x,int *n,int *p,double *time,int *id,int *yindx,double *lambda,int *nsplit,int *nboost,int *rf,int *nsamp,int *time_split)
{
  int i,j,ncol,k;
  struct data_options *d;
  double *v;
  int max_history_length;
  d=&daop;

  // THIS NEEDS TO SET AS PARAMETER LATER .... 
  daop.method=1;

  // start: testing data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  d->split_matrix=(double**) malloc(sizeof(double*)*(1000));
  v=(double*) malloc(sizeof(double)*1000*13);
  d->split_matrix[0]=v;
  for(i=1;i<1000;i++)
    d->split_matrix[i]=d->split_matrix[i-1]+13;
  // end: testing data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  
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
  
  // matrix of split_result structures
  d->sr=(struct split_result**) malloc(sizeof(struct split_result*)*nsamp[0]);
  for(k=0;k<nsamp[0];k++)
    d->sr[k]=(struct split_result*) malloc(sizeof(struct split_result)*(i+2)); // way to many...


  d->n_changed=(int**) malloc(sizeof(int*)*nsamp[0]);       // keeps track of which n_i(tau) change (for given tau)
  for(k=0;k<nsamp[0];k++)
    d->n_changed[k]=(int*) malloc(sizeof(int)*i*n[0]);
  
  d->counter_changed=(int*) malloc(sizeof(int)*nsamp[0]); // keeps track of number of n_i(tau) changed (given tau)

  d->n_row_number=(double**) malloc(sizeof(double*)*nsamp[0]);
  for(k=0;k<nsamp[0];k++)
    d->n_row_number[k]=(double*) malloc(sizeof(double)*i*n[0]);

  max_history_length=i;
  // -------------------------------------------------------
  ncol=13;

  d->splitvar_history=(int*) malloc(sizeof(int)*p[0]);
  d->splitvar_concurrent=(int*) malloc(sizeof(int)*p[0]);
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
    d->tree_matrix[i]=d->tree_matrix[i-1]+ncol;

  for(i=0;i<j;i++)
    for(k=0;k<ncol;k++)
      d->tree_matrix[i][k]=0;


  

  
  // Rprintf(" read in data \n");
  d->varimp=(double**) malloc(sizeof(double*)*p[0]);
  v=(double*) malloc(sizeof(double)*d->nboost*p[0]);
  d->varimp[0]=v;
  for(i=1;i<p[0];i++)
    d->varimp[i]=d->varimp[i-1]+d->nboost;
  
  d->row_counter=-1; // this is incremented as rows are added
  d->tree_counter=0;

  d->tree_start=(int*) malloc(sizeof(int)*d->nboost);
  d->tree_end=(int*) malloc(sizeof(int)*d->nboost);

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
  d->nsamp_var=d->nsamp*d->nsamp*d->nsamp;
  d->min_nodesize=5;
  
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


  daop.mean_y=0;
  for(i=0;i<d->n;i++)
    daop.mean_y+=d->response[i];

  daop.mean_y=daop.mean_y/((double) d->n);
  
   d->train=(int*) malloc(sizeof(int)*d->n);
     for(j=0;j<d->n;j++)
    d->train[j]=1;

   
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


    d->predictions_response=(double*) malloc(sizeof(double)*d->n);
  for(j=0;j<d->n;j++)
    d->predictions_response[j]=0;

// oob counter vector 
  d->oob=(int*) malloc(sizeof(int)*d->n);
  for(j=0;j<d->n;j++)
    d->oob[j]=0;
  
  // allocate and get range of predictor variables 
  d->lower=(double*) malloc(sizeof(double)*d->p);
  d->upper=(double*) malloc(sizeof(double)*d->p);

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
  d->node=(int*) malloc(sizeof(int)*d->n);
  for(i=0;i<d->n;i++)
    d->node[i]=0;

  
  // scratch space
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

 
}



void get_tree(int *tree_number,double *res,int *n)
{
  int i,j,k;
  i=0;
  for(k=0;k<=daop.row_counter;k++)
    {

      if(daop.tree_matrix[k][0]==tree_number[0])
	{
	  for(j=0;j<=11;j++)
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
	  for(j=0;j<=12;j++)
	    {
	      res[i]=daop.tree_matrix[k][j];
	      i++;
	    }
    }
  n[0]=i;
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

 d=&daop;

 for(i=0;i<d->nsamp;i++)
   free(d->n_row_number[i]);
 free(d->n_row_number);

 
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





void read_predict(double *x,int *n,int *p,double *time,int *id,int *yindx,double *trees,int *nrow_trees,int *nboost,int *time_split)
{
  // Read in data needed for predictions 

 
  
  int i,j,ncol,ctree;
  struct data_options *d;
  double *v;
  int ht;
  d=&daop;

  ncol=13;
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

  // loop through trees and form permuted oob predictions 

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

  if(n<.5)
    n=1;
  a=a/n;

  return(a);
}




void update_nodevector_time(int nodeid,int nodeL,int nodeR,struct split_result *sr)
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
  create_history_summary_uthash(nodeid,vindx,sr,x);

  cut=sr->point;


  //for(i=0;i<d->n;i++)
  HASH_FIND_INT(nodes,&nodeid,s);
  if(s==NULL){}else{
    HASH_ITER(hh,s->rows,r,rtemp)
  {

    i=r->row_number;
	  ff=x[i];
	  if((ff+EPSILON)<cut)
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
      // y[counter_obs]=d->x[d->yindx][i];
      //id[counter_obs]=d->id[i];
      //time[counter_obs]=d->time[i];
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




void initialize_split_result(struct split_result *spr)
{
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
  

}

void error_split_result(struct split_result *spr)
{

  // Compute error of split result (also predL and predR are formed) 
  spr->predL=spr->sumL/spr->nL;
  spr->predR=spr->sumR/spr->nR;
  spr->error=(spr->nL/(spr->nL+spr->nR))*(spr->sumsqL/spr->nL-spr->predL*spr->predL)+(spr->nR/(spr->nL+spr->nR))*(spr->sumsqR/spr->nR-spr->predR*spr->predR);



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

void get_split_matrix(double *res)
{

  int i,j;
  for(i=0;i<1000;i++)
    for(j=0;j<13;j++)
      res[i*13+j]=daop.split_matrix[i][j];

}

void record_split(struct split_result *sr,int counter)
{

  daop.split_matrix[counter][0]=sr->tau;
  daop.split_matrix[counter][1]=sr->error;
  daop.split_matrix[counter][2]=sr->delta;
  daop.split_matrix[counter][3]=sr->point;
  daop.split_matrix[counter][4]=sr->nL;
  daop.split_matrix[counter][5]=sr->nR;
  daop.split_matrix[counter][6]=sr->sumL;
  daop.split_matrix[counter][7]=sr->sumR;
  daop.split_matrix[counter][8]=sr->sumsqL;
  daop.split_matrix[counter][9]=sr->sumsqR;
  daop.split_matrix[counter][10]=sr->error_reduction;
  daop.split_matrix[counter][11]=((double) sr->varindx);
  daop.split_matrix[counter][12]=0;
  

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

      // Print split result ...
      /*  Rprintf(" ####   Print split result for delta=%lf ######################## \n",delta[k]);
      print_split_result(&sresult);
      */

      // ... WORKING HERE ... what the is initialize_sr() doing?? 
      // fill up daop.sr recording split-information for each value of n[delta,tau,k]
      // fills up daop.sr[j][k] giving split info for splitting on n_ih(theta_j)<k => left, and n_ih>=k => right
      initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);

      
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      
	      best_split=sresult;
	      best_split.delta=(delta[0]+100);
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
  
  for(k=0;k<daop.nsamp;k++)
    daop.counter_changed[k]=0;

 
  record_counter=0;
  for(k=0;k<npast;k++)
    {

      if(fabs(x[k]-xc)>EPSILON)
	{
	  /* if(unif_rand()<.01) */
	  /*   { */
	  
	  // new split point, evaluate previous one 
	  for(j=0;j<daop.nsamp;j++)
	    {
	      if(daop.counter_changed[j]>0)
		{
		  for(i=0;i<daop.counter_changed[j];i++)
		    {
		      n_i=daop.n_changed[j][i];
		      sp=&daop.sr[j][n_i];
		      if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
			{
			  
			  error_split_result(sp);
			  // ########### START: TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			  sp->tau=xc;
			  sp->delta=delta[j];
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
  

      
  } // matches if(0){.. 
  return_split[0]=best_split;

 
      }
     return init;
}


int split_attempt_history_testing(int nodeid,int vindx,double *delta,double tau_ext,int n_ext,struct split_result *return_split);


void split_test(int *vindx,double *delta,double *tau,int *n,double *result)
{

  // THIS and split_attempt_..._testing function are for extracting results of a single split 
  int i,j;
  struct split_result sr;
  
  // get result for split on vindx for settings:delta,tau,n => result
  

  // set nodeid==0 for all data

 // initialize root node 
  j=0;
  for(i=0;i<daop.n;i++)
    add_row(j,i);  

  daop.nsamp=1;
  
  split_attempt_history_testing(j,vindx[0],delta,tau[0],n[0],&sr);



  result[0]=sr.nL;
  result[1]=sr.nR;
  result[2]=sr.sumL;
  result[3]=sr.sumR;
  result[4]=sr.sumsqL;
  result[5]=sr.sumsqR;
  

  // delete nodes
  delete_all_nodes();
}


int split_attempt_history_testing(int nodeid,int vindx,double *delta,double tau_ext,int n_ext,struct split_result *return_split)
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
 
  return_split->success=0;
  init=0;
  
  x=daop.daux1;
  time_past=daop.daux2;
  y=daop.daux3;
  time_current=daop.daux4;

  ri_past=daop.iaux1;    // row index associated with response, attached to each historical value
  ri_node=daop.iaux4;    // unique values of ri_past

  extract_history_uthash(nodeid,vindx,x,time_past,ri_past,&npast,y,time_current,ri_node,delta[0],&n);
 
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


	// At this point, n_row_number[j][rnumber] should contain the number of observations in history defined by delta[j]
	// Looks pretty good, the first rows... 
	/*	for(k=0;k<10;k++)
	  {
	    for(j=0;j<daop.nsamp;j++)
	      Rprintf(" nrn[%d][%d]=%lf ",j,k,daop.n_row_number[j][k]);
	    Rprintf(" \n"); 
	  }
	Rprintf(" ### delta ### \n");
	for(j=0;j<daop.nsamp;j++)
	  Rprintf(" delta[%d]=%lf ",j,delta[j]);
	Rprintf("\n");
	*/
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

      // Print split result ...
      /*  Rprintf(" ####   Print split result for delta=%lf ######################## \n",delta[k]);
      print_split_result(&sresult);
      */

      // ... WORKING HERE ... what the is initialize_sr() doing?? 
      // fill up daop.sr recording split-information for each value of n[delta,tau,k]
      // fills up daop.sr[j][k] giving split info for splitting on n_ih(theta_j)<k => left, and n_ih>=k => right
      initialize_sr(k,vindx,daop.daux5,daop.daux6,daop.iaux5,n);

      
      
      if(sresult.success==1)
	{
	  if(init==0||(sresult.error<best_split.error))
	    {
	      
	      best_split=sresult;
	      best_split.delta=(delta[0]+100);
	      best_split.varindx=vindx;
	      init=1;
	      
	    }
	}
  

    }


     // Print best split result ...
  /*  Rprintf(" ####   Print best-split result for delta=%lf ######################## \n",delta[k]);
      print_split_result(&best_split);
  */
  
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
	  // new split point, evaluate previous one 
	  for(j=0;j<daop.nsamp;j++)
	    {
	      if(daop.counter_changed[j]>0)
		{
		  for(i=0;i<daop.counter_changed[j];i++)
		    {
		      n_i=daop.n_changed[j][i];
		      sp=&daop.sr[j][n_i];
		      if((sp->nL>=mnsd)&&(sp->nR>=mnsd))
			{
			  
			  error_split_result(sp);
			  // ########### START: TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			  sp->tau=xc;
			  sp->delta=delta[j];
			  sp->point=n_i;
			  if(record_counter<1000){
			    record_split(sp,record_counter);
			    record_counter++;
			  }
			  // ########### END: TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

	  
	    }
	}

      xc=x[k];
      dd=(time_current[daop.iaux3[k]]-time_past[daop.iaux3[k]]);
      for(j=0;j<daop.nsamp;j++)
	{

	  if(dd<=delta[j])
	    {


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

        if((fabs(x[k]-tau_ext)<EPSILON)&&(fabs(x[k+1]-tau_ext)>EPSILON))
	 {
	   // ... need fill in rest .. 
	   return_split[0]=daop.sr[0][n_ext];
	   break;
	 }

    }

      
  } // matches if(0){.. 

  // return tau,delta,n to return_split
  // return_split[0]=best_split;

 
      }
     return init;
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



void sample_delta(int nodeid,double *delta)
{

  struct node *s,*stemp;
  struct rnumber *r,*rtemp;
  int init,i,k;
  double max_time,a;
  
 
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






void update_treematrix_time(struct split_result *sr,int nodeid,int node_row,int node_counter)
{

  // update daop.tree_matrix according to split found (historical splitting) 

  // Information for parent node is filled as is information for offspring

  // sr: pointer to split_result containing split information of best split
  // nodeid: id of node just split
  // node_row: row in tree_matrix of node just split
  // node_counter: node_counter value just after split (incremented by two from whatever it was) 


  

  int m,ni;
  double a;
  struct data_options *d;
  d=&daop;
  
  // update split info for nodeid (node just split)
  d->tree_matrix[node_row][4]=((double) sr->varindx); // split variable
  d->tree_matrix[node_row][5]=sr->point; // fraction where split is performed 
  d->tree_matrix[node_row][6]=sr->tau; //  variable cut
  d->tree_matrix[node_row][7]=sr->delta; // delta
  d->tree_matrix[node_row][8]=(d->row_counter+1); // row number corresponding left node 
  d->tree_matrix[node_row][9]=(d->row_counter+2); // row number corresponding right node
  d->tree_matrix[node_row][10]=sr->error_reduction;
  d->tree_matrix[node_row][11]=((double) (sr->nL+sr->nR));
  d->tree_matrix[node_row][12]=-99;

						     
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
  d->tree_matrix[d->row_counter][11]=(double) sr->nL;
  d->tree_matrix[d->row_counter][12]=-99;
	  

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
  d->tree_matrix[d->row_counter][11]=sr->nR;
  d->tree_matrix[d->row_counter][12]=-99;


}



void update_treematrix(struct split_result *sr,int nodeid,int node_row,int node_counter)
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


}




void tree_historic(int *ptr_nsplit)
{

  int i,m,k,j;
  int nsplit,node_counter;
  double best_frac,best_cut,best_delta,best_sumsq;
  double best_frac_2,best_cut_2,best_delta_2,best_sumsq_2;
  int best_n_left,best_n_right,best_vtype;
  int best_n_left_2,best_n_right_2;
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

      success=0;
      if(d->time_split==1)
	{
	  if(daop.random_split==0){
	    split_historical(nodeid,&best_split_time);
	  }else{
	  split_historical_sample(nodeid,&best_split_time);
	  }
	  
	  success=best_split_time.success;
	  if(success==1)
	    error_reduction_split_result(&best_split_time);
	}


      split_standard(nodeid,&best_split);
      success_2=best_split.success;
      if(success_2==1)
	error_reduction_split_result(&best_split);
 				   
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
	  update_nodevector_time(nodeid,(node_counter-1),node_counter,&best_split_time);
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


void boost_uthash_dev(int *random_split,int *ptr_nboost,double *te,double *predictions,int *vimp,double *rvimp)
{
  int b,i,j;
  int B;

   GetRNGstate();
  B=ptr_nboost[0];

  daop.random_split=random_split[0];
 
  for(b=0;b<B;b++)
    {
     if(daop.subsample==1)
	subsample_oob(b);

      if(daop.rf==1&&daop.subsample==0)
	set_train_oob(b);
 
      tree_historic(&(daop.nsplit));

      update_residual_uthash(daop.tree_counter);
      daop.tree_counter++;

      te[b]=test_error();

      delete_all_nodes();      
    }

  if(vimp[0]==1){
  for(j=0;j<daop.p;j++)
    {
      for(i=0;i<B;i++)
	rvimp[j*B+i]=daop.varimp[j][i];
    }
  
  }

  for(i=0;i<daop.n;i++)
    predictions[i]=daop.predictions[i];


   PutRNGstate();
}










/* --------------- testing functions ------------------------------------------------------------ */

 

 void print_split_result(struct split_result *s)
 {
	Rprintf(" nL=%lf nR=%lf sumL=%lf sumR=%lf tau=%lf point=%lf delta=%lf \n",s->nL,s->nR,s->sumL,s->sumR,s->tau,s->point,s->delta);
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





void boost_wrapper(int *random_split,int *method,int *ptr_nboost,double *te,double *predictions,int *vimp,double *rvimp)
{
  int b,i,j;
  int B;

  GetRNGstate();
  B=ptr_nboost[0];

  
  // set some parameters 
  daop.method=method[0];                // type of regression 
  daop.random_split=random_split[0];    // random history splitting 

  delete_oob();
   
  set_id_rows(daop.id);
   
  for(b=0;b<B;b++)
    {
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
 
      tree_wrapper(&(daop.nsplit));

      
      update_terminal_nodes_wrapper(daop.tree_counter); // does 1-step NR for non-L2 loss functions
      
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
  
  }

  for(i=0;i<daop.n;i++)
    predictions[i]=daop.predictions[i];

   PutRNGstate();
}


double test_error_wrapper()
{

  double terr;

  terr=0;
  if(daop.method==1)
    {
      // daop.method=1: regression standard
      terr=test_error();
    }
  if(daop.method==2)
    {
      terr=test_error_logistic();
    }
  // for other methods .... 

  return(terr);
}


void update_residual_wrapper(int tree_counter)
{

  if(daop.method==1) // regression 
    update_residual_uthash(tree_counter);

   if(daop.method==2) // regression 
    update_residual_logistic(tree_counter);


    //  if(daop.method==3) // classification 
    //update_residual_classification(tree_counter);


  
}


int split_standard_wrapper(int nodeid,struct split_result *sr)
{
  int success;

  success=0;

  if(daop.method==1||daop.method==2)
    {
      // regression (standard) and logistic
      split_standard(nodeid,sr);
      success=sr->success;
      if(success==1)
	error_reduction_split_result(sr);

    }



  if(daop.method==3)
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
	      // classification trees (for random forest) 

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
  a=node_prediction_uthash(ni); // not  ..._wrapper function here, thats for the terminal nodes .. 
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
	  update_nodevector_time(nodeid,(node_counter-1),node_counter,&best_split_time);
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

  double prediction;
  prediction=0;
  if(daop.method==1)
    {
      // method=1: regression
      prediction=node_prediction_uthash(nodeid);
    }

  if(daop.method==2)
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
	  b+=(daop.predictions_response[i])*(1-(daop.predictions_response[i]));
	  n++;
	    }
       
    }
  }

  if(n<.5)
    n=1;
  if(b<.00001)
    b=1;
  //a=a/(b+1);

  // testing
  a=a/(n*(daop.mean_y*(1-daop.mean_y)));
  
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


}



void update_terminal_nodes_wrapper(int treeindx)
{

  // For loss functions other than squared error, a single newton step is performed to fit the terminal node fits.
  // This function performs the update.


  // if(daop.method==1) : do nothing. 

  // This function just decides whether to run update_terminal_nodes() function
  // (not needed for regression, but needed for logistic regression, fex)
  
  if(daop.method==2)
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
  //Rprintf(" Finished with node_predictions vector \n");
  //for(i=0;i<daop.n;i++)
  for(k=0;k<ntnodes;k++)
    {

      ni=tnodes[k];
      HASH_FIND_INT(nodes,&ni,s);
      if(s==NULL){}else{

	HASH_ITER(hh,s->rows,r,rtemp){
      //  Rprintf(" node[%d]=%d n=%d \n",i,daop.node[i],daop.n);
	  i=r->row_number;

	  daop.predictions[i]+=(node_predictions[ni]);
	  daop.predictions_response[i]=1/(1+exp(-daop.predictions[i]));
	  daop.response[i]=(daop.x[daop.yindx][i]-daop.predictions_response[i]);

      //k=find_tnode(tree_indx,i);
      //if(k!=daop.node[i])
      //Rprintf(" [tree_indx=%d] tnode=%d and daop.node=%d dont match for %d.  \n",tree_indx,k,daop.node[i],i);

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
