

get_tree<-function(ti,all=TRUE,max_size=20)
{
ncol=15
if(!all){
	h<-.C("get_tree",as.integer(ti),res=double(ncol),n=integer(1))
}else{
 # gets all trees upto and including ti 
	h<-.C("get_tree_all",as.integer(ti),res=double(ncol*(ti+1)*max_size),n=integer(1))

}
m<-matrix(h$res[1:h$n],ncol=ncol,byrow=TRUE)
m
}

set_ncat=function(ncat)
{

	h=.C("set_ncat",as.integer(ncat))
}


get_tree_classify<-function(ti,all=TRUE,max_size=20,ncat)
{
ncol=15+ncat
if(!all){
	#h<-.C("get_tree",as.integer(ti),res=double(ncol),n=integer(1))
}else{
 # gets all trees upto and including ti 
	h<-.C("get_tree_all_gini",as.integer(ti),as.integer(ncat),res=double(ncol*(ti+1)*max_size),n=integer(1))

}
m<-matrix(h$res[1:h$n],ncol=ncol,byrow=TRUE)
m
}

get_info<-function()
{
h<-.C("get_data_options_info",n=integer(1),p=integer(1),nboost=integer(1))
list(n=h$n,p=h$p,B=h$nboost)
}

predict_ltb<-function(ntrees,n)
{
h<-.C("predict_trees",as.integer(ntrees),pred=double(n))
h$pred
}


#find_tnode<-function(ti,ri)
#{
#h<-.C("find_tnode_wrapper",as.integer(ti),as.integer(ri),tnode=integer(1))
#h$tnode
#}

set_data<-function(x,time,id,yindx)
{
#void read_data(double *x,int *n,int *p,double *time,int *id,int *yindx)
n=nrow(x)
p=ncol(x)
h<-.C("set_data",as.double(x),as.integer(n),as.integer(p),as.double(time),as.integer(id),as.integer(yindx))

}


set_train<-function(train)
{
h<-.C("set_train",as.integer(train))
}

set_lambda<-function(lambda)
{
h<-.C("set_lambda",as.double(lambda))
}

free<-function()
{
h<-.C("free_daop")
}


read_in_data<-function(x,time,id,yindx,nboost=1000,lambda=.1,nsplit=1,rf=0,nsamp=5,time_split=1,window_summary=0)
{
#void read_data(double *x,int *n,int *p,double *time,int *id,int *yindx)

n=nrow(x)
p=ncol(x)
h<-.C("read_data",as.double(as.matrix(x)),as.integer(n),as.integer(p),as.double(time),as.integer(id),as.integer(yindx),
	as.double(lambda),as.integer(nsplit),as.integer(nboost),as.integer(rf),as.integer(nsamp),as.integer(time_split),as.integer(window_summary))

}


tsummary_wrapper<-function(n,cut,delta,vindx)
{
#void tsummary_wrapper(double *nrec,double *nrec_condition,double *cut,double *delta,int *vindx)

h<-.C("tsummary_wrapper",nrec=double(n),nrecc=double(n),as.double(cut),as.double(delta),as.integer(vindx))

dd=data.frame(nrec=h$nrec,nrecc=h$nrecc)
dd
}







# --- uthash based functions 

sample_indx=function(n)
{
	h=.C("sample_indx_wrapper",as.integer(n),res=integer(1))

h$res
}

permute_w=function(n)
{
	h=.C("permute_wrapper",as.integer(n),res=integer(n))

h$res

}

add_row<-function(n,r)
{
	h=.C("add_row_wrapper",as.integer(n),as.integer(r))
}

delete_node<-function(n)
{
	h=.C("delete_node_wrapper",as.integer(n))
}

print_node=function(n)
{
	h=.C("print_node_members",as.integer(n))
}


get_oob<-function()
{
	h=get_info()
	m=.C("get_oob",res=integer(h$B*h$n))

	oob=matrix(m$res,ncol=h$B,nrow=h$n,byrow=FALSE)
	oob
}


set_oob<-function(oob)
{
# This function sets the out-of-bag list, and turns off subsampling. 

	h=.C("set_oob",as.integer(oob))
}

set_mtry=function(mtry)
{
	h=.C("set_mtry",as.integer(mtry))
}



###################################################################################

set_variable_status=function(vh,vc)
{
#void set_variable_status(int *vhistory,int *nhistory,int *vconcurrent,int *nconcurrent)se

h=.C("set_variable_status",as.integer(vh),as.integer(length(vh)),as.integer(vc),as.integer(length(vc)))

}





