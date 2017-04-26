



uboost_dev<-function(B,vi=0,oob=FALSE,random_split=0)
{
i=get_info()
n=i$n
p=i$p
h<-.C("boost_uthash_dev",as.numeric(random_split),as.integer(B),testerror=double(B),pred=double(n),as.integer(vi),vimp=double(B*p))

#print(" R: Finished boost_uthash_dev ")

oobmatrix=NULL
if(oob)
{
	# get oob matrix 
	oobmatrix=get_oob()
}

if(vi==0){
	retlist=list(testerror=h$testerror,pred=h$pred,oob=oobmatrix)
}else{
	m=.C("delete_oob")
	retlist=list(testerror=h$testerror,pred=h$pred,vimp=matrix(h$vimp,ncol=p,nrow=B,byrow=FALSE),oob=oobmatrix)
}

retlist
}






uboost_dev_wrapper<-function(B,vi=0,oob=FALSE,random_split=0,method=1)
{
i=get_info()
n=i$n
p=i$p
h<-.C("boost_wrapper",as.integer(random_split),as.integer(method),as.integer(B),testerror=double(B),pred=double(n),as.integer(vi),vimp=double(B*p))

#print(" R: Finished boost_uthash_dev ")

oobmatrix=NULL
if(oob)
{
	# get oob matrix 
	oobmatrix=get_oob()
}

if(vi==0){
	retlist=list(testerror=h$testerror,pred=h$pred,oob=oobmatrix)
}else{
	m=.C("delete_oob")
	retlist=list(testerror=h$testerror,pred=h$pred,vimp=matrix(h$vimp,ncol=p,nrow=B,byrow=FALSE),oob=oobmatrix)
}

retlist
}



htree<-function(x,time,id,yindx,ntrees=100,lambda=1,rf=1,nsplit=200,nsamp=5,tf=.5,id_sampling=TRUE,rsplit=0,
		mtry=NULL,vi=0,time_split=1,oob=FALSE,oobmatrix=NULL,keep_data=TRUE,vh=NULL,vc=NULL,cv_train=NULL,method=1)
{
read_in_data(x=x,time=time,id=id,yindx=yindx,nboost=ntrees,lambda=lambda,nsplit=nsplit,rf=rf,nsamp=nsamp,time_split=time_split)

if(is.null(mtry))
{
	mtry=max(c(floor(ncol(x)/3),1))
}

set_mtry(mtry=mtry)


if(!is.null(vh))
{
if(is.null(vc))
	stop(" Both 'vh' and 'vc' need to be specified ....")
set_variable_status(vh=vh,vc=vc)
}


train=NULL
if(rf!=1)
{

if(is.null(cv_train))
{

	if(id_sampling){
	uid<-unique(id)
	tr_id<-sample(uid,round(length(uid)*tf))
	train<-rep(0,nrow(x))
	train[is.element(id,tr_id)]=1
	}else{
	train<-rep(0,nrow(x))
	train[sample(nrow(x),round(nrow(x)*tf))]=1
	}
}else{
	train=cv_train
}
	set_train(train=train)
}

if(!is.null(oobmatrix))
	set_oob(oobmatrix)

#print(1)
#h=uboost_dev(B=ntrees,vi=vi,oob=oob,random_split=rsplit)
h=uboost_dev_wrapper(B=ntrees,vi=vi,oob=oob,random_split=rsplit,method=method)

#print(" getting tree ")
m=get_tree(ti=(ntrees-1),all=TRUE,max_size=(nsplit*4))
#print(" freeing data")

######################## TESTING ############################################
#hh=.C("get_split_matrix",res=double(1000*13))
#split_matrix=matrix(hh$res,ncol=13,nrow=1000,byrow=TRUE)
######################## DONE: TESTING ############################################


free()
#print(" finished with data-free")



	retlist=list(trees=m,error=h$testerror,pred=h$pred,yindx=yindx,
		time=time,rf=rf,mtry=mtry,nsamp=nsamp,
		nsplit=nsplit,rf=rf,lambda=lambda,ntrees=ntrees,vimp=h$vimp,
		time_split=time_split,oob=h$oob,train=train,vi=vi,oob=oob,vh=vh,vc=vc,method=method,id_sampling=id_sampling)

if(keep_data){
	retlist$x=x
	retlist$id=id
	}

class(retlist)<-"htree"
retlist
}



seboot=function(B,M=10,object)
{
m=object

# start end for id (assume time sorted, within id) 

ii=diff(m$id)
start=1
end=length(m$id)
jj=c(2:length(m$id))[ii!=0]
start=c(start,jj)
end=c(jj-1,end)
uid=unique(m$id)
iList=list()
nList=list()

for(k in 1:length(uid))
{
	iList[[as.character(uid[k])]]=c(start[k]:end[k])
	nList[[as.character(uid[k])]]=end[k]-start[k]+1
}

mList=list()
for(b in 1:B)
{
boot_id=sample(uid,length(uid),replace=T)

new_id=NULL
for(j in 1:length(boot_id))
	new_id=c(new_id,rep(j,nList[[as.character(boot_id[j])]]))

ii=as.numeric(unlist(iList[as.character(boot_id)]))

x_boot=m$x[ii,]
id_boot=new_id
time_boot=m$time[ii]
hboot=htree(x=x_boot,time=time_boot,id=id_boot,yindx=m$yindx,ntrees=M,lambda=1,rf=1,nsplit=m$nsplit,nsamp=m$nsamp,tf=1,id_sampling=TRUE,rsplit=m$rsplit,
		mtry=m$mtry,vi=0,time_split=m$time_split,oob=FALSE,oobmatrix=NULL,keep_data=FALSE,vh=m$vh,vc=m$vc,cv_train=NULL,method=m$method)
mList[[b]]=list(bag=boot_id,fit=hboot)

}

mList
}



predict_se=function(m,object,x=NULL,id=NULL,time=NULL)
{

if(is.null(x))
{
	x=object$x
	id=object$id
	time=object$time

}

for(b in 1:length(m))
{
	pp=predict_htree(object=m[[b]]$fit,x=x,id=id,time=time,yindx=object$yindx,all.trees=TRUE)
	if(b==1){
		pa=apply(pp,1,mean)
		sa=apply(pp,1,var)
	}else{
		pa=cbind(pa,apply(pp,1,mean))
		sa=sa+apply(pp,1,var)
	}
}

B=length(m)
M=object$ntrees
R=m[[1]]$fit$ntrees
bias=sa/B*(1/R-1/M)
var_biased=apply(pa,1,var)
var_hat=var_biased-bias
mean_var=mean(var_hat,na.rm=T)
var_hat[var_hat<0]=mean_var
se=sqrt(var_hat)
se
}





predict_htree<-function(object,x=NULL,yindx=NULL,time=NULL,id=NULL,ntrees=NULL,time_split=1,all.trees=FALSE)
{

if(is.null(ntrees))
	ntrees=length(unique(object$trees[,1]))

ntrees_in_matrix=length(unique(object$trees[,1]))

if(is.null(x))
{
	x=object$x
	yindx=object$yindx
	time=object$time
	id=object$id
}

n=nrow(x)
p=ncol(x)

h=.C("read_predict",as.double(as.matrix(x)),as.integer(n),as.integer(p),as.double(time),as.integer(id),
		as.integer(yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),as.integer(ntrees_in_matrix),as.integer(time_split))

if(!all.trees){
# --- predict 
	h=.C("predict_trees_fast",as.integer(ntrees),pred=double(nrow(x)))

predictions=h$pred
if(object$rf==1)
	predictions=predictions/ntrees
}else{
	h=.C("predict_trees_all",as.integer(ntrees),pred=double(nrow(x)*ntrees))

	predictions=matrix(h$pred,ncol=ntrees,nrow=nrow(x),byrow=FALSE)

}

# -- free pred info  
h=.C("free_predict")


predictions
}



if(0)
{
p=partdep(object=h,xindx=4,xlim=NULL,ngrid=100,ntrees=NULL,subsample=1)


}

partdep=function(object,xindx,xlim=NULL,ngrid=100,ntrees=NULL,subsample=1)
{

# object=h;xindx=4;subsample=.1

uid=unique(object$id)
ns=round(subsample*length(uid))
sid=sample(uid,size=ns,replace=F)

ii=c(1:length(object$id))[is.element(object$id,sid)]
nii=length(ii)
mid=max(object$id)+1
idaux=object$id[ii]
id=NULL
for(k in 1:ngrid)
id=c(id,idaux+mid*k)

ii=rep(ii,ngrid)
x=object$x[ii,]
time=object$time[ii]



if(is.null(ntrees))
	ntrees=object$ntrees

fit=object  #$fit

uv=unique(object$x[,xindx])
nuv=length(uv)
if(nuv<ngrid)
{
	tt=sort(uv)
}else{
	if(is.null(xlim)){
		tt=seq(min(x[,xindx]),max(x[,xindx]),length=ngrid)
	}else{
		tt=seq(xlim[1],xlim[2],length=ngrid)
	}
}
 ### working here .... 
tt=sort(rep(tt,nii)) #expand.grid(tt,rep(1,nii))[,1]

x[,xindx]=tt

vv=NULL
xx=as.matrix(x)
yindx=fit$yindx

hh=predict_htree(object=fit,x=x,yindx=yindx,time=time,id=id,ntrees=ntrees,time_split=fit$time_split)

dd=data.frame(y=hh,tt=tt)
ds=aggregate(dd$y,list(tt),mean)


list(y=ds[,2],x=ds[,1])
}




partdep_se=function(object,m,xindx,xlim=NULL,ngrid=100,ntrees=NULL,subsample=.1)
{


	uv=unique(object$x[,xindx])
	nuv=length(uv)	
	ngrid_use=ngrid
	if(ngrid>nuv)
	{
		ngrid_use=nuv
	}


uid=unique(object$id)
ns=round(subsample*length(uid))
sid=sample(uid,size=ns,replace=F)

ii=c(1:length(object$id))[is.element(object$id,sid)]
nii=length(ii)
mid=max(object$id)+1
idaux=object$id[ii]
id=NULL
for(k in 1:ngrid_use)
id=c(id,idaux+mid*k)

ii=rep(ii,ngrid_use)
x=object$x[ii,]
time=object$time[ii]



if(is.null(ntrees))
	ntrees=object$ntrees

fit=object  #$fit
	if(nuv<ngrid)
	{
		tt=sort(uv)
	}else{

if(is.null(xlim)){
	tt=seq(min(x[,xindx]),max(x[,xindx]),length=ngrid)
}else{
	tt=seq(xlim[1],xlim[2],length=ngrid)
}
}
tt_orig=tt
 ### working here .... 
tt=sort(rep(tt,nii)) #expand.grid(tt,rep(1,nii))[,1]

x[,xindx]=tt

vv=NULL
xx=as.matrix(x)
yindx=fit$yindx
B=length(m)
bias=rep(0,ngrid_use)
sum_est=rep(0,ngrid_use)
ss=rep(0,ngrid_use)
hh=predict_htree(object=object,x=x,yindx=yindx,time=time,id=id,ntrees=object$ntrees,time_split=object$time_split,all.trees=FALSE)
ds=aggregate(hh,list(tt),mean)
pd_fit=ds[,-1]

for(b in 1:B){
hh=predict_htree(object=m[[b]]$fit,x=x,yindx=yindx,time=time,id=id,ntrees=m[[b]]$fit$ntrees,time_split=m[[b]]$fit$time_split,all.trees=TRUE)

dd=data.frame(y=hh,tt=tt)
ds=aggregate(hh,list(tt),mean)

bias_aux=apply(ds[,-1],1,var)
bias=bias+bias_aux
est_b=apply(ds[,-1],1,mean)
ss=ss+est_b^2
sum_est=sum_est+est_b

}
R=m[[1]]$fit$ntrees
bias=bias*(1/R-1/object$ntrees)/B
est_d=sum_est/B
v_bias=ss/B-(sum_est/B)^2
var_est=v_bias-bias
var_est[var_est<0]=median(var_est)
se=sqrt(var_est)
list(y=pd_fit,se=se,x=tt_orig)
}















