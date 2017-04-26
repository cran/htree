
# -------------------------------- #
# hrf: historical random forest
# -------------------------------- #

hrf=function(x,time,id,yindx,ntrees=100,mtry=NULL,se=FALSE,B=100,R=10,nsplit=NULL,nsamp=5,rsplit=FALSE,
		historical=TRUE,keep_data=TRUE,vh=NULL,vc=NULL,order.cat=FALSE)
{


 ## Check missing values
  if (any(is.na(x))) {
    stop("Missing data in 'x'.", call. = FALSE)
  }
  if (any(is.na(id))) {
    stop("Missing data in 'id'.", call. = FALSE)
  }
  if (any(is.na(time))) {
    stop("Missing data in 'time'.", call. = FALSE)
  }
  
	if(!is.data.frame(x))
		x=as.data.frame(x)

	tiny_number=.00001
	lambda=1
	id_sampling=TRUE
	rf=1

	if(is.null(nsplit))
		nsplit=nrow(x)

	tf=1
	cv_train=NULL
	method=1
	oobmatrix=NULL
	varimp=FALSE	
	vi=as.numeric(varimp)
	rsplit=as.numeric(rsplit)
	time_split=as.numeric(historical)
	oob=TRUE

	## --- set up historical and concurrent predictor vectors: vh and vc ----------  ##
	if(!is.null(vh))
	{

		if(is.null(vc))
	    		vc=c(1:ncol(x))[-yindx] #[-vh]
		vc=vc-1
		vh=vh-1

	}else{
	    # identify time-varying predictors 
	    xaux=aggregate(x,by=list(id=id),sd)[,-1]
	    xaux[is.na(xaux)]=0
	    sd_x_id=apply(xaux,2,sd)
	   # print(sd_x_id)
	    vh=c(1:ncol(x))[sd_x_id>tiny_number]
	    if(length(vh)==0)
		{
			print(" No across time variation within levels of 'id'. ")
			time_split=0
			vh=NULL
			vc=NULL
		}else{
			if(is.null(vc))
				vc=c(1:ncol(x))[-yindx] # [-vh]
	    		vc=vc-1
	   		vh=vh-1
		}
	
	}

	## cat_predictors goes here .... 
	##
	mm="unordered"
	if(order.cat)
		mm="ordered"
	cc=cat_predictors(x=x,yindx=yindx,method=mm)
 

	fit_rf=htree(x=cc$x,time=time,id=id,yindx=(yindx-1),ntrees=ntrees,lambda=lambda,rf=rf,nsplit=nsplit,nsamp=nsamp,tf=tf,id_sampling=id_sampling,rsplit=rsplit,
		mtry=mtry,vi=vi,time_split=time_split,oob=oob,oobmatrix=oobmatrix,keep_data=keep_data,vh=vh,vc=vc,cv_train=cv_train,method=method)

	# mapping of categorical variables to ordered variables (if any)
	fit_rf$mapping_categorical=cc$mapping

	fit_boot=NULL
	if(se)
	{
		# Compute standard errors (noisy bootstrap) 
		fit_boot=seboot(B=B,M=R,object=fit_rf)
	
	}

	fit_rf$boot=fit_boot
	fit_rf
}




cat_predictors=function(x,yindx,method="ordered")
{
	# method="one-hot"
	# 	one-hot encoding is used on all categorical variables  (beware of many categories!!!) 
	# method="ordered"
	#	ordered conversion. That is, each categorical predictor is converted into an ordered variable 
	#	
	# make sure you dont muck up yindx .... 

	#method="ordered" ## if you want one-hot then you need to do it yourself ....
	if(!is.element(method,c("ordered","unordered")))
		stop(" 'method' must be one of 'ordered' or 'unordered'.",.call=FALSE)	
	
	findx=c(1:ncol(x))[unlist(lapply(x,is.factor))]

	xc=data.frame(x)
	# xc=xc[,c(yindx,c(1:ncol(xc))[-yindx])]

	# convert factors to strings 
	if(length(findx)>0)
	{
		for(k in 1:length(findx))
			xc[,findx[k]]=as.character(xc[,findx[k]])	
	}

	y=xc[,yindx]
	#xc=xc[,-yindx]
	retList=list(mapping=NULL,x=xc)
	cindx=c(1:ncol(xc))[unlist(lapply(xc,is.character))]
	clist=list()  # 
	
	if(length(cindx)>0)
	{
		
		# there are categorical predictors 

		nx=names(xc)
		if(method=="one-hot")
		{

		# one-hot encoding 
	 	vm=NULL
		vm_names=NULL
		clist=list()
		for(k in 1:length(cindx))
		{
			xx=xc[,cindx[k]]
			ux=unique(xx)
			vv=NULL
			vv_names=paste(nx[cindx[k]],ux,sep="_")
			for(j in 1:length(ux))
			{
				vv=cbind(vv,as.numeric(xx==ux[k]))

			}
			vm=cbind(vm,vv)
			vm_names=c(vm_names,vv_names)
			clist[[nx[cindx[k]]]]=ux
		}

		xnew=xc[,-cindx]
		nxnew=names(xnew)
		xnew=data.frame(xnew,vm)
		names(xnew)[(length(nxnew)+1):ncol(xnew)]=vm_names
		retList=list(method="one-hot",mapping=clist,x=xnew,y=y)
		}
		if(method=="ordered")
		{
		# order conversion, factors/string variables are replaced by ordered variables according to response mean 

		if(length(cindx)>0)	
		{
		y=xc[,yindx]
		mapping=list()
		nx=names(xc)
		for(k in 1:length(cindx))
		{
			xx=xc[,cindx[k]]
			aa=aggregate(y,by=list(ff=xx),mean,na.rm=TRUE)
			aa=aa[order(aa[,2]),]

			ml=list()
			ml[aa[,1]]=c(1:nrow(aa))
			xnew=as.numeric(unlist(ml[xc[,cindx[k]]]))
			mapping[[nx[cindx[k]]]]=ml
			xc[,cindx[k]]=xnew

		}
		retList=list(method="ordered",x=xc,mapping=mapping)
		}
		}
		
		if(method=="unordered")
		{
		# arbitrary ordering of categorical variables  
		# .... each categorical is mapped into an ordered variable (based simply on order of occurence in predictor variable ) ....
 
		if(length(cindx)>0)	
		{
		y=xc[,yindx]
		mapping=list()
		nx=names(xc)
		for(k in 1:length(cindx))
		{
			xx=xc[,cindx[k]]
			ux=unique(xx)

			ml=list()
			ml[ux]=c(1:length(ux))
			xnew=as.numeric(unlist(ml[xc[,cindx[k]]]))
			mapping[[nx[cindx[k]]]]=ml
			xc[,cindx[k]]=xnew

		}
		retList=list(method="unordered",x=xc,mapping=mapping)
		}
		}
	}
	
retList
}


map_categorical=function(object,x)
{
	# Return x with categorical variables (if any) mapped to integers 
	xnew=x
	if(length(object$mapping)>0)
	{
		xnew=data.frame(x)
		nx=names(xnew)
		nc=names(object$mapping)
		for(k in 1:length(object$mapping))
		{
			xk=as.numeric(unlist(object$mapping[[nc[k]]][xnew[[nc[k]]]]))
			xnew[[nc[k]]]=xk
		}

	}

	return(xnew)
}






predict_hrf=function(object,x=NULL,yindx=NULL,time=NULL,id=NULL,all.trees=FALSE,se=FALSE)
{

	if(!is.null(yindx))
		yindx=yindx-1

	if(!is.null(x))
	{
		if(length(object$mapping)>0)
		{
			# Map categorical variables to integers 
			x=map_categorical(object=object,x=x)
		}

	}

	# predict_se=function(m,object,x=NULL,id=NULL,time=NULL)
	pred=predict_htree(object=object,x=x,yindx=yindx,time=time,id=id,ntrees=NULL,time_split=1,all.trees=all.trees)
	retlist=pred
	if(se)
	{
		if(!is.null(object$boot))
		{
			serr=predict_se(m=object$boot,object=object,x=x,id=id,time=time)
			retlist=cbind(pred,serr)
			colnames(retlist)=c("pred","se")
		}else{
			print(" No bootstrap replicates found, run 'hrf' with 'se=TRUE'.")
		}
	}
	
	retlist
}



partdep_hrf=function(object,xindx,xlim=NULL,ngrid=25,subsample=.1,plot.it=TRUE,cat.plot=FALSE)
{
	if(is.null(object$boot))
	{
		# No bootstrap samples, no standard error bars produced
		pd=partdep(object=object,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=NULL,subsample=subsample)

	}else{
		# partial dependence with standard errors 
		pd=partdep_se(object=object,m=object$boot,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=NULL,subsample=subsample)
	}

	if(plot.it)
	{
		nx=names(object$x)[xindx]
		if(!is.null(pd$se))
		{
			if(!cat.plot){
	
				upper=pd$y+2*pd$se;lower=pd$y-2*pd$se
				ylim=c(min(lower),max(upper))
				plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main="Partial dependence",ylim=ylim,cex.axis=1.5,cex.lab=2,cex.main=2)
				points(pd$x,upper,type="l",lwd=3,lty=2)
				points(pd$x,lower,type="l",lwd=3,lty=2)
			}else{

				plot_cat2(m=pd$y,l=(pd$y-2*pd$se),u=(pd$y+2*pd$se),cc=pd$x,ylim=NULL,mult=1,bf=4,lwd=3,main="Partial dependence")


			}
	
		}else{
			plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main="Partial dependence",cex.axis=1.5,cex.lab=2,cex.main=2)
		}
	}
pd
}





varimp_hrf=function(object,nperm=10)
{


	ntrees=length(unique(object$trees[,1]))
	ntrees_in_matrix=length(unique(object$trees[,1]))

	group="id"
	x=object$x
	yindx=object$yindx
	time=object$time
	id=object$id
	time_split=1
	y=x[,yindx+1]


	n=nrow(x)
	p=ncol(x)

	h=.C("read_predict",as.double(as.matrix(x)),as.integer(n),as.integer(p),as.double(time),as.integer(id),
		as.integer(yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),as.integer(ntrees_in_matrix),as.integer(time_split))
	
	# 
	set_oob(object$oob)

	h=.C("predict_trees_all",as.integer(ntrees),pred=double(nrow(x)*ntrees))	
	predictions=matrix(h$pred,ncol=ntrees,nrow=nrow(x),byrow=FALSE)
	noob=apply(object$oob==1,1,sum)
	pred_oob=apply(predictions*object$oob,1,sum)
	pred_oob=pred_oob/noob

 	#void varimp(int *nperm,double *res)
	h=.C("varimp",as.integer(nperm),res=double(n*p))
	perm_pred=matrix(h$res,ncol=p,nrow=n,byrow=FALSE)

	# p=5;id=hb$id;pred_oob=hz$pred_oob;perm_pred=hz$pred_perm
	zscore=NULL
	for(k in 1:p)
	{
		z=abs(y-perm_pred[,k])-abs(y-pred_oob)
		if(group=="id")
		{	# compute error change by subject/id 
			zid=aggregate(z,by=list(id=id),mean) 
			se_zid=sqrt(var(zid[,2])/length(unique(id)))
			zscore=c(zscore,mean(zid[,2])/se_zid)
		}else{
			# error change by row, assuming independence between 
			se_z=sqrt(var(z)/n)
			zscore=c(zscore,mean(z)/se_z)

		}
	}

	names(zscore)=colnames(x)

	# -- free pred info  
	h=.C("free_predict")



	#list(zscore=zscore,pred_oob=pred_oob,pred_perm=perm_pred)
	zscore
}






















quantile_hrf<-function(object,x=NULL,yindx=NULL,time=NULL,id=NULL,ntrees=NULL,prob=seq(.1,.9,length=10))
{

time_split=1

if(is.null(ntrees))
	ntrees=length(unique(object$trees[,1]))

ntrees_in_matrix=length(unique(object$trees[,1]))


n=nrow(object$x)
p=ncol(object$x)

h=.C("read_predict",as.double(as.matrix(object$x)),as.integer(n),as.integer(p),as.double(object$time),as.integer(object$id),
		as.integer(object$yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),as.integer(ntrees_in_matrix),as.integer(time_split))

h=.C("set_tnode_row")



# -- free pred info  
h=.C("free_predict")



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


h=.C("get_weights",as.integer(nrow(object$x)),res=integer(nrow(object$x)*nrow(x)))

ii=order(object$x[,yindx+1])
ysorted=object$x[ii,yindx+1]

weights=matrix(h$res,nrow=nrow(x),ncol=nrow(object$x),byrow=TRUE)
weights=weights[,ii]
weights=weights/apply(weights,1,sum)

# -- free pred info  
h=.C("free_predict")

np=length(prob)
n=length(ysorted)
m=nrow(x)
h=.C("quantile_R",as.double(prob),as.integer(np),as.double(ysorted),as.integer(n),as.double(as.numeric(t(weights))),as.integer(m),quant=double(m*np))
quants=matrix(h$quant,nrow=m,ncol=np,byrow=TRUE)

#list(weights=weights,y=ysorted,)
quants
}


quant_test=function(q,w,y)
{

h=.C("quantile_aux",as.double(q),as.double(y),as.integer(length(y)),as.double(w),quant=double(1))

h$quant
}


## TESTING FUNCTION 

get_tnode=function(ti,ni)
{

#void get_rows(int *treeindx,int *nodeindx,int *res,int *n)
h=.C("get_rows",as.integer(ti),as.integer(ni),res=integer(1000),n=integer(1))
h$res[1:h$n]
}

