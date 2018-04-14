# -------------------------------- #
# htb: historical tree boosting
# -------------------------------- #

htb=function(x,time=NULL,id=NULL,yindx,ntrees=100,method="freqw",nsplit=1,lambda=.1,family="gaussian",cv.fold=0,cv.rep=NULL,nsamp=5,
		historical=TRUE,keep_data=TRUE,vh=NULL,vc=NULL,delta=NULL)
{

	rsplit=FALSE

	if((is.null(id)|is.null(time)))
	{
		## if time/id then data assumed iid at same time point (ie standard random forest)
		id=c(1:nrow(x))
		time=rep(1,nrow(x))
	}
	

	if(!is.element(class(yindx),c("numeric","character")))
	{
		stop(" 'yindx' must give column number, or name, of response in 'x'.") 
	}else{
		if(class(yindx)=="character")
			yindx=which(names(x)==yindx)

		if(length(yindx)==0)
			stop(" 'yindx' not found ")

	}


  if (any(is.na(x))) {
    stop("Missing data in 'x'.", call. = FALSE)
  }
  if (any(is.na(id))) {
    stop("Missing data in 'id'.", call. = FALSE)
  }
  if (any(is.na(time))) {
    stop("Missing data in 'time'.", call. = FALSE)
  }
  

	# method code
	method=mcode(m=method)

	if(!is.data.frame(x))
		x=as.data.frame(x)

 	if(yindx>ncol(x)|yindx<1)
		stop(" 'yindx' out of range")

 	if(!is.numeric(time)){
		stop("'time' variable must be numeric. ")
	}  


	vtype=unlist(lapply(x,class))

	if(is.character(id))
		id=as.factor(id)

	if(is.factor(id))
		id=as.numeric(id)




	if(length(unique(c(nrow(x),length(id),length(time))))!=1) 
		stop(" length of 'id', 'time' and number of rows of 'x' must be the same.")

	classInfo=NULL

	rd=reformat_data(x=x)
	x=rd$x

	ii=order(id,time)
	x=x[ii,]
	id=id[ii]
	time=time[ii]

	flist=rd$flist


	tiny_number=.00001
	
	mtry=ncol(x)
	rf=0
	id_sampling=TRUE
	varimp=FALSE
	oobmatrix=NULL

	if(is.element(family,c("gaussian","bernoulli")))
	{
		if(family=="gaussian")
			family_code=1
		if(family=="bernoulli")
			family_code=2

	}else{
		# -- 
		family="gaussian"
		family_code=1
		print(" 'family' not specified, set to 'gaussian' ")
	}

	if(is.null(nsplit))
		nsplit=nrow(x)

	tf=1
	cv_train=NULL
	method_family=family_code

	vi=as.numeric(varimp)
	rsplit=as.numeric(rsplit)
	time_split=as.numeric(historical)
	oob=TRUE

	## --- set up historical and concurrent predictor vectors: vh and vc ----------  ##
	hc=histcon_aux(vh=vh,vc=vc,x=x,id=id)	
	vc=hc$vc
	vh=hc$vh


	# -- train model 
	fit_gb=htree(x=x,time=time,id=id,yindx=(yindx-1),ntrees=ntrees,lambda=lambda,rf=rf,nsplit=nsplit,nsamp=nsamp,tf=tf,id_sampling=id_sampling,rsplit=rsplit,
		mtry=mtry,vi=vi,time_split=time_split,oob=oob,oobmatrix=oobmatrix,keep_data=keep_data,vh=vh,vc=vc,cv_train=cv_train,method=method,
		delta=delta,control=list(method_family=method_family))

	fit_cv=list()
	if(cv.fold>1)
	{
		tf_cv=(1-1/cv.fold)
		
		if(is.null(cv.rep))
			cv.rep=cv.fold
		
		cv_error=rep(0,ntrees)
		for(k in 1:cv.rep)
		{
			#cat(paste(" cv-fold ",k," out-of ",cv.rep,".\n",sep=""))
			fit_cv[[k]]=htree(x=x,time=time,id=id,yindx=(yindx-1),ntrees=ntrees,lambda=lambda,rf=rf,
					nsplit=nsplit,nsamp=nsamp,tf=tf_cv,id_sampling=id_sampling,rsplit=rsplit,
					mtry=mtry,vi=vi,time_split=time_split,oob=oob,oobmatrix=oobmatrix,keep_data=FALSE,
					vh=vh,vc=vc,cv_train=cv_train,method=method,control=list(method_family=method_family))

			cv_error=cv_error+fit_cv[[k]]$error
		}
	
		
		fit_gb$cv_error=cv_error/cv.rep
		if(family=="bernoulli")
			fit_gb$cv_error=-fit_gb$cv_error
	}

	fit_gb$vtype=vtype
	fit_gb$indx_original=ii
	fit_gb$flist=flist
	fit_gb$family=family
	fit_gb$method_family=method_family

	fit_gb$cv.rep=cv.rep
	fit_gb$cv.fold=cv.fold
	fit_gb$cv=fit_cv
	fit_gb
}


varimp_htb=function(object,nperm=20,ntrees=NULL)
{
# object=ff;ntrees=100;nperm=10
	group="id";stat=1
	sumstat=stat
	if(length(object$cv)==0)
		stop(" 'object' does not contain cross-validation runs. No variable importance produced.")


	if(is.null(ntrees))
	{
		# select cross-validation error minimizing iteration
		ntrees=order(object$cv_error)[1]
		cat(paste("'ntrees' not specificed. Using ",ntrees," for prediction. \n",sep=""))

	}

	# Get out-of-sample predictions (non-permuted and permuted) 
	pred_orig=rep(0,nrow(object$x))
	pred_permuted=array(0,dim=dim(object$x))
	nout=rep(0,nrow(object$x))
	for(k in 1:length(object$cv))
	{

		train_k=object$cv[[k]]$train
		ii_out=c(1:length(train_k))[train_k==0]  # predict observations indexes 
		x_test=object$x[ii_out,];time_test=object$time[ii_out];id_test=object$id[ii_out];yindx=object$yindx
		vi=varimp_aux(object=object$cv[[k]],x=x_test,time=time_test,id=id_test,nperm=nperm,ntrees=ntrees)
		pred_orig[ii_out]=pred_orig[ii_out]+vi$pred
		pred_permuted[ii_out,]=pred_permuted[ii_out,]+vi$perm_pred
		nout[ii_out]=nout[ii_out]+1  # Note: since out-of-sample set is randomly sampled could have repetitions and non-inclusion.. 

	}

	ii_keep=c(1:length(nout))[nout>0]
	nout=nout[ii_keep]
	pred_orig=pred_orig[ii_keep]
	pred_orig=pred_orig/nout
	pred_permuted=pred_permuted[ii_keep,]
	pred_permuted=pred_permuted/nout

	# -- prediction error change ---- 
	method=object$method_family
	id=object$id[ii_keep]
	y=object$x[ii_keep,object$yindx+1]
	zscore=NULL
	pchange=NULL
	sdchange=NULL

	for(k in 1:ncol(object$x))
	{
			delta=pred_error_change(y=y,pred=pred_orig,pred_perm=pred_permuted[,k],method=method,sumstat=sumstat)
			delta_id=aggregate(delta,by=list(id=id),mean)[,2]
			se_delta=sqrt(var(delta_id)/length(delta_id))
			zscore=c(zscore,mean(delta_id)/se_delta)
			pchange=c(pchange,mean(delta_id))
			sdchange=c(sdchange,se_delta)
			

	}
	
	names(zscore)=colnames(object$x)
	err=min(object$cv_error)
	dz=data.frame(delta_rel=round(pchange/err,3),delta_error=round(pchange,3),stderr=round(sdchange,3),zscore=round(zscore,3),pvalue=pnorm(zscore,lower.tail=FALSE))
	names(dz)=c("Relative change","Mean change","SE","Z-value","P-value") 
	#rownames(dz)=NULL

	if(object$time_split==0)
		dz=dz[-which(rownames(dz)==names(object$x)[object$yindx+1]),,drop=FALSE]

	#names(zscore)=colnames(object$x)
	#zscore
	dz[order(dz[,4],decreasing=TRUE),]

}


pred_error_change=function(y,pred,pred_perm,method,sumstat=1)
{

	
	if(method==1)
	{
		# Regression
		z=abs(y-pred_perm)-abs(y-pred)

	}


	if(method==2)
	{
		if(sumstat==1)
		z=abs(y-pred_perm)-abs(y-pred)

		# -- below gave wierd answers ... 
		# Logistic regression
		#p=1/(1+exp(-pred))
		#p_perm=1/(1+exp(-pred_perm)) 
		p_perm=pred_perm
		p=pred
		z1=-(y*log(p)+(1-y)*log((1-p)))
		z2=-(y*log(p_perm)+(1-y)*log(1-p_perm))
		if(sumstat==2)
		z=sqrt(z2)-sqrt(z1)

		if(sumstat==3)
		z=(z2)-(z1)

	}

	z
}




varimp_aux=function(object,x,time,id,nperm=10,ntrees)
{

	# Helper function for varimp_htb
	# Returns prediction and average permuted prediction. (pred and perm_pred, resp.)

	
	ntrees_in_matrix=length(unique(object$trees[,1]))
	if(ntrees_in_matrix<ntrees)
		cat(" Fewer number of trees than specified by 'ntrees' argument. \n")

	yindx=object$yindx
	
	time_split=1
	y=x[,yindx+1]


	n=nrow(x)
	p=ncol(x)

	h=.C("read_predict",as.double(as.matrix(x)),as.integer(n),as.integer(p),as.double(time),as.integer(id),
		as.integer(yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),as.integer(ncol(object$trees)),as.integer(ntrees_in_matrix),
		as.integer(time_split),as.integer(object$method))
	
	# 
	
	h=.C("predict_trees_fast",as.integer(ntrees),pred=double(nrow(x)))	
	pred=h$pred

 	#void varimp(int *nperm,double *res)
	h=.C("varimp_boost",as.integer(nperm),res=double(n*p))
	perm_pred=matrix(h$res,ncol=p,nrow=n,byrow=FALSE)
	
	
	h=.C("free_predict")

	if(length(object$control$method_family)>0){
	if(object$control$method_family==2)
	{
		# logistic regression 
		pred=1/(1+exp(-pred))
		perm_pred=1/(1+exp(-perm_pred))
	}
	}
	#list(zscore=zscore,pred_oob=pred_oob,pred_perm=perm_pred)
	list(pred=pred,perm_pred=perm_pred)
}



partdep_htb=function(object,xindx,xlim=NULL,ngrid=25,subsample=.1,ntrees=NULL,plot.it=TRUE,cat.plot=FALSE)
{
	if(length(object$cv)==0)
	{
		# No cross-validation runs, no standard error bars produced
		pd=partdep(object=object,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=ntrees,subsample=subsample)

	}else{
		# partial dependence with standard errors 
		pd=partdep_seb(object=object,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=ntrees,subsample=subsample)
	}
	
	if((subsample>1)|(subsample<=0))
		stop(" Invalid 'subsample' value.")


	if(class(xindx)=="character")
	{
		xindx=which(names(object$x)==xindx)
	
		if(length(xindx)==0)
			stop(" 'xindx' not found.")
	}

	if(xindx>object$dimx[2])
		stop(" Invalid 'xindx' ")


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
			plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main="Partial dependence",cex.axis=1.5,cex.lab=2,cex.main=2,nx=nx)
		}
	}

pd
}


partdep_seb=function(object,xindx,xlim=NULL,ngrid=100,ntrees=NULL,subsample=.1)
{
	# compute approximate standard errors of partial dependence using leave-out-m jackknife 

	# object=ff;xindx=2;xlim=NULL;ngrid=25;ntrees=150;subsample=.1

	if(length(object$cv)==0)
		stop(" No cross-validation runs found. Standard errors cannot be computed.")


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
  
	tt=sort(rep(tt,nii)) #expand.grid(tt,rep(1,nii))[,1]

	x[,xindx]=tt

	vv=NULL
	xx=as.matrix(x)
	yindx=fit$yindx
	B=length(fit$cv)

	sum_est=rep(0,ngrid_use)
	ss=rep(0,ngrid_use)
	hh=predict_htree(object=object,x=x,yindx=yindx,time=time,id=id,ntrees=object$ntrees,time_split=object$time_split,all.trees=FALSE)
	if(object$method_family==2)
	{
		# logistic regression 
		hh=1/(1+exp(-hh))
	}
	ds=aggregate(hh,list(tt),mean)
	pd_fit=ds[,-1]

	for(b in 1:B){
		hh=predict_htree(object=fit$cv[[b]],x=x,yindx=yindx,time=time,id=id,ntrees=ntrees,time_split=fit$time_split,all.trees=FALSE)
		if(object$method_family==2)
		{
			# logistic regression 
			hh=1/(1+exp(-hh))
		}

		dd=data.frame(y=hh,tt=tt)
		ds=aggregate(hh,list(tt),mean)
		est_b=ds[,-1]
		ss=ss+est_b^2
		sum_est=sum_est+est_b

	}

	est_d=sum_est/B
	v_jack=ss/B-(sum_est/B)^2
	v_jack=v_jack*(fit$cv.fold-1) # if cv.rep!=cv.fold then B isn't equal to cv.fold


	se=sqrt(v_jack)
	list(y=pd_fit,se=se,x=tt_orig)
}





plot_cat2=function(m,l,u,cc,ylim=NULL,mult=1,bf=4,lwd=3,main="",nx="")
{
# Plot confidence intervals = m+/-2*se with categorical labels cc
#
#
#cc=rownames(t1)
#m=as.numeric(t1[,1])
#se=as.numeric(t1[,2])
upper=max(u)
lower=min(l)

if(is.null(ylim))
ylim=c(lower,upper)

x=seq(0,1,length=length(cc))
dd=mean(diff(x))

xlim=c(min(x-dd/bf),max(x+dd/bf))
plot(x,m, xaxt="n", xlab=nx, ylab="", ylim=ylim,lwd=lwd,main=main,xlim=xlim,cex.axis=1.5,cex.lab=2,cex.main=2)
axis(1,x, cc)
for(k in 1:length(cc))
{
xx=seq(x[k]-dd/bf,x[k]+dd/bf,length=2)
points(xx,rep(l[k],length(xx)),type="l",lwd=lwd)
points(xx,rep(u[k],length(xx)),type="l",lwd=lwd)
points(c(x[k],x[k]),c(l[k],u[k]),type="l",lwd=lwd)
}

}



predict_htb=function(object,x=NULL,time=NULL,id=NULL,ntrees=NULL,type="response",se=FALSE)
{

	# object=ff;x=NULL;yindx=NULL;time=NULL;id=NULL;ntrees=NULL;type="response";se=FALSE

	yindx=NULL;
	if((!is.null(x))&(is.null(id)|is.null(time)))
	{
		## if x supplied but not time/id then data assumed iid at same time point (ie standard random forest)
		id=c(1:nrow(x))
		time=rep(1,nrow(x))
	}


	if(class(object)!="htree")
		stop("'object' not of class 'htree' ")

	if(object$rf==1)
		stop("'object' fit by 'hrf', use 'predict_hrf' ")


	
	ii=NULL 
	if(!is.null(x))
	{
	
		  if (any(is.na(x))) {
		    	stop("Missing data in 'x'.", call. = FALSE)
  			}
		  if (any(is.na(id))) {
 			   stop("Missing data in 'id'.", call. = FALSE)
			  }
		  if (any(is.na(time))) {
 			   stop("Missing data in 'time'.", call. = FALSE)
 			 }


		if(sum(unlist(lapply(x,class))!=object$vtype)>0)
			stop("Variable class mismatch with training data.")

		if(ncol(x)!=object$dimx[2])
			stop("Number of columns in 'x' differs from training data.")

		if(nrow(x)==0)
			stop("Zero of rows in 'x'.") 
		
	
		if(length(object$flist)>0)
		{

			# map strings/factors into integers 
			x=format_data(x=x,flist=object$flist)
		}
		if(is.null(id)|is.null(time))
			stop(" Arguments 'id' and 'time' cannot be empty.")
		if(is.character(id))
			id=as.factor(id)

		if(is.factor(id))
			id=as.numeric(id)

		if(!is.numeric(time))
			stop(" 'time' must be numeric.")

		ii=order(id,time)
		x=x[ii,,drop=FALSE]
		id=id[ii]
		time=time[ii]

		
	
	}else{
		id=time=NULL
	}


	if(is.null(x))
	{
		x=object$x
		yindx=object$yindx
		time=object$time
		id=object$id
	}


	if(0){
	if(!is.element(class(yindx),c("numeric","character")))
	{
		stop(" 'yindx' must give column number, or name, of response in 'x'.") 
	}else{
		if(class(yindx)=="character")
			yindx=which(names(x)==yindx)

		if(length(yindx)==0)
			stop(" 'yindx' not found ")

	}

	if(!is.null(yindx))
		yindx=yindx-1  # 0-offset 

	if(is.null(yindx))
		yindx=object$yindx	


	 if(yindx>ncol(x)|yindx<0)
		stop(" 'yindx' out of range")

	}
	
	pred=predict_htree(object=object,x=x,yindx=object$yindx,time=time,id=id,ntrees=ntrees,time_split=1)
	
	if(object$method_family==2) # logistic regression
	{
		if(type=="response") 
			pred=1/(1+exp(-pred))
	}


	if((!se)|(length(object$cv)==0))
	{
	}else{

		B=length(object$cv)
		pm=NULL
		for(b in 1:B)
		{
			hh=predict_htree(object=object$cv[[b]],x=x,yindx=yindx,time=time,id=id,ntrees=ntrees,time_split=object$time_split,all.trees=FALSE)
			pm=cbind(pm,hh)
		}	
		
		# delete-d jackknife (d=object$cv.fold)
		v=apply(pm,1,var)*((B-1)/B)*(object$cv.fold-1)
		se_pred=sqrt(v)
		#pred=list(pred=pred,se=se_pred)
		pred=cbind(pred,se_pred)
		colnames(pred)=c("pred","se")
		if(object$method_family==2&type=="response")
			cat("Prediction standard errors apply to link, not probabilitites.\n")
	}
	pred
}


histcon_aux=function(vh,vc,x,id)
{
	# return concurrent and historical column indexes 

	tiny_number=.0000001

	if(!is.null(vh))
	{

		if(is.null(vc))
	    		vc=c(1:ncol(x)) #[-vh]
		vc=vc-1
		vh=vh-1

	}else{
	    # identify time-varying predictors 
	    xaux=aggregate(x,by=list(id=id),sd)[,-1]
	    xaux[is.na(xaux)]=0

	    sd_x_id=apply(xaux,2,sd)
	    vh=c(1:ncol(x))[sd_x_id>tiny_number]
	    if(length(vh)==0)
		{
			#cat("No across time variation within levels of 'id'. \n")
			time_split=0
			vh=NULL
			vc=NULL
		}else{
			if(is.null(vc))
				vc=c(1:ncol(x)) # [-vh]
	    		vc=vc-1
	   		vh=vh-1
		}
	
	}

	list(vh=vh,vc=vc)

}

