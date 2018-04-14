mcode=function(m)
{
	# translate m string to method integer code 
	clist=list()
	clist[["freq"]]=1
	clist[["frac"]]=4
	clist[["mean"]]=3
	clist[["sum"]]=5
	clist[["mean0"]]=6
	clist[["sum0"]]=7
	clist[["meanw0"]]=8
	clist[["sumw0"]]=9
	clist[["freqw"]]=10
	clist[["fracw"]]=11
	clist[["class"]]=12


	mc=clist[[as.character(m)]]
	if(is.null(mc))
	{
		cat(paste(" method: '",m,"' not matched, setting to 'freq'. \n",sep=""))
		mc=1
	}

	return(mc)
}

# -------------------------------- #
# hrf: historical random forest
# -------------------------------- #

hrf=function(x,time=NULL,id=NULL,yindx,ntrees=100,method="freqw",mtry=NULL,se=FALSE,B=100,R=10,nsplit=NULL,nsamp=5,
		historical=TRUE,keep_data=TRUE,vh=NULL,vc=NULL,delta=NULL,control=list(nmax=10,nodesize=1))
{

	rsplit=FALSE;

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

	if(!is.data.frame(x))
		x=as.data.frame(x)

	if(is.factor(x[,yindx]))
	{
		if(length(levels(x[,yindx]))>50)
			stop(" Number of classes too large (limit is 50) ")
	
		if(method!="freqw")
			cat("Setting 'method=freqw', only method implemented for classification.\n")
		method="class"

		if(se==TRUE)
		{
			se=FALSE
			cat("Standard error estimation not implemented for classification.\n")
		}
	}

	## -- convert factors/strings to numeric 
	rd=reformat_data(x=x)
	x=rd$x

	ii=order(id,time)
	x=x[ii,]
	id=id[ii]
	time=time[ii]

	flist=rd$flist


	# method code
	method=mcode(m=method)


	tiny_number=.00001
	lambda=1
	id_sampling=TRUE
	rf=1

	if(is.null(nsplit))
		nsplit=nrow(x)

	tf=1
	cv_train=NULL
	##method=1
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
			#cat(" No across time variation within levels of 'id'. \n")
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



	fit_rf=htree(x=x,time=time,id=id,yindx=(yindx-1),ntrees=ntrees,lambda=lambda,rf=rf,nsplit=nsplit,nsamp=nsamp,tf=tf,id_sampling=id_sampling,rsplit=rsplit,
		mtry=mtry,vi=vi,time_split=time_split,oob=oob,oobmatrix=oobmatrix,keep_data=keep_data,vh=vh,vc=vc,cv_train=cv_train,method=method,delta=delta,control=control)

	
	fit_rf$classInfo=flist[[names(x)[yindx]]]
	fit_rf$vtype=vtype
	fit_rf$indx_original=ii




	fit_boot=NULL
	if(se)
	{
		# Estimate standard errors (noisy bootstrap) 
		fit_boot=seboot(B=B,M=R,object=fit_rf)
	
	}

	#if(keep_data)
	#	fit_rf$x=deformat_data(fit_rf$x,flist)

	fit_rf$flist=flist


	fit_rf$boot=fit_boot
	fit_rf
}




factor_conversion=function(y)
{
	# convert factor-levels to integer: 0:(length(levels(y))-1)
	# 


	flist=list()
	flist[levels(y)]=c(0:(length(levels(y))-1))
	ilist=list()
	ilist[as.character(c(0:(length(levels(y))-1)))]=levels(y) 

	yint=as.numeric(unlist(flist[y]))
	
	classInfo=list()
	classInfo=list(levels=levels(y),flist=flist,ilist=ilist)
	rlist=list(classInfo=classInfo,y=yint)

	rlist

}


reformat_data=function(x)
{
	# Convert factors/character-variables to numeric 

	fvar=NULL
	flist=list()
	nx=names(x)
	for(k in nx)
	{
		if(is.character(x[[k]]))
			x[[k]]=as.factor(x[[k]])

		if(is.factor(x[[k]]))
		{
			h=factor_conversion(x[[k]])
			x[[k]]=h$y
			flist[[k]]=h$classInfo
		}
		
	} 

	list(x=x,flist=flist)
}


deformat_data=function(x,flist)
{

	if(!is.null(flist))
	{

		for(k in names(flist))
			x[[k]]=as.factor(as.character(unlist(flist[[k]]$ilist[as.character(x[[k]])])))
		
	}
	

	x
}


format_data=function(x,flist)
{
	# Apply data-formatting used in model fit 
	# .. convert strings/factors to integers 

	if(!is.null(flist))
	{

		for(k in names(flist))
		{
			uc=unique(as.character(x[[k]]))
			if(sum(!is.element(uc,flist[[k]]$levels))>0)
				stop(paste(" Factor levels in ",k," not in training data, exiting...",sep="")) 
			
			x[[k]]=(as.numeric(unlist(flist[[k]]$flist[as.character(x[[k]])])))
			
		}


	}


	x
}






predict_hrf=function(object,x=NULL,time=NULL,id=NULL,all.trees=FALSE,se=FALSE)
{

	yindx=NULL

	if((!is.null(x))&(is.null(id)|is.null(time)))
	{
		## if x supplied but not time/id then data assumed iid at same time point (ie standard random forest)
		id=c(1:nrow(x))
		time=rep(1,nrow(x))
	}


	if(class(object)!="htree")
		stop(" 'object' not of class 'htree' ")

	if(object$rf!=1)
		stop(" 'object' not fit by 'hrf' ")


	if(!is.null(yindx))
		yindx=yindx-1

	if(is.null(yindx))
		yindx=object$yindx

	
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

		  if(yindx>ncol(x)|yindx<0)
			stop(" 'yindx' out of range")
		
	
	}else{
		id=time=NULL
	}

		

	# predict_se=function(m,object,x=NULL,id=NULL,time=NULL)
	pred=predict_htree(object=object,x=x,yindx=yindx,time=time,id=id,ntrees=NULL,time_split=1,all.trees=all.trees)
	retlist=pred
	if(!is.null(ii))
		retlist[ii]=pred

	if(is.null(x))  # -- put predictions on training data in original order 
	{

		if(object$ncat<0)
			retlist[object$indx_original]=pred
		if(object$ncat>0){
			
			retlist[object$indx_original,]=pred
			colnames(retlist)=object$classInfo$levels
		}
	}
	if(se)
	{
		if(!is.null(object$boot))
		{
			serr=predict_se(m=object$boot,object=object,x=x,id=id,time=time)
			retlist=cbind(pred,serr)
			if(!is.null(ii))
				retlist[ii,]=retlist
			if(is.null(x))  # -- put predictions on training data in original order 
				retlist[object$indx_original,]=retlist

			colnames(retlist)=c("pred","se")

		}else{
			cat(" Set argument 'se=TRUE' in 'hrf' to get standard errors.\n")
		}
	}
	
	retlist
}



partdep_hrf=function(object,xindx,xlim=NULL,ngrid=25,subsample=.1,plot.it=TRUE,cat.plot=FALSE,which.class=1)
{
	main="Partial dependence"

	if(object$ncat>0)
	{
		if(class(which.class)=="character")
			which.class=which(object$classInfo$levels==which.class)

		if(length(which.class)==0)
			stop(" 'which.class' not found. ")

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

	if(is.null(object$boot))
	{
		# No bootstrap samples, no standard error bars produced
		pd=partdep(object=object,xindx=xindx,xlim=xlim,ngrid=ngrid,ntrees=NULL,subsample=subsample,which.class=which.class)

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
				plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main=main,ylim=ylim,cex.axis=1.5,cex.lab=2,cex.main=2)
				points(pd$x,upper,type="l",lwd=3,lty=2)
				points(pd$x,lower,type="l",lwd=3,lty=2)
			}else{

				plot_cat2(m=pd$y,l=(pd$y-2*pd$se),u=(pd$y+2*pd$se),cc=pd$x,ylim=NULL,mult=1,bf=4,lwd=3,main=main)


			}
	
		}else{
			plot(pd$x,pd$y,type="l",lwd=3,xlab=nx,ylab="",main="Partial dependence",cex.axis=1.5,cex.lab=2,cex.main=2)
		}
	}
pd
}





varimp_hrf=function(object,nperm=20)
{


	if(class(object)!="htree")
		stop(" 'object' not of class 'htree'.")

	if(object$rf!=1)
		stop(" 'object' not fit by 'hrf'.")

	if(nperm<1)
		stop(" 'nperm<1'.")  

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
		as.integer(yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),
		as.integer(ncol(object$trees)),as.integer(ntrees_in_matrix),as.integer(time_split),as.integer(object$method))
	
	# 
	set_oob(object$oobmatrix)


	## -- dont need to do this .... 
	h=.C("predict_trees_all",as.integer(ntrees),pred=double(nrow(x)*ntrees))	
	predictions=matrix(h$pred,ncol=ntrees,nrow=nrow(x),byrow=FALSE)
	noob=apply(object$oobmatrix==1,1,sum)
	pred_oob=apply(predictions*object$oobmatrix,1,sum)
	pred_oob=pred_oob/noob

 	#void varimp(int *nperm,double *res)
	h=.C("varimp",as.integer(nperm),res=double(n*p))
	perm_pred=matrix(h$res,ncol=p,nrow=n,byrow=FALSE)

	# p=5;id=hb$id;pred_oob=hz$pred_oob;perm_pred=hz$pred_perm
	zscore=NULL
	pchange=NULL
	sdchange=NULL
	for(k in 1:p)
	{
		z=abs(y-perm_pred[,k])-abs(y-pred_oob)
		if(group=="id")
		{	# compute error change by subject/id 
			zid=aggregate(z,by=list(id=id),mean) 
			se_zid=sqrt(var(zid[,2])/length(unique(id)))
			pchange=c(pchange,mean(zid[,2]))
			sdchange=c(sdchange,se_zid)
			zscore=c(zscore,mean(zid[,2])/se_zid)
		}else{
			# error change by row, assuming independence between 
			se_z=sqrt(var(z)/n)
			zscore=c(zscore,mean(z)/se_z)
			pchange=c(pchange,mean(z))
			sdchange=c(sdchange,se_z)

		}
	}

	names(zscore)=colnames(x)
	err=min(object$error)
	dz=data.frame(delta_rel=round(pchange/err,3),delta_error=round(pchange,3),stderr=round(sdchange,3),zscore=round(zscore,3),pvalue=pnorm(zscore,lower.tail=FALSE))
	names(dz)=c("Relative change","Mean change","SE","Z-value","P-value") 
	#rownames(dz)=NULL
	# -- free pred info  
	h=.C("free_predict")

	if(object$time_split==0)
		dz=dz[-which(rownames(dz)==names(object$x)[object$yindx+1]),,drop=FALSE]

	#list(zscore=zscore,pred_oob=pred_oob,pred_perm=perm_pred)
	dz[order(dz[,4],decreasing=TRUE),]
}






















quantile_hrf<-function(object,x=NULL,time=NULL,id=NULL,ntrees=NULL,prob=seq(.1,.9,length=10))
{


if(class(object)!="htree")
	stop(" 'object' must be of class 'htree'.")

if(object$rf!=1)
	stop(" 'object' is not a random forest.")

if(object$ncat>0)
	stop(" 'quantile_hrf' only for regression, not classification.")


if((!is.null(x))&(is.null(time)|is.null(id)))
{
	## rows assumed independent
	id=c(1:nrow(x))
	time=rep(1,nrow(x))

}




time_split=1

if(is.null(ntrees))
	ntrees=length(unique(object$trees[,1]))

ntrees_in_matrix=length(unique(object$trees[,1]))


n=nrow(object$x)
p=ncol(object$x)

h=.C("read_predict",as.double(as.matrix(object$x)),as.integer(n),as.integer(p),as.double(object$time),as.integer(object$id),
		as.integer(object$yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),
		as.integer(ncol(object$trees)),as.integer(ntrees_in_matrix),as.integer(time_split),as.integer(object$method))

h=.C("set_tnode_row")



# -- free pred info  
h=.C("free_predict")



if(!is.null(x))
{

	hc=check_data(object=object,x=x,id=id,time=time)

	x=hc$x
	yindx=hc$yindx
	time=hc$time
	id=hc$id
	indx_original=hc$ii

}else{
	
	x=object$x
	yindx=object$yindx
	time=object$time
	id=object$id
	indx_original=object$indx_original

}

n=nrow(x)
p=ncol(x)

h=.C("read_predict",as.double(as.matrix(x)),as.integer(n),as.integer(p),as.double(time),as.integer(id),
		as.integer(yindx),as.double(t(object$trees)),as.integer(nrow(object$trees)),
		as.integer(ncol(object$trees)),as.integer(ntrees_in_matrix),as.integer(time_split),as.integer(object$method))


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

## -- put quantiles in order of original data 
qu=quants

qu[indx_original,]=quants

qu
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




check_data=function(object,x,id,time)
{


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

	list(x=x,id=id,time=time,ii=ii)

}

