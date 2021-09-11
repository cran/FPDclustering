TPDC<- function(data=NULL,  k=2, method="kmedoids",nr=5,iter=100) {
  # Cluster the data whit pd-clustering Algoritmh
  #
  #
  #%%%%%%INPUT%%%%%%%%%
  #data=input data
  #k=number of cluster
  #
  #
  #%%%%%%OUTPUT%%%%%%%%
  #cnew=cluster's center  
  #l=class label
  #p=nxk matrix
  #probability to belong to each class 
  #JDF  join distance function
  #cont=number of iterations until convergence
  if((!is.double(k))&(!is.integer(k))){stop("The number of clusters (k) must take an integer value.")}
  if(k<2){stop("The number of clusters (k) must be greater than one.");}
  if((k-round(k)!=0)){stop("The number of clusters (k) must be a whole number.");}
  data=as.matrix(data)
  if(!is.double(data)){stop("All elements of data must have type double.");}
  n=nrow(data)
  J=ncol(data)
  temp.center<-list()
  JDFini=matrix(0,5,1)
  s=list()
  for(i in 1:k){
    s[[i]]=cov(data)#diag(J)
  }
  if(method=="random"){
    for(t in 1:nr){
      x<-vector()
      for(i in 1:J)
      {
        x<-c(x,runif(k,min(data[,i]), max(data[,i])))
      }
      center<-matrix(x,k,J)
      temp.center[[t]]<-center
      update=corePDT(data,k,center,s,n=nrow(data),J=ncol(data),iter=4)
      JDFini[t]=update$JDF}
    center= temp.center[[which(JDFini==min(JDFini))[1]]]}
  else{center=pam(data,k)$medoids}
  cnew=center
  update=corePDT(data,k,cnew,s,n=nrow(data),J=ncol(data),iter=iter)
  #check classification
  class<-apply(update$probability,1,which.max)
  #output 
  out<-list(label=class,centers=update$centers,sigma=update$sigma,df=update$df,probability=update$probability,JDF=update$JDF, iter=update$cont,data=data)
  class(out) <- "FPDclustering"
  return(out)
}

#km=kmeans(data, k)
#cnew=km$centers



corePDT<- function(data=NULL,k,cnew,s,n,J,iter=100){
  dfo<-df<-vector(mode="numeric",length=k)+20
  ver=100
  cont=2
  JDFv=matrix(0,iter,1)
  JDFv[1,1]=10000000000002
  JDFv[2,1]=10000000000000
  eps=.Machine$double.xmin
  rmax=.Machine$double.xmax
  ##Step 0
  dis=matrix(0,n,k)
  for( i in 1:k){
    sig=sweep(data,2,cnew[i,],FUN="-")
    sig=apply(sig^2,2,sum)
    sig=sqrt(sig/n)
    dataz=scale(data,cnew[i,],scale=sig)#/sqrt(diag(s[[i]]))#sqrt(diag(cov(data)))#,sqrt(diag(s[[i]])))
    #     meanden=dmvt(x=cnew[i,],delta=cnew[i,],sigma=cov2cor(s[[i]]),df=df[i],log=FALSE,type="shifted")
    #     dis[,i]<-log(meanden/(dmvt(x=data,delta=cnew[i,],sigma=cov2cor(s[[i]]),df=df[i],log=FALSE,type="shifted")))
    #      meanden=dmvt(x=rep(0,J),sigma=cov2cor(s[[i]]),df=df[i],log=TRUE)
    #     dis[,i]<-meanden-(dmvt(x=dataz,sigma=cov2cor(s[[i]]),df=df[i],log=TRUE))
    meanden=dmvt(x=rep(0,J),sigma=cov2cor(s[[i]]),df=df[i],log=FALSE)
    den=(dmvt(x=dataz,sigma=cov2cor(s[[i]]),df=df[i],log=FALSE))
    mm=min(den[which(den>0)])
    den[which(den==0)]=mm
    dis[,i]<-log(meanden/den)
    
    
  }
  
  dis[dis>rmax^(1/(k))]=rmax^(1/(k))
  dis[is.nan(dis)]=eps
  dis[is.na(dis)]=eps
  
  #computation of centers and probabilities
  t=matrix(1,n,k)
  p=matrix(0,n,k)
  for( i in 1:k){
    t2=as.matrix(dis[,-i])
    t[,i]=apply(t2,1,prod)}
  tot=apply(t,1,sum)
  p=t/tot
  while(ver>0.0000001 && cont<iter && JDFv[cont-1,1]!=JDFv[cont,1]){
    
    cont=cont+1
    c=cnew
    #STEP 1
    #computation of distance matrix
    dis=matrix(0,n,k)
    for( i in 1:k){
      sig=sweep(data,2,cnew[i,],FUN="-")
      sig=apply(sig^2,2,sum)
      sig=sqrt(sig/n)
      dataz=scale(data,cnew[i,],scale=sig)#/sqrt(diag(cov(data)))#,sqrt(diag(s[[i]])))
      
      #     meanden=dmvt(x=cnew[i,],delta=cnew[i,],sigma=cov2cor(s[[i]]),df=df[i],log=FALSE,type="shifted")
      #     dis[,i]<-log(meanden/(dmvt(x=data,delta=cnew[i,],sigma=cov2cor(s[[i]]),df=df[i],log=FALSE,type="shifted")))
      #      meanden=dmvt(x=rep(0,J),sigma=cov2cor(s[[i]]),df=df[i],log=TRUE)
      #     dis[,i]<-meanden-(dmvt(x=dataz,sigma=cov2cor(s[[i]]),df=df[i],log=TRUE))
      meanden=dmvt(x=rep(0,J),sigma=cov2cor(s[[i]]),df=df[i],log=FALSE)
      den=(dmvt(x=dataz,sigma=cov2cor(s[[i]]),df=df[i],log=FALSE))
      mm=min(den[which(den>0)])
      den[which(den==0)]=mm
      dis[,i]<-log(meanden/den)
      
     
    }
  
    ##Cluster size
    clus.size<-vector()
    for(l in 1:k)
    {
      clus.size[l]<-(n*sqrt(sum(dis[,l]*(p[,l]^2))))/sum(sqrt(colSums((dis*(p^2)))))
    }
    
    dis2=sweep(dis,2,clus.size,FUN="/")
    
    #STEP 2
    #computation of centers and probabilities
    t=matrix(1,n,k)
    p=matrix(0,n,k)
    for( i in 1:k){
      t2=as.matrix(dis2[,-i])
      t[,i]=apply(t2,1,prod)}
    tot=apply(t,1,sum)
    p=t/tot
   
    
    for(i in 1:k){
      #  mah<-mahalanobis(data, cnew[i,], s[[i]],inverted=FALSE)	
      #centers
      
      # b<-apply((data*p[,i]^2)/(df[i]+mah),2,sum)
      #cnew[i,]<-b/sum((p[,i]^2)/(df[i]+mah))			
      #sigma
      ##############
      #############
      mah<-mahalanobis(data, cnew[i,], s[[i]],inverted=FALSE)	
      diff<-sweep(data, 2, cnew[i,], FUN="-")
      
      w<-((df[i]+J)*p[,i]^2/(mah+df[i]))#/(sum(p[,i]^2))
      wm=w/sum(w)
      ws=w/sum(p[,i]^2)
      cnew[i,]=apply((data*wm),2,sum)
      #s[[i]]<-(t(diff)%*%(diff*ws))
      
      
      s[[i]]= cov.wt(data, center=cnew[i,],wt = ws, method="ML")$cov
      
      
      mah1<-mah2<-0
      mah<-mahalanobis(data, cnew[i,], s[[i]],inverted=FALSE)	
      mah1=sum(p[,i]^2*log(1+mah/df[i]))
      mah2=sum(p[,i]^2*mah/(df[i]+mah))
      
      
      
      df[i]= mle.nu(p[,i],J,mah,mah1,mah2,df[i],k)#,silent=TRUE
      
      
      df=as.vector(as.numeric(df))
      
    }
    #check if centers move
    ver1=matrix(0,1,k)
    for( j in 1:k){
      ver1[j]=sqrt(sum((cnew[j,]-c[j,])^2))
    }
    ver=sum(ver1)
    JDFv[cont,]= trunc(sum(dis*p),4)
    #  JDFv[cont,]= trunc(sum(p^2*(log(meanden)-log(dataden))),4)
  }
  if( cont==1000){
    print('Convergence not reached')}
  JDFv=JDFv[3:cont,]
  #memebership definition according to the maximum probability
  l=max.col(p)
  # computation of JDF
  JDF=sum(mean(dis*p))
  out=list(label=l, centers=cnew,sigma=s, df=df,probability=p, JDF=JDF,iter=cont,JDFv=JDFv)
  out
}



mle.nu <- function(p,J,mah,mah1,mah2,df,k){
  #p is the probability matrix, k is the cluster, J is the # of params, mah is the mahalanobis
  df=as.numeric(df)
  p2=sum(p^2)
  #	val <- (uniroot.all(function( z ) ((-digamma((z+J)/2) + digamma(z/2) +.5*J/z)*p2+ (.5*mah1-.5*((z+J)/z)*mah2)),lower=2,upper=200,maxiter=30))
  val <- (uniroot.all(function( z ) (1+(digamma((z+J)/2) - digamma(z/2) )*2+log(z)- sum(p^2*(log(df+mah)))/p2-(z+J)*sum(p^2*(df+mah)^(-1))/p2),lower=2,upper=200,maxiter=10))
  
  
  if(length(val)!=1)val=val[1]
  #if(!is.numeric(val))val=df
  #	if(val > 200) val <- 200
  #	if(val < 2) val <- 2
  return(val)
}

