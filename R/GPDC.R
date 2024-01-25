GPDC<- function(data=NULL,  k=2, ini="kmedoids",nr=5,iter=100) {
  # Cluster the data whit pd-Gaussian Algoritmh
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
#   if((!is.double(k))&(!is.integer(k))){stop("The number of clusters (k) must take an integer value.")}
#   if(k<2){stop("The number of clusters (k) must be greater than one.");}
#   if((k-round(k)!=0)){stop("The number of clusters (k) must be a whole number.");}
#   data=as.matrix(data)
#   if(!is.double(data)){stop("All elements of data must have type double.");}
  method=ini
  n=nrow(data)
  J=ncol(data)
 temp.center<-list()
 temp.l<-list()
 JDFini=matrix(0,nr,1)
 if(method=="random"){
   for(t in 1:nr){
     x<-matrix(0,k,J)
     for(i in 1:J)
     {
       x[,i]<-runif(k,min(data[,i]), max(data[,i]))
     }
     center<-x
     temp.center[[t]]<-center
     
     update=corePDGaus(data,k,center,n=n,J=J,iter=10)
     temp.l[[t]]=apply(update$probability,1,which.max)
     JDFini[t]=update$JDF[length(update$JDF)]}
   center= temp.center[[which(JDFini==min(JDFini))[1]]]
   l=temp.l[[which(JDFini==min(JDFini))[1]]]
   }
 else if(method=="PDclust"){
   ini=PDC(data,k)
   center=ini$centers
   l=ini$label
 }
 else{kmed=pam(data,k)
   center=kmed$medoids
 l=kmed$clustering}
 cnew=center
 update=corePDGaus(data,k,cnew,n=n,J=J,iter=iter,lab=l)
 #check classification
 class<-apply(update$probability,1,which.max)
 #output 
 out<-list(label=class,centers=update$centers,sigma=update$sigma,probability=update$probability,JDF=update$JDF, iter=update$iter,data=data)
 class(out) <- "FPDclustering"
 return(out)
}
  #km=kmeans(data, k)
  #cnew=km$centers
  
 
 corePDGaus<-function(data=NULL,k,cnew,n,J,iter=100,lab=NULL){
  ver=100
  cont=2
  sigma=list()
 #   matrix(diag(J),J*k,J,1)
  JDFv=matrix(0,iter,1)
  dis=matrix(0,n,k)
  den=matrix(0,n,k)
  p=matrix(0,n,k)
  for(j in 1:k){
    if(is.null(lab)){s=cov(data)}
   else{ s=cov(data[which(lab==j),])}#diag(J)#sigma[(1+(J*(j-1))):(J*j),]
    sigma[[j]]=s
       den[,j]=log(dmvnorm(cnew[j,],cnew[j,],s)/dmvnorm(data,cnew[j,],s))
      # den[,j]=(1/dmvnorm(data,cnew[j,],s))
       
         }
  dis=(den)#-matrix((1/dmax),n,k)
  eps=.Machine$double.xmin
  rmax=.Machine$double.xmax
  dis[dis>rmax^(1/(k))]=rmax^(1/(k))
  dis[is.nan(dis)]=eps
  dis[is.na(dis)]=eps
  t=matrix(1,n,k)
  
  for( i in 1:k){
    t2=as.matrix((dis[,-i]))
    t[,i]=apply(t2,1,prod)}
  tot=apply(t,1,sum)
  p=sweep(t,1,tot,FUN="/")
  JDFv[1,1]= sum(dis*p)+2
  JDFv[2,1]= sum(dis*p)


  while(ver>0.0000001 && cont<iter){#abs(JDFv[cont,1]-JDFv[cont-1,1])>0.01 &
    cont=cont+1
    c=cnew
    
    for( j in 1:k){
      
      
      #Update mu and sigma
      s=sigma[[j]]
    # w= p[,j]/sum( p[,j])
     w= p[,j]^2/sum( p[,j]^2)
      cnew[j,]=apply((data*w),2,sum)
      #dif=sweep(data, 2, cnew[j,], FUN="-") 
      sigma[[j]]= cov.wt(data, center=cnew[j,],wt = w, method="ML")$cov
      #  sigma[[j]]=t(dif*w)%*%(dif)
        s= sigma[[j]]
        
        ##Updates distances
        
        sig=sweep(data,2,cnew[j,],FUN="-")
        sig=apply(sig^2,2,sum)
        sig=sqrt(sig/n)
        dataz=scale(data,cnew[j,],scale=sig)#/sqrt(diag(cov(data)))#,sqrt(diag(s[[i]])))
        
        #dif=sweep(data, 2, cnew[j,], FUN="-")#/sqrt(diag(cov(data)))
   #   den[,j]=log(dmvnorm(cnew[j,],cnew[j,],s)/dmvnorm(data,cnew[j,],s))
       den[,j]=log(dmvnorm(rep(0,J),sigma=cov2cor(s))/dmvnorm(dataz,sigma=cov2cor(s)))
    #    den[,j]=log(dmvnorm(rep(0,J),sigma=(s))/dmvnorm(dif,sigma=(s)))
       den[which(den[,j]==Inf),j]=rmax/2

       
      }
      
    dis=den
    dis[dis>rmax^(1/(k))]=rmax^(1/(k))
    dis[is.nan(dis)]=eps
    dis[is.na(dis)]=eps
    
    ##Cluster size
    clus.size<-vector()
    pdsum=sum(sqrt(colSums((dis*(p^2)))))
    for(l in 1:k){
      clus.size[l]<-(n*sqrt(sum(dis[,l]*(p[,l]^2))))/pdsum
    # clus.size[l]<-(n*sqrt(sum(dis[,l]*(p[,l]))))/sum(sqrt(colSums((dis*(p)))))
    }
    
    ##check for empty clusters
    if(sum(clus.size<3)>=1){
      clus.size=rep(n/k,k)
      x<-matrix(0,k,J)
      for(i in 1:J){
        x[,i]<-runif(k,min(data[,i]), max(data[,i]))
      }
      cnew<-x
      
    }

  #  clus.size[k]=n-(sum(clus.size)-clus.size[k])
    dis2=sweep(dis,2,clus.size,FUN="/")
    for( i in 1:k){
      t2=as.matrix(dis2[,-i] )  
      t[,i]=apply(t2,1,prod)}
    tot=apply(t,1,sum)
    p=sweep(t,1,tot,FUN="/")
   # p[p<0]=0
  #  traces= (sapply(sigma,trace))
    el1=apply(dis*p,2,sum)
    JDFv[cont,]= trunc(sum(el1),4)#*traces
    
    #check if centers move
    ver1=matrix(0,1,k)
    for( j in 1:k){
      ver1[j]=sqrt(sum((cnew[j,]-c[j,])^2))
    }
    ver=sum(ver1)
    
 # par(new=TRUE)
  #    points(cnew,col=cont)
    
  }
 # if( cont==1000){
 #   print('Convergence not reached')}
  #JDFv=sum(log(dis)*p^2)#
  JDFv=JDFv[3:cont,]
  #memebership definition according to the maximum probability

  # computation of JDF
  out=list(centers=cnew,sigma=sigma, probability=p,iter=cont,JDF=JDFv,data=data)

  out
 }

 