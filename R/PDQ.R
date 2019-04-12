PDQ<-function(data=NULL,K=2,method="random", distance="euc", cent=NULL)
{
  #   %%%%%%%%%%%%%%Performs PDQ-Clustering%%%%%%%%%%%%%%%%%%
  #     %%%%%%%%INPUT%%%%%%%%%
  #   %data=data matrix  nxp
  #   %K= number of clusters
  #   %method=method to initialize centers
  #   %distance=dissimilarity/similarity measure
  #   %cent=user inputed center starts
  #   %%%%%%%OUTPUT%%%%%%%%%%
  #   %center=cluster's center  
  #   %class=class label
  #   %prob.m=nxk matrix-probability to belong to each class
  #   %JDF  join distance function
  #   %count=number of iterations until convergence

  # library(klaR)
 
  #   %input parameter check 
  if((!is.double(K))&(!is.integer(K))){stop("The number of clusters (K) must take an integer value.")}
  if(K<2){stop("The number of clusters (K) must be greater than one.");}
  if((K-round(K)!=0)){stop("The number of clusters (K) must be a whole number.");}
  if(!is.double(data)){stop("All elements of data must have type double.");}
  if(distance!="euc"&&distance!="chi"){stop("Distance measures available are euc and chi")}
  if(distance=="chi" && all(data>=0)!=TRUE){stop("When using chi distance all values must be positive.")}
  if(distance=="chi" && all(apply(data,1,sum)!=0)==FALSE){stop("When using chi distance no row entry can have all zero entries")}
  if(method!="random"&&method!="kmedoid"&&method!="center"){stop("Method options available are random, kmedoid and center")}
  if(method=="center" && is.null(cent)){stop("Error: If using method center must provide the starting centers")}
  if(method=="center"&&(!is.matrix(cent))){stop("Error: Centers must be in matrix form")}
  if(method=="center"&& K!=nrow(cent)){stop("Number of Centers must be the same as the number of clusters")}
  if(method=="center" && ncol(data)!=ncol(cent)){stop("The centers must have the same number of parameters as the dataset")}
  
  #  %initialize objects 
  data<-as.matrix(data)
  N=nrow(data)
  J=ncol(data)
  jdf<-vector()
  count<-1 
  epsilon=.001
  iter<-1
  dist=distance
 
  
  


  
 #    %initialize centers based on method 
  
 #    %method random, random 10 starts and choose best center based on JDF
  if(method=="random"){
    JDFini=vector()
   temp.center<-list()
  for(t in 1:10){
  clus.size<-c(rep(floor(N/K),each=K-1),(N-((K-1)*floor(N/K)))) 
  x<-vector()
  for(i in 1:J)
  {
    x<-c(x,runif(K,min(data[,i]), max(data[,i])))
  }
  center<-matrix(x,K,J)
  temp.center[[t]]<-center
  update=corePDQ(data,K,center,N,J,iter=10,clus.size,dist)
  JDFini[t]=update$jdfFinal}
  center= temp.center[[which(JDFini==min(JDFini))]]}

  #    %method user inputed center or kmedoid center
  if(method=="center"){center=cent}
  if(method=="kmedoid"){center=pam(data,K)$medoids}
  
  #    %re-initialize cluster sizes to restart algorithm 
  clus.size<-c(rep(floor(N/K),each=K-1),(N-((K-1)*floor(N/K))))
  
  #    %run main function  
  update=corePDQ(data,K,center,N,J,iter=100,clus.size,distance = dist)
  #check classification
  class<-apply(update$prob.m,1,which.max)
  #output 
  out<-list(label=class,centers=update$center,probability=update$prob.m,JDF=update$jdfFinal, iter=update$count,jdfvector=update$jdfseq)
  return(out)
}


  #    %main pdq function 

corePDQ<- function(data=NULL,K,center,N,J,iter=100,clus.size,distance){
  #    %initialize code objects 
  
  jdf<-vector()
  cent.tot=.2
  count=1
  epsilon=.01
  #    %run code until covergence by epsilon or iteration max 
  while(cent.tot>epsilon && count<iter)
  {
    
  #    %calculate distances by euc or chi 
    dist.m<-matrix(0,N,K)
    
if(distance=="euc")
{
    for(l in 1:K)
    {
      for(i in 1:N)
      {
        dist.m[i,l]<-sqrt(colSums(as.matrix(as.numeric(((data[i,]-center[l,])^2)))))
      }
    }}

  #    %calculate chi distances 
    if(distance=="chi")
    {
      for(l in 1:K)
      {
        for(i in 1:N)
        {
          dist.m[i,l]<-chi2Dist(rbind(data[i,],center[l,]))$D[2,1]
        }
      }}
  
  #    %check distance matrix for non numeric entries
    if(sum(!is.finite(dist.m))!=0){dist.m<-finite_check(dist.m)}
    
  #    %calculate probabilities 
    prob.num<-matrix(0,N,K)
    for(l in 1:K)
    {
      for(i in 1:N)
      {
        prob.num[i,l]<-prod(dist.m[i,-l])/prod(clus.size[-l])
      }
    }
    prob.den<-rowSums(prob.num)
    prob.m<-prob.num/prob.den
  #    %check matrix for non numeric entries
    if(sum(!is.finite(prob.m))!=0){prob.m<-finite_check(prob.m)}
    
  #    %update cluster sizes 
    clus.size<-vector()
    for(l in 1:K)
    {
      clus.size[l]<-(N*sqrt(sum(dist.m[,l]*(prob.m[,l]^2))))/sum(sqrt(colSums((dist.m*(prob.m^2)))))
    }
    
   #   %update centers 
    
    
   #   %calculate individual uK_i
   uk<-((prob.m)^2)/dist.m
   
   #    %check matrix for non numeric entries
   if(sum(!is.finite(uk))!=0){uk<-finite_check(uk)}
   
    coluk<-matrix(colSums(uk),N,K, byrow = T)
    uk_coluk<-uk/coluk
   #    %check matrix for non numeric entries  
    if(sum(!is.finite(uk_coluk))!=0){uk_coluk<-finite_check(uk_coluk)}
    
    cent_1<-matrix(0,K,J)
    for(l in 1:K)
    {
      cent_1[l,]<-colSums(data*uk_coluk[,l])
    }
    
  #    %check matrix for non numeric entries
    if(sum(!is.finite(cent_1))!=0){cent_1<-finite_check(cent_1)}
    
    
  #    % calculate joint distance fuction 
    dist.num<-apply(dist.m,1,prod)/prod(clus.size)
    dist.den<-matrix(0,N,K)
    for(l in 1:K)
    {
      for(i in 1:N)
      {
        dist.den[i,l]<-prod(dist.m[i,-l])/(prod(clus.size[-l]))
      }}
    
    jdf.temp<-sum(dist.num/rowSums(dist.den))
    jdf[count]<-jdf.temp
    count<-count+1
    
    
  #    %calculate difference in centers for convergence 
    
    cent.tot<-sum(rowSums(((center-cent_1)^2)))
    center<-cent_1
  }
return(list(center=center,prob.m=prob.m,jdfseq=jdf,jdfFinal= jdf[count-1], count=count))
}

#    %ancillary function to check matrix for non-numeric entries
finite_check<-function(mat)
{
  temp<-ifelse(!is.finite(mat),0,mat)
  
  for( i in 1:ncol(mat)){
    for( j in 1:nrow(mat)){
      if(!is.finite(mat[j,i]))
      {
        
        mat[j,i]=max(temp[1:(j-1),i])
      }
    }
  }
  return(mat)
}

