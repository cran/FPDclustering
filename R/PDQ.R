PDQ<-function(x=NULL,k=2,ini="kmd",dist="euc",cent=NULL,ord=NULL,cat=NULL,bin=NULL,cont=NULL,w=NULL)

#%%%%%%INPUT%%%%%%%%%
  
#data=input data
#K=number of cluster
#method=method to initialize centers
#distance=dissimilarity/similarity measure
#cent=user inputed center starts


#%%%%%%OUTPUT%%%%%%%%
#center=cluster's center
#class=class label
#prob.m=nxk matrix-probability to belong to each class
#JDF  join distance function
#count=number of iterations until convergence
#jdfvector  join distance function per iteration

{K=k
    data=x
    method=ini
  distance=dist
  ord.loc=ord
  cat.loc=cat
  bin.loc=bin
  cont.loc=cont   
  weights=w  
  #Check Input Parameters
    if((!is.double(K))&(!is.integer(K))){stop("The number of clusters (K) must take an integer value.")}
    if(K<2){stop("The number of clusters (K) must be greater than one.");}
    if((K-round(K)!=0)){stop("The number of clusters (K) must be a whole number.");}
   # if(!is.double(data)){stop("All elements of data must have type double.");}
    if(distance!="euc"&&distance!="chi"&&distance!="gower"){stop("Distance measures available are euc and chi and gower")}
    if(distance=="chi" && all(data>=0)!=TRUE){stop("When using chi distance all values must be positive.")}
    if(distance=="chi" && all(apply(data,1,sum)!=0)==FALSE){stop("When using chi distance no row entry can have all zero entries")}
    if(method=="center" && is.null(cent)){stop("Error: If using method center must provide the starting centers")}
    if(method=="center"&&(!is.matrix(cent))){stop("Error: Centers must be in matrix form")}
    if(method=="center"&& K!=nrow(cent)){stop("Number of Centers must be the same as the number of clusters")}
    if(method=="center" && ncol(data)!=ncol(cent)){stop("The centers must have the same number of parameters as the dataset")}
    if(distance!="gower"&&(!is.null(ord.loc)||(!is.null(cat.loc))||(!is.null(bin.loc))||(!is.null(cont.loc))))
    {stop("Are you sure you did not intend to use gower's distance. You provided labels for your variables. If you inted to use any other distance metric other then gowers please remove the label identifiers from the function parameters.")}
  N=nrow(data)
  J=ncol(data) 
   
  if(distance=="gower"){
    dataT=data
    
    if(!is.null(cont.loc)){
    dataN=as.matrix(data[,cont.loc])
    dataO=data[,-cont.loc]
    pcon=ncol(dataN)
    dataN=matrix(as.numeric(unlist(dataN)),N,pcon)}
    data=data.frame(dataT)
    if(!is.null(cont.loc)){data[,cont.loc]=dataN}
    
    dist=daisy(data,metric="gower")
    center=pam(data,K)$medoids
  }else{
    #Initialize Objects
    data<-as.matrix(data)
   
    jdf<-vector()
    count<-1
    epsilon=.001
    iter<-1
    dist=distance
 
    #Initialize centers based on method
  
    # method random, random 10 starts and choose best center based on JDF
    if(method=="random"){ #&& distance=="euc"){
        # JDFini=list()
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
            update=corePDQ(data,K,center,N,J,iter=10,clus.size,dist,ord.loc,cat.loc,bin.loc,cont.loc,weights)
            # JDFini[[t]]=update$jdf}
            
            JDFini[t]=update$jdfFinal}
        center= temp.center[[which(JDFini==min(JDFini))]]}

    #method user inputed center or kmedoid center OR kmode center
    if(method=="center"){center=cent}
    if(method=="kmd"){center=pam(data,K)$medoids}
    if(method=="kmode"){center=kmodes(data,K)$modes}
    
  }
    #re-initialize cluster sizes to restart algorithm
    clus.size<-c(rep(floor(N/K),each=K-1),(N-((K-1)*floor(N/K))))
    
    #run main function
    update=corePDQ(data,K,center,N,J,iter=100,clus.size, distance, ord.loc,cat.loc,bin.loc,cont.loc,weights)
        
    #check classification
    class<-apply(update$prob.m,1,which.max)
    
    #output
    out<-list(label=class,centers=update$center, clust_size=update$clust_size,probability=update$prob.m,JDF=update$jdfFinal, iter=update$count,jdfvector=update$jdfseq,data=data,cont.loc=cont.loc)
    class(out) <- "FPDclustering"
    return(out)
}



#Main PDQ Function
corePDQ<- function(data=NULL,K,center,N,J,iter=100,clus.size,distance,ord.loc=NULL,cat.loc=NULL,bin.loc=NULL,cont.loc=NULL,weights=NULL){
    
    #initialize code objects
    jdf<-vector()
    cent.tot=.2
    count=1
    epsilon=.01
    
    #run code until covergence by epsilon or iteration max
    while(cent.tot>epsilon && count<iter)
    {
        
        #calculate distances
        
        dist.m<-matrix(0,N,K)
        
        #calculate gower distance
        if(distance=="gower")
        {
            dist.m<- gowers_dist(data, center=center,K=K, ord.loc,cat.loc,bin.loc,cont.loc,weights)
        }
        
        #calculate eucledian distance
        if(distance=="euc")
        {
            for(l in 1:K)
            {
                for(i in 1:N)
                {
                    dist.m[i,l]<-sqrt(colSums(as.matrix(as.numeric(((data[i,]-center[l,])^2)))))
                }
            }}
        # fr=apply(data,1,sum)
        # fc=apply(data,2,sum)
        
        
        
        #calculate chi distance
        if(distance=="chi")
        {
            for(l in 1:K)
            {
                for(i in 1:N)
                {
                    dist.m[i,l]<-chi2Dist(rbind(data[i,],center[l,]))$D[2,1]
                    # dist.m[i,l]<-sqrt(sum((data[i,]-center[l,])^2/center[l,]))
                }
            }}
        
        
        
        #check distance matrix for non numeric entries
        if(sum(!is.finite(dist.m))!=0){dist.m<-finite_check(dist.m)}
        
        
        #calculate probabilities
        prob.num<-matrix(0,N,K)
        for(l in 1:K) {
            for(i in 1:N) {
                prob.num[i,l]<-prod(dist.m[i,-l])/prod(clus.size[-l])
            }
        }
        prob.den<-rowSums(prob.num)
        prob.m<-prob.num/prob.den
        
        #check matrix for non numeric entries
        if(sum(!is.finite(prob.m))!=0){prob.m<-finite_check(prob.m)}
        
        
        #update cluster sizes
        clus.size<-vector()
        for(l in 1:K) {
            clus.size[l]<-(N*sqrt(sum(dist.m[,l]*(prob.m[,l]^2))))/sum(sqrt(colSums((dist.m*(prob.m^2)))))
        }
        
        #update centers
     
        
        #calculate individual uK
        uk<-((prob.m)^2)/dist.m
        
        #check matrix for non numeric entries
        if(sum(!is.finite(uk))!=0){uk<-finite_check(uk)}
        # temp<-ifelse(!is.finite(uk),0,uk)
        
      
        
        
        coluk<-matrix(colSums(uk),N,K, byrow = T)
        uk_coluk<-uk/coluk
        
        #check matrix for non numeric entries
        if(sum(!is.finite(uk_coluk))!=0){uk_coluk<-finite_check(uk_coluk)}
        
        cent_1<-matrix(0,K,J)
        if(distance!="gower"){
        for(l in 1:K) {
            cent_1[l,]<-colSums(data*uk_coluk[,l])
        }
        }
        #check for gowers distance to update center of binary variable to mode:
        if(distance=="gower") {
          if(!is.null(cont.loc)){
            if(length(cont.loc)>1){
          for(l in 1:K) {
            cent_1[l,cont.loc]<-colSums(data[,cont.loc]*uk_coluk[,l])
          }
            }else{ 
              for(l in 1:K) {
              cent_1[l,cont.loc]<-sum(data[,cont.loc]*uk_coluk[,l])}}}
            mode_vector<-length(K)
            #x<-which(lapply((apply(data,2,unique)),length)==2)
            #classify cluster
            class<-apply(prob.m,1,which.max)
            # unique_class<-unique(class)
            if(!is.null(ord.loc)){
              for(i in ord.loc)
              {
                for(l in 1:K)
                {
                  cent_1[l,i]<-sum(data[,i]*(prob.m[,l])^2)/sum((prob.m[,l])^2)
                }

                
              }}
            if(!is.null(bin.loc))
            for(i in bin.loc)
            {
                for(l in 1:K)
                {
                    mode_data<-getMode(data[class==l,i])
                    mode_vector[l]=mode_data
                }
                cent_1[,i]<-mode_vector
                
            }
            
            #update the categorical centers to mode
            if(!is.null(cat.loc))
            {
                mode_cat_vect<-length(K)
                for(m in cat.loc)
                {
                    for(j in 1:K)
                    {
                        mode_dat<-getMode(data[class==j,m])
                        mode_cat_vect[j]=mode_dat
                    }
                    cent_1[,m]<-mode_cat_vect
                }
            }
            
        }
        #check if errors/check matrix for non numeric entries
        if(sum(!is.finite(cent_1))!=0){cent_1<-finite_check(cent_1)}
        
        
        # calculate joint distance fuction
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
        
        
        #calculate difference in centers for convergence
        cent.tot<-sum(rowSums(((center-cent_1)^2)))
        center<-cent_1
    }
    return(list(center=center,clust_size=clus.size,prob.m=prob.m,jdfseq=jdf,jdfFinal= jdf[count-1], count=count))
}



#ancillary function to check matrix for non-numeric entries
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


#gowers distance function

gowers_dist<-function(data,center, K,ord.loc=NULL,cat.loc=NULL,bin.loc=NULL,cont.loc=NULL, weights=NULL){
    
    if(is.null(ord.loc)&(is.null(cat.loc))&(is.null(bin.loc))&(is.null(cont.loc)))
    {
        stop("When Using Gower's Distance you must specify the labels for all the variables")
    }
    
    label<-sort(c(ord.loc,cat.loc,bin.loc,cont.loc))
    
    if(sum(!(label==(1:ncol(data))))!=0){stop("You have Not specified all the labels for the variables")}
    
    # #check data type
    # if(method=="nc")
    # {
    # data.check<-lapply(apply(data,2,unique),length)
    # bin.loc<-as.vector(which(unlist(data.check)==2))
    # cont.loc<-as.vector(which(unlist(data.check)>15))
    #
    # if(!is.null(ord))
    # {
    # cat.loc<-(1:ncol(data))[-c(ord,bin.loc,cont.loc)]
    # }
    # cat.loc<-(1:ncol(data))[-c(bin.loc,cont.loc)]
    # }
    #
    #
    #
    #
    # if(method=="cat")
    # {
    #   data.check<-lapply(apply(data,2,unique),length)
    #   bin.loc<-as.vector(which(unlist(data.check)==2))
    # }
    
    dist.mat<-matrix(0,nrow=nrow(data), ncol=K)
    for(l in 1:nrow(data))
    {
        for(m in 1:K)
        {temp.vec<-vector()
            for(i in 1:ncol(data))
            {
                if(sum(i==bin.loc)!=0)
                {
                    if(!is.null(weights))
                    {
                        temp.vec[i]<-(weights[i]*(ifelse(center[m,i]==data[l,i],0,1)))/sum(weights)
                    }
                    else{temp.vec[i]<-(ifelse(center[m,i]==data[l,i],0,1))}
                    #else{temp.vec[i]<-((data[l,i]-center[m,i])^2)}
                }
                if(sum(i==cont.loc)!=0)
                {
                    if(!is.null(weights))
                    {
                       # temp.vec[i]<-(weights[i]*(abs(data[l,i]-center[m,i])/(max(data[,i])-min(data[,i]))))/sum(weights)
                        temp.vec[i]<-(weights[i]*(data[l,i]-center[m,i])^2)/sum(weights)
                        

                        }
                    #else{temp.vec[i]<-(abs(data[l,i]-center[m,i])/(max(data[,i])-min(data[,i])))}
                    else{temp.vec[i]<-((data[l,i]-center[m,i])^2)}
                }
                if(sum(i==cat.loc)!=0){
                    if(!is.null(weights))
                    {
                        temp.vec[i]<- (weights[i]*(ifelse(center[m,i]==data[l,i],0,1)))/sum(weights)
                    }
                    else{temp.vec[i]<- (ifelse(center[m,i]==data[l,i],0,1))}
                    #else{temp.vec[i]<-((data[l,i]-center[m,i])^2)}
                }
                if(sum(i==ord.loc)!=0)
                {
                    if(!is.null(weights))
                    {
                        temp.vec[i]<-((weights[i]*(abs(center[m,i]-data[l,i]))/(max(data[,i])-min(data[,i])))/sum(weights))^2
                    }
                    #else{temp.vec[i]<-(abs(center[m,i]-data[l,i]))/(max(data[,i])-min(data[,i]))}
                    else{temp.vec[i]<-(center[m,i]-data[l,i])^2}
                }
            }
            #euc.dis<-sqrt(sum(temp.vec[c(cont.loc)]))*(length(cont.loc)/length(label))+sqrt(sum(temp.vec[c(ord.loc)]))*(length(ord.loc)/length(label))
            euc.dis<-sqrt(sum(temp.vec[cont.loc]))*(length(cont.loc)/length(label))+sum(temp.vec[ord.loc])*(length(ord.loc)/length(label))
            cat.dis<-sum(temp.vec[c(cat.loc,bin.loc)])/length(label)
            dist.mat[l,m]<-sum(euc.dis,cat.dis) }
    }
    
    return(dist.mat)
    
}
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

