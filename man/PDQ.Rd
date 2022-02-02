\name{PDQ}
\alias{PDQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probabilistic Distance Clustering Adjusted for Cluster Size 
}
\description{ An implementation of probabilistic distance clustering adjusted for cluster size (PDQ), a probabilistic distance clustering algorithm that involves optimizing the PD-clustering criterion. The algorithm can be used, on continous, count, or mixed type data setting Euclidean, Chi square, or Gower  as dissimilarity measurments. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PDQ(x=NULL,k=2,ini='kmd',dist='euc',cent=NULL,ord=NULL,cat=NULL,bin=NULL,cont=NULL,w=NULL)
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{data} here~~
}
  \item{k}{A numerical parameter giving the number of clusters.
%%     ~~Describe \code{k} here~~
}
  \item{ini}{A parameter that selects center starts. Options available are random ("random"), kmedoid ("kmd", by default"), center ("center", the user inputs the center), and kmode ("kmode", for categoriacal data sets).
%%     ~~Describe \code{method} here~~
}
  \item{dist}{A parameter that selects the distance measure used. Options available are Eucledean ("euc"), Gower ("gower") and chi square ("chi").
%%     ~~Describe \code{method} here~~
}
  \item{cent}{User inputed centers if method selected is "random".
%%     ~~Describe \code{cent} here~~
}  
\item{ord}{column indices of the x matrix indicating which columns are ordinal variables.
%%     ~~Describe \code{method} here~~
}
\item{cat}{column indices of the x matrix indicating which columns are  categorical variables.
%%     ~~Describe \code{cent} here~~
}  
\item{bin}{column indices of the x matrix indicating which columns are  binary variables.
%%     ~~Describe \code{cent} here~~
}  
\item{cont}{column indices of the x matrix indicating which columns are  continuous variables.
}
\item{w}{numerical vector same length as the columns of the data, ccontaining the variable weights when using Gower distance, equal weights  by default.
}
%%     ~~Describe \code{cent} here~~

}
\value{ A class FPDclustering list with components
%%  ~Describe the value returned
%%  If it is a LIST, use label=l, centers=c, probability=p, JDF=JDF, JDFIter=JDFv, iter=iter, explained
 \item{label }{A vector of integers indicating the cluster membership for each unit}
 \item{centers }{A matrix of cluster centers}
 \item{probability }{A matrix of probability of each point belonging to each cluster}
 \item{JDF }{The value of the Joint distance function}
 \item{iter}{The number of iterations}
 \item{jdfvector }{collection of all jdf calculations at each iteration}
  \item{data }{the data set}
%% ...
}
\references{
	Iyigun, Cem, and Adi Ben-Israel. \emph{Probabilistic distance clustering adjusted for cluster size.} Probability in the Engineering and Informational Sciences 22.4 (2008): 603-621.
doi.org/10.1017/S0269964808000351.
	
	Tortora and Palumbo. \emph{Clustering mixed-type data using a probabilistic distance algorithm.} submitted.
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora and Noe Vidales
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{PDC}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#Mixed type data

sig=matrix(0.7,4,4)
diag(sig)=1###creat a correlation matrix
x1=rmvnorm(200,c(0,0,3,3))##  cluster 1
x2=rmvnorm(200,c(4,4,6,6),sigma=sig)##  cluster 2
x=rbind(x1,x2)# data set with 2 clusters
l=c(rep(1,200),rep(2,200))#creating the labels
x1=cbind(x1,rbinom(200,4,0.2),rbinom(200,4,0.2))#categorical variables
x2=cbind(x2,rbinom(200,4,0.7),rbinom(200,4,0.7))
x=rbind(x1,x2) ##Data set

#### Performing PDQ
pdq_class<-PDQ(x=x,k=2, ini="random", dist="gower", cont= 1:4, cat = 5:6)

###Output
table(l,pdq_class$label)
plot(pdq_class)
summary(pdq_class)



###Continuous data example
# Gaussian Generated Data  no  overlap 
x<-rmvnorm(100, mean=c(1,5,10), sigma=diag(1,3))
y<-rmvnorm(100, mean=c(4,8,13), sigma=diag(1,3))
data<-rbind(x,y)

#### Performing PDQ
pdq1=PDQ(data,2,ini="random",dist="euc")
table(rep(c(2,1),each=100),pdq1$label)
Silh(pdq1$probability)
plot(pdq1)
summary(pdq1)


# Gaussian Generated Data with  overlap 
x2<-rmvnorm(100, mean=c(1,5,10), sigma=diag(1,3))
y2<-rmvnorm(100, mean=c(2,6,11), sigma=diag(1,3))
data2<-rbind(x2,y2)

#### Performing PDQ
pdq2=PDQ(data2,2,ini="random",dist="euc")
table(rep(c(1,2),each=100),pdq2$label)
plot(pdq2)
summary(pdq2)
}

