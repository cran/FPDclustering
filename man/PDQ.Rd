\name{PDQ}
\alias{PDQ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probabilistic Clustering Adjusted for Cluster Size 
}
\description{ An implementation of PDQ, a probabilistic distance clustering algorithm that involves optimizing the PD-clustering criterion with the option  of Euclidean and Chi as dissimilarity measurments. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PDQ(data=NULL,K=2,method="random", distance="euc", cent=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{data} here~~
}
  \item{K}{A numerical parameter giving the number of clusters
%%     ~~Describe \code{K} here~~
}
  \item{method}{A parameter that selects center starts. Options available are random, kmedoid, and center(user inputs center starts)
%%     ~~Describe \code{method} here~~
}
  \item{distance}{A parameter that selects the distance measure used. Options available are Eucledean euc  and chi square chi
%%     ~~Describe \code{method} here~~
}
  \item{cent}{User inputed centers if method selected is "random"
%%     ~~Describe \code{cent} here~~
}
}
\value{ A list with components
%%  ~Describe the value returned
%%  If it is a LIST, use label=l, centers=c, probability=p, JDF=JDF, JDFIter=JDFv, iter=iter, explained
 \item{label }{A vector of integers indicating the cluster membership for each unit}
 \item{centers }{A matrix of cluster centers}
 \item{probability }{A matrix of probability of each point belonging to each cluster}
 \item{JDF }{The value of the Joint distance function}
 \item{iter}{The number of iterations}
 \item{jdfvector }{collection of all jdf calculations at each iteration}
%% ...
}
\references{
	Iyigun, Cem, and Adi Ben-Israel. \emph{Probabilistic distance clustering adjusted for cluster size.} Probability in the Engineering and Informational Sciences 22.4 (2008): 603-621.
doi.org/10.1017/S0269964808000351.

%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora and Noe Vidales
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{PDclust}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Gaussian Generated Data  no  overlap 
x<-rmvnorm(100, mean=c(1,5,10), sigma=diag(1,3))
y<-rmvnorm(100, mean=c(4,8,13), sigma=diag(1,3))
data<-rbind(x,y)
pdq1=PDQ(data,2,method="random",distance="euc")
table(rep(c(2,1),each=100),pdq1$label)
Silh(pdq1$probability)



# Gaussian Generated Data with  overlap 
x2<-rmvnorm(100, mean=c(1,5,10), sigma=diag(1,3))
y2<-rmvnorm(100, mean=c(2,6,11), sigma=diag(1,3))
data2<-rbind(x2,y2)
pdq2=PDQ(data2,2,method="random",distance="euc")
table(rep(c(1,2),each=100),pdq2$label)
Silh(pdq2$probability)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
