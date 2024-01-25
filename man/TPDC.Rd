\name{TPDC}
\alias{TPDC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Student-t PD-Clustering
}
\description{ An implementation of Student-t PD-Clustering TPDC, an extention of PD-clustering adjusted for cluster size that uses a dissimilarity measure based on the multivariate Student-t density. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TPDC(data=NULL,k=2,ini="kmedoids", nr=5,iter=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{data} here~~
}
  \item{k}{A numerical parameter giving the number of clusters
%%     ~~Describe \code{K} here~~
}
  \item{ini}{A parameter that selects center starts. Options available are random ("random"), kmedoid ("kmedoid", by default), and PDC ("PDclust").
%%     ~~Describe \code{method} here~~
}
  \item{nr}{Number of random starts if ini is "random"
%%     ~~Describe \code{method} here~~
}
  \item{iter}{Maximum number of iterations
%%     ~~Describe \code{method} here~~
}
}
\value{ A class FPDclustering list with components
%%  ~Describe the value returned
%%  If it is a LIST, use label=l, centers=c, probability=p, JDF=JDF, JDFIter=JDFv, iter=iter, explained
 \item{label }{A vector of integers indicating the cluster membership for each unit}
 \item{centers }{A matrix of cluster means}
  \item{sigma }{A list of K elements, with the variance-covariance matrix per cluster}
    \item{df }{A vector of K degrees of freedom }
 \item{probability }{A matrix of probability of each point belonging to each cluster}
 \item{JDF }{The value of the Joint distance function}
 \item{iter}{The number of iterations}
   \item{data }{the data set}
%% ...
}
\references{
		Tortora C., McNicholas P.D.,  and Palumbo F. \emph{A probabilistic distance clustering algorithm using Gaussian and Student-t multivariate density distributions.} SN Computer Science, 1:65, 2020.

C. Rainey, C. Tortora and F.Palumbo.\emph{ A parametric version of probabilistic distance clustering}. In: Greselin F., Deldossi L., Bagnato L., Vichi M. (eds) Statistical Learning of Complex Data. CLADAG 2017. Studies in Classification, Data Analysis, and Knowledge Organization. Springer, Cham, 33-43 2019. doi.org/10.1007/978-3-030-21140-0_4
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora and Francesco Palumbo
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{PDC},\link{PDQ}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Load the data
data(ais)
dataSEL=ais[,c(10,3,5,8)]

#Clustering
res=TPDC(dataSEL,k=2,ini = "kmedoids")

#Results
table(res$label,ais$sex)
summary(res)
plot(res)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
