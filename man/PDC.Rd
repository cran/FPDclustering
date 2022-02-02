\name{PDC}
\alias{PDC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probabilistic Distance Clustering
%%  ~~function to do ... ~~
}
\description{Probabilistic distance clustering (PD-clustering) is an iterative, distribution free, probabilistic clustering method. PD clustering is based on the constraint that the product of the probability and the distance of each point to any cluster centre is a constant.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PDC(data = NULL, k = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{data} here~~
}
  \item{k}{A numerical parameter giving the number of clusters
%%     ~~Describe \code{k} here~~
}
}

\value{ A class FPDclustering list with components
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{label }{A vector of integers indicating the cluster membership for each unit}
 \item{centers }{A matrix of cluster centers}
 \item{probability }{A matrix of probability of each point belonging to each cluster}
 \item{JDF }{The value of the Joint distance function}
 \item{iter}{The number of iterations}
  \item{data }{the data set}
}
\references{Ben-Israel C. and  Iyigun C. Probabilistic D-Clustering.\emph{ Journal of Classification}, \bold{25}(1), 5-26, 2008.
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora and Paul D. McNicholas
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

#Normally generated clusters
c1 = c(+2,+2,2,2)
c2 = c(-2,-2,-2,-2)
c3 = c(-3,3,-3,3)
n=200
x1 = cbind(rnorm(n, c1[1]), rnorm(n, c1[2]), rnorm(n, c1[3]), rnorm(n, c1[4]) )
x2 = cbind(rnorm(n, c2[1]), rnorm(n, c2[2]),rnorm(n, c2[3]), rnorm(n, c2[4]) )
x3 = cbind(rnorm(n, c3[1]), rnorm(n, c3[2]),rnorm(n, c3[3]), rnorm(n, c3[4]) )
x = rbind(x1,x2,x3)

#Clustering
pdn=PDC(x,3)

#Results
plot(pdn)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
