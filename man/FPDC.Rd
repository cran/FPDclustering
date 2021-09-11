\name{FPDC}
\alias{FPDC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Factor probabilistic distance clustering
}
\description{ An implementation of FPDC, a probabilistic factor clustering algorithm that involves a linear transformation of variables and a cluster optimizing the PD-clustering criterion
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FPDC(data = NULL, k = 2, nf = 2, nu = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{data} here~~
}
  \item{k}{A numerical parameter giving the number of clusters
%%     ~~Describe \code{k} here~~
}
  \item{nf}{A numerical parameter giving the number of factors for variables
%%     ~~Describe \code{nf} here~~
}
  \item{nu}{A numerical parameter giving the number of factors for units
%%     ~~Describe \code{nu} here~~
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
 \item{explained }{The explained variability}
%% ...
}
\references{
	Tortora, C., M. Gettler Summa, M. Marino, and F. Palumbo. \emph{Factor probabilistic distance clustering
		(fpdc): a new clustering method for high dimensional data sets}. Advanced in Data Analysis and Classification,  10(4), 441-464, 2016. doi:10.1007/s11634-015-0219-5.
	
	Tortora C.,  Gettler Summa M., and  Palumbo F..
Factor pd-clustering. In Lausen et al., editor, \emph{Algorithms from and for Nature and Life, Studies in Classification}, Data Analysis, and Knowledge Organization DOI 10.1007/978-3-319-00035-011, 115-123, 2013.

Tortora C., \emph{Non-hierarchical clustering methods on factorial subspaces}, 2012.

%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora and Paul D. McNicholas
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{PDC}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
# Asymmetric data set clustering example (with shape=3).
data('asymmetric3')
x<-asymmetric3[,-1]

#Clustering
fpdas3=FPDC(x,4,3,3)

#Results
table(asymmetric3[,1],fpdas3$label)
Silh(fpdas3$probability)
summary(fpdas3)
plot(fpdas3)
}

\dontrun{
# Asymmetric data set clustering example (with shape=20).
data('asymmetric20')
x<-asymmetric20[,-1]

#Clustering
fpdas20=FPDC(x,4,3,3)

#Results
table(asymmetric20[,1],fpdas20$label)
Silh(fpdas20$probability)
summary(fpdas20)
plot(fpdas20)
}

\dontrun{
# Clustering example with outliers.
data('outliers')
x<-outliers[,-1]

#Clustering
fpdout=FPDC(x,4,5,4)

#Results
table(outliers[,1],fpdout$label)
Silh(fpdout$probability)
summary(fpdout)
plot(fpdout)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
