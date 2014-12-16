\name{outliers}
\alias{outliers}
\docType{data}
\title{Data set with outliers
%%   ~~ data name/kind ... ~~
}
\description{Each cluster has been generated according to a multivariate Gaussian distribution, with centers c randomly generated. For each cluster, 20\% of uniform distributed outliers have been generated at a distance included in max(x-c) and max(x-c)+5 form the center.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(outliers)}
\format{
  A data frame with 960 observations on the following 101 variables.
  The first variable corresponds to the membership
}

\source{generated with R
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(outliers)
 plot(outliers[,2:3]) 
}

