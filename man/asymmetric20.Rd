\name{asymmetric20}
\alias{asymmetric20}
\docType{data}
\title{Asymmetric data set shape=20
%%   ~~ data name/kind ... ~~
}
\description{Each cluster has been generated according to a multivariate asymmetric Gaussian distribution, with shape=20, covariance matrix equal to the identity matrix and randomly generated centres.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(asymmetric20)}
\format{
  A data frame with 800 observations on the following 101 variables.
 The first variable is the membership.
 
}

\source{Generated with R using the package sn (The skew-normal and skew-t distributions), function rsn
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(asymmetric20)
plot(asymmetric20[,2:3])
}

