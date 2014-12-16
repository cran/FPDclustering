\name{asymmetric3}
\alias{asymmetric3}
\docType{data}
\title{Asymmetric data set shape=3
%%   ~~ data name/kind ... ~~
}
\description{Each cluster has been generated according to a multivariate asymmetric Gaussian distribution, with shape=3, covariance matrix equal to the identity matrix and randomly generated centres.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(asymmetric3)}
\format{
  A data frame with 800 observations on 101 variables.
 The first variable is the membership labels.
}

\source{Generated with R using the package sn (The skew-normal and skew-t distributions), function rsn
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(asymmetric3)
plot(asymmetric3[,2:3])
}

