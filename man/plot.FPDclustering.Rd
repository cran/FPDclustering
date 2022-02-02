\name{plot.FPDclustering}
\alias{plot.FPDclustering}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots for FPDclusteringt Objects
}
\description{
Probability Silhouette plot, Scatterplot up to 10 varaibles, and  parallel coordinate plot up to 10 variables, for objects of class FPDclustering.
}
\usage{
\method{plot}{FPDclustering}(x, maxVar=30, ... ) 
}

\arguments{
  \item{x}{
an object of class FPDclustering
}
  \item{maxVar}{
a scalar indicating the maximum number of variables to display on the parallel plot, 30 by default
}
\item{...}{Additional parameters for the function paris}
 
}

\author{
Cristina Tortora
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

