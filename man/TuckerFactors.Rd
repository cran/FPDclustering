\name{TuckerFactors}
\alias{TuckerFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Choice of the number of Tucker 3 factors
%%  ~~function to do ... ~~
}
\description{ An empirical way of choosing the number of factors. The algorithm returns a graph and a table representing the explained variability varying the number of factors.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TuckerFactors(data = NULL, nc = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A matrix or data frame such that rows correspond to observations and columns correspond to variables.
%%     ~~Describe \code{data} here~~
}
  \item{nc}{A numerical parameter giving the number of clusters
%%     ~~Describe \code{nc} here~~
}
}

\value{A table containing the explained variability varying the number of factors for units (column) and for variables (row) and a plot
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ 
Kiers H, Kinderen A. A fast method for choosing the numbers of components in
Tucker3 analysis.\emph{British Journal of Mathematical and Statistical Psychology}, \bold{56}(1), 119-125, 2003.

Kroonenberg P. \emph{Applied Multiway Data Analysis}. Ebooks Corporation, Hoboken, New Jersey, 2008.

Tortora C.,  Gettler Summa M., and  Palumbo F..
Factor pd-clustering. In Lausen et al., editor, \emph{Algorithms from and for Nature and Life, Studies in Classification}, Data Analysis, and Knowledge Organization DOI 10.1007/978-3-319-00035-011, 115-123, 2013.
%% ~put references to the literature/web site here ~
}
\author{Cristina Tortora
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[ThreeWay]{T3}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
# Asymmetric data set example (with shape=3).
data('asymmetric3')
xp=TuckerFactors(asymmetric3[,-1], nc = 4)

}

\dontrun{
# Asymmetric data set example (with shape=20).
data('asymmetric20')
xp=TuckerFactors(asymmetric20[,-1], nc = 4)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

