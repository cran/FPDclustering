\name{Silh}
\alias{Silh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probabilistic silhouette plot
%%  ~~function to do ... ~~
}
\description{ Graphical tool to evaluate the clustering partition.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Silh(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ A matrix of probabilities such that rows correspond to observations and columns correspond to clusters.
%%     ~~Describe \code{p} here~~
}
}
\details{The probabilistic silhouettes are an adaptation of the ones proposed by Menardi(2011) according to the following formula:
\deqn{dbs_i =   (log(p_{im_k}/p_{im_1}))/max_i |log(p_{im_k}/p_{im_1})|}

where \eqn{m_k} is such that \eqn{x_i} belongs to cluster \eqn{k} and \eqn{m_1} is such that \eqn{p_{im_1}} is maximum for \eqn{m} different from\eqn{m_k}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Probabilistic silhouette plot
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{ Menardi G. Density-based Silhouette diagnostics for clustering methods.\emph{Statistics and Computing}, \bold{21}, 295-308, 2011.
%% ~put references to the literature/web site here ~
}
\author{ Cristina Tortora
%%  ~~who you are~~
}

\examples{
\dontrun{
# Asymmetric data set silhouette example (with shape=3).
data('asymmetric3')
x<-asymmetric3[,-1]
fpdas3=FPDC(x,4,3,3)
Silh(fpdas3$probability)
}

\dontrun{
# Asymmetric data set shiluette example (with shape=20).
data('asymmetric20')
x<-asymmetric20[,-1]
fpdas20=FPDC(x,4,3,3)
Silh(fpdas20$probability)
}

\dontrun{
# Shiluette example with outliers.
data('outliers')
x<-outliers[,-1]
fpdout=FPDC(x,4,4,3)
Silh(fpdout$probability)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

