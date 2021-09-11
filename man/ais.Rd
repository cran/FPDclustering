\name{ais}
\alias{ais}
\docType{data}
\title{Australian institute of sport data
%%   ~~ data name/kind ... ~~
}
\description{Data obtained to study sex, sport and body-size dependency of hematology in highly trained athletes.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(ais)}
\format{
 A data frame with 202 observations and 13 variable.
 \describe{
\item{rcc}{ red blood cell count, in}
\item{wcc}{ while blood cell count, in per liter}
\item{hc}{ hematocrit, percent}
\item{hg}{ hemaglobin concentration, in g per decaliter}
\item{ferr}{ plasma ferritins, ng}
\item{bmi}{ Body mass index, kg}
\item{ssf}{ sum of skin folds}
\item{pcBfat}{ percent Body fat}
\item{lbm}{ lean body mass, kg}
\item{ht}{ height, cm}
\item{wt}{ weight, kg}
\item{sex}{ a factor with levels f m}
\item{sport}{ a factor with levels B_Ball Field Gym Netball Row Swim T_400m T_Sprnt Tennis W_Polo}
}}
\source{R package DAAG
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\references{
Telford, R.D. and Cunningham, R.B. 1991. Sex, sport and body-size dependency of hematology in highly trained athletes. Medicine and Science in Sports and Exercise 23: 788-794.}

\examples{
data(ais)
pairs(ais[,1:11],col=ais$sex)
}

