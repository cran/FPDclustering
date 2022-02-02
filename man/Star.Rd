\name{Star}
\alias{Star}
\docType{data}
\title{Star dataset to predict star types
%%   ~~ data name/kind ... ~~
}
\description{A 6 class star dataset for star classification with Deep Learned approaches
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(ais)}
\format{
 A data frame with 202 observations and 13 variable.
 \describe{
\item{K}{ Absolute Temperature (in K)}
\item{Lum}{ Relative Luminosity (L/Lo)}
\item{Rad}{ Relative Radius (R/Ro)}
\item{Mag}{ Absolute Magnitude (Mv)}
\item{Col}{ Star Color (white,Red,Blue,Yellow,yellow-orange etc)}
\item{Spect}{ Spectral Class (O,B,A,F,G,K,,M)}
\item{Type}{ Star Type (Red Dwarf, Brown Dwarf, White Dwarf, Main Sequence , SuperGiants, HyperGiants)}
}
}


\source{ https://www.kaggle.com/deepu1109/star-dataset
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}


\examples{
data(Star)

}

