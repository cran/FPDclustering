\name{Country_data}
\alias{Country_data}
\docType{data}
\title{Unsupervised Learning on Country Data
%%   ~~ data name/kind ... ~~
}
\description{Ten vables recorded on 167 countries. The goal is to categorize the countries using socio-economic and health indicators that determine the country's overall development. The data set has been donated by the HELP International organization, an international humanitarian NGO that needs to
identify the countries that need aid and asked the analysts to categorize the countries.

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Country_data)}
\format{
 A data frame with 167 observations and 10 variables.
 \describe{
\item{country}{country name}
\item{child_mort}{Death of children under 5 years of age per 1000 live births}
\item{exports}{ Exports of goods and services per capita. Given as \%age of the GDP per capita}
\item{health}{ Total health spending per capita. Given as \%age of GDP per capita}
\item{imports}{ Imports of goods and services per capita. Given as \%age of the GDP per capita}
\item{income}{ Net income per person}
\item{inflation}{The measurement of the annual growth rate of the Total GDP}
\item{life_expec}{The average number of years a new born child would live if the current mortality patterns are to remain the same}
\item{total_fer}{ The number of children that would be born to each woman if the current age-fertility rates remain the same.}
\item{gdpp}{ The GDP per capita. Calculated as the Total GDP divided by the total population.}
}}
\source{https://www.kaggle.com/datasets/rohan0301/unsupervised-learning-on-country-data/metadata?resource=download
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\references{
R. Kokkula. Unsupervised 
learning on country data. kaggle, 2022. URL
https://www.kaggle.com/datasets/rohan0301/unsupervised-learning-on-country-data/metadata?resource=download
}

\examples{
data(Country_data)
pairs(Country_data[,2:10])
}

