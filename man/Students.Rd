\name{Students}
\alias{Students}
\docType{data}
\title{Statistics 1 students 
%%   ~~ data name/kind ... ~~
}
\description{Data set collected in 2022 that contains 10 variables recorded on a convenience sample of 253 students enrolled in the first year at the University od Naples FedericoII and attending an introductory Statistics course.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data(Students)}
\format{
 A data frame with 253 observations and 10 variable.
 \describe{
\item{Sex}{ gender, binary}
\item{HS_qual}{ high school type, categorical}
\item{Stud_stat}{ prior knowladge of statistics, binary}
\item{Course_modality}{ course modality of attendance (in presence, online, mixed), categorical}
\item{HE_Parents}{parents' education degree, categorical}
\item{PMP}{ mathematical prerequisits for psychometric, continuous}
\item{SAS}{ statistical anxiety sale, continuous}
\item{RAI}{ relative authonomy index, continuous}
\item{S_EFF}{ self-efficacy, continuous}
\item{COG}{ cognitive competence, continuous}
}
}


\references{ R. Fabbricatore. Latent class analysis for proficiency assessment in higher education: integrating multidimensional latent traits and learning topics. Ph.D. thesis, University of Naples Federico II, 2023
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}


\examples{
data(Students)

}

