\name{tsd_IPD}
\alias{tsd_IPD}
\docType{data}
\title{Set of Individual Patient Data from Appendix D of NICE TSD 18}
\description{A set of IPD created using the \code{wakefield} package and consists of two patient characteristics: age and gender. Patients are randomised to one of two treatments and a binary outcome is reported.}
\usage{data("tsd_IPD")}
\format{
  Patient data set; 500 subjects randomised 50:50 to treatment A or treatment B. Dataset includes patient ID, age, gender and a binary outcome.
  \describe{
    \item{\code{ID}}{a character vector, giving a unique ID number allocated to each patient.}
    \item{\code{age}}{a numeric vector, giving patient age within the range 45-75.}
    \item{\code{gender}}{a factor with levels \code{Male}, \code{Female}, 64\% patients classed as female.}
    \item{\code{trt}}{a character vector, defining the treatment to which the patient was assigned, A or B.}
    \item{\code{y}}{a binary variable indicating whether the patient experienced the event or not.}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Phillippo, D.M., Ades, A.E., Dias, S., Palmer, S., Abrams, K.R., Welton, N.J. NICE DSU Technical Support Document 18: Methods for population-adjusted indirect comparisons in submission to NICE. 2016. Available from http://www.nicedsu.org.uk
}
\examples{
data(tsd_IPD)
}
\keyword{datasets}
