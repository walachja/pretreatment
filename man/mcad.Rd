\docType{data}
\name{mcad}
\alias{mcad}
\title{metabolomics mcad data set}
\format{A data frame with 50 observations and 279 variables}
\usage{
data(mcad)
}
\description{
The aim of the experiment was to ascertain novel biomarkers of
MCAD (Medium chain acyl-CoA dehydrogenase) deficiency.
The data consists of 25 patients and 25 controls and the analysis was done by LC-MS.
Rows represent patients and controls and columns represent chemical
entities with their quantity.
}
\details{
\itemize{
 \item{\code{group }}{patient group}
 \item{\code{... }}{the remaining variables columns are represented by m/z which are chemical characterizations of individual chemical components on exact mass measurements..}
}
}
\examples{

data(mcad)
str(mcad)
}
\references{
Najdekr L., Gardlo A., Madrova L., Friedeckyy D., Janeckova H., Correa E.S., Goodacre R., Adam T., Oxidized phosphatidylcholines suggest oxidative stress in patients with medium-chain acyl-CoA dehydrogenase deficiency, \emph{Talanta} 139, 2015, 62-66.
}
\keyword{data}
