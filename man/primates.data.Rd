% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/document-data.R
\docType{data}
\name{primates.data}
\alias{primates.data}
\title{Primates dataset: Example dataset for the package sensiPhy}
\format{A data frame with 95 rows and 3 variables:
\itemize{
  \item family: Taxonomic family
  \item adultMass: Mean adult body mass (g)
  \item sexMaturity: Age when individuals are first physically capable of reproducing (days)
  \item homeRange: Mean home range (km)
  }}
\usage{
data(primates)
}
\description{
A comparative dataset containing traits for 95 Primates species 
(primates.data) and a multiphylo object with 101 phylogenies matching the 
data (primates.phy). Tip labels are the binomial species names and match 
with data rownames. Data was taken from (Jones et al. 2009) and phylogenies from
(Fritz et al 2009) and (Kuhn et al 2011).
}
\references{
Data downloaded from: http://esapubs.org/archive/ecol/E090/184/

Jones, K. E., Bielby, J., Cardillo, M., Fritz, S. A., O'Dell, J., Orme, C. D. L., Safi, K., Sechrest,
 W., Boakes, E. H., Carbone, C., Connolly, C., Cutts, M. J., Foster, J. K., Grenyer, R., Habib, M., Plaster, 
 C. A., Price, S. A., Rigby, E. A., Rist, J., Teacher, A., Bininda-Emonds, O. R. P., Gittleman, J. L., 
 Mace, G. M., Purvis, A. (2009), PanTHERIA: a species-level database of life history, ecology, 
 and geography of extant and recently extinct mammals. Ecology, 90: 2648. doi: 10.1890/08-1494.1

Phylogeny: Kuhn, Tyler S., Arne O. Mooers, and Gavin H. Thomas. "A simple polytomy resolver for 
dated phylogenies." Methods in Ecology and Evolution 2.5 (2011): 427-436.

Fritz, Susanne A., Olaf RP Bininda-Emonds, and Andy Purvis. "Geographical variation in predictors 
of mammalian extinction risk: big is bad, but only in the tropics." Ecology letters 12.6 (2009): 538-549.
}
\keyword{datasets}
