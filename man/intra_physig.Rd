% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intra_physig.R
\name{intra_physig}
\alias{intra_physig}
\title{Intraspecific variability - Phylogenetic signal}
\usage{
intra_physig(
  trait.col,
  data,
  phy,
  V = NULL,
  n.intra = 100,
  distrib = "normal",
  method = "K",
  track = TRUE
)
}
\arguments{
\item{trait.col}{The name of a column in the provided data frame with trait 
to be analyzed  (e.g. "Body_mass").}

\item{data}{Data frame containing species traits with row names matching tips
in \code{phy}.}

\item{phy}{A phylogeny (class 'phylo', see ?\code{ape}).}

\item{V}{Name of the column containing the standard deviation or the standard error of the trait 
variable. When information is not available for one taxon, the value can be 0 or \code{NA}.}

\item{n.intra}{Number of times to repeat the analysis generating a random trait value.
If NULL, \code{n.intra} = 30}

\item{distrib}{A character string indicating which distribution to use to generate a random value for the response 
and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
Uniform distribution: "uniform" (\code{\link{runif}})
Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.}

\item{method}{Method to compute signal: can be "K" or "lambda".}

\item{track}{Print a report tracking function progress (default = TRUE)}
}
\value{
The function \code{intra_physig} returns a list with the following
components:

\code{Trait}: Column name of the trait analysed

\code{data}: Original full dataset

\code{intra.physig.estimates}: Run number, phylogenetic signal estimate 
(lambda or K) and the p-value for each run with a different simulated datset.

\code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.

\code{stats}: Main statistics for signal estimate\code{CI_low} and \code{CI_high} are the lower 
and upper limits of the 95% confidence interval.
}
\description{
Performs Phylogenetic signal estimates evaluating
trait intraspecific variability
}
\details{
This function estimates phylogenetic signal using \code{\link[phytools]{phylosig}}.
The analysis is repeated \code{n.intra} times. At each iteration the function generates a random value
for each row in the dataset using the standard deviation or errors supplied and assuming a normal or uniform distribution.
To calculate means and se for your raw data, you can use the \code{summarySE} function from the 
package \code{Rmisc}. 

Output can be visualised using \code{sensi_plot}.
}
\note{
The argument "se" from \code{\link[phytools]{phylosig}} is not available in this function. Use the 
argument "V" instead with \code{\link{intra_physig}} to indicate the name of the column containing the standard 
deviation or the standard error of the trait variable instead.
}
\examples{
\dontrun{
data(alien)
alien.data<-alien$data
alien.phy<-alien$phy
# Run sensitivity analysis:
intra <- intra_physig(trait.col = "gestaLen", V = "SD_gesta" ,
                     data = alien.data, phy = alien.phy[[1]])
summary(intra)
sensi_plot(intra)
sensi_plot(intra, graphs = 1)
sensi_plot(intra, graphs = 2)
}

\dontshow{
data(alien)
# Run sensitivity analysis:
intra <- intra_physig(trait.col = "gestaLen", V = "SD_gesta" ,
                     data = alien.data, n.intra = 5,
                     phy = alien.phy[[1]])
summary(intra)
}
}
\references{
Paterno, G. B., Penone, C. Werner, G. D. A. 
\href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
An r-package for sensitivity analysis in phylogenetic 
comparative methods.} Methods in Ecology and Evolution 
2018, 9(6):1461-1467

Martinez, P. a., Zurano, J.P., Amado, T.F., Penone, C., Betancur-R, R., 
Bidau, C.J. & Jacobina, U.P. (2015). Chromosomal diversity in tropical reef 
fishes is related to body size and depth range. Molecular Phylogenetics and 
Evolution, 93, 1-4

Blomberg, S. P., T. Garland Jr., A. R. Ives (2003) 
Testing for phylogenetic signal in comparative data: 
Behavioral traits are more labile. Evolution, 57, 717-745.

Pagel, M. (1999) Inferring the historical patterns of biological evolution. 
Nature, 401, 877-884.

Kamilar, J. M., & Cooper, N. (2013). Phylogenetic signal in primate behaviour,
 ecology and life history. Philosophical Transactions of the Royal Society 
 B: Biological Sciences, 368: 20120341.
}
\seealso{
\code{\link[phytools]{phylosig}}, \code{\link{sensi_plot}}
}
\author{
Caterina Penone & Pablo Ariel Martinez & Gustavo Paterno
}
