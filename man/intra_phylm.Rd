% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intra_phylm.R
\name{intra_phylm}
\alias{intra_phylm}
\title{Intraspecific variability - Phylogenetic Linear Regression}
\usage{
intra_phylm(
  formula,
  data,
  phy,
  Vy = NULL,
  Vx = NULL,
  y.transf = NULL,
  x.transf = NULL,
  n.intra = 30,
  distrib = "normal",
  model = "lambda",
  track = TRUE,
  ...
)
}
\arguments{
\item{formula}{The model formula: \code{response~predictor}.}

\item{data}{Data frame containing species traits and species names as row names.}

\item{phy}{A phylogeny (class 'phylo', see ?\code{ape}).}

\item{Vy}{Name of the column containing the standard deviation or the standard error of the response 
variable. When information is not available for one taxon, the value can be 0 or \code{NA}.}

\item{Vx}{Name of the column containing the standard deviation or the standard error of the predictor 
variable. When information is not available for one taxon, the value can be 0 or \code{NA}}

\item{y.transf}{Transformation for the response variable (e.g. \code{log} or \code{sqrt}). Please use this 
argument instead of transforming data in the formula directly (see also details below).}

\item{x.transf}{Transformation for the predictor variable (e.g. \code{log} or \code{sqrt}). Please use this 
argument instead of transforming data in the formula directly (see also details below).}

\item{n.intra}{Number of times to repeat the analysis generating a random value for response and/or predictor variables.
If NULL, \code{n.intra} = 30}

\item{distrib}{A character string indicating which distribution to use to generate a random value for the response 
and/or predictor variables. Default is normal distribution: "normal" (function \code{\link{rnorm}}).
Uniform distribution: "uniform" (\code{\link{runif}})
Warning: we recommend to use normal distribution with Vx or Vy = standard deviation of the mean.}

\item{model}{The phylogenetic model to use (see Details). Default is \code{lambda}.}

\item{track}{Print a report tracking function progress (default = TRUE)}

\item{...}{Further arguments to be passed to \code{phylolm}}
}
\value{
The function \code{intra_phylm} returns a list with the following
components:

\code{formula}: The formula

\code{data}: Original full dataset

\code{sensi.estimates}: Coefficients, aic and the optimised
value of the phylogenetic parameter (e.g. \code{lambda}) for each regression.

\code{N.obs}: Size of the dataset after matching it with tree tips and removing NA's.

\code{stats}: Main statistics for model parameters.\code{CI_low} and \code{CI_high} are the lower 
and upper limits of the 95% confidence interval.

\code{all.stats}: Complete statistics for model parameters. \code{sd_intra} is the standard deviation 
due to intraspecific variation. \code{CI_low} and \code{CI_high} are the lower and upper limits 
of the 95% confidence interval.

\code{sp.pb}: Species that caused problems with data transformation (see details above).
}
\description{
Performs Phylogenetic linear regression evaluating
intraspecific variability in response and/or predictor variables.
}
\details{
This function fits a phylogenetic linear regression model using \code{\link[phylolm]{phylolm}}.
The regression is repeated \code{n.intra} times. At each iteration the function generates a random value
for each row in the dataset using the standard deviation or errors supplied and assuming a normal or uniform distribution.
To calculate means and se for your raw data, you can use the \code{summarySE} function from the 
package \code{Rmisc}.

#' All phylogenetic models from \code{phylolm} can be used, i.e. \code{BM},
\code{OUfixedRoot}, \code{OUrandomRoot}, \code{lambda}, \code{kappa},
\code{delta}, \code{EB} and \code{trend}. See ?\code{phylolm} for details.

Currently, this function can only implement simple linear models (i.e. \eqn{trait~
predictor}). In the future we will implement more complex models.

Output can be visualised using \code{sensi_plot}.
}
\section{Warning}{
  
When Vy or Vx exceed Y or X, respectively, negative (or null) values can be generated, this might cause problems
for data transformation (e.g. log-transformation). In these cases, the function will skip the simulation. This problem can
be solved by increasing \code{n.intra}, changing the transformation type and/or checking the target species in output$sp.pb.
}

\examples{
# Load data:
data(alien)
# run PGLS accounting for intraspecific variation:
intra <- intra_phylm(gestaLen ~ adultMass, y.transf = log, x.transf = log, 
phy = alien$phy[[1]], data = alien$data, Vy = "SD_gesta", n.intra = 30)
# To check summary results:
summary(intra)
# Visual diagnostics
sensi_plot(intra)

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

Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
}
\seealso{
\code{\link[phylolm]{phylolm}}, \code{\link{sensi_plot}}
}
\author{
Caterina Penone & Pablo Ariel Martinez
}
