% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sensi_plot.sensiInflu.R
\name{sensi_plot.sensiInflu}
\alias{sensi_plot.sensiInflu}
\title{Graphical diagnostics for class 'sensiInflu'}
\usage{
\method{sensi_plot}{sensiInflu}(x, graphs = "all", param = "estimate", ...)
}
\arguments{
\item{x}{output from \code{influ_phylm}}

\item{graphs}{choose which graph should be printed on the output ("all", 1,2,3 or 4)}

\item{param}{choose which parameter ("intercept" or "estimate" should be printed)}

\item{...}{further arguments to methods}
}
\description{
\code{plot_influ_phylm} Plot results from \code{influ_phylm} and 
\code{influ_phyglm}
}
\details{
For 'x' from influ_phylm or influ_phyglm:

\strong{Graph 1:} Distribution of estimated slopes (estimates) or intercepts for each 
simulation (leave-one-out deletion). Red vertical line represents the original
slope or intercept from the full model (with all species). 

\strong{Graph 2:} Original regression plot (\eqn{trait~predictor}). Standardized 
difference in slope or intercept is represented by a continuous size scale. 
The names of the most influential species (sDF > cutoff) are ploted in the
graph. 

\strong{Graph 3:} Distribution of standardized difference in slope or intercept. Red 
colour indicates influential species (with a standardised difference above 
the value of \code{cutoff}).

\strong{Graph 4:} Distribution of the percentage of change in slope or intercept.
}
\seealso{
\code{\link[ggplot2]{ggplot}}
}
\author{
Gustavo Paterno
}
