#' Interaction between phylogenetic uncertainty and influential clade detection - Phylogenetic Logistic Regression
#'
#' Estimate the impact on model estimates of phylogenetic logistic regression after 
#' removing clades from the analysis and evaluating uncertainty in trees topology.
#'
#' @param formula The model formula
#' @param data Data frame containing species traits with row names matching tips
#' in \code{phy}.
#' @param phy A phylogeny (class 'multiPhylo', see ?\code{ape}).
#' @param clade.col The column in the provided data frame which specifies the
#' clades (a character vector with clade names).
#' @param n.species Minimum number of species in a clade for the clade to be
#' included in the leave-one-out deletion analysis. Default is \code{5}.
#' @param n.sim Number of simulations for the randomization test.
#' @param n.tree Number of times to repeat the analysis with n different trees picked 
#' randomly in the multiPhylo file.
#' If NULL, \code{n.tree} = 2
#' @param btol Bound on searching space. For details see \code{phyloglm}
#' @param track Print a report tracking function progress (default = TRUE)
#' @param ... Further arguments to be passed to \code{phyloglm}
#' @details
#' Currently only logistic regression using the "logistic_MPLE"-method from
#' \code{phyloglm} is implemented.
#' 
#' This function sequentially removes one clade at a time, fits a phylogenetic
#' logistic regression model using \code{\link[phylolm]{phyloglm}} and stores the
#' results. The impact of of a specific clade on model estimates is calculated by the
#' comparison between the full model (with all species) and the model without 
#' the species belonging to a clade. It repeats this operation using n trees, 
#' randomly picked in a multiPhylo file.
#' 
#' Additionally, to account for the influence of the number of species on each 
#' clade (clade sample size), this function also estimates a null distribution of slopes
#' expected for the number of species in a given clade. This is done by fitting
#' models without the same number of species in the given clade. 
#'  The number of simulations to be performed is set by 'n.sim'. To test if the 
#'  clade influence differs from the null expectation for a clade of that size, 
#'  a randomization test can be performed using 'summary(x)'. 
#'
#'
#' \code{clade_phyglm} detects influential clades based on
#' difference in intercept and/or slope when removing a given clade compared
#' to the full model including all species. This is done for n trees in the multiphylo file.
#' 
#' Currently, this function can only implements simple logistic models (i.e. \eqn{trait~
#' predictor}). In the future we will implement more complex models.
#'
#' Output can be visualised using \code{sensi_plot}.
#'
#' @return The function \code{clade_phyglm} returns a list with the following
#' components:
#' @return \code{formula}: The formula
#' @return \code{full.model.estimates}: Coefficients, aic and the optimised
#' value of the phylogenetic parameter (e.g. \code{lambda}) for the full model
#' without deleted species.
#' @return \code{sensi.estimates}: A data frame with all simulation
#' estimates. Each row represents a deleted clade for a tree iteration. Columns report the calculated
#' regression intercept (\code{intercept}), difference between simulation
#' intercept and full model intercept (\code{DIFintercept}), the percentage of change
#' in intercept compared to the full model (\code{intercept.perc}) and intercept
#' p-value (\code{pval.intercept}). All these parameters are also reported for the regression
#' slope (\code{DIFestimate} etc.). Additionally, model aic value (\code{AIC}) and
#' the optimised value (\code{optpar}) of the phylogenetic parameter 
#' (e.g. \code{kappa} or \code{lambda}, depending on the phylogenetic model used) 
#' are reported.
#' @return \code{null.dist}: A data frame with estimates for the null distributions
#' for all clades analysed.
#' @return \code{data}: Original full dataset.
#' @return \code{errors}: Clades and/or trees where deletion resulted in errors.
#' @author Gustavo Paterno, Caterina Penone & Gijsbert D.A. Werner
#' @seealso \code{\link[phylolm]{phyloglm}}, \code{\link[sensiPhy]{tree_phyglm}},
#'  \code{\link{clade_phyglm}}, \code{\link{tree_clade_phylm}},
#' \code{\link{sensi_plot}}
#' @references 
#' 
#' Paterno, G. B., Penone, C. Werner, G. D. A. 
#' \href{http://doi.wiley.com/10.1111/2041-210X.12990}{sensiPhy: 
#' An r-package for sensitivity analysis in phylogenetic 
#' comparative methods.} Methods in Ecology and Evolution 
#' 2018, 9(6):1461-1467
#'
#' Ho, L. S. T. and Ane, C. 2014. "A linear-time algorithm for 
#' Gaussian and non-Gaussian trait evolution models". Systematic Biology 63(3):397-408.
#' @examples 
#' \dontrun{
#' # Simulate Data:
#' set.seed(6987)
#' mphy = rmtree(150, N = 30)
#' x = rTrait(n=1,phy=mphy[[1]])
#' X = cbind(rep(1,150),x)
#' y = rbinTrait(n=1,phy=mphy[[1]], beta=c(-1,0.5), alpha=.7 ,X=X)
#' cla <- rep(c("A","B","C","D","E"), each = 30)
#' dat = data.frame(y, x, cla)
#' # Run sensitivity analysis:
#' tree_clade <- tree_clade_phyglm(y ~ x, phy = mphy, data = dat, 
#' n.tree = 10, n.sim = 10, clade.col = "cla")
#'# To check summary results and most influential clades:
#'summary(tree_clade)
#'# Visual diagnostics for clade removal:
#'sensi_plot(tree_clade)
#'# Specify which clade removal to plot:
#'sensi_plot(tree_clade, "B")
#'sensi_plot(tree_clade, "C", graphs = 2)
#'sensi_plot(tree_clade, "D", graphs = 2) 
#'}
#'\dontshow{
#'set.seed(6987)
#'mphy = rmtree(150, N = 30)
#'x = rTrait(n=1,phy=mphy[[1]])
#'X = cbind(rep(1,150),x)
#'y = rbinTrait(n=1,phy=mphy[[1]], beta=c(-1,0.5), alpha=.7 ,X=X)
#'cla <- rep(c("A","B","C","D","E"), each = 30)
#'dat = data.frame(y, x, cla)
#'# Run sensitivity analysis:
#'tree_clade <- tree_clade_phyglm(y ~ x, phy = mphy, data = dat, 
#'                                n.tree = 2, n.sim = 1, clade.col = "cla")
#'# To check summary results and most influential clades:
#'summary(tree_clade)
#'# Visual diagnostics for clade removal:
#'sensi_plot(tree_clade)
#'}
#' @export

tree_clade_phyglm <-
  function(formula,
           data,
           phy,
           clade.col,
           n.species = 5,
           n.sim = 100,
           n.tree = 2,
           btol = 50,
           track = TRUE,
           ...) {
    # Error checking:
    if (!inherits(data, "data.frame"))
      stop("data must be class 'data.frame'")
    if (missing(clade.col))
      stop("clade.col not defined. Please, define the",
           " column with clade names.")
    if (!inherits(formula, "formula"))
      stop("formula must be class 'formula'")
    if (!inherits(phy, "multiPhylo"))
      stop("phy must be class 'multiPhylo'")
    if (length(phy) < n.tree)
      stop("'times' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
    
    #Match data and phy
    data_phy <- match_dataphy(formula, data, phy, ...)
    phy <- data_phy$phy
    full.data <- data_phy$data
    if (is.na(match(clade.col, names(full.data)))) {
      stop("Names column '", clade.col, "' not found in data frame'")
    }
    
    
    # If the class of tree is multiphylo pick n=times random trees
    trees <- sample(length(phy), n.tree, replace = F)
    
    
    # Identify CLADES to use and their sample size
    wc <- table(full.data[, clade.col]) > n.species
    uc <- table(full.data[, clade.col])[wc]
    
    if (length(uc) == 0)
      stop(
        paste(
          "There is no clade with more than ",
          n.species,
          " species. Change 'n.species' to fix this
                                  problem",
          sep = ""
        )
      )
    
    #List to store information
    tree.clade <- list ()
    
    #Start tree loop here
    errors <- NULL
    if (track == TRUE)
      pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
    counter = 1
    
    for (j in trees) {
      #Match data order to tip order
      full.data <- full.data[phy[[j]]$tip.label, ]
      
      #Select tree
      tree <- phy[[j]]
      
      tree.clade[[counter]] <-
        clade_phyglm(
          formula,
          data = full.data,
          phy = tree,
          btol,
          track = FALSE,
          clade.col,
          n.species,
          n.sim,
          verbose = FALSE,
          ...
        )
      
      if (track == TRUE)
        utils::setTxtProgressBar(pb, counter)
      counter = counter + 1
    }
    
    if (track == TRUE)
      close(pb)
    names(tree.clade) <- trees
    
    # Merge lists into data.frames between iterations:
    full.estimates  <-
      suppressWarnings(recombine(tree.clade, slot1 = 3, slot2 = 1))
    clade.estimates <- recombine(tree.clade, slot1 = 4)
    clade.estimates$info <- NULL
    null.dist       <- recombine(tree.clade, slot1 = 5)
    null.dist$info <- NULL
    
    #Generate output:
    res <- list(
      call = match.call(),
      formula = formula,
      full.model.estimates = full.estimates,
      sensi.estimates = clade.estimates,
      null.dist = null.dist,
      data = full.data,
      errors = errors,
      clade.col = clade.col
    )
    
    class(res) <- c("sensiTree_Clade", "sensiTree_CladeL")
    
    ### Warnings:
    if (length(res$errors) > 0) {
      warning("Some clades deletion presented errors, please check: output$errors")
    }
    else {
      res$errors <- "No errors found."
    }
    return(res)
  }

