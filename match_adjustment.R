#' Matching-Adjustment
#'
#' \code{match_adjustment} is used to estimate the weights needed to make the necessary adjustment for MAIC.
#'
#' Specifying the covariates to match on begins with creating a "matching set"
#' which is made up of the individual level patient data and the summary
#' information of that covariate. This requires the type of the summary
#' statistic from the aggregate data to which the covariate needs to match to.
#' @export
match_adjustment <- function(matching_object, outcomes = NULL, arm_identifier = NULL)
{
  centred <- centre_covariates(matching_object)
  beta <- find_beta(centred)
  weights <- find_weights(centred, beta)
  rescaled_weights <- rescale_weights(weights, nrow(centred)) # Need to be careful using nrow - might perform test for verification of sample size?
  ess <- approximate_ess(weights)
  essprop <- ess/length(weights)
  rescaled_dist <- plot_weights(rescaled_weights)

  match_adjustment_object <- list()
  match_adjustment_object$centred <- centred
  match_adjustment_object$beta <- beta
  match_adjustment_object$ess <- ess
  match_adjustment_object$essprop <- essprop
  match_adjustment_object$weights <- weights
  match_adjustment_object$rescaled <- rescaled_weights
  match_adjustment_object$rescaled_weights_hist <- rescaled_dist
  match_adjustment_object$matching_set <- matching_object$matching_set
  match_adjustment_object$targets <- matching_object$targets
  match_adjustment_object$match_adjusted_covariates <- weight_data(match_adjustment_object$matching_set, match_adjustment_object$rescaled)

  match_adjustment_object$verification <- matched_check(match_adjustment_object$matching_set, match_adjustment_object$match_adjusted_covariates, match_adjustment_object$targets, match_adjustment_object$weights)

  match_adjustment_object$verification_plot <- plot_verify(match_adjustment_object$verification)

  if (!is.null(outcomes) & is.null(arm_identifier))
  {
    stop("Require arm_identifier argument if entering outcome information.")
  }
  if (is.null(outcomes & is.null(arm_identifier)))
  {
    stop("Require outcome argument if entering arm_identifier information.")
  }
  if (!is.null(outcomes) & !is.null(arm_identifier))
  {
    if(length(outcomes) != nrow(match_adjustment_object$matching_set))
    {
      stop("Number of outcomes not equal to number of patients.")
    }
    if(length(outcomes) != length(arm_identifier))
    {
      stop("Number of outcomes different to treatment assignment information entered.")
    }

    oalist <- split(outcomes, arm_identifier)
    if(length(oalist) != 2)
    {
      warning("Patient data has not been assigned to 2 arms.")
    }

    group1 <- oalist[1]
    group2 <- oalist[2]

    match_adjustment_object$results <- list()
    match_adjustment_object$results$average_outcome_group1 <- lapply(oalist[1], mean)
    match_adjustment_object$results$average_outcome_group2 <- lapply(oalist[2], mean)
    match_adjustment_object$results$match_adjusted_means <- colMeans(match_adjustment_object$match_adjusted_covariates)

    new_oalist <- split(weight_data(outcomes, match_adjustment_object$rescaled), arm_identifier)
    # match_adjustment_object$results$new_outcomes_group1 <- new_oalist[1]
    # match_adjustment_object$results$new_outcomes_group2 <- new_oalist[2]
    match_adjustment_object$results$new_average_outcome_group1 <- lapply(new_oalist[1], mean)
    match_adjustment_object$results$new_average_outcome_group2 <- lapply(new_oalist[2], mean)
  }

  class(match_adjustment_object) <- "maic.match_adjustment"
  return(match_adjustment_object)
}

#' @export
print.maic.match_adjustment <- function(x, ...)
{
  cat("Weighting IPD for MAIC resulted in the following:\n")
  cat("\nESS:", x$ess, "\t(ESS/N):", x$essprop, "\n")
  cat("\nSummary of weights generated, scaled to original weight:\n(1 = no change, <1 = downweighted, >1 = upweighted)\n\n")
  print(summary(x$rescaled))

  if(!is.null(x$results))
  {
    cat("\nOriginal average outcomes by arm prior to weighting:\n\n")
    print(x$results$average_outcome_group1)
    print(x$results$average_outcome_group2)
    cat("New average outcomes by arm after weighting:\n\n")
    print(x$results$new_average_outcome_group1)
    print(x$results$new_average_outcome_group2)
    cat("See $results for:\n")
    print(ls(x$results))
  }

  cat("\nThe rescaled weights can be accessed by calling $rescaled from the match_adjustment x.\nSee summary() for more information. See plot() for a histogram of the rescaled weights.\n\n")
  print(ls(x))
  cat("\nFor verification of successful matching, see $verification and $verification_plot.\n")
}

#' @export
plot.maic.match_adjustment <- function(x, ...)
{
  plot(plot_weights(object$rescaled))
  # plot_verify(object$verification)
}

#' @export
summary.maic.object <- function(object, ...)
{
  if (!is.null(object$results))
  {
    cat("Based on the IPD outcomes entered:\n")
    cat("\nOriginal average outcomes by arm prior to weighting:\n")
    print(object$results$average_outcome_group1)
    print(object$results$average_outcome_group2)
    cat("New average outcomes by arm after weighting:\n")
    print(object$results$new_average_outcome_group1)
    print(object$results$new_average_outcome_group2)
    cat("Access these results and the means of the weighted matching covariates via $results and:\n")
    print(ls(object$results))
  }

  cat("\nWeighting IPD for MAIC resulted in the following:\n")
  cat("\nCovariates used in matching (centred):\n\n")
  print(head(object$centred))
  cat("(Showing up to the first 6 patients)\n")
  cat("\nSummary of weights generated, scaled to original weight:\n(1 = no change, <1 = downweighted, >1 = upweighted)\n\n")
  print(summary(object$rescaled))
  cat("\nThe Effective Sample Size (ESS) is:\n")
  cat(object$ess, "\n")
  cat("\nThe ESS as a proportion of the original sample size is:\n")
  cat(object$essprop, "\n")
  cat("\nVerification of newly weighted IPD summary statistics matching the target AgD:\n\n")
  print(object$verification)
}
