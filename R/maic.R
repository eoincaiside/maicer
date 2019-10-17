# Core S3 functions and definition of maic object

#' Performing a Matching-Adjusted Indirect Comparison
#'
#' \code{maic} is used to perform a Matching_Adjusted Indirect Comparison. Given a set of individual patient data (IPD), a set of covariates to match on, \code{maic} will generate weights to be applied to the IPD in order to perform a Matching-Adjusted Indirect Comparison. Currently the package can match continuous covariates to a mean and standard deviation as well as binary covariates to a target proportion.
#'
#' It is important that the IPD is in a particular format.
#'
#' @section IPD Format:
#' A patient ID column is optional. The maicer package anticipates that IPD will be in a data frame format with each covariate listed in a separate column. For an anchored MAIC, a variable denoting the arm to which the patient was assigned is necessary and should follow the patient characteristic columns. The final column should give the outcome, y.
#' @export
maic <- function(ipd, matching_covariates, outcome_identifier, arm_identifier, outcome_type = NULL, binary_event_type = NULL, comparator_to_anchor = NULL, se_comparator_to_anchor = NULL, name_intervention = "Intervention", name_anchor = "Anchor", name_comparator = "Comparator")
{
  maic_object <- match_adjustment(matching_covariates)
  maic_object$outcomes <- ipd[,outcome_identifier]
  maic_object$average_outcome <- mean(maic_object$outcomes)
  maic_object$arms <- ipd[,arm_identifier]
  if (is.null(comparator_to_anchor))
  {
    warning("Require relative treatment effect between comparator and anchor in order to calculate relative treatment effect between intervention and comparator.")
  }
  maic_object$comparator_to_anchor <- comparator_to_anchor
  if (is.null(se_comparator_to_anchor))
  {
    warning("Standard error of relative effect between comparator and anchor required for standard error calculation of relative effect between intervention and comparator.")
  }
  maic_object$se_comparator_to_anchor <- se_comparator_to_anchor
  class(maic_object) <- "maic"
  weights <- maic_object$weights

  if (outcome_type == "binary")
  {
    if (is.null(binary_event_type))
    {
      stop("Warning: For a binary outcome, the outcome value as specified by the event_identifier_name must be defined as either a positive event for the patient (enter event_type = \"positive\"), or a negative event (event_type = \"negative\").")
    }
    adjusted_model <- anchored_relative_effect_binary_outcome(ipd, outcome_identifier, event_type = binary_event_type, arm_identifier, maic_object$weights)
    maic_object$population_adjusted_effect <- coef(adjusted_model)[2]
    names(maic_object$population_adjusted_effect) <- NULL
    maic_object$population_adjusted_se_estimate <- sandwich::vcovHC(adjusted_model)[4]
    maic_object$relative_adjusted_effect <- maic_object$comparator_to_anchor - maic_object$population_adjusted_effect
    maic_object$relative_adjusted_effect_se <- maic_object$population_adjusted_se_estimate + maic_object$se_comparator_to_anchor

    unadjusted_model <- anchored_relative_effect_binary_outcome(ipd, outcome_identifier, event_type = binary_event_type, arm_identifier, weights = NULL)
    maic_object$unadjusted_intervention_to_anchor <- coef(unadjusted_model)[2]
    names(maic_object$unadjusted_intervention_to_anchor) <- NULL
    maic_object$unadjusted_intervention_to_anchor_se <- sandwich::vcovHC(unadjusted_model)[4]

    maic_object$relative_unadjusted_effect <- maic_object$comparator_to_anchor - maic_object$unadjusted_intervention_to_anchor
    maic_object$relative_unadjusted_effect_se <- maic_object$se_comparator_to_anchor + maic_object$unadjusted_intervention_to_anchor_se

    maic_object$results <- data_frame(id = 1:5,
                                      Comparison = c(paste(name_intervention,"\nvs\n", name_anchor), paste(name_intervention,"\nvs\n", name_anchor), paste(name_comparator,"\nvs\n", name_anchor), paste(name_comparator,"\nvs\n", name_intervention), paste(name_comparator,"\nvs\n", name_intervention)),

                                      Estimate = c(maic_object$population_adjusted_effect, maic_object$unadjusted_intervention_to_anchor, maic_object$comparator_to_anchor, maic_object$relative_adjusted_effect, maic_object$relative_unadjusted_effect),

                                      var = c(maic_object$population_adjusted_se_estimate, maic_object$unadjusted_intervention_to_anchor_se, maic_object$se_comparator_to_anchor, maic_object$relative_adjusted_effect_se, maic_object$relative_unadjusted_effect_se),

                                      lo = Estimate + qnorm(0.025)*sqrt(var),

                                      hi = Estimate + qnorm(0.975)*sqrt(var),

                                      type = c("MAIC", "Unadjusted", "Unadjusted", "MAIC", "Unadjusted"))

    return (maic_object)
  }

  # Something to go here.
}

#' @export
print.maic <- function(x, ...)
{
  # maic_print <- list()
  # maic_print$ess <- object$ess
  # maic_print$essprop <- object$essprop
  # maic_print$comparator_to_anchor <- object$comparator_to_anchor
  # maic_print$population_adjusted_effect <- object$population_adjusted_effect
  # maic_print$population_adjusted_se_estimate <- object$population_adjusted_se_estimate
  #
  # maic_print$relative_treatment_effect <- object$relative_treatment_effect

  cat("A Matching-Adjusted Indirect Comparison resulted in the following:\n\nPopulation-Adjusted Relative Treatment effect: ", x$population_adjusted_effect)
  cat("\nStandard error:\t\t\t\t\t", x$population_adjusted_se_estimate)
  cat("\n\nThe MAIC had an Effective Sample Size (ESS) of:\t", x$ess, "\nESS as a proportion of sample size:\t\t", x$essprop)
  cat("\n\nSee summary() for more information. If comparator to anchor trial information has been entered, see plot() for a forest plot of adjusted vs unadjusted results.")
}

#' @export
plot.maic <- function(x, ...)
{

  maic_plot <- ggplot2::ggplot(aes(x = Estimate, y = id, col = type, shape = type), data = x$results) + geom_vline(xintercept = 0, lty = 2) + geom_point(size = 2) + geom_segment(aes(y = id, yend = id, x = lo, xend = hi), na.rm = TRUE) + xlab("Estimate (Log OR)") + facet_grid(Comparison~., switch = "y", scales = "free_y", space = "free_y") + scale_y_reverse(name = "", breaks = NULL, expand = c(0, 0.6))

  return(maic_plot)
}

#' @export
summary.maic <- function(object, ...)
{
  cat("A Matching-Adjusted Indirect Comparison resulted in the following:\n\nPopulation-Adjusted Relative Treatment effect: ", object$population_adjusted_effect)
  cat("\nStandard error:\t\t\t\t\t", object$population_adjusted_se_estimate)
  cat("\n\nThe MAIC had an Effective Sample Size (ESS) of:\t", object$ess, "\nESS as a proportion of sample size:\t\t", object$essprop)
  cat("\n\nThe overall results:")
  print(object$results)
  cat("\n")
  print(ls(object))
}

