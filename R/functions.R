# devtools::load_all()

## suggestion -- set defaults when arguments take specific values

# It may also be worth building a "Summary" function to output everything together e.g. head of rescaled weights, summary of weight distribution, ESS, percentage of original sample size etc. Function might be applied to "MAIC" object which consists of ipd, matching set etc.

# summary <- function(results) ??

# This function should take in the IPD
# This function takes in a data frame and returns the variables and their type; quantitative, categorical or ordinal.

# Important for the user to note that if the IPD contains binary variables, these must be coded as {0, 1} with 1 corresponding to the level which is examined in the AgD e.g. if AgD reports male proportion, then gender must be encoded {female = 0, male = 1}.
summarise_ipd <- function(data_file)
{
  variables <- colnames(data_file)
  cat("R recognised the following variables and their type as follows:\n")
  print(lapply(data_file, recognise_variable_type))
}

# This function below may be better off being placed in a separate script with the summarise_ipd() function.
recognise_variable_type <- function(variable)
{
  type <- "Quantitative"
  if(is.factor(variable))
  {
    type <- "Categorical"
  }
  if(is.ordered(variable))
  {
    type <- "Ordinal"
  }
  return(type)
}

#' Initialise the Set of Matching Covariates
#'
#' \code{set_matching} creates a data frame of one covariate from the individual
#' patient data to be used for matching purposes and to calculate the weights.
#'
#' Specifying the covariates to match on begins with creating a "matching set"
#' which is made up of the individual level patient data and the summary
#' information of that covariate. This requires the type of the summary
#' statistic from the aggregate data to which the covariate needs to match to.
#' @export
set_matching <- function(ipd, covariate_name, target_mean = NULL, target_sd = NULL, target_proportion = NULL)
{
  matching_information <- list()
  ipd <- as.data.frame(ipd)
  matching_information$matching_set <- as.data.frame(ipd[,covariate_name])
  names(matching_information$matching_set) <- covariate_name
  # matching_information$matching_set <- matching_set

  if (!is.null(target_mean))
  {
    mean_statistic <- mean(ipd[,covariate_name])
    matching_information$mean_statistic[covariate_name] <- mean_statistic
    # names(matching_information$mean_statistic) <- covariate_name
    matching_information$targets[paste(covariate_name, ".mean", sep = "")] <- target_mean
    # names(matching_information$targets) <- paste(covariate_name, ".mean", sep = "")
  }
  if (!is.null(target_sd))
  {
    matching_information$matching_set[,paste(covariate_name, "^2", sep = "")] <- (ipd[,covariate_name])^2
    sd <- sd(ipd[,covariate_name])
    matching_information$sd_statistic[covariate_name] <- sd
    # names(matching_information$sd_statistic) <- covariate_name
    matching_information$targets[paste(covariate_name, ".sd", sep = "")] <- target_sd
    # names(matching_information$targets)[2] <- paste(covariate_name, ".sd", sep = "")
  }
  if (!is.null(target_proportion))
  {
    matching_information$matching_set[,covariate_name] <- ipd[,covariate_name]
    counts <- table(ipd[,covariate_name])
    matching_information$proportion[covariate_name] <- as.numeric(prop.table(counts)[1])
    # names(matching_information$proportion) <- covariate_name
    matching_information$targets[paste(covariate_name, ".proportion", sep = "")] <- target_proportion
    # names(matching_information$targets) <- paste(covariate_name, ".proportion", sep = "")
  }
  matching_information
}

#' Add a Covariate to the Matching Set
#'
#' \code{add_matching_covariate} adds a covariate to the matching set initialised
#' by \code{set_matching}.
#'
#' When the covariate is added the patient information will be added to the
#' \code{matching_set} object, the target will be added to $targets and the
#' relevant summary statistic of the covariate prior to matching will be
#' calculated and appended to the \code{matching_set}.
#' @export
add_matching_covariate <- function(matching_information, ipd, covariate_name, target_mean = NULL, target_sd = NULL, target_proportion = NULL)
{
  ipd <- as.data.frame(ipd)
  # relevant_index <- ncol(matching_information$matching_set)+1
  matching_information$matching_set <- as.data.frame(matching_information$matching_set)
  matching_information$matching_set[,covariate_name] <- ipd[,covariate_name]
  # names(matching_information$matching_set[,covariate_name]) <- covariate_name

  # The order of these if-statements is important. Need some sort of error handling for ensuring user doesn't enter mean, sd and proportion together.
  if (!is.null(target_mean))
  {
    mean_statistic <- mean(ipd[,covariate_name])
    matching_information$mean_statistic[covariate_name] <- mean_statistic
    # names(matching_information$mean_statistic)[length(names(matching_information$mean_statistic))] <- covariate_name
    matching_information$targets[paste(covariate_name, ".mean", sep = "")] <- target_mean
    # names(matching_information$targets)[length(names(matching_information$targets))] <- paste(covariate_name, ".mean", sep = "")
  }
  if (!is.null(target_sd))
  {
    matching_information$matching_set[,paste(covariate_name, "^2", sep = "")] <- (ipd[,covariate_name])^2
    sd <- sd(ipd[,covariate_name])
    matching_information$sd_statistic[covariate_name] <- sd
    # names(matching_information$sd_statistic)[length(names(matching_information$sd_statistic))] <- covariate_name
    matching_information$targets[paste(covariate_name, ".sd", sep = "")] <- target_sd
    # names(matching_information$targets)[length(names(matching_information$targets))] <- paste(covariate_name, ".sd", sep = "")
  }
  if (!is.null(target_proportion))
  {
    counts <- table(ipd[,covariate_name])
    # Note that this [2] index means that the proportion recorded for a binary covariate will relate to the event (not the non-event)
    matching_information$proportion[covariate_name] <- as.numeric(prop.table(counts)[2])
    # names(matching_information$proportion)[length(names(matching_information$proportion))] <- covariate_name

    matching_information$targets[paste(covariate_name, ".proportion", sep = "")] <- target_proportion
    # names(matching_information$targets)[length(names(matching_information$targets))] <- paste(covariate_name, ".proportion", sep = "")
  }
  matching_information
}


# THIS function is legacy and needs updating
#' @export
remove_matching_covariate <- function(matching_information, covariate_name)
{
  if (length(matching_information$matching_set) > 0)
  {
    for(removal in 1:length(matching_information$matching_set))
    {
      if ( grepl(covariate_name, names(matching_information$matching_set)[removal] ))
      {
        matching_information$matching_set <- matching_information$matching_set[,-removal]
      }
    }
  }
  if (length(matching_information$mean_statistic) > 0)
  {
    for(removal in 1:length(matching_information$mean_statistic))
    {
      if ( grepl(covariate_name, names(matching_information$mean_statistic)[removal] ))
      {
        matching_information$mean_statistic <- matching_information$mean_statistic[-removal]
      }
    }
  }
  if (length(matching_information$sd_statistic) > 0)
  {
    for(removal in 1:length(matching_information$sd_statistic))
    {
      if ( grepl(covariate_name, names(matching_information$sd_statistic)[removal] ))
      {
        matching_information$sd_statistic <- matching_information$sd_statistic[-removal]
      }
    }
  }
  if (length(matching_information$proportion) > 0)
  {
    for(removal in 1:length(matching_information$proportion))
    {
      if ( grepl(covariate_name, names(matching_information$proportion)[removal] ))
      {
        matching_information$proportion <- matching_information$proportion[-removal]
      }
    }
  }
  if (length(matching_information$targets) > 0)
  {
    for(removal in 1:length(matching_information$targets))
    {
      if ( grepl(covariate_name, names(matching_information$targets)[removal] ))
      {
        matching_information$targets <- matching_information$targets[-removal]
      }
    }
  }

  # Note that this will not remove the column with a squared version of a covariate, this currently has to be entered as a separate remove call.
  # relevant_index <- which(colnames(matching_information$matching_set) == covariate_name)
  # matching_information$matching_set[c(relevant_index, relevant_index+1)] <- NULL
  # matching_information$mean_statistic <- matching_information$mean_statistic[-c(relevant_index, relevant_index+1)]
  # matching_information$sd_statistic <- matching_information$sd_statistic[-c(relevant_index, relevant_index+1)]
  return(matching_information)
}

# Here ipd is our patient data, matching_vector is the vector of covariates to be used for matching and targets is the AgD or summary data in our competitor trial which we need to match to.

# May need to eliminate the 'for' loop.


centre_covariates <- function(matching_object)
{
  centred_ipd <- matrix(NA, nrow = nrow(matching_object$matching_set), ncol = ncol(matching_object$matching_set), dimnames = list(rownames(matching_object$matching_set), colnames(matching_object$matching_set)))

  for (covariate in 1: ncol(matching_object$matching_set))
  {
    covariate_name <- names(matching_object$matching_set)[covariate]
    if (grepl("mean",names(matching_object$targets)[covariate]))
    {
      centred_ipd[,covariate] <- matching_object$matching_set[,covariate_name] - matching_object$targets[paste(covariate_name, ".mean", sep = "")]
    }
    if (grepl("sd",names(matching_object$targets)[covariate]))
    {
      covariate_name <- sub(pattern = "^2", replacement = "", x = covariate_name, fixed = TRUE)
      centred_ipd[,covariate] <- matching_object$matching_set[,paste(covariate_name, "^2", sep = "")] - ( matching_object$targets[paste(covariate_name, ".mean", sep = "")]^2 + matching_object$targets[paste(covariate_name, ".sd", sep = "")]^2 )
    }
    if (grepl("proportion",names(matching_object$targets)[covariate]))
    {
      centred_ipd[,covariate] <- matching_object$matching_set[,covariate_name] - matching_object$targets[paste(covariate_name, ".proportion", sep = "")]
    }
  }
  centred_ipd
}


# For the optimisation function in find_beta
objective_function <- function(beta, x)
{
  sum( exp(x %*% beta) )
}
gradient_function <- function(beta, x)
{
  colSums(sweep(x, 1, exp(x %*% beta), "*"))
}

# Using BFGS minimisation optimisation algorithm
find_beta <- function(centered_data)
{
  # For the optimising function, the data needs to be of matrix type.
  optimising <- optim(par = rep(0, ncol(centered_data)), fn = objective_function, gr = gradient_function, x = centered_data, method = "BFGS")

  # Error handling:
  if (optimising$convergence != 0)
  {
    warning("Convergence not reached. Inadequate patient population overlap.")
    # stop()
  }
  beta_parameter <- optimising$par
  return(beta_parameter)
}

find_weights <- function(centered_data, beta)
{
  exp(centered_data %*% beta)
}

rescale_weights <- function(weights, ipd_size)
{
  rescaled_weights <- (weights/ sum(weights)) * ipd_size
  return(rescaled_weights)
}

approximate_ess <- function(weights)
{
  sum(weights)^2/sum(weights^2)
}

# Provide multiple functions to plot weight distribution?
plot_weights <- function(weights)
{
  ggplot2::qplot(weights, geom="histogram", xlab = "Rescaled weight (multiple of original unit weight)", xlim = c(0, max(weights)+1), ylab = "Frequency", binwidth = 0.25) + geom_vline(xintercept = 1, lty = 2)
}

# We need to consider whether this function is needed and what the scope of this tool is. This would be a form of error-handling. Can check that centred data actually matches to the input targets.

#' Verify matching of weighted IPD to target AgD.
#'
#' \code{matched_check}
#'
#' @export
matched_check <- function(original_data, matched_data, targets, weights, ...)
{
  verification <- as.data.frame(targets)
  for (target in 1:length(targets))
  {
    recognised <- FALSE
    if (grepl("mean",names(targets[target])))
    {
      recognised <- TRUE
      covariate_name <- unlist(strsplit(names(targets)[target], ".", fixed = TRUE))[1]
      verification[target, 2] <- mean(matched_data[,covariate_name])
    }
    if (grepl("sd",names(targets[target])))
    {
      recognised <- TRUE
      covariate_name <- unlist(strsplit(names(targets)[target], ".", fixed = TRUE))[1]

      verification[target, 2] <- sqrt( sum(weights/sum(weights)*(original_data[,covariate_name] - mean(matched_data[,covariate_name]))^2 ) )

    }
    if (grepl("proportion",names(targets[target])))
    {
      recognised <- TRUE
      covariate_name <- unlist(strsplit(names(targets)[target], ".", fixed = TRUE))[1]
      verification[target, 2] <- mean(matched_data[,covariate_name])
    }
    if (!recognised)
    {
      warning("Target type not recognised! Inspect targets entered.")
    }
  }
  colnames(verification)[2] <- "weighted_data"
  divergence <- verification$weighted_data - verification$targets
  percent_divergence <- divergence/verification$targets
  verification <- cbind(summary = row.names(verification), verification, divergence, percent_divergence)
  verification$percent_divergence <- paste(round(verification$percent_divergence*100, digits = 5), "%", sep="")
  return (verification)
}

weight_data <- function(to_weight, weights)
{
  return (to_weight*weights)
}

# For plotting the distribution of the new covariates (verification exercise).
plot_covariate <- function(match_adjusted_covariate)
{
  ggplot2::ggplot(match_adjusted_covariate) + geom_histogram(aes(x = match_adjusted_covariate))
  ggplot(match_adjusted_data, aes(x = Age)) + geom_density()
}

plot_verify <- function(verification)
{
  # divergence <- verification$weighted_data - verification$targets
  # percent_divergence <- divergence/verification$targets
  # new_verification <- cbind(summary = row.names(verification), verification, divergence, percent_divergence)
  # new_verification$percent_divergence <- paste(round(new_verification$percent_divergence*100, digits = 5), "%", sep="")
  ggplot2::ggplot(verification, aes(summary, percent_divergence)) + geom_bar(stat = "identity") + xlab("Targets") + ylab("% Error from Target Value") + ggtitle("Deviation of match-adjusted data as a proportion of target value") + expand_limits(x = 0, y = 0)
}

# For binary outcomes, we can generate the population-adjusted relative treatment effect for  by using a binomial GLM.
# Note that it is important to observe in the glm output that the coefficient is labelled "as.factor(treatment_identifer)" and then TREATMENT or CONTROL, this tells us that the predicted Control takes the intercept value whereas Treatment takes: intercept + coefficient.
anchored_relative_effect_binary_outcome <- function(ipd, event_identifier_name, event_type = NULL, treatment_identifier_name, weights = NULL)
{
  if (is.null(event_type))
  {
    stop("Warning: For a binary outcome, the outcome value as specified by the event_identifier_name must be defined as either a positive event for the patient (enter event_type = \"positive\"), or a negative event (event_type = \"negative\").")
  }
  ipd <- as.data.frame(ipd)
  if (event_type == "negative")
  {
    event_identifier <- ipd[,event_identifier_name]
    non_event_identifier <- 1 - event_identifier
    non_event_identifier_name <- paste(event_identifier_name, "positive", sep="_")
    ipd[,non_event_identifier_name] <- non_event_identifier
  }
  if (event_type == "positive")
  {
    non_event_identifier_name <- event_identifier_name
    non_event_identifier <- ipd[,non_event_identifier_name]
    event_identifier <- 1 - non_event_identifier
    event_identifier_name <- paste(event_identifier_name, "negative", sep="_")
    ipd[,non_event_identifier_name] <- non_event_identifier
  }
  treatment_identifier <- ipd[,treatment_identifier_name]

  binary_model <- glm( cbind(event_identifier, non_event_identifier)~as.factor(treatment_identifier), weights = weights, family = "binomial" )
  return (binary_model)
}

# The final indirect comparison
indirect_comparison <- function(relative_effect_comparator_to_placebo, population_adjusted_intervention_to_placebo, outcome_scale = "Log Odds")
{
  if (outcome_scale == "Log Odds")
  {
    return (indirect_comparison <- relative_effect_comparator_to_placebo - population_adjusted_intervention_to_placebo)
  }
}

# We'll take the raw outcomes into the function and calculate the mean.
# unanchored_relative_effect_binary_outcome <- function(intervention_individual_outcomes, maic, comparator_summary_outcomes, comparator_n)
# {
#   intervention_n <- length(maic$weights) # Is this acceptable?
#   comparator_outcome_sum <- sum(comparator_summary_outcomes) # SUM(Y_C_(C))
#   # sum_intervention_outcomes <- sum(intervention_individual_outcomes) # SUM(Y_B_(B))
#   intervention_summary_outcomes <- mean(intervention_individual_outcomes) # Ybar_B_(B)
#   match_adjusted_intervention_sum <- weighted.mean(intervention_individual_outcomes, maic$weights) # Y_B_(C)
#   match_adjusted_intervention_fails <- intervention_n - match_adjusted_intervention_sum
#   comparator_fails <- comparator_n - comparator_outcome_sum
#
#   unanchored_relative_effect <- log( (match_adjusted_intervention_sum / comparator_outcome_sum) / (match_adjusted_intervention_fails / comparator_fails) )
#   return (unanchored_relative_effect)
# }
