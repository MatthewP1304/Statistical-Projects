#' Matthew Phillips, s2190456
#' Add your own function definitions on this file.

#' Log-Exponential density
#'
#' Compute the density or log-density for a Log-Exponential (LogExp)
#' distribution
#'
#' @param x vector of quantiles
#' @param rate vector of rates
#' @param log logical; if TRUE, the log-density is returned

dlogexp <- function(x, rate = 1, log = FALSE) {
  result <- log(rate) + x - rate * exp(x)
  if (!log) {
    exp(result)
  }
  result
}

#' Log-Sum-Exp
#'
#' Convenience function for computing log(sum(exp(x))) in a
#' numerically stable manner
#'
#' @param x numerical vector

log_sum_exp <- function(x) {
  max_x <- max(x, na.rm = TRUE)
  max_x + log(sum(exp(x - max_x)))
}


#' wquantile 
#'
#' Calculates empirical sample quantiles with optional weights, for given probabilities. 
#' Like in quantile(), the smallest observation corresponds to a probability of 0 and the largest to a probability of 1. 
#' Interpolation between discrete values is done when type=7, as in quantile(). 
#' Use type=1 to only generate quantile values from the raw input samples.
#'
#' @param x numeric vector whose sample quantiles are wanted
#' NA and NaN values are not allowed in numeric vectors unless na.rm is TRUE
#' @param probs numeric vector of probabilities with values in [0,1]
#' @param na.rm logical; if true, any NA and NaN's are removed from x before the quantiles are computed
#' @param type numeric, 1 for no interpolation, or 7, for interpolated quantiles. Default is 7
#' @param weights	 numeric vector of non-negative weights, the same length as x, or NULL. The weights are normalised to sum to 1. If NULL, then wquantile(x) behaves the same as quantile(x), with equal weight for each sample value

wquantile <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, type = 7, 
                       weights = NULL, ...) 
{
  if (is.null(weights) || (length(weights) == 1)) {
    weights <- rep(1, length(x))
  }
  stopifnot(all(weights >= 0))
  stopifnot(length(weights) == length(x))
  if (length(x) == 1) {
    return(rep(x, length(probs)))
  }
  n <- length(x)
  q <- numeric(length(probs))
  reorder <- order(x)
  weights <- weights[reorder]
  x <- x[reorder]
  wecdf <- pmin(1, cumsum(weights)/sum(weights))
  if (type == 1) {
  }
  else {
    weights2 <- (weights[-n] + weights[-1])/2
    wecdf2 <- pmin(1, cumsum(weights2)/sum(weights2))
  }
  for (pr_idx in seq_along(probs)) {
    pr <- probs[pr_idx]
    if (pr <= 0) {
      q[pr_idx] <- x[1]
    }
    else if (pr >= 1) {
      q[pr_idx] <- x[n]
    }
    else {
      if (type == 1) {
        j <- 1 + pmax(0, pmin(n - 1, sum(wecdf <= pr)))
        q[pr_idx] <- x[j]
      }
      else {
        j <- 1 + pmax(0, pmin(n - 2, sum(wecdf2 <= pr)))
        g <- (pr - c(0, wecdf2)[j])/(wecdf2[j] - c(0, 
                                                   wecdf2)[j])
        q[pr_idx] <- (1 - g) * x[j] + g * x[j + 1]
      }
    }
  }
  q
}

#' Compute empirical weighted cumulative distribution
#'
#' Version of `ggplot2::stat_ecdf` that adds a `weights` property for each
#' observation, to produce an empirical weighted cumulative distribution function.
#' The empirical cumulative distribution function (ECDF) provides an alternative
#' visualisation of distribution. Compared to other visualisations that rely on
#' density (like [geom_histogram()]), the ECDF doesn't require any
#' tuning parameters and handles both continuous and discrete variables.
#' The downside is that it requires more training to accurately interpret,
#' and the underlying visual tasks are somewhat more challenging.
#'
# @inheritParams layer
# @inheritParams geom_point
#' @param na.rm If `FALSE` (the default), removes missing values with
#'    a warning.  If `TRUE` silently removes missing values.
#' @param n if NULL, do not interpolate. If not NULL, this is the number
#'   of points to interpolate with.
#' @param pad If `TRUE`, pad the ecdf with additional points (-Inf, 0)
#'   and (Inf, 1)
#' @section Computed variables:
#' \describe{
#'   \item{x}{x in data}
#'   \item{y}{cumulative density corresponding x}
#' }
#' @seealso wquantile
#' @export
#' @examples
#' library(ggplot2)
#'
#' n <- 100
#' df <- data.frame(
#'   x = c(rnorm(n, 0, 10), rnorm(n, 0, 10)),
#'   g = gl(2, n),
#'   w = c(rep(1/n, n), sort(runif(n))^sqrt(n))
#' )
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step")
#'
#' # Don't go to positive/negative infinity
#' ggplot(df, aes(x, weights = w)) + stat_ewcdf(geom = "step", pad = FALSE)
#'
#' # Multiple ECDFs
#' ggplot(df, aes(x, colour = g, weights = w)) + stat_ewcdf()
#' ggplot(df, aes(x, colour = g, weights = w)) +
#'   stat_ewcdf() +
#'   facet_wrap(vars(g), ncol = 1)

stat_ewcdf <- function(mapping = NULL, data = NULL,
                       geom = "step", position = "identity",
                       ...,
                       n = NULL,
                       pad = TRUE,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatEwcdf,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      n = n,
      pad = pad,
      na.rm = na.rm,
      ...
    )
  )
}


#' @title StatEwcdf ggproto object
#' @name StatEwcdf
#' @rdname StatEwcdf
#' @aliases StatEwcdf
#' @format NULL
#' @usage NULL
#' @export
#' @importFrom ggplot2 aes after_stat has_flipped_aes Stat
NULL

StatEwcdf <- ggplot2::ggproto(
  "StatEwcdf", ggplot2::Stat,
  required_aes = c("x|y", "weights"),
  dropped_aes = c("weights"),     
  
  default_aes = ggplot2::aes(y = ggplot2::after_stat(y)),
  
  setup_params = function(data, params) {
    params$flipped_aes <-
      ggplot2::has_flipped_aes(data,
                               params,
                               main_is_orthogonal = FALSE,
                               main_is_continuous = TRUE)
    
    has_x <- !(is.null(data$x) && is.null(params$x))
    has_y <- !(is.null(data$y) && is.null(params$y))
    if (!has_x && !has_y) {
      rlang::abort("stat_ewcdf() requires an x or y aesthetic.")
    }
    has_weights <- !(is.null(data$weights) && is.null(params$weights))
    #    if (!has_weights) {
    #      rlang::abort("stat_ewcdf() requires a weights aesthetic.")
    #    }
    
    params
  },
  
  compute_group = function(data, scales, n = NULL, pad = TRUE, flipped_aes = FALSE) {
    data <- flip_data(data, flipped_aes)
    # If n is NULL, use raw values; otherwise interpolate
    if (is.null(n)) {
      x <- unique(data$x)
    } else {
      x <- seq(min(data$x), max(data$x), length.out = n)
    }
    
    if (pad) {
      x <- c(-Inf, x, Inf)
    }
    if (is.null(data$weights)) {
      data_ecdf <- ecdf(data$x)(x)
    } else {
      data_ecdf <-
        spatstat.geom::ewcdf(
          data$x,
          weights = data$weights / sum(abs(data$weights)) 
        )(x)
    }
    
    df_ecdf <- vctrs::new_data_frame(list(x = x, y = data_ecdf), n = length(x))
    df_ecdf$flipped_aes <- flipped_aes
    ggplot2::flip_data(df_ecdf, flipped_aes)
  }
)

# 2 Classical Estimation

#' Negative Log-Likelihood Function
#'
#' This function computes the negative log-likelihood for a given set of parameters and data, based on specified models 'A' or 'B'.
#'
#' @param beta Numeric vector of parameters.
#' @param data Data frame containing the columns 'CAD_Weight' and 'Actual_Weight'.
#' @param model Character specifying the model to use ('A' or 'B').
#'
#' @return The negative log-likelihood value for the given parameters and model.

neg_log_like <- function(beta, data, model) {

  x = data$CAD_Weight
  y = data$Actual_Weight

  if (model == 'A') {
    mu <- beta[1] + beta[2] * x
    var <- exp(beta[3] + beta[4] * x)
    neg_log_lik <- -sum(dnorm(y, mean = mu, sd = sqrt(var), log = TRUE))
  }
  else if (model == 'B') {
    mu <- beta[1] + beta[2] * x
    var <- exp(beta[3]) + (exp(beta[4]) * x^2)
    neg_log_lik <- -sum(dnorm(y, mean = mu, sd = sqrt(var), log = TRUE))
  }
  else {
    stop("Invalid model. Use 'A' or 'B'.")
  }
  return(neg_log_lik)
}

#' Parameter Estimation for Filament Model
#'
#' This function estimates the parameters for a filament model based on the specified model ('A' or 'B').
#'
#' @param data Data frame containing the columns 'CAD_Weight' and 'Actual_Weight'.
#' @param model Character specifying the model to use ('A' or 'B').
#'
#' @return A list containing the estimated parameters and Hessian matrix.
#' 
#' @seealso neg_log_like

filament1_estimate <- function(data, model) {
  if (model == 'A') {
    beta_init = c(-0.1, 1.07, -2, 0.05)
  }
  if (model == 'B') {
    beta_init = c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(
    par = beta_init,
    fn = neg_log_like,
    data = data,
    model = model,
    method = "Nelder-Mead",
    hessian = TRUE
  )

  return(list(par = opt$par, hessian = opt$hessian))
}


#' Confidence Interval Calculation
#'
#' This function calculates confidence intervals for estimated parameters based on their standard errors and a specified significance level.
#'
#' @param estimates A list containing estimated parameters and the Hessian matrix.
#' @param alpha Numeric value specifying the significance level (default is 0.1).
#'
#' @return A matrix containing the lower and upper bounds of the confidence intervals for each parameter.

CI <- function(estimates, alpha = 0.1) {
  se <- sqrt(diag(solve(estimates$hessian)))
  beta <- estimates$par
  lower <- estimates$par - se * qnorm(1-alpha/2)
  upper <- estimates$par + se * qnorm(1-alpha/2)
  intervals <- cbind(lower, upper)
}



# 3 Bayesian Estimation

# 3.1

#' Log Prior Density Calculation
#'
#' This function calculates the log prior density for a set of parameters (theta) based on specified hyperparameters.
#'
#' @param theta Numeric vector of parameters.
#' @param params Numeric vector of hyperparameters.
#'
#' @return The log prior density value for the given parameters and hyperparameters.
#'
#' @examples
#' theta <- c(0.5, 1.2, 3, 0.8)
#' params <- c(0.2, 0.5, 1, 0.3)
#' log_prior_density(theta, params)
#'
#' @seealso dlogexp

log_prior_density <- function(theta, params) {

  theta1 <- theta[1]
  theta2 <- theta[2]
  theta3 <- theta[3]
  theta4 <- theta[4]

  gamma1 <- params[1]
  gamma2 <- params[2]
  gamma3 <- params[3]
  gamma4 <- params[4]

  log_prior_theta1 <- dnorm(theta1, mean = 0, sd = gamma1, log = TRUE)

  log_prior_theta2 <- dnorm(theta2, mean = 1, sd = gamma2, log = TRUE)

  log_prior_theta3 <- dlogexp(theta3, rate = gamma3, log = TRUE)

  log_prior_theta4 <- dlogexp(theta4, rate = gamma4, log = TRUE)

  log_prior <- log_prior_theta1 + log_prior_theta2 + log_prior_theta3 + log_prior_theta4

  return(log_prior)
}

# 3.2

#' Log Likelihood Calculation
#'
#' This function calculates the log likelihood for a set of parameters (theta) given observed data (x, y).
#'
#' @param theta Numeric vector of parameters.
#' @param x Numeric vector of predictor variable values.
#' @param y Numeric vector of response variable values.
#'
#' @return The total log likelihood for the given parameters and data.
#'
#' @examples
#' theta <- c(0.5, 1.2, 3, 0.8)
#' x <- c(1, 2, 3)
#' y <- c(2, 4, 6)
#' log_like(theta, x, y)

log_like <- function(theta, x, y) {

  beta1 <- theta[1]
  beta2 <- theta[2]
  beta3 <- exp(theta[3])
  beta4 <- exp(theta[4])

  mean <- beta1 + beta2 * x
  var <- beta3 + beta4 * x^2

  log_likelihood <- dnorm(y, mean = mean, sd = sqrt(var), log = TRUE)

  total_log_likelihood <- sum(log_likelihood)

  return(total_log_likelihood)
}

# 3.3

#' Log Posterior Density Calculation
#'
#' This function calculates the log posterior density for a set of parameters (theta) given observed data (x, y) and hyperparameters.
#'
#' @param theta Numeric vector of parameters.
#' @param x Numeric vector of predictor variable values.
#' @param y Numeric vector of response variable values.
#' @param params Numeric vector of hyperparameters.
#'
#' @return The log posterior density value for the given parameters, data, and hyperparameters.
#' 
#' @examples
#' theta <- c(0.5, 1.2, 3, 0.8)
#' x <- c(1, 2, 3)
#' y <- c(2, 4, 6)
#' params <- c(0.2, 0.5, 1, 0.3)
#' log_posterior_density(theta, x, y, params)
#'
#' @seealso log_like, log_prior_density

log_posterior_density <- function(theta, x, y, params) {

  log_likelihood <- log_like(theta, x, y)

  log_prior <- log_prior_density(theta, params)

  log_posterior <- log_likelihood + log_prior

  return(log_posterior)
}

# 3.4

#' Posterior Mode Estimation
#'
#' This function estimates the posterior mode, Hessian matrix, and the inverse of the negated Hessian matrix (S) for a set of parameters (theta) given observed data (x, y) and hyperparameters.
#'
#' @param theta_start Numeric vector of initial values for the parameters.
#' @param x Numeric vector of predictor variable values.
#' @param y Numeric vector of response variable values.
#' @param params Numeric vector of hyperparameters.
#'
#' @return A list containing the posterior mode, Hessian matrix, and the inverse of the negated Hessian matrix (S).
#' 
#' @examples
#' theta_start <- c(0, 1, 2, 0.5)
#' x <- c(1, 2, 3)
#' y <- c(2, 4, 6)
#' params <- c(0.2, 0.5, 1, 0.3)
#' posterior_mode(theta_start, x, y, params)
#'
#' @seealso log_posterior_density

posterior_mode <- function(theta_start, x, y, params) {

  optim_result <- optim(par = theta_start,
                        fn = log_posterior_density,
                        x = x,
                        y = y,
                        params = params,
                        method = 'BFGS',
                        control = list(fnscale = -1),
                        hessian = TRUE)

  mode <- optim_result$par
  hessian <- optim_result$hessian
  S <- solve(-hessian)

  result_list <- list("mode" = mode, "hessian" = hessian, "S" = S)

  return(result_list)
}

# 3.6

#' Importance Sampling
#'
#' This function performs importance sampling by sampling from a multivariate normal distribution, calculating log-importance weights, and normalizing them.
#'
#' @param N Number of samples to draw.
#' @param mu Mean vector for the multivariate normal distribution.
#' @param S Covariance matrix for the multivariate normal distribution.
#' @param x Numeric vector of predictor variable values.
#' @param y Numeric vector of response variable values.
#' @param params Numeric vector of hyperparameters.
#'
#' @return A data frame containing samples and normalized log-importance weights.
#'
#' @examples
#' N <- 1000
#' mu <- c(0, 1, 2, 0.5)
#' S <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' x <- c(1, 2, 3)
#' y <- c(2, 4, 6)
#' params <- c(0.2, 0.5, 1, 0.3)
#' do_importance(N, mu, S, x, y, params)
#'
#' @seealso log_posterior_density, log_sum_exp

do_importance <- function(N, mu, S, x, y, params) {

  samples <- rmvnorm(n = N, mean = mu, sigma = S)

  log_weights <- log_posterior_density(samples, x, y, params) - dmvnorm(samples, mean = mu, sigma = S, log=TRUE)

  log_weights_normalised <- log_weights - log_sum_exp(log_weights)

  β1 = samples[, 1]
  β2 = samples[, 2]
  β3 = exp(samples[, 3])
  β4 = exp(samples[, 4])
  
  result_df <- cbind(β1, β2, β3, β4, log_weights_normalised)

  return(data.frame(result_df))
}


# 3.7

make_CI <- function(x, weights, prob) {
  cred_interval <- wquantile(x, probs = c((1 - prob)/2, 1 - (1 - prob)/2), weights = weights)
  CI <- data.frame(Lower = cred_interval[1], Upper = cred_interval[2])
  return(CI)
}



