#' Matthew Phillips, s2190456
#' Add your own function definitions on this file.

#' neg_log_lik
#
#' @description Evaluate the negated log-likelihood for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model

neg_log_lik <- function(beta, data, model){
  
  mu <- beta[1] + beta[2]*data[["CAD_Weight"]]
  
  # distinguish between the two models to find the particular standard deviation for the betas
  if(model == "A") {
    sigma <- sqrt(exp(beta[3] + beta[4]*data[["CAD_Weight"]]))
  }else{
    sigma <- sqrt(exp(beta[3])+exp(beta[4]) * (data[["CAD_Weight"]]^2))
  }
  - sum(dnorm(data[["Actual_Weight"]],
              mean = mu,
              sd=sigma,
              log = TRUE))
  
}

#' filament_estimate
#
#' @description Estimate filament models with different variance structure
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @return An estimation object suitable for use with [filament1_predict()]

filament1_estimate <- function(data, model) {
  model <- match.arg(model, c("A", "B"))
  if (model == "A") {
    beta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    beta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(beta_start,
               neg_log_lik,
               data = data,
               model = model,
               hessian = TRUE,
               method = "Nelder-Mead",
               control = list(maxit = 5000)
  )
  fit <- list(
    model = model,
    par = opt$par,
    hessian = opt$hessian
  )
  class(fit) <- c("filament1_estimate", "list")
  fit
}

#' filament1_aux_EV
#' 
#' @description Evaluate the expectation and variance for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` containing the required predictors, including `CAD_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param Sigma_beta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated betas
#' @return A list with four elements:
#     E : E(y|beta,x)
#     V : Var(y|beta,x)
#     VE : Var(E(y|beta,x)|x) or NULL
#     EV : E(Var(y|beta,x)|x) or NULL

filament1_aux_EV <- function(beta, data, model = c("A", "B"),
                             Sigma_beta = NULL) {
  
  model <- match.arg(model)
  if (model == "A") {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      EV <- exp(ZV %*% beta + rowSums(ZV * (ZV %*% Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = exp(ZV %*% beta),
      VE = VE,
      EV = EV
    )
  } else {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + I(CAD_Weight^2), data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      # (pmin: Ignore large Sigma_beta values)
      EV <- ZV %*% exp(beta + pmin(0.5^2, diag(Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = ZV %*% exp(beta),
      VE = VE,
      EV = EV
    )
  }
  out
}

#' filament1_predict
#' 
#' @description Compute predictive distributions and prediction intervals for model A and B
#' @param new_data A `data.frame` used to compute the predictive mean and sd
#' @param model_estimate A `list` outputted from `filament1_estimate` containing model 
#' parameters and the hessian
#' @param alpha The level of confidence required in the prediction intervals, if not
#' specified, is set to 0.05
#' @return A `data.frame` with four columns:
#     mean : predictive mean
#     sd : predictive mean
#     uppr : upper value of prediction interval
#     lwr : upper value of prediction interval

filament1_predict <- function(new_data, model_estimate, alpha=0.05) {
  
  # Extract model information
  model <- model_estimate$model
  beta <- model_estimate$par
  sig_beta <- solve(model_estimate$hessian)
  
  aux <- filament1_aux_EV(beta, new_data, model, sig_beta)
  
  mean <- aux$E
  sd <- sqrt(aux$EV + aux$VE)
  
  z_alpha = qnorm(1 - alpha/2)
  
  n <- nrow(new_data)

  lwr <- mean - z_alpha * sd
  uppr <- mean + z_alpha * sd

  predictions_df <- data.frame(mean = mean, 
                               sd = sd,
                               lwr = lwr,
                               uppr = uppr)
  return(predictions_df)
}

#' squared_error
#' 
#' @description Compute squared error score 
#' @param predicted A single predicted mean value or `list` of predicted mean values 
#' @param actual A single observed value or `list` of observed values
#' @return The single squared error score or `list` of squared error scores

squared_error <- function(predicted, actual) {
  (predicted - actual)^2
}

#' squared_error
#' 
#' @description Compute Dawid-Sebastiani score 
#' @param predicted_mean A single predicted mean value or `list` of predicted mean values 
#' @param predicted_sd A single predicted standard deviation value or `list` of predicted 
#' standard deviation values 
#' @param actual_mean A single observed mean value or `list` of observed mean values
#' @return The single Dawid-Sebastiani score or `list` of Dawid-Sebastiani scores

dawid_sebastiani <- function(predicted_mean, predicted_sd, actual_mean) {
  (predicted_mean - actual_mean)^2 / predicted_sd^2 + 2 * log(predicted_sd)
}

#' leave1out
#' 
#' @description Compute squared error and Dawid-Sebastiani scores for a given model
#' and dataset
#' @param data A `data.frame` from which to calculate the required scores 
#' @param model A `str` defining the model for which to calculate the score
#' @return The original `data.frame` with columns added for the computed predictive 
#' mean and standard deviation and the ES and DS scores

leave1out <- function(data, model) {
  n <- nrow(data)
  
  data <- data  %>% 
    mutate(mean = NA_real_, sd = NA_real_)
  
  for (i in 1:n) {
    # Exclude ith observation
    train_data <- data[-i, , drop=FALSE]
    
    # Estimate model parameters using train_data
    model_fit <- filament1_estimate(train_data, model)
    
    # Predict for the excluded observation
    prediction <- filament1_predict(data[i, , drop=FALSE], model_fit)
    
    data[i, "mean"] <- prediction$mean
    data[i, "sd"] <- prediction$sd
  }
  
  data <- data %>%
    mutate(se = squared_error(mean, data$Actual_Weight),
           ds = dawid_sebastiani(mean, sd, data$Actual_Weight))

  return(data)
  
}

#' arch_loglike
#' 
#' @description Compute combined log-likelihood for a collection of observations and
#' and parameters
#' @param y A `list` of observations or data frame containing columns for y, N and phi
#' @param N Model parameter defining the number of buried individuals specified in model
#' @param phi Model parameter defining the probability of finding a femur, left or right, 
#' specified in model
#' @return If y is a `data.frame` return the original `data.frame` with column added 
#' to contain the combined log-likelihood for each row, else return the combined log-
#' likelihood for the given observations and parameters   

arch_loglike <- function(y, N = NULL, phi = NULL) {
  
  # Function to compute combined log-likelihood
  combined_log_likelihood <- function(y, N, phi) {
    
    log_gamma_y <- sum(lgamma(y + 1))
    log_gamma_Ny <- sum(lgamma(N - y + 1))
    
    log_likelihood <- -log_gamma_y - log_gamma_Ny  +
      2 * lgamma(N + 1) + sum(y) * log(phi) + (2 * N - sum(y)) * log(1 - phi)
    
    return(log_likelihood)
  }
  
  if (is.data.frame(y)) {
    if (!all(c("N", "phi") %in% colnames(y))) {
      stop("Data frame must contain columns named 'N' and 'phi'.")
    }
    log_likelihoods <- apply(y, 1, function(row) {
      combined_log_likelihood(row["y"], row["N"], row["phi"])
    })
    return(log_likelihoods)
  } else {
    if (is.null(N) || is.null(phi)) {
      stop("If 'y' is not a data frame, 'N' and 'phi' must be provided.")
    }
    return(combined_log_likelihood(y, N, phi))
  }
}

#' estimate
#' 
#' @description Compute Monte Carlo estimates for p_y, E_Ny and E_phi_y
#' @param y A `list` of observations 
#' @param xi Model parameter used to generate values for N
#' @param a,b Model parameters used to generate values for phi
#' @param K Number of values for N and phi to generate for Monte Carlo integration
#' @return A `list` of 3 values:
#     p_y : probability of observing given observations under provided model 
#     E_Ny : expected value of N under the provided model and given observations
#     E_phi_y : expected value of phi under the provided model and given 
#               observations

estimate <- function(y, xi, a, b, K) {
  
  N_samples <- rgeom(K, xi)
  phi_samples <- rbeta(K, a, b)
  
  # Compute likelihood for each sample
  likelihoods <- (sapply(1:K, function(k) exp(arch_loglike(y, N_samples[k], phi_samples[k]))))
  
  # Compute Monte Carlo estimates
  p_y <- mean(likelihoods)
  E_Ny <- mean(N_samples * likelihoods) / p_y
  E_phi_y <- mean(phi_samples * likelihoods) / p_y
  
  return(list(p_y = p_y, E_Ny = E_Ny, E_phi_y = E_phi_y))
}
























