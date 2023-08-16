

change_varname <- function(data, old_name, new_name){
  sel <- names(data) %in% old_name
  if (sum(sel) == 0){
    stop(old_name, " not found")
  } else {
    names(data)[sel] <- new_name
  }
  data
}

# for package R2jags
lm_leftcensored <- function(data, var_x, var_y, var_threshold, var_is_over_threshold,
                            quiet = FALSE,
                            init_alpha = NULL, init_beta = 0,
                            sd_alpha = 100, sd_beta = 100){
  
  data <- as.data.frame(data)
  
  # Internally, assumes variables x, y_censored, not_censored, threshold
  data <- change_varname(data, var_x, "x")
  data <- change_varname(data, var_y, "y_censored")
  data <- change_varname(data, var_threshold, "threshold")
  data <- change_varname(data, var_is_over_threshold, "not_censored")
  
  # Default alpha (= intercept) is the median of the observed values  
  if (is.null(init_alpha)){
    init_alpha <- median(data$y_censored, na.rm = TRUE)
  }
  
  # Set all cencored data to NA (if not already done)
  # Important! Otherwise all LOD stuff is ignored
  data$y_censored[data$not_censored == 0] <- NA
  
  # Jags code to fit the model to the simulated data
  # Jags code to fit the model to the simulated data
  
  model_code <- glue::glue('
model
{
  # Likelihood
  for (i in 1:n) {
    not_censored[i] ~ dinterval(y_censored[i], threshold[i])
    y_censored[i] ~ dnorm(alpha + beta * x[i], sigma^-2)
  }

  # Priors - very diffuse  
  alpha ~ dnorm({{init_alpha}}, {{sd_alpha}}^-2)
  beta ~ dnorm({{init_beta}}, {{sd_beta}}^-2)
  sigma ~ dunif(0, 10)
}
', .open = "{{", .close = "}}")  # must use double braces to avoid confusion with code

  ### Set up data and parameters
  # Set up the data
  model_data = list(n = nrow(data), 
                    y_censored = data$y_censored, 
                    not_censored = data$not_censored,
                    threshold = data$threshold,
                    x = data$x)
  # Choose the parameters to watch
  model_parameters =  c("alpha", "beta", "sigma")
  
  ### Run model
  # Run the model
  model_run = jags(data = model_data,
                   parameters.to.save = model_parameters,
                   model.file = textConnection(model_code),
                   n.chains = 4, # Number of different starting positions
                   n.iter = 5000, # Number of iterations
                   n.burnin = 1000, # Number of iterations to remove at start
                   n.thin = 2,
                   quiet = quiet) # Amount of thinning
  
  # model_run
  model_mcmc <- as.mcmc(model_run)
  # summary(model_mcmc) %>% str()
  summary(model_mcmc)
  
}


lm_leftcensored_estimates <- function(modelsummary){
  result <- cbind(modelsummary$statistics[1:2,1:2], 
                  modelsummary$quantiles[1:2, c("2.5%", "50%", "97.5%")]) %>%
    as.data.frame()
  # Change column names
  names(result)[names(result) == "50%"] <- "Median" 
  names(result)[names(result) == "2.5%"] <- "CI_lower" 
  names(result)[names(result) == "97.5%"] <- "CI_upper" 
  # Instead of rownames we want a variable 
  result$term <- rownames(result)
  rownames(result) <- NULL
  # Order columns
  result[c("term", "Mean", "Median", "SD", "CI_lower", "CI_upper")]
}

if (FALSE){
  
  sim <- sim_data()

  res <- lm_leftcensored(data = sim$data, 
                         var_x = "x", var_y = "y_cens", 
                         var_threshold = "y_LOD", var_is_over_threshold = "y_aboveLOD")
  
  # Test with other variable names
  # names(sim$data) <- c("Year", "Conc", "LOQvalue", "Flag")
  # res <- lm_leftcensored(data = sim$data, 
  #                        var_x = "Year", var_y = "Conc", 
  #                        var_threshold = "LOQvalue", var_is_over_threshold = "Flag")
  
  plot_data(sim$data)
  # True slope of simulation
  abline(sim$alpha, sim$beta,
         col = "red3",
         lwd = 2)
  # Regression based on LOD values
  naive_model <- lm(y_cens~x, sim$data)
  abline(coef(naive_model),
         col = "blue2",
         lwd = 2, 
         lty = "dashed")
  # Regression based on JAGS
  abline(res$quantiles["alpha", "50%"],
         res$quantiles["beta", "50%"],
         col = "red3",
         lwd = 2, 
         lty = "dashed")
  
  
  
}

#
# Function for simulating data from the above model
#
sim_data <- function(n = 100,
                     alpha = 30,
                     beta = -3,
                     sigma = 4,
                     c1 = 20, c2 = 5,  # Censoring limits
                     c_change = 4){
  x <- sort(runif(n, 0, 10)) # Sort as it makes the plotted lines neater
  y <- rnorm(n, mean = alpha + beta * x, sd = sigma)
  # Censoring
  y_LOD <- ifelse(x <= c_change, c1, c2)
  y_aboveLOD <- ifelse(y > y_LOD, 1, 0)
  y_cens <- pmax(y_LOD, y)
  # Alt. 1. We only have the LOD value if we are under LOD
  # y_LOD <- ifelse(y > y_LOD, NA, y_LOD)
  # Alt. 2. LOD = value
  # y_LOD <- ifelse(y > y_LOD, y, y_LOD)
  # Make data frame
  data <- data.frame(x, y_cens, y_LOD, y_aboveLOD)
  list(data = data, alpha = alpha, beta = beta)
}

#
# Function for plotting simulated data
# Nite: variable names hard-coded
#
plot_data <- function(df, version = 1){   # Version 1 is without LOD; Version 2 is with LOD
  if (version == 1){
    plot(df$x, df$y_cens, pch = 19, col = "blue2")
    # Concentrations below LOD
    sel <- df$y_aboveLOD == 0
    points(df$x[sel], df$y_cens[sel], pch = 20, col = "red3")
  } else if (version == 2){
    plot(df$x, df$y_cens, pch = 19, col = "blue2", cex = 1.8)
    # Concentrations below LOD
    sel <- df$y_aboveLOD == 0
    points(df$x[sel], df$y_cens[sel], pch = 20, col = "red3", cex = 1.8)
    # Add LOD as well (as small green dots)
    points(df$x, df$y_LOD, pch = 20, col = "green")
  }
}

# Test:
if (FALSE){
  sim <- sim_data()
  plot_data(sim$data)
  plot_data(sim$data, version = 2)
}




#
# In contrast to the new version (lm_leftcensored),
#    this demands that variables x, y_cens, y_aboveLOD, y_LOD exist
#
lm_leftcensored_OLD <- function(data){
  
  # Set all cencored data to NA (if not already done)
  # Imortant! Otherwise all LOD stuff is ignored
  data$y_cens[data$y_aboveLOD == 0] <- NA
  
  # Jags code to fit the model to the simulated data
  # Jags code to fit the model to the simulated data
  
  model_code = '
model
{
  # Likelihood
  for (i in 1:n) {
    y_aboveLOD[i] ~ dinterval(y_cens[i], y_LOD[i])
    y_cens[i] ~ dnorm(alpha + beta * x[i], sigma^-2)
  }

  # Priors
  alpha ~ dnorm(0, 100^-2)
  beta ~ dnorm(0, 100^-2)
  sigma ~ dunif(0, 10)
}
'
  ### Set up data and parameters
  # Set up the data
  model_data = list(n = nrow(data), 
                    y_cens = data$y_cens, 
                    y_aboveLOD = data$y_aboveLOD,
                    y_LOD = data$y_LOD,
                    x = data$x)
  # Choose the parameters to watch
  model_parameters =  c("alpha", "beta", "sigma")
  
  ### Run model
  # Run the model
  model_run = jags(data = model_data,
                   parameters.to.save = model_parameters,
                   model.file=textConnection(model_code),
                   n.chains=4, # Number of different starting positions
                   n.iter = 5000, # Number of iterations
                   n.burnin = 1000, # Number of iterations to remove at start
                   n.thin=2) # Amount of thinning
  
  # model_run
  model_mcmc <- as.mcmc(model_run)
  # summary(model_mcmc) %>% str()
  summary(model_mcmc)
  
}
