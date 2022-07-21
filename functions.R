
##' The objective function f constructed by 3/2 Matern
##' 
##' @param x
##' 
##' @return f(x)
f <- function(x){
  l <- length(x)
  Sigma <- Matern_k(x,x)
  R <- chol(Sigma)
  return( t(R)%*%c(rnorm(n=l,mean=0,sd = 1)) )
}

##' The objective function g constructed by Exp
##' 
##' @param x
##' 
##' @return g(x)
g <- function(x){
  l <- length(x)
  Sigma <- Exp_k(x,x)
  R <- chol(Sigma)
  return( t(R)%*%c(rnorm(n=l,mean=0,sd = 1))  )
}


##' Matern kernel function to calculate the covariance
##' 
##' @param X1 a vector 
##' @param X2 a vector 
##' 
##' @return the covariance matrix of X1 and X2
Matern_k <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      r <- abs(X1[i]-X2[j])
      Sigma[i,j] <- (1+sqrt(3)*r/l)*exp(-sqrt(3)*r/l)
    }
  }
  return (Sigma)
}

##' Squared Exponential kernel 
##' @param X1 a vector 
##' @param X2 a vector 
##' 
##' @return the covariance matrix of X1 and X2
SE_k <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      r <- abs(X1[i]-X2[j])
      Sigma[i,j] <- exp(-1/2*(r/l)^2)
    }
  }
  return (Sigma)
}

##' Exponential kernel 
##' @param X1 a vector 
##' @param X2 a vector 
##' 
##' @return the covariance matrix of X1 and X2
Exp_k <- function(X1,X2,l=1){
  Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow=length(X1))
  for (i in 1:nrow(Sigma)) {
    for (j in 1:ncol(Sigma)) {
      r <- abs(X1[i]-X2[j])
      Sigma[i,j] <- exp(-r/l)
    }
  }
  return (Sigma)
}

##' Criterion 1--MMSE
##' Choosing the next point where the variance is the largest
##' 
##' @param Sigma the covariance matrix
##' @param sigma2_noise variance of observation noise
##' 
##' @return the index of the next point to be chosen
criterion_MMSE <- function(Sigma,sigma2_noise=0){
  return (which.max(diag(Sigma)))
}

##' Criterion 2--IMSE
##' Choosing the next point at which the overall
##'  variance can be reduced the most
##'
##' @param Sigma the covariance matrix
##' @param sigma2_noise variance of observation noise
##' 
##' @return the index of the next point to be chosen
criterion_IMSE <- function(Sigma,sigma2_noise=0){
  int_var <- rowSums(Sigma^2)/(diag(Sigma) + sigma2_noise)
  return (which.max(int_var))
}

##' Criterion 3--Dawid-Sebastiani (IMDS = Integrated mean D-S)
##' Choosing the next point at which the overall
##'  expected Dawid-Sebastiani score is minimised
##'
##' @param Sigma the covariance matrix
##' @param sigma2_noise variance of observation noise
##'
##' @return the index of the next point to be chosen
criterion_IMDS <- function(Sigma,sigma2_noise=0){
  # Try to handle the zero noise case:
  sigma2_noise <- max(1e-6, sigma2_noise)
  Var <- diag(Sigma)
  ok <- is.finite(Var)
  ok[ok] <- Var[ok] > 0
  n <- sum(ok)
  V <- matrix(Var[ok], n, n)
  int_ds <- colSums(log(V - Sigma[ok, ok]^2 / (t(V) + sigma2_noise)))
  return (which(ok)[which.min(int_ds)])
}
criterion_0 <- function(Sigma,sigma2_noise=0){
  # Try to handle the zero noise case:
  sigma2_noise <- max(1e-6, sigma2_noise)
  Var <- diag(Sigma)
  ok <- is.finite(Var)
  ok[ok] <- Var[ok] > 0
  pred_v <- Var[ok] - colSums(Sigma[ok, ok]^2) / (Var[ok] + sigma2_noise)
  return (which(ok)[which.min(pred_v)])
}


#' Compute an observation update
#'
#' @param mu Previous expectation vector
#' @param Sigma Previous covariance matrix
#' @param sigma2_noise variance of observation noise
#' @param y_new Vector of the new observation(s)
#' @param obs logical vector or index vector (into mu and Sigma) for
#'   the new observation(s)
update_estimate <- function(mu, Sigma, sigma2_noise,
                            y, obs) {
  if (sigma2_noise == 0) {
    # Only involve elements with nonzero variance
    obs_prev <- diag(Sigma) == 0
    if (any(obs_prev[obs])) {
      warning("Element already observed, and noise variance is zero. 
              Expect numerical problems.")
    }
    obs_all <- obs_prev
    obs_all[obs] <- TRUE
    # No observation noise
    Sigma_obs <- Sigma[obs, obs, drop=FALSE]
    # Update unobserved covariances and expectations:
    Sigma_new <- Sigma
    Sigma_new[!obs_all,!obs_all] <- Sigma[!obs_all,!obs_all] -
      Sigma[!obs_all, obs, drop=FALSE] %*%
      solve(Sigma_obs, Sigma[obs, !obs_all,drop=FALSE])
    Sigma_new[obs, ] <- 0
    Sigma_new[, obs] <- 0
    # Symmetrise:
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    
    mu_new <- mu
    mu_new[obs] <- y
    mu_new[!obs_all] <- mu[!obs_all] +
      Sigma[!obs_all, obs, drop=FALSE] %*%
      solve(Sigma_obs, y - mu[obs])
  } else {
    # Find out how many new observations there are:
    if (is.logical(obs)) {
      n <- sum(obs)
    } else {
      n <- length(obs)
    }
    # Add observation noise to the observation submatrix:
    Sigma_obs <- (Sigma[obs, obs, drop=FALSE] +
                    diag(x = sigma2_noise, nrow = n, ncol = n))
    # Update unobserved covariances and expectations:
    Sigma_new <- Sigma -
      Sigma[, obs, drop=FALSE] %*%
      solve(Sigma_obs, Sigma[obs, ,drop=FALSE])
    # Symmetrise:
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    
    mu_new <- mu +
      Sigma[, obs,drop=FALSE] %*%
      solve(Sigma_obs, y - mu[obs])
  }
  list(mu = mu_new, Sigma = Sigma_new)
}




##' Fit a GP model with sequentially chosen points
##' 
##' @param x a vector 
##' @param y a vector, true value of f over the domain
##' @param initial_design vector indices, indicating the initial designs
##' @param criterion MMSE or IMSE 
##' @param kernel a function, to calculate initial covariance matrix
##' @param sigma2_noise a scalar, observation variance
##' @param max_iterations the number of total iterations 
##' 
##' @return a dataframe of 5 veriables, 
##' including x, y, prediction expectation, prediction variance,
##' and set of all the observations.
seq_DoE <- function(x,y,initial_design, criterion,
                    kernel=Matern_k, sigma2_noise=0, max_iterations=30){
  
  y_hat <- rep(NA, length(y))
  
  # Initial mean
  mu <- rep(0,length(y))
  mu_new <-mu
  
  # Initial covariance:
  Sigma <- kernel(x,x)
  Sigma_new <- Sigma
  
  # Initial observation
  obs_new <- rep(FALSE, nrow(Sigma))
  x_idx_new <- initial_design
  obs_new[x_idx_new] <- TRUE
  
  new_values <- update_estimate(mu, Sigma, sigma2_noise = sigma2_noise,
                                y = y[x_idx_new], obs = x_idx_new)
  mu <- new_values$mu
  Sigma <- new_values$Sigma
  obs <- obs_new
  
  # Update current "best estimate"
  y_hat <- mu
  
  for (i in 1:max_iterations){
    # Choose the next point
    x_idx_new <- criterion(Sigma, sigma2_noise = sigma2_noise)
    
    # Now condition on the new point, at index x_idx_new
    obs_new <- rep(FALSE, nrow(Sigma))
    obs_new[x_idx_new] <- TRUE
    
    new_values <- update_estimate(mu, Sigma, sigma2_noise = sigma2_noise,
                                  y = y[x_idx_new], obs = x_idx_new)
    mu <- new_values$mu
    Sigma <- new_values$Sigma
    
    # Update current "best estimate"
    y_hat <- mu
    
    # Update:
    obs <- obs | obs_new
    
  }
  # Calculate proper scores
  #se <- SE_score(y_pred=y_hat,y_true=y)
  #ds <- DS_score(y_pred=y_hat,y_true=y,sigma2=pmax(0,diag(Sigma)))
  
  #Return data frame
  data.frame(x=x, y=y, obs=obs, 
             y_hat=y_hat, var=pmax(0,diag(Sigma)))
}

##' Plot prediction
##' 
##' @param df dataframe of 5 variables (x, y, y_hat, var and obs)
##'    usually the output of seq_DoE
##' 
##' @return ggplot of prediction mean and confidence interval
plt_pred <- function(df){
  ggplot(df)+
    geom_ribbon(aes(x = x, ymin = y_hat-1.97*sqrt(var), ymax = y_hat+1.97*sqrt(var)),
                alpha = 0.2) +
    geom_line(aes(x = x, y = y_hat), linetype = "dashed") +
    geom_point(data = data.frame(x=x[df$obs],y=y[df$obs]),
               aes(x = x, y = y),size = 2) +
    theme_bw() +
    geom_line(aes(x=x,y=y), colour= "red")+
    ggtitle("Fitted GP model")
}


##' Plot prediction variance 
##' 
##' @param var prediction variance, diagonal of covariance matrix
##' 
##' @return ggplot of variance at each x 
plt_var <- function(var){
  ggplot(data.frame(x=x, var=var))+
    geom_line(aes(x = x, y = var))+
    ylab("variance") +
    xlab("x") +
    ggtitle("Prediction variance")
}


##' Calculate SE scores for each x
##' 
##' @param y_pred vector of prediction values (\mu(x) in the case)
##' @param y_true vector of actual values (y(x))
##' 
##' @return vector of SE scores
SE_score <- function(y_pred, y_true){
  (y_true-y_pred)^2
}

##' Calculate DS scores for each x
##' 
##' @param y_pred vector of prediction values (\mu(x) in the case)
##' @param y_true vector of actual values (y(x))
##' @param sigma2 vector of the prediction varainces
##' 
##' @return vector of DS scores
DS_score <- function(y_pred, y_true,sigma2){
  ds<- (y_true-y_pred)^2/sigma2+log(sigma2)
  ds
}













