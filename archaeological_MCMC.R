library(readxl)
library(mvtnorm)

# load data
all.dat <- read_xlsx("Copy of HazelnutsAllData_forStats_240105.xlsx")
nr.all <- dim(all.dat)[1]; nc.all <- dim(all.dat)[2]
nr.old <- 63 # old archaeological data
nr.new <- 192
mod.dat <- all.dat[1:192,]
old.dat <- all.dat[193:255,] 

# discretize p
S <- 83 # assign p by every 100 years (from=-7100, to=1100), there are S=83 p's in total.

y_M <- mod.dat$D13C
y_A <- old.dat$D13C
lM <- mod.dat$LAI
lower_tau <- old.dat$lower
upper_tau <- old.dat$upper
mu_mean <- mean(lM[lM>0])
b <- sqrt(var(lM[lM>0]))

# MCMC sampling
## prior 

### tilde_p ~ AR(1)
logprior_tildep <- function(tilde_p, theta, v){
  logdensity <- dnorm(tilde_p[1], 0, sqrt(v/(1-theta^2)), log=TRUE)
  for(i in 1:(S-1)){
    logdensity <- logdensity + dnorm(tilde_p[i+1], theta*tilde_p[i], v^0.5, log=TRUE)
  }
  return(logdensity)
}

logprior_alpha <- function(alpha){
  return(dnorm(alpha, 20, 1, log=TRUE))
}

logprior_beta <- function(beta){
  return(dnorm(beta, 0.5, 0.1, log=TRUE))
}

logprior_sigma2_0 <- function(alpha, beta, lA, sigma2_0){
  
  fisher_info <- 127.5/(sigma2_0^2) 
  return(log(sqrt(fisher_info)))
}

logprior_lA <- function(lA, tau, tilde_p, mu, s){
  p <- exp(tilde_p)/(1+exp(tilde_p))
  valid_index <- c()
  for(i in 1:nr.old){
    if((p[i]>0) & (p[i]<1)){
      valid_index <- c(valid_index, i)
    }
  }
  llk <- 0
  for(i in valid_index){
    p_index <- (tau[i]+7100 - (tau[i]+7100)%%100)/100 + 1
    if(lA[i] == 0){
      llk <- llk + log(p[p_index]) 
    }
    else{
      llk <- llk + log(1-p[p_index]) + dgamma(lA[i], shape=(mu^2/s^2), rate=(mu/s^2), log=TRUE)
    }
  }
  return(llk)
}

logprior_tau <- function(tau){
  llk <- 0
  for(i in 1:nr.old){
    llk <- llk - log(upper_tau[i]-lower_tau[i])
  }
  return(llk)
}


logprior_mu <- function(mu){
  return(dexp(mu, rate=1/mu_mean, log=TRUE))
}

logprior_s <- function(s){
  return(dexp(s, rate=1/b, log=TRUE))
}

# prior(v) exp(1)
logprior_v <- function(v){
  return(dexp(v, rate=1, log=TRUE))
}

# prior(\theta) = Uniform (0,1)
logprior_theta <- function(theta){
  return(0)
}

# likelihood
loglkd_yM <- function(alpha, beta, sigma2_0){
  llk <- 0
  for (i in 1:nr.mod) {
    llk <- llk + dnorm(y_M[i], mean=alpha+beta*lM[i], sd=sqrt(sigma2_0), log=TRUE)
  }
  return(llk)
}

loglkd_yA <- function(alpha, beta, sigma2_0, lA){
  llk <- 0
  for (i in 1:nr.old) {
    llk <- llk + dnorm(y_A[i], mean=alpha+beta*lA[i], sd=sqrt(sigma2_0), log=TRUE)
  }
  return(llk)
}

# posterior
logposterior <- function(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v){
  return(loglkd_yM(alpha, beta, sigma2_0) + loglkd_yA(alpha, beta, sigma2_0, lA) +
           logprior_alpha(alpha) + logprior_beta(beta) + logprior_sigma2_0(alpha, beta, lA, sigma2_0) +
           logprior_lA(lA, tau, tilde_p, mu, s) + logprior_mu(mu) + logprior_s(s) +
           logprior_tau(tau) + logprior_tildep(tilde_p, theta, v) + logprior_v(v) + logprior_theta(theta))
}

update_all <- function(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v,
                       rw_alpha, rw_beta, rw_sigma2_0, rw_tildep, rw_mu, rw_s, rw_theta, rw_v){
  
  logdensity0 <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v)
  
  # randomly choose a variable to update
  n <- sample(1:206, 1)
  
  # update alpha
  if(n == 1){ 
    alpha_new <- rnorm(1, alpha, rw_alpha)
    logdensity1 <- logposterior(alpha_new, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v)
    r <- logdensity1 - logdensity0
    if(log(runif(1)) < r){
      alpha <- alpha_new
    }
  }
  
  # update beta
  if(n == 2){
    beta_new <- rnorm(1, beta, rw_beta)
    logdensity1 <- logposterior(alpha, beta_new, sigma2_0, lA, tau, tilde_p, mu, s, theta, v)
    r <- logdensity1 - logdensity0
    if(log(runif(1)) < r){
      beta <- beta_new
    }
  }
  
  # update sigma0
  if(n == 3){
    sigma2_0_new <- rnorm(1, sigma2_0, rw_sigma2_0)
    if(sigma2_0_new >= 0){
      logdensity1 <- logposterior(alpha, beta, sigma2_0_new, lA, tau, tilde_p, mu, s, theta, v)
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r){
        sigma2_0 <- sigma2_0_new
      }
    }
  }
  
  # update lA
  if((n >= 4) & (n <= 66)){
    lA_new <- lA
    k <- n-3
    #p <- exp(tilde_p[k])/(1+exp(tilde_p[k]))
    if(lA_new[k] == 0){
      #lA_new[k] <- rgamma(1, shape=mu^2/(s^2), rate=mu/(s^2))
      lA_new[k] <- sample(c(0, rgamma(1, shape=mu^2/(s^2), rate=mu/(s^2))), 1, prob=c(0.5, 0.5))
      logdensity1 <- logposterior(alpha, beta, sigma2_0, lA_new, tau, tilde_p, mu, s, theta, v)
      if(lA_new[k] == 0){
        r <- min(1, logdensity1-logdensity0)
      }
      else{
        r <- min(1, logdensity1 - logdensity0 - dgamma(lA_new[k], shape=mu^2/(s^2), rate=mu/(s^2), log=TRUE)) 
      }
    }
    else{
      lA_new[k] <- sample(c(0, rgamma(1, shape=mu^2/(s^2), rate=mu/(s^2))), 1, prob=c(0.5, 0.5))
      logdensity1 <- logposterior(alpha, beta, sigma2_0, lA_new, tau, tilde_p, mu, s, theta, v)
      if(lA_new[k] == 0){
        r <- min(1, logdensity1 + dgamma(lA[k], shape=mu^2/(s^2), rate=mu/(s^2), log=TRUE) - logdensity0)
      }
      else{
        r <- min(1, logdensity1 + dgamma(lA[k], shape=mu^2/(s^2), rate=mu/(s^2), log=TRUE) - logdensity0 - dgamma(lA_new[k], shape=mu^2/(s^2), rate=mu/(s^2), log=TRUE))
      }
    }
    if(log(runif(1)) < r){
      lA <- lA_new
    }
  }
  
  # update tau
  if((n >= 67) & (n <= 129)){
    tau_new <- tau
    k <- n-66
    tau_new[k] <- runif(1, lower_tau[k], upper_tau[k])
    logdensity1 <- logposterior(alpha, beta, sigma2_0, lA, tau_new, tilde_p, mu, s, theta, v)
    r <- logdensity1 - logdensity0
    if(log(runif(1)) < r){
      tau <- tau_new
    }
  }
  
  # update tilde_p
  if((n >= 130) & (n <= 202)){
    tilde_p_new <- tilde_p
    k <- n-129
    tilde_p_new[k] <- rnorm(1, tilde_p_new[k], rw_tildep)
    logdensity1 <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p_new, mu, s, theta, v)
    r <- logdensity1 - logdensity0
    if(log(runif(1)) < r){
      tilde_p <- tilde_p_new
    }
  }
  
  # update mu
  if(n == 203){
    mu_new <- rnorm(1, mu, rw_mu)
    if(mu_new > 0){
      logdensity1 <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p, mu_new, s, theta, v)
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r){
        mu <- mu_new
      }
    }
  }
  
  # update s
  if(n == 204){
    s_new <- rnorm(1, s, rw_s)
    if(s_new > 0){
      logdensity1 <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s_new, theta, v)
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r){
        s <- s_new
      }
    }
  }
  
  # update theta
  if(n == 205){
    theta_new <- rnorm(1, theta, rw_theta)
    if((theta_new >= 0) & (theta_new <= 1)){
      logdensity1 <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta_new, v)
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r){
        theta <- theta_new
      }
    }
  }
  
  # update v
  if(n == 206){
    v_new <- rnorm(1, v, rw_v)
    if(v > 0){
      logdensity1 <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v_new)
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r){
        v <- v_new
      }
    }
  }
  
  llk <- logposterior(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v)
  list("alpha"=alpha, "beta"=beta, "sigma2_0"=sigma2_0, "lA"=lA, "tau"=tau, 
       "tilde_p"=tilde_p, "mu"=mu, "s"=s, "theta"=theta, "v"=v, "llk"=llk)
}

sampling_all <- function(iter, thin, rw_alpha, rw_beta, rw_sigma2_0, rw_tildep, rw_mu, rw_s, rw_theta, rw_v){
  
  alpha.samples <- beta.samples <- sigma2_0.samples <- mu.samples <- s.samples <- theta.samples <- v.samples <- llk.samples <- numeric(iter/thin)
  lA.samples <- tau.samples <- matrix(NA, nrow=iter/thin, ncol=nr.old)
  tilde_p.samples <- matrix(NA, nrow=iter/thin, ncol=S)
  
  alpha <- alpha.samples[1] <- rnorm(1, 20, 1)
  beta <- beta.samples[1] <- rnorm(1, 0.5, 0.1)
  sigma2_0 <- sigma2_0.samples[1] <- runif(1, 3, 4)
  mu <- mu.samples[1] <- rexp(1, rate=1/mu_mean)
  s <- s.samples[1] <- rexp(1, rate=1/b)
  theta <- theta.samples[1] <- runif(1, 0, 1)
  v <- v.samples[1] <- runif(1, 0, 10)
  
  for(i in 1:nr.old){
    lA.samples[1,i] <- sample(c(0, rgamma(1, mu^2/(s^2), mu/(s^2))), 1)
    tau.samples[1,i] <- runif(1, lower_tau[i], upper_tau[i])
  }
  
  lA <- lA.samples[1,]
  tau <- tau.samples[1,]
  
  for(i in 1:S){
    tilde_p.samples[1,i] <- runif(1, -10, 10)
  }
  tilde_p <- tilde_p.samples[1,]
  
  result <- update_all(alpha, beta, sigma2_0, lA, tau, tilde_p, mu, s, theta, v,
                       rw_alpha, rw_beta, rw_sigma2_0, rw_tildep, rw_mu, rw_s, rw_theta, rw_v)
  llk.samples[1] <- logposterior(result$alpha, result$beta, result$sigma2_0, result$lA, result$tau, result$tilde_p, result$mu, result$s, result$theta, result$v)
  k <- 1
  thin_round <- 0
  
  for(i in 2:iter){
    k <- k+1
    result <- update_all(result$alpha, result$beta, result$sigma2_0, result$lA, result$tau, result$tilde_p, result$mu, result$s, result$theta, result$v,
                         rw_alpha, rw_beta, rw_sigma2_0, rw_tildep, rw_mu, rw_s, rw_theta, rw_v)
    if(k == thin){
      k <- 0
      thin_round <- thin_round + 1
      alpha.samples[thin_round] <- result$alpha
      beta.samples[thin_round] <- result$beta
      sigma2_0.samples[thin_round] <- result$sigma2_0
      lA.samples[thin_round,] <- result$lA
      tau.samples[thin_round,] <- result$tau
      tilde_p.samples[thin_round,] <- result$tilde_p
      mu.samples[thin_round] <- result$mu
      s.samples[thin_round] <- result$s
      theta.samples[thin_round] <- result$theta
      v.samples[thin_round] <- result$v
      llk.samples[thin_round] <- result$llk
    }
    
    # save data in progress
    if((thin_round %in% c(10, 50, 100, 200, 400, 600, 800, 1000)) & (k == 0)){
      print(thin_round)
      write.csv(data.frame(alpha.samples, beta.samples, sigma2_0.samples, 
                           mu.samples, s.samples, theta.samples, v.samples), file="old_single_inprogress.csv")
      write.csv(tilde_p.samples, file="tildep_old_inprogress.csv")
      write.csv(lA.samples, file="lA_old_inprogress.csv")
      write.csv(tau.samples, file="tau_old_inprogess.csv")
    }
    
  } 
  list("alpha"=alpha.samples, "beta"=beta.samples, "sigma2_0"=sigma2_0.samples,
       "lA"=lA.samples, "tau"=tau.samples, "tilde_p"=tilde_p.samples, "mu"=mu.samples, "s"=s.samples,
       "theta"=theta.samples, "v"=v.samples, "llk"=llk.samples)
}