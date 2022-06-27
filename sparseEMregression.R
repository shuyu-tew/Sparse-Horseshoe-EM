# ==============================================================================
# Estimate tau2 & lambda2, depending on prior
# ==============================================================================

est.tau2.lambda2.sigma2 <- function(prior, model, Eb2, n = NULL, E.RSS = NULL, Enu.inv = NULL, Exi.inv = NULL, log.tau2.max = 20, tau2=1, sigma2=1)
{
  # ----------------------------
  # Gaussian regression
  if (model == "gaussian")
  {
    # Find tau2/lambda2 from Eb2 depending on the prior
    if (prior == "hs")
    {
      rv = find.tau2.sigma2.hs(Eb2,n,E.RSS, log.tau2.max=log.tau2.max,tau2=tau2,sigma2=sigma2)
      tau2 = rv$tau2
      sigma2 = rv$sigma2
      lambda2 = est.lambda2.hs(Eb2, sigma2, tau2)
    }  
    else if (prior == "hs.lambda2")
    {
      rv = find.tau2.sigma2.hs.lambda2(Eb2,n,E.RSS, log.tau2.max=log.tau2.max,tau2=tau2,sigma2=sigma2)
      tau2 = rv$tau2
      sigma2 = rv$sigma2
      lambda2 = est.lambda2.hs.lambda2(Eb2, sigma2, tau2)      
    }
    else if (prior == "lasso")
    {
      rv = find.tau2.sigma2.lasso(Eb2,n,E.RSS, log.tau2.max=log.tau2.max,tau2=tau2,sigma2=sigma2)
      tau2 = rv$tau2
      sigma2 = rv$sigma2
      lambda2 = est.lambda2.lasso(Eb2, sigma2, tau2)      
    }
    else if (prior == "ridge.tau2")
    {
      sigma2 = find.sigma2.ridge.tau2(Eb2,n,E.RSS)
      tau2 = find.tau2.ridge.tau2(Eb2, sigma2)
      lambda2 = rep(1,length(Eb2))
    }
  }
  
  # ----------------------------
  # GLMs -- no sigma2  
  else
  {
    # Find tau2/lambda2 from Eb2 depending on the prior
    if (prior == "hs")
    {
      tau2 = find.tau2.hs(Eb2, sigma2=1, log.tau2.max=log.tau2.max)
      lambda2 <- est.lambda2.hs(Eb2, sigma2=1, tau2)
    }
    else if (prior == "hs.lambda2")
    {
      tau2 = find.tau2.hs.lambda2(Eb2, sigma2=1, log.tau2.max=log.tau2.max)
      lambda2 = est.lambda2.hs.lambda2(Eb2, sigma=1, tau2)
    }
    else if (prior == "lasso")
    {
      tau2 = find.tau2.lasso(Eb2, sigma2=1, log.tau2.max=log.tau2.max)
      lambda2 <- est.lambda2.lasso(Eb2, sigma2=1, tau2)
    }
    else if (prior == "lasso.lambda")
    {
      tau2 = find.tau2.lasso.lambda(Eb2, sigma2=1, log.tau2.max=log.tau2.max)
      lambda2 <- est.lambda2.lasso.lambda(Eb2, sigma2=1, tau2)
    }      
    else if (prior == "1.on.lambda2")
    {
      tau2 = find.tau2.lasso.1.on.lambda2(Eb2, sigma2=1, log.tau2.max=log.tau2.max)
      lambda2 <- est.lambda2.lasso.1.on.lambda2(Eb2, sigma2=1, tau2)
    }      
    else if (prior == "hs.igig")
    {
      tau2 = find.tau2.hs.igig(Eb2, Enu.inv, Exi.inv, sigma2=1, log.tau2.max=log.tau2.max)
      lambda2 = est.lambda2.hs.igig(Eb2, Enu.inv, sigma2=1, tau2)
    }
    else if (prior == "ridge")
    {
      tau2 = find.tau2.ridge.tau(Eb2, sigma2=1)
      lambda2 = rep(1,length(Eb2))
    }
    else if (prior == "ridge.tau2")
    {
      tau2 = find.tau2.ridge.tau2(Eb2, sigma2=1)
      lambda2 = rep(1,length(Eb2))
    }
    
    sigma2 = NULL
  }
  
  lambda2 = pmax(lambda2, 1e-20)
  
  return(list(lambda2=lambda2,tau2=tau2,sigma2=sigma2))
}

# ------------------------------------------------------------------------------
# Horseshoe (lambda)
tau2.f.hs <- function(log.tau2, Eb2, sigma2=1)
{
  p = length(Eb2)
  tau2 = exp(log.tau2)
  lambda2 = est.lambda2.hs(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + log(1+tau2) + sum(log(1+lambda2))
  f
}

tau2.sigma2.f.hs <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  sigma2 = exp(x[2])
  lambda2 = est.lambda2.hs(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + log(1+tau2) + sum(log(1+lambda2))
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2
  f
}

tau2.sigma2.f.hs.2 <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  sigma2 = Erss/n
  lambda2 = est.lambda2.hs(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + log(1+tau2) + sum(log(1+lambda2))
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2 + log(sigma2)
  f
}

est.lambda2.hs <- function(Eb2,sigma2,tau2)
{
  W = Eb2/2/sigma2/tau2
  
  pmax((sqrt(-1 + 2*W + sqrt(1 + 20*W + 4*W^2))/sqrt(6))^2, 1e-10)
}

find.tau2.hs <- function(Eb2,sigma2,log.tau2.max=20)
{
  exp((optimise(function(log.tau2){tau2.f.hs(log.tau2,Eb2,sigma2)}, c(-10,log.tau2.max)))$minimum)
}

#find.tau2.sigma2.hs <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1)
#{
#  min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.hs(x, Eb2, n, E.RSS) },
#lower = c(-10,-10), upper = c(log.tau2.max,10), method="L-BFGS-B"))$par)
#min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.hs(x, Eb2, n, E.RSS) }))$par)

#  tau2 = min.par[1]
#  sigma2 = min.par[2]

#  return(list(tau2=tau2,sigma2=sigma2))
#}

find.tau2.sigma2.hs <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1)
{
  #min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.hs(x, Eb2, n, E.RSS) },
  #lower = c(-10,-10), upper = c(log.tau2.max,10), method="L-BFGS-B"))$par)
  #min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.hs(x, Eb2, n, E.RSS) }))$par)
  
  #tau2 = exp((optimise(function(log.tau2){tau2.sigma2.f.hs.2(log.tau2,Eb2, n, Erss)}, c(-10,log.tau2.max)))$minimum)
  
  #tau2 = min.par[1]
  p = length(Eb2)
  
  tau2 = exp((optimise(function(log.tau2){tau2.sigma2.f.hs.2(log.tau2,Eb2,n,E.RSS)}, c(-10,log.tau2.max)))$minimum)
  #sigma2 = E.RSS/n
  #lambda2 = est.lambda2.hs.lambda2(Eb2,sigma2,tau2)
  #sigma2 = (E.RSS + sum(Eb2/lambda2/tau2))/(n+p)
  #sigma2 = (E.RSS + sum(Eb2/lambda2/tau2))/(n+p)
  sigma2 = E.RSS/n
  
  return(list(tau2=tau2,sigma2=sigma2))  
}


# ------------------------------------------------------------------------------
# Horseshoe (lambda2)
tau2.f.hs.lambda2 <- function(log.tau2, Eb2, sigma2=1)
{
  p = length(Eb2)
  tau2 = exp(log.tau2)
  lambda2 = est.lambda2.hs.lambda2(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum( (1/2)*log(lambda2) + log(lambda2+1) )
  f = f + (1/2)*log(tau2) + log(1+tau2)
  f
}

tau2.sigma2.f.hs.lambda2 <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  sigma2 = exp(x[2])
  lambda2 = est.lambda2.hs.lambda2(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum( (1/2)*log(lambda2) + log(lambda2+1) )
  f = f + (1/2)*log(tau2) + log(1+tau2)
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2
  f
}

tau2.sigma2.f.hs.lambda2.2 <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  sigma2 = Erss/n
  lambda2 = est.lambda2.hs.lambda2(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum( (1/2)*log(lambda2) + log(lambda2+1) )
  f = f + (1/2)*log(tau2) + log(1+tau2)
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2 + log(sigma2)
  if (is.nan(f))
  {
    f
  }
  f
}

est.lambda2.hs.lambda2 <- function(Eb2,sigma2,tau2)
{
  W = Eb2/2/sigma2/tau2
  
  pmax( (sqrt(W^2 + 6*W + 1) + W - 1) / 4, 1e-50  )
}

find.tau2.hs.lambda2 <- function(Eb2,sigma2,log.tau2.max=20)
{
  exp((optimise(function(log.tau2){tau2.f.hs.lambda2(log.tau2,Eb2,sigma2)}, c(-10,log.tau2.max)))$minimum)
}

# find.tau2.sigma2.hs.lambda2 <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1)
# {
#  min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.hs.lambda2(x, Eb2, n, E.RSS) },
#                        lower = c(-10,-10), upper = c(log.tau2.max,10), method="L-BFGS-B"))$par)
#  #min.par = exp((optim( (c(0,0)), function(x){ tau2.sigma2.f.hs.lambda2(x, Eb2, n, E.RSS) }                       ))$par)
# 
#  tau2 = min.par[1]
#  sigma2 = min.par[2]
# 
#  return(list(tau2=tau2,sigma2=sigma2))
# }

find.tau2.sigma2.hs.lambda2 <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1)
{
  p = length(Eb2)
  tau2 = exp((optimise(function(log.tau2){tau2.sigma2.f.hs.lambda2.2(log.tau2,Eb2, n, E.RSS)}, c(-10,log.tau2.max)))$minimum)
  #sigma2 = (E.RSS+p)/(n+p)
  sigma2 = E.RSS/n
  #lambda2 = est.lambda2.hs.lambda2(Eb2,sigma2,tau2)
  #sigma2 = (E.RSS + sum(Eb2/lambda2/tau2))/(n+p)
  
  return(list(tau2=tau2,sigma2=sigma2))
}

# ------------------------------------------------------------------------------
# Horseshoe (IGIG)
tau2.f.hs.igig <- function(log.tau2, Eb2, Enu.inv, Exi.inv, sigma2=1)
{
  p = length(Eb2)
  tau2 = exp(log.tau2)
  lambda2 = est.lambda2.hs.igig(Eb2,Enu.inv,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum(3/2*log(lambda2) + Enu.inv/lambda2) + 3/2*log(tau2) + Exi.inv/tau2
  if (is.infinite(f) || is.na(f))
  {
    f
  }
  f
}

est.lambda2.hs.igig <- function(Eb2, Enu.inv, sigma2, tau2)
{
  W = Eb2/2/sigma2/tau2
  pmax( (W + Enu.inv)/2, 1e-50 )
}

find.tau2.hs.igig <- function(Eb2,Enu.inv,Exi.inv,sigma2,log.tau2.max=20)
{
  exp((optimise(function(log.tau2){tau2.f.hs.igig(log.tau2,Eb2,Enu.inv,Exi.inv,sigma2)}, c(-10,log.tau2.max)))$minimum)
}

# ------------------------------------------------------------------------------
# Lasso (lambda2)
tau2.f.lasso <- function(log.tau2, Eb2, sigma2=1)
{
  p = length(Eb2)
  tau2 = exp(log.tau2)
  lambda2 = est.lambda2.lasso(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum(lambda2/2) 
  f = f + 2*log(tau2) + 1/tau2
  #f = f + log(1+tau2) + (1/2)*log(tau2)
  f
}

tau2.sigma2.f.lasso <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  sigma2 = exp(x[2])
  lambda2 = est.lambda2.lasso(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum(lambda2/2) + 2*log(tau2) + 1/tau2
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2
  f
}

tau2.sigma2.f.lasso.2 <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  sigma2 = Erss/n
  lambda2 = est.lambda2.lasso(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum(lambda2/2) + 2*log(tau2) + 1/tau2
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2
  f
}

tau2.sigma2.f.lasso.3 <- function(x, Eb2, n, Erss)
{
  p = length(Eb2)
  tau2 = exp(x[1])
  lambda2 = est.lambda2.lasso(Eb2,sigma2=1,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/tau2*sum(Eb2/lambda2) 
  f = f + sum(lambda2/2) + 2*log(tau2) + 1/tau2
  f
}

est.lambda2.lasso <- function(Eb2, sigma2, tau2)
{
  W = Eb2/2/sigma2/tau2
  pmax((sqrt(1 + 8*W)-1)/2, 1e-50)
}

find.tau2.lasso <- function(Eb2,sigma2,log.tau2.max=20)
{
  exp((optimise(function(log.tau2){tau2.f.lasso(log.tau2,Eb2,sigma2)}, c(-10,log.tau2.max)))$minimum)
}

#find.tau2.sigma2.lasso <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1)
#{
#  min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.lasso(x, Eb2, n, E.RSS) },
#                        lower = c(-10,-10), upper = c(log.tau2.max,10), method="L-BFGS-B"))$par)
#  tau2 = min.par[1]
#  sigma2 = min.par[2]
#  
#  return(list(tau2=tau2,sigma2=sigma2))
#}

find.tau2.sigma2.lasso <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1)
{
  #min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f.lasso(x, Eb2, n, E.RSS) },
  #                    lower = c(-10,-10), upper = c(log.tau2.max,10), method="L-BFGS-B"))$par)
  
  
  #tau2 = exp((optimise(function(log.tau2){tau2.sigma2.f.lasso.3(log.tau2,Eb2, n, E.RSS)}, c(-10,log.tau2.max)))$minimum)
  tau2 = exp((optimise(function(log.tau2){tau2.sigma2.f.lasso.2(log.tau2,Eb2, n, E.RSS)}, c(-10,log.tau2.max)))$minimum)
  
  #tau2 = min.par[1]
  sigma2 = E.RSS/n
  
  return(list(tau2=tau2,sigma2=sigma2))
}


# ------------------------------------------------------------------------------
# Lasso (lambda)
tau2.f.lasso.lambda <- function(log.tau2, Eb2, sigma2=1)
{
  p = length(Eb2)
  tau2 = exp(log.tau2)
  lambda2 = est.lambda2.lasso.lambda(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum(lambda2/2) - sum(log(lambda2)/2) + 2*log(tau2) + 1/tau2
  if (is.infinite(f) || is.na(f))
  {
    f
  }
  f
}

est.lambda2.lasso.lambda <- function(Eb2, sigma2, tau2)
{
  W = Eb2/2/sigma2/tau2
  pmax(sqrt(2)*sqrt(W), 1e-50)
}

find.tau2.lasso.lambda <- function(Eb2,sigma2,log.tau2.max=20)
{
  exp((optimise(function(log.tau2){tau2.f.lasso.lambda(log.tau2,Eb2,sigma2)}, c(-10,log.tau2.max)))$minimum)
}


# ------------------------------------------------------------------------------
# Lasso (1.on.lambda)
tau2.f.lasso.1.on.lambda2 <- function(log.tau2, Eb2, sigma2=1)
{
  p = length(Eb2)
  tau2 = exp(log.tau2)
  lambda2 = est.lambda2.lasso.1.on.lambda2(Eb2,sigma2,tau2)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) 
  f = f + sum(2*log(1/lambda2) + 1/2/(1/lambda2) ) + 2*log(tau2) + 1/tau2
  if (is.infinite(f) || is.na(f))
  {
    f
  }
  f
}

est.lambda2.lasso.1.on.lambda2 <- function(Eb2, sigma2, tau2)
{
  W = Eb2/2/sigma2/tau2
  pmax( 1/( (sqrt(8*W+9)-3)/(4*W) ), 1e-20 )
}

find.tau2.lasso.1.on.lambda2 <- function(Eb2,sigma2,log.tau2.max=20)
{
  exp((optimise(function(log.tau2){tau2.f.lasso.1.on.lambda2(log.tau2,Eb2,sigma2)}, c(-10,log.tau2.max)))$minimum)
}

# ------------------------------------------------------------------------------
# Ridge (tau)
find.tau2.ridge.tau <- function(Eb2,sigma2)
{
  p = length(Eb2)
  W = sum(Eb2/sigma2/2)
  
  (sqrt(4*W^2+4*p*W+16*W+p^2)/(p+2)+(2*W)/(p+2)-p/(p+2))/2
}

# ------------------------------------------------------------------------------
# Ridge (tau2)
find.tau2.ridge.tau2 <- function(Eb2,sigma2)
{
  p = length(Eb2)
  W = sum(Eb2/sigma2/2)
  
  (sqrt(4*W^2+(4*p+20)*W+p^2+2*p+1)+2*W-p-1)/(2*p+6)
}


sigma2.f.ridge.tau2 <- function(log.sigma2, Eb2, n, E.RSS)
{
  p = length(Eb2)
  sigma2 = exp(log.sigma2)
  tau2 = find.tau2.ridge.tau2(Eb2,sigma2)
  
  f = p/2*log(tau2) + 1/2/sigma2/tau2*sum(Eb2)
  f = f + (1/2)*log(tau2) + log(1+tau2)
  f = f + (n+p)/2*log(sigma2) + E.RSS/2/sigma2
  f
}

find.sigma2.ridge.tau2 <- function(Eb2,n,E.RSS)
{
  exp((optimise(function(log.sigma2){sigma2.f.ridge.tau2(log.sigma2,Eb2,n,E.RSS)}, c(-10,20)))$minimum)
}


# coord.wise.opt <- function(Eb2,sigma2,max.iter=1000)
# {
#   tau2 = 1
#   for (i in 1:max.iter)
#   {
#     lambda2 = est.lambda2(Eb2,sigma2,tau2)
#     tau2 = est.tau2(Eb2,lambda2,sigma2)
#     print(log(tau2))
#   }
#   list(lambda2=lambda2,tau2=tau2)
# }


# ==============================================================================
# Multiple means ...
# ==============================================================================
#

mu.est <- function(y, sigma2, tau2)
{
  # Initialise
  Eb2 = y^2 + sigma2
  
  # EM
  done = F
  i = 1;
  while (!done)
  {
    # M-step
    #tau2 <- find.tau2(Eb2, sigma2)
    lambda2 <- est.lambda2.hs(Eb2, sigma2, tau2)
    lambda2 = pmax(lambda2, 1e-6)
    
    # E-step
    kappa = 1/(1+tau2*lambda2)
    beta.hat = y*(1-kappa)
    
    Eb2 = beta.hat^2 + (1-kappa)^2*sigma2
    
    # Termination conditions
    i = i+1
    if (i == 3000)
    {
      break
    }
  }
  
  beta.hat[abs(beta.hat)<1e-4] = 0
  
  return(list(beta=beta.hat,lambda2=lambda2,tau2=tau2))
}


mu.est.kappa <- function(y, sigma2, tau2)
{
  # Initialise
  Eb2 = y^2 + sigma2
  #Eb2 = 0
  
  Eb2 = 0.5
  
  # EM
  done = F
  i = 1;
  while (!done)
  {
    # M-step
    #tau2 <- find.tau2(Eb2, sigma2)
    #lambda2 <- est.lambda2(Eb2, sigma2, tau2)
    #lambda2 = pmax(lambda2, 1e-6)
    k = pmax(1 - Eb2, 1e-10)
    lambda2 = sqrt((1-k)/k)
    
    # E-step
    kappa = 1/(1+tau2*lambda2)
    beta.hat = y*(1-kappa)
    
    Eb2 = beta.hat^2 + (1-kappa)^2*sigma2
    
    # Termination conditions
    i = i+1
    if (i == 3000)
    {
      break
    }
  }
  
  beta.hat[abs(beta.hat)<1e-4] = 0
  
  return(list(beta=beta.hat,lambda2=lambda2,tau2=tau2))
}

# ==============================================================================
# EM algorithm
# ==============================================================================

beta.EM <- function(formula, data, model, prior = "hs", approx = F, e = 1e-6, 
                    stochastic = F, PG = F, rv.br = NULL, log.tau2.max=15, start.values = NULL, 
                    use.glmnet = F, Eb2.start = NULL)
{
  rv = list()
  
  # -------------------------------------------------------------------    
  # Process and set up the data from the model formula
  rv$terms <- stats::terms(x = formula, data = data)
  
  mf = stats::model.frame(formula = formula, data = data)
  rv$target.var = names(mf)[1]
  
  if (model == 'binomial')
  {
    # Check to ensure target is a factor
    if (!is.factor(mf[,1]) || (!is.factor(mf[,1]) && length(levels(mf[,1])) != 2))
    {
      stop('Target variable must be a factor with two levels for logistic regression')
    }
    
    rv$ylevels <- levels(mf[,1])
    mf[,1] = as.numeric(mf[,1])
    mf[,1] = (as.numeric(mf[,1]) - 1)
  }  
  
  y = mf[,1]
  X = stats::model.matrix(formula, data=data)
  
  # Convert to a numeric matrix and drop the target variable
  X = as.matrix(X)
  X = X[,-1,drop=FALSE]
  X0 = X
  
  std.X  = my.standardise(X)
  I.keep = std.X$std.X!=0
  X      = std.X$X[,I.keep]
  
  y.mu  = 0
  y.sd  = 1  
  if (model == "gaussian")
  {
    y.mu  = mean(y)
    y     = y - y.mu
    y.sd  = sqrt( mean( (y-mean(y))^2 ) )
    y     = y / y.sd
  }
  
  ## EM-stuff
  if (model == "gaussian")
  {
    rv.em = beta.est(X, y, prior, approx = approx, 
                     e = e, stochastic = stochastic, log.tau2.max = log.tau2.max,
                     use.glmnet = use.glmnet, start.values = start.values)
  }
  else
  {
    if (!is.null(rv.br))
    {
      #rv.start = find.start.glm(rv.br, X0, y)
    }
    else
    {
      rv.em = beta.est.glm(X, y, model, prior, approx = approx, 
                           e = e, stochastic = stochastic, PG = PG, log.tau2.max = log.tau2.max, 
                           Eb2.start=Eb2.start, use.glmnet = use.glmnet)
    }
  }
  
  # Done
  rv$beta  = rv.em$beta / (std.X$std.X[I.keep]) * y.sd
  rv$beta0 = rv.em$beta0 + y.mu
  rv$beta0 = rv$beta0 - std.X$mean.X[I.keep] %*% rv$beta
  rv$rv.em = rv.em
  rv$model = model
  rv$I.keep = I.keep
  
  rv
}


# ==============================================================================
# Find the posterior mode in a Gaussian regression
# ==============================================================================
#
beta.est <- function(X, y, prior, approx = F, use.glmnet = F, e = 1e-6, stochastic = F, 
                     start.values = NULL, log.tau2.max = 20)
{
  p = ncol(X)
  n = nrow(X)
  
  # Initialise
  XtX = t(X) %*% X
  Xty = t(X) %*% y
  
  Enu.inv = rep(1,p)
  Exi.inv = 1  
  
  tau2 = 1
  if (is.null(start.values))
  {
    #rv = linreg.post(XtX, Xty, X, y, rep(1/tau2,p), 1)
    #sigma2 = sum( (y-X %*% rv$beta.mu - mean(y))^2 ) / (n-1) 
    
    beta.mu = Xty/diag(XtX)
    Eb2 = beta.mu^2
    rv = list(beta.mu = beta.mu)
    
    #E.RSS = n*sigma2
    E.RSS = 1e10
    
  }
  else
  {
    Eb2 = start.values$Eb2
    E.RSS = start.values$E.RSS
    rv = list(beta.mu = rep(0,p))
  }
  
  t = 1/sqrt(n)/5
  
  
  # EM
  done = F
  i = 1
  lambda2 = matrix(0,p,1)
  while (!done && i < 1e4)
  {
    beta.old = rv$beta.mu
    lambda2.old = lambda2
    tau2.old = tau2
    
    # M-step
    rv.hyp = est.tau2.lambda2.sigma2(prior, "gaussian", Eb2, n, E.RSS, 
                                     Enu.inv, Exi.inv, log.tau2.max=log.tau2.max,
                                     tau2=tau2, sigma2=sigma2)
    
    tau2 = rv.hyp$tau2
    sigma2 = rv.hyp$sigma2
    lambda2 = rv.hyp$lambda2
    
    lambda2 = pmax(lambda2, 1e-20)
    
    d = rep(1/(tau2*lambda2*sigma2))
    if (sum(is.nan(d)))
    {
      d
    }
    
    # E-step
    if (!stochastic)
    {
      if (i < 10 || (approx == F && use.glmnet == F))
      {
        rv = linreg.post(XtX, Xty, X, y, rep(1/(tau2*lambda2*sigma2)), sigma2, approx = F, sample = F)
      }
      else
      {
        if (!use.glmnet)
        {
          rv = linreg.post(XtX, Xty, X, y, rep(1/(tau2*lambda2*sigma2)), sigma2, approx = T, sample = F)
        }
        else
        {
          rv = linreg.post.glmnet(X, y, lambda2, tau2, sigma2)
        }
      }
      
      Eb2 = rv$beta.mu^2 + rv$beta.v2
      E.RSS = rv$E.RSS
      if(E.RSS<0)
      {
        E.RSS
      }
    }
    # Stochastic E-Step
    else
    {
      rv = linreg.post(XtX, Xty, X, y, rep(1/(tau2*lambda2*sigma2)), sigma2, approx = T, sample = T)
      
      rv.stem = update.Eb2.E.RSS(Eb2, E.RSS, rv$b, X, y, i, chi=100)
      Eb2 = rv.stem$Eb2
      E.RSS = rv.stem$E.RSS
    }
    
    
    # Threshold
    if (prior != "ridge.tau" && prior != "ridge.tau2")
    {
      rv$beta.mu[abs(rv$beta.mu)<t] = 0
    }
    I = (abs(beta.old - rv$beta.mu) / (1+abs(rv$beta.mu))) < 1e-5
    
    # Termination conditions
    if (sum(abs(rv$beta.mu-beta.old)) / (1+sum(abs((rv$beta.mu)))) < e || sum(I) == p)
    {
      break
    }
    i = i+1
  }
  
  # Done -- final sparsification
  t = 1/sqrt(n)/5
  
  beta.hat = rv$beta.mu
  if (prior != "ridge.tau" && prior != "ridge.tau2")
  {
    beta.hat[abs(beta.hat)<t] = 0
  }
  
  # Compute the negative log-posterior
  L = E.neg.log.posterior(Eb2, rv$E.RSS, lambda2, tau2, sigma2, n, p)
  
  # Return
  return(list(beta=beta.hat,
              beta0=mean(y),
              lambda2=lambda2,
              tau2=tau2,
              sigma2=sigma2,
              num.iter = i,
              Eb2 = Eb2, 
              E.RSS = E.RSS,
              L = L
  ))
}

# ==============================================================================
# Estimate the beta's in a GLM
# ==============================================================================
#
beta.est.glm <- function(X, y, model, prior, approx = F, e = 1e-6, stochastic = F, Eb2.start = NULL, 
                         PG = F, log.tau2.max=20, use.glmnet = F)
{
  p = ncol(X)
  n = nrow(X)
  
  # Initialise
  X = cbind(matrix(1,n,1),X)
  lambda2 = matrix(0,p,1)
  
  Enu.inv = rep(1,p)
  Exi.inv = 1  
  
  # If no initial Eb2 provided
  tau2 = min(1, exp(log.tau2.max))
  if (is.null(Eb2.start))
  {
    #beta = coefficients(glmnet(X[,2:(p+1)], y, alpha=0, family="binomial", lambda=1/1e3))
    
    #rv = glmreg.post.glmnet(X[,2:(p+1)], y, model, rep(1,p), 1/n)
    #beta = rv$beta.mu
    
    rv = glm.fit(X, y, c(0, rep(1/1,p)), model)
    Eb2 = rv$beta.mu^2 + rv$beta.v2
    Eb2 = Eb2[2:(p+1)]
    beta = rv$beta.mu
  }
  # Starting point passed
  else
  {
    # Find tau2/lambda2 from starting Eb2
    rv.hyp = est.tau2.lambda2.sigma2(prior, model, Eb2.start, Enu.inv=Enu.inv, Exi.inv=Exi.inv, log.tau2.max=log.tau2.max)
    
    lambda2 = rv.hyp$lambda2
    tau2    = rv.hyp$tau2    
    
    # Find posterior mode
    rv = glm.fit(X, y, c(0, 1/(tau2*lambda2)), model)
    Eb2 = Eb2.start
  }
  
  # Initialise Stochastic elements  
  #beta.mu.stochastic = beta
  #Eb2.stochastic = beta^2
  
  # Initial linear predictor and conditional mean
  eta = X %*% beta
  if (model == "binomial")
  {
    mu = 1/(1+exp(-eta))
  }
  else
  {
    mu = exp(eta)
  }
  
  # Sparsity threshold  
  t = 1/sqrt(n)/5
  
  # EM algorithm
  done = F
  i = 1
  while (!done && i < 5e4)
  {
    # Previous values
    beta.old = beta
    lambda2.old = lambda2
    tau2.old = tau2
    
    # M-step
    if (stochastic)
    {
      # ...
    }
    else
    {
      rv.hyp = est.tau2.lambda2.sigma2(prior, model, Eb2, Enu.inv=Enu.inv, Exi.inv=Exi.inv, log.tau2.max=log.tau2.max)
      
      lambda2 = rv.hyp$lambda2
      tau2    = rv.hyp$tau2
    }
    lambda2 = pmin(lambda2,1e10)
    lambda2 = pmax(lambda2,1e-10)
    
    # E-step
    # Quadratic approximation
    if (!use.glmnet)
    {
      #omega2 = pgdraw.moments(1, eta)$mu
      #for (j in 1:100)
      #{
      #  omega2 = pgdraw(1, eta)
      #  z = (y-1/2)/omega2
      #  
      #  rv = linreg.w2.post(X, z, 1/omega2, c(0,1/(tau2*lambda2)), approx = approx, sample = T)
      #}
      
      rv.omega2 = compute.omega2(model,mu,eta,y,PG)
      omega2 = rv.omega2$omega2
      z = rv.omega2$z
      
      rv = linreg.w2.post(X, z, 1/omega2, c(0,1/(tau2*lambda2)), approx = approx, sample = F)
    }
    else
    {
      rv = glmreg.post.glmnet(X[,2:(p+1)], y, model, lambda2, tau2)
    }
    
    beta = rv$beta.mu
    eta = X %*% beta
    
    # Conditional means
    if (model == "binomial")
    {
      mu = 1/(1+exp(-eta))    
    }
    else if (model == "poisson")
    {
      mu = exp(eta)
    }
    
    # Expectations
    Eb2 = rv$beta.mu^2 + rv$beta.v2
    Eb2 = Eb2[2:(p+1)]
    if (prior != "ridge" && prior != "ridge.tau2")
    {
      #rv$beta.mu[abs(rv$beta.mu)<t] = 0
    }
    
    I = (abs(beta.old - rv$beta.mu) / (1+abs(rv$beta.mu))) < 1e-4
    
    #cat(sum(I),": ", sum(rv$beta.mu==0))
    #cat(rv$beta.mu)
    #cat("\n")
    
    # IG-IG
    if (prior == "hs.igig")
    {
      Enu.inv = lambda2/(1+lambda2)
      Exi.inv = tau2/(1+tau2)
    }
    
    if (sum(abs((beta.old)-(rv$beta.mu))) / (1+sum(abs((rv$beta.mu)))) < e)
    {
      break
    }
    if (sum(I) == p)
    {
      break
    }
    i = i+1
  }
  
  beta.hat = rv$beta.mu
  # Done -- final sparsification
  if (prior != "ridge" && prior != "ridge.tau2")
  {
    t = 1/sqrt(n)/5
    beta.hat[abs(beta.hat)<t] = 0
  }
  
  #beta.mu.stochastic[abs(beta.mu.stochastic) < 1e-2] = 0
  
  # Compute the negative log-posterior
  beta = beta.hat[2:(p+1)]
  beta0 = beta.hat[1]
  
  #NLL = neg.log.likelihood(X[,-1], y, beta, beta0, model)  
  #L = E.neg.log.posterior.glm(NLL, Eb2, lambda2, tau2)
  
  # Return
  return(list(beta=beta.hat[2:(p+1)],
              beta0=beta.hat[1],
              lambda2=lambda2,
              tau2=tau2,
              num.iter = i,
              Eb2 = Eb2
              #Eb2.stochastic = Eb2.stochastic, 
              #beta.mu.stochastic = beta.mu.stochastic
              #NLL = NLL,
              #L = L
  ))
}

compute.omega2 <- function(model, mu, eta, y, PG)
{
  if (model == "binomial")
  {
    if (!PG)
    {
      deta = 1/mu/(1-mu)
      z = eta + (y-mu)*deta
      omega2 = mu*(1-mu)
    }
    else
    {
      omega2 = pgdraw.moments(1, eta)$mu
      z = (y-1/2)/omega2
    }
  }
  else if (model == "poisson")
  {
    deta = 1/mu
    z = eta + (y-mu)*deta
    omega2 = mu      
  }  
  return(list(omega2=omega2,z=z))
}

update.Eb2.E.RSS <- function(Eb2, E.RSS, b, X, y, i, chi)
{
  alpha = (i/(i+1))^chi
  
  Eb2 = alpha*Eb2 + (1-alpha)*b^2
  E.RSS = alpha*E.RSS + (1-alpha)*(sum((y - X %*% b - mean(y))^2))
  
  return(list(Eb2=Eb2, E.RSS=E.RSS))
}

# ============================================================================================================================
# Linear regression conditional posterior statistics
linreg.post <- function(XtX, Xty, X, y, d, sigma2, approx = F, sample = F, return.E.RSS = T, bhat = F)
{
  p = ncol(X)
  n = nrow(X)
  
  bhat = F
  
  ## Rue's algorithm (p < 2n)
  if (!bhat)
  {
    L = chol(XtX/sigma2 + diag(d))
    beta.mu = backsolve(L,forwardsolve(t(L), Xty/sigma2))
    E.RSS = NULL
    
    # If a sample is requested
    if (sample)
    {
      b = beta.mu + backsolve(L, stats::rnorm(p,0,1))
    }
    else
    {
      b = NULL
    }
    
    # Exact first and second moments
    if (!approx)
    {
      Li       = forwardsolve(t(L), diag(1,p))
      beta.var = colSums(Li^2)
      
      if (return.E.RSS)
      {
        E.RSS    = sum((y - X%*%beta.mu - mean(y))^2) + sum(diag(XtX %*% Li%*%t(Li)))
      }
    }
    # Approximate first and second moments
    else
    {
      beta.var = 1/diag(XtX/sigma2+diag(d))
      if (return.E.RSS)
      {
        E.RSS = sum((y - X%*%beta.mu - mean(y))^2) + sum( diag(XtX)*beta.var )
      }
    }
  }
  
  ## Bhattarchaya's algorithm
  else
  {
    E.RSS = NULL
    
    #bayesLREM.beta.cond.posterior.bhat <- function(x, y, Lambda, quick = F, xtx){
    u     = as.matrix(stats::rnorm(p,0,1)) / sqrt(d)
    delta = as.matrix(stats::rnorm(n,0,1))
    
    v     = (X/sqrt(sigma2)) %*% u + delta
    
    XL    = X
    
    #performing x%*%diag(as.vector(Lambda))
    for (i in 1:length(d))
    {
      XL[,i] = XL[,i]/d[i]/sqrt(sigma2)
    }
    
    LXt = t(XL)
    
    W   = diag(n) + X %*% LXt/sqrt(sigma2)
    
    Up  = chol(W)
    tUp = t(Up)
    
    #sample
    b = NULL
    if (sample)
    {
      vv  = forwardsolve(tUp, (y/sqrt(sigma2)-v))
      w   = backsolve(Up, vv)
      b   = u + LXt %*% w
    }
    
    #beta mean
    vv_  = forwardsolve(tUp, y/sqrt(sigma2))
    w_   = backsolve(Up, vv_)
    beta.mu = LXt %*% w_
    
    #beta variance
    if(approx)
    {
      #beta.var = 1/(diag(XtX/sigma2) + d)
      XtX = colSums(X^2)
      beta.var = 1/(XtX/sigma2 + d)
      if (return.E.RSS)
      {
        E.RSS = sum((y - X %*% beta.mu - mean(y))^2) + sum( diag(XtX)*beta.var )
      }      
    }
    else
    {
      vvv      = forwardsolve(tUp, XL)
      omega    = backsolve(Up, vvv)
      beta.var = 1/d - colSums(t(LXt) * omega)
    }
  }
  
  #cat("E.RSS: ")
  #cat(E.RSS,"\n")
  #cat("bv2: ")
  #cat(beta.var,"\n")  
  
  return(list(beta.mu = beta.mu, beta.v2 = beta.var, E.RSS = E.RSS, b = b))
}

# ============================================================================================================================
# Linear regression conditional posterior statistics with heteroskedasticity
linreg.w2.post <- function(X, y, w2, d, approx = F, sample = F)
{
  p = ncol(X)
  X0 = X
  
  # Create weighted sufficient statistics
  for (i in 1:p)
  {
    X0[,i] = as.vector(X[,i]/w2)
  }  
  
  # Precompute if using Rue's algorithm
  #if (p < 2*n)
  {
    XtX = t(X0)%*%X
    Xty = t(X0)%*%y
  }
  
  # 
  linreg.post(XtX, Xty, X0, y, d, sigma2 = 1, approx = approx, sample = sample, return.E.RSS = F)
}

# ============================================================================================================================
# Simple GLM ridge fit
glm.fit <- function(X, y, d, model)
{
  n = nrow(X)
  p = ncol(X)
  
  # Initialise
  if (model == "binomial")
  {
    mu = (y + 1/2) / 2
    eta = log(mu/(1-mu))
  }
  else if (model == "poisson")
  {
    mu = y+1/2
    eta = log(mu)
  }
  
  beta = matrix(0, p, 1)
  beta.v2 = matrix(0, p, 1)
  chi = 2
  
  i = 1
  while (1)
  {
    beta.old = beta
    
    # Get weights
    if (model == "binomial")
    {
      #deta = 1/mu/(1-mu)
      #z = eta + (y-mu)*deta
      #w2 = mu*(1-mu)
      w2 = pgdraw.moments(1, eta)$mu
      z = (y-1/2)/w2
    }
    else if (model == "poisson")
    {
      deta = 1/mu
      z = eta + (y-mu)*deta
      w2 = mu
    }
    
    # Update
    rv = linreg.w2.post(X, z, 1/w2, d, approx = F, sample = F)
    beta = rv$beta.mu
    
    eta = X %*% beta
    
    # Conditional mean
    if (model == "binomial")
    {
      mu = 1/(1+exp(-eta))
    }
    else if (model == "poisson")
    {
      mu = exp(eta)
    }
    
    i = i+1
    
    QQ=( sum(abs(beta-beta.old))/(1+sum(abs(beta))) )
    if (is.nan(QQ))
    {
      QQ
    }
    
    #
    e = 1e-3
    if (sum(abs(beta-beta.old))/(1+sum(abs(beta))) < e)
    {
      break
    }
  }
  
  return(list(beta.mu=beta, beta.v2=rv$beta.v2, i = i))
}

# ============================================================================================================================
# function to standardise columns of X to have mean zero and unit length
my.standardise <- function(X)
{
  n = nrow(X)
  p = ncol(X)
  
  # 
  r       = list()
  r$X     = X
  if (p > 1)
  {
    r$mean.X = colMeans(X)
  } else
  {
    r$mean.X = mean(X)
  }
  r$std.X  = sqrt( t(apply(X,2,stats::sd))^2*(n-1)/n )
  
  # Perform the standardisation
  if (p == 1)
  {
    r$X <- as.matrix(apply(X,1,function(x)(x - r$mean.X)))
    r$X <- as.matrix(apply(r$X,1,function(x)(x / r$std.X)))
  } else
  {
    r$X <- t(as.matrix(apply(X,1,function(x)(x - r$mean.X))))
    r$X <- t(as.matrix(apply(r$X,1,function(x)(x / r$std.X))))
  }
  
  return(r)
}

# ============================================================================================================================
# Find starting point
find.start <- function(rv, y, X)
{
  n.samples = ncol(rv.br$beta)
  p = nrow(rv.br$beta)
  n = length(y)
  
  alpha = 0.9
  
  Eb2 = matrix(rv$beta[,1]^2, p, 1)
  E.RSS = sum( (y - X%*%rv$beta[,1] - rv$beta0[1])^2 )
  
  Eb2.store = matrix(0,p,n.samples)
  E.RSS.store = matrix(0,1,n.samples)
  
  tau2 = rv$tau2[1]
  sigma2 = rv$sigma2[1]
  
  L = matrix(0,n.samples,1)
  
  for (i in 2:n.samples)
  {
    Eb2 = Eb2*alpha + rv.br$beta[,i]^2*(1-alpha)
    E.RSS = E.RSS*alpha + sum( (y - X%*%rv$beta[,i] - rv$beta0[i])^2 )*(1-alpha)
    
    Eb2.store[,i] = Eb2
    E.RSS.store[i] = E.RSS
    
    #tau2 = tau2*alpha + rv$tau2[i]*(1-alpha)
    #sigma2 = sigma2*alpha + rv$sigma2[i]*(1-alpha)
    
    min.par = exp((optim( log(c(tau2,sigma2)), function(x){ tau2.sigma2.f(x, Eb2, n, E.RSS) } ))$par)
    tau2 = min.par[1]
    sigma2 = min.par[2]
    lambda2 <- est.lambda2(Eb2, sigma2, tau2)
    
    # Neg-log-posterior
    L[i] = E.neg.log.posterior(Eb2, E.RSS, lambda2, tau2, sigma2, n, p)
  }
  
  return(list(Eb2 = Eb2.store[,100:n.samples], E.RSS = E.RSS.store[100:n.samples], L = L[100:n.samples]))
}

# ============================================================================================================================
# Find starting point
find.start.glm <- function(rv, X, y)
{
  n.samples = ncol(rv.br$beta)
  p = nrow(rv.br$beta)
  n = length(y)
  
  alpha = 0
  
  # Starting point
  Eb2 = matrix(rv$beta[,1]^2, p, 1)
  E.NLL = neg.log.likelihood(X,y,rv$beta[,1],rv$beta0[1],rv$model)
  
  Eb2.store = matrix(0,p,n.samples)
  E.NLL.store = matrix(0,1,n.samples)
  
  tau2 = rv$tau2[1]
  sigma2 = rv$sigma2[1]
  
  L = matrix(0,n.samples,1)
  for (i in 2:n.samples)
  {
    Eb2 = Eb2*alpha + rv.br$beta[,i]^2*(1-alpha)
    E.NLL = E.NLL*alpha + (1-alpha)*neg.log.likelihood(X,y,rv$beta[,i],rv$beta0[i],rv$model)
    
    Eb2.store[,i] = Eb2
    E.NLL.store[i] = E.NLL
    
    #tau2 = tau2*alpha + rv$tau2[i]*(1-alpha)
    #sigma2 = sigma2*alpha + rv$sigma2[i]*(1-alpha)
    
    tau2 = find.tau2(Eb2,sigma2=1)
    lambda2 <- est.lambda2(Eb2, sigma2=1, tau2)
    
    # Neg-log-posterior
    L[i] = E.neg.log.posterior.glm(E.NLL, Eb2, lambda2, tau2)
  }
  
  return(list(Eb2 = Eb2.store[,100:n.samples], E.NLL = E.NLL.store[100:n.samples], L = L[100:n.samples]))
}

# ============================================================================================================================
# Expected negative log-posterior for linear regression model
E.neg.log.posterior <- function(Eb2, E.RSS, lambda2, tau2, sigma2, n, p)
{
  # Neg-log-posterior
  L = (n+p)/2*log(sigma2) + E.RSS/2/sigma2 + p/2*log(tau2) + (1/2)*sum(log(lambda2)) + 1/sigma2/tau2*sum(Eb2/lambda2)
  L + sum(log(1 + lambda2)) + log(1+tau2)
}

# ============================================================================================================================
# Expected negative log-posterior for GLM
E.neg.log.posterior.glm <- function(L, Eb2, lambda2, tau2)
{
  p = length(Eb2)
  
  #
  L = L + p/2*log(tau2) + (1/2)*sum(log(lambda2)) + 1/tau2*sum(Eb2/lambda2)
  L + sum(log(1 + lambda2)) + log(1+tau2)
}

# ============================================================================================================================
# Negative log-likelihood
neg.log.likelihood <- function(X, y, beta, beta0, model, sigma2 = NULL)
{
  n = nrow(X)
  eta = X %*% beta + beta0
  
  if (model == "gaussian")
  {
    L = (n/2)*log(sigma2) + (y-eta)^2/2/sigma2
  }
  else if (model == "binomial")
  {
    mu = 1/(1+exp(-eta))
    L = -sum(log(mu[y==1])) - sum(log(1-mu[y==0]))
  }
  else if (model == "poisson")
  {
    mu = exp(eta)
    L = sum( -y*eta + mu + lgamma(y+1) )
  }
  L
}

em.predict <- function(object, newdata)
{
  # Build the fully specified formula using the covariates that were fitted
  f <- stats::as.formula(paste("~",paste(attr(object$terms,"term.labels"),collapse="+")))
  
  # Extract the design matrix
  X = stats::model.matrix(f, data=newdata)
  X = as.matrix(X[,-1])
  X = X[,object$I.keep]
  n = nrow(X)
  p = ncol(X)
  
  # Make predictions
  yp = X %*% object$beta + as.numeric(object$beta0)
  
  # If GLM
  if (object$model == "binomial")
  {
    yp = 1/(1+exp(-yp))
  }
  else if (object$model == "poisson")
  {
    yp = exp(yp)
  }
  
  yp
}

# ==============================================================================
# Horseshoe-like 
# ==============================================================================

HS.like.EM <- function(formula, data, n.iter = 1000, standardize = T, a = 1e5, tol = 1e-4)
{
  r = list()
  
  # -------------------------------------------------------------------    
  # Process and set up the data from the model formula
  r$terms <- stats::terms(x = formula, data = data)
  
  mf = stats::model.frame(formula = formula, data = data)
  r$target.var = names(mf)[1]
  
  y = mf[,1]
  X = stats::model.matrix(formula, data=data)
  
  # Convert to a numeric matrix and drop the target variable
  X = as.matrix(X)
  X = X[,-1,drop=FALSE]
  X0 = X
  
  I.keep = rep(T, ncol(X))
  if(standardize)
  {
    std.X = my.standardise(X)
    X = std.X$X
    
    I.keep = std.X$std.X!=0
    X      = std.X$X[,I.keep]
  }
  
  #initialize parameters
  n = nrow(X)
  p = ncol(X)
  theta = rep(1, p)
  
  # ------
  xty = crossprod(X,y)
  xtx = crossprod(X,X)
  
  theta = solve(xtx + diag(1,p),xty)
  
  condition = T
  
  counter = 0
  while (condition) 
  {
    theta_old = theta
    
    u     = (1/(2*pi*sqrt(a))) * ((a/theta^2) - (a/(theta^2 + a)))
    a     = ((a^(3/2))/(p*pi)) * sum(1/(a + theta^2))
    
    u = pmax(u, 1e-10)
    
    #d = pmin(a/(2*u),1e10)
    #d = pmax(a/(2*u),1e-30)
    d = a/(2*u)
    
    #D     = diag(as.vector(a/(2*u)))
    D     = diag(as.vector(d))
    X0    = X
    for (i in 1:p)
    {
      X0[,i] = X0[,i]*d[i]
    }
    W     = X0 %*% t(X) + diag(n)
    
    #W     = X%*%D%*%t(X) + diag(n)
    
    L     = chol(W)
    vv    = forwardsolve(t(L), y)
    w     = backsolve(L, vv)
    
    theta   = D%*%t(X)%*% w
    
    d = abs(theta_old - theta)
    #difference_small = (max(d/abs(theta_old), na.rm = T) < tol)
    
    
    difference_small = (sum(abs(theta_old-theta)) / (1+sum(abs(theta))) < tol) # might end too early
    # difference_small = difference_small | (sum(abs(a_old-a)) / (1+sum(abs(a))) < (tol^2))
    
    #cat(sum(abs(theta_old-theta)) / (1+sum(abs(theta))),",")
    
    
    # update condition
    if (difference_small && counter > 20)
    {
      condition = FALSE
    }
    if (counter == 1000)
    {
      break
    }
    
    counter=counter+1
  }
  
  # Done
  r$a = a
  r$beta0 = mean(y)
  r$beta = theta
  r$converge = counter
  
  if(standardize)
  {
    r$beta  <- r$beta / (std.X$std.X[I.keep])
    r$beta0 <- r$beta0 - std.X$mean.X[I.keep] %*% r$beta
    
    r$std.X = std.X
  }
  r$model = "gaussian"
  r$I.keep = I.keep
  
  return(r)
}



# ==============================================================================
# Horseshoe-like 
# ==============================================================================

HS.like.EM.2 <- function(formula, data, n.iter = 1000, standardize = T, a = 1e5, tol = 1e-5)
{
  r = list()
  
  # -------------------------------------------------------------------    
  # Process and set up the data from the model formula
  r$terms <- stats::terms(x = formula, data = data)
  
  mf = stats::model.frame(formula = formula, data = data)
  r$target.var = names(mf)[1]
  
  y = mf[,1]
  X = stats::model.matrix(formula, data=data)
  
  # Convert to a numeric matrix and drop the target variable
  X = as.matrix(X)
  X = X[,-1,drop=FALSE]
  X0 = X
  
  if(standardize)
  {
    std.X = my.standardise(X)
    X = std.X$X
  }
  
  #initialize parameters
  n = nrow(X)
  p = ncol(X)
  theta = rep(1, p)
  
  # ------
  xty = crossprod(X,y)
  xtx = crossprod(X,X)
  
  condition = T
  
  counter = 0
  while (condition) 
  {
    theta_old = theta
    
    u     = (1/(2*pi*sqrt(a))) * ((a/theta^2) - (a/(theta^2 + a)))
    a     = pmax( ((a^(3/2))/(p*pi)) * sum(1/(a + theta^2)), 1e-20)
    
    #d = pmin(a/(2*u),1e10)
    
    d = pmax(2*u/a, 1e-5)
    d = pmin(d, 1e10)
    
    theta = solve(xtx+diag(as.vector(d)),xty)
    
    #D     = diag(as.vector(a/(2*u)))
    #D     = diag(as.vector(d))
    #W     = X%*%D%*%t(X) + diag(n)
    
    #L     = chol(W)
    #vv    = forwardsolve(t(L), y)
    #w     = backsolve(L, vv)
    
    #theta   = D%*%t(X)%*% w
    
    d = abs(theta_old - theta)
    difference_small = (max(d/abs(theta_old), na.rm = T) < tol)
    
    
    difference_small = (sum(abs(theta_old-theta)) / (1+sum(abs(theta))) < tol) # might end too early
    #difference_small = difference_small | (sum(abs(a_old-a)) / (1+sum(abs(a))) < (tol^2))
    
    #cat(sum(abs(theta_old-theta)) / (1+sum(abs(theta))),",")
    
    
    # update condition
    if (counter == 500 || difference_small)
    {
      condition = FALSE
    }
    
    counter=counter+1
  }
  
  t = 1/sqrt(n)/5
  theta[abs(theta)<t] = 0
  
  # Done
  r$a = a
  r$beta0 = mean(y)
  r$beta = theta
  r$converge = counter
  
  if(standardize)
  {
    r$beta  <- r$beta / t(std.X$std.X)
    r$beta0 <- r$beta0 - std.X$mean.X %*% r$beta
    
    r$std.X = std.X
  }
  r$model = "gaussian"
  
  return(r)
}

tr <- function(X)
{
  sum(diag(X))
}

linreg.post.glmnet <- function(X, y, lambda2, tau2, sigma2)
{
  n = nrow(X)
  p = ncol(X)
  
  s = sum(1/lambda2)
  bhat = as.matrix(coefficients(glmnet(X,y,"gaussian",
                                       alpha=0,lambda=s/(p*tau2*n),
                                       penalty.factor=1/lambda2,
                                       thresh = 1e-10)))
  
  beta.var = 1/(n/sigma2 + 1/(lambda2*tau2*sigma2))
  
  E.RSS = sum( (y - X %*% bhat[2:(p+1)] - bhat[1])^2 ) + sum( n*beta.var )
  
  return(list(beta.mu = bhat[2:(p+1)], beta0 = bhat[1], beta.v2 = beta.var, E.RSS = E.RSS))
}

glmreg.post.glmnet <- function(X, y, model, lambda2, tau2)
{
  n = nrow(X)
  p = ncol(X)
  
  # Fit
  s = sum(1/lambda2)
  bhat = as.matrix(coefficients(glmnet(X,y,model,
                                       alpha=0,lambda=s/(p*tau2*n),
                                       penalty.factor=1/lambda2,
                                       thresh = 1e-10)))
  b0 = bhat[1]
  beta = bhat[2:(p+1)]
  
  
  eta = X %*% beta + b0
  mu = 1/(1+exp(-eta))
  rv.omega2 = compute.omega2(model,mu,eta,y,T)
  omega2 = as.vector(rv.omega2$omega2)
  
  #X0 = X
  XtX = rep(0, p)
  for (i in 1:p)
  {
    XtX[i] = sum(X[,i]^2*omega2)
  }
  
  # Approximate variances  
  beta.var = 1/(XtX + 1/(lambda2*tau2))
  
  return(list(beta.mu = bhat, beta.v2 = c(0,beta.var)))
}
