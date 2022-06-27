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
# Estimate tau2 & lambda2, depending on prior
# ==============================================================================

est.tau2.lambda2.sigma2 <- function(prior, model, Eb2, n = NULL, E.RSS = NULL, Enu.inv = NULL, Exi.inv = NULL, log.tau2.max = 20, tau2=1, sigma2=1)
{
  # ----------------------------
  # Gaussian regression
  if (model == "gaussian")
  {
    # Find tau2/lambda2 from Eb2 depending on the prior
    rv = find.tau2.sigma2(Eb2,n,E.RSS, log.tau2.max=log.tau2.max,tau2=tau2,sigma2=sigma2, prior = prior)
    tau2 = rv$tau2
    sigma2 = rv$sigma2
    lambda2 = est.lambda2(Eb2, sigma2, tau2, prior)
  }
  
  # ----------------------------
  # GLMs -- no sigma2  
  else
  {
    
    tau2 = find.tau2(Eb2, sigma2=1, log.tau2.max=log.tau2.max, prior = prior)
    lambda2 = est.lambda2(Eb2, sigma2=1, tau2, prior = prior)
    
    sigma2 = NULL
  }
  
  return(list(lambda2=lambda2,tau2=tau2,sigma2=sigma2))
}


find.tau2.sigma2 <- function(Eb2,n,E.RSS,log.tau2.max=20,tau2=1,sigma2=1, prior)
{
  p = length(Eb2)
  
  tau2 = exp((optimise(function(log.tau2){tau2.sigma2.f(log.tau2,Eb2,n,E.RSS,prior)}, c(-10,log.tau2.max)))$minimum)
  sigma2 = E.RSS/n
  
  return(list(tau2=tau2,sigma2=sigma2))  
}



tau2.sigma2.f <- function(x, Eb2, n, Erss, prior)
{
  p = length(Eb2)
  tau2 = exp(x)
  sigma2 = Erss/n
  lambda2 = est.lambda2(Eb2,sigma2,tau2,prior)
  
  f = p/2*log(tau2) + 1/2*sum(log(lambda2)) + 1/2/sigma2/tau2*sum(Eb2/lambda2) # beta 
  f = f + (n+p)/2*log(sigma2) + Erss/2/sigma2 + log(sigma2) # data model
  
  if(prior == "hs"){
    f = f + log(1+tau2) + sum(log(1+lambda2))
  }else if (prior == "hs.lambda2"){
    f = f + sum( (1/2)*log(lambda2) + log(lambda2+1) )
    f = f + (1/2)*log(tau2) + log(1+tau2)
  }else if (prior == "lasso"){
    f = f + sum(lambda2/2) + 2*log(tau2) + 1/tau2
  }else if (prior == "lasso.lambda"){
    f = f + sum(lambda2/2) - sum(log(lambda2)/2) + 2*log(tau2) + 1/tau2
  }
  
  f
}

est.lambda2 <- function(Eb2,sigma2,tau2,prior)
{
  W = Eb2/2/sigma2/tau2
  
  if(prior == "hs"){
    lambda2 = (sqrt(-1 + 2*W + sqrt(1 + 20*W + 4*W^2))/sqrt(6))^2
  }else if (prior == "hs.lambda2"){
    lambda2 = (sqrt(W^2 + 6*W + 1) + W - 1) / 4
  }else if (prior == "lasso"){
    lambda2 = (sqrt(1 + 8*W)-1)/2
  }else if (prior == "lasso.lambda"){
    lambda2 = sqrt(2)*sqrt(W)
  }
  
  return(pmax(lambda2,1e-50))
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
# Expected negative log-posterior for linear regression model
E.neg.log.posterior <- function(Eb2, E.RSS, lambda2, tau2, sigma2, n, p)
{
  # Neg-log-posterior
  L = (n+p)/2*log(sigma2) + E.RSS/2/sigma2 + p/2*log(tau2) + (1/2)*sum(log(lambda2)) + 1/sigma2/tau2*sum(Eb2/lambda2)
  L + sum(log(1 + lambda2)) + log(1+tau2)
}


