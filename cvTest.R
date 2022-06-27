# ==============================================================================
# Simple CV split with added noise variables for real data tests
# ==============================================================================
source("sparseEMregression.R")

EM.real.data.CV <- function(data, target, n.train, p.extra=0, 
                            use.interactions=F, use.logs=F, use.squares=F, use.cubics=F, seed = NULL,
                            rho = 0)
{
  # Set the seed if passed
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  
  # Get the CV indices
  n = nrow(data)
  ix = sample(n)
  
  tr.ix = ix[1:n.train]
  tst.ix = ix[(n.train+1):n]
  
  # Extra noise variables, if requested
  if (p.extra > 0)
  {
    #X.extra = matrix(rnorm(n*p.extra),n,p.extra)
    C = toeplitz(rho^(0:(p.extra-1)))
    X.extra = mvrnorm(n, rep(0,p.extra), C)    
    
    var.names = c()
    for (j in 1:p.extra)
    {
      var.names = c(var.names, sprintf("X.%d",j))
    }
    colnames(X.extra) = var.names
    
    data = cbind(data, X.extra)
  }
  
  my.form = my.make.formula(target,data,use.interactions=use.interactions,use.logs=use.logs,use.squares=use.squares,use.cubics=use.cubics)
  my.term = stats::terms(x = my.form, data = data)
  
  return(list(train.data=data[tr.ix,], test.data=data[tst.ix,], my.form = my.form, tr.ix=tr.ix, tst.ix=tst.ix, my.term = my.term))
}

# ==============================================================================
# Prediction stats for binary models
# ==============================================================================

my.pred.stats <- function(prob, target, display = TRUE)
{
  rv = list()
  
  classes = levels(target)
  
  # Convert probabilities to best guesses at classes
  pred = factor(prob > 1/2, c(F,T), classes)
  
  # Compute statistics
  T = table(pred,target)
  roc.obj = roc(response=as.numeric(target)-1, as.vector(prob), quiet=TRUE)
  rv$ca   = mean(pred==target)
  rv$sens = T[2,2]/(T[1,2]+T[2,2])
  rv$spec = T[1,1]/(T[1,1]+T[2,1])
  rv$auc  = as.numeric(roc.obj$auc)
  
  # Prob is probability of success, so if the target is not a success, flip the probability
  # to get probability of failure
  prob[target==classes[1]] = 1 - prob[target==classes[1]]
  # Also make sure we never get exactly zero or one for probabilities due to numerical rounding
  prob = (prob+1e-10)/(1+2e-10)
  
  rv$log.loss = -sum(log(prob))
  
  # Display, if requested    
  if (display == TRUE)
  {
    cat("---------------------------------------------------------------------------\n")
    cat("Performance statistics:\n")
    cat("\n")
    cat("Confusion matrix:\n\n")
    print(T)
    cat("\n")
    cat("Classification accuracy =", rv$ca, "\n")
    cat("Sensitivity             =", rv$sens, "\n")
    cat("Specificity             =", rv$spec, "\n")
    cat("Area-under-curve        =", rv$auc, "\n")
    
    cat("Logarithmic loss        =", rv$log.loss, "\n")
    cat("\n")
    
    #plot(roc.obj)
    
    cat("---------------------------------------------------------------------------\n")
  }
  else 
  {
    #
    return(rv)
  }
}

# ==============================================================================
# Binary real data test script
# ==============================================================================

binary.real.data.test <- function(data, target, n.train, p.extra, n.iter, methods, log.tau2.max, use.approx = F, seed = NULL)
{
  # CV test
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  
  n.methods = length(methods)
  
  CA       = matrix(nrow = n.iter, ncol = n.methods)
  colnames(CA) = methods
  log.loss = CA
  AUC      = CA
  run.time = CA
  n.coefs  = CA
  
  for (j in 1:n.iter)
  {
    cat(j,", ")
    rv.cv = EM.real.data.CV(data, target, n.train, p.extra, use.interactions = T, use.logs = T, use.squares = T, use.cubics = T)
    
    # GLMNET
    for (i in 1:n.methods)
    {
      # GLMNET
      if (methods[i] == "glmnet")
      {
        s = system.time({rv.glmnet = cv.glmnet.f(rv.cv$my.form, data=rv.cv$train.data, family="binomial", alpha = 1)})
        rv.ps = my.pred.stats(predict.glmnet.f(rv.glmnet, rv.cv$test.data, type="response"), rv.cv$test.data[,target], display=F)
        
        run.time[j,i] = s[[3]]
        log.loss[j,i] = rv.ps$log.loss
        AUC[j,i] = rv.ps$auc
        CA[j,i] = rv.ps$ca
        n.coefs[j,i] = sum(coefficients(rv.glmnet)!=0)      
      }
      else
      {
        if (i == 3)
        {
          Eb2.start = rv.em$rv.em$Eb2
          # EM
          s = system.time({rv.em = beta.EM(rv.cv$my.form, data=rv.cv$train.data, "logistic", prior = methods[i], approx = use.approx, 
                                           e = 1e-4, stochastic = F, PG = T, 
                                           log.tau2.max = log.tau2.max[i],
                                           Eb2.start = Eb2.start)})
        }
        else
        {
          s = system.time({rv.em = beta.EM(rv.cv$my.form, data=rv.cv$train.data, "logistic", prior = methods[i], approx = use.approx, 
                                           e = 1e-4, stochastic = F, PG = T, 
                                           log.tau2.max = log.tau2.max[i])})
        }
        
        rv.ps = my.pred.stats(em.predict(rv.em, rv.cv$test.data), rv.cv$test.data[,target], display=F)
        
        run.time[j,i] = s[[3]]
        log.loss[j,i] = rv.ps$log.loss
        AUC[j,i] = rv.ps$auc
        CA[j,i] = rv.ps$ca
        n.coefs[j,i] = sum(rv.em$beta!=0)        
      }
    }
  }
  
  #
  scores = matrix(nrow = 5, ncol = n.methods)
  colnames(scores) = methods
  rownames(scores) = c("CA","AUC","log.loss","run.time","n.coefs")
  
  scores[1,] = colMeans(CA)
  scores[2,] = colMeans(AUC)
  scores[3,] = colMeans(log.loss)
  scores[4,] = colMeans(run.time)
  scores[5,] = colMeans(n.coefs)
  
  cat("\n")
  
  return(scores)
}


# ==============================================================================
# Continuous real data test script
# ==============================================================================

continuous.real.data.test <- function(data, target, n.train, p.extra, n.iter, 
                                      methods, log.tau2.max, e = 1e-4, use.approx = NULL, seed = NULL, 
                                      use.interactions = T, use.logs = T, use.squares = T, use.cubics = T, 
                                      use.glmnet = NULL, rho = 0)
{
  # CV test
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  
  n.methods = length(methods)
  
  if (is.null(use.approx))
  {
    use.approx = rep(F, n.methods)
  }
  if (is.null(use.glmnet))
  {
    use.glmnet = rep(F, n.methods)
  }
  
  MSE      = matrix(nrow = n.iter, ncol = n.methods)
  colnames(MSE) = methods
  run.time = MSE
  n.coefs  = MSE
  corr.zero = MSE
  
  rv.cv.pre = list()
  for (j in 1:n.iter)
  {
    rv.cv.pre[[j]] = EM.real.data.CV(data, target, n.train, p.extra, 
                                     use.interactions = use.interactions, use.logs = use.logs, 
                                     use.squares = use.squares, use.cubics = use.cubics,
                                     rho = rho)
  }
  
  for (j in 1:n.iter)
  {
    cat(j,", ")
    
    rv.cv = rv.cv.pre[[j]]
    
    noisy_coef = grep("X.", attr(rv.cv$my.term, "term.labels"))
    
    # GLMNET
    for (i in 1:n.methods)
    {
      print(methods[i])
      # GLMNET
      if (methods[i] == "glmnet")
      {
        s = system.time({rv.glmnet = cv.glmnet.f(rv.cv$my.form, data=rv.cv$train.data, family="gaussian", alpha = 1)})
        mse = mean( (predict.glmnet.f(rv.glmnet, rv.cv$test.data, type="response") - rv.cv$test.data[,target])^2 )
        
        run.time[j,i] = s[[3]]
        MSE[j,i] = mse
        n.coefs[j,i] = sum(coefficients(rv.glmnet)[-1]!=0)
        corr.zero[j,i] = sum(coefficients(rv.glmnet)[noisy_coef]==0)/length(noisy_coef)
      }
      # GLMNET-Ridge
      else if (methods[i] == "glmnet.ridge")
      {
        s = system.time({rv.glmnet = cv.glmnet.f(rv.cv$my.form, data=rv.cv$train.data, family="gaussian", alpha = 0)})
        mse = mean( (predict.glmnet.f(rv.glmnet, rv.cv$test.data, type="response") - rv.cv$test.data[,target])^2 )
        
        run.time[j,i] = s[[3]]
        MSE[j,i] = mse
        n.coefs[j,i] = sum(coefficients(rv.glmnet)[-1]!=0)
        corr.zero[j,i] = sum(coefficients(rv.glmnet)[noisy_coef]==0)/length(noisy_coef)
      }
      else if (methods[i] == "MCP" || methods[i] == "SCAD")
      {
        s = system.time({rv.ncv = cv.ncvreg.f(rv.cv$my.form, data=rv.cv$train.data, family="gaussian", penalty = methods[i])})
        mse = mean( (predict.ncvreg.f(rv.ncv, data=rv.cv$test.data) - rv.cv$test.data[,target])^2 )
        
        run.time[j,i] = s[[3]]
        MSE[j,i] = mse
        n.coefs[j,i] = sum(coefficients(rv.ncv)[-1]!=0)
        corr.zero[j,i] = sum(coefficients(rv.ncv)[noisy_coef]==0)/length(noisy_coef)
      }
      else if (methods[i] == "hs.like")
      {
        s = system.time({rv.hs.like = HS.like.EM(rv.cv$my.form, data=rv.cv$train.data)})
        mse = mean( (em.predict(rv.hs.like, rv.cv$test.data) - rv.cv$test.data[,target])^2 )
        
        run.time[j,i] = s[[3]]
        MSE[j,i] = mse
        n.coefs[j,i] = sum(rv.hs.like$beta!=0)
        corr.zero[j,i] = sum(rv.hs.like$beta[noisy_coef]==0)/length(noisy_coef)
      }
      else if (methods[i] == "bayesreg")
      {
        s = system.time({rv.br = bayesreg(rv.cv$my.form, data=rv.cv$train.data, model="gaussian", prior="hs", n.samples=5e3)})
        
        mse = mean( (predict(rv.br, rv.cv$test.data) - rv.cv$test.data[,target])^2 )
        
        run.time[j,i] = s[[3]]
        MSE[j,i] = mse
        n.coefs[j,i] = sum(rv.br$mu.beta!=0)
        corr.zero[j,i] = sum(rv.br$mu.beta[noisy_coef]==0)/length(noisy_coef)
      }
      else
      {
        # EM
        s = system.time({rv.em = beta.EM(rv.cv$my.form, data=rv.cv$train.data, "gaussian", 
                                         prior = methods[i], approx = use.approx[i], e = 1e-4, 
                                         stochastic = F, PG = T, log.tau2.max = log.tau2.max[i],
                                         use.glmnet = use.glmnet[i])})
        
        mse = mean( (em.predict(rv.em, rv.cv$test.data) - rv.cv$test.data[,target])^2 )
        
        run.time[j,i] = s[[3]]
        MSE[j,i] = mse
        n.coefs[j,i] = sum(rv.em$beta!=0)
        corr.zero[j,i] = sum(rv.em$beta[noisy_coef]==0)/length(noisy_coef)
      }
    }
  }
  
  #
  scores = matrix(nrow = 4, ncol = n.methods)
  colnames(scores) = methods
  rownames(scores) = c("MSE","run.time","n.coefs", "corr.zero")
  
  scores[1,] = colMeans(MSE)
  scores[2,] = colMeans(run.time)
  scores[3,] = colMeans(n.coefs)
  scores[4,] = colMeans(corr.zero)
  
  
  cat("\n")
  
  return(list(p=length(rv.em$beta),scores=scores,MSE=MSE,run.time=run.time,n.coefs=n.coefs, corr.zero = corr.zero))
}


# ==============================================================================
# Continuous synthetic data test script
# ==============================================================================

continuous.synthetic.data.test <- function(beta, cov.X, SNR, n.train, n.iter, 
                                           methods, log.tau2.max, e = 1e-4, use.approx = NULL, seed = NULL, 
                                           use.glmnet = NULL)
{
  beta = as.vector(beta)
  
  # Setup
  if (!is.null(seed))
  {
    set.seed(seed)
  }
  
  n.methods = length(methods)
  
  if (is.null(use.approx))
  {
    use.approx = rep(F, n.methods)
  }
  if (is.null(use.glmnet))
  {
    use.glmnet = rep(F, n.methods)
  }
  
  # Determine noise level from SNR
  # SNR = sqrt( E[FSS/n] / sigma2) => 
  FSS = t(beta) %*% cov.X %*% beta
  sigma2 = FSS/SNR^2
  
  # Results
  MSE      = matrix(nrow = n.iter, ncol = n.methods)
  colnames(MSE) = methods
  SSE       = MSE
  run.time  = MSE
  n.coefs   = MSE
  TZ = MSE
  FZ = MSE
  TNZ = MSE
  FNZ = MSE
  
  rv.cv.pre = list()
  
  p    = length(beta)
  
  X.tr = list()
  y.tr = list()
  for (j in 1:n.iter)
  {
    #rv.cv.pre[[j]] = EM.real.data.CV(data, target, n.train, p.extra, 
    #                                 use.interactions = use.interactions, use.logs = use.logs, 
    #                                 use.squares = use.squares, use.cubics = use.cubics,
    #                                 rho = rho)
    X.tr[[j]] = mvrnorm(n.train, rep(0,p), cov.X)
    y.tr[[j]] = X.tr[[j]] %*% beta + rnorm(n.train, 0, sqrt(sigma2))
  }
  
  for (j in 1:n.iter)
  {
    cat(j,", ")
    
    df.tr = data.frame(y=y.tr[[j]],X=X.tr[[j]])
    
    # GLMNET
    for (i in 1:n.methods)
    {
      # GLMNET
      if (methods[i] == "glmnet")
      {
        s = system.time({rv.glmnet = cv.glmnet.f(y ~ ., data=df.tr, family="gaussian", alpha = 1)})
        b.hat = coefficients(rv.glmnet)[2:(p+1)]
      }
      # GLMNET-Ridge
      else if (methods[i] == "glmnet.ridge")
      {
        s = system.time({rv.glmnet = cv.glmnet.f(y ~ ., data=df.tr, family="gaussian", alpha = 0)})
        b.hat = coefficients(rv.glmnet)[2:(p+1)]
      }
      else if (methods[i] == "MCP" || methods[i] == "SCAD")
      {
        s = system.time({rv.ncv = cv.ncvreg.f(y ~ ., data=df.tr, family="gaussian", penalty = methods[i])})
        b.hat = coefficients(rv.ncv)[2:(p+1)]
      }
      else if (methods[i] == "hs.like")
      {
        s = system.time({rv.hs.like = HS.like.EM(y ~ ., data=df.tr)})
        b.hat = rv.hs.like$beta
      }
      else if (methods[i] == "bayesreg")
      {
        s = system.time({rv.br = bayesreg(y ~ ., data=df.tr, model="gaussian", prior="hs", n.samples=5e3)})
        b.hat = rv.br$mu.beta
      }
      else
      {
        # EM
        s = system.time({rv.em = beta.EM(y ~ ., data=df.tr, model = "gaussian", 
                                         prior = methods[i], approx = use.approx[i], e = e, 
                                         stochastic = F, PG = T, log.tau2.max = log.tau2.max[i],
                                         use.glmnet = use.glmnet[i])})
        b.hat = rv.em$beta        
      }
      
      # Compute scores
      run.time[j,i] = s[[3]]
      MSE[j,i] = t(b.hat - beta) %*% cov.X %*% (b.hat - beta)
      SSE[j,i] = sum( (b.hat - beta)^2 )
      
      #n.z.corr[j,i] = sum( (b.hat == 0) & (beta == 0) )
      #n.nz.corr[j,i] = sum( (b.hat != 0) & (beta != 0) )
      
      TZ[j,i] = sum( (b.hat == 0) & (beta == 0) )
      FZ[j,i] = sum( (b.hat == 0) & (beta != 0) )
      FNZ[j,i] = sum( (b.hat != 0) & (beta == 0) )
      TNZ[j,i] = sum( (b.hat != 0) & (beta != 0) )
      
      n.coefs[j,i] = sum(b.hat != 0)      
    }
  }
  
  #
  scores = matrix(nrow = 8, ncol = n.methods)
  colnames(scores) = methods
  rownames(scores) = c("MSE","SSE","run.time","n.coefs", "TZ", "FZ", "TNZ", "FNZ")
  
  scores[1,] = colMeans(MSE)
  scores[2,] = colMeans(SSE)
  scores[3,] = colMeans(run.time)
  scores[4,] = colMeans(n.coefs)
  scores[5,] = colMeans(TZ)
  scores[6,] = colMeans(FZ)
  scores[7,] = colMeans(TNZ)
  scores[8,] = colMeans(FNZ)
  
  cat("\n")
  
  return(list(p=length(rv.em$beta),scores=scores,MSE=MSE,SSE=SSE,TZ=TZ,FZ=FZ,TNZ=TNZ,FNZ=FNZ,run.time=run.time,n.coefs=n.coefs))
}



compute_se <- function(my_list){
  
  dimension = dim(my_list$scores)
  n.iter = nrow(my_list[[3]])
  
  final      = matrix(nrow = dimension[1], ncol = dimension[2])
  rownames(final) = names(my_list)[-(1:2)]
  colnames(final) = colnames(my_list[[3]])
  
  for(i in 3:length(my_list)){
    final[i-2,] = apply(my_list[[i]],2,sd)/sqrt(n.iter)
  }
  
  return(final)
  
}
