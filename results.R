source("cvTest.R")

####################################################################################
# DIABETES
####################################################################################

diabetes = read.csv("diabetes.train.csv")
diabetes = rbind(diabetes,read.csv("diabetes.test.csv"))

diabetes.scores = continuous.real.data.test(diabetes, "Y", 100, 15, 50, c("glmnet","hs.lambda2","hs.lambda2","glmnet.ridge","MCP","SCAD", "hs.like"), 
                                   log.tau2.max = c(0,0,0,0,0), use.approx = c(F,F,T,F,T,F,F), seed = 1000, use.glmnet = c(F,F,F,F,F), 
                                   rho = 0.8)


####################################################################################
# BOSTON HOUSING
####################################################################################

bh.data = read.csv("boshouse.csv",header=F)

colnames(bh.data) <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","y")

bh.scores = continuous.real.data.test(bh.data, "y", 100, 15, 50, c("glmnet","hs.lambda2","hs.lambda2","glmnet.ridge","MCP","SCAD", "hs.like"), 
                                   log.tau2.max = c(0,0,0,0,0), use.approx = c(F,F,T,F,T,F,F), seed = 1000, use.glmnet = c(F,F,F,F,F), 
                                   rho = 0.8)


####################################################################################
# CONCRETE
####################################################################################

concrete = read.csv("Concrete_Data.csv",header=T)

concrete.scores = continuous.real.data.test(concrete, "Strength", 100, 15, 50, c("glmnet","hs.lambda2","hs.lambda2","glmnet.ridge","MCP","SCAD", "hs.like"), 
                                   log.tau2.max = c(0,0,0,0,0), use.approx = c(F,F,T,F,T,F,F), seed = 1000, use.glmnet = c(F,F,F,F,F), 
                                   use.interactions = T, use.logs = T, use.squares = T, use.cubics = T,
                                   rho = 0.8)


####################################################################################
# EYE DATA
####################################################################################

eye.data = read.csv("eyedata.csv")

eye.data[,1] <- NULL


eye.scores = continuous.real.data.test(eye.data, "y", 100, 0, 50, c("glmnet","hs.lambda2","hs.lambda2","glmnet.ridge","MCP","SCAD", "hs.like"), 
                                   log.tau2.max = c(0,0,0,0,0), use.approx = c(F,F,T,F,T,F,F), seed = 1000, use.glmnet = c(F,F,F,F,F),
                                   use.interactions = F, use.logs = F, use.squares = F, use.cubics = F)


####################################################################################
# SIMULATED DATA
####################################################################################

## Setup some synthetic data
p = 350
I = 1:20
coef = c(3, -3)
beta = rep(0, p)
beta[I[1:10]] = coef[1]
beta[I[11:20]] = coef[2]

rho = 0.7
C = toeplitz(rho^(0:(p-1)))
beta = as.vector(beta)
SNR = sqrt( t(beta)%*%C%*%beta / 1 ) # Calibrate SNR so sigma2 = 1

#sigma2 = 1
scores_071 = continuous.synthetic.data.test(beta, cov.X=C, SNR, 70, 100, c("glmnet","hs.lambda2","hs.lambda2","MCP","SCAD", "hs.like", "glmnet.ridge"), 
                                                log.tau2.max=c(0,0,0,0,0,0,0), e = 1e-5, use.approx = c(F,F,T,F,F,F,F), seed = 1000, 
                                                use.glmnet = c(F,F,F,F,F,F,F))
#sigma2 = 9
scores_079 = continuous.synthetic.data.test(beta, cov.X=C, SNR/3, 70, 100, c("glmnet","hs.lambda2","hs.lambda2","MCP","SCAD", "hs.like", "glmnet.ridge"), 
                                                log.tau2.max=c(0,0,0,0,0,0,0), e = 1e-5, use.approx = c(F,F,T,F,T,F,F), seed = 1000, 
                                                use.glmnet = c(F,F,F,F,F,F,F))



rho = 0
C = toeplitz(rho^(0:(p-1)))
SNR = sqrt( t(beta)%*%C%*%beta / 1 )

#sigma2 = 1
scores_001 = continuous.synthetic.data.test(beta, cov.X=C, SNR, 70, 100, c("glmnet","hs.lambda2","hs.lambda2","MCP","SCAD", "hs.like", "glmnet.ridge"), 
                                                log.tau2.max=c(0,0,0,0,0,0,0), e = 1e-5, use.approx = c(F,F,T,F,F,F,F), seed = 1000, 
                                                use.glmnet = c(F,F,F,F,F,F,F))
#sigma2 = 9
scores_009 = continuous.synthetic.data.test(beta, cov.X=C, SNR/3, 70, 100, c("glmnet","hs.lambda2","hs.lambda2","MCP","SCAD", "hs.like", "glmnet.ridge"), 
                                                log.tau2.max=c(0,0,0,0,0,0,0), e = 1e-5, use.approx = c(F,F,T,F,F,F,F), seed = 1000, 
                                                use.glmnet = c(F,F,F,F,F,F,F))


