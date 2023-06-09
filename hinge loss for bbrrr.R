logit = function(x)1*(exp(x)/(1+exp(x))>=0.5)
library(Rcpp)
sourceCpp('test.cpp')
library(rrpack)
# random data generation
n = 100  # samples
l = 8   # response
p = 12   # predictors
r = 2    # true rank
simdata <- rrr.sim3(n = n, p = p, q.mix = c(0, l, 0), intercept = rep(0,l), nrank = r, mis.prop = 0)
family <- simdata$family 
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000, trace =FALSE,gammaC0 = 1.1, plot.cv = F,conv.obj=TRUE)

X = matrix(rnorm(n*p),nc=p); tX = t(X)
B = matrix(rnorm(p*r),nr=p)%*%t(matrix(rnorm(l*r),nr=l))# *2 + rnorm(p*l,0,1)
Z = sample(c(-1,1),size = n*l,replace = T,prob = c(0.1 , 0.9))
# hinge model
Ytrue = sign( X%*%B ) 
Y = sign( X%*%B + rnorm(n*l)) 

#switching model
Ytrue = sign( X%*%B + rnorm(n*l)) 
Y = Ytrue*Z 

# binomial model
MU0 = X%*%B  + rnorm(n*l)
prob = exp(MU0)/(1+ exp(MU0) )
Y = apply(prob, 2, function(a) rbinom(n = n, size = 1, a))
Ytrue =Y

### mRRR
C = B
Yna <- Y; Yna[Y==-1] <-0 ; library(rrpack)
fit.cv.mrrr <- cv.mrrr(Yna , X, family = family,control = control, penstr = list(penaltySVD = "rankCon"))
hatc = coef(fit.cv.mrrr$fit)[-1,]
mean((hatc - C)^2 ) 

Y[Y == 0] = -1
tau = 1  # in the prior
ystar = diag(p)*tau^2
Bm_hinge = matrix( 0 ,nr=p,nc=l)
Iters = 15000
burnin = 1000
h = 1/(p*l)^1.9
a = 0  
M =  hatc
lambda = 10
for(s in 1:Iters){
  MtM = tcrossprod(M)
  tam1 = solve(ystar + MtM, M )
  hxtyxm = h*tX%*%(Y*((1-Y*(X%*%M))>0))
  tam = M + hxtyxm*lambda + h*(p+l+2)*tam1 + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
  
  tattam = tcrossprod(tam)
  pro.tam = - sum((1 - Y*(X%*%tam))*(1 - Y*(X%*%tam)>0) )*lambda  - 0.5*(p+l+2)*determinant(ystar + tattam)[[1]][1]
  pro.M = - sum((1 - Y*(X%*%M))*(1 - Y*(X%*%M)>0) )*lambda  - 0.5*(p+l+2)*determinant(ystar + MtM)[[1]][1]
  
  tam2 = solve(ystar + tattam,tam)
  tran.m = -sum((M-tam -h*tX%*%(Y*((1-Y*(X%*%tam))>0))*lambda  - h*(p+l+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M - hxtyxm*lambda  - h*(p+l+2)*tam1 )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) <= pro.trans){
    M = tam;  print(s);   a = a+1
  } 
  if (s>burnin)Bm_hinge = Bm_hinge + M/(Iters-burnin)
}

#Langevin MC for BRRR
M1 = hatc
h_lmc =  h/n/n
Mlmc_hinge = matrix( 0 ,nr=p,nc=l)
for(s in 1:Iters){
  tam = solve(a = ystar + tcrossprod(M1),M1)
  M1 = M1 + h_lmc*tX%*%(Y*((1-Y*(X%*%M1))>0))*lambda + h_lmc*(p+l+2)*tam + sqrt(2*h_lmc)*rnorm(p*l)
  if(s>burnin) Mlmc_hinge =  Mlmc_hinge + M1/(Iters-burnin) }
c(Matrix::rankMatrix(hatc,.1)[1],  mean((hatc - C )^2)/mean(C^2) )
c(Matrix::rankMatrix(Mlmc_hinge,.1)[1],mean((Mlmc_hinge - C )^2)/mean(C^2) )
c(Matrix::rankMatrix(Bm_hinge,.1)[1],  mean((Bm_hinge - C )^2)/mean(C^2) )
a/Iters
source('logit model.R')
a/Iters
c(mean(sign(X%*%C) != Y) ,mean(sign(X%*%hatc) != Y) )
c(mean(sign(X%*%Bm_hinge) != Y) , mean(sign(X%*%Bm_logit) != Y) )
c(mean(sign(X%*%Mlmc_hinge) != Y) , mean(sign(X%*%Blmc_logit) != Y))
mean((hatc - C)^2 ) ; mean((Bm_hinge - C)^2 ) ;mean((Bm_logit - C)^2 ) ;mean((Mlmc_hinge - C)^2 )  ;mean((Blmc_logit - C)^2 ) 
