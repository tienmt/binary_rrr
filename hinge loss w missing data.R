logit = function(x)1*(exp(x)/(1+exp(x))>=0.5)
P_Omega = function(a,entri){ a[entri] = 0; return(a) }; P_Omega <- compiler::cmpfun(P_Omega)
library(rrpack)
# random data generation
n = 100  # samples
l = 8   # response
p = 12   # predictors
r = 2    # true rank
missfrac = 0.3
simdata <- rrr.sim3(n = n, p = p, q.mix = c(0, l, 0), intercept = rep(0,l), nrank = r, mis.prop = 0)
family <- simdata$family 
control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 1000, trace =FALSE,gammaC0 = 1.1, plot.cv = F,conv.obj=TRUE)

X = matrix(rnorm(n*p),nc=p); tX = t(X)
B = matrix(rnorm(p*r),nr=p)%*%t(matrix(rnorm(l*r),nr=l)) #*2 + rnorm(p*l,0,1)
#Z = sample(c(-1,1),size = n*l,replace = T,prob = c(0.1,0.9))
# hinge model
Ytrue = sign( X%*%B ) # *Z #  + rnorm(n*l)
#MU0 = X%*%B  + rnorm(n*l)
# binomial model
#prob = exp(MU0)/(1+ exp(MU0) )
#Ytrue = apply(prob, 2, function(a) rbinom(n = n, size = 1, a))

C = B
Yna <- Ytrue; Yna[Ytrue==-1] <-0 ; library(rrpack)
imiss = sample(n*l,size = n*l*missfrac,replace = F)
Yna[ imiss ] <- NA
fit.cv.mrrr <- cv.mrrr(Yna , X, family = family,control = control, 
                       penstr = list(penaltySVD = "rankCon"))
hatc = coef(fit.cv.mrrr$fit)[-1,]
mean((hatc - C)^2 ) ; mean((fit.cv.mrrr$fit$mu - X%*%B)^2 )

Y = Yna; Y[Y == 0] = -1
lam = 1  # in the prior
ystar = diag(p)*lam^2
Bm_hinge = matrix( 0 ,nr=p,nc=l)
Iters = 15000
burnin = 1000
h = 1/(p*l)^1.65
a = 0  
M =  hatc
for(s in 1:Iters){
  MtM = tcrossprod(M)
  tam1 = solve(ystar + MtM, M )
  hxtyxm = h*tX %*% P_Omega(Y*((1-Y*(X%*%M))>0),imiss)
  tam = M + hxtyxm + h*(p+l+2)*tam1 + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
  
  tattam = tcrossprod(tam)
  pro.tam = -sum(P_Omega((1 - Y*(X%*%tam))*(1 - Y*(X%*%tam)>0),imiss) ) - 0.5*(p+l+2)*determinant(ystar + tattam)[[1]][1]
  pro.M = - sum(P_Omega((1 - Y*(X%*%M))*(1 - Y*(X%*%M)>0),imiss) ) - 0.5*(p+l+2)*determinant(ystar + MtM)[[1]][1]
  
  tam2 = solve(ystar + tattam,tam)
  tran.m = -sum((M-tam -h*tX%*%P_Omega(Y*((1-Y*(X%*%tam))>0),imiss) - h*(p+l+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M - hxtyxm - h*(p+l+2)*tam1 )^2)/(4*h)
  
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
  M1 = M1 + h_lmc*tX %*% P_Omega(Y*((1-Y*(X%*%M1))>0),imiss) + h_lmc*(p+l+2)*tam + sqrt(2*h_lmc)*rnorm(p*l)
  if(s>burnin) Mlmc_hinge =  Mlmc_hinge + M1/(Iters-burnin) }
c(Matrix::rankMatrix(hatc,.1)[1],  mean((hatc - C )^2)/mean(C^2) )
c(Matrix::rankMatrix(Mlmc_hinge,.1)[1],mean((Mlmc_hinge - C )^2)/mean(C^2) )
c(Matrix::rankMatrix(Bm_hinge,.1)[1],  mean((Bm_hinge - C )^2)/mean(C^2) )
a/Iters
source('logit_model_missing.R')
a/Iters;  Ytrue[Ytrue == 0] = -1
c(mean(sign(X%*%C) != Ytrue) ,mean(sign(X%*%hatc) != Ytrue ) )
c(mean(sign(X%*%Bm_hinge) != Ytrue) , mean(sign(X%*%Bm_logit) != Ytrue) )
c(mean(sign(X%*%Mlmc_hinge) != Ytrue) , mean(sign(X%*%Blmc_logit) != Ytrue))
