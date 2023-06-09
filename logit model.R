library(Rcpp)

Y[Y == 0] = -1
Bm_logit = matrix(data=0,nr=p,nc=l)
h = 1/(p*l)^1.7
a = 0  
M =  hatc
for(s in 1:Iters){
  MtM = tcrossprod(M)
  tam1 = solve(ystar + MtM, M )
  XM = eigenMapMatMult(X , M)
  hxtyxm = h * eigenMapMatMult(tX, Y* exp(-Y*XM)/(1+exp(-Y*XM)) )
  tam = M + hxtyxm + h*(p+l+2)*tam1 + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
  
  tattam = tcrossprod(tam)
  Xtam = eigenMapMatMult(X, tam)
  pro.tam = - sum(log(1+exp(-Y*Xtam)) ) - 0.5*(p+l+2)*determinant(ystar + tattam)[[1]][1]
  pro.M = - sum(log(1+exp(-Y*XM))) - 0.5*(p+l+2)*determinant(ystar + MtM)[[1]][1]
  
  tam2 = solve(ystar + tattam,tam)
  tran.m = -sum((M-tam -h* eigenMapMatMult(tX, Y* exp(-Y*Xtam)/(1+exp(-Y*Xtam)) ) - h*(p+l+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M - hxtyxm - h*(p+l+2)*tam1 )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) <= pro.trans){
    M = tam; a = a+1  } 
  if (s>burnin)Bm_logit = Bm_logit+ M/(Iters-burnin)
} 
M1 = hatc ; h_lmc =  h/n/n
Blmc_logit = matrix(data=0,nr=p,nc=l)
for(s in 1:Iters){
  tam = solve(a = ystar + tcrossprod(M1),M1)
  XM = eigenMapMatMult(X , M1)
  M1 = M1 + h_lmc * eigenMapMatMult(tX, Y* exp(-Y*XM)/(1+exp(-Y*XM)) )+ h_lmc *(p+l+2)*tam + sqrt(2*h_lmc)*rnorm(p*l)
  if(s>burnin) Blmc_logit =  Blmc_logit + M1/(Iters-burnin)
}


