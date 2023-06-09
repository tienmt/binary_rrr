
Y = Yna; Y[Y == 0] = -1
Bm_logit = matrix(data=0,nr=p,nc=l)
h = 1/(p*l)^1.65
a = 0  
M =  hatc
for(s in 1:Iters){
  MtM = tcrossprod(M)
  tam1 = solve(ystar + MtM, M )
  Y_o_XM = Y * (X%*%M)
  hxtyxm = h * tX%*%P_Omega(Y* exp(-Y_o_XM)/(1+exp(-Y_o_XM)),imiss)
  tam = M + hxtyxm + h*(p+l+2)*tam1 + sqrt(2*h)*matrix(rnorm(p*l),nr = p)
  
  tattam = tcrossprod(tam)
  Xtam = X%*%tam
  pro.tam = - sum(P_Omega(log(1+exp(-Y*Xtam)),imiss ) ) - 0.5*(p+l+2)*determinant(ystar + tattam)[[1]][1]
  pro.M = - sum(P_Omega(log(1+exp(-Y_o_XM)),imiss ) ) - 0.5*(p+l+2)*determinant(ystar + MtM)[[1]][1]
  
  tam2 = solve(ystar + tattam,tam)
  tran.m = -sum((M-tam -h* tX%*%P_Omega(Y* exp(-Y*Xtam)/(1+exp(-Y*Xtam)),imiss) - h*(p+l+2)*tam2 )^2)/(4*h)
  tran.tam = -sum((tam-M - hxtyxm - h*(p+l+2)*tam1 )^2)/(4*h)
  
  pro.trans = pro.tam+tran.m-pro.M-tran.tam
  if(log(runif(1)) <= pro.trans){
    M = tam; a = a+1  } 
  if (s>burnin)Bm_logit = Bm_logit+ M/(Iters-burnin)
};  

#Langevin MC for BRRR
M1 = hatc ; h_lmc =  h/n/n
Blmc_logit = matrix(data=0,nr=p,nc=l)
Y_Omega = P_Omega(Y,imiss)
for(s in 1:Iters){
  tam = solve(a = ystar + tcrossprod(M1),M1)
  Y_o_XM = Y_Omega * (X%*%M1)
  M1 = M1 + h_lmc * tX%*%(Y_Omega* exp(-Y_o_XM)/(1+exp(-Y_o_XM)))+ h_lmc *(p+l+2)*tam + sqrt(2*h_lmc)*rnorm(p*l)
  if(s>burnin) Blmc_logit =  Blmc_logit + M1/(Iters-burnin) }

