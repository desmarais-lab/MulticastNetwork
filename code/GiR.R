source("/Users/bomin8319/Desktop/MulticastNetwork/code/Multicast.R")

D = 100
A = 5
P = 5
#Q = 3
X = array(rnorm(D*A*A*P), dim = c(D,A,A,P))
X[,,,1] = 1
#Y = array(rnorm(D*A*Q), dim = c(D,A,Q))
Q = 3
Y = array(1, dim = c(D,A,Q))
for (a in 1:A) {
  Y[,a,2] = a
  if (a != 1) {
    Y[,a,3] = 0
  }
}
support = gibbs.measure.support(A-1)
prior.beta = list(mean = c(-3, rep(0, P-1)), var = diag(P))
prior.eta = list(mean = rep(0, Q), var = diag(Q))
prior.sigma2 = list(a = 3, b = 1)
Nsamp = 5000
outer = 50
inner = c(5, 5, 1)
burn = 0
#Schein test
result = matrix(NA, Nsamp, 2*(2+P+Q))
for (n in 1:Nsamp) {
  if (n %% 100 == 0) print(n)
  beta = rmvnorm_arma(1, prior.beta$mean, prior.beta$var)
  eta = rmvnorm_arma(1, prior.eta$mean, prior.eta$var)
  sigma2 = rinvgamma(1, prior.sigma2$a, prior.sigma2$b)
  initial = Generate(D, A, beta, eta, sigma2, X, Y, support)
  infer = Inference(initial$data, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, 
                    initial = initial, proposal.var = c(0.01, 0.1, 0.01, 0.5), lasttime = 0)
  result[n, ] = c(mean(vapply(initial$u, function(x) rowSums(x), rep(0, A))),
                  initial$beta, initial$eta, initial$sigma2, 
                  mean(vapply(infer$u, function(x) rowSums(x), rep(0, A))),
                  infer$beta[outer,], infer$eta[outer,], infer$sigma2[outer,])			 
}
par(mfrow=c(2,5))
GiR_PP_Plots(result[,c(1:(2+P+Q))], result[,c((3+P+Q):(2*(2+P+Q)))])
