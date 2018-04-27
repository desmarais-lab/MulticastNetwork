#install the package from Github
library(devtools)
install_github("desmarais-lab/MulticastNetwork/pkg")
library(MulticastNetwork)
#load our Montgomery county email data
load("/Users/bomin8319/Box/gainlab_example/Bomin/Montgomery.RData")
names(Montgomery)

#data including sender a_d, reciever vector r_d, timestamp t_d
edge = Montgomery$edge
head(edge)

#covariates affecting "who sends to whom"
X = Montgomery$X
head(X[100, 1, , ])
X = X[,,,c(1:5)]    #use few covariates for this application

#covariates affecting "when to send"
Y = Montgomery$Y
head(Y[100, , ])
Y = Y[,,c(1:5)]    #use few covariates for this application


P = dim(X)[4]
Q = dim(Y)[3]
#run inference to estimate beta, eta, u, and sigma2
prior.beta = list(mean = rep(0, P), var = 2*diag(P))
prior.eta = list(mean = rep(0, Q), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)

outer = 500
inner = c(1, 1, 1)
burn = 0

#this initial values are my results after convergence
initialval = list()
initialval$beta = c(-3.548488159, -0.108898033,  0.086928281,  0.274102964,  0.029273815)
initialval$eta = c(7.38197540,  0.12883606, -1.07065227, -0.20698218, -0.05894448)
initialval$sigma2 = 14.09288
initialval$u = lapply(1:dim(X)[1], function(d) matrix(0, dim(X)[2], dim(X)[2]))

#run infernece
Montgomery_infer = Inference(edge, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initialval = initialval,
		  proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, lasttime = Montgomery$lasttime, timedist = "lognormal")
names(Montgomery_infer)

# check convergence
plot(Montgomery_infer$loglike, type = 'l')

# generate data from the model estimates
Montgomery_PPC = PPC(length(edge), beta = colMeans(Montgomery_infer$beta), eta = colMeans(Montgomery_infer$eta), 
                     sigma2 = mean(Montgomery_infer$sigma2), X, Y, timeunit = 3600, u = Montgomery_infer$u, 
        			 lasttime = Montgomery$lasttime, timedist = "lognormal")
