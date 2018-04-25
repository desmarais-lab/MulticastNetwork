#install the package from Github
library(devtools)
install_github("desmarais-lab/MulticastNetwork/pkg")

#generate synthetic data




#load our Montgomery county email data
load("/Users/bomin8319/Box/gainlab_example/Bomin/Montgomery.RData")

#data including sender a_d, reciever vector r_d, timestamp t_d
edge = Montgomery$edge
head(edge)

#covariates affecting "who sends to whom"
X = Montgomery$X
head(X[100, 1, , ])
 
#covariates affecting "when to send"
Y = Montgomery$Y
head(Y[100, , ])

#run inference to estimate beta, eta, u, and sigma2
prior.beta = list(mean = rep(0, P), var = 2*diag(P))
prior.eta = list(mean = rep(0, Q), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)

outer = 100
inner = c(1, 1, 1)
burn = 0

Montgomery_infer = Inference(edge, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
		  proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = , timedist = "lognormal")
