library(MulticastNetwork)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
set.seed(526113322)
nDocs = 100
node = 1:4
netstat = c("dyadic")
#timestat = c("timeofday", "dayofweek")
timestat = c()
L = 3
P = 6
prior.b = list(rep(0.5, P), 0.5* diag(P))
prior.delta = c(-3, 0.1)
#prior.eta = list(rep(2.5, length(node) + length(timestat)), 0.5*diag(length(node) +length(timestat)))
prior.eta = list(rep(2.5, 1), 0.5*diag(1))
prior.tau =  c(2,1)
sigma.Q = c(0.075, 0.005, 0.25, 0.05)

b = prior.b[[1]]
eta = prior.eta[[1]]
delta = prior.delta[1]
sigma_tau = prior.tau[2]/(prior.tau[1]-1)
support = gibbs.measure.support(length(node)-1)
base.data = GenerateDocs(500, node, b, eta, delta, sigma_tau, support, netstat, timestat,
                        base.data= NULL, backward = FALSE, base = TRUE)


Outer = 5
Inner = c(10, 10, 10)
Outer = 1
Inner = c(5, 5, 5)

Schein <- Schein(10000, nDocs, node, prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner,
               netstat = c("dyadic"), timestat = timestat,
              base.data = base.data, generate_PP_plots = TRUE)

