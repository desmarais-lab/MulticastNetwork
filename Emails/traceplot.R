load('~/Desktop/MulticastNetwork/Emails/Montgomery_infer.RData')
load('~/Desktop/MulticastNetwork/Emails/Montgomery_infer2.RData')

library(ggplot2)
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(latex2exp)
library(MCMCpack)
beta = t(sapply(1:5000, function(x) Montgomery_infer2$beta[8*x, c(2,5,10)]))
colnames(beta) = c("beta2","beta5", "beta10")
S <-  ggs(mcmc(beta))
plot = ggs_traceplot(S)
density = ggs_density(S)
grid.arrange(plot, density, ncol = 2, nrow = 1)


beta = t(sapply(1:5000, function(x) Montgomery_infer2$beta[8*x, ]))
colnames(beta) = c(sapply(1:11, function(x) paste0("beta",x)))
S <-  ggs(mcmc(beta))
ggs_geweke(S)



eta = t(sapply(1:5000, function(x) Montgomery_infer$eta[8*x, c(2,7)]))
eta = cbind(eta, sapply(1:5000, function(x) Montgomery_infer$sigma2[8*x]))
colnames(eta) = c("eta2", "eta7", "sigma2")
S2 <-  ggs(mcmc(eta))
plot = ggs_traceplot(S)
density = ggs_density(S)
grid.arrange(plot, density, ncol = 2, nrow = 1)

eta = t(sapply(1:5000, function(x) Montgomery_infer2$eta[8*x, ]))
colnames(eta) = c(sapply(1:7, function(x) paste0("eta",x)))
S <-  ggs(mcmc(eta))
ggs_geweke(S)



