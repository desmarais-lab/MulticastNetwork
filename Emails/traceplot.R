load('~/Desktop/MulticastNetwork/Emails/Montgomery_infer.RData')
load('~/Desktop/MulticastNetwork/Emails/Montgomery_infer2.RData')

library(ggplot2)
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(latex2exp)
library(MCMCpack)
beta = t(sapply(5000+c(1:5000), function(x) Montgomery_infer$beta[8*x, c(1,5,10)]))
beta = t(sapply(1:5000, function(x) Montgomery_infer$beta[8*x,c(3,5,8)]))
colnames(beta) = c("beta3","beta5", "beta8")
S <-  ggs(mcmc(beta))
plot = ggs_traceplot(S)+theme(plot.title = element_blank(),text = element_text(size = rel(4.25)),legend.position="none")+ scale_x_continuous(breaks=seq(0,5000,2500))
density = ggs_density(S)+theme(plot.title = element_blank(),text = element_text(size = rel(4.5)),legend.position="none")
grid.arrange(plot, density, ncol = 2, nrow = 1)


beta = t(sapply(1:5000, function(x) Montgomery_infer$beta[8*x, ]))
beta = t(sapply(1:800, function(x) Montgomery_infer$beta[50*x,]))
colnames(beta) = c(sapply(1:14, function(x) paste0("beta",x)))
S <-  ggs(mcmc(beta))
S$Parameter <- factor(S$Parameter, levels = c(sapply(1:14, function(x) paste0("beta",x))))
ggs_geweke(S)+scale_colour_manual(values = c("black","black"))+theme(plot.title = element_blank(),text = element_text(size = rel(4.25)),legend.position="none")



eta = t(sapply(1:5000, function(x) Montgomery_infer$eta[8*x, c(2,7)]))
eta = cbind(eta, sapply(1:5000, function(x) Montgomery_infer$sigma2[8*x]))
colnames(eta) = c("eta2", "eta7", "sigma2")
S <-  ggs(mcmc(eta))
plot = ggs_traceplot(S)+theme(plot.title = element_blank(),text = element_text(size = rel(4.25)),legend.position="none")+ scale_x_continuous(breaks=seq(0,5000,2500))
density = ggs_density(S)+theme(plot.title = element_blank(),text = element_text(size = rel(4.25)),legend.position="none")

grid.arrange(plot, density, ncol = 2, nrow = 1)

eta = t(sapply(2000:5000, function(x) Montgomery_infer$eta[8*x, ]))
colnames(eta) = c(sapply(1:7, function(x) paste0("eta",x)))
S <-  ggs(mcmc(eta))
ggs_geweke(S)+scale_colour_manual(values = c("black","black"))+theme(plot.title = element_blank(),text = element_text(size = rel(4.25)),legend.position="none")




