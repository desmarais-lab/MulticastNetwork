source("/Users/bomin8319/Desktop/MulticastNetwork/code/Multicast.R")

D = 100
A = 5
P = 4
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
Nsamp = 10^5
outer = 10
inner = c(5, 5, 5)
burn = 0
#Schein test
result = matrix(NA, Nsamp, 2*(5+P+Q))
for (n in 1:Nsamp) {
  if (n %% 100 == 0) print(n)
  beta = rmvnorm_arma(1, prior.beta$mean, prior.beta$var)
  eta = rmvnorm_arma(1, prior.eta$mean, prior.eta$var)
  sigma2 = rinvgamma(1, prior.sigma2$a, prior.sigma2$b)
  initial = Generate(D, A, beta, eta, sigma2, X, Y, support)
  infer = Inference(initial$data, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, 
                    initial = initial, proposal.var = c(0.01, 0.1, 0.01, 0.5), lasttime = 0)
  initial2 = Generate(D, A, infer$beta[outer,], infer$eta[outer,], infer$sigma2[outer,], X, Y, support)                  
  result[n, ] = c(mean(vapply(initial$data, function(x) sum(x$r_d), c(1))),
                  var(vapply(initial$data, function(x) sum(x$r_d), c(1))),
 				 mean(vapply(2:D, function(d) initial$data[[d]]$t_d-initial$data[[d-1]]$t_d, c(1))),
 				var(vapply(2:D, function(d) initial$data[[d]]$t_d-initial$data[[d-1]]$t_d, c(1))),
                  initial$beta, initial$eta, initial$sigma2, 
                 mean(vapply(initial2$data, function(x) sum(x$r_d), c(1))),
 			var(vapply(initial2$data, function(x) sum(x$r_d), c(1))),
                  mean(vapply(2:D, function(d) initial2$data[[d]]$t_d-initial2$data[[d-1]]$t_d, c(1))),
 				var(vapply(2:D, function(d) initial2$data[[d]]$t_d-initial2$data[[d-1]]$t_d, c(1))),
                  initial2$beta, initial2$eta, initial2$sigma2)			 
}
par(mfrow=c(3,4))
GiR_PP_Plots(result[,c(1:(5+P+Q))], result[,c((6+P+Q):(2*(5+P+Q)))])
save(result, file = "/Users/bomin8319/Desktop/result.RData")

setwd("/Users/bomin8319/Desktop/")
Forward_stats = result[,c(1:(5+P+Q))]
Backward_stats = result[,c((6+P+Q):(2*(5+P+Q)))]
colnames(Forward_stats) = c("Mean Recipient Size", "Mean Time-increments", sapply(1:P, function(p){paste0("b",p," Estimates")}), sapply(1:Q, function(q){paste0(expression(eta),q," Estimates")}), "Variance Estimate")
plot = list()
library(ggplot2)
 nms = colnames(Forward_stats)
  
  for (i in 1:ncol(Forward_stats)) {
    all = c(Backward_stats[, i], Forward_stats[, i])
    
    quantiles = 1000
    
    uniqueValues = quantile(all,seq(0, 1, length = quantiles))
    qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
  		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(Forward_stats[, i] <= uniqueValues[j])
  		qx2[j] = mean(Backward_stats[, i] <= uniqueValues[j])
  	}
  	data = data.frame(qx1 = qx1, qx2 = qx2)
	plot[[i]] =qplot(data = data, x = qx1, y = qx2, colour = I("blue"), size = I(2)) + labs(x = "Forward", y = "Backward")+geom_abline(intercept = 0, col = 'red', size = I(1)) + theme(text = element_text(size = rel(5)), plot.title= element_text(size=rel(3.5)))
	#+labs(title = expression(nms[i]))
mname = paste0("plot", i, ".png")
print(plot[[i]])
ggsave(filename = mname)
}	
library(gridExtra)
library(latex2exp)
grid.arrange(arrangeGrob(plot[[1]]+labs(title = "Mean of Recipient Sizes"),	plot[[2]]+labs(title = "Var of Recipient Sizes"),
	plot[[3]]+labs(title = "Mean of Time-increments"),	plot[[4]]+labs(title = "Var of Time-increments"),
	plot[[5]]+labs(title = expression(b*" Estimates for p=1")),plot[[6]]+labs(title = expression(b*" Estimates for p=2")),plot[[7]]+labs(title = expression(b*" Estimates for p=3")),
plot[[8]]+labs(title = expression(b*" Estimates for p=4")),plot[[9]]+labs(title = expression(eta*" Estimates for q=1")),plot[[10]]+labs(title = expression(eta*" Estimates for q=2")),
plot[[11]]+labs(title = expression(eta*" Estimates for q=3")),plot[[12]]+labs(title = expression(sigma^2*" Estimates")), nrow=3))
