library(MulticastNetwork)
load("/Users/bomin8319/Box/gainlab_example/Bomin/Montgomery.RData")
edge = Montgomery$edge
X = Montgomery$X
Y = Montgomery$Y
P = dim(X)[4]
Q = dim(Y)[3]
A = dim(Y)[2]

#run inference to estimate beta, eta, u, and sigma2
prior.beta = list(mean = c(-3.5, rep(0, P-1)), var = 2*diag(P))
prior.eta = list(mean = c(7, rep(0, Q-1)), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)

outer = 500
inner = c(1, 1, 1)
burn = 0

#run infernece
Montgomery_infer = Inference(edge, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
		  proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, lasttime = Montgomery$lasttime, timedist = "lognormal")

# generate data from the model estimates
#Montgomery_PPC = PPC(length(edge), beta = colMeans(Montgomery_infer$beta), eta = colMeans(Montgomery_infer$eta), 
#                     sigma2 = mean(Montgomery_infer$sigma2), X, Y, timeunit = 3600, u = Montgomery_infer$u, timedist = "lognormal")


set.seed(1)
missing = list()
#missingness of senders
missing[[1]] = matrix(0, nrow = dim(Y)[1], 1)    
missing[[1]][sample(1:dim(Y)[1], 62, replace = FALSE), ] = 1
#missingness of receivers
missing[[2]] = matrix(0, nrow = dim(Y)[1], A)    
missing[[2]][sample(1:(dim(Y)[1]*A), 1118, replace = FALSE)] = 1
#missingness of timestamps
missing[[3]] = matrix(0, nrow = dim(Y)[1], 1)
missing[[3]][sample(1:dim(Y)[1], 62, replace = FALSE), ] = 1


for (d in 1:dim(Y)[1]) {
	if (missing[[1]][d,1] == 1) {
		edge[[d]]$a_d = NA
	}
	if (sum(missing[[2]][d,]) > 0) {
		edge[[d]]$r_d[which(missing[[2]][d,]==1)] = NA
	}
	if (missing[[3]][d,1] == 1) {
		edge[[d]]$t_d = NA
	}
}


initial = list()
initial$beta = colMeans(Montgomery_infer$beta)
initial$eta =  colMeans(Montgomery_infer$eta)
initial$u = Montgomery_infer$u
initial$sigma2 = mean(Montgomery_infer$sigma2)

#will generate 10 predictions (iterate two steps: imputation -> inference)
Montgomery_PPE = PPE(edge, X, Y, 550, c(5,5,1), 50, prior.beta, prior.eta, prior.sigma2, 
                     initial = initial, proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, 
                     lasttime = Montgomery$lasttime, MHprop.var = 0.15, timedist = "lognormal")

save(Montgomery_PPE, file = "/Users/bomin8319/Desktop/MulticastNetwork/Emails/Montgomery_PPE.RData")

initial = list()
initial$beta = colMeans(Montgomery_infer2$beta)
initial$eta =  colMeans(Montgomery_infer2$eta)
initial$u = Montgomery_infer2$u
initial$sigma2 = mean(Montgomery_infer2$sigma2)

Montgomery_PPE2 = PPE(edge, X, Y, 550, c(5,5,1), 50, prior.beta, prior.eta, prior.sigma2, 
                     initial = initial, proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, 
                     lasttime = Montgomery$lasttime, MHprop.var = 0.15, timedist = "exponential")

save(Montgomery_PPE2, file = "/Users/bomin8319/Desktop/MulticastNetwork/Emails/Montgomery_PPE2.RData")
save(Montgomery_PPE2, file = "/Users/bomin8319/Desktop/Montgomery_PPE2.RData")


names(Montgomery_PPE)

library(gridExtra)
library(latex2exp)
p = list()
#################
truesender = sapply(Montgomery_PPE$sendermissing, function(d) edge[[d]]$a_d)
predprob = Montgomery_PPE$senderprob/rowSums(Montgomery_PPE$senderprob)
predprob2 = Montgomery_PPE2$senderprob/rowSums(Montgomery_PPE2$senderprob)

sender = data.frame(probtrue = sapply(1:62, function(d) predprob[d,truesender[d]]), dist = rep("lognormal", 62))
sender = rbind(sender, data.frame(probtrue = sapply(1:62, function(d) predprob2[d,truesender[d]]), dist = rep("exponential", 62)))
library(ggplot2)
library(reshape)
sender = melt(sender)
colnames(sender)[3] = "correct"
p[[1]]=ggplot(data = sender, aes(x = dist, y = correct, fill = dist))+geom_boxplot()+ylab("Correct sender posterior probability")+geom_hline(yintercept=1/18, col = 'blue', lty=2, lwd=1)
boxplot(sapply(1:62, function(d) predprob[d,truesender[d]]), sapply(1:62, function(d) predprob2[d,truesender[d]]))

###############################################
truereceiver = unlist(sapply(Montgomery_PPE$receivermissing, function(d) edge[[d]]$r_d[which(missing[[2]][d,]==1)]))
predprob = Montgomery_PPE$receiverprob
predprob2 = Montgomery_PPE2$receiverprob

probtrue = sapply(1:1118, function(d) predprob[d,truereceiver[d]+1])
#probtrue[probtrue==1] = 0.999
probtrue2 = sapply(1:1118, function(d) predprob2[d,truereceiver[d]+1])
#probtrue2[probtrue2==1] = 0.999

receiver = data.frame(probtrue = log(probtrue /(1-probtrue)), dist = rep("lognormal", 1118))
receiver = rbind(receiver, data.frame(probtrue = log(probtrue2/(1-probtrue2)), dist = rep("exponential", 1118)))
receiver = melt(receiver)
colnames(receiver)[3] = "logit"
ggplot(data = receiver, aes(x = dist, y = logit, fill = dist))+geom_boxplot()
p[[2]]=ggplot(data = receiver, aes(x = logit, fill = dist))+geom_histogram(position = "dodge")
#boxplot(sapply(1:62, function(d) predprob[d,truesender[d]]), sapply(1:62, function(d) predprob2[d,truesender[d]]))

######################
library(MLmetrics)

truereceiver = unlist(sapply(Montgomery_PPE$receivermissing, function(d) edge[[d]]$r_d[which(missing[[2]][d,]==1)]))
predreceiver = Montgomery_PPE$receiverpredict
predreceiver2 = Montgomery_PPE2$receiverpredict

probtrue = sapply(1:500, function(d) F1_Score(truereceiver, predreceiver[,d], 1))
probtrue2 = sapply(1:500, function(d) F1_Score(truereceiver, predreceiver2[,d], 1))

receiver = data.frame(F1score = probtrue, dist = rep("lognormal", 500))
receiver = rbind(receiver, data.frame(F1score = probtrue2, dist = rep("exponential", 500)))
receiver = melt(receiver)
colnames(receiver)[3] = "F1score"
p[[2]]=ggplot(data = receiver, aes(x = dist, y = F1score, fill = dist))+geom_boxplot()


#################################################
truetime = sapply(Montgomery_PPE$timemissing, function(d) edge[[d]]$t_d-edge[[d-1]]$t_d)/3600
predtime = Montgomery_PPE$timepredict
predtime2 = Montgomery_PPE2$timepredict

time = data.frame(MdAPE = sapply(1:62, function(d) median(abs((predtime[d,]-truetime[d])/truetime[d]))), dist = rep("lognormal", 62))
time = rbind(time, data.frame(MdAPE = sapply(1:62, function(d) median(abs((predtime2[d,]-truetime[d])/truetime[d]))), dist = rep("exponential", 62)))
time$MdAPE = log(time$MdAPE)
time = melt(time)
stats = boxplot.stats(time$value)$stats

colnames(time)[3] = "MdAPE"
p[[3]]=ggplot(data = time, aes(x = dist, y = MdAPE, fill = dist))+geom_boxplot()+theme(legend.position = "bottom")+ylab("log(MdAPE)")
ggplot(data = time, aes(x = MdAPE, fill = dist))+geom_histogram(position = "dodge")

####################################
time = data.frame(MdAE = sapply(1:62, function(d) mean(abs((predtime[d,]-truetime[d])/truetime[d]))), dist = rep("lognormal", 62))
time = rbind(time, data.frame(MdAE = sapply(1:62, function(d) mean(abs((predtime2[d,]-truetime[d])/truetime[d]))), dist = rep("exponential", 62)))
time$MdAE = log(time$MdAE)

time = melt(time)
colnames(time)[3] = "MdAE"
ggplot(data = time, aes(x = dist, y = MdAE, fill = dist))+geom_boxplot()
#####################################
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(p[[3]])
p3 <- grid.arrange(arrangeGrob(
				p[[1]] + theme(legend.position="none"),
                        p[[2]] + theme(legend.position="none"),
                         p[[3]] + theme(legend.position="none"), 
                       nrow=1), mylegend,
        heights=c(10, 1))




