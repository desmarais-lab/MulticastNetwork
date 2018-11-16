#source("/Users/bomin8319/Desktop/MulticastNetwork/code/Multicast.R")
library(MulticastNetwork)
library(lubridate)
load('~/Desktop/MulticastNetwork/code/Temporal_Email_Data.Rdata')
Montgomery = Temporal_Email_Data$Montgomery
email = Montgomery$email_data
email = unique(email)
email$timepoints =  as.numeric(as.POSIXct(strptime(email[,1], "%d %b %Y %H:%M:%S")))
email = email[order(email$timepoints), ]
edge = list()
initialtime =  as.numeric(as.POSIXct(strptime("01 Mar 2012 00:00:00", "%d %b %Y %H:%M:%S")))
for (d in 1:nrow(email)) {
	t_d = email[d, 21] - initialtime
	edge[[d]] = list(a_d = email[d,2], r_d = as.numeric(email[d,-c(1:2, 21)]), t_d = t_d)
}

uniqtime = unique(email$timepoints)
# construct time covariates Y
D = length(edge)
A = length(Montgomery$manager_gender)
Q = 7
Y = array(1, dim = c(D,A,Q))
for (a in 1:A) {
	Y[,a,2] = 1* (Montgomery$manager_gender[a]=="Female")
	Y[,a,3] = 1* (Montgomery$manager_department[a]=="County Manager")
}
timeunit = 3600
Y[1,,6] = rep(as.numeric(wday(as.POSIXct(strptime("01 Mar 2012 00:00:00", "%d %b %Y %H:%M:%S"))) %in% c(1, 7)), A)
Y[1,,7] = rep(pm(as.POSIXct(strptime("01 Mar 2012 00:00:00", "%d %b %Y %H:%M:%S"))), A)
for (d in 2:D) {
	index = which(email$timepoints >= uniqtime[which(uniqtime==email$timepoints[d])-1]-7*24*timeunit & email$timepoints < email$timepoints[d])
	sent = email[index, 2]
	received = email[index, 3:(2+A)]
	Y[d, ,4] = tabulate(sent, A) 
	Y[d, ,5] = colSums(received)
	Y[d,,6] = rep(as.numeric(wday(as.POSIXct(strptime(email[d-1,1], "%d %b %Y %H:%M:%S"))) %in% c(1, 7)), A)
	Y[d,,7] = rep(as.numeric(pm(as.POSIXct(strptime(email[d-1,1], "%d %b %Y %H:%M:%S")))), A)
}


sendraw = function(data, a, r) {
	sum(data[,2] == a & data[, 2+r]==1)
}

# construct recipient covariates X
D = length(edge)
A = length(Montgomery$manager_gender)
P = 11
X = array(0, dim = c(D,A,A,P))
X[,,,1] = 1
timeunit = 3600
for (d in 2:D) {
	index = which(email$timepoints >= uniqtime[which(uniqtime==email$timepoints[d])-1]-7*24*timeunit & email$timepoints < email$timepoints[d])
	data = email[index, ]
	sent = data[, 2]
	received = data[, 3:(2+A)]
	outdegree = tabulate(sent, A)
	indegree = colSums(received)
	for (a in 1:A) {
		for (r in c(1:A)[-a]) {
			X[d, a, r, 2] = outdegree[a]  
			X[d, a, r, 3] = indegree[r]	
			X[d, a, r, 4] = sendraw(data, a, r)
			X[d, a, r, 5] = sendraw(data, r, a)
			X[d, a, r, 6] = sum(sapply(c(1:A)[-c(a,r)], function(h) {
				sendraw(data, a, h) * sendraw(data, h, r) / 10
				}))
			X[d, a, r, 7] = sum(sapply(c(1:A)[-c(a,r)], function(h) {
				sendraw(data, h, a) * sendraw(data, r, h)
				})) / 10
			X[d, a, r, 8] = sum(sapply(c(1:A)[-c(a,r)], function(h) {
				sendraw(data, h, a) * sendraw(data, h, r)
				})) / 10
			X[d, a, r, 9] = sum(sapply(c(1:A)[-c(a,r)], function(h) {
				sendraw(data, a, h) * sendraw(data, r, h)
				}))	/10	
		}
	  X[d, a, , 10] = ifelse(outdegree[a] > 0, sum(X[d,a,,4]), 0)
	  X[d, a, , 11] = X[d, a, , 2] * X[d, a, , 10] / 10
	}
}

prior.beta = list(mean = c(-3.5, rep(0, P-1)), var = 2*diag(P))
prior.eta = list(mean = c(7, rep(0, Q-1)), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)
email$timepoints =  as.numeric(as.POSIXct(strptime(email[,1], "%d %b %Y %H:%M:%S")))
trim = which(email$timepoints >=7*24*timeunit+email$timepoints[1])
edge = edge[trim]
X = X[trim,,,]
Y = Y[trim,,]

set.seed(1)
missing = list()
missing[[1]] = matrix(0, nrow = dim(Y)[1], 1)
missing[[1]][sample(1:dim(Y)[1], 62, replace = FALSE), ] = 1
missing[[2]] = matrix(0, nrow = dim(Y)[1], A)
missing[[2]][sample(1:(dim(Y)[1]*A), 1118, replace = FALSE)] = 1
missing[[3]] = matrix(0, nrow = dim(Y)[1], 1)
missing[[3]][sample(1:dim(Y)[1], 62, replace = FALSE), ] = 1

senderprob = tabulate(sapply(which(missing[[1]][,1]==0), function(x) edge[[x]]$a_d), 18)

library(MCMCpack)
posterior <- MCmultinomdirichlet(senderprob, rep(1, 18), mc=10000)
senderpred = rmultinom(sum(missing[[1]][,1]), 1, colMeans(posterior))

library(gridExtra)
library(latex2exp)
p = list()
#################
load("/Users/bomin8319/Desktop/MulticastNetwork/Emails/Montgomery_PPE.RData")
load("/Users/bomin8319/Desktop/MulticastNetwork/Emails/Montgomery_PPE2.RData")
truesender = sapply(which(missing[[1]][,1]==1), function(d) edge[[d]]$a_d)
predprob = Montgomery_PPE$senderprob/rowSums(Montgomery_PPE$senderprob)
predprob2 = Montgomery_PPE2$senderprob/rowSums(Montgomery_PPE2$senderprob)
predprob3 = colMeans(posterior)

sender = data.frame(probtrue = sapply(1:62, function(d) predprob[d,truesender[d]]), Dist = rep("log-normal", 62))
sender = rbind(sender, data.frame(probtrue = sapply(1:62, function(d) predprob2[d,truesender[d]]), Dist = rep("exponential", 62)))
sender = rbind(sender, data.frame(probtrue = sapply(1:62, function(d) predprob3[truesender[d]]), Dist = rep("baseline", 62)))

library(ggplot2)
library(reshape)
sender = melt(sender)
colnames(sender)[3] = "correct"
p[[1]]=ggplot(data = sender, aes(x = Dist, y = correct, fill = Dist))+geom_boxplot()+ylab(expression(pi))+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")

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

receiver = data.frame(F1 = probtrue, Dist = rep("log-normal", 500))
receiver = rbind(receiver, data.frame(F1 = probtrue2, Dist = rep("exponential", 500)))
receiver = melt(receiver)
colnames(receiver)[3] = "F1"
p[[2]]=ggplot(data = receiver, aes(x = Dist, y = F1, fill = Dist))+geom_boxplot()+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")


#################################################
truetime = sapply(Montgomery_PPE$timemissing, function(d) edge[[d]]$t_d-edge[[d-1]]$t_d)/3600
predtime = Montgomery_PPE$timepredict
predtime2 = Montgomery_PPE2$timepredict

time = data.frame(MdAPE = sapply(1:62, function(d) median(abs((predtime[d,]-truetime[d])/truetime[d]))), Dist = rep("log-normal", 62))
time = rbind(time, data.frame(MdAPE = sapply(1:62, function(d) median(abs((predtime2[d,]-truetime[d])/truetime[d]))), Dist = rep("exponential", 62)))
time$MdAPE = log(time$MdAPE)
time = melt(time)
stats = boxplot.stats(time$value)$stats

colnames(time)[3] = "MdAPE"
p[[3]]=ggplot(data = time, aes(x = Dist, y = MdAPE, fill = Dist))+geom_boxplot()+theme(legend.position = "bottom")+ylab("log(MdAPE)")+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")

#######################
truetime = sapply(Montgomery_PPE$timemissing, function(d) edge[[d]]$t_d-edge[[d-1]]$t_d)/3600
predtime = Montgomery_PPE$timepredict
predtime2 = Montgomery_PPE2$timepredict

time = data.frame(MdAPE = sapply(1:62, function(d) median(abs((predtime[d,]-truetime[d])/truetime[d]))), Dist = rep("log-normal", 62))
time = rbind(time, data.frame(MdAPE = sapply(1:62, function(d) median(abs((predtime2[d,]-truetime[d])/truetime[d]))), Dist = rep("exponential", 62)))
time = melt(time)
stats = boxplot.stats(time$value)$stats

colnames(time)[3] = "MdAPE"
p[[3]]=ggplot(data = time, aes(x = Dist, y = MdAPE, fill = Dist))+geom_boxplot()+theme(legend.position = "bottom")+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")+scale_y_continuous(trans = "log", breaks = trans_breaks("log", function(x) exp(x)), labels = scientific)



ggplot(data = time, aes(x = MdAPE, fill = dist))+geom_histogram(position = "dodge")+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")

####################################
time = data.frame(MdAE = sapply(1:62, function(d) mean(abs((predtime[d,]-truetime[d])/truetime[d]))), Dist = rep("log-normal", 62))
time = rbind(time, data.frame(MdAE = sapply(1:62, function(d) mean(abs((predtime2[d,]-truetime[d])/truetime[d]))), Dist = rep("exponential", 62)))
time$MdAE = log(time$MdAE)

time = melt(time)
colnames(time)[3] = "MdAE"
ggplot(data = time, aes(x = Dist, y = MdAE, fill = Dist))+geom_boxplot()+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")
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



