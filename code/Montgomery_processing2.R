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
P = 14
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
		for (r in c(1:A)) {
		  if (r != a) {
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
		    X[d, a, , 12] = 1* (Montgomery$manager_gender[a]=="Female")
			  X[d, a, r, 13] =1* (Montgomery$manager_gender[r]=="Female")
			  X[d, a, r, 14] =1* (Montgomery$manager_gender[a]==Montgomery$manager_gender[r])
		}
	  X[d, a, , 10] = ifelse(outdegree[a] <sum(X[d,a,,4]),sum(X[d,a,,4])/outdegree[a] , 1)
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

prior.beta = list(mean = c(-3.5, rep(0, P-1)), var = 1*diag(P))
prior.eta = list(mean = c(7, rep(0, Q-1)), var = 1*diag(Q))
Montgomery_infer = Inference(edge, X, Y, 55000, c(20,10,1), 15000, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
                             proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = email[min(trim-1), 21] - initialtime, timedist = "lognormal")
save(Montgomery_infer, file= "/Users/bomin8319/Desktop/Montgomery_infer.RData")

prior.eta = list(mean = c(7, rep(0, Q-1)), var = 4*diag(Q))
Montgomery_infer4 = Inference(edge, X, Y, 55000, c(20,10,1), 15000, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
                             proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = email[min(trim-1), 21] - initialtime, timedist = "lognormal")
save(Montgomery_infer4, file= "/Users/bomin8319/Desktop/Montgomery_infer4.RData")

############################################################################
Montgomery_infer2 = Inference(edge, X, Y, 55000, c(20,10,1), 15000, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
                             proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = email[min(trim-1), 21] - initialtime, timedist = "exponential")
save(Montgomery_infer2, file= "/Users/bomin8319/Desktop/Montgomery_infer2.RData")
###############################################


initialval = list()
initialval$u = Montgomery_infer$u
initialval$beta = colMeans(Montgomery_infer$beta)
initialval$eta = colMeans(Montgomery_infer$eta)
initialval$sigma2 = mean(Montgomery_infer$sigma2)

Montgomery_infer_new = Inference(edge, X, Y, 40000, c(40,10,1), 0, prior.beta, prior.eta, prior.sigma2, initialval = initialval,
                             proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = email[min(trim-1), 21] - initialtime, timedist = "lognormal")


load("/Users/bomin8319/Box/gainlab_example/Bomin/Montgomery.RData")
edge = Montgomery$edge
X = Montgomery$X
Y = Montgomery$Y
P = dim(X)[4]
Q = dim(Y)[3]
A = dim(Y)[2]
prior.beta = list(mean = c(-3.5, rep(0, P-1)), var = 2*diag(P))
prior.eta = list(mean = c(7, rep(0, Q-1)), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)
Montgomery_infer = Inference(edge, X, Y, 55000, c(20,10,1), 15000, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
		  proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = email[min(trim-1), 21] - initialtime, timedist = "lognormal")
save(Montgomery_infer, file = "/Users/bomin8319/Desktop/Montgomery_infer.RData")

dimnames(X)[[4]] = c("intercept", "outdegree", "indegree", "send", "receive", "2send", "2receive", "sibling", "cosibling", "hyperedge_size", "outdegree*hyperedge_size")
dimnames(Y)[[3]] = c("intercept", "female", "manager", "outdegree", "indegree", "weekend", "pm")

Montgomery = list(edge = edge, X = X, Y = Y, lasttime = email[min(trim-1), 21] - initialtime )

save(Montgomery, file = "Montgomery.RData")
###################################################

load("/Users/bomin8319/Desktop/Montgomery_infer.RData")
load("/Users/bomin8319/Desktop/MulticastNetwork/Emails/Montgomery_infer.RData")
setwd("/Users/bomin8319/Desktop/PPC")
initial = list()
initial$sender = email[1:(min(trim)-1), 2]
initial$receiver = email[1:(min(trim)-1), 3:20]
initial$time = email[1:(min(trim)-1),1]
for (n in 1:500) {
    Montgomery_PPC = PPC(Montgomery$email_data, length(edge), beta = colMeans(Montgomery_infer$beta), eta = colMeans(Montgomery_infer$eta),
    sigma2 = mean(Montgomery_infer$sigma2), X, Y, timeunit = 3600, u = Montgomery_infer$u, timedist = "lognormal")
  filename = paste0("Montgomery_PPCnew", n,".RData")
  save(Montgomery_PPC, file = filename)
}



Montgomery_infer2 = Inference(edge, X, Y, 55000, c(10,1,1), 15000, prior.beta, prior.eta, prior.sigma2, initialval = NULL,
proposal.var = c(0.00001, 0.001, 0.1), timeunit = 3600, lasttime = Montgomery$lasttime, timedist = "exponential")

save(Montgomery_infer2, file = "/Users/bomin8319/Desktop/Montgomery_infer2.RData")
load("/Users/bomin8319/Desktop/MulticastNetwork/Montgomery_infer2.RData")
setwd("/Users/bomin8319/Desktop/MulticastNetwork/Emails/PPC2")
for (n in 1:500) {
	print(n)
    Montgomery_PPC2 = PPC(length(edge), beta = colMeans(Montgomery_infer2$beta), eta = colMeans(Montgomery_infer2$eta),
    sigma2 = mean(Montgomery_infer2$sigma2), X, Y, timeunit = 3600, u = Montgomery_infer2$u, timedist = "exponential")
  filename = paste0("Montgomery_PPCnew2", n,".RData")
  save(Montgomery_PPC2, file = filename)
}

