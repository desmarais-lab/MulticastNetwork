source("/Users/bomin8319/Desktop/MulticastNetwork/code/Multicast.R")
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
	index = which(email$timepoints >= email$timepoints[d-1]-7*24*timeunit & email$timepoints < email$timepoints[d])
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
P = 9
X = array(0, dim = c(D,A,A,P))
X[,,,1] = 1
timeunit = 3600
for (d in 2:D) {
	index = which(email$timepoints >= email$timepoints[d-1]-7*24*timeunit & email$timepoints < email$timepoints[d])
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
	}
}

prior.beta = list(mean = c(-3, rep(0, P-1)), var = 2*diag(P))
prior.eta = list(mean = c(5, rep(0, Q-1)), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)
email$timepoints =  as.numeric(as.POSIXct(strptime(email[,1], "%d %b %Y %H:%M:%S")))
trim = which(email$timepoints >=7*24*timeunit+email$timepoints[1])
edge = edge[trim]
X = X[trim,,,]
Y = Y[trim,,]
Montgomery_infer = Inference(edge, X, Y, 55000, c(5,1,1), 15000, prior.beta, prior.eta, prior.sigma2, initial = NULL,
		  proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, lasttime = email[min(trim-1), 21] - initialtime)

save(Montgomery_infer, file = "Montgomery_infer.RData")

load("/Users/bomin8319/Desktop/MulticastNetwork/code/Montgomery_infer.RData")
initial = list()
initial$sender = email[1:(min(trim)-1), 2]
initial$receiver = email[1:(min(trim)-1), 3:20]
initial$time = email[1:(min(trim)-1),1]
for (n in 1:100) {
  Montgomery_PPC = PPC(length(edge), A, colMeans(Montgomery_infer$beta), colMeans(Montgomery_infer$eta), 
                       mean(Montgomery_infer$sigma2), X, Y, timeunit = 3600, lasttime = email[min(trim-1), 21], 
                       Montgomery_infer$u, initial =initial)
  filename = paste0("Montgomery_PPC", n,".RData")
  save(Montgomery_PPC, file = filename)
}


