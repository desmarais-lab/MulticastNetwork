for (i in 1:17) {
 print(dim(Temporal_Email_Data[[i]]$email_data))
}


#Caldwell email data
Caldwell_email = unique(Temporal_Email_Data[[2]]$email_data)
time = as.numeric(as.POSIXct(strptime(Caldwell_email[,1], "%d %b %Y %H:%M:%S")))
Caldwell_email = Caldwell_email[order(time),]

email = Caldwell_email

edge = list()
initialtime =  as.numeric(as.POSIXct(strptime("01 Jan 2012 00:00:00", "%d %b %Y %H:%M:%S")))
for (d in 1:nrow(email)) {
	t_d = as.numeric(as.POSIXct(strptime(Caldwell_email[d,1], "%d %b %Y %H:%M:%S")))
	edge[[d]] = list(a_d = email[d,2], r_d = as.numeric(email[d,-c(1:2)]), t_d = t_d)
}

library(lubridate)
uniqtime = unique(time[order(time)])

# construct time covariates Y
D = length(edge)
A = length(Caldwell$manager_gender)
Q = 7
Y = array(1, dim = c(D,A,Q))
for (a in 1:A) {
	Y[,a,2] = 1* (Caldwell$manager_gender[a]=="Female")
	Y[,a,3] = 1* (Caldwell$manager_department[a]=="County Manager")
}
timeunit = 3600
Y[1,,6] = rep(as.numeric(wday(as.POSIXct(strptime("01 Jan 2012 00:00:00", "%d %b %Y %H:%M:%S"))) %in% c(1, 7)), A)
Y[1,,7] = rep(pm(as.POSIXct(strptime("01 Jan 2012 00:00:00", "%d %b %Y %H:%M:%S"))), A)
for (d in 2:D) {
	index = which(uniqtime >= uniqtime[which(uniqtime==uniqtime[d])-1]-7*24*timeunit & uniqtime < uniqtime[d])
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
A = length(Caldwell$manager_gender)
P = 7
X = array(0, dim = c(D,A,A,P))
X[,,,1] = 1
timeunit = 3600
for (d in 2:D) {
	index = which(uniqtime  >= uniqtime[which(uniqtime==uniqtime [d])-1]-7*24*timeunit & uniqtime  < uniqtime[d])
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
			}
	  		X[d, a, , 6] = ifelse(outdegree[a] > 0, sum(X[d,a,,4]), 0)
	  		X[d, a, , 7] = X[d, a, , 2] * X[d, a, ,6] / 10
	}
}

trim = which(uniqtime >=7*24*timeunit+uniqtime[1])

Caldwell = list()
Caldwell$email = edge[trim]
Caldwell$X = X[trim,,,]
Caldwell$Y = Y[trim,,]

save(Caldwell, file = "Caldwell.RData")