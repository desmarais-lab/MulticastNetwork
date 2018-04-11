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
A = length(Montgomery$manager_gender)

#observed statistics
timeunit = 3600
trim = which(email$timepoints >=3*24*timeunit+email$timepoints[1])

indegree = tabulate(email[trim,2], A)
outdegree = colSums(email[trim,3:20])
recipients = tabulate(rowSums(email[trim,3:20]), A-1)
timeinc = diff(email$timepoints[21:max(trim)])

indegreedist = matrix(NA, 50, A)
outdegreedist = matrix(NA, 50, A)
recipientsdist = matrix(NA, 50, A-1)
setwd("/Users/bomin8319/Desktop/MulticastNetwork/code/PPC")
for (n in 1:50) {
	filename = paste0("Montgomery_PPC", n,".RData")
	load(filename)
	indegreedist[n, ] = tabulate(vapply(1:642, function(x) Montgomery_PPC[[x]]$a_d, c(1)), A)
	outdegreedist[n, ] = rowSums(sapply(1:642, function(x) Montgomery_PPC[[x]]$r_d))
	recipientsdist[n, ] = tabulate(vapply(1:642, function(x) sum(Montgomery_PPC[[x]]$r_d), c(1)), A-1)
} 

par(mfrow = c(3,1))

boxplot(indegreedist)
lines(indegree, col = 2)

boxplot(outdegreedist)
lines(outdegree, col = 2)

boxplot(recipientsdist)
lines(recipients, col = 2)



load("/Users/bomin8319/Desktop/MulticastNetwork/code/Montgomery_infer.RData")
initial = list()
initial$sender = email[1:21, 2]
initial$receiver = email[1:21, 3:20]
initial$time = email[1:21,1]
