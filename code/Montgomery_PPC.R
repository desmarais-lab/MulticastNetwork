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
trim = which(email$timepoints >=7*24*timeunit+email$timepoints[1])

outdegree = tabulate(email[trim,2], A)
indegree = colSums(email[trim,3:20])
recipients = tabulate(rowSums(email[trim,3:20]), A-1)
timeinc = diff(sort(email$timepoints)[43:max(trim)])/3600


indegreedist = matrix(NA, 500, A)
outdegreedist = matrix(NA,  500, A)
recipientsdist = matrix(NA,  500, A-1)
timedist = matrix(NA, 500, 621)
setwd("/Users/bomin8319/Desktop/MulticastNetwork/code/PPC")
for (n in 1:500) {
	filename = paste0("Montgomery_PPCnew", n,".RData")
	load(filename)
	outdegreedist[n, ] = tabulate(vapply(1:621, function(x) Montgomery_PPC[[x]]$a_d, c(1)), A)
	indegreedist[n, ] = rowSums(sapply(1:621, function(x) Montgomery_PPC[[x]]$r_d))
	recipientsdist[n, ] = tabulate(vapply(1:621, function(x) sum(Montgomery_PPC[[x]]$r_d), c(1)), A-1)
	timedist[n, ] = c(email$timepoints[42], diff(vapply(1:621, function(x) sum(Montgomery_PPC[[x]]$t_d), c(1)))) / timeunit
} 

par(mfrow = c(3,1))

# indegreesum = table(floor(indegreedist/10))
# boxplot(floor(indegreedist/10))

boxplot(outdegreedist, ylim = c(0, 175), main = "outdegree")
lines(outdegree, col = 2)

boxplot(indegreedist,ylim = c(0, 275),  main = "indegree")
lines(indegree, col = 2)

boxplot(recipientsdist, ylim = c(0, 515),  main = "receiver size")
lines(recipients, col = 2)


hi = quantile(c(timedist[,-1]), c(.025, .975 ))
qqplot(c(timedist[,-1])[c(timedist[,-1])>=hi[1] & c(timedist[,-1])<=hi[2]], timeinc, xlab = "post", ylab = "obs", main = "timeinc")
abline(0, 1, col = 2)

uniqueValues = quantile(c(timedist[,-1], timeinc), seq(0, 1, length = 1000))
  qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
 		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(c(timedist[,-1]) <= uniqueValues[j])
  		qx2[j] = mean(c(timeinc) <= uniqueValues[j])
  	}


GiR_PP_Plots = function(Forward_stats, Backward_stats) {
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
    
    qqplot(x = qx1,
           y = qx2,
           ylim = c(0, 1),
           xlim = c(0, 1),
           ylab = "Backward",
           xlab = "Forward",
           col = "blue",
           pch = 19,
           cex = 0.25,
           main = nms[i],
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    
      if (nrow(Forward_stats) > 10000) {
       thinning2 = seq(from = floor(nrow(Forward_stats) / 10), to = nrow(Forward_stats), length.out = 10000)
       Forward_test2 = Forward_stats[thinning2, i]
       Backward_test2 = Backward_stats[thinning2, i]
       } else {
        Forward_test2 = Forward_stats[, i]
        Backward_test2 = Backward_stats[, i]    	
      }
    text(paste("Backward Mean:", round(mean(Backward_stats[, i]), 4),
                "\nForward Mean:", round(mean(Forward_stats[, i]), 4),
                "\nt-test p-value:", round(t.test(Backward_test2, Forward_test2)$p.value, 4),
                "\nMann-Whitney p-value:", round(wilcox.test(Backward_test2, Forward_test2)$p.value, 4)),
                x = 0.65, y = 0.15, cex = 0.4)
  }
}      


load("/Users/bomin8319/Desktop/MulticastNetwork/code/Montgomery_infer.RData")
initial = list()
initial$sender = email[1:21, 2]
initial$receiver = email[1:21, 3:20]
initial$time = email[1:21,1]
