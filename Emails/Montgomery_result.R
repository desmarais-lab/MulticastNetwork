#source("/Users/bomin8319/Desktop/MulticastNetwork/code/Multicast.R")
library(MulticastNetwork)
library(lubridate)

#####################################
#inference results
load("/Users/bomin8319/Desktop/MulticastNetwork/Montgomery_infer.RData")
names(Montgomery_infer)

library(anytime)
library(ggplot2)
library(MCMCpack)
library(reshape2)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
number_ticks <- function(n) {function(limits) pretty(limits, n)}
colorsbeta = ggplotColours(13)

colnames(Montgomery_infer$beta) = c("intercept", "outdegree", "indegree", "send", "receive",
"two_send", "two_receive", "sibling", "cosibling", "hyperedge_size", "interaction", "gender_sender", "gender_receiver",
"gender_homophily")
beta = data.frame(t(sapply(1:5000, function(i) Montgomery_infer$beta[8*i,])))
# colnames(beta) = dimnames(Montgomery$X)[[4]][1:10]
# colnames(beta)[6:7] = c("twosend", "tworeceive")
# colnames(beta) = dimnames(Montgomery$X)[[4]][1:10]
beta.est = melt(beta[-1])
colnames(beta.est)[1] = "covariates"
colnames(beta.est)[2] = "b"
beta.est[,1] = as.factor(beta.est[,1])
beta.est$covariates <- factor(beta.est$covariates , levels = c("gender_sender", "gender_receiver",
"gender_homophily", "outdegree", "indegree","hyperedge_size", "interaction", "send", "receive",
"two_send", "two_receive", "sibling", "cosibling" ))

bplot = list()
bplot[[1]] = ggplot(data = beta.est[which(beta.est$covariates %in% c("gender_sender", "gender_receiver",
"gender_homophily")),], aes(x = reorder(covariates,-as.numeric(covariates)) , y = b))+ geom_boxplot(fill=colorsbeta[3:1])+coord_flip()+ geom_hline(yintercept = 0.0, colour = "red", size = 1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(4)),legend.position="none")
bplot[[2]] = ggplot(data = beta.est[which(beta.est$covariates %in% c("outdegree", "indegree","hyperedge_size", "interaction")),], aes(x = reorder(covariates,-as.numeric(covariates)), y = b)) + geom_boxplot(fill=colorsbeta[7:4])+coord_flip()+ geom_hline(yintercept = 0.0, colour = "red", size = 1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(4)),legend.position="none")
bplot[[3]] = ggplot(data = beta.est[which(beta.est$covariates %in% c("send", "receive",
"two_send", "two_receive", "sibling", "cosibling")),], aes(x = reorder(covariates,-as.numeric(covariates)), y = b)) + geom_boxplot(fill=colorsbeta[13:8])+coord_flip()+ geom_hline(yintercept = 0.0, colour = "red", size =1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(4)),legend.position="none")
grid.arrange(bplot[[1]], bplot[[2]],bplot[[3]], ncol = 1, nrow = 3)

ggplot(data = beta.est, aes(x = reorder(covariates,-as.numeric(covariates)) , y = b))+ geom_boxplot(fill=colorsbeta[13:1])+coord_flip()+ geom_hline(yintercept = 0.0, colour = "red", size = 1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(4.25)),legend.position="none")

bplot = list()
bplot[[1]] = ggplot(data = beta.est[which(beta.est$covariates %in% c("gender_sender", "gender_receiver",
"gender_homophily", "outdegree", "indegree","hyperedge_size", "interaction")),], aes(x = reorder(covariates,-as.numeric(covariates)) , y = b))+ geom_boxplot(fill=colorsbeta[7:1])+coord_flip()+ geom_hline(yintercept = 0.0, colour = "blue", size = 1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")
bplot[[2]] = ggplot(data = beta.est[which(beta.est$covariates %in% c("send", "receive",
"two_send", "two_receive", "sibling", "cosibling")),], aes(x = reorder(covariates,-as.numeric(covariates)), y = b)) + geom_boxplot(fill=colorsbeta[13:8])+coord_flip()+ geom_hline(yintercept = 0.0, colour = "blue", size =1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")
grid.arrange(bplot[[1]], bplot[[2]], ncol =2, nrow = 1)

##################################################
colnames(Montgomery_infer$eta) = c("intercept","gender", "manager", "outdegree", "indegree", "weekend", "PM")
eta = data.frame(t(sapply(1:5000, function(i) Montgomery_infer$eta[8*i,])))
eta.est = melt(eta[-1])
colnames(eta.est)[1] = "covariates"
colnames(eta.est)[2] = "eta"
eta.est[,1] = as.factor(eta.est[,1])
ggplot(data = eta.est, aes(x = reorder(covariates,-as.numeric(covariates)), y = eta, fill =covariates)) + geom_boxplot()+coord_flip()+ geom_hline(yintercept = 0.0, colour = "blue", size =1, linetype = "dashed")+labs(x = NULL, fill = "Covariates")+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")+ylab(expression(eta))

#####################################
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


indegreedist1 = matrix(NA, 500, A)
outdegreedist1 = matrix(NA,  500, A)
recipientsdist1 = matrix(NA,  500, A-1)
timedist1 = matrix(NA, 500, 621)
#setwd("/Users/bomin8319/Desktop/MulticastNetwork/Emails/PPC")
#setwd("/Users/bomin8319/")
setwd("/Users/bomin8319/Desktop/MulticastNetwork/Emails/PPC4")

for (n in 1:500) {
	filename = paste0("Montgomery_PPCnew", n,".RData")
	load(filename)
	outdegreedist1[n, ] = tabulate(vapply(1:621, function(x) Montgomery_PPC[[x]]$a_d, c(1)), A)
	indegreedist1[n, ] = rowSums(sapply(1:621, function(x) Montgomery_PPC[[x]]$r_d))
	recipientsdist1[n, ] = tabulate(vapply(1:621, function(x) sum(Montgomery_PPC[[x]]$r_d), c(1)), A-1)
	timedist1[n, ] = vapply(1:621, function(x) Montgomery_PPC[[x]]$t_d, c(1))/ timeunit
	} 
###
outdegreedist1 = data.frame(outdegreedist1)
colnames(outdegreedist1)=1:18
library(ggplot2)
library(reshape)
data = melt(outdegreedist1)
colnames(data) = c("Node", "Outdegree")

outdegree = data.frame(Node = 1:18, Outdegree = outdegree)
ggplot(data = data, aes(x = Node, y = Outdegree, fill = as.factor("hi"))) +geom_boxplot() + geom_line(data = outdegree, colour = "blue", size =0.5)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")
###
indegreedist1 = data.frame(indegreedist1)
colnames(indegreedist1)=1:18
library(ggplot2)
library(reshape)
data = melt(indegreedist1)
colnames(data) = c("Node", "Indegree")

indegree = data.frame(Node = 1:18, Indegree = indegree)
ggplot(data = data, aes(x = Node, y = Indegree, fill = as.factor("hi"))) +geom_boxplot() + geom_line(data = indegree, colour = "blue", size =0.5)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")
###
recipientsdist1 = data.frame(recipientsdist1[,1:14])
colnames(recipientsdist1)=1:14
library(ggplot2)
library(reshape)
data = melt(recipientsdist1)
colnames(data) = c("RecipientSize", "Documents")

recipients = data.frame(RecipientSize = 1:14, Documents = recipients[1:14])
ggplot(data = data, aes(x = RecipientSize, y = Documents, fill = as.factor("hi")) )+geom_boxplot() + geom_line(data = recipients, colour = "blue", size =0.5)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")


uniqueValues = quantile(c(timedist1[,-1], timeinc), seq(0, 1, length = 500))
  qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
 		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(c(timedist1[,-1]) <= uniqueValues[j])
  		qx2[j] = mean(c(timeinc) <= uniqueValues[j])
}

time = data.frame(Simulated = qx1, Observed = qx2)
ggplot(data = time, aes(x = Simulated, y = Observed, colour = as.factor("hi"))) + geom_point() + geom_abline(intercept = 0, slope = 1,colour = "blue", size =0.5)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")

##########################################
#################################

indegreedist2 = matrix(NA, 500, A)
outdegreedist2 = matrix(NA,  500, A)
recipientsdist2 = matrix(NA,  500, A-1)
timedist2 = matrix(NA, 500, 621)
#setwd("/Users/bomin8319/Desktop/MulticastNetwork/Emails/PPC")
#setwd("/Users/bomin8319/")
setwd("/Users/bomin8319/Desktop/MulticastNetwork/Emails/PPC3")

for (n in 1:500) {
	filename = paste0("Montgomery_PPCnew2", n,".RData")
	load(filename)
	outdegreedist2[n, ] = tabulate(vapply(1:621, function(x) Montgomery_PPC2[[x]]$a_d, c(1)), A)
	indegreedist2[n, ] = rowSums(sapply(1:621, function(x) Montgomery_PPC2[[x]]$r_d))
	recipientsdist2[n, ] = tabulate(vapply(1:621, function(x) sum(Montgomery_PPC2[[x]]$r_d), c(1)), A-1)
	timedist2[n, ] = vapply(1:621, function(x) Montgomery_PPC2[[x]]$t_d, c(1))/ timeunit
} 

outdegreedist = data.frame(outdegreedist)
colnames(outdegreedist)=1:18
library(ggplot2)
library(reshape)
data = melt(outdegreedist)
colnames(data) = c("Node", "Outdegree")

outdegree = data.frame(Node = 1:18, Outdegree = outdegree)
ggplot(data = data, aes(x = Node, y = Outdegree)) +geom_boxplot() + geom_line(data = outdegree, col = 'red')

indegreedist = data.frame(indegreedist)
colnames(indegreedist)=1:18
library(ggplot2)
library(reshape)
data = melt(indegreedist)
colnames(data) = c("Node", "Indegree")

indegree = data.frame(Node = 1:18, Indegree = indegree)
ggplot(data = data, aes(x = Node, y = Indegree)) +geom_boxplot() + geom_line(data = indegree, col = 'red')

recipientsdist = data.frame(recipientsdist[,1:14])
colnames(recipientsdist)=1:14
library(ggplot2)
library(reshape)
data = melt(recipientsdist)
colnames(data) = c("RecipientSize", "Documents")

recipients = data.frame(RecipientSize = 1:14, Documents = recipients[1:14])
ggplot(data = data, aes(x = RecipientSize, y = Documents) )+geom_boxplot() + geom_line(data = recipients, col = 'red')

uniqueValues = quantile(c(timedist[,-1], timeinc), seq(0, 1, length = 500))
  qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
 		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(c(timedist[,-1]) <= uniqueValues[j])
  		qx2[j] = mean(c(timeinc) <= uniqueValues[j])
}

time = data.frame(Simulated = qx1, Observed = qx2)
ggplot(data = time, aes(x = Simulated, y = Observed)) + geom_point() + geom_abline(intercept = 0, slope = 1, col = 'red')
###################################

outdegreedist = data.frame(outdegreedist1, model = rep("lognormal", 500))
outdegreedist = rbind(outdegreedist, data.frame(outdegreedist2, model = rep("exponential", 500)))
colnames(outdegreedist)[1:18]=1:18
library(ggplot2)
library(reshape)
data = melt(outdegreedist)
colnames(data) = c("Model","Node", "Outdegree")

outdegreeobs = data.frame(Node = 1:18, Outdegree = outdegree, Model = rep("lognormal", 18))
ggplot(data = data, aes(x = Node, y = Outdegree, fill= Model)) +geom_boxplot() + geom_line(data = outdegreeobs, colour = "blue", size =0.5, group = 1)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")

indegreedist = data.frame(indegreedist1, model = rep("lognormal", 500))
indegreedist = rbind(indegreedist, data.frame(indegreedist2, model = rep("exponential", 500)))
colnames(indegreedist)[1:18]=1:18
data = melt(indegreedist)
colnames(data) = c("Model","Node", "Indegree")

indegreeobs = data.frame(Node = 1:18, Indegree = indegree, Model = rep("lognormal", 18))
ggplot(data = data, aes(x = Node, y =Indegree, fill= Model)) +geom_boxplot() + geom_line(data = indegreeobs, colour = "blue", size =0.5, group = 1)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")

recipientsdist = data.frame(recipientsdist1[,1:14], model = rep("lognormal", 500))
recipientsdist = rbind(recipientsdist, data.frame(recipientsdist2[,1:14], model = rep("exponential", 500)))
colnames(recipientsdist)[1:14]=1:14
data = melt(recipientsdist)
colnames(data) = c("Model","RecipientSize","Documents")

recipientsobs = data.frame(RecipientSize = 1:14, Documents = recipients[1:14], Model = rep("lognormal", 14))
ggplot(data = data, aes(x =RecipientSize, y = Documents, fill= Model)) +geom_boxplot()+ geom_line(data = recipientsobs, colour = "blue", size =0.5, group = 1)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none") 




uniqueValues = quantile(c(timedist1[,-1], timeinc), seq(0, 1, length = 500))
  qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
 		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(c(timedist1[,-1]) <= uniqueValues[j])
  		qx2[j] = mean(c(timeinc) <= uniqueValues[j])
}


time = data.frame(Simulated = qx1, Observed = qx2, Model = rep("log-normal", 500))

uniqueValues = quantile(c(timedist2[,-1], timeinc), seq(0, 1, length = 500))
  qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
 		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(c(timedist2[,-1]) <= uniqueValues[j])
  		qx2[j] = mean(c(timeinc) <= uniqueValues[j])
}

time = rbind(time, data.frame(Simulated = qx1, Observed = qx2, Model = rep("exponential", 500)))
ggplot(data = time, aes(x = Simulated, y = Observed, colour = Model)) + geom_point()+ geom_abline(intercept = 0, slope = 1,colour = "blue", size =0.5)+theme(plot.title = element_blank(),text = element_text(size = rel(5.5)),legend.position="none")



#################################
set.seed(1)
missing = list()
missing[[1]] = matrix(0, nrow = 621, 1)
missing[[1]][sample(1:621[1], 62, replace = FALSE), ] = 1
missing[[2]] = matrix(0, nrow =621, A)
missing[[2]][sample(1:(621*A), 1118, replace = FALSE)] = 1
missing[[3]] = matrix(0, nrow = 621, 1)
missing[[3]][sample(1:621, 62, replace = FALSE), ] = 1

truesenders = sapply(trim[which(missing[[1]]==1)], function(x) edge[[x]]$a_d)
truereceivers = sapply(trim[which(rowSums(missing[[2]])>=1)], function(x) edge[[x]]$r_d[which(missing[[2]][x-42,]==1)])
truetime = sapply(trim[which(missing[[3]]==1)], function(x) (edge[[x]]$t_d - edge[[x-1]]$t_d)/3600)

setwd("/Users/bomin8319/Desktop/MulticastNetwork/Emails/PPE")
sender = matrix(NA, nrow = 500, ncol = 62)
sender2 = matrix(NA, nrow = 500, ncol = 62)
time = matrix(NA, nrow = 500, ncol = 62)
time2 = matrix(NA, nrow = 500, ncol = 62)
for (n in 1:500) {
	filename = paste0("Montgomery_PPE", n,".RData")
	load(filename)
	sender[n, ] = sapply(1:62, function(x) sum(Montgomery_PPE$senderpredict[x,]==truesenders[x])/50)
	time[n, ] = sapply(1:62, function(x) sqrt(mean((Montgomery_PPE$timepredict[x,]-truetime[x])^2)))
    filename = paste0("Montgomery_PPE2", n,"exp.RData")
    load(filename)
    sender2[n, ] = sapply(1:62, function(x) sum(Montgomery_PPE2$senderpredict[x,]==truesenders[x])/50)
    time2[n, ] = sapply(1:62, function(x) sqrt(mean((Montgomery_PPE2$timepredict[x,]-truetime[x])^2)))
} 

data = data.frame(melt(sender))
data = rbind(data, melt(sender2))
data$dist =c(rep("lognormal", dim(data)[1]/2), rep("exponential", dim(data)[1]/2))
data = data[,-1]
colnames(data) = c("Document", "Precision", "Dist")
data[,1] = as.factor(data[,1])
data[,3] = as.factor(data[,3])
ggplot(data = data, aes(x = Document, y = Precision, fill = Dist))+geom_boxplot()


data = data.frame(melt(time))
data = rbind(data, melt(time2))
data$dist =c(rep("lognormal", dim(data)[1]/2), rep("exponential", dim(data)[1]/2))
data = data[,-1]
colnames(data) = c("Document", "RMSE", "Dist")
data[,1] = as.factor(data[,1])
data[,3] = as.factor(data[,3])
ggplot(data = data, aes(x = Document, y = RMSE, fill = Dist))+geom_boxplot()+coord_cartesian(ylim = c(0,1.0e+05))

