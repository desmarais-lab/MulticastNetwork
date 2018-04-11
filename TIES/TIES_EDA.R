library(readxl)
TIES_raw <-as.data.frame(read_excel("/Users/bomin8319/Desktop/MulticastNetwork/TIES/TIESv4.xls", col_names = TRUE))


######Below generates node (Line 6 - 22)#####
senders = cbind(TIES_raw$sender1, TIES_raw$sender2, TIES_raw$sender3, TIES_raw$sender4, TIES_raw$sender5)
senders[which(is.na(senders[,1])),1] = 2000
sender = c(senders)
sender = sender[-which(is.na(sender))] #106 unique senders
receiver = TIES_raw$targetstate #169 unique targets
node = sort(as.numeric(unique(union(sender, receiver)))) #176 unique nodes 
COW = read.csv("/Users/bomin8319/Desktop/MulticastNetwork/TIES/COW country codes.csv")
COW = unique(COW)
node2 = cbind(sort(as.numeric(node[-c(175, 176)])), COW[which(COW$CCode %in% node), c(1,3)])
node2 = as.matrix(node2)
node2 = rbind(node2, c(1000, "EEC/EU", "European Economic Community/European Union"))
node2 = rbind(node2, c(2000, "II", "International Institution"))
node = as.data.frame(node2)
colnames(node)[1] = "COWcode"
rownames(node) = NULL
node = data.frame(ID = 1:176, node)
node[,2] =  as.numeric(as.character(node[,2]))
############################################
head(node)


#####Below generates edges (Line 27 - 35)########
TIES <- cbind(TIES_raw$caseid, TIES_raw$sender1, TIES_raw$sender2, TIES_raw$sender3, TIES_raw$sender4, TIES_raw$sender5, TIES_raw$targetstate)
TIES[,1] = sapply(TIES[,1], function(x) substr(x, 1, nchar(x)-2))
TIES[which(is.na(TIES[,2])),2] = rep(2000, sum(is.na(TIES[,2])))
TIES[209, 1] = "19650630"  #one weird date which has June 31th -> fix to 30th
TIES = as.data.frame(TIES)
TIES[,1] = as.character(TIES[,1])
TIES$time = vapply(TIES[,1], function(x) as.numeric(strptime(x, format = "%Y%m%d")), c(1))
TIES = TIES[-which(is.na(TIES$time)),] #delete unknown dates
colnames(TIES) = c("date", "sender1", "sender2", "sender3", "sender4", "sender5", "target", "unixtime")
##################################################
head(TIES)


#####Below generate edges that happened at the same date#####
TIES_redundant = matrix(NA, nrow = 0, ncol = 8)
for (d in unique(TIES$unixtime)) {
	if (length(which(TIES$unixtime == d)) > 1) {
		TIES_redundant = rbind(TIES_redundant, TIES[which(TIES$unixtime == d),])
	}
}
#########################################################
head(TIES_redundant)