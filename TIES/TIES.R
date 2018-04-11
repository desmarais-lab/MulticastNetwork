library(readxl)
TIES<-as.data.frame(read_excel("/Users/bomin8319/Desktop/MulticastNetwork/TIES/TIESv4.xls", col_names = TRUE))
time = as.numeric(TIES[,1])
time = sapply(time, function(x) substr(x, 1, nchar(x)-2))

senders = cbind(TIES$sender1, TIES$sender2, TIES$sender3, TIES$sender4, TIES$sender5)
senders[which(is.na(senders[,1])),1] = 2000
sender = c(senders)

sender = sender[-which(is.na(sender))] #106 unique senders
receiver = TIES$targetstate #169 unique targets

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
TIES <- cbind(TIES$caseid, TIES$sender1, TIES$sender2, TIES$sender3, TIES$sender4, TIES$sender5, TIES$targetstate)
TIES[,1] = sapply(TIES[,1], function(x) substr(x, 1, nchar(x)-2))
TIES[which(is.na(TIES[,2])),2] = rep(2000, sum(is.na(TIES[,2])))

#delete edges with missing senders (zero known sender) -> discontinuity issue?
#TIES = TIES[-which(is.na(TIES[,2])),]
TIES[209, 1] = "19650630"
TIES = as.data.frame(TIES)
TIES[,1] = as.character(TIES[,1])
TIES$time = vapply(TIES[,1], function(x) as.numeric(strptime(x, format = "%Y%m%d")), c(1))

TIES_reduced = TIES[-which(is.na(TIES$time)),]
save(TIES_reduced, file = "TIES_reduced.RData")
#impute 113 cases for missing dates. 
#day missing, plug in 01 unless previous event has same month 
#day missing, plug in 30 or 31 if previous event has same month 
#month & day missing, plug in next month (from previous) with 01
# for (d in 1:nrow(TIES)) {
	# if (is.na(TIES$time[d])) {
		# browser()
	# }
# }
# d = 37
# TIES[d, 1] = "19490401"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 94
# TIES[d, 1] = "19540501"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 100
# TIES[d, 1] = "19550630"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 133
# TIES[d, 1] = "19581001"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 143
# TIES[d, 1] = "19600331"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 148
# TIES[d, 1] = "19600731"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 155
# TIES[d, 1] = "19610601"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 190
# TIES[d, 1] = "19630401"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 191
# TIES[d, 1] = "19630801"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 198
# TIES[d, 1] = "19650101"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 199
# TIES[d, 1] = "19650131"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 223
# TIES[d, 1] = "19670531"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 241
# TIES[d, 1] = "19680131"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 284
# TIES[d, 1] = "19721001"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 300
# TIES[d, 1] = "19731031"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 312
# TIES[d, 1] = "19741101"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 316
# TIES[d, 1] = "19750630"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 317
# TIES[d, 1] = "19750701"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 318
# TIES[d, 1] = "19751201"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 319
# TIES[d, 1] = "19760101"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 321
# TIES[d, 1] = "19760801"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 325
# TIES[d, 1] = "19770101"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 326
# TIES[d, 1] = "19770115"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 343
# TIES[d, 1] = "19771225"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))
# d = 344
# TIES[d, 1] = "19771231"
# TIES$time[d] = as.numeric(strptime(TIES[d,1], format = "%Y%m%d"))


TIES = TIES_reduced
for (i in 1:7) {
	TIES[,i] = as.character(TIES[,i])
}
edge = list()
init = 1
for (d in unique(TIES$time)) {
	data = TIES[which(TIES$time==d),]
	for (m in 1:nrow(data)) {
	sender = as.numeric(data[m,2:6])
	if (sum(is.na(sender)) > 0) {
	sender = sender[-which(is.na(sender))]
	} 
	for (s in 1:length(sender)) {
		hi = node[node[,2] == sender[s],1]
		sender[s] = node[node[,2] == sender[s],1]
	}
	target = node[node[,2]==data[m,7],1]
	t_d = data[m,8]
	edge[[init]] = list(a_d = sender, r_d = target, t_d = t_d)
	init = init+1	
	}
}
save(edge, file = "TIES_edge.RData")
save(node, file = "TIES_node.RData")

timediff = sapply(2:length(edge), function(d) edge[[d]]$t_d - edge[[d-1]]$t_d)
tieevents = edge[which(timediff==0)]


#112 either day missing or month missing
edge = list()
initial = as.numeric(strptime("1945-01-01", format = "%Y-%m-%d"))
for (d in 1:nrow(TIES)) {
	sender = as.numeric(TIES[d,2:6])
	if (sum(is.na(sender)) > 0) {
	sender = sender[-which(is.na(sender))]
	}
	for (s in 1:length(sender)) {
		sender[s] = node[node[,2] == sender[s],1]
	}
	target = node[node[,2]==TIES[d,7],1]
	t_d = TIES[d,8]
	# if (d >= 2 && t_d == edge[[d-1]]$t_d) {
		# browser()
		# t_d = mean(t_d, TIES[d+1,8])
	# }
	edge[[d]] = list(a_d = sender, r_d = target, t_d = t_d)
}
timediff = sapply(2:length(edge), function(d) edge[[d]]$t_d - edge[[d-1]]$t_d)