# only common mean
nDocs = 50
s2 = 1
node = 10
X = matrix(1, nDocs, node)
	eta.true = rnorm(1, 0, 1)
	mu = eta.true[1] *  X
	obs = rep(NA, nDocs)
	sender = rep(NA, nDocs)
	for (d in 1:nDocs) {
		times = rlnorm(node, mu[d,], s2)
		obs[d] = min(times)
		sender[d] = which(times == obs[d])
	}
post.old = function(eta.old, obs, sender) {
	 mu = eta.old * X
	 sum(sapply(1:nDocs, function(d) dlnorm(obs[d], mu[d,sender[d]], s2, log = TRUE) + sum(sapply(c(1:node)[-sender[d]], function(x) plnorm(obs[d], mu[d, x], s2, lower.tail = FALSE, log.p = TRUE)))))
	}
eta.grid = seq(-10, 10, by = 0.01)
y = sapply(eta.grid, function(x) post.old(x, obs, sender))
plot(eta.grid, y)
c(eta.true,eta.grid[which(y==max(y))])

#mean --- test of s2
nDocs = 100
s2 = rhalfcauchy(1, 0.5)
node = 4
X = matrix(1, nDocs, node)
	eta.true = rnorm(1, 0, 1)
	mu = eta.true[1]* X 
	obs = rep(NA, nDocs)
	sender = rep(NA, nDocs)
	for (d in 1:nDocs) {
		times = rlnorm(node, mu[d,], s2)
		obs[d] = min(times)
		sender[d] = which(times == obs[d])
	}
post.old = function(eta.old, obs, sender, s2) {
	 mu = eta.old[1]* X
	sum(sapply(1:nDocs, function(d) dlnorm(obs[d], mu[d,sender[d]], s2, log = TRUE) + sum(sapply(c(1:node)[-sender[d]], function(x) plnorm(obs[d], mu[d, x], s2, lower.tail = FALSE, log.p = TRUE)))))
	}
	
s2.grid = seq(0.1, 10, by = 0.01)
y = sapply(s2.grid, function(x) post.old(eta.true, obs, sender, x))
plot(s2.grid, y)
c(s2,s2.grid[which(y==max(y))])



#mean + one covariate
nDocs = 10
s2 = 1
node = 4
X = matrix(sample(1:5, nDocs * node, replace = TRUE), nDocs)
	eta.true = rnorm(2, 0, 1)
	mu = eta.true[1] + X * eta.true[2]
	obs = rep(NA, nDocs)
	sender = rep(NA, nDocs)
	for (d in 1:nDocs) {
		times = rlnorm(node, mu[d,], s2)
		obs[d] = min(times)
		sender[d] = which(times == obs[d])
	}
post.old = function(eta.old, obs, sender) {
	 mu = eta.old[1] + X * eta.old[2]
	 sum(sapply(1:nDocs, function(d) dlnorm(obs[d], mu[d,sender[d]], s2, log = TRUE) + sum(sapply(c(1:node)[-sender[d]], function(x) plnorm(obs[d], mu[d, x], s2, lower.tail = FALSE, log.p = TRUE)))))
	}
eta.grid = seq(-10, 10, by = 0.01)
y = sapply(eta.grid, function(x) post.old(c(eta.true[1],x), obs, sender))
plot(eta.grid, y)
c(eta.true[2],eta.grid[which(y==max(y))])


#mean + one covariate --- test of s2
nDocs = 100
s2 = rhalfcauchy(1, 0.1)
node = 4
X = matrix(sample(0:3, nDocs * node, replace = TRUE), nDocs)
	eta.true = rnorm(2, 0, 1)
	mu = eta.true[1] + X * eta.true[2]
	obs = rep(NA, nDocs)
	sender = rep(NA, nDocs)
	for (d in 1:nDocs) {
		times = rlnorm(node, mu[d,], s2)
		obs[d] = min(times)
		sender[d] = which(times == obs[d])
	}
post.old = function(eta.old, obs, sender, s2) {
	 mu = eta.old[1] + X * eta.old[2]
	sum(sapply(1:nDocs, function(d) dlnorm(obs[d], mu[d,sender[d]], s2, log = TRUE) + sum(sapply(c(1:node)[-sender[d]], function(x) plnorm(obs[d], mu[d, x], s2, lower.tail = FALSE, log.p = TRUE)))))
	}
	
s2.grid = seq(0.1, 10, by = 0.01)
y = sapply(s2.grid, function(x) post.old(eta.true, obs, sender, x))
plot(s2.grid, y)
c(s2,s2.grid[which(y==max(y))])


# node-specific intercepts
nDocs =10
s2 = 1
node = 5
X = matrix(1, nDocs, node)
	eta.true = rnorm(node, 0, 1)
	mu = matrix(eta.true, nDocs, node, byrow = TRUE)
	obs = rep(NA, nDocs)
	sender = rep(NA, nDocs)
	for (d in 1:nDocs) {
		times = rlnorm(node, mu[d,], s2)
		obs[d] = min(times)
		sender[d] = which(times == obs[d])
	}

post.old = function(eta.old, obs, sender) {
	 mu = matrix(eta.old, nDocs, node, byrow = TRUE)
	 sum(sapply(1:nDocs, function(d) dlnorm(obs[d], mu[d,sender[d]], s2, log = TRUE) + sum( plnorm(obs[d], mu[d,-sender[d]], s2, lower.tail = FALSE, log.p = TRUE))))
	}
eta.grid = seq(-10, 10, by = 0.01)
y = sapply(eta.grid, function(x) post.old(c(eta.true[1:4],x), obs, sender))
plot(eta.grid, y)
c(eta.true[5],eta.grid[which(y==max(y))])






#MLE using log-normal

log.lik <- function(par,val,index){
    mu <- par[1:3]
    sig <- par[4]
    ll <- 0
    for(i in 1:length(index)){
        ll <- ll + dlnorm(val[i],mu[index[i]],sig, log = TRUE) + sum(plnorm(val[i],mu[-index[i]],sig, FALSE, TRUE))
    }
    -ll
}

y <- cbind(rlnorm(1000,0,5),rlnorm(1000,1,5),rlnorm(1000,2,5))

vals <- apply(y,1,min)
index <- apply(y,1,which.min)

optim(par = rep(1,4),fn=log.lik,val=vals,index=index)

s2grid = seq(3, 10, by = 0.1)
y = sapply(s2grid, function(x) log.lik(c(0,1,2,x), vals, index))
plot(s2grid, y)
