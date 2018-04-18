# Generative process
library(combinat)
library(mvtnorm)
library(MCMCpack)
library(Rcpp)
library(lubridate)
sourceCpp('/Users/bomin8319/Desktop/MulticastNetwork/code/Multicast_rcpp.cpp')

send = function(sender, receiver, a, r) {
  D = length(sender)
  sums = 0
  for (d in 1:D) {
    sums = sums + as.numeric(sender[d] == a & receiver[d, r]==1)
  }
  return(sums)
}

gibbs.measure.support = function(n) {
	gibbs.support = rbind(rep(1, n))
	for(i in 1:(n-1)){
		gibbs.mat.i = do.call('rbind',permn(c(rep(1, i), rep(0, n-i))))
		gibbs.support = rbind(gibbs.support, gibbs.mat.i)
	}
	out = as.matrix(unique(gibbs.support))
	return(out)
}

r.gibbs.measure <- function(lambda.i, support) {
	logitNumerator = vapply(1:nrow(support), function(s) {
		sum(lambda.i*support[s,])
		}, c(1))		
	samp = multinom_vec(exp(logitNumerator))
	return(support[samp,])	
}

dinvgamma = function(x, shape, scale) {
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
        1) * log(x) - (beta/x)
    return(log.density)
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

Generate = function(D, A, beta, eta, sigma2, X, Y, support, timeunit = 3600) {
	P = length(beta)
	Q = length(eta)
	u = list()
	data = list()
	t_d = 0
	lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta))
	mu = mu_cpp(Y, eta)
	d = 1
	while (d <= D) {
		u[[d]] = matrix(0, A, A)
		for (a in 1:A) {
			u[[d]][a, -a] = r.gibbs.measure(lambda[[d]][a, -a], support) 
		}
		tau = rlnorm(A, mu[d, ], sqrt(sigma2))
		for (n in 1:length(which(tau == min(tau)))) {
		a_d = which(tau == min(tau))
		r_d = u[[d]][a_d,]
		t_d = t_d + min(tau) * timeunit
		data[[d]] = list(a_d = a_d, r_d = r_d, t_d = t_d)
		d = d+1
		}
	}
	return(list(data = data, u = u, beta = beta, eta = eta, sigma2 = sigma2))
}

PPC = function(D, A, beta, eta, sigma2, X, Y, timeunit = 3600, lasttime, u, initial = NULL) {
  P = length(beta)
  Q = length(eta)
  X = X
  Y = Y
  u = u
  data = list()
  t_d = lasttime
  lambda = list()
  mu = matrix(NA, D, A)
  sender = initial$sender
  receiver = initial$receiver
  time = as.numeric(as.POSIXct(strptime(initial$time, "%d %b %Y %H:%M:%S")))
  d = 1
  while (d <= D) {
    # if (d %% 10 == 0) print(d)
    # index = which(time >= t_d-7*24*timeunit & time <= t_d)
    # sent = sender[index]
    # received = receiver[index, ]
    # outdegree = tabulate(sent, A)
    # indegree = colSums(received)
    # Y[d,,4] = outdegree
    # Y[d,,5] = indegree
    # Y[d,,6] = rep(as.numeric(wday(as.POSIXct(t_d, tz = getOption("tz"), origin = "1970-01-01")) %in% c(1, 7)), A)
    # Y[d,,7] = rep(as.numeric(pm(as.POSIXct(t_d, tz = getOption("tz"), origin = "1970-01-01"))), A)
    # for (a in 1:A) {
    #   for (r in c(1:A)[-a]) {
    #     X[d, a, r, 2] = outdegree[a]  
    #     X[d, a, r, 3] = indegree[r]	
    #     X[d, a, r, 4] = send(sent, received, a, r)
    #     X[d, a, r, 5] = send(sent, received, r, a)
    #     atoh = vapply(c(1:A)[-c(a,r)], function(h) {send(sent, received, a, h)}, c(1))
    #     htoa = vapply(c(1:A)[-c(a,r)], function(h) {send(sent, received, h, a)}, c(1))
    #     rtoh = vapply(c(1:A)[-c(a,r)], function(h) {send(sent, received, r, h)}, c(1))
    #     htor = vapply(c(1:A)[-c(a,r)], function(h) {send(sent, received, h, r)}, c(1))
    #     X[d, a, r, 6] = sum(atoh * htor) / 10
    #     X[d, a, r, 7] = sum(htoa * rtoh) / 10
    #     X[d, a, r, 8] = sum(htoa * htor) / 10
    #     X[d, a, r, 9] = sum(atoh * rtoh) / 10		
    #   }
    # }
    
    lambda[[d]] = lambda_cpp(X[d,,,], beta)
    u[[d]] = u_cpp_d(lambda[[d]], u[[d]])
    mu[d, ] = mu_cpp_d(Y[d,,], eta)
    tau = rlnorm(A, mu[d, ], sqrt(sigma2))
    for (n in 1:length(which(tau == min(tau)))) {
      a_d = which(tau == min(tau))
      r_d = u[[d]][a_d,]
      t_d = t_d + min(tau) * timeunit
      data[[d]] = list(a_d = a_d, r_d = r_d, t_d = t_d)
      d = d+1
      sender = c(sender, a_d)
      receiver = rbind(receiver, r_d)
      time = c(time, t_d)
    }
  }
  return(data)
}

Inference = function(data, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initial = initial,
		proposal.var, timeunit = 3600, lasttime) {
	D = dim(X)[1]
	A = dim(X)[2]
	P = dim(X)[4]
	Q = dim(Y)[3]
	
	if (length(initial) > 0) {
		u = initial$u
		beta = initial$beta
		eta = initial$eta
		sigma2 = initial$sigma2
	} else {
		u = lapply(1:D, function(d) matrix(0, A, A))
		beta = matrix(prior.beta$mean, nrow = 1)
		eta = matrix(prior.eta$mean, nrow = 1)
		sigma2 = prior.sigma2$b / (prior.sigma2$a-1)
	}
	#output matrix
	betamat = matrix(beta, nrow = outer-burn, ncol = P)
	etamat = matrix(eta, nrow = outer-burn, ncol = Q)
	sigma2mat = matrix(sigma2, nrow = outer-burn, ncol = 1)
	loglike = matrix(NA, nrow = outer-burn, ncol = 1)
	senders = vapply(data, function(d) { d[[1]] }, c(1))
	timestamps = vapply(data, function(d) { d[[3]] }, c(1))
	timeinc = c(timestamps[1]-lasttime, timestamps[-1]-timestamps[-length(timestamps)]) / timeunit
	timeinc[timeinc == 0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
	for (o in 1:outer) {
		if (o %% 100 == 0) print(o)
		lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta))
		u = u_cpp(lambda, u)
		for (d in 1:D) {
		  u[[d]][senders[d],] = data[[d]][[2]]
		}
		  prior.old1 = dmvnorm_arma(beta, prior.beta$mean, prior.beta$var)
    	post.old1 = Edgepartsum(lambda, u)
    	for (i1 in 1:inner[1]) {
			beta.new = rmvnorm_arma(1, beta, proposal.var[1]*diag(P))
     	prior.new1 = dmvnorm_arma(beta.new, prior.beta$mean, prior.beta$var)
			lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta.new))
			post.new1 = Edgepartsum(lambda, u)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
			if (log(runif(1, 0, 1)) < loglike.diff) {
        			beta = beta.new
        			prior.old1 = prior.new1
        			post.old1 = post.new1
	      	}
		}
		prior.old2 = dmvnorm_arma(eta, prior.eta$mean, prior.eta$var) 
    mu = mu_cpp(Y, eta)
    post.old2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
		for (i2 in 1:inner[2]) {
			eta.new = rmvnorm_arma(1, eta, proposal.var[2]*diag(Q))
      prior.new2 = dmvnorm_arma(eta.new, prior.eta$mean, prior.eta$var) 	
      mu = mu_cpp(Y, eta.new)
    	post.new2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
    	loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      		if (log(runif(1, 0, 1)) < loglike.diff) {
        			eta = eta.new
        			prior.old2 = prior.new2
        			post.old2 = post.new2
	      	}
		}
		prior.old3 = dinvgamma(sigma2, prior.sigma2$a, prior.sigma2$b) 
    post.old3 = post.old2
    mu = mu_cpp(Y, eta)

		for (i3 in 1:inner[3]) {
			sigma2.new = exp(rnorm(1, log(sigma2), proposal.var[3]))
      prior.new3 = dinvgamma(sigma2.new, prior.sigma2$a, prior.sigma2$b)
    	post.new3 = Timepartsum(mu, sqrt(sigma2.new), senders, timeinc)
    	loglike.diff = prior.new3+post.new3-prior.old3-post.old3
    			if (log(runif(1, 0, 1)) < loglike.diff) {
        			sigma2 = sigma2.new
        			prior.old3 = prior.new3
        			post.old3 = post.new3
	      	}
		}
		if (o > burn) {
			betamat[o-burn, ] = beta
			etamat[o-burn, ] = eta
			sigma2mat[o-burn, ] = sigma2	
			loglike[o-burn, ] = post.old1 + post.old3
		}		
	}
	return(list(u = u, beta = betamat, eta = etamat, sigma2 = sigma2mat, loglike = loglike))
}

# missing is a list object
PPE = function(data, missing, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initial = initial,
		proposal.var, timeunit = 3600, lasttime) {
	D = dim(X)[1]
	A = dim(X)[2]
	P = dim(X)[4]
	Q = dim(Y)[3]
	
	if (length(initial) > 0) {
		u = initial$u
		beta = initial$beta
		eta = initial$eta
		sigma2 = initial$sigma2
	} else {
		u = lapply(1:D, function(d) matrix(0, A, A))
		beta = matrix(prior.beta$mean, nrow = 1)
		eta = matrix(prior.eta$mean, nrow = 1)
		sigma2 = prior.sigma2$b / (prior.sigma2$a-1)
	}
	#output matrix
	betamat = matrix(beta, nrow = outer-burn, ncol = P)
	etamat = matrix(eta, nrow = outer-burn, ncol = Q)
	sigma2mat = matrix(sigma2, nrow = outer-burn, ncol = 1)
	loglike = matrix(NA, nrow = outer-burn, ncol = 1)
	senders = vapply(data, function(d) { d[[1]] }, c(1))
	timestamps = vapply(data, function(d) { d[[3]] }, c(1))
	timeinc = c(timestamps[1]-lasttime, timestamps[-1]-timestamps[-length(timestamps)]) / timeunit
	timeinc[timeinc == 0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
		
	senderpredict = matrix(NA, nrow = sum(missing[[1]]), ncol = outer)
    receiverpredict = lapply(1:sum(missing[[2]]), function(d) {c()})
    timepredict = matrix(NA, nrow = sum(missing[[3]]), ncol = outer)
    sendermissing = which(missing[[1]]==1)
    receivermissing = which(missing[[2]]==1)
    timemissing = which(missing[[3]]==1)
	
	for (o in 1:outer) {
		
	#imputation
    iter1 = 1
    iter2 = 1
    iter3 = 1
    for (d in sendermissing) {
        probi = Timepartindiv(mu[d,], sigma_tau, timeinc[d])
        senders[d] = lmultinom(probi)
        senderpredict[iter1, o] = senders[d]
        iter1 = iter1+1
    }
    for (d in timemissing) {
        timeinc[d] = rlnorm(1, mu[d, senders[d]], sigma_tau)
        timepredict[iter2, o] = timeinc[d]
        iter2 = iter2+1
    }
    
    timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
    timestamps[-1] = timestamps[1]+cumsum(timeinc[-1])*timeunit
    for (d in edge.trim) {
        history.t = History(edge, timestamps, p.d, node, d, timeunit)
        X[[d]] = Netstats_cpp(history.t, node, netstat)
    }
    for (d in receivermissing) {
        vu = MultiplyXB(X[[d]], b.old)
        lambda = lambda_cpp(p.d[d,], vu)
        i = senders[d]
        for (j in sample(node[-i], A-1)) {
            probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
            u[[d]][i, j] = lmultinom(probij)-1
        }
        receiverpredict[[iter3]] = rbind(receiverpredict[[iter3]], u[[d]][i, ])
        iter3 = iter3+1
    }  

	#run inference
		if (o %% 100 == 0) print(o)
		lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta))
		u = u_cpp(lambda, u)
		for (d in 1:D) {
		  u[[d]][senders[d],] = data[[d]][[2]]
		}
		prior.old1 = dmvnorm_arma(beta, prior.beta$mean, prior.beta$var)
    	post.old1 = Edgepartsum(lambda, u)
    	for (i1 in 1:inner[1]) {
			beta.new = rmvnorm_arma(1, beta, proposal.var[1]*diag(P))
     		prior.new1 = dmvnorm_arma(beta.new, prior.beta$mean, prior.beta$var)
			lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta.new))
			post.new1 = Edgepartsum(lambda, u)
      		loglike.diff = prior.new1+post.new1-prior.old1-post.old1
			if (log(runif(1, 0, 1)) < loglike.diff) {
        			beta = beta.new
        			prior.old1 = prior.new1
        			post.old1 = post.new1
	      	}
		}
		prior.old2 = dmvnorm_arma(eta, prior.eta$mean, prior.eta$var) 
    	mu = mu_cpp(Y, eta)
   		post.old2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
		for (i2 in 1:inner[2]) {
			eta.new = rmvnorm_arma(1, eta, proposal.var[2]*diag(Q))
      		prior.new2 = dmvnorm_arma(eta.new, prior.eta$mean, prior.eta$var) 	
      		mu = mu_cpp(Y, eta.new)
    		post.new2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
    		loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      		if (log(runif(1, 0, 1)) < loglike.diff) {
        			eta = eta.new
        			prior.old2 = prior.new2
        			post.old2 = post.new2
	      	}
		}
		prior.old3 = dinvgamma(sigma2, prior.sigma2$a, prior.sigma2$b) 
    	post.old3 = post.old2
   	 	mu = mu_cpp(Y, eta)

		for (i3 in 1:inner[3]) {
			sigma2.new = exp(rnorm(1, log(sigma2), proposal.var[3]))
     	 	prior.new3 = dinvgamma(sigma2.new, prior.sigma2$a, prior.sigma2$b)
    		post.new3 = Timepartsum(mu, sqrt(sigma2.new), senders, timeinc)
    		loglike.diff = prior.new3+post.new3-prior.old3-post.old3
    			if (log(runif(1, 0, 1)) < loglike.diff) {
        			sigma2 = sigma2.new
        			prior.old3 = prior.new3
        			post.old3 = post.new3
	      	}
		}
		if (o > burn) {
			betamat[o-burn, ] = beta
			etamat[o-burn, ] = eta
			sigma2mat[o-burn, ] = sigma2	
			loglike[o-burn, ] = post.old1 + post.old3
		}		
	}
	return(list(u = u, beta = betamat, eta = etamat, sigma2 = sigma2mat, loglike = loglike))
}


Inference_ties = function(data, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initial = initial,
                     proposal.var, timeunit = 3600, lasttime) {
  D = dim(X)[1]
  A = dim(X)[2]
  P = dim(X)[4]
  Q = dim(Y)[3]
  
  if (length(initial) > 0) {
    u = initial$u
    beta = initial$beta
    eta = initial$eta
    sigma2 = initial$sigma2
  } else {
    u = lapply(1:D, function(d) matrix(0, A, A))
    beta = matrix(prior.beta$mean, nrow = 1)
    eta = matrix(prior.eta$mean, nrow = 1)
    sigma2 = prior.sigma2$b / (prior.sigma2$a-1)
  }
  #output matrix
  betamat = matrix(beta, nrow = outer-burn, ncol = P)
  etamat = matrix(eta, nrow = outer-burn, ncol = Q)
  sigma2mat = matrix(sigma2, nrow = outer-burn, ncol = 1)
  loglike = matrix(NA, nrow = outer-burn, ncol = 1)
  senders = vapply(data, function(d) { d[[1]] }, c(1))
  timestamps = vapply(data, function(d) { d[[3]] }, c(1))
  uniqtime = unique(timestamps)
  timeinc = c(timestamps[1]-lasttime, timestamps[-1]-timestamps[-length(timestamps)]) / timeunit
  for (o in 1:outer) {
    if (o %% 100 == 0) print(o)
    lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta))
    u = u_cpp(lambda, u)
    for (d in 1:D) {
      u[[d]][senders[d],] = data[[d]][[2]]
    }
    prior.old1 = dmvnorm_arma(beta, prior.beta$mean, prior.beta$var)
    post.old1 = Edgepartsum(lambda, u)
    for (i1 in 1:inner[1]) {
      beta.new = rmvnorm_arma(1, beta, proposal.var[1]*diag(P))
      prior.new1 = dmvnorm_arma(beta.new, prior.beta$mean, prior.beta$var)
      lambda = lapply(1:D, function(d) lambda_cpp(X[d,,,], beta.new))
      post.new1 = Edgepartsum(lambda, u)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        beta = beta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
      }
    }
    prior.old2 = dmvnorm_arma(eta, prior.eta$mean, prior.eta$var) 
    mu = mu_cpp(Y, eta)
    post.old2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
    for (i2 in 1:inner[2]) {
      eta.new = rmvnorm_arma(1, eta, proposal.var[2]*diag(Q))
      prior.new2 = dmvnorm_arma(eta.new, prior.eta$mean, prior.eta$var) 	
      mu = mu_cpp(Y, eta.new)
      post.new2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      }
    }
    prior.old3 = dinvgamma(sigma2, prior.sigma2$a, prior.sigma2$b) 
    post.old3 = post.old2
    mu = mu_cpp(Y, eta)
    
    for (i3 in 1:inner[3]) {
      sigma2.new = exp(rnorm(1, log(sigma2), proposal.var[3]))
      prior.new3 = dinvgamma(sigma2.new, prior.sigma2$a, prior.sigma2$b)
      post.new3 = Timepartsum(mu, sqrt(sigma2.new), senders, timeinc)
      loglike.diff = prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma2 = sigma2.new
        prior.old3 = prior.new3
        post.old3 = post.new3
      }
    }
    if (o > burn) {
      betamat[o-burn, ] = beta
      etamat[o-burn, ] = eta
      sigma2mat[o-burn, ] = sigma2	
      loglike[o-burn, ] = post.old1 + post.old3
    }		
  }
  return(list(u = u, beta = betamat, eta = etamat, sigma2 = sigma2mat, loglike = loglike))
}

