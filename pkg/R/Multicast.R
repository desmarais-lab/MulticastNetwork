#' @useDynLib MulticastNetwork
#' @import stats
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom combinat permn
#' @importFrom lubridate wday hour
#' @importFrom MCMCpack rinvgamma dinvgamma


#' @title gibbs.measure.support
#' @description List out the support of Gibbs measure
#'
#' @param n length of the vector to be sampled
#'
#' @return a 2^n x n binary matrix representing the support of the binary Gibbs measure in n elements
#'
#' @export
gibbs.measure.support = function(n) {
	gibbs.support = rbind(rep(1, n))
	for(i in 1:(n-1)){
		gibbs.mat.i = do.call('rbind',permn(c(rep(1, i), rep(0, n-i))))
		gibbs.support = rbind(gibbs.support, gibbs.mat.i)
	}
	out = as.matrix(unique(gibbs.support))
	return(out)
}

#' @title r.gibbs.measure
#' @description Sample a binary vector from Gibbs measure
#'
#' @param lambda.i lambda of sender i
#' @param support support of Gibbs measure
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
r.gibbs.measure <- function(lambda.i, support) {
	logitNumerator = vapply(1:nrow(support), function(s) {
		sum(lambda.i*support[s,])
		}, c(1))		
	samp = multinom_vec(exp(logitNumerator))
	return(support[samp,])	
}


#' @title dinvgamma
#' @description log of inverse-gamma pdf
#'
#' @param x Scalar location to evaluate density
#' @param shape Scalar shape parameter
#' @param scale Scalar shape parameter
#'
#' @return evaluate the density at x
#'
#' @export
dinvgamma = function(x, shape, scale) {
    alpha <- shape
    beta <- scale
    log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 
        1) * log(x) - (beta/x)
    return(log.density)
}


#' @title GiR_PP_Plots
#' @description Generate PP-plots (probability-probability plot) for Gettig it Right test
#'
#' @param Forward_stats statistics obtained from forward sampling
#' @param Backward_stats statistics obtained from backward sampling
#'
#' @return PP-plots for different GiR statistics of interest
#'
#' @export
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

#' @title Generate
#' @description Generate a collection of events according to the generative process
#'
#' @param D number of events to be generated
#' @param A vector of node id's (ID starting from 1)
#' @param beta P-length vector of coefficients for recipients
#' @param eta Q-length vector of coefficients for timestamps
#' @param sigma2 variance parameter for the timestamps
#' @param X an array of dimension D x A x A x P for covariates used for Gibbs measure
#' @param Y an array of dimension D x A x Q for covariates used for timestamps GLM
#' @param support support of latent recipients from Gibbs measure
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param timedist lognormal or exponential (will include others)
#'
#' @return generated data including (sender, recipients, timestamp)
#'
#' @export
Generate = function(D, A, beta, eta, sigma2, X, Y, support, timeunit = 3600, timedist = "lognormal") {
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
		if (timedist == "lognormal") {
            tau = rlnorm(A, mu[d, ], sqrt(sigma2))
        } else {
            tau = rexp(A, 1/exp(mu[d, ]))
        }
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

#' @title PPC
#' @description Generate a collection of events according to the posterior predictive checks
#'
#' @param D number of events to be generated
#' @param beta P-length vector of coefficients for recipients
#' @param eta Q-length vector of coefficients for timestamps
#' @param sigma2 variance parameter for the timestamps
#' @param X an array of dimension D x A x A x P for covariates used for Gibbs measure
#' @param Y an array of dimension D x A x Q for covariates used for timestamps GLM
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param u D-length list of latent receiver vectors
#' @param lasttime last timestamp of the event used as initial history (in unix.time format)
#' @param timedist lognormal or exponential (will include others)
#'
#' @return generated data including (sender, recipients, timestamp)
#'
#' @export
PPC = function(D, beta, eta, sigma2, X, Y, timeunit = 3600, u, lasttime, timedist = "lognormal") {
  P = length(beta)
  Q = length(eta)
  A = dim(X)[2]
  u = u
  data = list()
  lambda = list()
  mu = matrix(NA, D, A)
  t_d = lasttime
  d = 1
  while (d <= D) {
    lambda[[d]] = lambda_cpp(X[d,,,], beta)
    u[[d]] = u_cpp_d(lambda[[d]], u[[d]])
    mu[d, ] = mu_cpp_d(Y[d,,], eta)
    if (timedist == "lognormal") {
        tau = rlnorm(A, mu[d, ], sqrt(sigma2))
    } else {
        tau = rexp(A, 1/exp(mu[d, ]))
    }
    for (n in 1:length(which(tau == min(tau)))) {
      a_d = which(tau == min(tau))
      r_d = u[[d]][a_d,]
      t_d = t_d + min(tau) * timeunit
      data[[d]] = list(a_d = a_d, r_d = r_d, t_d = t_d)
      d = d+1
    }
  }
  return(data)
}

#' @title Inference
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm to infer the parameters
#'
#' @param data list of tie data with 3 elements (1: sender, 2: recipient, 3: timestamp in unix.time format)
#' @param X an array of dimension D x A x A x P for covariates used for Gibbs measure
#' @param Y an array of dimension D x A x Q for covariates used for timestamps GLM
#' @param outer size of outer iterations
#' @param inner size of inner iteration for Metropolis-Hastings updates
#' @param burn size of burn-in
#' @param prior.beta prior mean and covariance of beta in multivariate Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.sigma2 prior shape and scale parameter of sigma2 in inverse-Gamma distribution
#' @param initialval initial value of the parameters (if speficied)
#' @param proposal.var proposal variance for beta, eta, and sigma2
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param lasttime last timestamp of the event used as initial history (in unix.time format)
#' @param timedist lognormal or exponential (will include others)
#'
#' @return generated data including (sender, recipients, timestamp)
#'
#' @export
Inference = function(data, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initialval = NULL, proposal.var, timeunit = 3600, lasttime, timedist = "lognormal") {
	D = dim(X)[1]
	A = dim(X)[2]
	P = dim(X)[4]
	Q = dim(Y)[3]
	
	if (length(initialval) > 0) {
		u = initialval$u
		beta = matrix(initialval$beta, nrow = 1)
		eta = matrix(initialval$eta, nrow = 1)
		sigma2 = initialval$sigma2
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
        if (timedist == "lognormal") {
            post.old2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
        } else {
            post.old2 = Timepartsum2(mu, senders, timeinc)
        }
		for (i2 in 1:inner[2]) {
			eta.new = rmvnorm_arma(1, eta, proposal.var[2]*diag(Q))
     	 	prior.new2 = dmvnorm_arma(eta.new, prior.eta$mean, prior.eta$var) 	
      		mu = mu_cpp(Y, eta.new)
            if (timedist == "lognormal") {
                post.new2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
            } else {
                post.new2 = Timepartsum2(mu, senders, timeinc)
            }
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
        
        if (timedist == "lognormal") {
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


#' @title PPE
#' @description Posterior predictive experiments for
#'
#' @param data list of tie data with 3 elements (1: sender, 2: recipient, 3: timestamp in unix.time format)
#' @param missing list of indicators denoting missingness of sender, receiver, and timestamps
#' @param X an array of dimension D x A x A x P for covariates used for Gibbs measure
#' @param Y an array of dimension D x A x Q for covariates used for timestamps GLM
#' @param outer size of outer iterations
#' @param inner size of inner iteration for Metropolis-Hastings updates
#' @param burn size of burn-in
#' @param prior.beta prior mean and covariance of beta in multivariate Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.sigma2 prior shape and scale parameter of sigma2 in inverse-Gamma distribution
#' @param initial initial value of the parameters (if speficied)
#' @param proposal.var proposal variance for beta, eta, and sigma2
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param lasttime last timestamp of the event used as initial history (in unix.time format)
#' @param MHprop.var proposal variance for time predictions
#' @param timedist lognormal or exponential (will include others)
#'
#' @return generated data including (sender, recipients, timestamp)
#'
#' @export
PPE = function(data, missing, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initial = NULL,
		proposal.var, timeunit = 3600, lasttime, MHprop.var, timedist = "lognormal") {
	D = dim(X)[1]
	A = dim(X)[2]
	P = dim(X)[4]
	Q = dim(Y)[3]
	
	if (length(initial) > 0) {
		u = initial$u
		beta = matrix(initial$beta, nrow = 1)
		eta = matrix(initial$eta, nrow = 1)
		sigma2 = initial$sigma2
	} else {
		u = lapply(1:D, function(d) matrix(0, A, A))
		beta = matrix(prior.beta$mean, nrow = 1)
		eta = matrix(prior.eta$mean, nrow = 1)
		sigma2 = prior.sigma2$b / (prior.sigma2$a-1)
	}
	mu = mu_cpp(Y, eta)
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
    receiverpredict = matrix(NA, nrow = sum(missing[[2]]), ncol = outer)
    timepredict = matrix(NA, nrow = sum(missing[[3]]), ncol = outer)
    sendermissing = which(missing[[1]]==1)
    receivermissing = which(rowSums(missing[[2]]) > 0)
    timemissing = which(missing[[3]]==1)
	timeinc[timemissing] = vapply(timemissing, function(d) rlnorm(1, mu[d, senders[d]], sqrt(sigma2)), c(1))
	for (o in 1:outer) {
		
	#imputation
    iter1 = 1
    iter2 = 1
    iter3 = 1
    for (d in sendermissing) {
        probi = Timepartindiv(mu[d,], sqrt(sigma2), timeinc[d])
        senders[d] = multinom_vec(exp(probi))
        senderpredict[iter1, o] = senders[d]
        iter1 = iter1+1
    }
    for (d in receivermissing) {
    	missingr = which(missing[[2]][d,]==1)
    	data[[d]][[2]][missingr] = u[[d]][senders[d], missingr]
        for (it in 1:length(missingr)) {
        	receiverpredict[iter3, o] = data[[d]][[2]][missingr[it]]
        	iter2 = iter2+1
        }
    }  
    for (d in timemissing) {
    	tau_new = rlnorm(1, mu[d, senders[d]], MHprop.var)
    	prior.new0 = dlnorm(tau_new, mu[d, senders[d]], sqrt(sigma2), TRUE)
    	prior.old0 = dlnorm(timeinc[d], mu[d, senders[d]], sqrt(sigma2), TRUE)
    	post.new0 = Timepartindiv(mu[d,], sqrt(sigma2), tau_new)[senders[d]]
    	post.old0 = Timepartindiv(mu[d,], sqrt(sigma2), timeinc[d])[senders[d]]
		loglike.diff = prior.new0+post.new0-prior.old0-post.old0
    	if (log(runif(1, 0, 1)) < loglike.diff) {
        	timeinc[d] = tau_new
	    }
        timepredict[iter3, o] = timeinc[d]
        iter3 = iter3+1
    }    
    timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))

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
   		if (timedist == "lognormal") {
            post.old2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
        } else {
            post.old2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
        }
		for (i2 in 1:inner[2]) {
			eta.new = rmvnorm_arma(1, eta, proposal.var[2]*diag(Q))
      		prior.new2 = dmvnorm_arma(eta.new, prior.eta$mean, prior.eta$var) 	
      		mu = mu_cpp(Y, eta.new)
    		if (timedist == "lognormal") {
                post.new2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
            } else {
                post.new2 = Timepartsum(mu, sqrt(sigma2), senders, timeinc)
            }
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

        if (timedist == "lognormal") {
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
        }
		if (o > burn) {
			betamat[o-burn, ] = beta
			etamat[o-burn, ] = eta
			sigma2mat[o-burn, ] = sigma2	
			loglike[o-burn, ] = post.old1 + post.old3
		}		
	}
	return(list(senderpredict = senderpredict, receiverpredict = receiverpredict, timepredict = timepredict))
}
