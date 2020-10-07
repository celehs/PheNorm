#' Fit the phenotyping algorithm with EHR features. The function requires a surrogate (ICD) and
#' the health utilization as its input and can leverage other EHR features (optional) to assist
#' risk prediction.
#' @param nm.logS.ori name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
#' @param nm.utl name of healthcare utlization (e.g. note count, encounter_num etc)
#' @param dat all data columns need to be log-transformed and need column names
#' @param nm.X additional features other than the main ICD and NLP
#' @param corrupt.rate
#' @param train.size
#' @return
#' @export
PheNorm.Prob = function(nm.logS.ori,nm.utl,dat, nm.X=NULL,corrupt.rate=0.3,train.size=10000){
  #browser()
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]; utl = dat[,nm.utl]
  a.hat = apply(S.ori, 2, function(S){findMagicNumber(S,utl)$coef})
  S.norm = S.ori - VTM(a.hat,nrow(dat))*utl
  if(!is.null(nm.X)){
    X = as.matrix(dat[,nm.X])
    SX.norm = cbind(S.norm,X,utl)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    b.all[is.na(b.all)] = 0
    S.norm = as.matrix(SX.norm)%*%b.all
    b.all = b.all[-dim(b.all)[1],]
  }
  else{
    b.all = NULL
  }
  if(length(nm.logS.ori)>1){
    postprob = apply(S.norm,2,function(x){fit = normalmixEM2comp2(x, lambda=0.5, mu=quantile(x,probs=c(1/3,2/3)), sigsqrd=1);fit$posterior[,2]})
    list("probs"=rowMeans(postprob,na.rm = T), "betas"=b.all)

  }else{
    fit = normalmixEM2comp2(unlist(S.norm), lambda=0.5, mu=quantile(S.norm,probs=c(1/3,2/3)), sigsqrd=1)
    list("probs"=fit$posterior[,2], "betas"=b.all)
  }
}

PheNorm.Prob_noUTL = function(nm.logS.ori,dat, nm.X=NULL,corrupt.rate=0.3,train.size=10000){
  #browser()
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]
  S.norm = S.ori
  if(!is.null(nm.X)){
    X = as.matrix(dat[,nm.X])
    SX.norm = cbind(S.norm,X)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    b.all[is.na(b.all)] = 0
    S.norm = as.matrix(SX.norm)%*%b.all
    b.all = b.all[-dim(b.all)[1],]
  }
  else{
    b.all = NULL
  }
  if(length(nm.logS.ori)>1){
    postprob = apply(S.norm,2,function(x){fit = normalmixEM2comp2(x, lambda=0.5, mu=quantile(x,probs=c(1/3,2/3)), sigsqrd=1);fit$posterior[,2]})
    list("probs"=rowMeans(postprob,na.rm = T), "betas"=b.all)

  }else{
    fit = normalmixEM2comp2(unlist(S.norm), lambda=0.5, mu=quantile(S.norm,probs=c(1/3,2/3)), sigsqrd=1)
    list("probs"=fit$posterior[,2], "betas"=b.all)
  }
}

PheNorm = function(nm.logS.ori,nm.utl,dat, nm.X=NULL,corrupt.rate=0.3,train.size=100000){
  ## dat: all data columns need to be log-transformed and need column names; ##
  ## nm.logS.ori is the name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
  ## nm.utl: is the name of healthcare utlization (e.g. note count, encounter_num etc)
  ## nm.X: additional features other than the main ICD and NLP
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]; utl = dat[,nm.utl]
  a.hat = apply(as.matrix(S.ori), 2, function(S){findMagicNumber(S,utl)$coef})
  S.norm = S.ori - VTM(a.hat,nrow(dat))*utl
  if(!is.null(nm.X)){
    #ZH:apply instead of sapply
    X = as.matrix(dat[,nm.X]); a.X = apply(X,2,function(xx){findMagicNumber(xx,utl)$coef})
    X.norm = X - utl %*% t(a.X); SX.norm = cbind(S.norm,X.norm)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    ## SX.norm insead of X.norm
    S.norm = as.matrix(SX.norm)%*%b.all
  }
  if(length(nm.logS.ori)>1){
    postprob = apply(S.norm,2,function(x){fit = normalmixEM2comp2(x, lambda=0.5, mu=quantile(x,probs=c(1/3,2/3)), sigsqrd=1);fit$posterior[,2]})
    keep = apply(postprob,1,function(x){if(sum(x>0.5)>=2) x[which(x<0.5)]=NA else x[which(x>0.5)]=NA; x})
    keep = as.matrix(1*(keep>=0)); if(nrow(keep)!=nrow(dat)){keep=t(keep)}
    rowMeans(S.norm*keep,na.rm = T)
  }else{
    unlist(S.norm)
  }
}

PheNorm_noUTL = function(nm.logS.ori,dat, nm.X=NULL,corrupt.rate=0.3,train.size=100000){
  ## dat: all data columns need to be log-transformed and need column names; ##
  ## nm.logS.ori is the name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
  ## nm.utl: is the name of healthcare utlization (e.g. note count, encounter_num etc)
  ## nm.X: additional features other than the main ICD and NLP
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]
  S.norm = S.ori
  if(!is.null(nm.X)){
    #ZH:apply instead of sapply
    X = as.matrix(dat[,nm.X])
    X.norm = X; SX.norm = cbind(S.norm,X.norm)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    ## SX.norm insead of X.norm
    S.norm = as.matrix(SX.norm)%*%b.all
  }
  if(length(nm.logS.ori)>1){
    postprob = apply(S.norm,2,function(x){fit = normalmixEM2comp2(x, lambda=0.5, mu=quantile(x,probs=c(1/3,2/3)), sigsqrd=1);fit$posterior[,2]})
    keep = apply(postprob,1,function(x){if(sum(x>0.5)>=2) x[which(x<0.5)]=NA else x[which(x>0.5)]=NA; x})
    keep = as.matrix(1*(keep>=0)); if(nrow(keep)!=nrow(dat)){keep=t(keep)}
    rowMeans(S.norm*keep,na.rm = T)
  }else{
    unlist(S.norm)
  }
}

PheNorm_noUTL_beta = function(nm.logS.ori,dat, nm.X=NULL,corrupt.rate=0.3,train.size=100000){
  ## dat: all data columns need to be log-transformed and need column names; ##
  ## nm.logS.ori is the name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
  ## nm.utl: is the name of healthcare utlization (e.g. note count, encounter_num etc)
  ## nm.X: additional features other than the main ICD and NLP
  dat = as.matrix(dat)
  S.ori = dat[,nm.logS.ori,drop=F]
  S.norm = S.ori
  if(!is.null(nm.X)){
    #ZH:apply instead of sapply
    X = as.matrix(dat[,nm.X])
    X.norm = X; SX.norm = cbind(S.norm,X.norm)
    id = sample(1:nrow(dat), train.size, replace=T)
    SX.norm.corrupt = apply(SX.norm[id,],2,function(x){ifelse(rbinom(length(id),1,corrupt.rate),mean(x),x)})
    b.all = apply(S.norm, 2, function(ss){lm(ss[id]~SX.norm.corrupt-1)$coef})
    ## SX.norm insead of X.norm
    S.norm = as.matrix(SX.norm)%*%b.all
  }
  return(b.all)
}


normalmixEM2comp2 <- function (x, lambda, mu, sigsqrd, eps = 1e-08, maxit = 1000, verb = FALSE) {
  arbvar <- (length(sigsqrd) == 2)
  mu1 <- mu[1]
  mu2 <- mu[2]
  sigsqrd1 <- sigsqrd[1]
  sigsqrd2 <- sigsqrd[arbvar + 1]
  mx <- mean(x)
  const <- length(x) * 0.918938533204673
  dl <- 1 + eps
  iter <- 0
  ll <- rep(0, maxit + 1)
  a1 <- (x - mu1)^2
  b1 <- (lambda/sqrt(sigsqrd1)) * exp(-a1/2/sigsqrd1)
  a2 <- (x - mu2)^2
  b2 <- ((1 - lambda)/sqrt(sigsqrd2)) * exp(-a2/2/sigsqrd2)
  l <- sum(log(b1 + b2))
  while (dl > eps && iter < maxit) {
    iter <- iter + 1
    ll[iter] <- l
    postprobs <- b1/(b1 + b2)
    lambda <- mean(postprobs)
    mu1 <- mean(postprobs * x)/lambda
    mu2 <- (mx - lambda * mu1)/(1 - lambda)
    if (arbvar) {
      sigsqrd1 <- mean(postprobs * a1)/lambda
      sigsqrd2 <- mean((1 - postprobs) * a2)/(1 - lambda)
    }
    else {
      sigsqrd1 <- sigsqrd2 <- mean(postprobs * a1 + (1 -
                                                       postprobs) * a2)
    }
    a1 <- (x - mu1)^2
    b1 <- (lambda/sqrt(sigsqrd1)) * exp(-a1/2/sigsqrd1)
    a2 <- (x - mu2)^2
    b2 <- ((1 - lambda)/sqrt(sigsqrd2)) * exp(-a2/2/sigsqrd2)
    oldl <- l
    l <- sum(log(b1 + b2))
    dl <- l - oldl
    if (verb) {
      cat("iteration =", iter, " log-lik diff =", dl, " log-lik =",
          l - const, "\n")
    }
  }
  ##cat("number of iterations=", iter, "\n")
  iter <- iter + 1
  ll[iter] <- l
  postprobs <- cbind(postprobs, 1 - postprobs)
  colnames(postprobs) <- c(paste("comp", ".", 1:2, sep = ""))
  out <- list(x = x, lambda = c(lambda, 1 - lambda), mu = c(mu1,
                                                            mu2), sigma = sqrt(c(sigsqrd1, sigsqrd2)[1:(1 + arbvar)]),
              loglik = l - const, posterior = postprobs, all.loglik = ll[1:iter] -
                const, restarts = 0, ft = "normalmixEM")
  class(out) <- "mixEM"
  out
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


findMagicNumber = function(surrogate, log_note_count, n.boot=10) {
  a.values = rep(1,n.boot)
  err.values = rep(Inf,n.boot)
  for(k in 1:n.boot) {
    idx = sample(1:length(log_note_count), replace=n.boot>1) # bootstrap to make the process more robust
    coefs = seq(0,1.2,0.05)
    err = rep(0, length(coefs))
    err.prev = Inf
    for (j in 1:length(coefs)) {
      score = surrogate[idx] - coefs[j] * log_note_count[idx]
      fit = normalmixEM2comp2(score, lambda=0.5, mu=quantile(score,probs=c(1/3,2/3)), sigsqrd=1)
      # the following approximates \int|Fn(x)-Fmix(x)|dx
      Fn = ecdf(score) # empirical cdf
      Fmix = function(x) {fit$lambda[1]*pnorm(x,fit$mu[1],fit$sigma) + fit$lambda[2]*pnorm(x,fit$mu[2],fit$sigma)}
      id.small = which.min(fit$mu); id.large = which.max(fit$mu)
      fit.range = seq(fit$mu[id.small]-5*fit$sigma, fit$mu[id.large]+5*fit$sigma, length.out=1000)
      x.step = (10*fit$sigma + fit$mu[id.large] - fit$mu[id.small])/length(fit.range)
      fit.err = function(x){abs(Fn(x)-Fmix(x))}
      err[j] = sum(fit.err(fit.range))*x.step
      if (err[j] > err.prev)
        break
      else
        err.prev = err[j]
    }
    a.values[k] = coefs[j-1]
    err.values[k] = err[j-1]
  }
  return(list("coef"=mean(a.values), "error"=mean(err.values)))
}

