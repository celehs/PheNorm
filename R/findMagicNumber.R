findMagicNumber = function(surrogate, log_note_count, n.boot = 10) {
  a.values <- rep(1, n.boot)
  err.values <- rep(Inf, n.boot)
  for (k in 1:n.boot) {
    idx <- sample(1:length(log_note_count), replace = n.boot > 1) # bootstrap to make the process more robust
    coefs <- seq(0, 1.2, 0.05)
    err <- rep(0, length(coefs))
    err.prev <- Inf
    for (j in 1:length(coefs)) {
      score <- surrogate[idx] - coefs[j] * log_note_count[idx]
      fit <- normalmixEM2comp2(score, lambda = 0.5, mu = quantile(score, probs=c(1/3, 2/3)), sigsqrd = 1)
      # the following approximates \int|Fn(x)-Fmix(x)|dx
      Fn <- ecdf(score) # empirical cdf
      Fmix = function(x) {fit$lambda[1] * pnorm(x, fit$mu[1], fit$sigma) + fit$lambda[2]*pnorm(x, fit$mu[2], fit$sigma)}
      id.small <- which.min(fit$mu)
      id.large <- which.max(fit$mu)
      fit.range <- seq(fit$mu[id.small] - 5*fit$sigma, fit$mu[id.large] + 5*fit$sigma, length.out = 1000)
      x.step <- (10*fit$sigma + fit$mu[id.large] - fit$mu[id.small])/length(fit.range)
      fit.err = function(x){abs(Fn(x)-Fmix(x))}
      err[j] <- sum(fit.err(fit.range))*x.step
      if (err[j] > err.prev)
        break
      else
        err.prev <- err[j]
    }
    a.values[k] <- coefs[j-1]
    err.values[k] <- err[j-1]
  }
  return(list("coef" = mean(a.values), "error" = mean(err.values)))
}
