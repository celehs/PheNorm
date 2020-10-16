
#' Fit the phenotyping algorithm with EHR features. The function requires a surrogate (ICD) and
#' the health utilization as its input and can leverage other EHR features (optional) to assist
#' risk prediction.
#' @param nm.logS.ori name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
#' @param dat all data columns need to be log-transformed and need column names
#' @param nm.X additional features other than the main ICD and NLP
#' @param corrupt.rate rate for random corruption denoising, between 0 and 1
#' @param train.size size of training sample, default value 10 * nrow(dat)
#' @return beta coefficient
#' @export
PheNorm_noUTL_beta <- function(nm.logS.ori, dat, nm.X = NULL, corrupt.rate = 0.3, train.size = 10 * nrow(dat)){
  dat <- as.matrix(dat)
  S.ori <- dat[, nm.logS.ori, drop = FALSE]
  S.norm <- S.ori
  if (!is.null(nm.X)) {
    X <- as.matrix(dat[, nm.X])
    X.norm <- X
    SX.norm <- cbind(S.norm, X.norm)
    id <- sample(1:nrow(dat), train.size, replace=TRUE)
    SX.norm.corrupt <- apply(SX.norm[id, ], 2,
                             function(x) {ifelse(rbinom(length(id), 1, corrupt.rate), mean(x), x)}
                             )
    b.all <- apply(S.norm, 2, function(ss) {lm(ss[id] ~ SX.norm.corrupt - 1)$coef})
    S.norm <- as.matrix(SX.norm) %*% b.all
  }
  return(b.all)
}
