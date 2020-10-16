#' Fit the phenotyping algorithm PheNorm using EHR features
#'
#' @description
#' The function requires as input:
#' * a surrogate, such as the ICD code
#' * the healthcare utilization
#' It can leverage other EHR features (optional) to assist risk prediction.
#'
#' @param nm.logS.ori name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
#' @param nm.utl name of healthcare utilization (e.g. note count, encounter_num etc)
#' @param dat all data columns need to be log-transformed and need column names
#' @param nm.X additional features other than the main ICD and NLP
#' @param corrupt.rate rate for random corruption denoising, between 0 and 1
#' @param train.size size of training sample, default value 10 * nrow(dat)
#' @return S.norm
#' @export
PheNorm <- function(nm.logS.ori, nm.utl, dat, nm.X = NULL, corrupt.rate = 0.3, train.size = 10 * nrow(dat)){
  dat <- as.matrix(dat)
  S.ori <- dat[, nm.logS.ori, drop = FALSE]
  utl <- dat[, nm.utl]
  a.hat <- apply(as.matrix(S.ori), 2, function(S) {findMagicNumber(S, utl)$coef})
  S.norm <- S.ori - VTM(a.hat, nrow(dat)) * utl
  if (!is.null(nm.X)) {
    X <- as.matrix(dat[, nm.X])
    a.X <- apply(X, 2, function(xx) {findMagicNumber(xx, utl)$coef})
    X.norm <- X - utl %*% t(a.X)
    SX.norm <- cbind(S.norm, X.norm)
    id <- sample(1:nrow(dat), train.size, replace = TRUE)
    SX.norm.corrupt <- apply(SX.norm[id, ], 2,
                             function(x) {
                               ifelse(rbinom(length(id), 1, corrupt.rate), mean(x), x)
                               }
                             )
    b.all <- apply(S.norm, 2, function(ss) {lm(ss[id] ~ SX.norm.corrupt - 1)$coef})
    S.norm <- as.matrix(SX.norm) %*% b.all
  }
  if (length(nm.logS.ori) > 1) {
    postprob <- apply(S.norm, 2,
                      function(x) {
                        fit = normalmixEM2comp2(x, lambda = 0.5,
                                                mu = quantile(x, probs = c(1/3, 2/3)), sigsqrd = 1
                                                )
                        fit$posterior[, 2]
                        }
                      )
    keep <- apply(postprob, 1,
                  function(x) {
                    if (sum(x > 0.5) >= 2) x[which(x < 0.5)] = NA else x[which(x > 0.5)] = NA
                    x
                    }
                  )
    keep <- as.matrix(1 * (keep >= 0))
    if (nrow(keep) != nrow(dat)) {keep = t(keep)}
    rowMeans(S.norm * keep, na.rm = TRUE)
  } else {
    unlist(S.norm)
  }
}

