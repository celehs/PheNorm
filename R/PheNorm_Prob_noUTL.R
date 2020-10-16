#' Fit the phenotyping algorithm PheNorm using EHR features
#'
#' @description
#' The function requires as input a surrogate, such as the ICD code
#' It can leverage other EHR features (optional) to assist risk prediction.
#'
#' @param nm.logS.ori name of the surrogates (log(ICD+1), log(NLP+1) and log(ICD+NLP+1)
#' @param dat all data columns need to be log-transformed and need column names
#' @param nm.X additional features other than the main ICD and NLP
#' @param corrupt.rate rate for random corruption denoising, between 0 and 1
#' @param train.size size of training sample, default value 10 * nrow(dat)
#' @return list containing probability and beta coefficient
#' @export
PheNorm.Prob_noUTL <- function(nm.logS.ori, dat, nm.X = NULL, corrupt.rate = 0.3, train.size = 10 * nrow(dat)) {
  dat <- as.matrix(dat)
  S.ori <- dat[, nm.logS.ori, drop = FALSE]
  S.norm <- S.ori
  if (!is.null(nm.X)) {
    X <- as.matrix(dat[, nm.X])
    SX.norm <- cbind(S.norm, X)
    id <- sample(1:nrow(dat), train.size, replace = TRUE)
    SX.norm.corrupt <- apply(SX.norm[id,], 2,
                             function(x) {ifelse(rbinom(length(id), 1, corrupt.rate), mean(x), x)}
                             )
    b.all <- apply(S.norm, 2, function(ss) {lm(ss[id] ~ SX.norm.corrupt - 1)$coef})
    b.all[is.na(b.all)] <- 0
    S.norm <- as.matrix(SX.norm) %*% b.all
    b.all <- b.all[-dim(b.all)[1], ]
  } else {
    b.all = NULL
  }
  if (length(nm.logS.ori) > 1) {
    postprob <- apply(S.norm, 2,
                      function(x) {
                        fit <- normalmixEM2comp2(x, lambda = 0.5,
                                                 mu = quantile(x, probs = c(1/3, 2/3)), sigsqrd = 1
                                                 )
                        fit$posterior[, 2]
                        }
                      )
    list("probs" = rowMeans(postprob, na.rm = TRUE), "betas" = b.all)

  } else {
    fit <- normalmixEM2comp2(unlist(S.norm), lambda = 0.5,
                             mu = quantile(S.norm, probs = c(1/3, 2/3)), sigsqrd = 1
                             )
    list("probs" = fit$posterior[, 2], "betas" = b.all)
  }
}

