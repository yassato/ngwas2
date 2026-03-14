simPhenoWithCov <- function(ZETA, covList, herit,
                            weights, rhosMat,
                            X = NULL, beta = NULL) {
  weights <- weights / sum(weights)
  nKernels <- length(weights)
  nDims <- sapply(
    X = ZETA,
    FUN = function(x) {
      nrow(x$K)
    }, simplify = TRUE
  )

  kernelNamesEach <- unlist(sapply(
    X = ZETA,
    FUN = function(x) {
      rownames(x$K)
    }, simplify = FALSE
  ), use.names = FALSE)

  nDimAll <- sum(nDims)
  cumDims <- c(0, cumsum(nDims))
  allKernels <- matrix(
    data = 0,
    nrow = nDimAll,
    ncol = nDimAll
  )
  covKernels <- rhosMat

  for (i in 1:nKernels) {
    rangeI <- (cumDims[i] + 1):cumDims[i + 1]
    for (j in 1:nKernels) {
      rangeJ <- (cumDims[j] + 1):cumDims[j + 1]
      covNow <- sqrt(weights[i] * weights[j]) * covKernels[i, j]
      covKernels[i, j] <- covNow

      if (i == j) {
        allKernels[rangeI, rangeJ] <- covNow * ZETA[[i]]$K
      } else {
        allKernels[rangeI, rangeJ] <- covNow * covList[[i]][[j]]
      }
    }
  }

  ranEffs <- MASS::mvrnorm(n = 1, mu = rep(0, nDimAll),
                           Sigma = allKernels)
  ranEffs <- as.matrix(ranEffs)
  rownames(ranEffs) <- kernelNamesEach

  Zs <- do.call(
    what = cbind,
    args = lapply(
      X = ZETA,
      FUN = function(x) {
        x$Z
      })
  )
  uAll <- Zs %*% ranEffs

  varU <- var(uAll)
  varE <- varU * (1 - herit) / herit

  e <- rnorm(n = nrow(uAll), mean = 0, sd = sqrt(varE))

  if (is.null(X)) {
    y <- uAll + e
  } else {
    y <- X %*% beta + uAll + e
  }

  return(list(y = y,
              uAll = uAll,
              uEach = ranEffs,
              varE = varE,
              varU = varU,
              e = e,
              herit = herit,
              weights = weights,
              rhosMat = rhosMat,
              beta = beta))
}
