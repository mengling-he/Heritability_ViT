# Choose the optimal tuning parameter when there is no fixed effects
#
# This uses a simple GCV criterion for penalized least squares:
#   bhat = (X'X + lambda * Omega)^{-1} X'Y
#
# Notes:
# - Omega should be a KxK penalty matrix.
# - Prior versions mistakenly used adiag(Omega*lambda). The correct penalty is
#   (lambda * Omega).

tuning_nointer <- function(lower, upper, Omega, Xmat, Y_vec) {
  lam.list <- exp(seq(lower, upper, 1))
  XtX <- t(Xmat) %*% Xmat
  XtY <- t(Xmat) %*% Y_vec
  nobs <- nrow(Xmat)
  gcv <- vapply(lam.list, function(lam) {
    A <- solve(XtX + lam * Omega)
    trH <- sum(diag(XtX %*% A))
    mean((Y_vec - Xmat %*% (A %*% XtY))^2) / (1 - trH / nobs)^2
  }, numeric(1))
  lam.list[which.min(gcv)]
}
