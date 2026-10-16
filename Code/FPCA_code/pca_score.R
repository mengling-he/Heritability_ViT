# Principal component scores prediction (BLUP)
#
# Key fixes vs older versions:
# - Accept Y_vec explicitly (no global dependency)
# - Define n before allocating PC
# - Add numerical safeguards for sigma_e2

pca_score <- function(Xmat, bhat, Y_vec, T_list, K.est, v1.est, Phi1.est, sigma.e2.hat,
                      eps_sigma = 1e-8) {

  n <- length(T_list)
  N <- nrow(Xmat)

  sigma2 <- max(as.numeric(sigma.e2.hat), eps_sigma)
  g <- rep(v1.est, n)

  # Eigenfunction matrix Phi (block-diagonal by subject)
  Z <- matrix(0, N, n * K.est)
  start.temp <- 1
  for (i in seq_len(n)) {
    n.i <- length(T_list[[i]])
    Z[start.temp:(start.temp + n.i - 1), ((i - 1) * K.est + 1):(i * K.est)] <- Phi1.est[[i]]
    start.temp <- start.temp + n.i
  }

  # Sigma^{-1} using Woodbury identity: (sigma^2 I + Z G Z')^{-1}
  M <- crossprod(Z) / sigma2
  diag(M) <- diag(M) + rep(1 / v1.est, n)
  mid <- solve(M)

  # BLUP: avoid forming the full N x N matrix S.inv
  e <- as.numeric(Y_vec - Xmat %*% bhat)
  Zt_e <- t(Z) %*% e
  Se <- Zt_e / sigma2 - t(Z) %*% (Z %*% (mid %*% Zt_e)) / sigma2^2
  xi <- g * Se
  PC <- t(matrix(xi, K.est, n))

  colnames(PC) <- paste0('fpca', seq_len(K.est))
  rownames(PC) <- names(T_list)
  PC
}
