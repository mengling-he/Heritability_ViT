# Uni-level FPCA with irregular time points
#
# This function estimates:
# 1) The smooth covariance surface of residuals (G1)
# 2) The smooth variance function (for sigma_e^2)
# 3) Eigenvalues/eigenfunctions via basis-space eigen-decomposition
#
# Key fix vs older versions:
# - Use (X'X + lambda * Omega) rather than adiag(Omega*lambda)

pca_fun <- function(Y_list, Xmat, bhat, K, J, BS, basis, Omega, Omega2,
                    lower_g1 = -10, upper_g1 = 10,
                    lower_var = -10, upper_var = 15,
                    eps_sigma = 1e-8) {

  n <- length(Y_list)
  Y_vec <- unlist(Y_list)
  e <- as.numeric(Y_vec - Xmat %*% bhat)

  # --- Estimate G1 (off-diagonal) ---
  ee.list <- vector('list', n)
  Xmat.tensor.list <- vector('list', n)
  start.temp <- 1
  for (i in seq_len(n)) {
    n.i <- length(Y_list[[i]])
    e.i <- e[start.temp:(start.temp + n.i - 1)]
    BS.i <- Xmat[start.temp:(start.temp + n.i - 1), ]
    start.temp <- start.temp + n.i

    # response: vec(e_i e_i^T) excluding diagonal
    ee.i <- kronecker(e.i, e.i)
    idx_diag <- (0:(n.i - 1)) * (n.i + 1) + 1
    ee.list[[i]] <- ee.i[-idx_diag]

    # covariates: kron(BS_i, BS_i) excluding diagonal rows
    BS.i.tensor <- kronecker(BS.i, BS.i)
    Xmat.tensor.list[[i]] <- BS.i.tensor[-idx_diag, , drop = FALSE]
  }
  ee.vec <- do.call(c, ee.list)
  Xmat.tensor <- do.call(rbind, Xmat.tensor.list)

  lam.G1 <- tuning_nointer(lower_g1, upper_g1, Omega2, Xmat.tensor, ee.vec)
  bhat.G1 <- solve(t(Xmat.tensor) %*% Xmat.tensor + lam.G1 * Omega2) %*%
    (t(Xmat.tensor) %*% ee.vec)

  # --- Estimate variance function (for sigma_e^2) ---
  e2.vec <- e^2
  lam.var <- tuning_nointer(lower_var, upper_var, Omega, Xmat, e2.vec)
  bhat.e2 <- solve(t(Xmat) %*% Xmat + lam.var * Omega) %*% (t(Xmat) %*% e2.vec)

  # --- Eigenvalues/eigenfunctions ---
  Jmat <- inprod(basis, basis)
  Jmat.sqrt <- chol(Jmat)

  A1 <- matrix(bhat.G1, K, K)
  fpca1 <- eigen(Jmat.sqrt %*% A1 %*% t(Jmat.sqrt), symmetric = TRUE)
  v1 <- fpca1$values
  V1 <- solve(Jmat.sqrt) %*% fpca1$vectors

  # estimated covariance surface on grid (used for sigma_e2)
  BS.tensor <- kronecker(BS, BS)
  K.hat.vec <- BS.tensor %*% bhat.G1
  K.hat <- matrix(K.hat.vec, J, J)

  sigma.e2.hat <- mean(as.numeric(BS %*% bhat.e2) - diag(K.hat))
  sigma.e2.hat <- max(sigma.e2.hat, eps_sigma)

  list(
    v1 = v1,
    V1 = V1,
    sigma_e2 = sigma.e2.hat,
    lam_G1 = lam.G1,
    lam_var = lam.var,
    bhat_G1 = bhat.G1,
    bhat_e2 = bhat.e2
  )
}
