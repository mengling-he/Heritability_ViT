# Code Update Notes

**Updated:** 2026-06-01

The following files were cleaned up for simplification and efficiency. Logic and results are unchanged.

## tuning_nointer.R
- Replaced `for` loop over lambda grid with `vapply`.

## pca_fun.R
- `Reduce(c, Y_list)` → `unlist(Y_list)`.
- `ee.vec` and `Xmat.tensor` no longer grown element-by-element in a loop; collected in lists and bound once after the loop.
- Per-subject `eval.basis(T_list[[i]], ...)` call removed; rows are sliced directly from the already-built `Xmat`. `T_list` parameter removed from the function signature.

## pca_score.R
- Removed `K.pc <- K.est` alias; use `K.est` directly.
- Replaced dense `(n·K) × (n·K)` matrix `G` with a vector `g = rep(v1.est, n)`; `solve(G)` replaced by in-place diagonal update.
- Eliminated the N×N matrix `S.inv`; `S.inv %*% e` is now computed directly via the Woodbury formula as a vector operation.

## fpca_analysis_single.R
- Removed unused variables: `knots`, `ylimit`, `df_mu`, `df_phi`, `label_map`.
- `Reduce(c, Tlist)` / `Reduce(c, Ylist)` → `unlist(...)`.
- `eval.basis` results cached per subject during the `Xmat` build loop (`BS_by_id`); the separate second loop for `Phi_by_id_orig` is replaced by `lapply(BS_by_id, ...)`.
- `phi_grid` column truncation moved before the matrix multiply: `BS %*% V1[, seq_len(k_phi)]`.
- Redundant `paste0` labels in `df_mean_pm$curve` replaced with already-computed `lbl_pos`/`lbl_neg`.
- `pca_fun` call updated to drop the now-unused `T_list` argument.
