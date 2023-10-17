## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 6, warning = FALSE)

## ----load_packages------------------------------------------------------------
library(flashier)
library(ggplot2)
library(cowplot)

## ----gtex.const---------------------------------------------------------------
data(gtex, gtex_colors)
fit_default <- flash(gtex, greedy_Kmax = 5)

## ----gtex.output--------------------------------------------------------------
fit_default$elbo

## ----plot.fn------------------------------------------------------------------
plot(fit_default, pw_which = "factors", pm_colors = gtex_colors)

## ----gtex.bycol---------------------------------------------------------------
fit_var_bycol <- flash(gtex, greedy_Kmax = 5, var_type = 2, verbose = 0)
c(const_elbo = fit_default$elbo, bycol_elbo = fit_var_bycol$elbo)

## ----gtex.kron----------------------------------------------------------------
fit_var_kronecker <- flash(gtex, greedy_Kmax = 5, var_type = c(1, 2), verbose = 0)
c(bycol_elbo = fit_var_bycol$elbo, kron_elbo = fit_var_kronecker$elbo)

## ----gtex.bycol.plot----------------------------------------------------------
p1 <- plot(
  fit_var_bycol, 
  pm_which = "factors", 
  pm_colors = gtex_colors, 
  order_by_pve = FALSE,
  include_scree = FALSE
) + ggtitle("By-column variance assumption")

p2 <- plot(
  fit_var_kronecker, 
  pm_which = "factors", 
  pm_colors = gtex_colors, 
  order_by_pve = FALSE,
  include_scree = FALSE
) + ggtitle("Kronecker variance assumption")

plot_grid(p1, p2, nrow = 2)

## ----var.compare--------------------------------------------------------------
comparison_df <- data.frame(
  var_by_col = as.vector(ldf(fit_var_bycol, type = "m")$F),
  var_kronecker = as.vector(ldf(fit_var_kronecker, type = "m")$F),
  k = factor(rep(1:5, each = ncol(gtex)))
)
ggplot(comparison_df, aes(x = var_by_col, y = var_kronecker, col = k)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "By-column variance assumption", y = "Kronecker variance assumption") +
  theme_minimal() 

## ----fit_var_bycol_plus1------------------------------------------------------
fit_with_meas_err <- flash(gtex, S = 1, greedy_Kmax = 5, var_type = 2, verbose = 0)
c(bycol_elbo = fit_var_bycol$elbo, meas_err_elbo = fit_with_meas_err$elbo)

## ----gtex.ebnm.fn-------------------------------------------------------------
fit_normalmix <- flash(
  gtex, 
  greedy_Kmax = 5, 
  ebnm_fn = ebnm_normal_scale_mixture, 
  verbose = 0
)
fit_unimodal <- flash(
  gtex, 
  greedy_Kmax = 5, 
  ebnm_fn = ebnm_unimodal_symmetric, 
  verbose = 0
)
c(pn_elbo = fit_default$elbo, 
  normalmix_elbo = fit_normalmix$elbo, 
  unimodal_elbo = fit_unimodal$elbo)

## ----gtex.snn-----------------------------------------------------------------
fit_seminonnegative <- flash(
  gtex,
  greedy_Kmax = 5, 
  var_type = 2,
  ebnm_fn = c(ebnm_point_normal, ebnm_point_exponential), 
  verbose = 0
)

## -----------------------------------------------------------------------------
plot(
  fit_seminonnegative,
  pm_which = "factors", 
  pm_colors = gtex_colors, 
  include_scree = FALSE
)

## ----backfit------------------------------------------------------------------
fit_greedy <- flash(
  gtex, 
  greedy_Kmax = 5, 
  var_type = 2,
  ebnm_fn = c(ebnm_point_normal, ebnm_point_exponential), 
  backfit = FALSE, 
  verbose = 0
)
fit_backfit <- flash(
  gtex, 
  greedy_Kmax = 5, 
  var_type = 2,
  ebnm_fn = c(ebnm_point_normal, ebnm_point_exponential), 
  backfit = TRUE, 
  verbose = 0
)

c(greedy_elbo = fit_greedy$elbo, bf_elbo = fit_backfit$elbo)

## ----backfit.plot-------------------------------------------------------------
p1 <- plot(
  fit_greedy, 
  pm_which = "factors",
  pm_colors = gtex_colors, 
  order_by_pve = FALSE,
  include_scree = FALSE
) + ggtitle("Greedy fit")
p2 <- plot(
  fit_backfit, 
  pm_which = "factors",
  pm_colors = gtex_colors, 
  order_by_pve = FALSE,
  include_scree = FALSE
) + ggtitle("Backfit")

plot_grid(p1, p2, nrow = 2)

## ----bf.compare---------------------------------------------------------------
comparison_df <- data.frame(
  fit_greedy = as.vector(ldf(fit_greedy, type = "m")$F),
  fit_backfit = as.vector(ldf(fit_backfit, type = "m")$F),
  k = factor(rep(1:5, each = ncol(gtex)))
)
ggplot(comparison_df, aes(x = fit_greedy, y = fit_backfit, col = k)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(x = "Greedy fit", y = "Backfit") +
  theme_minimal() 

## ----final.bf.ci--------------------------------------------------------------
# Set seed for reproducibility.
set.seed(1)
# Use returned sampler to sample from posterior.
samp <- fit_backfit$sampler(nsamp = 200)
# Only keep factor 3.
factor3_samp <- lapply(samp, function(x) x[[2]][, 3])
# Normalize the loadings.
factor3_samp <- sapply(factor3_samp, function(x) x / max(abs(x)))
# Get 95% confidence intervals.
factor3_ci <- apply(factor3_samp, 1, quantile, c(0.025, 0.975))

## ----final.bf.ci.plot---------------------------------------------------------
plot(fit_backfit, kset = 3, pm_colors = gtex_colors, include_scree = FALSE) +
  geom_errorbar(aes(ymin = factor3_ci[1, ], ymax = factor3_ci[2, ]))

## -----------------------------------------------------------------------------
sessionInfo()

