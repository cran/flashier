## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 6, warning = FALSE)

## ----load_packages, message=FALSE---------------------------------------------
library(flashier)
library(ggplot2)
library(cowplot)
library(dplyr)

## ----adv_interface------------------------------------------------------------
# # Basic interface (not run):
# fit_backfit <- flash(
#     gtex,
#     greedy_Kmax = 5,
#     var_type = 2,
#     backfit = TRUE,
#     verbose = 0
#   )

# Pipeable interface:
t_backfit <- system.time(
  fit_backfit <- flash_init(gtex, var_type = 2) %>%
    flash_set_verbose(verbose = 0) %>%
    flash_greedy(Kmax = 5) %>%
    flash_backfit() %>%
    flash_nullcheck()
)

## ----cust_order_ops-----------------------------------------------------------
# Pipeable interface:
fit_multiple_backfits <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(verbose = 0) %>%
  flash_greedy(Kmax = 3) %>%
  flash_backfit() %>%
  flash_nullcheck() %>%
  flash_greedy(Kmax = 2) %>%
  flash_backfit() %>%
  flash_nullcheck()

c(one_bf_elbo = fit_backfit$elbo, two_bf_elbo = fit_multiple_backfits$elbo)

## ----init.factors-------------------------------------------------------------
fit_alternative_backfit <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(verbose = 0) %>%
  flash_factors_init(svd(gtex, nu = 5, nv = 5)) %>%
  flash_backfit(verbose = 0)
c(bf_elbo = fit_backfit$elbo, alt_bf_elbo = fit_alternative_backfit$elbo)

## ----no.extrap----------------------------------------------------------------
t_no_extrapolate <- system.time(
  fit_no_extrapolate <- flash_init(gtex, var_type = 2) %>%
    flash_set_verbose(verbose = 0) %>%
    flash_greedy(Kmax = 5) %>%
    flash_backfit(extrapolate = FALSE) %>%
    flash_nullcheck()
)
c(extrapolate_elbo = fit_backfit$elbo, no_extrapolate_elbo = fit_no_extrapolate$elbo)

## ----no.extrap.time-----------------------------------------------------------
c(t_extrapolate = t_backfit[3], t_no_extrapolate = t_no_extrapolate[3])

## ----intercept----------------------------------------------------------------
fit_with_intercept <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(verbose = 0) %>%
  flash_add_intercept(rowwise = FALSE) %>%
  flash_greedy(Kmax = 4) %>%
  flash_backfit() %>%
  flash_nullcheck()

p1 <- plot(
  fit_backfit, 
  pm_which = "factors", 
  pm_colors = gtex_colors,
  include_scree = FALSE
) + ggtitle("No intercept")
p2 <- plot(
  fit_with_intercept, 
  pm_which = "factors", 
  pm_colors = gtex_colors,
  include_scree = FALSE
) + ggtitle("With intercept")
plot_grid(p1, p2, nrow = 2)

## ----fix.mean, eval = FALSE---------------------------------------------------
#  ones <- matrix(1, nrow = ncol(gtex), ncol = 1)
#  init_loadings <- matrix(rowMeans(gtex), ncol = 1)
#  
#  fit_with_intercept <- flash_init(gtex, var_type = 2) %>%
#    flash_set_verbose(0) %>%
#    flash_factors_init(list(init_loadings, ones)) %>%
#    flash_factors_fix(kset = 1, which_dim = "factors") %>%
#    flash_greedy(Kmax = 4) %>%
#    flash_backfit()

## ----fixed.sprs---------------------------------------------------------------
is_brain <- grepl("Brain", colnames(gtex))
init_loadings <- rowMeans(gtex[, is_brain]) - rowMeans(gtex[, !is_brain])

fit_fixed_pattern <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(0) %>%
  flash_add_intercept(rowwise = FALSE) %>%
  flash_factors_init(list(matrix(init_loadings, ncol = 1),
                          matrix(is_brain, ncol = 1))) %>%
  flash_factors_fix(kset = 2, 
                    which_dim = "factors", 
                    fixed_idx = !is_brain) %>%
  flash_greedy(3) %>%
  flash_backfit()

plot(
  fit_fixed_pattern, 
  pm_which = "factors",
  pm_colors = gtex_colors, 
  include_scree = FALSE,
  order_by_pve = FALSE
)

## ----conv.crit----------------------------------------------------------------
gtex_conv_crit <- flash_init(gtex, var_type = 2) %>%
  flash_set_conv_crit(fn = flash_conv_crit_max_chg_F, tol = .001) %>%
  flash_set_verbose(
    fns = c(flash_verbose_elbo, flash_verbose_max_chg_F),
    colnames = c("ELBO", "Max.Chg.Factors"),
    colwidths = c(18, 18)
  ) %>%
  flash_greedy(Kmax = 3) %>%
  flash_backfit()

## ----custom-------------------------------------------------------------------
verbose_sparsity <- function(new, old, k, f_idx) {
  g <- flash_fit_get_g(new, n = 2) # setting n = 2 gets g_f (n = 1 would get g_\ell)
  pi0 <- g[[f_idx]]$pi[1] # return point mass weight
  return(formatC(pi0, format = "f", digits = 3)) 
}
verbose_sprs2 <- function(new, old, k) verbose_sparsity(new, old, k, 2)
verbose_sprs3 <- function(new, old, k) verbose_sparsity(new, old, k, 3)
verbose_sprs4 <- function(new, old, k) verbose_sparsity(new, old, k, 4)
verbose_sprs5 <- function(new, old, k) verbose_sparsity(new, old, k, 5)

fit_monitor_sparsity <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(0) %>%
  flash_greedy(Kmax = 5) %>%
  flash_set_verbose(
    verbose = 3,
    fns = c(flash_verbose_elbo, verbose_sprs2, verbose_sprs3, verbose_sprs4, verbose_sprs5),
    colnames = c("ELBO", paste0("Sparsity (", 2:5, ")")),
    colwidths = rep(14, 5)
  ) %>%
  flash_backfit()

## ----normal.est.mode----------------------------------------------------------
fit_flash_ebnm <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(0) %>%
  flash_greedy(ebnm_fn = flash_ebnm(prior_family = "normal", mode = "estimate")) %>%
  flash_greedy(Kmax = 4, ebnm_fn = ebnm_point_normal)

fit_flash_ebnm$F_ghat[[1]]

## ----custom.ebnm--------------------------------------------------------------
ebnm_custom <- function(x, s, g_init, fix_g, output) {
  if (fix_g) {
    ebnm_res <- ebnm_ash(
      x, s, g_init = g_init, fix_g = TRUE, output = output,
      mixcompdist = "normal"
    )
  } else {
    # Parameters are:
    #   1. mean of normal component
    #   2. sd of normal component
    neg_llik <- function(par) {
      g <- ashr::normalmix(c(0.5, 0.5), c(0, par[1]), c(0, par[2]))
      ebnm_res <- ebnm_ash(
        x, s, g_init = g, fix_g = FALSE, mixcompdist = "normal"
      )
      return(-ebnm_res$log_likelihood)
    }
    
    # Run optim to get mean and sd of normal component:
    opt_res <- optim(
      par = c(0, 1), # Initial values
      fn = neg_llik, 
      method = "L-BFGS-B", 
      lower = c(-Inf, 0.01), 
      upper = c(Inf, Inf)
    )
    
    # Now re-run ash to get mixture weights:
    opt_par <- opt_res$par
    g <- ashr::normalmix(c(0.5, 0.5), c(0, opt_par[1]), c(0, opt_par[2]))
    ebnm_res <- ebnm_ash(
        x, s, g_init = g, fix_g = FALSE, output = output,
        mixcompdist = "normal"
    )
  } 
  
  return(ebnm_res)
}

fit_custom <- flash_init(gtex, var_type = 2) %>%
  flash_set_verbose(0) %>%
  flash_greedy(
    Kmax = 2,
    ebnm_fn = c(ebnm_point_normal, ebnm_custom)
  )

fit_custom$F_ghat

## ----plot.history, eval = FALSE-----------------------------------------------
#  sink("zz.tsv")
#  tmp <- flash_init(gtex, var_type = 2) %>%
#    flash_set_verbose(-1) %>%
#    flash_factors_init(svd(gtex, nu = 5, nv = 5)) %>%
#    flash_backfit()
#  progress_extrapolate <- read.delim("zz.tsv")
#  sink()
#  
#  sink("zz.tsv")
#  tmp <- flash_init(gtex, var_type = 2) %>%
#    flash_set_verbose(-1) %>%
#    flash_factors_init(svd(gtex, nu = 5, nv = 5)) %>%
#    flash_backfit(extrapolate = FALSE)
#  progress_no_extrapolate <- read.delim("zz.tsv")
#  sink()
#  
#  rm(tmp)
#  file.remove("zz.tsv")
#  
#  progress_extrapolate <- progress_extrapolate %>%
#    mutate(Extrapolate = TRUE) %>%
#    select(Iter, ELBO, Extrapolate)
#  
#  progress_no_extrapolate <- progress_no_extrapolate %>%
#    group_by(Iter) %>%
#    summarize(ELBO = max(ELBO, na.rm = TRUE)) %>%
#    ungroup() %>%
#    mutate(Extrapolate = FALSE)
#  
#  tib <- progress_extrapolate %>%
#    bind_rows(progress_no_extrapolate) %>%
#    mutate(Iter = as.numeric(Iter),
#           ELBO = as.numeric(ELBO))
#  
#  ggplot(tib, aes(x = Iter, y = ELBO, col = Extrapolate)) +
#    geom_line() +
#    theme_minimal()

## -----------------------------------------------------------------------------
sessionInfo()

