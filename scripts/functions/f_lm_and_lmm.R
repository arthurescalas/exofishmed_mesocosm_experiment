# ------------------------------------------------------------------------------
# fit the different models on OBSERVED data
#   We use random models and the function lmer() from the package lme4
#   site and sample are used as random effects
#  r2 is obtained using the r.squaredGLMM {MuMIn} function and provide two r2:
#   - marginal (m) : variance explained by the fixed effects
#   - conditional (c) : variance explained by the entire model (fixed + rdm effects)
# ------------------------------------------------------------------------------

fit_lmm <- function(input_data, model_list, nm_var_Y = NULL, reml = TRUE) {
  
  res <- list()
  
  if (! is.null(nm_var_Y)) {
    res$mods <- purrr::map(model_list, 
                    ~ lmer(formula(paste0(nm_var_Y , .x)), data = input_data,
                           REML = reml)
                    )
  } else {
    res$mods <- map(model_list, 
                    ~ lmer(formula(.x), data = input_data, REML = reml)
    )
  }
  
  for (i in names(res$mods)) {
    eval(parse(text = paste0(i, " <- res$mods[[i]]")))
  }

  res$summaries <- lapply(res$mods, summary)
  res$anova     <- lapply(res$mods, anova)
  res$rdm_anova <- lapply(res$mods, function(m) ranova(m))
  res$mod_comp  <- eval(parse(text = paste0("anova(", 
                                            paste0(names(res$mods), 
                                                   collapse = ","), ")")))
  res$aic_r2   <- cbind(do.call(rbind, lapply(res$mods, extractAIC)),
                        data.frame(do.call(rbind, lapply(res$mods, r.squaredGLMM)),
                                   row.names = names(res$mods))) %>% 
    setNames(c("df", "AIC","r2m","r2c"))
  
  res$mod_comp <- res$mod_comp[names(model_list), ]
  res$aic_r2   <- res$aic_r2[names(model_list), ]
  
  return(res)
}


extract_best_lmm <- function(X, force_linear = FALSE) {
  
  # keep only the models different from model 1
  mask <- X$mod_comp[, "Pr(>Chisq)"] < 0.05
  mask[is.na(mask)] <- TRUE
  
  if (sum(mask, na.rm = T) != 0) {
    # the best model is the one with the lowest AIC
    tmp  <- X$aic_r2[mask, ]
    bm <- row.names(tmp[which(tmp[, 2] == min(tmp[, 2], na.rm = T)), ])
  } else {
    bm <- "mod1"
  }
  
  if (force_linear) {bm <- "mod1"}
  
  # get the model parameters
  
  # fixed effects
  pars <- rep(NA, 3) %>% setNames(c("(Intercept)", "Rank", "I(Rank^2)"))
  coeff <- fixef(X$mods[[bm]])
  pars[names(coeff)] <- coeff
  names(pars) <- c("intercept", "rank", "rank_squ")
  
  # model summary
  p_pars <- rep(NA, 3) %>% setNames(c("(Intercept)", "Rank", "I(Rank^2)"))
  tmp <- t(X$summaries[[bm]]$coefficients)
  pvals <- tmp["Pr(>|t|)", ]
  p_pars[names(pvals)] <- pvals
  p_pars <- p_pars %>% setNames(c("pars_Pval_intercept", "pars_Pval_rank",
                                   "pars_Pval_rank_squ"))

  # result of ANOVA on fixed effects
  tmp <- t(X$anova[[bm]])
  f <- tmp["F value", ] %>% setNames(paste0("Fval_", names(coeff)[-1]))
  p <- tmp["Pr(>F)",] %>% setNames(paste0("Pval_", names(coeff)[-1]))
  f_and_p <- rep(NA, 4) %>% setNames(c(names(f), names(p)))
  f_and_p[names(f)] <- f
  f_and_p[names(p)] <- p
  f_and_p <- f_and_p %>% setNames(c("Fval_rank", "Fval_rank_squ", "Pval_rank", 
                                    "Pval_rank_squ"))
  
  # random effects
  # ranef(X$mods[[bm]])
  
  # does the random effect improve the model ?
  # if this value is negative then NO, the régular model has a lower AIC
  delta_AIC_lm_lmm <- diff(abs(X$rdm_anova[[bm]][, "AIC"]))
  
  # r square
  mod_eval <- X$aic_r2[bm, ] %>% setNames(c("df", "AIC","r2m","r2c"))
  
  # difference in AIC between the two models: model2 - model1
  delta_models <- apply(X$aic_r2[c("mod1", bm), 2:4], 2, diff) %>% 
    setNames(paste0("delta_", c("AIC","r2m","r2c")))
  
  # final results
  out <- c(best_model = bm, pars, p_pars, f_and_p, 
           delta_AIC_lm_lmm = delta_AIC_lm_lmm, 
           delta_models, mod_eval)
  out <- data.frame(out)
  
  return(out)
}


just_fit_lmm <- function(input_data, model_list, nm_var_Y = NULL, reml = TRUE,
                         wght = NULL) {
  
  res <- list()
  
  if (! is.null(nm_var_Y)) {
    res$mods <- purrr::map(model_list, 
                    ~ lmer(formula(paste0(nm_var_Y , .x)), data = input_data,
                           REML = reml, weights = wght)
    )
  } else {
    res$mods <- map(model_list, 
                    ~ lmer(formula(.x), data = input_data, REML = reml,
                           weights = wght)
    )
  }
  
  for (i in names(res$mods)) {
    eval(parse(text = paste0(i, " <- res$mods[[i]]")))
  }
  
  res$summaries <- lapply(res$mods, summary)
  res$anova     <- lapply(res$mods, anova)
  res$rdm_anova <- lapply(res$mods, function(m) ranova(m))
  res$aic_r2   <- cbind(do.call(rbind, lapply(res$mods, extractAIC)),
                        data.frame(do.call(rbind, lapply(res$mods, r.squaredGLMM)),
                                   row.names = names(res$mods))) %>% 
    setNames(c("df", "AIC","r2m","r2c"))
  
  res$aic_r2   <- res$aic_r2[names(model_list), ]
  
  return(res)
}

# ------------------------------------------------------------------------------
# fit classic regression models on OBSERVED data
# The best mixed model for process and traits is the model 5 and in each case the
#  inclusion of random effect does not improve the model.
# Hence the selected model will be the model 5 without random effect
# ------------------------------------------------------------------------------


fit_lm <- function(input_data, model_list, nm_var_Y) {
  
  res <- list()
  
  if (! is.null(nm_var_Y)) {
    res$mods <- purrr::map(model_list, 
                    ~ lm(formula(paste0(nm_var_Y , .x)), data = input_data)
    )
  } else {
    res$mods <- map(model_list, 
                    ~ lm(formula(.x), data = input_data)
    )
  }

  for (i in names(res$mods)) {
    eval(parse(text = paste0(i, " <- res$mods[[i]]")))
  }
  
  res$summaries <- lapply(res$mods, summary)
  res$anova     <- lapply(res$mods, anova)
  res$mod_comp  <- eval(parse(text = paste0("anova(", 
                                            paste0(names(res$mods), 
                                                   collapse = ","), ")")))
  res$aic_r2   <- cbind(do.call(rbind, lapply(res$mods, extractAIC)),
                        data.frame(do.call(rbind, lapply(res$mods, function(m) {
                          c(summary(m)$r.squared, summary(m)$adj.r.squared)})),
                          row.names = names(res$mods))) %>% 
    setNames(c("df", "AIC","r2m","r2c"))
  
  res$mod_comp <- res$mod_comp[names(model_list), ]
  res$aic_r2   <- res$aic_r2[names(model_list), ]
  
  return(res)
}


extract_lm_results <- function(X, force_linear = FALSE) {
  
  # model names
  nms <- names(X$mods)
  
  # best model
  tmp <- X$aic_r2
  bm <- nms[which(tmp[, "AIC"] == min(tmp[, "AIC"], na.rm = T))]
  
  if (force_linear) {bm <- "mod1"}
  
  # model coefficients
  coeffs <- do.call(c, lapply(nms, function(x) {
    tmp <- coefficients(X$mods[[x]])
    names(tmp) <- paste0(x, "_", names(tmp))
    tmp
  }))
  
  # pvalue of model parameters
  pval_coeff <- do.call(c, lapply(nms, function(x) {
    tmp <- X$summaries[[x]]$coefficients[, "Pr(>|t|)"]
    names(tmp) <- paste0(x, "_pval_", names(tmp))
    tmp
  }))
  
  # result of ANOVA : pvalues of model variables
  aov_pval_coeff <- do.call(c, lapply(nms, function(x) {
    tmp <- X$anova[[x]]
    tmp <- tmp[! rownames(tmp) %in% "Residuals", ]
    res <- tmp[, "Pr(>F)"]
    names(res) <- paste0(x, "_aov_pval_", rownames(tmp))
    res
  }))
  
  # result of ANOVA : F values of model variables
  aov_fval_coeff <- do.call(c, lapply(nms, function(x) {
    tmp <- X$anova[[x]]
    tmp <- tmp[! rownames(tmp) %in% "Residuals", ]
    res <- tmp[, "F value"]
    names(res) <- paste0(x, "_aov_fval_", rownames(tmp))
    res
  }))
  
  # model evaluation
  mod_eval <- do.call(cbind, lapply(nms, function(x) {
    tmp <- X$aic_r2[x, ]
    names(tmp) <- paste0(x, "_", names(tmp))
    tmp
  }))
  
  # model comparison
  # mod_compar <- do.call(c, lapply(combn(nms, 2, simplify = FALSE), function(xx) {
  #   tmp <- apply(X$aic_r2[xx,c("AIC","r2m","r2c")], 2, diff)
  #   names(tmp) <- paste0("delta_", names(tmp), "_", paste0(xx, collapse = "_"))
  #   tmp
  # }))
  
  # final results
  out <- c(best_model = bm, coeffs, pval_coeff, aov_pval_coeff, aov_fval_coeff,
           mod_eval)
  names(out) <- gsub("I\\(", "", names(out))
  names(out) <- gsub("\\(", "", names(out))
  names(out) <- gsub("\\)", "", names(out))
  names(out) <- gsub("\\^2", "sq", names(out))
  out <- data.frame(out)

return(out)
}

