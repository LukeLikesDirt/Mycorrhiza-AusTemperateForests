# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers for the alpha-diversity bootstrap scripts
# (01a_alpha_diversity_amf.R, 01b_alpha_diversity_g_amf.R,
#  01c_alpha_diversity_m_amf.R, 01d_alpha_diversity_emf.R)
#
# Confirmed byte-identical (bar one blank-line whitespace difference) across
# the AM and ECM scripts before being centralised here — see
# rope_implementation_plan.md. AME helpers (ame_pct(), ame_pct_fast(), etc.)
# were added for the naming-migration + AME spec in edit.pdf, and are used
# only by 01a/01d — see the "Average marginal effects" section below.
# ─────────────────────────────────────────────────────────────────────────────

# Plot theme
theme_custom <- theme_minimal() +
  theme(
    panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
    panel.grid = element_blank(),
    axis.ticks = element_line(colour = "grey80", linewidth = 0.25),
    axis.ticks.length = unit(-0.1, "cm"),
    axis.title = element_markdown(size = 10),
    axis.text = element_markdown(size = 8),
    strip.text = element_markdown(size = 10),
    plot.tag = element_markdown(size = 12),
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1),
    aspect.ratio = 1
  )

# Function to std covariates
std <- function(x) {
  x <- as.numeric(x)
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Function to unstandardise covariates
unstd <- function(z, original_mean, original_sd) {
  (z * original_sd) + original_mean
}

# ─────────────────────────────────────────────────────────────────────────────
# BCa confidence intervals
#
# Refactored per rope_implementation_plan.md §6: the jackknife pass (n model
# refits) used to run once per parameter inside calculate_bca_ci(), which for
# ~9 parameters across 7 models meant ~7,700 redundant refits. One jackknife
# pass yields the estimates for every parameter, so jackknife_coefs() is now
# called once per model in bootstrap_model() and the resulting column handed
# to calculate_bca_ci() per parameter. The BCa maths themselves are unchanged.
# ─────────────────────────────────────────────────────────────────────────────

# One jackknife pass (n model refits) yielding every parameter's estimates.
jackknife_coefs <- function(data, formula_obj, family_obj = NULL, is_negbin = FALSE) {
  n <- nrow(data)
  out <- NULL

  for (i in seq_len(n)) {
    jack_data <- data[-i, ]

    jack_model <- if (is_negbin) {
      MASS::glm.nb(formula_obj, data = jack_data)
    } else if (!is.null(family_obj)) {
      glm(formula_obj, data = jack_data, family = family_obj)
    } else {
      lm(formula_obj, data = jack_data)
    }

    cf <- coef(jack_model)
    if (is.null(out)) {
      out <- matrix(NA_real_, n, length(cf), dimnames = list(NULL, names(cf)))
    }
    out[i, ] <- cf
  }

  out
}

# Function to calculate BCa confidence intervals for one parameter, given a
# pre-computed jackknife column (see jackknife_coefs()).
calculate_bca_ci <- function(boot_estimates, original_estimate, jack_estimates,
                             alpha = 0.05) {

  # Step 1: Calculate bias correction (z0)
  n_boot <- length(boot_estimates)
  prop_less <- sum(boot_estimates < original_estimate) / n_boot
  z0 <- qnorm(prop_less)

  # Step 2: Calculate acceleration (a) from the jackknife estimates
  jack_mean <- mean(jack_estimates)
  numerator <- sum((jack_mean - jack_estimates)^3)
  denominator <- 6 * (sum((jack_mean - jack_estimates)^2))^(3/2)
  a <- numerator / denominator

  # Step 3: Calculate adjusted percentiles
  z_alpha_lower <- qnorm(alpha / 2)
  z_alpha_upper <- qnorm(1 - alpha / 2)

  p_lower <- pnorm(z0 + (z0 + z_alpha_lower) / (1 - a * (z0 + z_alpha_lower)))
  p_upper <- pnorm(z0 + (z0 + z_alpha_upper) / (1 - a * (z0 + z_alpha_upper)))

  # Step 4: Get BCa confidence limits
  ci_lower <- quantile(boot_estimates, p_lower)
  ci_upper <- quantile(boot_estimates, p_upper)

  return(c(ci_lower, ci_upper))
}

# Function to perform bootstrap with both Percentile and BCa CIs.
#
# Returns coefficients$Original alongside coefficients$Mean: Original is
# coef(model)[j] (the point estimate to report — see plan §0 and §9), Mean is
# the bootstrap mean (kept only as an extra diagnostic column, since BCa's
# bias correction z0 is computed against the original estimate, not the
# bootstrap mean, so reporting Mean alongside a BCa interval is inconsistent).
bootstrap_model <- function(model, data, formula_obj, family_obj = NULL,
                            n_iter = 10000, is_negbin = FALSE, ci_level = 0.95,
                            ame_focal = NULL, ame_step = NULL, ame_backtf = NULL,
                            ame_p10 = NULL, ame_p90 = NULL) {

  alpha <- 1 - ci_level
  n <- nrow(data)

  # Average marginal effect (percentage change) is opt-in: when ame_focal is
  # supplied, every bootstrap replicate's own refitted model and resampled
  # data (already in hand for this iteration) are also used to compute the
  # closed-form AME and the two gradient-position predictions, so the interval
  # comes from the same case-resampling bootstrap as every other statistic in
  # this function (rope_implementation_plan.md's sibling spec, "edit.pdf" §3).
  # 01b/01c never pass ame_focal, so they are completely unaffected.
  track_ame <- !is.null(ame_focal)
  if (track_ame) {
    stopifnot(!is.null(ame_step), !is.null(ame_backtf), !is.null(ame_p10), !is.null(ame_p90))
    boot_ame_mean      <- numeric(n_iter)
    boot_ame_n_invalid <- numeric(n_iter)
    boot_pred_p10      <- numeric(n_iter)
    boot_pred_p90      <- numeric(n_iter)
  }

  # Storage for bootstrap coefficients
  coef_names <- names(coef(model))
  n_params <- length(coef_names)
  boot_coefs <- matrix(NA, nrow = n_iter, ncol = n_params)
  colnames(boot_coefs) <- coef_names

  # Storage for R-squared (or pseudo R-squared for GLM)
  boot_r2 <- numeric(n_iter)

  # Storage for partial R² per predictor (r2glmm, lm models only)
  boot_partial_r2    <- NULL   # matrix, initialised on first successful r2beta call
  partial_r2_effects <- NULL   # predictor names, set on first successful call

  # Track convergence issues
  n_failed <- 0

  # Bootstrap loop
  cat("Running bootstrap iterations...\n")
  i <- 1
  attempts <- 0
  max_attempts <- n_iter * 3  # Allow up to 3x iterations to get n_iter successful fits

  while (i <= n_iter && attempts < max_attempts) {
    attempts <- attempts + 1

    # Resample with replacement
    boot_indices <- sample(1:n, size = n, replace = TRUE)
    boot_data <- data[boot_indices, ]

    # Check for NA/NaN/Inf in predictors and response
    predictor_cols <- all.vars(formula_obj)[-1]  # Exclude response variable
    response_col <- all.vars(formula_obj)[1]

    # Check if bootstrap sample has valid data
    # Convert to data.frame for consistent subsetting
    boot_df <- as.data.frame(boot_data)

    if (any(!is.finite(as.matrix(boot_df[, predictor_cols, drop = FALSE]))) ||
        any(!is.finite(boot_df[[response_col]]))) {
      n_failed <- n_failed + 1
      next
    }

    # Check for sufficient variation in predictors
    predictor_sds <- sapply(boot_df[, predictor_cols, drop = FALSE], sd, na.rm = TRUE)
    if (any(predictor_sds == 0 | is.na(predictor_sds))) {
      n_failed <- n_failed + 1
      next
    }

    # Fit model to bootstrap sample with error handling
    converged <- TRUE
    boot_model <- NULL

    boot_model <- tryCatch({
      if (is_negbin) {
        # Negative binomial GLM with stricter controls
        m <- suppressWarnings(
          glm.nb(formula_obj, data = boot_data,
                 control = glm.control(maxit = 100, epsilon = 1e-8, trace = FALSE),
                 init.theta = model$theta,  # Use original theta as starting value
                 link = log)  # Explicit link specification
        )
        # Check convergence and coefficient validity
        if (!m$converged || any(is.na(coef(m))) || any(is.infinite(coef(m))) ||
            any(!is.finite(coef(m)))) {
          converged <- FALSE
          NULL
        } else {
          m
        }
      } else if (!is.null(family_obj)) {
        # GLM (Poisson for richness)
        m <- glm(formula_obj, data = boot_data, family = family_obj,
                 control = glm.control(maxit = 100, trace = FALSE))
        if (!m$converged || any(is.na(coef(m))) || any(is.infinite(coef(m))) ||
            any(!is.finite(coef(m)))) {
          converged <- FALSE
          NULL
        } else {
          m
        }
      } else {
        # LM (for evenness, Shannon, Faith)
        m <- lm(formula_obj, data = boot_data)
        if (any(is.na(coef(m))) || any(is.infinite(coef(m))) ||
            any(!is.finite(coef(m)))) {
          converged <- FALSE
          NULL
        } else {
          m
        }
      }
    }, error = function(e) {
      # Silently catch any errors (including NA/NaN/Inf in x, theta estimation, etc.)
      converged <<- FALSE
      NULL
    })

    # Skip this iteration if model didn't converge or has invalid coefficients
    if (!converged || is.null(boot_model)) {
      n_failed <- n_failed + 1
      next
    }

    # Calculate R-squared with error handling
    r2_calculated <- FALSE
    if (is_negbin) {
      null_model <- tryCatch({
        suppressWarnings(
          glm.nb(update(formula_obj, . ~ 1), data = boot_data,
                 control = glm.control(maxit = 100, epsilon = 1e-8))
        )
      }, error = function(e) NULL)

      if (!is.null(null_model) && null_model$converged &&
          !any(is.na(coef(null_model))) && !any(is.infinite(coef(null_model)))) {
        ll_full <- tryCatch(logLik(boot_model), error = function(e) NULL)
        ll_null <- tryCatch(logLik(null_model), error = function(e) NULL)
        if (!is.null(ll_full) && !is.null(ll_null) && is.finite(ll_full) && is.finite(ll_null)) {
          boot_r2[i] <- as.numeric(1 - (ll_full / ll_null))
          r2_calculated <- TRUE
        }
      }
    } else if (!is.null(family_obj)) {
      null_model <- tryCatch({
        glm(update(formula_obj, . ~ 1), data = boot_data, family = family_obj)
      }, error = function(e) NULL)

      if (!is.null(null_model) && null_model$converged) {
        ll_full <- tryCatch(logLik(boot_model), error = function(e) NULL)
        ll_null <- tryCatch(logLik(null_model), error = function(e) NULL)
        if (!is.null(ll_full) && !is.null(ll_null) && is.finite(ll_full) && is.finite(ll_null)) {
          boot_r2[i] <- as.numeric(1 - (ll_full / ll_null))
          r2_calculated <- TRUE
        }
      }
    } else {
      r2_val <- tryCatch(summary(boot_model)$r.squared, error = function(e) NULL)
      if (!is.null(r2_val) && is.finite(r2_val)) {
        boot_r2[i] <- r2_val
        r2_calculated <- TRUE
      }
    }

    # Skip if R2 calculation failed
    if (!r2_calculated) {
      n_failed <- n_failed + 1
      next
    }

    # Store coefficients
    boot_coefs[i, ] <- coef(boot_model)

    # Average marginal effect (percentage change) for this replicate, using
    # this replicate's own refitted model and resampled data -- see ame_pct_fast()
    if (track_ame) {
      ame_i <- ame_pct_fast(boot_model, boot_data, ame_focal, ame_step, ame_backtf)
      boot_ame_mean[i]      <- ame_i$mean
      boot_ame_n_invalid[i] <- ame_i$n_invalid
      boot_pred_p10[i]      <- avg_pred_backtf(boot_model, boot_data, ame_focal, ame_p10, ame_backtf)
      boot_pred_p90[i]      <- avg_pred_backtf(boot_model, boot_data, ame_focal, ame_p90, ame_backtf)
    }

    # Partial R² via r2glmm (lm models only)
    if (!is_negbin && is.null(family_obj)) {
      r2b <- tryCatch(
        r2glmm::r2beta(boot_model, method = "lm", partial = TRUE),
        error = function(e) NULL
      )
      if (!is.null(r2b)) {
        if (is.null(boot_partial_r2)) {
          partial_r2_effects <- as.character(r2b$Effect)
          boot_partial_r2    <- matrix(NA_real_, nrow = n_iter,
                                       ncol = length(partial_r2_effects))
          colnames(boot_partial_r2) <- partial_r2_effects
        }
        boot_partial_r2[i, ] <- r2b$Rsq
      }
    }

    # Progress indicator
    if (i %% 1000 == 0) {
      cat(sprintf("  %d/%d iterations complete", i, n_iter))
      if (n_failed > 0) {
        cat(sprintf(" (%d failed, %d attempts)\n", n_failed, attempts))
      } else {
        cat("\n")
      }
    }

    # Increment successful iteration counter
    i <- i + 1
  }

  # Check if we got enough successful iterations
  if (attempts >= max_attempts) {
    warning(sprintf("Reached maximum attempts (%d). Only %d/%d successful iterations.",
                    max_attempts, i - 1, n_iter))
  }

  if (n_failed > 0) {
    cat(sprintf("\nBootstrap summary: %d failed iterations out of %d attempts (%.1f%% failure rate)\n",
                n_failed, attempts, 100 * n_failed / attempts))
  }

  # Get original model coefficients — the point estimate to report (plan §0, §9)
  original_coefs <- coef(model)

  # Calculate percentile CIs
  cat("Calculating percentile confidence intervals...\n")
  coef_mean <- apply(boot_coefs, 2, mean)
  coef_ci_lower_pct <- apply(boot_coefs, 2, quantile, probs = alpha / 2)
  coef_ci_upper_pct <- apply(boot_coefs, 2, quantile, probs = 1 - alpha / 2)

  # Calculate BCa CIs — one jackknife pass for the whole model (plan §6),
  # then one calculate_bca_ci() call per parameter using its jackknife column
  cat("Calculating jackknife estimates for BCa...\n")
  jack <- jackknife_coefs(data, formula_obj, family_obj = family_obj, is_negbin = is_negbin)
  stopifnot(setequal(colnames(jack), coef_names))

  cat("Calculating BCa confidence intervals...\n")
  bca_cis <- matrix(NA, nrow = n_params, ncol = 2)

  for (j in 1:n_params) {
    bca_cis[j, ] <- calculate_bca_ci(
      boot_estimates = boot_coefs[, j],
      original_estimate = original_coefs[j],
      jack_estimates = jack[, coef_names[j]],
      alpha = alpha
    )
  }

  # Create summary dataframe with both CI methods
  coef_summary <- data.frame(
    Parameter = coef_names,
    Original = as.numeric(original_coefs),
    Mean = coef_mean,
    Percentile_Lower = coef_ci_lower_pct,
    Percentile_Upper = coef_ci_upper_pct,
    BCa_Lower = bca_cis[, 1],
    BCa_Upper = bca_cis[, 2],
    stringsAsFactors = FALSE
  )

  # Format for display (BCa only)
  coef_summary$CI_95 <- sprintf("[%.2f, %.2f]",
                                coef_summary$BCa_Lower,
                                coef_summary$BCa_Upper)

  # R-squared summary (percentile method for R2)
  r2_mean <- mean(boot_r2)
  r2_ci_lower <- quantile(boot_r2, alpha / 2)
  r2_ci_upper <- quantile(boot_r2, 1 - alpha / 2)

  r2_summary <- data.frame(
    Metric = "R-squared",
    Mean = r2_mean,
    CI_Lower = r2_ci_lower,
    CI_Upper = r2_ci_upper,
    stringsAsFactors = FALSE
  )

  r2_summary$CI_95 <- sprintf("[%.2f, %.2f]",
                              r2_summary$CI_Lower,
                              r2_summary$CI_Upper)

  # Partial R² summary (lm models only)
  partial_r2_summary <- NULL
  if (!is.null(boot_partial_r2)) {
    pr2_mean  <- apply(boot_partial_r2, 2, mean,     na.rm = TRUE)
    pr2_lower <- apply(boot_partial_r2, 2, quantile, probs = alpha / 2, na.rm = TRUE)
    pr2_upper <- apply(boot_partial_r2, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
    partial_r2_summary <- data.frame(
      Effect   = partial_r2_effects,
      Mean_R2  = round(pr2_mean,  3),
      CI_Lower = round(pr2_lower, 3),
      CI_Upper = round(pr2_upper, 3),
      stringsAsFactors = FALSE
    )
    partial_r2_summary$CI_95 <- sprintf("[%.3f, %.3f]",
                                         partial_r2_summary$CI_Lower,
                                         partial_r2_summary$CI_Upper)
  }

  # AME (percentage-change) bootstrap summary: only the interval comes from
  # here -- the point estimate (AME_pct, median, IQR, n_invalid) is computed
  # separately on the *observed* fit by the caller, via ame_pct()/ame_pct_fast()
  # directly ("edit.pdf" §5 wants those columns on the observed fit, not
  # averaged across replicates).
  ame_summary <- NULL
  if (track_ame) {
    ame_ci   <- quantile(boot_ame_mean, c(alpha / 2, 1 - alpha / 2))
    p10_ci   <- quantile(boot_pred_p10, c(alpha / 2, 1 - alpha / 2))
    p90_ci   <- quantile(boot_pred_p90, c(alpha / 2, 1 - alpha / 2))
    pct_change_boot <- (boot_pred_p90 - boot_pred_p10) / boot_pred_p10 * 100
    pct_change_ci   <- quantile(pct_change_boot, c(alpha / 2, 1 - alpha / 2))

    ame_summary <- list(
      ame_pct_ci_lower           = unname(ame_ci[1]),
      ame_pct_ci_upper           = unname(ame_ci[2]),
      n_invalid_total_bootstrap  = sum(boot_ame_n_invalid),
      pred_p10_ci_lower          = unname(p10_ci[1]),
      pred_p10_ci_upper          = unname(p10_ci[2]),
      pred_p90_ci_lower          = unname(p90_ci[1]),
      pred_p90_ci_upper          = unname(p90_ci[2]),
      pct_change_p10_p90_ci_lower = unname(pct_change_ci[1]),
      pct_change_p10_p90_ci_upper = unname(pct_change_ci[2])
    )
  }

  cat("Bootstrap complete!\n\n")

  return(list(
    coefficients    = coef_summary,
    r2              = r2_summary,
    partial_r2      = partial_r2_summary,
    boot_coefs      = boot_coefs,
    boot_r2         = boot_r2,
    boot_partial_r2 = boot_partial_r2,
    ame             = ame_summary
  ))
}

# ─────────────────────────────────────────────────────────────────────────────
# Average marginal effects: percentage change on the response's original
# scale, averaged across the observed data ("edit.pdf" §0-3)
#
# The models are lm(sqrt(y) ~ ...) in most responses, but not all -- EMF's
# richness and Faith's PD are modelled on the raw scale. A constant slope on
# sqrt(y) is not a constant proportional change on y (y = eta^2 is convex), so
# there is no single "% per SD"; only an average over the observed data. The
# back-transform must exactly invert each response's own LHS, so
# response_backtf_for() derives it per response instead of assuming sqrt()
# everywhere.
# ─────────────────────────────────────────────────────────────────────────────

# Back-transform for a response's modelled scale, derived from the formula's
# own LHS rather than assumed: sqrt(y) -> y = eta^2; a bare response column ->
# y = eta (no transform needed -- EMF's richness and Faith's PD). Anything
# else stops loudly rather than silently applying the wrong inverse.
response_backtf_for <- function(formula_obj) {
  lhs <- deparse(formula_obj[[2]])
  if (grepl("^sqrt\\(", lhs)) {
    return(function(eta) eta^2)
  }
  if (grepl("^[A-Za-z_][A-Za-z0-9_.]*$", lhs)) {
    return(function(eta) eta)
  }
  stop(sprintf(
    "response_backtf_for(): unrecognised response transform '%s' in formula %s -- add explicit handling before computing AME for this response.",
    lhs, paste(deparse(formula_obj), collapse = " ")
  ))
}

# Point estimate: average percentage change in the response per `step` units
# of `focal`, on the response's original scale, averaged across the observed
# data via marginaleffects::avg_comparisons() (holds every other covariate at
# its own observed value for each row -- do not evaluate at covariate means).
# hi/lo are predictions on the modelled scale at focal+step and focal;
# back-transforming inside the comparison returns each observation to the
# response scale before the ratio is taken, and avg_comparisons() rebuilds the
# model frame when it perturbs focal, so I(x^2) terms are handled correctly
# (the same reason ame_pct_fast() uses predict() rather than hand-rolled
# model-matrix arithmetic).
#
# newdata is pre-filtered to the same "ok" rows ame_pct_fast() averages over
# (non-positive eta at focal or focal+step, on the modelled scale, is dropped
# there before backtf() is applied -- see its comment). Without this filter,
# avg_comparisons() would silently include those rows (squaring a negative eta
# is still a finite number, so nothing here would error), and its mean would
# disagree with ame_pct_fast()'s whenever n_invalid > 0 -- caught by
# acceptance test 1 in the calling scripts, which requires the two to match to
# 6 significant figures on the observed fit.
ame_pct <- function(model, data, focal, step, backtf) {
  d_lo <- data
  d_hi <- data
  d_hi[[focal]] <- d_hi[[focal]] + step

  eta_lo <- predict(model, newdata = d_lo)
  eta_hi <- predict(model, newdata = d_hi)
  ok <- is.finite(eta_lo) & is.finite(eta_hi) & eta_lo > 0 & eta_hi > 0

  marginaleffects::avg_comparisons(
    model,
    newdata    = data[ok, , drop = FALSE],
    variables  = setNames(list(step), focal),
    comparison = function(hi, lo) (backtf(hi) - backtf(lo)) / backtf(lo) * 100
  )
}

# Closed-form equivalent of ame_pct(), fast enough to call inside a bootstrap
# loop (avg_comparisons() is not -- see acceptance test 1 in the calling
# script, which licenses using this in place of ame_pct() for every bootstrap
# replicate). ok guards against a non-positive fitted value on the *modelled*
# scale (eta) -- for a sqrt() response a negative eta is not a valid square
# root regardless of what its square looks like, so the guard must run before
# backtf(), not after.
ame_pct_fast <- function(fit, data, focal, step, backtf) {
  d_lo <- data
  d_hi <- data
  d_hi[[focal]] <- d_hi[[focal]] + step

  eta_lo <- predict(fit, newdata = d_lo)
  eta_hi <- predict(fit, newdata = d_hi)

  ok <- is.finite(eta_lo) & is.finite(eta_hi) & eta_lo > 0 & eta_hi > 0
  y_lo <- backtf(eta_lo[ok])
  y_hi <- backtf(eta_hi[ok])
  pct  <- (y_hi - y_lo) / y_lo * 100

  list(mean = mean(pct), median = stats::median(pct),
       iqr = stats::IQR(pct), n_invalid = sum(!ok))
}

# Predicted response at a fixed value of `focal` (e.g. its 10th or 90th
# percentile), averaged across observations with every other covariate held
# at its own observed value -- the same "hold covariates at observed values"
# rule as ame_pct(), used for edit.pdf §6's pred_p10/pred_p90. Back-transforming
# the *averaged* eta estimates the median of y on the response scale, not the
# mean (note this when reporting, not when computing).
avg_pred_backtf <- function(fit, data, focal, x_val, backtf) {
  d <- data
  d[[focal]] <- x_val
  eta <- predict(fit, newdata = d)
  ok  <- is.finite(eta) & eta > 0
  mean(backtf(eta[ok]))
}

# ─────────────────────────────────────────────────────────────────────────────
# ROPE-based effect classification (rope_implementation_plan.md §3–5)
# ─────────────────────────────────────────────────────────────────────────────

# SD of the response on the scale it enters the model (e.g. sqrt(richness),
# not richness), read from the model's own frame so NA-dropped rows are
# handled automatically. Do not compute sd(data$richness) directly — wrong
# scale, and may include rows the model itself dropped.
sd_response <- function(model) {
  stats::sd(stats::model.frame(model)[[1]])
}

# Classify a bootstrap CI against a region of practical equivalence (ROPE).
# The four outcomes are mutually exclusive and exhaustive: "unresolved"
# absorbs both a wide interval spanning the ROPE and a narrow interval
# straddling one ROPE boundary.
classify_effect <- function(lower, upper, rope_lower, rope_upper) {
  dplyr::case_when(
    is.na(lower) | is.na(upper) ~ NA_character_,
    upper < rope_lower           ~ "negative",
    lower > rope_upper           ~ "positive",
    lower >= rope_lower & upper <= rope_upper ~ "neutral",
    TRUE                         ~ "unresolved"
  )
}

# Function to extract formatted coefficient for annotation (using BCa),
# fully standardised (SD of modelled response per SD of predictor — plan §0,
# §9). Uses the original model coefficient, not the bootstrap mean (see
# bootstrap_model()'s Original column).
# Linear terms:    β = value [lower, upper]
# Quadratic terms: β₁ = value [lower, upper]  (linear component)
#                  β₂ = value [lower, upper]  (quadratic component)
get_boot_annotation <- function(boot_results, param_name, sd_y) {

  lin_row  <- boot_results$coefficients %>% dplyr::filter(Parameter == param_name)
  quad_row <- boot_results$coefficients %>%
    dplyr::filter(Parameter == paste0("I(", param_name, "^2)"))

  if (nrow(lin_row) > 0 && nrow(quad_row) > 0) {
    b1 <- lin_row$Original / sd_y;  l1 <- lin_row$BCa_Lower / sd_y;  u1 <- lin_row$BCa_Upper / sd_y
    b2 <- quad_row$Original / sd_y; l2 <- quad_row$BCa_Lower / sd_y; u2 <- quad_row$BCa_Upper / sd_y
    return(sprintf("β<sub>1</sub> = %.2f [%.2f, %.2f]<br>β<sub>2</sub> = %.2f [%.2f, %.2f]",
                   b1, l1, u1, b2, l2, u2))
  }

  if (nrow(lin_row) > 0) {
    m   <- lin_row$Original / sd_y
    lcl <- lin_row$BCa_Lower / sd_y
    ucl <- lin_row$BCa_Upper / sd_y
    return(sprintf("β = %.2f [%.2f, %.2f]", m, lcl, ucl))
  }

  return("")
}

# Build a poly_sel data frame from a named logical vector
make_poly_sel <- function(covariates, quad_covs = character(0)) {
  data.frame(
    covariate     = covariates,
    delta_aicc    = NA_real_,
    use_quadratic = covariates %in% quad_covs,
    stringsAsFactors = FALSE
  )
}

# Build GLM/LM formula: linear or quadratic (x + I(x²)) per covariate
build_formula_mixed <- function(response, selection_results) {
  terms <- ifelse(
    selection_results$use_quadratic,
    paste0(selection_results$covariate, " + I(", selection_results$covariate, "^2)"),
    selection_results$covariate
  )
  as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

# Print selection table
print_term_selection <- function(sel, response_name) {
  cat(sprintf("\n--- %s ---\n", response_name))
  print(
    sel %>%
      mutate(
        term = ifelse(use_quadratic, "quadratic [x + I(x²)]", "linear")
      ) %>%
      dplyr::select(covariate, term),
    row.names = FALSE
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# Excel assembly (rope_implementation_plan.md §7-8)
# ─────────────────────────────────────────────────────────────────────────────

# One per-response sheet: every coefficient on both scales plus its ROPE
# classification, with the existing R² and partial-R² blocks appended below.
# Built with bind_rows() on tibbles rather than rbind() on a fixed-width
# character matrix, so a spacer row (tibble(Parameter = "")) fills missing
# columns with NA instead of misaligning them once the table is wider than
# three columns (plan §7.4 — this was the concrete breakage point named in
# the plan).
build_response_sheet <- function(boot, sd_y, rope_std) {

  coefs <- boot$coefficients %>%
    dplyr::transmute(
      Parameter,
      Beta_partial = round(Original, 3),
      CI_partial   = sprintf("[%.3f, %.3f]", BCa_Lower, BCa_Upper),
      Beta_std     = round(Original / sd_y, 3),
      CI_std       = sprintf("[%.3f, %.3f]", BCa_Lower / sd_y, BCa_Upper / sd_y),
      ROPE_std     = sprintf("[%.2f, %.2f]", -rope_std, rope_std),
      Classification = classify_effect(BCa_Lower / sd_y, BCa_Upper / sd_y, -rope_std, rope_std)
    ) %>%
    dplyr::mutate(
      # (Intercept) has no ROPE interpretation (plan §5)
      Classification = dplyr::if_else(Parameter == "(Intercept)", NA_character_, Classification)
    )

  spacer <- tibble::tibble(Parameter = "")

  r2_block <- boot$r2 %>%
    dplyr::transmute(Parameter = Metric, Beta_partial = round(Mean, 3), CI_partial = CI_95)

  pr2_block <- boot$partial_r2 %>%
    dplyr::transmute(Parameter = Effect, Beta_partial = Mean_R2, CI_partial = CI_95)

  dplyr::bind_rows(coefs, spacer, r2_block, spacer, pr2_block)
}

# One row per response, for the focal term only — the table that feeds the
# results text (plan §7.2).
#
# If the focal term carries a quadratic component in this response (plan
# §10), a single linear beta is not a valid summary of "the effect" — neither
# beta1 nor beta2 alone describes it, and the plan explicitly forbids
# reporting either as if it were. The plan gestures at a fix ("a derived
# scalar... for example the predicted change across the observed range,
# computed inside each bootstrap replicate") without fully specifying the
# scalar or what ROPE it should be judged against — both are real
# methodological choices, not mechanical ones, so this reports the gap
# explicitly (NA + a note) rather than inventing an unspecified statistic.
# The per-response sheet still classifies beta1 and beta2 separately for
# transparency (build_response_sheet() always does this).
build_summary_row <- function(label, boot, sd_y, rope_std, focal_term) {
  fr   <- boot$coefficients %>% dplyr::filter(Parameter == focal_term)
  quad <- boot$coefficients %>% dplyr::filter(Parameter == paste0("I(", focal_term, "^2)"))

  if (nrow(quad) > 0) {
    return(tibble::tibble(
      Response       = label,
      Beta_std       = NA_real_,
      CI_std         = NA_character_,
      Classification = "quadratic — see per-response sheet",
      sd_y           = round(sd_y, 4),
      ROPE_partial   = sprintf("[%.3f, %.3f]", -rope_std * sd_y, rope_std * sd_y)
    ))
  }

  tibble::tibble(
    Response       = label,
    Beta_std       = round(fr$Original / sd_y, 3),
    CI_std         = sprintf("[%.3f, %.3f]", fr$BCa_Lower / sd_y, fr$BCa_Upper / sd_y),
    Classification = classify_effect(fr$BCa_Lower / sd_y, fr$BCa_Upper / sd_y, -rope_std, rope_std),
    sd_y           = round(sd_y, 4),
    ROPE_partial   = sprintf("[%.3f, %.3f]", -rope_std * sd_y, rope_std * sd_y)
  )
}

# One row per response, average marginal effect columns (edit.pdf §5-6) --
# joined onto the Summary sheet by Response rather than folded into
# build_summary_row(), so 01b/01c (which never compute AME) are unaffected and
# the ROPE-based Summary sheet stays identical across all four scripts save
# for these extra columns where they apply.
build_ame_row <- function(label, ame_pct, ame_ci_lower, ame_ci_upper,
                          ame_median, ame_iqr, n_invalid, ame_flag,
                          pred_p10 = NA_real_, pred_p10_ci = NA_character_,
                          pred_p90 = NA_real_, pred_p90_ci = NA_character_,
                          pct_change_p10_p90 = NA_real_,
                          pct_change_p10_p90_ci = NA_character_) {
  tibble::tibble(
    Response               = label,
    AME_pct                = round(ame_pct, 1),
    AME_pct_CI              = sprintf("[%.1f, %.1f]", ame_ci_lower, ame_ci_upper),
    AME_pct_median          = round(ame_median, 1),
    AME_pct_IQR             = round(ame_iqr, 1),
    n_invalid_pred          = n_invalid,
    AME_flag                = ame_flag,
    pred_p10                = round(pred_p10, 2),
    pred_p10_CI              = pred_p10_ci,
    pred_p90                = round(pred_p90, 2),
    pred_p90_CI              = pred_p90_ci,
    pct_change_p10_p90      = round(pct_change_p10_p90, 1),
    pct_change_p10_p90_CI    = pct_change_p10_p90_ci
  )
}

# One row per response, full run provenance (plan §7.3).
build_metadata_row <- function(label, formula_obj, n_obs, n_boot, seed, sd_y,
                               rope_std, ci_level, session_info_summary) {
  tibble::tibble(
    response        = label,
    formula         = deparse(formula_obj),
    n_obs           = n_obs,
    n_boot          = n_boot,
    seed            = seed,
    sd_y            = round(sd_y, 4),
    rope_std        = rope_std,
    rope_partial    = round(rope_std * sd_y, 4),
    ci_level        = ci_level,
    interval_method = "BCa",
    date_run        = as.character(Sys.Date()),
    session_info    = session_info_summary
  )
}

# Sensitivity of the focal-term classification to the ROPE half-width and to
# BCa vs percentile intervals (plan §8) — cheap, no refitting, since Beta_std
# is an exact rescaling of the already-bootstrapped partial coefficient.
#
# See build_summary_row() for why a quadratic focal term short-circuits to a
# placeholder instead of a computed sweep — the sweep is specifically over
# ROPE half-widths for a single linear beta, which doesn't exist here.
build_sensitivity_rows <- function(label, boot, sd_y, rope_sens, rope_std, focal_term) {
  fr   <- boot$coefficients %>% dplyr::filter(Parameter == focal_term)
  quad <- boot$coefficients %>% dplyr::filter(Parameter == paste0("I(", focal_term, "^2)"))

  if (nrow(quad) > 0) {
    return(tibble::tibble(
      Response = label,
      Setting  = c(sprintf("BCa, ROPE = %.2f", rope_sens),
                   sprintf("Percentile, ROPE = %.2f", rope_std)),
      Classification = "quadratic — not applicable"
    ))
  }

  bca_rows <- purrr::map_df(rope_sens, function(rp) {
    tibble::tibble(
      Response = label,
      Setting  = sprintf("BCa, ROPE = %.2f", rp),
      Classification = classify_effect(fr$BCa_Lower / sd_y, fr$BCa_Upper / sd_y, -rp, rp)
    )
  })

  pct_row <- tibble::tibble(
    Response = label,
    Setting  = sprintf("Percentile, ROPE = %.2f", rope_std),
    Classification = classify_effect(
      fr$Percentile_Lower / sd_y, fr$Percentile_Upper / sd_y, -rope_std, rope_std
    )
  )

  dplyr::bind_rows(bca_rows, pct_row)
}
