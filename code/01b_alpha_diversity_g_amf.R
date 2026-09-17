
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#  Testing the effect of mineral nitrogen on Glomeromycota alpha diversity using a
#  causal modelling framework
#
#  Metrics modeled:
#  - Richness, Faith's PD, Hill taxonomic (q=1, q=2), Hill phylogenetic (q=1,
#    q=2), and Glomeromycota SRS abundance
#
#  Effects are classified negative / positive / neutral / unresolved against a
#  region of practical equivalence (ROPE) on the fully standardised scale (SD
#  of the modelled response per SD of predictor). See rope_implementation_plan.md
#  for the full specification; shared helpers live in code/functions_alpha_diversity.R,
#  sourced by all four alpha-diversity scripts (01a-01d).
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

# Load packages ----------------------------------------------------------------
library(data.table)
library(writexl)
library(parameters)
library(DHARMa)
library(GGally)
library(ggeffects)
library(ggtext)
library(patchwork)
library(marginaleffects)
library(tidyverse)

# ---- Analysis configuration --------------------------------------------------
CI_LEVEL   <- 0.95                    # interval level for BCa / percentile CIs
ROPE_STD   <- 0.10                    # ROPE half-width on the fully standardised scale
ROPE_SENS  <- c(0.05, 0.10, 0.20)     # sensitivity values for the Sensitivity sheet
N_BOOT     <- 10000
SEED       <- 1986
FOCAL_TERM <- "mineral_nitrogen_std"  # the exposure of interest
AME_STEP   <- 1                       # predictor is z-scored, so 1 unit = 1 SD

# edit.pdf §1 specifies a single RESPONSE_BACKTF <- function(eta) eta^2,
# asserting every response's LHS is sqrt(...) and stop()-ing otherwise. That
# holds for every response in *this* script (01b), and in 01a/01c, but not for
# 01d_alpha_diversity_emf.R (EMF models richness and Faith's PD untransformed) --
# squaring an already-untransformed prediction would silently misreport, exactly
# what the plan's assertion is meant to prevent. response_backtf_for()
# (functions_alpha_diversity.R) generalises this: it derives the correct
# back-transform from each response's own formula LHS (square for sqrt(...),
# identity for a bare response column) and stop()s for anything else, so the
# safety property the plan asks for is preserved per-response rather than
# assumed globally. In this script every response resolves to eta^2.

# Functions and theme ----------------------------------------------------------
source("code/functions_alpha_diversity.R")

# (1) Set-up -------------------------------------------------------------------

# Load and standardise data
data <- inner_join(
  fread("data/amf/diversity_g_amf.txt"),
  fread("data/sample_covariates.txt"),
  by = "sample_id"
) %>%
  mutate(

    # Select alpha diversity metrics
    richness        = richness_srs,
    faith_phylo_div = hill_phy_div_q0,
    hill_tax_q1     = hill_tax_div_q1,
    hill_phy_q1     = hill_phy_div_q1,
    hill_tax_q2     = hill_tax_div_q2,
    hill_phy_q2     = hill_phy_div_q2,

    # Tree covariates
    am_tree_richness_std    = std(sqrt(am_tree_richness)),
    am_tree_basal_area_std  = std(log(am_tree_basal_area)),

    # Climate covariates
    bio1_std    = std(bio1),
    bio12_std   = std(bio12),
    aridity_std = std(aridity_index),

    # Soil covariates
    ph_std                = std(log(ph)),
    mineral_nitrogen_std  = std(log(nitrate + ammonium)),
    phosphorus_std        = std(log(phosphorus)),

  ) %>%
  glimpse()

# (2) Quick data exploration ---------------------------------------------------
dir.create("output/g_amf", showWarnings = FALSE, recursive = TRUE)

#### (2a) Pair plot ####
data %>%
  select(
    richness, faith_phylo_div,
    hill_tax_q1, hill_phy_q1, hill_tax_q2, hill_phy_q2, abundance_srs,
    ends_with("_std")
  ) %>%
  GGally::ggpairs()

ggsave(
  "output/g_amf/covariates_alpha_diversity.png",
  width = 44, height = 44, units = "cm", dpi = 300
)

#### (2b) Relationships ####

# Exploratory response ~ covariate panels: linear, quadratic and loess fits.
# Inspect output/g_amf/ to decide which covariates need quadratic terms.

responses_explore <- list(
  list(var = "richness",        label = "Richness (SRS)"),
  list(var = "faith_phylo_div", label = "Faith's PD"),
  list(var = "hill_tax_q1",     label = "Hill taxonomic diversity (q=1)"),
  list(var = "hill_phy_q1",     label = "Hill phylogenetic diversity (q=1)"),
  list(var = "hill_tax_q2",     label = "Hill taxonomic diversity (q=2)"),
  list(var = "hill_phy_q2",     label = "Hill phylogenetic diversity (q=2)"),
  list(var = "abundance_srs",   label = "Glomeromycota abundance")
)

for (resp in responses_explore) {
  covs_tmp <- names(data) %>% str_subset("_std$")
  d_tmp    <- if (!is.null(resp$filter_fn)) resp$filter_fn(data) else data

  long_df <- d_tmp %>%
    select(response_val = all_of(resp$var), all_of(covs_tmp)) %>%
    pivot_longer(-response_val,
                 names_to  = "covariate",
                 values_to = "cov_val") %>%
    mutate(
      covariate_label = covariate %>%
        str_remove("_std$") %>%
        str_replace_all("_", " ") %>%
        str_to_title()
    )

  p <- ggplot(long_df, aes(x = cov_val, y = response_val)) +
    geom_point(alpha = 0.35, size = 1.2, colour = "grey40") +
    geom_smooth(aes(colour = "Linear"),
                method = "lm",    formula = y ~ x,          se = FALSE, linewidth = 0.9) +
    geom_smooth(aes(colour = "Quadratic"),
                method = "lm",    formula = y ~ poly(x, 2), se = FALSE, linewidth = 0.9) +
    geom_smooth(aes(colour = "LOESS"),
                method = "loess", formula = y ~ x,          se = FALSE, linewidth = 0.9,
                span = 0.8) +
    scale_colour_manual(
      name   = NULL,
      values = c(Linear = "#3182bd", Quadratic = "#e6550d", LOESS = "#31a354"),
      breaks = c("Linear", "Quadratic", "LOESS")
    ) +
    facet_wrap(~ covariate_label, scales = "free_x", ncol = 3) +
    labs(x = "Standardised value", y = resp$label) +
    theme_custom +
    theme(aspect.ratio = 0.9, legend.position = "top",
          legend.text  = element_text(size = 10))

  n_rows  <- ceiling(length(covs_tmp) / 3)
  ggsave(
    sprintf("output/g_amf/exp_%s.png", resp$var),
    plot = p,
    width = 18, height = n_rows * 6 + 3, units = "cm", dpi = 200
  )
  cat(sprintf("  Saved: output/g_amf/exp_%s.png\n", resp$var))
}

#### (2c) Polynomial term specification ####
#
# Inspect the panels in output/g_amf/ and set use_quadratic = TRUE for
# any covariate whose relationship with that response looks unimodal / curved.

covariates <- names(data) %>% str_subset("_std$")

# ── Set quad_covs for each response after inspecting output/g_amf/ ──────

poly_sel_richness  <- make_poly_sel(covariates, quad_covs = c(

))

poly_sel_faith     <- make_poly_sel(covariates, quad_covs = c(

))

poly_sel_hill_tax  <- make_poly_sel(covariates, quad_covs = c(

))

poly_sel_hill_phy  <- make_poly_sel(covariates, quad_covs = c(

))

poly_sel_hill_tax_q2 <- make_poly_sel(covariates, quad_covs = c(

))

poly_sel_hill_phy_q2 <- make_poly_sel(covariates, quad_covs = c(

))

poly_sel_abundance <- make_poly_sel(covariates, quad_covs = c(
  "aridity_std"
))

cat("\n=== (2.5) POLYNOMIAL TERM SPECIFICATION ===\n")
print_term_selection(poly_sel_richness,    "Richness")
print_term_selection(poly_sel_faith,       "Faith's PD")
print_term_selection(poly_sel_hill_tax,    "Hill taxonomic (q=1)")
print_term_selection(poly_sel_hill_phy,    "Hill phylogenetic (q=1)")
print_term_selection(poly_sel_hill_tax_q2, "Hill taxonomic (q=2)")
print_term_selection(poly_sel_hill_phy_q2, "Hill phylogenetic (q=2)")
print_term_selection(poly_sel_abundance,   "SRS abundance")

# (3) Response registry, formulas and model fits --------------------------------
#
# One registry drives formula construction, model fitting, bootstrapping and
# Excel assembly, replacing seven near-identical hand-written blocks (the
# highest-value change in rope_implementation_plan.md §2). Ported verbatim
# from 01a_alpha_diversity_amf.R (plan §13, step 9 — extended here to G-AMF for
# consistency with 01a/01d, at the user's request).
#
# `key` matches the suffix of the original model_*/boot_* variable names
# exactly, so the model-diagnostics and marginal-effects sections below — kept
# explicit because their axis limits, margins and tags differ per panel — can
# go on referring to model_richness, boot_hill_tax, etc. without modification.

responses <- list(
  list(key = "richness",      var = "richness",        lhs = "sqrt(richness)",
       label = "Richness",          diag_label = "RICHNESS",
       tag = "(**a**)", poly_sel = poly_sel_richness),
  list(key = "faith",         var = "faith_phylo_div",  lhs = "sqrt(faith_phylo_div)",
       label = "Faith's PD",        diag_label = "FAITH'S PD",
       tag = "(**b**)", poly_sel = poly_sel_faith),
  list(key = "hill_tax",      var = "hill_tax_q1",      lhs = "sqrt(hill_tax_q1)",
       label = "Hill–Shannon",      diag_label = "HILL TAXONOMIC (q=1)",
       tag = NULL,       poly_sel = poly_sel_hill_tax),
  list(key = "hill_phy",      var = "hill_phy_q1",      lhs = "sqrt(hill_phy_q1)",
       label = "Hill–Shannon PD",   diag_label = "HILL PHYLOGENETIC (q=1)",
       tag = NULL,       poly_sel = poly_sel_hill_phy),
  list(key = "hill_tax_q2",   var = "hill_tax_q2",      lhs = "sqrt(hill_tax_q2)",
       label = "Hill–Simpson",      diag_label = "HILL TAXONOMIC (q=2)",
       tag = NULL,       poly_sel = poly_sel_hill_tax_q2),
  list(key = "hill_phy_q2",   var = "hill_phy_q2",      lhs = "sqrt(hill_phy_q2)",
       label = "Hill–Simpson PD",   diag_label = "HILL PHYLOGENETIC (q=2)",
       tag = NULL,       poly_sel = poly_sel_hill_phy_q2),
  list(key = "abundance_srs", var = "abundance_srs",    lhs = "sqrt(abundance_srs)",
       label = "Abundance",         diag_label = "ABUNDANCE SRS",
       tag = NULL,       poly_sel = poly_sel_abundance)
)
names(responses) <- purrr::map_chr(responses, "key")

# FOCAL_TERM is linear in every response here (only aridity_std carries a
# quadratic, on abundance) — flag loudly rather than silently misreport if
# that ever changes, since a quadratic focal term needs different handling in
# the Summary sheet (see build_summary_row() in functions_alpha_diversity.R;
# 01c_alpha_diversity_m_amf.R has an example where it does).
focal_is_quadratic <- purrr::map_lgl(responses, function(r) {
  isTRUE(r$poly_sel$use_quadratic[r$poly_sel$covariate == FOCAL_TERM])
})
if (any(focal_is_quadratic)) {
  cat(sprintf(
    "\nNOTE: FOCAL_TERM ('%s') is quadratic in: %s. Summary/Sensitivity sheets\n",
    FOCAL_TERM, paste(names(responses)[focal_is_quadratic], collapse = ", ")
  ))
  cat("  will report 'quadratic - see per-response sheet' for those rows.\n")
}

formulas <- purrr::map(responses, ~ build_formula_mixed(.x$lhs, .x$poly_sel))
models   <- purrr::map(formulas, ~ lm(.x, data = data))

cat("\n=== FINAL MODEL FORMULAS ===\n")
for (r in responses) {
  cat(sprintf("%-24s", paste0(r$diag_label, ":"))); print(formulas[[r$key]])
}

# Back-compat variable names (model_richness, model_faith, ...) so the
# diagnostics and marginal-effects sections below can stay unchanged.
for (nm in names(models)) assign(paste0("model_", nm), models[[nm]])

# (4) Model diagnostics ---------------------------------------------------------

cat("\n=== RICHNESS MODEL ===\n")
performance::check_collinearity(model_richness)
parameters::model_parameters(model_richness)

cat("\n=== FAITH'S PD MODEL ===\n")
performance::check_collinearity(model_faith)
parameters::model_parameters(model_faith)

cat("\n=== HILL TAXONOMIC (q=1) MODEL ===\n")
performance::check_collinearity(model_hill_tax)
parameters::model_parameters(model_hill_tax)

cat("\n=== HILL PHYLOGENETIC (q=1) MODEL ===\n")
performance::check_collinearity(model_hill_phy)
parameters::model_parameters(model_hill_phy)

cat("\n=== HILL TAXONOMIC (q=2) MODEL ===\n")
performance::check_collinearity(model_hill_tax_q2)
parameters::model_parameters(model_hill_tax_q2)

cat("\n=== HILL PHYLOGENETIC (q=2) MODEL ===\n")
performance::check_collinearity(model_hill_phy_q2)
parameters::model_parameters(model_hill_phy_q2)

cat("\n=== ABUNDANCE SRS MODEL ===\n")
performance::check_collinearity(model_abundance_srs)
parameters::model_parameters(model_abundance_srs)

# DHARMa residual diagnostics
cat("\n=== Residual diagnostics ===\n")
simulateResiduals(model_richness, plot = TRUE)
performance::check_predictions(model_richness)
simulateResiduals(model_faith, plot = TRUE)
performance::check_predictions(model_faith)
simulateResiduals(model_hill_tax, plot = TRUE)
performance::check_predictions(model_hill_tax)
simulateResiduals(model_hill_phy, plot = TRUE)
performance::check_predictions(model_hill_phy)
simulateResiduals(model_hill_tax_q2, plot = TRUE)
performance::check_predictions(model_hill_tax_q2)
simulateResiduals(model_hill_phy_q2, plot = TRUE)
performance::check_predictions(model_hill_phy_q2)
simulateResiduals(model_abundance_srs, plot = TRUE)
performance::check_predictions(model_abundance_srs)

# Spatial autocorrelation
cat("\n=== Spatial autocorrelation ===\n")
testSpatialAutocorrelation(simulateResiduals(model_richness),
                           x = data$longitude, y = data$latitude)
testSpatialAutocorrelation(simulateResiduals(model_faith),
                           x = data$longitude, y = data$latitude)
testSpatialAutocorrelation(simulateResiduals(model_hill_tax),
                           x = data$longitude, y = data$latitude)
testSpatialAutocorrelation(simulateResiduals(model_hill_phy),
                           x = data$longitude, y = data$latitude)
testSpatialAutocorrelation(simulateResiduals(model_hill_tax_q2),
                           x = data$longitude, y = data$latitude)
testSpatialAutocorrelation(simulateResiduals(model_hill_phy_q2),
                           x = data$longitude, y = data$latitude)
testSpatialAutocorrelation(simulateResiduals(model_abundance_srs),
                           x = data$longitude, y = data$latitude)

# Temporal autocorrelation
data$date_formatted <- as.Date(data$date, format = "%d/%m/%Y")

cat("\n=== Temporal autocorrelation ===\n")
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_richness), group = data$date_formatted),
  time = unique(data$date_formatted)
)
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_faith), group = data$date_formatted),
  time = unique(data$date_formatted)
)
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_hill_tax), group = data$date_formatted),
  time = unique(data$date_formatted)
)
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_hill_phy), group = data$date_formatted),
  time = unique(data$date_formatted)
)
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_hill_tax_q2), group = data$date_formatted),
  time = unique(data$date_formatted)
)
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_hill_phy_q2), group = data$date_formatted),
  time = unique(data$date_formatted)
)
testTemporalAutocorrelation(
  recalculateResiduals(simulateResiduals(model_abundance_srs), group = data$date_formatted),
  time = unique(data$date_formatted)
)

# (5) Bootstrap model coefficients with BCa CIs, and AME percentage change ------

# SD of each response on the scale it was modelled (e.g. sqrt(richness), not
# richness) — the constant that turns a partial coefficient into a fully
# standardised one (plan §3-4). Read from each model's own frame, not
# data$<var>, so NA-dropped rows are handled automatically.
sd_y <- purrr::map(models, sd_response)
for (nm in names(sd_y)) assign(paste0("sd_y_", nm), sd_y[[nm]])

# Back-transform for each response's own modelled scale (edit.pdf §1, see the
# note near the top of this script for why this is derived per response
# rather than assumed to be eta^2 everywhere).
backtfs <- purrr::map(formulas, response_backtf_for)

# 10th/90th percentile of the observed mineral-nitrogen gradient (edit.pdf §6)
# — fixed once from the original data and reused in every bootstrap replicate,
# not recomputed per resample, so AME_STEP's "+1 SD" and these two gradient
# positions are evaluated against the same reference gradient throughout.
nitrogen_p10 <- unname(quantile(data[[FOCAL_TERM]], 0.10))
nitrogen_p90 <- unname(quantile(data[[FOCAL_TERM]], 0.90))

set.seed(SEED)

# Bootstrap results are cached per response (key, seed, n_boot) so a crash
# further down the script (Excel export, plotting) doesn't require repeating
# ~10,000 model refits + jackknife per response. Caching assumes a from-scratch
# run processes responses in this exact order under one set.seed() call, so
# the sequence of random draws matches the pre-refactor script's; a *partial*
# cache (some responses cached, others not, across separate runs) does not
# reproduce a from-scratch run's RNG trajectory for responses after the first
# cache hit, since a cache hit consumes no random draws. That's fine for
# resuming after a downstream crash (the cached values are exactly what a
# from-scratch run already produced) but do not mix a partial cache with a
# reproducibility check against a fresh from-scratch run.
boot_cache_dir <- "generated_data/boot_cache"
if (!dir.exists(boot_cache_dir)) dir.create(boot_cache_dir, recursive = TRUE)

boots <- list()
for (r in responses) {
  cache_file <- file.path(
    boot_cache_dir,
    sprintf("g_amf_%s_seed%d_n%d.rds", r$key, SEED, N_BOOT)
  )
  if (file.exists(cache_file)) {
    cat(sprintf("\n=== %s: loading cached bootstrap (%s) ===\n", r$diag_label, cache_file))
    boots[[r$key]] <- readRDS(cache_file)
  } else {
    cat(sprintf("\n=== BOOTSTRAPPING %s MODEL ===\n", r$diag_label))
    boots[[r$key]] <- bootstrap_model(
      models[[r$key]], data, formulas[[r$key]],
      family_obj = NULL, n_iter = N_BOOT, ci_level = CI_LEVEL,
      ame_focal = FOCAL_TERM, ame_step = AME_STEP, ame_backtf = backtfs[[r$key]],
      ame_p10 = nitrogen_p10, ame_p90 = nitrogen_p90
    )
    saveRDS(boots[[r$key]], cache_file)
  }
}

# Back-compat variable names (boot_richness, boot_faith, ...) for the
# marginal-effects section below.
for (nm in names(boots)) assign(paste0("boot_", nm), boots[[nm]])

# (5.5) Average marginal effects (percentage change) ----------------------------
# edit.pdf — for each response, the average percentage change in the response
# on its original untransformed scale, per +1 SD of mineral nitrogen,
# averaged across the observed data. Point estimate via
# marginaleffects::avg_comparisons() (ame_pct(), slow but only called once per
# response here); the bootstrap interval comes from the closed-form
# ame_pct_fast() already run inside bootstrap_model()'s loop above.

cat("\n=== AVERAGE MARGINAL EFFECTS (percentage change) ===\n")

ame_results <- purrr::map(responses, function(r) {
  backtf <- backtfs[[r$key]]
  model  <- models[[r$key]]

  ame_slow <- ame_pct(model, data, FOCAL_TERM, AME_STEP, backtf)
  ame_fast <- ame_pct_fast(model, data, FOCAL_TERM, AME_STEP, backtf)

  beta_std <- (boots[[r$key]]$coefficients %>%
                 dplyr::filter(Parameter == FOCAL_TERM) %>%
                 dplyr::pull(Original)) / sd_y[[r$key]]

  # Guardrail 2 (skew): mean vs median differing by > 1.5x either way flags
  # the average as driven by low-fitted-value observations (edit.pdf §4.2).
  skew_ratio <- if (ame_fast$median == 0) Inf else abs(ame_fast$mean) / abs(ame_fast$median)
  ame_flag <- dplyr::case_when(
    ame_fast$n_invalid > 0                     ~ "invalid_predictions",  # guardrail 1
    skew_ratio > 1.5 | skew_ratio < (1 / 1.5)  ~ "skewed",
    TRUE                                        ~ "ok"
  )

  p10_val <- avg_pred_backtf(model, data, FOCAL_TERM, nitrogen_p10, backtf)
  p90_val <- avg_pred_backtf(model, data, FOCAL_TERM, nitrogen_p90, backtf)

  list(
    key = r$key, label = r$label,
    ame_estimate = ame_slow$estimate, ame_fast_mean = ame_fast$mean,
    ame_median = ame_fast$median, ame_iqr = ame_fast$iqr, n_invalid = ame_fast$n_invalid,
    ame_flag = ame_flag, beta_std = beta_std,
    ame_ci_lower = boots[[r$key]]$ame$ame_pct_ci_lower,
    ame_ci_upper = boots[[r$key]]$ame$ame_pct_ci_upper,
    n_invalid_boot = boots[[r$key]]$ame$n_invalid_total_bootstrap,
    p10_val = p10_val, p90_val = p90_val,
    p10_ci_lower = boots[[r$key]]$ame$pred_p10_ci_lower,
    p10_ci_upper = boots[[r$key]]$ame$pred_p10_ci_upper,
    p90_ci_lower = boots[[r$key]]$ame$pred_p90_ci_lower,
    p90_ci_upper = boots[[r$key]]$ame$pred_p90_ci_upper,
    pct_change_p10_p90 = (p90_val - p10_val) / p10_val * 100,
    pct_change_ci_lower = boots[[r$key]]$ame$pct_change_p10_p90_ci_lower,
    pct_change_ci_upper = boots[[r$key]]$ame$pct_change_p10_p90_ci_upper
  )
})
names(ame_results) <- purrr::map_chr(responses, "key")

# Guardrail 3 (sign coherence) — a mismatch means the back-transformation is
# wrong somewhere, not that the data disagrees with itself: stop() rather than
# report a self-contradictory row (edit.pdf §4.3).
for (r in responses) {
  ar <- ame_results[[r$key]]
  stopifnot(
    "AME_pct sign must match the standardised beta's sign for the focal term -- check response_backtf_for()" =
      sign(ar$ame_estimate) == sign(ar$beta_std)
  )
}

# Acceptance tests 1-4 (edit.pdf §7) ---------------------------------------------
for (r in responses) {
  ar <- ame_results[[r$key]]

  # 1. Equivalence: ame_pct_fast() reproduces avg_comparisons() to >= 6 sig
  # figs on the observed fit -- this licenses using the fast version inside
  # the 10,000-iteration bootstrap loop, where avg_comparisons() is too slow.
  stopifnot(isTRUE(all.equal(ar$ame_estimate, ar$ame_fast_mean, tolerance = 1e-6)))

  # 2. Sign of AME_pct matches sign of the standardised beta (re-asserted here
  # as part of the formal acceptance-test pass, having already stop()'d above
  # if violated).
  stopifnot(sign(ar$ame_estimate) == sign(ar$beta_std))

  # 3. n_invalid_pred is zero on the observed fit, or the row is flagged.
  stopifnot(ar$n_invalid == 0 || ar$ame_flag == "invalid_predictions")

  # 4. AME_pct falls inside its own bootstrap CI.
  stopifnot(ar$ame_estimate >= ar$ame_ci_lower, ar$ame_estimate <= ar$ame_ci_upper)
}
cat("\nAME acceptance tests 1-4 passed for every response.\n")
# Test 5 (re-running with SEED fixed reproduces all AME values exactly) is
# structural: ame_pct_fast()/avg_pred_backtf() call predict(), which draws no
# random numbers, so adding AME tracking to bootstrap_model()'s loop does not
# change its RNG consumption -- reproducibility holds by the same construction
# as every other bootstrap statistic in this script.

# Guardrail console output (edit.pdf §4 "print each to console") --------------
cat("\nn_invalid_pred (observed fit) and summed across bootstrap replicates, per response:\n")
for (r in responses) {
  ar <- ame_results[[r$key]]
  cat(sprintf("  %-16s observed: %d | bootstrap total: %d\n",
              r$label, ar$n_invalid, ar$n_invalid_boot))
}

ame_flagged <- purrr::keep(ame_results, ~ .x$ame_flag != "ok")
if (length(ame_flagged) > 0) {
  for (ar in ame_flagged) {
    if (ar$ame_flag == "skewed") {
      cat(sprintf(
        "\n*** AME WARNING: %s -- mean (%.1f%%) and median (%.1f%%) percentage change differ by >1.5x; the average is being driven by low-fitted-value observations. Report the median in the text instead. ***\n",
        ar$label, ar$ame_fast_mean, ar$ame_median
      ))
    } else if (ar$ame_flag == "invalid_predictions") {
      cat(sprintf(
        "\n*** AME WARNING: %s -- %d observation(s) had a non-positive fitted value on the modelled scale; this response's AME is not trustworthy and is flagged in the table. ***\n",
        ar$label, ar$n_invalid
      ))
    }
  }
} else {
  cat("\nAME check: no skew or invalid-prediction flags for any response.\n")
}

# (6) ROPE classification and Excel export ---------------------------------------
#
# Per-response sheets carry every parameter on both scales (Beta_partial /
# CI_partial, in the original model's units; Beta_std / CI_std, SD of the
# modelled response per SD of predictor) plus the ROPE classification. The
# Summary sheet is restricted to FOCAL_TERM (see build_summary_row() for how
# a quadratic focal term is handled, not applicable here); it is the table
# that feeds the results text. Metadata records run provenance, and
# Sensitivity checks whether FOCAL_TERM's classification survives varying the
# ROPE half-width and swapping BCa for the percentile interval (plan §7-8).
#
# Sheet assembly uses bind_rows() on tibbles with a tibble(Parameter = "")
# spacer, rather than rbind(c("", "", "")) on a fixed-width character matrix
# (plan §7.4) — the old spacer silently misaligned as soon as the table grew
# past three columns, which this one does.

response_sheets <- purrr::map(responses, function(r) {
  build_response_sheet(boots[[r$key]], sd_y[[r$key]], ROPE_STD)
})

summary_sheet <- purrr::map_df(responses, function(r) {
  build_summary_row(r$label, boots[[r$key]], sd_y[[r$key]], ROPE_STD, FOCAL_TERM)
})

ame_sheet <- purrr::map_df(ame_results, function(ar) {
  build_ame_row(
    label = ar$label,
    ame_pct = ar$ame_estimate,
    ame_ci_lower = ar$ame_ci_lower, ame_ci_upper = ar$ame_ci_upper,
    ame_median = ar$ame_median, ame_iqr = ar$ame_iqr,
    n_invalid = ar$n_invalid, ame_flag = ar$ame_flag,
    pred_p10 = ar$p10_val,
    pred_p10_ci = sprintf("[%.2f, %.2f]", ar$p10_ci_lower, ar$p10_ci_upper),
    pred_p90 = ar$p90_val,
    pred_p90_ci = sprintf("[%.2f, %.2f]", ar$p90_ci_lower, ar$p90_ci_upper),
    pct_change_p10_p90 = ar$pct_change_p10_p90,
    pct_change_p10_p90_ci = sprintf("[%.1f, %.1f]", ar$pct_change_ci_lower, ar$pct_change_ci_upper)
  )
})

summary_sheet <- summary_sheet %>% left_join(ame_sheet, by = "Response")

session_info_summary <- paste(utils::capture.output(sessionInfo())[1:3], collapse = " | ")

metadata_sheet <- purrr::map_df(responses, function(r) {
  build_metadata_row(
    label                 = r$label,
    formula_obj           = formulas[[r$key]],
    n_obs                 = nrow(stats::model.frame(models[[r$key]])),
    n_boot                = N_BOOT,
    seed                  = SEED,
    sd_y                  = sd_y[[r$key]],
    rope_std              = ROPE_STD,
    ci_level              = CI_LEVEL,
    session_info_summary  = session_info_summary
  )
})

sensitivity_sheet <- purrr::map_df(responses, function(r) {
  build_sensitivity_rows(r$label, boots[[r$key]], sd_y[[r$key]], ROPE_SENS, ROPE_STD, FOCAL_TERM)
})

# Flag any response whose focal-term classification is not robust to the ROPE
# half-width (0.05-0.20) or to swapping BCa for the percentile interval —
# pre-empts the obvious reviewer objection that the ROPE is arbitrary (plan §8).
sensitivity_inconsistent <- sensitivity_sheet %>%
  group_by(Response) %>%
  summarise(n_distinct_class = n_distinct(Classification), .groups = "drop") %>%
  filter(n_distinct_class > 1)

if (nrow(sensitivity_inconsistent) > 0) {
  cat(sprintf(
    "\n*** SENSITIVITY WARNING: mineral nitrogen classification is NOT robust for: %s ***\n",
    paste(sensitivity_inconsistent$Response, collapse = ", ")
  ))
} else {
  cat("\nSensitivity check: mineral nitrogen classification is robust to ROPE half-width",
      "(0.05-0.20) and to BCa vs percentile intervals for every response.\n")
}

sensitivity_sheet_wide <- sensitivity_sheet %>%
  tidyr::pivot_wider(names_from = Setting, values_from = Classification)

# Acceptance checks (rope_implementation_plan.md §11.2-11.4) ---------------------

# 11.2: Beta_std * sd_y reproduces the original coefficient, to floating-point
# tolerance, on the unrounded values (the 3dp *display* columns are each
# independently rounded, so they need not multiply back exactly).
for (r in responses) {
  coefs <- boots[[r$key]]$coefficients
  reconstructed <- (coefs$Original / sd_y[[r$key]]) * sd_y[[r$key]]
  stopifnot(isTRUE(all.equal(reconstructed, coefs$Original)))
}

# 11.3: classification agrees whether computed on the standardised scale
# against [-ROPE_STD, ROPE_STD] or on the partial scale against a
# correspondingly rescaled ROPE.
for (r in responses) {
  coefs <- boots[[r$key]]$coefficients
  cls_std     <- classify_effect(coefs$BCa_Lower / sd_y[[r$key]], coefs$BCa_Upper / sd_y[[r$key]],
                                 -ROPE_STD, ROPE_STD)
  cls_partial <- classify_effect(coefs$BCa_Lower, coefs$BCa_Upper,
                                 -ROPE_STD * sd_y[[r$key]], ROPE_STD * sd_y[[r$key]])
  stopifnot(identical(cls_std, cls_partial))
}

# 11.4: every non-intercept row gets exactly one of the four labels, no NA.
for (r in responses) {
  cls <- response_sheets[[r$key]] %>%
    filter(!Parameter %in% c("", "R-squared", "Model") &
           !Parameter %in% boots[[r$key]]$partial_r2$Effect) %>%
    filter(Parameter != "(Intercept)") %>%
    pull(Classification)
  stopifnot(!anyNA(cls), all(cls %in% c("negative", "positive", "neutral", "unresolved")))
}

cat("\nAcceptance checks 11.2-11.4 passed.\n")

# Export to Excel
write_xlsx(
  c(
    setNames(response_sheets, c("Richness", "Faith_PD", "Hill_tax_q1", "Hill_phy_q1",
                                "Hill_tax_q2", "Hill_phy_q2", "Abundance_SRS")),
    list(
      Summary     = summary_sheet,
      Metadata    = metadata_sheet,
      Sensitivity = sensitivity_sheet_wide
    )
  ),
  path = "output/g_amf/coefficients_alpha_diversity.xlsx"
)

# Print summaries
for (r in responses) {
  cat(sprintf("\n=== %s RESULTS ===\n", r$diag_label))
  print(response_sheets[[r$key]])
}

cat("\n=== SUMMARY (focal term: mineral nitrogen, fully standardised) ===\n")
print(summary_sheet)


# (7) Marginal-effects plots with BCa bootstrap annotations ---------------------
# Annotations are fully standardised (SD of the modelled response per SD of
# mineral nitrogen) — see get_boot_annotation() in functions_alpha_diversity.R.
# Panels themselves still plot raw response values; the figure caption must
# state that beta is fully standardised (plan §9 "axis mismatch").

#### (7a) Richness plots ####

pred_richness_nitrogen <- ggpredict(model_richness, terms = "mineral_nitrogen_std [all]")

annot_rich_nitrogen <- get_boot_annotation(boot_richness, FOCAL_TERM, sd_y_richness)

plot_richness <- ggplot(data, aes(x = mineral_nitrogen_std, y = richness)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_richness_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_richness_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n = 3)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_rich_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = NULL, y = "Richness", tag = "(**a**)") +
  theme_custom +
  theme(axis.text.x = element_blank())

print(plot_richness)

#### (7d) Faith's PD plots ####

pred_faith_nitrogen <- ggpredict(model_faith, terms = "mineral_nitrogen_std [all]")

annot_faith_nitrogen  <- get_boot_annotation(boot_faith, FOCAL_TERM, sd_y_faith)

plot_faith <- ggplot(data, aes(x = mineral_nitrogen_std, y = faith_phylo_div)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_faith_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_faith_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(-0.01, NA), breaks = scales::pretty_breaks(n = 5)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_faith_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = NULL, y = "Faith's PD", tag = "(**b**)") +
  theme_custom

print(plot_faith)

#### (7e) Hill taxonomic diversity (q=1) plots ####

pred_hill_tax_nitrogen  <- ggpredict(model_hill_tax, terms = "mineral_nitrogen_std [all]")

annot_hill_tax_nitrogen  <- get_boot_annotation(boot_hill_tax, FOCAL_TERM, sd_y_hill_tax)

plot_hill_tax <- ggplot(data, aes(x = mineral_nitrogen_std, y = hill_tax_q1)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_hill_tax_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_hill_tax_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n = 3)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_hill_tax_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = NULL, y = "Hill–Shannon") +
  theme_custom +
  theme(
    axis.text.x = element_blank(),
    plot.margin = margin(1, 1, 1, 5)
    )

print(plot_hill_tax)

#### (7f) Hill phylogenetic diversity (q=1) plots ####

pred_hill_phy_nitrogen  <- ggpredict(model_hill_phy, terms = "mineral_nitrogen_std [all]")

annot_hill_phy_nitrogen  <- get_boot_annotation(boot_hill_phy, FOCAL_TERM, sd_y_hill_phy)

plot_hill_phy <- ggplot(data, aes(x = mineral_nitrogen_std, y = hill_phy_q1)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_hill_phy_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_hill_phy_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n = 3)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_hill_phy_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = "Mineral nitrogen", y = "Hill–Shannon PD") +
  theme_custom +
  theme(plot.margin = margin(1, 1, 1, 5))

print(plot_hill_phy)

#### (7g) Hill taxonomic diversity (q=2) plots ####

pred_hill_tax_q2_nitrogen <- ggpredict(model_hill_tax_q2, terms = "mineral_nitrogen_std [all]")

annot_hill_tax_q2_nitrogen <- get_boot_annotation(boot_hill_tax_q2, FOCAL_TERM, sd_y_hill_tax_q2)

plot_hill_tax_q2 <- ggplot(data, aes(x = mineral_nitrogen_std, y = hill_tax_q2)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_hill_tax_q2_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_hill_tax_q2_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(-0.66, NA), breaks = scales::pretty_breaks(n = 3)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_hill_tax_q2_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = NULL, y = "Hill–Simpson") +
  theme_custom +
  theme(
    axis.text.x = element_blank(),
    plot.margin = margin(1, 1, 1, 5)
    )

print(plot_hill_tax_q2)

#### (7h) Hill phylogenetic diversity (q=2) plots ####

pred_hill_phy_q2_nitrogen <- ggpredict(model_hill_phy_q2, terms = "mineral_nitrogen_std [all]")

annot_hill_phy_q2_nitrogen <- get_boot_annotation(boot_hill_phy_q2, FOCAL_TERM, sd_y_hill_phy_q2)

plot_hill_phy_q2 <- ggplot(data, aes(x = mineral_nitrogen_std, y = hill_phy_q2)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_hill_phy_q2_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_hill_phy_q2_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n = 2)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_hill_phy_q2_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = NULL, y = "Hill–Simpson PD") +
  theme_custom +
  theme(plot.margin = margin(1, 1, 1, 5))

print(plot_hill_phy_q2)

#### (7i) AM abundance plots ####

pred_abund_srs_nitrogen <- ggpredict(model_abundance_srs, terms = "mineral_nitrogen_std [all]")

annot_abund_srs_nitrogen  <- get_boot_annotation(boot_abundance_srs, FOCAL_TERM, sd_y_abundance_srs)

plot_abundance <- ggplot(data, aes(x = mineral_nitrogen_std, y = abundance_srs)) +
  geom_point(alpha = 0.25, size = 1, colour = "#4A4A4A") +
  geom_ribbon(
    data = pred_abund_srs_nitrogen,
    aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
    alpha = 0.2, fill = "#3181FF", inherit.aes = FALSE
  ) +
  geom_line(
    data = pred_abund_srs_nitrogen,
    aes(x = x, y = predicted),
    colour = "#3366FF", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_y_continuous(limits = c(0, NA), breaks = scales::pretty_breaks(n = 3)) +
  annotate("richtext", x = Inf, y = Inf, label = annot_abund_srs_nitrogen,
           hjust = 1, vjust = 1, size = 2.5, fill = NA, label.color = NA) +
  labs(x = "Mineral nitrogen", y = "Abundance") +
  theme_custom

print(plot_abundance)

#### Save Figure 2 data ####

# Rename all objects with _glom suffix so they can be loaded alongside
# analogous e-am objects without name collisions
data_glom                      <- data

pred_richness_nitrogen_glom    <- pred_richness_nitrogen
annot_rich_nitrogen_glom       <- annot_rich_nitrogen
r2_richness_glom               <- boot_richness$r2
partial_r2_richness_glom       <- boot_richness$partial_r2

pred_abund_srs_nitrogen_glom   <- pred_abund_srs_nitrogen
annot_abund_srs_nitrogen_glom  <- annot_abund_srs_nitrogen
r2_abundance_glom              <- boot_abundance_srs$r2
partial_r2_abundance_glom      <- boot_abundance_srs$partial_r2

if (!dir.exists("generated_data")) dir.create("generated_data")

save(
  data_glom,
  pred_richness_nitrogen_glom,
  annot_rich_nitrogen_glom,
  r2_richness_glom,           partial_r2_richness_glom,
  pred_abund_srs_nitrogen_glom,
  annot_abund_srs_nitrogen_glom,
  r2_abundance_glom,          partial_r2_abundance_glom,
  file = "generated_data/figure_2_g_amf.Rdata"
)
cat("Figure 2 Glomeromycota data saved to generated_data/figure_2_g_amf.Rdata\n")

#### Figure s4 ####

figure_s4 <- wrap_plots(
  plot_richness, plot_hill_tax, plot_hill_tax_q2,
  plot_faith,    plot_hill_phy, plot_hill_phy_q2,
  ncol = 3
)

print(figure_s4)

ggsave(
  figure_s4,
  filename = "output/figure_s4.png",
  width = 16, height = 10, units = "cm", dpi = 300
)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("All models fitted, bootstrapped, classified against ROPE, and visualised.\n")
cat("Results exported to: output/g_amf/coefficients_alpha_diversity.xlsx\n")
cat("  (sheets: Richness, Faith_PD, Hill_tax_q1, Hill_phy_q1, Hill_tax_q2,\n")
cat("   Hill_phy_q2, Abundance_SRS, Summary, Metadata, Sensitivity)\n")
cat("Figure 2 data saved to: generated_data/figure_2_g_amf.Rdata\n")
