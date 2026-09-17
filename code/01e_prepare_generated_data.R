# ─────────────────────────────────────────────────────────────────────────────
# 01e: Assemble the Figure 2 panel data for the three mycorrhizal guilds
#
# Runs after 01a-01d (the four alpha-diversity modelling scripts).
#
#   G-AMF  Glomeromycota          (01b_alpha_diversity_g_amf.R,  data_glom)
#   M-AMF  Densosporales          (01c_alpha_diversity_m_amf.R,  data_endo)
#   EMF    Ectomycorrhizal fungi  (01d_alpha_diversity_emf.R,   data_emf)
#
# across three responses: Hill–Shannon (taxonomic, q=1), Hill–Shannon PD
# (phylogenetic, q=1) and SRS abundance.
#
# The alpha-diversity scripts only persisted abundance's pred_* object to
# generated_data/figure_2_*.Rdata — hill_tax_q1/hill_phy_q1 were modelled and
# bootstrapped but never saved. Rather than re-run three 10,000-iteration
# bootstraps, this script reconstructs every response uniformly from two
# sources that already exist:
#
#   1. Fitted curves  — the models are plain lm() with hard-coded quadratic
#      specifications, so refitting from the saved `data_*` frames is exact and
#      deterministic. Verified below: the refitted abundance curve reproduces
#      the saved pred_* object to the last bit (see `stopifnot`).
#
#   2. Coefficients and partial R2 — read back from the bootstrap results already
#      written to output/<group>/coefficients_alpha_diversity_*.xlsx. These are
#      the same numbers the saved annot_* / partial_r2_* objects were built from.
#
# Output: generated_data/figure_2_panel.Rdata
# ─────────────────────────────────────────────────────────────────────────────

library(readxl)
library(ggeffects)
library(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# Model specification, mirrored from the alpha-diversity scripts
# ─────────────────────────────────────────────────────────────────────────────

# Reproduced from 01b/01c/01d so this script stands alone.
make_poly_sel <- function(covariates, quad_covs = character(0)) {
  data.frame(
    covariate     = covariates,
    delta_aicc    = NA_real_,
    use_quadratic = covariates %in% quad_covs,
    stringsAsFactors = FALSE
  )
}

build_formula_mixed <- function(response, selection_results) {
  terms <- ifelse(
    selection_results$use_quadratic,
    paste0(selection_results$covariate, " + I(", selection_results$covariate, "^2)"),
    selection_results$covariate
  )
  as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

# Response transforms and quadratic terms exactly as specified in the source
# scripts. hill_tax_q1/hill_phy_q1 are sqrt-transformed in all three guilds
# (unlike richness/Faith's PD, where EMF is the untransformed exception) --
# see the "key = "hill_tax"" / "key = "hill_phy"" entries in 01b/01c/01d.
guilds <- list(
  `G-AMF` = list(
    rdata      = "generated_data/figure_2_g_amf.Rdata",
    data_obj   = "data_glom",
    suffix     = "glom",
    excel      = "output/g_amf/coefficients_alpha_diversity.xlsx",
    responses  = list(
      hill_tax  = list(var = "hill_tax_q1",     lhs = "sqrt(hill_tax_q1)",     quad = character(0), sheet = "Hill_tax_q1"),
      hill_phy  = list(var = "hill_phy_q1",     lhs = "sqrt(hill_phy_q1)",     quad = character(0), sheet = "Hill_phy_q1"),
      abundance = list(var = "abundance_srs",   lhs = "sqrt(abundance_srs)",   quad = "aridity_std", sheet = "Abundance_SRS")
    )
  ),
  `M-AMF` = list(
    rdata      = "generated_data/figure_2_m_amf.Rdata",
    data_obj   = "data_endo",
    suffix     = "endo",
    excel      = "output/m_amf/coefficients_alpha_diversity.xlsx",
    responses  = list(
      hill_tax  = list(var = "hill_tax_q1",     lhs = "sqrt(hill_tax_q1)",     quad = "ph_std", sheet = "Hill_tax_q1"),
      hill_phy  = list(var = "hill_phy_q1",     lhs = "sqrt(hill_phy_q1)",     quad = c("ph_std", "mineral_nitrogen_std"), sheet = "Hill_phy_q1"),
      abundance = list(var = "abundance_srs",   lhs = "sqrt(abundance_srs)",   quad = "ph_std", sheet = "Abundance_SRS")
    )
  ),
  `EMF` = list(
    rdata      = "generated_data/figure_2_emf.Rdata",
    data_obj   = "data_emf",
    suffix     = "emf",
    excel      = "output/emf/coefficients_alpha_diversity.xlsx",
    responses  = list(
      hill_tax  = list(var = "hill_tax_q1",     lhs = "sqrt(hill_tax_q1)",    quad = character(0), sheet = "Hill_tax_q1"),
      hill_phy  = list(var = "hill_phy_q1",     lhs = "sqrt(hill_phy_q1)",    quad = character(0), sheet = "Hill_phy_q1"),
      abundance = list(var = "abundance_srs",   lhs = "sqrt(abundance_srs)",  quad = "bio1_std",   sheet = "Abundance_SRS")
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# Bootstrap results are read back out of the exported workbooks
# ─────────────────────────────────────────────────────────────────────────────

# Each sheet is three stacked blocks separated by blank rows:
#   coefficients | R-squared | partial R2.
# Split on the blank rows and return the coefficient and partial-R2 blocks.
#
# The sheet carries both the partial coefficient (Beta_partial/CI_partial, in
# the original model's units) and the fully standardised one (Beta_std/CI_std,
# SD of the modelled response per SD of predictor -- rope_implementation_plan.md's
# reporting scale, and what get_boot_annotation() uses live in 01a-01d for their
# own figures). The R-squared/partial-R2 rows only ever populate the
# Beta_partial/CI_partial column position (Beta_std/CI_std is NA there -- R2
# has no "standardised" counterpart), so partial_r2 below still reads
# coef_partial/ci_partial; coefficients keeps both so callers can pick.
read_boot_sheet <- function(file, sheet) {
  d <- suppressMessages(read_excel(file, sheet = sheet)) %>%
    rename(param = Parameter,
           coef_partial = Beta_partial, ci_partial = CI_partial,
           coef_std     = Beta_std,     ci_std     = CI_std) %>%
    mutate(across(everything(), as.character))

  blanks <- which(is.na(d$param))
  stopifnot(length(blanks) >= 2)

  list(
    coefficients = d[seq_len(blanks[1] - 1), ],
    # +2 skips the blank row and the whole-model R2 row
    partial_r2   = d[(blanks[2] + 2):nrow(d), ]
  )
}

# "[-0.17, 0.54]" -> c(lower = -0.17, upper = 0.54)
parse_ci <- function(x) {
  m <- str_match(x, "\\[\\s*(-?[0-9.]+)\\s*,\\s*(-?[0-9.]+)\\s*\\]")
  c(lower = as.numeric(m[, 2]), upper = as.numeric(m[, 3]))
}

# Rebuild the beta annotation in the same format get_boot_annotation() used in
# the alpha-diversity scripts: beta1/beta2 when a quadratic term is present.
# Reads coef_std/ci_std (fully standardised, SD of response per SD of
# predictor) to match get_boot_annotation()'s convention -- coef_partial/
# ci_partial is in the original model's units and must not be shown here.
make_annotation <- function(coefs, param_name = "mineral_nitrogen_std") {
  lin  <- coefs %>% filter(param == param_name)
  quad <- coefs %>% filter(param == paste0("I(", param_name, "^2)"))

  if (nrow(lin) == 0) return("")

  l_ci <- parse_ci(lin$ci_std)
  if (nrow(quad) > 0) {
    q_ci <- parse_ci(quad$ci_std)
    return(sprintf(
      "β<sub>1</sub> = %.2f [%.2f, %.2f]<br>β<sub>2</sub> = %.2f [%.2f, %.2f]",
      as.numeric(lin$coef_std), l_ci[["lower"]], l_ci[["upper"]],
      as.numeric(quad$coef_std), q_ci[["lower"]], q_ci[["upper"]]
    ))
  }
  sprintf("β = %.2f [%.2f, %.2f]",
          as.numeric(lin$coef_std), l_ci[["lower"]], l_ci[["upper"]])
}

# ─────────────────────────────────────────────────────────────────────────────
# Build the long-format panel tables
# ─────────────────────────────────────────────────────────────────────────────

response_levels <- c("Hill–Shannon", "Hill–Shannon PD", "Abundance")
guild_levels    <- names(guilds)

response_labels <- c(hill_tax = "Hill–Shannon", hill_phy = "Hill–Shannon PD", abundance = "Abundance")

panel_points <- list()
panel_pred   <- list()
panel_annot  <- list()
panel_r2     <- list()

for (g in guild_levels) {

  spec <- guilds[[g]]
  env  <- new.env()
  load(spec$rdata, envir = env)
  dat  <- env[[spec$data_obj]]

  # Covariate set is taken from the data in the same way the source scripts do,
  # so any covariate present in the saved frame enters the model.
  covs <- names(dat) %>% str_subset("_std$")

  for (r in names(spec$responses)) {

    rs    <- spec$responses[[r]]
    model <- lm(build_formula_mixed(rs$lhs, make_poly_sel(covs, rs$quad)), data = dat)

    pred <- suppressMessages(
      ggpredict(model, terms = "mineral_nitrogen_std [all]")
    ) %>% as.data.frame()

    # Guard: abundance must reproduce the published curve exactly (its pred_*
    # object is the only one 01a-01d persist to the Rdata file). hill_tax_q1
    # and hill_phy_q1 have no saved counterpart -- like Faith's PD before them
    # -- so they are only checked structurally.
    saved_name <- switch(
      r,
      abundance = paste0("pred_abund_srs_nitrogen_", spec$suffix),
      NULL
    )
    if (!is.null(saved_name)) {
      saved <- as.data.frame(env[[saved_name]])
      stopifnot(nrow(saved) == nrow(pred))
      if (!isTRUE(all.equal(pred$predicted, saved$predicted, tolerance = 0))) {
        stop(sprintf(
          "Refitted %s / %s does not reproduce the saved curve (max diff %.6g)",
          g, r, max(abs(pred$predicted - saved$predicted))
        ))
      }
    }

    boot <- read_boot_sheet(spec$excel, rs$sheet)

    panel_points[[paste(g, r)]] <- tibble(
      guild    = g,
      response = response_labels[[r]],
      x        = dat$mineral_nitrogen_std,
      y        = dat[[rs$var]]
    )

    panel_pred[[paste(g, r)]] <- tibble(
      guild     = g,
      response  = response_labels[[r]],
      x         = pred$x,
      predicted = pred$predicted,
      conf.low  = pred$conf.low,
      conf.high = pred$conf.high
    )

    panel_annot[[paste(g, r)]] <- tibble(
      guild    = g,
      response = response_labels[[r]],
      label    = make_annotation(boot$coefficients)
    )

    # Partial R2 for the heatmap (its value only exists in the
    # coef_partial/ci_partial column position -- R2 has no standardised
    # counterpart); the coefficient block supplies the sign and CI used to
    # derive the direction-of-effect symbol, read on the fully standardised
    # scale for consistency with the rest of the reporting (sign is identical
    # either way, since Beta_std = Beta_partial / sd_y and sd_y > 0).
    panel_r2[[paste(g, r)]] <- boot$partial_r2 %>%
      transmute(
        guild     = g,
        response  = response_labels[[r]],
        effect    = param,
        Mean_R2   = as.numeric(coef_partial),
        CI_Lower  = map_dbl(ci_partial, ~ parse_ci(.x)[["lower"]]),
        CI_Upper  = map_dbl(ci_partial, ~ parse_ci(.x)[["upper"]])
      ) %>%
      left_join(
        boot$coefficients %>%
          transmute(
            effect    = param,
            coef_mean = as.numeric(coef_std),
            coef_lo   = map_dbl(ci_std, ~ parse_ci(.x)[["lower"]]),
            coef_hi   = map_dbl(ci_std, ~ parse_ci(.x)[["upper"]])
          ),
        by = "effect"
      )
  }
  cat("Assembled:", g, "\n")
}

panel_points <- bind_rows(panel_points)
panel_pred   <- bind_rows(panel_pred)
panel_annot  <- bind_rows(panel_annot)
panel_r2     <- bind_rows(panel_r2)

# Factor ordering shared by both sub-panels
for (nm in c("panel_points", "panel_pred", "panel_annot", "panel_r2")) {
  assign(nm, get(nm) %>%
    mutate(
      guild    = factor(guild,    levels = guild_levels),
      response = factor(response, levels = response_levels)
    ))
}

if (!dir.exists("generated_data")) dir.create("generated_data")

save(
  panel_points, panel_pred, panel_annot, panel_r2,
  guild_levels, response_levels,
  file = "generated_data/figure_2_panel.Rdata"
)

cat("\nSaved generated_data/figure_2_panel.Rdata\n")
cat("  points:", nrow(panel_points), "rows | pred:", nrow(panel_pred),
    "rows | r2:", nrow(panel_r2), "rows\n")

