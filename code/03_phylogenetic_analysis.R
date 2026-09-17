library(picante)
library(ape)
library(readxl)
library(writexl)
library(data.table)
library(tidyverse)

# ─────────────────────────────────────────────────────────────────────────────
# Phylogenetic analysis of nitrogen-response indicators (analysis step)
#
# Standardised effect sizes (SES-MPD and SES-MNTD) for the declining (z-) and
# increasing (z+) TITAN2 indicator sets, tested against the candidate-pool null
# (taxa entering that guild's TITAN run) separately within G-AMF, M-AMF and EMF.
# G-AMF and M-AMF are analysed apart, never pooled: they are sister clades
# separated by a single deep branch that would otherwise dominate a combined MPD
# null.
#
# Each effect size is computed twice on the same matched candidate pool: once on
# the pruned phylogram, and once on an integer taxonomic distance (1 = same
# genus, 2 = same lineage/family, 3 = different lineage/family). The tree-free
# run tests whether the effect sizes depend on the topology — relevant because
# EMF ITS branch lengths are unreliable across deep divergences and negative NJ
# branches were zeroed. It is independent of alignment and tree inference, but
# NOT of ITS itself, since the taxonomy is also ITS-derived. Because the
# taxonomic distance takes only three values its null is discrete, so the
# taxonomic SES is descriptive and inference rests on the rank p.
#
# REPORTING: SES is the effect size and the null randomisation characterises
# its uncertainty. We do NOT report SES +/- 1.96 as if it were a confidence
# interval — SES is a z-score locating the observation within the null, so its
# denominator is the null SD, not a standard error. SES +/- 1.96 has constant
# width 3.92 for every set, carries no information beyond the point estimate,
# and excludes zero exactly when |SES| > 1.96, i.e. it is a two-tailed test at
# alpha = 0.05 wearing interval clothing.
#
# ESTIMATION FRAMEWORK (Part B): p_rank and p_norm disagree substantially for
# several sets (see the "Calls differing..." block below), which means the null
# distributions are not Gaussian and |SES| = 1.96 does not mark their 5% tails.
# Classification is made against an empirical ROPE (region of practical
# equivalence) = the central 95% of the null distribution, rather than a fixed
# standardised threshold. See classify_phylo() below and the pre-specification
# block (B1).
#
# ROPE SCALE FIX (ses_rope_scale_fix_spec.md, 2026-08-19): classification and
# pct_outside_rope are computed on the SES scale, with the ROPE evaluated
# separately at EACH bootstrap replicate's own indicator-set size k, not once
# at the observed k on the raw metric scale. Replicate sets typically run
# 2-4x larger than the observed set (see the indicator-set definition mismatch
# note in the BOOT_SES block below), and because MNTD falls sharply as k rises
# while MPD is roughly size-invariant, comparing a single observed-k ROPE
# against replicates of every size produced an artefactual split where every
# MPD row had ~0-10% of replicates outside the ROPE and every MNTD row had
# ~7-93%, regardless of biology. rope_ses_at(guild, basis, metric, k) computes
# the ROPE fresh (cached) at whatever k is asked of it, on the SES scale, so
# replicates of differing size are directly comparable. See B4.
#
# Reported per set:
#   Observed / Null mean / Null 95%  the null band on the observed (MPD or MNTD)
#                                    scale, AT THE OBSERVED k — informational
#                                    (judge magnitude in tree units), not the
#                                    basis for classification since the fix.
#   SES                              standardised effect size, in null SDs
#   p rank                           rank of the observed value in the null; the
#                                    reported inference (continuous, no cutoff)
#   p norm                           pnorm(SES), the normal approximation; carried
#                                    for completeness and because it is what a
#                                    fixed SES +/- threshold would encode
#   Jackknife                        range of SES across leave-one-out sets;
#                                    if it straddles zero, one taxon drives the
#                                    result
#   Bootstrap 95% (raw + SES scale)  membership-uncertainty interval from
#                                    TITAN2's site-resampling replicates (below)
#   ROPE_ses / k_boot_range          the SES-scale ROPE at the observed k, and
#                                    the 2.5th/97.5th percentile of replicate
#                                    set size (so the size spread is visible)
#   pct_outside_rope                 continuous overlap measure (B5): % of
#                                    bootstrap replicates whose SES falls
#                                    outside the ROPE evaluated at THAT
#                                    replicate's own k
#   Classification                   clustered / overdispersed / random /
#                                    unresolved, from classify_phylo() (B4),
#                                    on the SES scale, k-matched
#
# A true sampling CI would require resampling sites, re-running TITAN2 and
# re-deriving the indicator sets — the uncertainty that matters is which taxa
# become indicators. Not attempted from scratch: the nested cost is prohibitive.
# Bootstrapping taxa within a fixed set is worse than useless here, because
# duplicated taxa contribute zero pairwise distance and bias MPD towards
# apparent clustering. Instead we reuse TITAN2's own site-resampling replicates
# (below) as an approximate, membership-uncertainty-only interval.
#
# SES convention: negative = clustered, positive = overdispersed. Mapping to the
# manuscript's four categories (B4): clustered -> negative, overdispersed ->
# positive, random -> neutral, unresolved -> unresolved.
#
# Outputs:
#   output/table_6_ses.xlsx — primary / sensitivity_basis /
#     sensitivity_metric / sensitivity_rope / resolution / metadata /
#     session_info sheets (B6)
#   generated_data/figure_6.Rdata      — plotting data for code/figure_6.R
#
# Indicator sets and change points come from the production TITAN2 runs
# (output/titan_amf.xlsx, output/titan_emf.xlsx; code/02a_titan_amf.R,
# code/02b_titan_emf.R). Earlier signal analyses (Blomberg's K, Fritz & Purvis D)
# are archived in code/drafts/xx_03_phylogenetic_dispersion.R.
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# B1. Pre-specified analysis choices — visible and auditable, not inferred from
# the manuscript.
# ─────────────────────────────────────────────────────────────────────────────
# Guild naming — single lookup, do not scatter string literals (naming migration)
GUILD_LABELS <- c(gamf = "G-AMF", mamf = "M-AMF", emf = "EMF", amf = "AMF")

SEED <- 1986

# Basis labels used throughout this script: "Phylogeny" = the pruned phylogram,
# "Lineage" = the integer taxonomic distance (tree-free). PRIMARY_BASIS is
# reported as the main result; the other basis is reported as a sensitivity
# analysis (sensitivity_basis sheet).
PRIMARY_BASIS <- "Phylogeny"

# MPD for G-AMF/M-AMF (SSU marker: reliable at deeper nodes), MNTD for EMF
# (full ITS: alignment unreliable at depth, so tip-level clustering is the
# defensible signal). State this split explicitly in the methods -- otherwise a
# reviewer will read a guild-specific metric choice as a fishing expedition.
PRIMARY_METRIC <- setNames(c("MPD", "MPD", "MNTD"),
                           c(GUILD_LABELS[["gamf"]], GUILD_LABELS[["mamf"]], GUILD_LABELS[["emf"]]))

ROPE_METHOD  <- "empirical_null_ses_scale_k_matched"  # see B4; "ses_fixed" available as a sensitivity (sensitivity_rope sheet)
ROPE_SES     <- 1.96              # used only for the ses_fixed sensitivity variants below (also see SES_FIXED_THRESHOLDS)
ROPE_NULL_Q  <- c(0.025, 0.975)   # null quantiles defining the empirical ROPE

# Number of TITAN2 site-resampling bootstrap replicates reused from 02a/02b
# (their nBoot argument) -- not a new bootstrap run here, just a fixed
# expectation asserted below so drift between scripts is caught, not silent.
N_BOOT <- 1000

# Sets smaller than this are forced to "unresolved" regardless of point
# estimate (B8): a 3-taxon set has 3 pairwise distances and a degenerate
# jackknife, so a clustered/overdispersed call would be read as a result when
# it is not one.
MIN_TAXA_FOR_INFERENCE <- 5

set.seed(SEED)

SES_RUNS  <- as.integer(Sys.getenv("FIG6_RUNS",  "9999"))  # null draws per statistic
JACK_RUNS <- as.integer(Sys.getenv("FIG6_JACK",  "999"))   # null draws per jackknife fold
VERIFY_PICANTE <- as.logical(Sys.getenv("FIG6_VERIFY", "TRUE"))

# Two-tailed alpha used only to compare the rank and normal-approximation calls.
ALPHA <- 0.05

# Approximate bootstrap CI from TITAN2's site-resampling replicates (below).
# FIG6_BOOTSES=FALSE reproduces the non-bootstrap output unchanged, and skips
# B3-B8 (raw-metric storage, empirical ROPE, classification, sensitivity_rope),
# which all depend on it.
BOOT_SES   <- as.logical(Sys.getenv("FIG6_BOOTSES", "TRUE"))
BOOT_ALPHA <- 0.05   # bootstrap quantiles at BOOT_ALPHA/2 and 1 - BOOT_ALPHA/2

# sensitivity_rope: SES thresholds compared against the empirical-ROPE call (B7)
SES_FIXED_THRESHOLDS <- c(0, 1, 2, 2.5)

# Guild display order (G-AMF and M-AMF analysed separately, then EMF)
guild_levels <- unname(GUILD_LABELS[c("gamf", "mamf", "emf")])

# ─────────────────────────────────────────────────────────────────────────────
# -- Indicator sets and change points from the production TITAN2 runs ---------
# ─────────────────────────────────────────────────────────────────────────────
# change_point   = zenv.cp, the taxon-specific IndVal maximum on the log-nitrogen
#                  gradient (drives the supplementary change-point bars)
# threshold_mgkg = exp(median bootstrap threshold), the back-transformed value
#                  reported in figures 4 and 5
# indicator      = neg (z-, nitrophobic) / pos (z+, nitrophilic) / none
# genus is carried through for the two-rank taxonomic distance (genus in lineage).
build_frame <- function(sppmax, tax) {
  sppmax %>%
    left_join(tax, by = "species") %>%
    transmute(
      otu_id, family_or_lineage, genus, species,
      change_point    = zenv.cp,
      threshold_mgkg  = exp(`50%`),
      indicator       = case_when(filter == 1 ~ "neg", filter == 2 ~ "pos", TRUE ~ "none")
    )
}

sppmax_amf <- read_excel("output/titan_amf.xlsx", sheet = "species")
sppmax_emf <- read_excel("output/titan_emf.xlsx", sheet = "species")

tax_amf <- fread("data/amf/classification.txt") %>%
  select(otu_id, am_group, family, genus, species) %>%
  rename(family_or_lineage = family)

dat_amf_full <- build_frame(sppmax_amf, tax_amf)

# am_group == "Glomeromycota" is G-AMF; "Endogonomycetes" is the Densosporales
# (Mucoromycota) lineage this project treats as M-AMF.
am_group_lookup <- fread("data/amf/classification.txt") %>% select(otu_id, am_group) %>% deframe()
dat_gamf <- dat_amf_full %>% filter(am_group_lookup[otu_id] == "Glomeromycota")
dat_mamf <- dat_amf_full %>% filter(am_group_lookup[otu_id] == "Endogonomycetes")

tax_emf <- fread("data/emf/classification.txt") %>%
  left_join(
    fread("data/emf/otu_table_srs_emf_genus.txt") %>% select(genus, lineage) %>% distinct(),
    by = "genus"
  ) %>%
  select(otu_id, species, genus, lineage) %>%
  rename(family_or_lineage = lineage)

dat_emf <- build_frame(sppmax_emf, tax_emf)

# ─────────────────────────────────────────────────────────────────────────────
# -- Prune each guild's phylogeny to its candidate pool ------------------------
# ─────────────────────────────────────────────────────────────────────────────
# Pool = the taxa entering that guild's TITAN run, i.e. every row of sppmax
# (already restricted to the n_samples >= 5 prevalence filter upstream in
# 02a/02b -- see the "Prevalence filter" row in the metadata sheet, B1/C.3).
# Pruning G-AMF and M-AMF to their own clades drops the deep branch between them.
prepare_tree <- function(tr_full, pool_ids) {
  tr <- keep.tip(tr_full, pool_ids)
  n_neg <- sum(tr$edge.length < 0)
  if (n_neg > 0) tr$edge.length[tr$edge.length < 0] <- 0  # NJ artefact
  if (!is.binary(tr)) tr <- multi2di(tr)
  tr
}

tree_amf_full <- read.tree("data/amf/tree.newick")
tree_emf_full <- read.tree("data/emf/tree.newick")

tree_gamf <- prepare_tree(tree_amf_full, dat_gamf$otu_id)
tree_mamf <- prepare_tree(tree_amf_full, dat_mamf$otu_id)
tree_emf  <- prepare_tree(tree_emf_full, dat_emf$otu_id)

cat("Candidate pools (z- / z+ indicators):\n")
cat(sprintf("  %-5s: %d OTUs (%d / %d)\n", GUILD_LABELS[["gamf"]], nrow(dat_gamf),
            sum(dat_gamf$indicator == "neg"), sum(dat_gamf$indicator == "pos")))
cat(sprintf("  %-5s: %d OTUs (%d / %d)\n", GUILD_LABELS[["mamf"]], nrow(dat_mamf),
            sum(dat_mamf$indicator == "neg"), sum(dat_mamf$indicator == "pos")))
cat(sprintf("  %-5s: %d OTUs (%d / %d)\n", GUILD_LABELS[["emf"]], nrow(dat_emf),
            sum(dat_emf$indicator == "neg"), sum(dat_emf$indicator == "pos")))

# ─────────────────────────────────────────────────────────────────────────────
# -- Matched pool: same taxa for the phylogenetic and taxonomic runs -----------
# ─────────────────────────────────────────────────────────────────────────────
# Taxa without a lineage/family assignment cannot be placed taxonomically. If
# they were retained in the phylogram but dropped from the taxonomic distance, a
# discrepancy between the two analyses could reflect the pool rather than the
# tree. Restrict BOTH to the same taxa and report how many were lost.
matched_pool <- function(dat, tr_full, group_col = "family_or_lineage") {
  keep <- dat %>% filter(!is.na(.data[[group_col]]), .data[[group_col]] != "")
  list(dat     = keep,
       tr      = prepare_tree(tr_full, keep$otu_id),
       n_drop  = nrow(dat) - nrow(keep),
       n_total = nrow(dat))
}

mp_gamf <- matched_pool(dat_gamf, tree_amf_full)
mp_mamf <- matched_pool(dat_mamf, tree_amf_full)
mp_emf  <- matched_pool(dat_emf,  tree_emf_full)

cat("\nMatched pool (taxa dropped for missing lineage/family):\n")
for (x in list(list(GUILD_LABELS[["gamf"]], mp_gamf),
               list(GUILD_LABELS[["mamf"]], mp_mamf),
               list(GUILD_LABELS[["emf"]],  mp_emf))) {
  cat(sprintf("  %-6s %d of %d dropped (%.1f%%)\n", x[[1]],
              x[[2]]$n_drop, x[[2]]$n_total, 100 * x[[2]]$n_drop / x[[2]]$n_total))
}

# ─────────────────────────────────────────────────────────────────────────────
# -- Integer taxonomic distance (tree-free relatedness) ------------------------
# ─────────────────────────────────────────────────────────────────────────────
# 1 = same genus; 2 = same lineage/family, different genus; 3 = different
# lineage/family. Fed to the same SES code path as the phylogram, so the only
# thing that differs between the two analyses is how relatedness is measured.
make_tax_dist <- function(dat, ranks = c("family_or_lineage", "genus")) {
  ids <- dat$otu_id
  d   <- matrix(length(ranks) + 1L, length(ids), length(ids),
                dimnames = list(ids, ids))
  for (k in seq_along(ranks)) {                 # coarse -> fine, finer overwrites
    v    <- dat[[ranks[k]]]
    same <- outer(v, v, "==")
    same[is.na(same)] <- FALSE
    d[same] <- length(ranks) + 1L - k
  }
  diag(d) <- 0
  d
}

# ─────────────────────────────────────────────────────────────────────────────
# -- Diagnostics to report before believing a null result ----------------------
# ─────────────────────────────────────────────────────────────────────────────
# If nearly every genus is a pool singleton, taxonomic MNTD is near-constant and
# the test is degenerate — that is an underpowered test, not evidence of lability.
# Same logic for lineage.
tax_freedom <- function(dat, label) {
  g <- table(dat$genus[!is.na(dat$genus)])
  l <- table(dat$family_or_lineage[!is.na(dat$family_or_lineage)])
  tibble(guild = label,
         n_genera = length(g), genera_ge2 = sum(g >= 2),
         n_lineages = length(l), lineages_ge2 = sum(l >= 2),
         pct_taxa_in_multi_genus = round(100 * sum(g[g >= 2]) / sum(g), 1))
}

freedom <- bind_rows(tax_freedom(mp_gamf$dat, GUILD_LABELS[["gamf"]]),
                     tax_freedom(mp_mamf$dat, GUILD_LABELS[["mamf"]]),
                     tax_freedom(mp_emf$dat,  GUILD_LABELS[["emf"]]))
cat("\nTaxonomic resolution of the pool:\n"); print(as.data.frame(freedom))

# Zeroing negative NJ branches sets some tip-pair cophenetic distances to ~0,
# inflating apparent MNTD clustering. A large count is a further argument for
# leaning on the tree-free result. Report these counts in the Methods.
count_neg_branches <- function(tr_full, ids) {
  tr <- keep.tip(tr_full, ids)
  sum(tr$edge.length < 0)
}
cat(sprintf("\nNegative branches zeroed — %s %d, %s %d, %s %d\n",
            GUILD_LABELS[["gamf"]], count_neg_branches(tree_amf_full, mp_gamf$dat$otu_id),
            GUILD_LABELS[["mamf"]], count_neg_branches(tree_amf_full, mp_mamf$dat$otu_id),
            GUILD_LABELS[["emf"]],  count_neg_branches(tree_emf_full, mp_emf$dat$otu_id)))

# ─────────────────────────────────────────────────────────────────────────────
# -- SES machinery: statistic, null, effect size, rank p, jackknife -------------
# ─────────────────────────────────────────────────────────────────────────────
# Written against a distance matrix rather than a tree so the phylogenetic and
# taxonomic runs share one code path (same statistic, null, pool, runs, seed).
#
# The null is a random draw of k taxa from the candidate pool, which is identical
# to picante's "taxa.labels" model for a single set of fixed size. Drawing
# directly gives access to the null values themselves, which ses.mpd/ses.mntd do
# not return — that is what makes the observed-scale band (and the empirical
# ROPE, B4) possible.

stat_set <- function(dist, ids_set, fun) {
  k <- length(ids_set)
  if (k < 2) return(NA_real_)
  d <- dist[ids_set, ids_set, drop = FALSE]
  if (fun == "mpd") mean(d[lower.tri(d)]) else mean(apply(d + diag(Inf, k), 1, min))
}

null_draws <- function(dist, k, runs, fun) {
  ids <- colnames(dist)
  vapply(seq_len(runs), function(i) stat_set(dist, sample(ids, k), fun), numeric(1))
}

ses_full <- function(dist, ids_set, fun, runs = SES_RUNS, jack_runs = JACK_RUNS) {
  k <- length(ids_set)
  if (k < 2) {
    return(tibble(n = k, obs = NA_real_, null_mean = NA_real_, null_lo = NA_real_,
                  null_hi = NA_real_, ses = NA_real_,
                  p_rank = NA_real_, p_norm = NA_real_,
                  jack_min = NA_real_, jack_max = NA_real_))
  }

  obs <- stat_set(dist, ids_set, fun)
  nd  <- null_draws(dist, k, runs, fun)
  qs  <- unname(quantile(nd, ROPE_NULL_Q))

  # Every leave-one-out fold has size k-1, so one null distribution serves all of
  # them. Fewer runs here than for the main statistic: the jackknife is an
  # influence diagnostic, not a precision estimate.
  jk <- NA_real_
  if (k >= 3 && jack_runs >= 1) {
    nj <- null_draws(dist, k - 1L, jack_runs, fun)
    jk <- vapply(seq_len(k), function(i)
            (stat_set(dist, ids_set[-i], fun) - mean(nj)) / sd(nj), numeric(1))
  }

  # Both p values are lower-tail, so they are directly comparable: small = the
  # observed value sits low in the null (clustered), large = high (overdispersed).
  ses_val <- (obs - mean(nd)) / sd(nd)

  tibble(n = k, obs = obs, null_mean = mean(nd), null_lo = qs[1], null_hi = qs[2],
         ses    = ses_val,
         p_rank = (sum(nd <= obs) + 1) / (runs + 1),  # rank in the null; reported
         p_norm = pnorm(ses_val),                     # normal approx
         jack_min = suppressWarnings(min(jk)), jack_max = suppressWarnings(max(jk)))
}

compute_ses_dist <- function(dist, dat, label, basis) {
  sets <- list("Nitrophobic (z-)" = dat$otu_id[dat$indicator == "neg"],
               "Nitrophilic (z+)" = dat$otu_id[dat$indicator == "pos"])
  expand_grid(set = names(sets), fun = c("mpd", "mntd")) %>%
    mutate(res = map2(set, fun, ~ { set.seed(SEED); ses_full(dist, sets[[.x]], .y) })) %>%
    unnest(res) %>%
    mutate(guild  = label,
           basis  = basis,
           metric = recode(fun, mpd = "SES-MPD", mntd = "SES-MNTD")) %>%
    select(-fun)
}

dist_gamf_phy <- cophenetic(mp_gamf$tr); dist_gamf_tax <- make_tax_dist(mp_gamf$dat)
dist_mamf_phy <- cophenetic(mp_mamf$tr); dist_mamf_tax <- make_tax_dist(mp_mamf$dat)
dist_emf_phy  <- cophenetic(mp_emf$tr);  dist_emf_tax  <- make_tax_dist(mp_emf$dat)

results <- bind_rows(
  # Phylogeny — matched pool, so directly comparable to the taxonomic run
  compute_ses_dist(dist_gamf_phy, mp_gamf$dat, GUILD_LABELS[["gamf"]], "Phylogeny"),
  compute_ses_dist(dist_mamf_phy, mp_mamf$dat, GUILD_LABELS[["mamf"]], "Phylogeny"),
  compute_ses_dist(dist_emf_phy,  mp_emf$dat,  GUILD_LABELS[["emf"]],  "Phylogeny"),
  # Tree-free — identical statistic, null and pool
  compute_ses_dist(dist_gamf_tax, mp_gamf$dat, GUILD_LABELS[["gamf"]], "Lineage"),
  compute_ses_dist(dist_mamf_tax, mp_mamf$dat, GUILD_LABELS[["mamf"]], "Lineage"),
  compute_ses_dist(dist_emf_tax,  mp_emf$dat,  GUILD_LABELS[["emf"]],  "Lineage")
) %>%
  mutate(basis = factor(basis, levels = c("Phylogeny", "Lineage")))

# ─────────────────────────────────────────────────────────────────────────────
# -- B2/B3. Approximate bootstrap CI + raw-metric storage, from TITAN2's -------
# -- site-resampling replicates -------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
# TITAN2's nBoot = 1000 replicates each re-run the IndVal on a bootstrap sample
# of sites. We saved the per-replicate metricArray (02a/02b) and here
# reconstruct each replicate's indicator sets, recompute the raw metric (MPD or
# MNTD) and SES against the same taxa.labels null (cached per set size), and
# take quantiles across replicates. This propagates the one uncertainty a fixed
# SES threshold cannot: WHICH taxa are indicators.
#
# APPROXIMATION, stated in the Methods: purity and reliability are cross-
# replicate summaries, so there is no per-replicate "pure and reliable"
# membership. We use the within-replicate rule TITAN2 itself uses to build
# reliability — direction from maxgrp, IndVal p (obsiv.prob) <= 0.05 — which is
# noisier than the cross-replicate filter, so replicate sets typically run
# larger than the reported (observed) set and the interval errs wide rather
# than narrow. This indicator-set definition mismatch (B2) is why the
# bootstrap interval need not be centred on the observed point estimate; it is
# not evidence of a null-calibration defect (see B2 verification note below).
#
# metricArray was saved subset to two slices: [,1,] = maxgrp (1 = z-, 2 = z+),
# [,2,] = obsiv.prob (the IndVal p-value). Rows are taxa (= `taxa`, TITAN species
# names), mapped to otu_id here so the guild split matches the distance matrices.
if (BOOT_SES) {

  spp2otu <- function(cls) { m <- fread(cls); setNames(m$otu_id, m$species) }
  s2o_amf <- spp2otu("data/amf/classification.txt")
  s2o_emf <- spp2otu("data/emf/classification.txt")

  # Per-replicate neg / pos otu_id sets. The guild split (G-AMF vs M-AMF within
  # the AMF run) happens later, by intersecting each set with the guild's own pool.
  boot_otu_sets <- function(path, s2o, p_cut = 0.05) {
    r <- readRDS(path); MA <- r$metricArray; taxa <- r$taxa
    stopifnot(dim(MA)[3] == N_BOOT)  # B1: catch drift from 02a/02b's nBoot silently
    lapply(seq_len(dim(MA)[3]), function(b) {
      grp <- MA[, 1, b]; p <- MA[, 2, b]
      list(neg = unname(s2o[taxa[which(grp == 1 & p <= p_cut)]]),
           pos = unname(s2o[taxa[which(grp == 2 & p <= p_cut)]]))
    })
  }
  boot_amf <- boot_otu_sets("generated_data/titan_boot_amf.rds", s2o_amf)
  boot_emf <- boot_otu_sets("generated_data/titan_boot_emf.rds", s2o_emf)

  # ─────────────────────────────────────────────────────────────────────────
  # Deterministic, order-independent null cache (ses_rope_scale_fix_spec.md
  # §2.1/§2.2). A null is needed at EVERY k that appears across the bootstrap
  # replicates, not only the observed k -- generated on demand and cached by
  # the full (guild, basis, metric, k) key. The RNG seed for a given key is
  # derived from a hash of the key itself, not from wherever the global stream
  # happens to be when that key is first encountered, so two runs (or the same
  # run after an unrelated code change reorders iteration) always draw the
  # same null for the same key. withr::with_seed() saves/restores the global
  # .Random.seed around each on-demand draw, so this cannot perturb any other
  # random draw in the script (e.g. the main ses_full() calls above, or the
  # jackknife).
  # ─────────────────────────────────────────────────────────────────────────
  dist_lookup <- list(
    "G-AMF" = list(Phylogeny = dist_gamf_phy, Lineage = dist_gamf_tax),
    "M-AMF" = list(Phylogeny = dist_mamf_phy, Lineage = dist_mamf_tax),
    "EMF"   = list(Phylogeny = dist_emf_phy,  Lineage = dist_emf_tax)
  )
  metric_fun_lookup <- c("SES-MPD" = "mpd", "SES-MNTD" = "mntd")

  null_seed <- function(guild, basis, metric, k) {
    SEED + as.integer(strtoi(substr(digest::digest(
      paste(guild, basis, metric, k, sep = "|")), 1, 7), 16L) %% 1e6)
  }

  null_cache <- new.env(parent = emptyenv())      # key -> raw null draws (numeric vector)
  null_requested <- new.env(parent = emptyenv())  # key -> TRUE, every distinct key ever asked for (acceptance test 6)
  get_null <- function(guild, basis, metric, k, runs = SES_RUNS) {
    key <- paste(guild, basis, metric, k, sep = "|")
    assign(key, TRUE, envir = null_requested)
    if (!exists(key, envir = null_cache)) {
      dist <- dist_lookup[[guild]][[basis]]
      fun  <- metric_fun_lookup[[metric]]
      if (is.null(dist) || is.null(fun))
        stop(sprintf("get_null(): no distance matrix / statistic found for key '%s'", key))
      withr::with_seed(null_seed(guild, basis, metric, k),
                        assign(key, null_draws(dist, k, runs, fun), envir = null_cache))
    }
    get(key, envir = null_cache)
  }
  null_stats <- function(guild, basis, metric, k, runs = SES_RUNS) {
    nd <- get_null(guild, basis, metric, k, runs)
    c(mean = mean(nd), sd = sd(nd))
  }

  # B4 fix: the ROPE on the SES scale, evaluated at a given k. Cached per key --
  # called once per bootstrap replicate otherwise, and every replicate sharing
  # a k gets the identical answer. Bounds are expected to be asymmetric where
  # the null is skewed (see acceptance test 3).
  rope_cache <- new.env(parent = emptyenv())
  rope_ses_at <- function(guild, basis, metric, k) {
    key <- paste(guild, basis, metric, k, sep = "|")
    if (!exists(key, envir = rope_cache)) {
      nd   <- get_null(guild, basis, metric, k)
      mu   <- mean(nd); sdev <- sd(nd)
      q    <- stats::quantile(nd, ROPE_NULL_Q, names = FALSE)
      assign(key, (q - mu) / sdev, envir = rope_cache)
    }
    get(key, envir = rope_cache)
  }

  # B3: per-replicate raw metric, null mean, n_boot_taxa and SES -- not just SES.
  # n_boot_taxa is the diagnostic B2.4 asks for (its distribution belongs in the
  # supplement; the per-guild x direction median, and now the 95% range
  # (k_boot_range), are exported below).
  boot_replicates <- function(guild, basis, metric, sets, which, runs = SES_RUNS) {
    dist <- dist_lookup[[guild]][[basis]]
    fun  <- metric_fun_lookup[[metric]]
    out <- lapply(sets, function(s) {
      ids <- intersect(s[[which]], colnames(dist))   # guild split + drop unpooled
      n_k <- length(ids)
      if (n_k < 2) return(c(n_boot_taxa = n_k, obs = NA_real_, null_mean = NA_real_, ses = NA_real_))
      ns  <- null_stats(guild, basis, metric, n_k, runs)
      obs <- stat_set(dist, ids, fun)
      c(n_boot_taxa = n_k, obs = obs, null_mean = unname(ns["mean"]),
        ses = unname((obs - ns["mean"]) / ns["sd"]))
    })
    as_tibble(do.call(rbind, out))
  }

  # B4/B5 fix: summarise the per-replicate SES against the ROPE evaluated at
  # THAT replicate's own k (rope_ses_at, cached -- the size-mismatch fix) into
  # a bootstrap 95% interval (raw + SES scale) and the continuous overlap
  # measures, now computed entirely on the SES scale. k_boot_lo/hi (2.5th/
  # 97.5th percentile of replicate set size) makes the size spread visible.
  summarise_boot <- function(guild, basis, metric, rep_tbl, alpha = BOOT_ALPHA) {
    v <- rep_tbl %>% filter(!is.na(obs))
    if (!nrow(v)) {
      return(tibble(metric_boot_lo = NA_real_, metric_boot_hi = NA_real_,
                    ses_boot_lo = NA_real_, ses_boot_hi = NA_real_, ses_boot_median = NA_real_,
                    pct_below_rope = NA_real_, pct_above_rope = NA_real_, pct_outside_rope = NA_real_,
                    n_boot_taxa_median = NA_real_, k_boot_lo = NA_real_, k_boot_hi = NA_real_,
                    boot_n_valid = 0L, boot_n_total = nrow(rep_tbl)))
    }
    rope_b <- t(vapply(v$n_boot_taxa, function(k) rope_ses_at(guild, basis, metric, k), numeric(2)))
    tibble(
      metric_boot_lo   = unname(quantile(v$obs, alpha / 2)),
      metric_boot_hi   = unname(quantile(v$obs, 1 - alpha / 2)),
      ses_boot_lo      = unname(quantile(v$ses, alpha / 2)),
      ses_boot_hi      = unname(quantile(v$ses, 1 - alpha / 2)),
      ses_boot_median  = median(v$ses),
      pct_below_rope   = 100 * mean(v$ses < rope_b[, 1]),
      pct_above_rope   = 100 * mean(v$ses > rope_b[, 2]),
      pct_outside_rope = 100 * mean(v$ses < rope_b[, 1] | v$ses > rope_b[, 2]),
      n_boot_taxa_median = median(v$n_boot_taxa),
      k_boot_lo = unname(quantile(v$n_boot_taxa, 0.025)),
      k_boot_hi = unname(quantile(v$n_boot_taxa, 0.975)),
      boot_n_valid = nrow(v), boot_n_total = nrow(rep_tbl)
    )
  }

  boot_cfg <- list(
    list(guild = GUILD_LABELS[["gamf"]], sets = boot_amf),
    list(guild = GUILD_LABELS[["mamf"]], sets = boot_amf),
    list(guild = GUILD_LABELS[["emf"]],  sets = boot_emf)
  )
  boot_tbl <- map_dfr(boot_cfg, function(g) {
    map_dfr(c("Phylogeny", "Lineage"), function(basis) {
      map_dfr(c("neg", "pos"), function(dir) {
        set_name <- if (dir == "neg") "Nitrophobic (z-)" else "Nitrophilic (z+)"
        map_dfr(c("SES-MPD", "SES-MNTD"), function(metric_label) {
          k_obs    <- results$n[results$guild == g$guild & results$basis == basis &
                                   results$set == set_name & results$metric == metric_label]
          rope_obs <- rope_ses_at(g$guild, basis, metric_label, k_obs)
          rep_tbl  <- boot_replicates(g$guild, basis, metric_label, g$sets, dir)
          summarise_boot(g$guild, basis, metric_label, rep_tbl) %>%
            mutate(guild = g$guild, basis = basis, set = set_name, metric = metric_label,
                   rope_ses_lo = rope_obs[1], rope_ses_hi = rope_obs[2])
        })
      })
    })
  }) %>% mutate(basis = factor(basis, levels = c("Phylogeny", "Lineage")))

  results <- results %>%
    left_join(boot_tbl, by = c("guild", "basis", "set", "metric"))

  # ─────────────────────────────────────────────────────────────────────────
  # B4. Empirical ROPE classification, SES scale, k-matched + B8 taxa floor --
  # ─────────────────────────────────────────────────────────────────────────
  # Mutually exclusive and exhaustive. Mapping to the manuscript's four
  # categories, stated once: clustered -> negative, overdispersed -> positive,
  # random -> neutral, unresolved -> unresolved.
  classify_phylo <- function(boot_lower, boot_upper, rope_lower, rope_upper) {
    dplyr::case_when(
      is.na(boot_lower) | is.na(boot_upper) ~ NA_character_,
      boot_upper < rope_lower ~ "clustered",
      boot_lower > rope_upper ~ "overdispersed",
      boot_lower >= rope_lower & boot_upper <= rope_upper ~ "random",
      TRUE ~ "unresolved"
    )
  }
  CLASS_TO_MANUSCRIPT <- c(clustered = "negative", overdispersed = "positive",
                           random = "neutral", unresolved = "unresolved")

  results <- results %>%
    mutate(
      classification = classify_phylo(ses_boot_lo, ses_boot_hi, rope_ses_lo, rope_ses_hi),
      # classification_fixed_k: the pre-fix behaviour (raw metric scale, one
      # ROPE at the observed k applied to every replicate regardless of its
      # own size) -- retained ONLY so sensitivity_rope can show the effect of
      # the fix (Classification_empirical_null_fixed_k); not used elsewhere.
      classification_fixed_k = classify_phylo(metric_boot_lo, metric_boot_hi, null_lo, null_hi),
      note = if_else(n < MIN_TAXA_FOR_INFERENCE,
                     sprintf("n = %d taxa (< MIN_TAXA_FOR_INFERENCE = %d): classification forced to unresolved",
                             n, MIN_TAXA_FOR_INFERENCE),
                     NA_character_),
      classification = if_else(n < MIN_TAXA_FOR_INFERENCE, "unresolved", classification),
      classification_fixed_k = if_else(n < MIN_TAXA_FOR_INFERENCE, "unresolved", classification_fixed_k)
    )

  # -- Phase 5 validation: is the interval trustworthy? ------------------------
  cat("\n=== Bootstrap CI validation ===\n")

  cat("1. ses_boot_median vs observed SES (should track; large offset => bias):\n")
  results %>% filter(!is.na(ses), !is.na(ses_boot_median)) %>%
    transmute(guild, basis, set, metric, ses = round(ses, 2),
              ses_boot_median = round(ses_boot_median, 2),
              offset = round(ses_boot_median - ses, 2)) %>%
    as.data.frame() %>% print(row.names = FALSE)

  cat(sprintf("\n2. Degeneracy: worst set dropped %.1f%% of replicates for k < 2.\n",
              100 * max(1 - results$boot_n_valid / results$boot_n_total, na.rm = TRUE)))
  deg <- results %>%
    mutate(drop = round(100 * (1 - boot_n_valid / boot_n_total), 1)) %>%
    filter(drop > 10) %>% distinct(guild, basis, set, drop)
  if (nrow(deg)) { cat("   sets losing > 10% of replicates:\n")
    print(as.data.frame(deg), row.names = FALSE)
  } else cat("   none above 10% — every set well populated.\n")

  cat("\n3. Per-replicate indicator-set size vs observed n (B2.4; distribution belongs in the supplement):\n")
  results %>% filter(!is.na(n_boot_taxa_median)) %>%
    distinct(guild, basis, set, n, n_boot_taxa_median, k_boot_lo, k_boot_hi) %>%
    as.data.frame() %>% print(row.names = FALSE)

  cat("\n4. Empirical-ROPE classification (SES scale, k-matched) vs rank-p call (should broadly agree):\n")
  results %>% filter(!is.na(ses_boot_lo)) %>%
    transmute(guild, basis, set, metric,
              ses_boot = sprintf("[%.2f, %.2f]", ses_boot_lo, ses_boot_hi),
              rope_ses = sprintf("[%.2f, %.2f]", rope_ses_lo, rope_ses_hi),
              classification,
              p_rank = round(p_rank, 3),
              rank_call = case_when(p_rank <= ALPHA / 2     ~ "clustered",
                                    p_rank >= 1 - ALPHA / 2 ~ "overdispersed",
                                    TRUE                    ~ "unresolved")) %>%
    as.data.frame() %>% print(row.names = FALSE)

  cat(sprintf("\n5. null_cache: %d distinct (guild x basis x metric x k) nulls drawn.\n",
              length(ls(null_cache))))

  # ─────────────────────────────────────────────────────────────────────────
  # Acceptance tests (ses_rope_scale_fix_spec.md §4) --------------------------
  # ─────────────────────────────────────────────────────────────────────────
  cat("\n=== Acceptance tests (ROPE scale fix) ===\n")

  # Test 1: the metric split must disappear. Under the pre-fix code every MPD
  # row's pct_outside_rope fell in ~0-10% and every MNTD row in ~7-93%, purely
  # from the size mismatch. Printed for every result row that has a
  # pct_outside_rope value (primary + both sensitivity bases/metrics).
  split_summary <- results %>% filter(!is.na(pct_outside_rope)) %>%
    group_by(metric) %>%
    summarise(min = min(pct_outside_rope), max = max(pct_outside_rope),
              median = median(pct_outside_rope), .groups = "drop")
  cat("Test 1 (metric split must disappear) -- pct_outside_rope range by metric:\n")
  print(as.data.frame(split_summary), row.names = FALSE)
  if (nrow(split_summary) == 2 && max(split_summary$min) > min(split_summary$max)) {
    cat("Test 1 FAIL: SES-MPD and SES-MNTD pct_outside_rope ranges do not overlap at all --\n")
    cat("             something is still comparing across mismatched sizes.\n")
  } else {
    cat("Test 1 PASS: no clean range separation between SES-MPD and SES-MNTD.\n")
  }

  # Test 2: ses_obs must fall within SES_boot_95 for every row. Hard invariant.
  t2 <- results %>% filter(!is.na(ses), !is.na(ses_boot_lo)) %>%
    mutate(inside = ses >= ses_boot_lo & ses <= ses_boot_hi)
  if (!all(t2$inside)) {
    print(as.data.frame(filter(t2, !inside) %>%
                          select(guild, basis, set, metric, ses, ses_boot_lo, ses_boot_hi)))
  }
  stopifnot("Test 2 FAILED: ses_obs falls outside SES_boot_95 for at least one row (see printed rows above)" = all(t2$inside))
  cat(sprintf("Test 2 PASS: ses_obs falls within SES_boot_95 for all %d rows.\n", nrow(t2)))

  # Test 3: for a k whose null is approximately normal (Shapiro-Wilk on a
  # 5000-draw subsample, p > 0.1), rope_ses_at() should return ~= +/-1.96.
  keys  <- ls(null_cache)
  found <- FALSE
  for (key in keys) {
    nd  <- get(key, envir = null_cache)
    sub <- if (length(nd) > 5000) nd[seq_len(5000)] else nd
    sw  <- tryCatch(stats::shapiro.test(sub), error = function(e) NULL)
    if (!is.null(sw) && sw$p.value > 0.1) {
      parts <- strsplit(key, "\\|")[[1]]
      r <- rope_ses_at(parts[1], parts[2], parts[3], as.integer(parts[4]))
      cat(sprintf("Test 3: key '%s' (Shapiro-Wilk p = %.3f) -> rope_ses_at = [%.2f, %.2f] (expect ~= +/-1.96)\n",
                  key, sw$p.value, r[1], r[2]))
      found <- TRUE
      break
    }
  }
  if (!found) cat("Test 3: no cached null passed the normality check (p > 0.1) among those drawn -- no worked example available.\n")

  # Test 4: every row carries exactly one classification; NA only permitted
  # where MIN_TAXA_FOR_INFERENCE forces "unresolved" (which is a string, not
  # NA, in this design) -- i.e. NA should never occur for n >= the floor.
  bad_na <- results %>% filter(n >= MIN_TAXA_FOR_INFERENCE, is.na(classification))
  if (nrow(bad_na) > 0) print(as.data.frame(bad_na %>% select(guild, basis, set, metric, n, classification)))
  stopifnot("Test 4 FAILED: NA classification for a row with n >= MIN_TAXA_FOR_INFERENCE (see printed rows above)" = nrow(bad_na) == 0)
  cat("Test 4 PASS: no unexplained NA classification.\n")

  # Test 6: every k present in the bootstrap has a cached (generated-on-demand)
  # null -- no silent misses. get_null() has no code path that returns without
  # either a cache hit or a successful generate-and-cache (it stop()s instead),
  # so this checks that invariant held for every key actually requested.
  missed <- setdiff(ls(null_requested), ls(null_cache))
  if (length(missed) > 0) cat("Missing keys:", paste(head(missed, 5), collapse = ", "), "\n")
  stopifnot("Test 6 FAILED: null(s) requested during the bootstrap are missing from the cache" = length(missed) == 0)
  cat(sprintf("Test 6 PASS: %d distinct nulls requested and cached; no misses.\n", length(ls(null_cache))))
}

# ─────────────────────────────────────────────────────────────────────────────
# -- Does the call change under the normal approximation? ----------------------
# ─────────────────────────────────────────────────────────────────────────────
# The rank p is the reported inference; the normal approximation is what a
# fixed SES threshold would encode. They agree only while the null is roughly
# symmetric and Gaussian, which the discrete taxonomic distance (values 1, 2, 3)
# and the small indicator sets can both break.
calls <- results %>%
  filter(!is.na(ses)) %>%
  mutate(
    rank_call = case_when(p_rank <= ALPHA / 2     ~ "clustered",
                          p_rank >= 1 - ALPHA / 2 ~ "overdispersed",
                          TRUE                    ~ "unresolved"),
    norm_call = case_when(p_norm <= ALPHA / 2     ~ "clustered",
                          p_norm >= 1 - ALPHA / 2 ~ "overdispersed",
                          TRUE                    ~ "unresolved")
  )

changed <- filter(calls, rank_call != norm_call)
if (nrow(changed)) {
  cat(sprintf("\nCalls differing between rank and normal-approximation p (alpha = %.2f):\n", ALPHA))
  print(as.data.frame(select(changed, guild, basis, set, metric, n, ses,
                             p_rank, p_norm, rank_call, norm_call)))
  cat("Rank p is the reported inference; the null is non-Gaussian for these sets, which is\n")
  cat("exactly why classification (B4) is made against an empirical ROPE (evaluated on the\n")
  cat("SES scale, at each replicate's own indicator-set size) rather than a fixed SES threshold.\n")
} else {
  cat(sprintf("\nRank and normal-approximation p agree for all sets (alpha = %.2f).\n", ALPHA))
}

# ─────────────────────────────────────────────────────────────────────────────
# -- One-off check that the rewritten SES reproduces picante -------------------
# ─────────────────────────────────────────────────────────────────────────────
# Drawing k taxa at random from the pool and shuffling taxa labels are the same
# null for a single fixed-size set, so the two SES values should agree to Monte
# Carlo error. Set FIG6_VERIFY=FALSE once satisfied.
if (VERIFY_PICANTE) {
  ids     <- colnames(dist_gamf_phy)
  neg_ids <- mp_gamf$dat$otu_id[mp_gamf$dat$indicator == "neg"]
  pos_ids <- mp_gamf$dat$otu_id[mp_gamf$dat$indicator == "pos"]
  # Both rows, not just z-. ses.mpd errors on a single-row community: it passes
  # row.names(samp) to data.frame(), and a length-1 character there is read as a
  # column name rather than as the row name.
  chk_comm <- rbind("Nitrophobic (z-)" = as.integer(ids %in% neg_ids),
                    "Nitrophilic (z+)" = as.integer(ids %in% pos_ids))
  colnames(chk_comm) <- ids
  set.seed(SEED)
  chk_pic <- ses.mpd(chk_comm, dist_gamf_phy, null.model = "taxa.labels", runs = SES_RUNS)
  set.seed(SEED)
  chk_new <- ses_full(dist_gamf_phy, neg_ids, "mpd", runs = SES_RUNS, jack_runs = 0)
  cat(sprintf("\nVerification (%s z-, phylogeny): picante SES = %.3f, rewritten SES = %.3f\n",
              GUILD_LABELS[["gamf"]], chk_pic["Nitrophobic (z-)", "mpd.obs.z"], chk_new$ses))
}

# ─────────────────────────────────────────────────────────────────────────────
# -- B6. Output specification: primary / sensitivity_basis / sensitivity_metric
# -- / sensitivity_rope / metadata (+ resolution, session_info) ----------------
# ─────────────────────────────────────────────────────────────────────────────
if (BOOT_SES) {

  base_cols <- function(df) {
    df %>%
      transmute(
        Guild = guild, Direction = set, k_observed = n, Metric = metric,
        Metric_observed = round(obs, 3),
        Null_mean       = round(null_mean, 3),
        `ROPE (raw scale)`     = sprintf("[%.3f, %.3f]", null_lo, null_hi),
        `Boot_95 (raw scale)`  = sprintf("[%.3f, %.3f]", metric_boot_lo, metric_boot_hi),
        SES             = round(ses, 2),
        SES_boot_95     = sprintf("[%.2f, %.2f]", ses_boot_lo, ses_boot_hi),
        ROPE_ses        = sprintf("[%.2f, %.2f]", rope_ses_lo, rope_ses_hi),
        k_boot_range    = sprintf("[%d, %d]", as.integer(round(k_boot_lo)), as.integer(round(k_boot_hi))),
        p_rank = round(p_rank, 4), p_norm = round(p_norm, 4),
        Jackknife = sprintf("[%.2f, %.2f]", jack_min, jack_max),
        Classification = classification,
        pct_outside_rope = round(pct_outside_rope, 1),
        pct_below_rope   = round(pct_below_rope, 1),
        pct_above_rope   = round(pct_above_rope, 1),
        n_boot_taxa_median = n_boot_taxa_median,
        note = note
      ) %>%
      arrange(factor(Guild, levels = guild_levels), desc(Direction))
  }

  results <- results %>%
    mutate(m = recode(metric, "SES-MPD" = "MPD", "SES-MNTD" = "MNTD"),
           is_primary_metric = m == PRIMARY_METRIC[guild])

  # primary_raw: PRIMARY_BASIS x each guild's PRIMARY_METRIC (B1), full numeric
  # columns -- this is what code/figure_6.R plots (via figure_6.Rdata below).
  # `primary` (the B6 Excel sheet) is a rounded/formatted display version of
  # the same rows, built from it just below.
  primary_raw <- results %>% filter(basis == PRIMARY_BASIS, is_primary_metric)
  primary     <- primary_raw %>% base_cols()

  # sensitivity_basis: the other basis, same (pre-specified) metric per guild
  sensitivity_basis <- results %>%
    filter(basis != PRIMARY_BASIS, is_primary_metric) %>% base_cols()

  # sensitivity_metric: PRIMARY_BASIS, the non-pre-specified metric per guild --
  # what MNTD says for G-AMF/M-AMF, what MPD says for EMF
  sensitivity_metric <- results %>%
    filter(basis == PRIMARY_BASIS, !is_primary_metric) %>% base_cols()

  # sensitivity_rope (B7): classification under ROPE_METHOD = "ses_fixed" at
  # SES thresholds 0/1/2/2.5, alongside the empirical-null (primary) result.
  # Reuses the SES-scale bootstrap CI (ses_boot_lo/hi) already computed above,
  # against symmetric bounds +/- threshold, via the same classify_phylo() (B4)
  # defined inside the BOOT_SES block earlier in this script.
  #
  # Classification_empirical_null is now the SES-scale, k-matched result
  # (post ROPE-scale-fix). Classification_empirical_null_fixed_k holds the
  # PRE-fix behaviour (raw metric scale, single ROPE at the observed k applied
  # to every replicate) so the effect of the fix is visible in the file, not
  # only in a commit message (ses_rope_scale_fix_spec.md §3.3).
  primary_rows <- primary_raw
  sensitivity_rope <- primary_rows %>%
    transmute(Guild = guild, Direction = set, n_taxa = n, Metric = metric, SES = round(ses, 2),
              Classification_empirical_null = if_else(n < MIN_TAXA_FOR_INFERENCE, "unresolved", classification),
              Classification_empirical_null_fixed_k = if_else(n < MIN_TAXA_FOR_INFERENCE, "unresolved", classification_fixed_k))
  for (thr in SES_FIXED_THRESHOLDS) {
    call <- classify_phylo(primary_rows$ses_boot_lo, primary_rows$ses_boot_hi, -thr, thr)
    call <- if_else(primary_rows$n < MIN_TAXA_FOR_INFERENCE, "unresolved", call)
    sensitivity_rope[[sprintf("Classification_ses_fixed_%s", format(thr))]] <- call
  }
  class_mat <- as.matrix(select(sensitivity_rope, starts_with("Classification_ses_fixed_")))
  sensitivity_rope$any_change <- apply(class_mat != sensitivity_rope$Classification_empirical_null, 1, any)
  sensitivity_rope <- sensitivity_rope %>%
    arrange(factor(Guild, levels = guild_levels), desc(Direction))

  if (any(sensitivity_rope$any_change)) {
    cat("\nB7: classification is threshold-dependent -- rows whose label changes under a fixed SES ROPE:\n")
    print(as.data.frame(filter(sensitivity_rope, any_change)), row.names = FALSE)
  } else {
    cat("\nB7: classification is unchanged across all SES-fixed thresholds tested.\n")
  }

  # metadata (B6): pool definition, prevalence filter, tree source, null model,
  # N_BOOT, SEED, per-replicate indicator-set definition, date. sessionInfo() is
  # a separate sheet (session_info) since it does not fit the Field/Value shape.
  metadata_tbl <- tibble(
    Field = c(
      "Pool definition", "Prevalence filter",
      "Tree source (AMF)", "Tree source (EMF)", "Null model",
      "Observed indicator-set definition", "Per-replicate indicator-set definition",
      "PRIMARY_BASIS", sprintf("PRIMARY_METRIC (%s)", names(PRIMARY_METRIC)),
      "ROPE_METHOD", "ROPE definition", "Rationale",
      "ROPE_NULL_Q", "MIN_TAXA_FOR_INFERENCE", "N_BOOT", "SEED",
      "SES_RUNS (null draws per statistic)", "JACK_RUNS (null draws per jackknife fold)",
      "Residual limitation", "Date generated"
    ),
    Value = c(
      "Taxa entering that guild's TITAN2 run (sppmax rows), further restricted to taxa with a non-missing family/lineage assignment (matched pool, see console output above for drop counts)",
      "n_samples >= 5 (out of ~123-126 sites), applied in 02a_titan_amf.R / 02b_titan_emf.R before TITAN2 -- this filter defines the SES null pool and is not re-tested at an alternative threshold here (C.3)",
      "data/amf/tree.newick",
      "data/emf/tree.newick",
      "Random draw of k taxa without replacement from the guild- and basis-specific candidate pool (equivalent to picante's taxa.labels null for one fixed k); drawn independently per (guild, basis, metric, k) with a deterministic per-key seed (hash of the key, not stream order) and cached",
      "purity >= 0.95 & reliability >= 0.90, both computed across all TITAN2 bootstraps",
      "direction = maxgrp, filtered on within-replicate IndVal p (obsiv.prob) <= 0.05 -- noisier than the cross-replicate purity/reliability filter, so replicate sets typically run 2-4x larger than the observed set",
      PRIMARY_BASIS, unname(PRIMARY_METRIC),
      ROPE_METHOD,
      "Central 95% of the null distribution at each replicate's own indicator-set size, expressed in SES units. Asymmetric where the null is skewed.",
      "Null distributions depart from normality (p_rank vs p_norm disagree by 2-3x), and MNTD is strongly size-dependent while bootstrap replicate sets are larger than the observed set. A single ROPE at the observed k is therefore not comparable across replicates.",
      paste(ROPE_NULL_Q, collapse = " - "), as.character(MIN_TAXA_FOR_INFERENCE),
      as.character(N_BOOT), as.character(SEED),
      as.character(SES_RUNS), as.character(JACK_RUNS),
      "The k-matched ROPE corrects the SIZE mismatch in the comparison, but the bootstrap distribution still describes a more permissively defined taxon set than the observed one (see 'Per-replicate indicator-set definition' above) -- replicate sets run 2-4x larger. This is a real property of the estimand, not an artefact the ROPE fix removes.",
      as.character(Sys.Date())
    )
  )
  session_info_tbl <- tibble(session_info = capture.output(sessionInfo()))

  write_xlsx(
    list(primary = primary,
         sensitivity_basis = sensitivity_basis,
         sensitivity_metric = sensitivity_metric,
         sensitivity_rope = sensitivity_rope,
         resolution = freedom,
         metadata = metadata_tbl,
         session_info = session_info_tbl),
    "output/table_6_ses.xlsx"
  )
  cat("\nprimary sheet (SES, null band / ROPE on the observed scale, classification):\n")
  print(as.data.frame(primary))
  cat("\nSaved output/table_6_ses.xlsx\n")

} else {
  cat("\nFIG6_BOOTSES=FALSE: skipping B3-B8 (raw-metric storage, empirical ROPE,\n")
  cat("classification, sensitivity_rope) and the primary/sensitivity output workbook,\n")
  cat("since all of them depend on the bootstrap replicates.\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# -- Save plotting data for code/figure_6.R ------------------------------------
# ─────────────────────────────────────────────────────────────────────────────
# results        : SES, null band/ROPE, jackknife range, bootstrap CIs and
#                  classification (both bases)
# primary_raw    : results filtered to PRIMARY_BASIS x each guild's
#                  PRIMARY_METRIC (B1), full numeric columns -- what
#                  figure_6.R plots. `primary` (the B6 Excel sheet) is a
#                  rounded/formatted version of the same rows and is not
#                  useful for plotting. Only exists when BOOT_SES = TRUE (it
#                  depends on the bootstrap CIs).
# dat_* / tree_* : full candidate-pool frames and pruned trees for the
#                  supplementary change-point tree figures
if (!dir.exists("generated_data")) dir.create("generated_data")
save_objs <- c("results", "freedom", "guild_levels", "GUILD_LABELS",
               "dat_gamf", "dat_mamf", "dat_emf",
               "tree_gamf", "tree_mamf", "tree_emf")
if (BOOT_SES) save_objs <- c(save_objs, "primary_raw")
save(list = save_objs, file = "generated_data/figure_6.Rdata")
cat("Saved generated_data/figure_6.Rdata\n")

cat("\n=== PHYLOGENETIC ANALYSIS COMPLETE ===\n")
