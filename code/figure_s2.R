library(dagitty)
library(ggdag)
library(ggplot2)
library(dplyr)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────
# DAG — MYCORRHIZAL FUNGI
#
# This DAG represents our hypothesised causal structure governing arbuscular
# mycorrhizal (AM) and ectomycorrhizal (ECM) diversity across environmental
# gradients.
#
# Mineral nitrogen is the exposure of interest. Mycorrhizal fungal diversity
# is the outcome. Aridity remains in the DAG as an upstream driver and
# potential confounder, but is not treated as an exposure here.
#
# LIMITATIONS
# Two feedback loops exist that cannot be represented in an acyclic graph and
# cannot be resolved with cross-sectional observational data:
#
#   (1) Carbon <-> Saprotrophs: Soil organic carbon is food for saprotrophic
#       fungi, but saprotrophic activity in turn governs carbon pool dynamics
#       through decomposition. We retained the Saprotrophs -> Carbon direction
#       as we are primarily concerned with how the decomposer community mediates
#       downstream effects on mycorrhizal diversity.
#
#   (2) ECM -> Saprotrophs (Gadgil effect): ECM fungi can suppress
#       saprotrophic communities through competitive exclusion and reduced
#       carbon availability. This creates cycles back to ECM via soil nutrient
#       paths. Both ECM and Saprotrophs are therefore treated as sharing
#       common upstream drivers rather than being directly causally linked.
#       This limitation applies to ECM only — AM fungi do not exhibit the
#       Gadgil effect.
#
# PATHS AND MECHANISMS
#
# ── CLIMATE DRIVERS ──────────────────────────────────────────────────────────
# PATH:       Temperature / Precipitation -> Aridity
# MECHANISM:  Higher temperatures increase evapotranspiration and lower
#             precipitation increases water deficit, together driving aridity.
#
# PATH:       Temperature / Precipitation -> Trees
# MECHANISM:  Climate directly filters tree species tolerances, growth rates,
#             and recruitment success independently of soil-mediated pathways.
#
# PATH:       Temperature / Precipitation -> Nitrogen, Carbon, pH, Phosphorus
# MECHANISM:  Climate governs decomposition rates, weathering, leaching, and
#             organic matter accumulation, directly altering soil chemistry.
#
# PATH:       Temperature / Precipitation -> Fire
# MECHANISM:  Temperature and precipitation deficit determine fuel moisture and
#             ignition probability, driving fire frequency and severity.
#
# PATH:       Temperature / Precipitation -> Pathogens
# MECHANISM:  Climatic conditions govern pathogen life cycles, dispersal,
#             host susceptibility, and epidemic potential.
#
# PATH:       Temperature / Precipitation -> Saprotrophs
# MECHANISM:  Saprotrophic fungal activity is strongly governed by temperature
#             and soil moisture, which control enzymatic decomposition rates
#             and hyphal growth.
#
# PATH:       Temperature / Precipitation -> Mycorrhizal
# MECHANISM:  Climate directly filters mycorrhizal community composition
#             through physiological tolerances, fruiting phenology, and
#             hyphal growth rates.
#
# ── ARIDITY ──────────────────────────────────────────────────────────────────
# PATH:       Aridity -> Nitrogen
# MECHANISM:  Aridity reduces decomposition and mineralisation rates, lowering
#             plant-available nitrogen in drier soils.
#
# PATH:       Aridity -> Carbon
# MECHANISM:  Aridity slows decomposition, allowing organic carbon to
#             accumulate, and reduces plant productivity, lowering carbon
#             inputs.
#
# PATH:       Aridity -> pH
# MECHANISM:  Aridity reduces leaching of base cations, typically raising soil
#             pH in arid systems relative to humid ones.
#
# PATH:       Aridity -> Phosphorus
# MECHANISM:  Aridity reduces weathering and mineralisation, limiting
#             plant-available phosphorus.
#
# PATH:       Aridity -> Trees
# MECHANISM:  Aridity imposes direct water stress on trees, filtering species
#             composition and reducing productivity and abundance.
#
# PATH:       Aridity -> Fire
# MECHANISM:  Aridity increases fuel dryness and fire weather, driving fire
#             frequency and severity.
#
# PATH:       Aridity -> Pathogens
# MECHANISM:  Aridity alters host tree stress susceptibility and pathogen
#             environmental tolerances, shifting pest and pathogen pressure.
#
# PATH:       Aridity -> Saprotrophs
# MECHANISM:  Soil moisture is a primary control on saprotrophic fungal
#             activity. Aridity reduces decomposer biomass and activity
#             through desiccation stress.
#
# PATH:       Aridity -> Mycorrhizal
# MECHANISM:  Aridity directly influences mycorrhizal community composition
#             through hyphal drought tolerance and host tree water stress
#             responses that alter root exudate allocation to mycorrhizal
#             partners.
#
# ── FIRE ─────────────────────────────────────────────────────────────────────
# PATH:       Fire -> Trees
# MECHANISM:  Fire causes direct tree mortality, alters species composition
#             through differential fire tolerance, and resets stand structure.
#
# PATH:       Fire -> Carbon
# MECHANISM:  Fire volatilises organic carbon and converts litter and biomass
#             to char, substantially reducing soil organic carbon pools.
#
# PATH:       Fire -> Nitrogen
# MECHANISM:  Fire volatilises nitrogen from biomass and litter, causing large
#             pulse losses of soil nitrogen post-fire.
#
# PATH:       Fire -> pH
# MECHANISM:  Ash deposition following fire temporarily raises soil pH through
#             addition of base cations.
#
# PATH:       Fire -> Saprotrophs
# MECHANISM:  Fire causes large shifts in saprotrophic communities. Post-fire
#             char and altered substrate chemistry drives succession of
#             specialist saprotrophic taxa and suppresses others. ECM
#             communities include well-documented post-fire specialist taxa
#             (e.g. Rhizopogon, Morchella) that are particularly responsive
#             to fire history.
#
# ── PATHOGENS ────────────────────────────────────────────────────────────────
# PATH:       Pathogens -> Trees
# MECHANISM:  Pest and pathogen pressure causes direct tree mortality and
#             suppresses recruitment, reducing host diversity and abundance
#             for mycorrhizal fungi.
#
# PATH:       Pathogens -> Carbon
# MECHANISM:  Tree mortality from pests and pathogens increases necromass and
#             deadwood inputs, elevating soil organic carbon pools.
#
# NOTE:       Pathogens -> Mycorrhizal direct path NOT included. Mycorrhizal
#             communities are governed by host plant physiology; pathogen
#             effects are mediated through tree health and abundance and soil
#             carbon inputs, not through direct pathogen-fungus interactions.
#
# ── SAPROTROPHS ──────────────────────────────────────────────────────────────
# PATH:       Saprotrophs -> Carbon
# MECHANISM:  Saprotrophic fungi are the primary decomposers of organic matter.
#             High saprotrophic activity reduces soil organic carbon pools
#             through mineralisation; low activity allows carbon accumulation.
#
# PATH:       Saprotrophs -> Nitrogen
# MECHANISM:  Saprotrophic fungi mineralise organic nitrogen during
#             decomposition, releasing plant-available inorganic nitrogen.
#             They can also immobilise nitrogen in high C:N substrates,
#             reducing availability.
#
# PATH:       Saprotrophs -> Phosphorus
# MECHANISM:  Saprotrophic fungi mineralise organic phosphorus through
#             phosphatase enzyme activity, releasing plant-available inorganic
#             phosphorus.
#
# NOTE:       Carbon -> Saprotrophs NOT INCLUDED. A bidirectional relationship
#             between Carbon and Saprotrophs exists ecologically but retaining
#             both directions violates acyclicity. Saprotrophs -> Carbon is
#             retained. See limitations above.
#
# NOTE:       For ECM specifically, the Gadgil effect also exists but creates
#             unresolvable cycles — see limitations above.
#
# ── SOIL CHEMISTRY ───────────────────────────────────────────────────────────
# PATH:       Carbon -> Nitrogen
# MECHANISM:  Soil organic carbon drives microbial mineralisation of organic
#             nitrogen into plant-available forms. C:N ratio governs net
#             mineralisation versus immobilisation. A separate C:N ratio node
#             is not included as its effect is fully captured by this path.
#
# PATH:       Carbon -> Phosphorus
# MECHANISM:  In highly weathered Australian soils, organic matter is a
#             primary driver of phosphorus mineralisation through microbial
#             phosphatase activity.
#
# PATH:       pH -> Nitrogen
# MECHANISM:  Soil pH governs nitrification and mineralisation rates. Acidic
#             soils suppress nitrifying bacteria, reducing plant-available
#             nitrogen.
#
# PATH:       Nitrogen -> Trees
# MECHANISM:  Nitrogen is a primary limiting nutrient for tree growth.
#             Nitrogen availability directly determines tree productivity and
#             species composition.
#
# PATH:       pH -> Trees
# MECHANISM:  Soil pH directly filters tree species tolerances and influences
#             nutrient availability and aluminium toxicity.
#
# PATH:       Phosphorus -> Trees
# MECHANISM:  Phosphorus is co-limiting with nitrogen in many systems,
#             particularly highly weathered Australian soils.
#
# ── MYCORRHIZAL FUNGI ────────────────────────────────────────────────────────
# PATH:       Trees -> Mycorrhizal
# MECHANISM:  Mycorrhizal fungi are obligate symbionts dependent on host plant
#             carbon. Tree species identity, diversity, and abundance are
#             primary determinants of mycorrhizal community composition,
#             richness, and abundance.
#
# PATH:       Nitrogen -> Mycorrhizal
# MECHANISM:  Mycorrhizal fungi are sensitive to nitrogen availability.
#             Elevated nitrogen suppresses diversity through reduced host
#             carbon allocation to mycorrhizal partners and direct competitive
#             shifts among taxa. This effect is particularly strong for ECM
#             fungi which are highly sensitive to nitrogen deposition.
#
# PATH:       pH -> Mycorrhizal
# MECHANISM:  Mycorrhizal community composition is strongly filtered by soil
#             pH through direct physiological constraints on hyphal growth
#             and indirect effects on nutrient availability.
#
# PATH:       Phosphorus -> Mycorrhizal
# MECHANISM:  AM fungi are the classical phosphorus acquisition symbiont and
#             respond strongly to phosphorus availability. ECM fungi also
#             respond to phosphorus — some taxa produce phosphatase enzymes
#             and actively acquire phosphorus, particularly in phosphorus-
#             limited systems such as highly weathered Australian soils.
#
# ─────────────────────────────────────────────────────────────────────────────

g_myc <- dagitty('
dag {

  Temperature -> Mycorrhizal
  Temperature -> Aridity
  Temperature -> Trees
  Temperature -> Nitrogen
  Temperature -> Carbon
  Temperature -> pH
  Temperature -> Phosphorus
  Temperature -> Fire
  Temperature -> Pathogens
  Temperature -> Saprotrophs

  Precipitation -> Mycorrhizal
  Precipitation -> Aridity
  Precipitation -> Trees
  Precipitation -> Nitrogen
  Precipitation -> Carbon
  Precipitation -> pH
  Precipitation -> Phosphorus
  Precipitation -> Fire
  Precipitation -> Pathogens
  Precipitation -> Saprotrophs

  Aridity -> Nitrogen
  Aridity -> Carbon
  Aridity -> pH
  Aridity -> Phosphorus
  Aridity -> Trees
  Aridity -> Mycorrhizal
  Aridity -> Fire
  Aridity -> Pathogens
  Aridity -> Saprotrophs

  Fire -> Trees
  Fire -> Carbon
  Fire -> Nitrogen
  Fire -> pH

  Pathogens -> Trees

  Carbon -> Nitrogen
  Carbon -> Phosphorus

  Saprotrophs -> Carbon
  Saprotrophs -> Nitrogen
  Saprotrophs -> Phosphorus

  pH -> Nitrogen

  Nitrogen -> Trees
  pH -> Trees
  Phosphorus -> Trees

  Nitrogen -> Mycorrhizal
  pH -> Mycorrhizal
  Phosphorus -> Mycorrhizal

  Trees -> Mycorrhizal

  Nitrogen    [exposure=TRUE]
  Mycorrhizal [outcome=TRUE]
}
')

# ─────────────────────────────────────────────────────────────────────────────
# 2. Check DAG is acyclic
# ─────────────────────────────────────────────────────────────────────────────
stopifnot(isAcyclic(g_myc))
cat("DAG is acyclic: OK\n")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Adjustment sets — direct and total effects of Nitrogen -> Mycorrhizal
# ─────────────────────────────────────────────────────────────────────────────
adj_dir_nitrogen <- adjustmentSets(g_myc, exposure = "Nitrogen", outcome = "Mycorrhizal", effect = "direct")[[1]]
adj_tot_nitrogen <- adjustmentSets(g_myc, exposure = "Nitrogen", outcome = "Mycorrhizal", effect = "total")[[1]]

cat("\nDirect effect adjustment set — Nitrogen -> Mycorrhizal:\n")
print(adj_dir_nitrogen)

cat("\nTotal effect adjustment set — Nitrogen -> Mycorrhizal:\n")
print(adj_tot_nitrogen)

# ─────────────────────────────────────────────────────────────────────────────
# 4. Implied conditional independencies — sanity check
# Cross-reference against ecological knowledge. Any surprising independency
# may indicate a missing path or incorrectly specified path direction.
# ─────────────────────────────────────────────────────────────────────────────
cat("\nImplied conditional independencies:\n")
print(impliedConditionalIndependencies(g_myc))

# ─────────────────────────────────────────────────────────────────────────────
# 5. Helper function — label nodes by role for plotting
# ─────────────────────────────────────────────────────────────────────────────
label_nodes <- function(graph, exposure, outcome, adj_set) {
  tidy_dagitty(graph) %>%
    mutate(role = case_when(
      name == exposure   ~ "Exposure",
      name == outcome    ~ "Outcome",
      name %in% adj_set  ~ "Adjust",
      TRUE               ~ "Other"
    ))
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. Shared colour scale and theme
# ─────────────────────────────────────────────────────────────────────────────
myc_colours <- scale_colour_manual(
  name   = NULL,
  values = c(
    "Exposure"           = "#8da0cb",
    "Outcome"            = "#fc8d62",
    "Adjustment set"     = "#66c2a5",
    "Other covariate"    = "grey90"
  ),
  breaks = c(
    "Exposure",
    "Outcome",
    "Adjustment set",
    "Other covariate"
  ),
  guide = guide_legend(
    nrow         = 2,
    override.aes = list(size = 5)
  )
)

myc_theme <- list(
  theme_dag(),
  theme(
    plot.caption    = element_text(size = 8, colour = "grey40", hjust = 0.5),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    legend.title    = element_text(size = 10, face = "bold")
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Full DAG plot
# ─────────────────────────────────────────────────────────────────────────────
plot_myc_full <- ggdag(g_myc, layout = "sugiyama") +
  myc_theme +
  ggtitle("Full Mycorrhizal DAG") +
  labs(
    caption = paste(
      "Known feedback limitations: Carbon <-> Saprotrophs and ECM-specific Gadgil effect.",
      "\nBoth cycles cannot be represented in an acyclic graph with cross-sectional data."
    )
  )

plot_myc_full

# ─────────────────────────────────────────────────────────────────────────────
# 8. Direct effect adjustment set plot — Nitrogen -> Mycorrhizal
# ─────────────────────────────────────────────────────────────────────────────

# Nodes in the adjustment set (excluding the exposure itself, defensively)
nitrogen_adj <- setdiff(adj_dir_nitrogen, "Nitrogen")

# Build node labels
# Rename "Mycorrhizal" to "Mycorrhizal\nFungi" for the plot label only
set.seed(1986)
nitrogen_roles <- tidy_dagitty(g_myc) %>%
  mutate(
    role = case_when(
      name == "Nitrogen"     ~ "Exposure",
      name == "Mycorrhizal"  ~ "Outcome",
      name %in% nitrogen_adj ~ "Adjustment set",
      TRUE                   ~ "Other covariate"
    ),
    # Line break in outcome label for plot legibility
    label = if_else(name == "Mycorrhizal", "Mycorrhizal\nFungi", name)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 9. Plot direct effect adjustment set
# ─────────────────────────────────────────────────────────────────────────────
plot_nitrogen_direct <- nitrogen_roles %>%
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges(edge_colour = "grey70") +
  geom_dag_point(aes(colour = role), shape = 19, size = 20) +
  geom_dag_text(
    aes(label = label),
    colour     = "black",
    fontface   = "bold",
    size       = 3,
    lineheight = 0.85   # tighten line spacing for two-line label
  ) +
  myc_colours +
  myc_theme

plot_nitrogen_direct

# ─────────────────────────────────────────────────────────────────────────────
# 10. Save
# ─────────────────────────────────────────────────────────────────────────────
ggsave(
  plot_nitrogen_direct,
  filename = "output/figure_s2.png",
  width    = 16,
  height   = 16,
  units    = "cm",
  dpi      = 300
)

