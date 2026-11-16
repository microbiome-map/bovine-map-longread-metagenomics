# 08_Farm_Adjusted_Alpha_Diversity.R
# Farm-adjusted sensitivity analyses for alpha diversity (03).
#
# Usage in RStudio: set the working directory to the folder that CONTAINS the
# extracted `data/` folder (Session > Set Working Directory > Choose Directory),
# or edit `data_dir` below, then run the whole script (Ctrl/Cmd + Shift + Enter).
#
# Alpha values are calculated exactly as in 03_Alpha_Diversity.R (MetaPhlAn GTDB
# normalized profiles; observed richness, Shannon and Gini-Simpson; genus and
# species). As in 07_Farm_Adjusted_Analyses.R, Swiss MAP-negative animals share
# no farm with MAP-positive animals, so the analyses are:
#   Reference: Kruskal-Wallis across the three Swiss MAP groups (as in 03).
#   A. Linear model value ~ Farm + MAP_group; the MAP_group term can only use
#      subclinical vs clinical within farms 6 and 7 (check its Df).
#   B. MAP-negative vs MAP-positive, Kruskal-Wallis H on animals with MAP status
#      permuted between whole farms (all labellings enumerated; exact).
#   B2. MAP-negative vs MAP-positive on farm means (one value per farm),
#      exact Wilcoxon rank-sum test.
#   B3. Optional (if lme4 is installed): linear mixed model with a random farm
#      intercept, likelihood-ratio test for MAP-negative vs MAP-positive.
#   C. Subclinical vs clinical within farms holding both, Kruskal-Wallis H with
#      permutations restricted within farm (exact).
#   South Africa: MAP-negative = farm 32, subclinical = farm 33; not separable.

# ---------------------------------------------------------------- settings ----
data_dir   <- "data"
output_dir <- "Farm_adjusted_alpha_diversity"
ranks      <- c("Genus", "Species")
metrics    <- c("Richness", "Shannon", "Gini_Simpson")

# ---------------------------------------------------------------- packages ----
packages <- c("dplyr", "readr", "tidyr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1),
                                     quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install first: install.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))")
}
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })
use_lme4 <- requireNamespace("lme4", quietly = TRUE)
if (!use_lme4) message("lme4 not installed: skipping the optional mixed model (B3).")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

taxonomy_file <- function(rank) {
  file.path(data_dir, "taxonomy", "MetaPhlAn_GTDB",
            paste0("MetaPhlAn_GTDB_", rank, "_normalized.tsv"))
}
for (f in c(file.path(data_dir, "metadata.txt"), taxonomy_file(ranks))) {
  if (!file.exists(f)) stop("Missing input: ", normalizePath(f, mustWork = FALSE),
                            "\nCheck the working directory / data_dir.")
}

# ---------------------------------------------------------------- metadata ----
eligible_groups <- c("Healthy" = "SW_NEG", "Positive, Subclinical" = "SW_SUB",
                     "Positive, Clinical" = "SW_POS", "ZA Healthy" = "ZA_NEG",
                     "ZA Positive, Subclinical" = "ZA_SUB")
excluded_ids <- c("SW_NEG_01", "ZA_NEG_01", "ZA_NEG_08")

meta <- read_tsv(file.path(data_dir, "metadata.txt"),
                 col_types = cols(.default = col_character())) %>%
  transmute(SAMPLE_ID = trimws(SAMPLE_ID),
            GROUP     = trimws(GROUP),
            Farm      = trimws(`Farm number`)) %>%
  filter(GROUP %in% names(eligible_groups), !SAMPLE_ID %in% excluded_ids) %>%
  mutate(Cohort_group = unname(eligible_groups[GROUP]),
         Country = if_else(startsWith(SAMPLE_ID, "SW_"), "Switzerland", "South Africa"),
         MAP_group = case_when(grepl("_NEG_", SAMPLE_ID) ~ "MAP-negative",
                               grepl("_SUB_", SAMPLE_ID) ~ "Subclinical",
                               grepl("_POS_", SAMPLE_ID) ~ "Clinical"),
         MAP_binary = if_else(MAP_group == "MAP-negative", "MAP_negative", "MAP_positive")) %>%
  arrange(SAMPLE_ID)

counts <- table(factor(meta$Cohort_group,
                       levels = c("SW_NEG", "SW_SUB", "SW_POS", "ZA_NEG", "ZA_SUB")))
if (!identical(as.integer(counts), c(16L, 13L, 7L, 6L, 8L))) {
  print(counts); stop("Expected 16/13/7/6/8 samples, as in scripts 03-05.")
}

swiss <- meta %>% filter(Country == "Switzerland") %>%
  mutate(Farm = factor(Farm, levels = unique(Farm[order(as.numeric(Farm))])),
         MAP_group = factor(MAP_group, levels = c("MAP-negative", "Subclinical", "Clinical")),
         MAP_binary = factor(MAP_binary, levels = c("MAP_negative", "MAP_positive")))

cat("\nSwiss animals per farm and MAP group:\n")
print(table(swiss$Farm, swiss$MAP_group))

# Farm-level labellings for B: every Swiss farm is entirely negative or positive.
farm_status <- swiss %>% distinct(Farm, MAP_binary)
if (anyDuplicated(farm_status$Farm)) {
  stop("A farm holds both MAP-negative and MAP-positive animals; ",
       "the farm-level permutation in B assumes it does not.")
}
farm_assignments <- combn(nrow(farm_status), sum(farm_status$MAP_binary == "MAP_negative"))
farm_index <- match(swiss$Farm, farm_status$Farm)
permuted_labels <- apply(farm_assignments, 2, function(neg_farms)
  ifelse(farm_index %in% neg_farms, "MAP_negative", "MAP_positive"))
cat("\nFarm-level permutation: ", nrow(farm_status), " farms; ",
    ncol(farm_assignments), " possible labellings; smallest attainable p = ",
    signif(1 / ncol(farm_assignments), 3), "\n", sep = "")

# Within-farm labellings for C: subclinical vs clinical on farms holding both.
within <- swiss %>% filter(MAP_group != "MAP-negative") %>%
  group_by(Farm) %>% filter(n_distinct(MAP_group) == 2) %>% ungroup() %>%
  mutate(Farm = droplevels(Farm), MAP_group = droplevels(MAP_group))
within_options <- lapply(split(seq_len(nrow(within)), within$Farm), function(rows) {
  n_clinical <- sum(within$MAP_group[rows] == "Clinical")
  combn(rows, n_clinical, simplify = FALSE)
})
within_grid <- expand.grid(lapply(within_options, seq_along))
within_labels <- apply(within_grid, 1, function(choice) {
  clinical_rows <- unlist(Map(function(opts, k) opts[[k]], within_options, choice))
  ifelse(seq_len(nrow(within)) %in% clinical_rows, "Clinical", "Subclinical")
})

# ------------------------------------------------------------- alpha values ----
# Same steps as read_profile() and calculate_alpha() in 03_Alpha_Diversity.R.
read_alpha <- function(rank, ids) {
  raw <- read_tsv(taxonomy_file(rank), show_col_types = FALSE)
  missing <- setdiff(ids, names(raw))
  if (length(missing)) stop(rank, " table lacks: ", paste(missing, collapse = ", "))
  raw %>%
    select(clade_name, all_of(ids)) %>%
    pivot_longer(-clade_name, names_to = "SAMPLE_ID", values_to = "Abundance") %>%
    mutate(Abundance = as.numeric(Abundance),
           Taxon = trimws(sub("^[dkpcofgst]__", "", sub("^.*[;|]", "", clade_name)))) %>%
    filter(!is.na(Taxon), nzchar(Taxon),
           !grepl("^(unclassified|unassigned|unknown|unidentified|na|none)$", tolower(Taxon)),
           !grepl("(^|[;|])t__", clade_name)) %>%
    group_by(SAMPLE_ID, Taxon) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(SAMPLE_ID) %>%
    summarise(Richness = sum(Abundance > 0),
              Shannon = { p <- Abundance[Abundance > 0] / sum(Abundance); -sum(p * log(p)) },
              Gini_Simpson = { p <- Abundance / sum(Abundance); 1 - sum(p^2) },
              Total = sum(Abundance), .groups = "drop") %>%
    mutate(Rank = rank)
}

alpha <- bind_rows(lapply(ranks, read_alpha, ids = swiss$SAMPLE_ID))
if (any(alpha$Total <= 0)) stop("Zero-abundance profile among Swiss samples.")
alpha_long <- alpha %>%
  select(-Total) %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Value") %>%
  left_join(swiss, by = "SAMPLE_ID")
write_csv(alpha_long, file.path(output_dir, "Alpha_values_Swiss.csv"))

# --------------------------------------------------------------- statistics ----
kw_h <- function(value, groups) {
  unname(kruskal.test(value, factor(groups))$statistic)
}

results <- list()
for (rank in ranks) for (metric in metrics) {
  d <- alpha_long %>% filter(Rank == rank, Metric == metric) %>% arrange(SAMPLE_ID)
  stopifnot(identical(d$SAMPLE_ID, swiss$SAMPLE_ID))

  reference <- kruskal.test(Value ~ MAP_group, data = d)

  # A. Farm as a fixed covariate.
  fit <- lm(Value ~ Farm + MAP_group, data = d)
  a_test <- drop1(fit, test = "F")["MAP_group", ]

  # B. Animal-level H, farm-level exact permutation.
  h_obs <- kw_h(d$Value, d$MAP_binary)
  h_perm <- apply(permuted_labels, 2, function(g) kw_h(d$Value, g))
  p_b <- mean(h_perm >= h_obs - 1e-12)

  # B2. Farm means, exact Wilcoxon rank-sum test (13 farms).
  farm_means <- d %>% group_by(Farm, MAP_binary) %>%
    summarise(Value = mean(Value), .groups = "drop")
  b2 <- suppressWarnings(wilcox.test(Value ~ MAP_binary, data = farm_means,
                                     exact = TRUE))

  # B3. Optional mixed model, random farm intercept, likelihood-ratio test.
  p_b3 <- NA_real_; est_b3 <- NA_real_
  if (use_lme4) {
    m1 <- lme4::lmer(Value ~ MAP_binary + (1 | Farm), data = d, REML = FALSE)
    m0 <- lme4::lmer(Value ~ 1 + (1 | Farm), data = d, REML = FALSE)
    p_b3 <- anova(m0, m1)$`Pr(>Chisq)`[2]
    est_b3 <- unname(lme4::fixef(m1)["MAP_binaryMAP_positive"])
  }

  # C. Subclinical vs clinical within mixed farms, within-farm exact permutation.
  dw <- d %>% filter(SAMPLE_ID %in% within$SAMPLE_ID) %>%
    arrange(match(SAMPLE_ID, within$SAMPLE_ID))
  h_within <- kw_h(dw$Value, within$MAP_group)
  h_within_perm <- apply(within_labels, 2, function(g) kw_h(dw$Value, g))
  p_c <- mean(h_within_perm >= h_within - 1e-12)

  medians <- d %>% group_by(MAP_binary) %>% summarise(m = median(Value))
  results[[paste(rank, metric)]] <- tibble(
    Rank = rank, Metric = metric,
    Median_MAP_negative = medians$m[medians$MAP_binary == "MAP_negative"],
    Median_MAP_positive = medians$m[medians$MAP_binary == "MAP_positive"],
    p_reference_KW_3groups = reference$p.value,
    A_MAP_group_Df = a_test$Df,
    p_A_lm_Farm_plus_MAP_group = a_test$`Pr(>F)`,
    p_unrestricted_KW_binary = kruskal.test(Value ~ MAP_binary, data = d)$p.value,
    p_B_farm_level_permutation = p_b,
    p_B2_farm_means_Wilcoxon = b2$p.value,
    B3_mixed_model_estimate_positive_minus_negative = est_b3,
    p_B3_mixed_model_LRT = p_b3,
    N_C_within_farm = nrow(dw),
    p_C_within_farm_permutation = p_c
  )
}

results <- bind_rows(results) %>%
  mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "BH"),
                .names = "{.col}_BH"))
write_csv(results, file.path(output_dir, "Alpha_farm_adjusted.csv"))

writeLines(c(
  "Farm-adjusted alpha-diversity sensitivity analyses (MetaPhlAn GTDB; Swiss n = 36).",
  "Alpha values calculated as in 03_Alpha_Diversity.R.",
  "Reference: Kruskal-Wallis across MAP-negative / subclinical / clinical (as in 03).",
  "A: lm(Value ~ Farm + MAP_group), F test for MAP_group (drop1). MAP-negative animals share",
  "   no farm with MAP-positive animals, so this term only reflects subclinical vs clinical",
  "   within mixed farms (Df = A_MAP_group_Df).",
  "B: MAP-negative vs MAP-positive, Kruskal-Wallis H on animals; MAP status permuted",
  paste0("   between whole farms, all ", ncol(farm_assignments),
         " labellings enumerated (smallest attainable p = ",
         signif(1 / ncol(farm_assignments), 3), ")."),
  "B2: same comparison on farm means (one value per farm), exact Wilcoxon rank-sum test.",
  if (use_lme4) "B3: lme4::lmer(Value ~ MAP_binary + (1 | Farm)), ML fit, likelihood-ratio test." else
    "B3: not run (lme4 not installed).",
  paste0("C: subclinical vs clinical within farms ", paste(levels(within$Farm), collapse = " & "),
         ", Kruskal-Wallis H, labels permuted within farm, all ", ncol(within_labels),
         " arrangements enumerated."),
  "*_BH columns: Benjamini-Hochberg across the 6 rank x metric tests within each analysis.",
  "South Africa: MAP-negative = farm 32, subclinical = farm 33; not separable, not tested.",
  paste("R", getRversion(), "| dplyr", packageVersion("dplyr"),
        if (use_lme4) paste("| lme4", packageVersion("lme4")) else "")
), file.path(output_dir, "README_farm_adjusted_alpha.txt"))

cat("\n==================== SUMMARY ====================\n")
print(as.data.frame(results %>% select(Rank, Metric, A_MAP_group_Df, !ends_with("_BH") & starts_with("p_"))),
      row.names = FALSE, digits = 3)
cat("\nOutputs written to: ", normalizePath(output_dir), "\n", sep = "")
