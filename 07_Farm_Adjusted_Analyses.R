# 07_Farm_Adjusted_Analyses.R
# Farm-adjusted sensitivity analyses for beta diversity (04) and indicator taxa (05).
#
# Usage in RStudio: set the working directory to the folder that CONTAINS the
# extracted `data/` folder (Session > Set Working Directory > Choose Directory),
# or edit `data_dir` below, then run the whole script (Ctrl/Cmd + Shift + Enter).
#
# Design note (why this is not only "Farm + MAP_group"):
#   Swiss MAP-negative animals come from farms 1-5 only, subclinical animals from
#   farms 6-8, and clinical animals from farms 6 and 7 (one each) or from
#   single-animal farms. MAP-negative vs MAP-positive is therefore completely
#   nested in farm. In `adonis2(bray ~ Farm + MAP_group, by = "margin")` the
#   MAP_group term can only use the subclinical-vs-clinical contrast inside
#   farms 6 and 7 (1 df). It cannot test negative vs positive at all.
#   This script therefore runs:
#     A. the requested covariate model (reported, with its effective df);
#     B. a farm-level exact permutation test of MAP-negative vs MAP-positive,
#        treating the farm as the unit of replication (13 Swiss farms);
#     C. a within-farm test of subclinical vs clinical (farms 6 and 7);
#     D. indicator taxa (IndVal.g) for MAP-negative vs MAP-positive with the
#        same farm-level exact permutation.
#   South Africa: all MAP-negative animals are on farm 32 and all subclinical
#   animals on farm 33, so farm and MAP status cannot be separated there.

# ---------------------------------------------------------------- settings ----
data_dir   <- "data"
output_dir <- "Farm_adjusted_analyses"
ranks      <- c("Genus", "Species")
n_perm     <- 9999            # for the adonis2 models (A, C)
minimum_prevalence_samples <- 2  # as in 05_IndVal.R
seed       <- 12345

# ---------------------------------------------------------------- packages ----
packages <- c("vegan", "indicspecies", "permute", "dplyr", "readr", "tidyr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1),
                                     quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install first: install.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))")
}
suppressPackageStartupMessages({
  library(vegan); library(indicspecies); library(permute)
  library(dplyr); library(readr); library(tidyr)
})
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

taxonomy_file <- function(rank, kind) {
  file.path(data_dir, "taxonomy", "MetaPhlAn_GTDB",
            paste0("MetaPhlAn_GTDB_", rank, "_", kind, ".tsv"))
}
for (f in c(file.path(data_dir, "metadata.txt"),
            taxonomy_file(ranks, "normalized"), taxonomy_file(ranks, "abundance"))) {
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

farm_table <- as.data.frame.matrix(table(swiss$Farm, swiss$MAP_group))
farm_table <- cbind(Farm = rownames(farm_table), farm_table)
cat("\nSwiss animals per farm and MAP group:\n"); print(farm_table, row.names = FALSE)
cat("\nSouth African animals per farm and MAP group:\n")
print(table(meta$Farm[meta$Country == "South Africa"],
            meta$MAP_group[meta$Country == "South Africa"]))
write_csv(farm_table, file.path(output_dir, "Swiss_farm_by_MAP_group.csv"))

# Farm-level labels: every Swiss farm is entirely MAP-negative or MAP-positive.
farm_status <- swiss %>% distinct(Farm, MAP_binary)
if (anyDuplicated(farm_status$Farm)) {
  stop("A farm holds both MAP-negative and MAP-positive animals; ",
       "the farm-level permutation in B/D assumes it does not.")
}
n_farms     <- nrow(farm_status)
n_neg_farms <- sum(farm_status$MAP_binary == "MAP_negative")
# All ways of labelling n_neg_farms of the n_farms as MAP-negative (exact test).
farm_assignments <- combn(n_farms, n_neg_farms)
farm_index <- match(swiss$Farm, farm_status$Farm)
observed_column <- which(apply(farm_assignments, 2, function(ix)
  setequal(ix, which(farm_status$MAP_binary == "MAP_negative"))))
stopifnot(length(observed_column) == 1)
permuted_labels <- apply(farm_assignments, 2, function(neg_farms)
  ifelse(farm_index %in% neg_farms, "MAP_negative", "MAP_positive"))
cat("\nFarm-level permutation: ", n_farms, " farms (", n_neg_farms,
    " MAP-negative); ", ncol(farm_assignments),
    " possible labellings; smallest attainable p = ",
    signif(1 / ncol(farm_assignments), 3), "\n", sep = "")

# ---------------------------------------------------------------- profiles ----
# Beta diversity input as in 04: normalized table, unclassified rows removed,
# renormalized per sample.
read_normalized <- function(rank, ids) {
  raw <- read_tsv(taxonomy_file(rank, "normalized"), show_col_types = FALSE,
                  name_repair = "minimal")
  terminal <- sub("^[a-z]__", "", basename(gsub("[|;]", "/", raw$clade_name)))
  raw <- raw[!tolower(trimws(terminal)) %in% c("unclassified", "unassigned", "unknown"), ]
  missing <- setdiff(ids, names(raw))
  if (length(missing)) stop(rank, " normalized table lacks: ", paste(missing, collapse = ", "))
  x <- t(as.matrix(raw[, ids]))
  colnames(x) <- make.unique(raw$clade_name)
  storage.mode(x) <- "double"
  x <- x[, colSums(x) > 0, drop = FALSE]
  if (any(rowSums(x) <= 0)) stop("Empty profile in ", rank)
  sweep(x, 1, rowSums(x), "/")
}

# Indicator input as in 05: original abundance table, named taxa at this rank,
# prevalence filter within the cohort.
read_abundance <- function(rank, ids) {
  raw <- read_tsv(taxonomy_file(rank, "abundance"), show_col_types = FALSE,
                  col_types = cols(.default = col_character()))
  pattern <- if (rank == "Genus") "(?:^|;)g__([^;]*)" else "(?:^|;)s__([^;]*)"
  taxon <- regmatches(raw$clade_name, regexec(pattern, raw$clade_name, perl = TRUE))
  taxon <- vapply(taxon, function(m) if (length(m) == 2) m[2] else NA_character_, "")
  keep <- !is.na(taxon) & nzchar(trimws(taxon))
  x <- t(sapply(raw[keep, ids], as.numeric))
  colnames(x) <- raw$clade_name[keep]
  rownames(x) <- ids
  if (anyNA(x) || any(x < 0)) stop("Invalid abundances in ", rank)
  prevalent <- colSums(x > 0) >= minimum_prevalence_samples &
    apply(x, 2, function(v) length(unique(v)) > 1)
  list(x = x[, prevalent, drop = FALSE],
       taxon = setNames(trimws(taxon[keep][prevalent]), colnames(x)[prevalent]))
}

# Pseudo-F (as in adonis2) for a one-factor design, from squared distances.
pseudo_f <- function(d2, groups) {
  n <- length(groups); lev <- unique(groups); g <- length(lev)
  ss_total <- sum(d2) / (2 * n)
  ss_within <- sum(vapply(lev, function(l) {
    i <- groups == l; sum(d2[i, i]) / (2 * sum(i))
  }, numeric(1)))
  c(F = ((ss_total - ss_within) / (g - 1)) / (ss_within / (n - g)),
    R2 = (ss_total - ss_within) / ss_total)
}

# IndVal.g (square root, as reported by multipatt) for each taxon and group;
# returns the maximum over groups and the group attaining it.
indval_g <- function(x, groups) {
  lev <- sort(unique(groups))
  means <- sapply(lev, function(l) colMeans(x[groups == l, , drop = FALSE]))
  freq  <- sapply(lev, function(l) colMeans(x[groups == l, , drop = FALSE] > 0))
  a <- means / rowSums(means)
  iv <- sqrt(a * freq)
  list(stat = apply(iv, 1, max), group = lev[apply(iv, 1, which.max)])
}

# ------------------------------------------------------------------ analyses ----
beta_results <- list()
indval_results <- list()
indval_summary <- list()

for (rank in ranks) {
  cat("\n==================== ", rank, " ====================\n", sep = "")
  x <- read_normalized(rank, swiss$SAMPLE_ID)
  bray <- vegdist(x, method = "bray")
  d2 <- as.matrix(bray)^2

  # Reference: unadjusted MAP-group PERMANOVA, as in 04 (Swiss_MAP).
  set.seed(seed)
  unadjusted <- adonis2(bray ~ MAP_group, data = swiss, permutations = n_perm)

  # A. Requested covariate model.
  set.seed(seed)
  covariate <- adonis2(bray ~ Farm + MAP_group, data = swiss,
                       permutations = n_perm, by = "margin")
  cat("\nA. adonis2(bray ~ Farm + MAP_group, by = 'margin'):\n"); print(covariate)

  # B. Farm-level exact permutation, MAP-negative vs MAP-positive.
  observed <- pseudo_f(d2, as.character(swiss$MAP_binary))
  set.seed(seed)
  check <- adonis2(bray ~ MAP_binary, data = swiss, permutations = 0)
  stopifnot(isTRUE(all.equal(unname(observed["F"]), check$F[1], tolerance = 1e-8)))
  perm_f <- apply(permuted_labels, 2, function(g) pseudo_f(d2, g)["F"])
  p_farm <- mean(perm_f >= observed["F"] - 1e-12)
  set.seed(seed)
  binary_unrestricted <- adonis2(bray ~ MAP_binary, data = swiss, permutations = n_perm)

  # C. Subclinical vs clinical within farms 6 and 7 (permutations within farm).
  within <- swiss %>% filter(MAP_group != "MAP-negative") %>%
    group_by(Farm) %>% filter(n_distinct(MAP_group) == 2) %>% ungroup() %>%
    mutate(Farm = droplevels(Farm), MAP_group = droplevels(MAP_group))
  bray_within <- vegdist(x[within$SAMPLE_ID, , drop = FALSE], method = "bray")
  set.seed(seed)
  within_test <- adonis2(bray_within ~ MAP_group, data = within,
                         permutations = how(blocks = within$Farm, nperm = n_perm))

  beta_results[[rank]] <- bind_rows(
    tibble(Analysis = "Reference (04): MAP_group, 3 groups, unrestricted permutations",
           N = nrow(swiss), Df = unadjusted$Df[1], R2 = unadjusted$R2[1],
           F = unadjusted$F[1], p = unadjusted$`Pr(>F)`[1]),
    tibble(Analysis = "A: Farm + MAP_group (margin); MAP_group term",
           N = nrow(swiss), Df = covariate["MAP_group", "Df"],
           R2 = covariate["MAP_group", "R2"], F = covariate["MAP_group", "F"],
           p = covariate["MAP_group", "Pr(>F)"]),
    tibble(Analysis = "A: Farm + MAP_group (margin); Farm term",
           N = nrow(swiss), Df = covariate["Farm", "Df"],
           R2 = covariate["Farm", "R2"], F = covariate["Farm", "F"],
           p = covariate["Farm", "Pr(>F)"]),
    tibble(Analysis = "MAP-negative vs MAP-positive, unrestricted permutations",
           N = nrow(swiss), Df = binary_unrestricted$Df[1],
           R2 = binary_unrestricted$R2[1], F = binary_unrestricted$F[1],
           p = binary_unrestricted$`Pr(>F)`[1]),
    tibble(Analysis = "B: MAP-negative vs MAP-positive, farm-level exact permutation",
           N = nrow(swiss), Df = 1, R2 = unname(observed["R2"]),
           F = unname(observed["F"]), p = p_farm),
    tibble(Analysis = paste0("C: Subclinical vs clinical within farms ",
                             paste(levels(within$Farm), collapse = " & "),
                             ", permutations within farm"),
           N = nrow(within), Df = within_test$Df[1], R2 = within_test$R2[1],
           F = within_test$F[1], p = within_test$`Pr(>F)`[1])
  ) %>% mutate(Rank = rank, .before = 1)
  cat("\nBeta-diversity summary:\n")
  print(as.data.frame(beta_results[[rank]]), row.names = FALSE, digits = 4)

  # D. Indicator taxa, MAP-negative vs MAP-positive, farm-level exact permutation.
  ab <- read_abundance(rank, swiss$SAMPLE_ID)
  groups <- as.character(swiss$MAP_binary)
  obs <- indval_g(ab$x, groups)
  set.seed(123)
  mp <- multipatt(ab$x, swiss$MAP_binary, func = "IndVal.g", duleg = TRUE,
                  control = how(nperm = n_perm))
  stopifnot(isTRUE(all.equal(unname(obs$stat), mp$sign[colnames(ab$x), "stat"],
                             tolerance = 1e-8)))
  perm_stats <- apply(permuted_labels, 2, function(g) indval_g(ab$x, g)$stat)
  p_taxa <- rowMeans(perm_stats >= obs$stat - 1e-12)

  indval_results[[rank]] <- tibble(
    Rank = rank, taxon = ab$taxon[colnames(ab$x)], taxon_id = colnames(ab$x),
    indicator_group = obs$group, IndVal_g_sqrt = obs$stat,
    p_farm_level = p_taxa,
    p_adj_BH_farm_level = p.adjust(p_taxa, method = "BH"),
    p_unrestricted_multipatt = mp$sign[colnames(ab$x), "p.value"],
    p_adj_BH_unrestricted = p.adjust(mp$sign[colnames(ab$x), "p.value"], method = "BH")
  ) %>% arrange(p_farm_level, desc(IndVal_g_sqrt))

  indval_summary[[rank]] <- indval_results[[rank]] %>% summarise(
    Rank = rank, Taxa_tested = n(),
    Unrestricted_raw_p05 = sum(p_unrestricted_multipatt < 0.05),
    Unrestricted_BH05 = sum(p_adj_BH_unrestricted < 0.05),
    Farm_level_raw_p05 = sum(p_farm_level < 0.05),
    Farm_level_BH05 = sum(p_adj_BH_farm_level < 0.05),
    Smallest_attainable_p = 1 / ncol(farm_assignments))
  cat("\nD. Indicator taxa (MAP-negative vs MAP-positive):\n")
  print(as.data.frame(indval_summary[[rank]]), row.names = FALSE)
}

# ------------------------------------------------------------------- outputs ----
beta_table   <- bind_rows(beta_results)
indval_table <- bind_rows(indval_results)
indval_overview <- bind_rows(indval_summary)
write_csv(beta_table, file.path(output_dir, "Beta_farm_adjusted.csv"))
write_csv(indval_table, file.path(output_dir, "IndVal_binary_farm_level.csv"))
write_csv(indval_overview, file.path(output_dir, "IndVal_binary_farm_level_overview.csv"))
writeLines(c(
  "Farm-adjusted sensitivity analyses (MetaPhlAn GTDB; Swiss n = 36).",
  "A: adonis2(bray ~ Farm + MAP_group, by = 'margin'). MAP-negative animals share no farm",
  "   with MAP-positive animals, so the MAP_group term only reflects subclinical vs",
  "   clinical within mixed farms; check its Df.",
  "B: MAP-negative vs MAP-positive with MAP status permuted between whole farms",
  "   (all labellings enumerated; farm = unit of replication).",
  "C: subclinical vs clinical within farms holding both, permutations restricted within farm.",
  "D: IndVal.g (sqrt), duleg, MAP-negative vs MAP-positive; p from the same farm-level",
  "   permutation as B; BH within rank. Unrestricted multipatt p shown for comparison.",
  "South Africa: MAP-negative = farm 32, subclinical = farm 33; not separable, not tested.",
  paste("R", getRversion(), "| vegan", packageVersion("vegan"),
        "| indicspecies", packageVersion("indicspecies"),
        "| permute", packageVersion("permute"))
), file.path(output_dir, "README_farm_adjusted.txt"))

cat("\n==================== SUMMARY ====================\n")
print(as.data.frame(beta_table), row.names = FALSE, digits = 4)
cat("\n"); print(as.data.frame(indval_overview), row.names = FALSE)
cat("\nOutputs written to: ", normalizePath(output_dir), "\n", sep = "")
