# 05_IndVal.R: indicator taxon analysis

# Map analysis methods to the original filenames in the Zenodo archive.
taxonomy_file <- function(method, rank, kind) {
  prefix <- ifelse(method == "MetaPhlAn_default", "MetaPhlAn", method)
  file.path("data", "taxonomy", method,
            paste0(prefix, "_", rank, "_", kind, ".tsv"))
}
# See README.md for inputs, provenance and validation status.

# Indicator taxa
local({
  input_dir <- file.path("data", "taxonomy")
  metadata_file <- file.path("data", "metadata.txt")
  output_dir <- "DVM_Indicator_Taxa_FINAL"

  n_permutations <- 9999
  minimum_prevalence_samples <- 2
  run_south_africa <- TRUE
  run_group_combinations <- TRUE
  set.seed(123)

  required_packages <- c(
    "readr", "dplyr", "tidyr", "stringr", "ggplot2",
    "indicspecies", "permute", "writexl"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Install missing packages in the console BEFORE running the script: ",
         paste(missing_packages, collapse = ", "),
         "\nFor example: install.packages(c(\"",
         paste(missing_packages, collapse = "\", \""), "\"))")
  }

  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ggplot2)
    library(indicspecies)
    library(permute)
    library(writexl)
  })

  if (!dir.exists(input_dir)) {
    stop("Missing standardized folder: ", input_dir,
         "\nExtract the Zenodo taxonomy archive to data/taxonomy/ first.")
  }
  if (!file.exists(metadata_file)) {
    stop("Metadata file missing: ", metadata_file)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sample_metadata <- dplyr::bind_rows(
    tibble(SAMPLE_ID = sprintf("SW_NEG_%02d", 2:17), Group = "SW_NEG"),
    tibble(SAMPLE_ID = sprintf("SW_SUB_%02d", 1:13), Group = "SW_SUB"),
    tibble(SAMPLE_ID = sprintf("SW_POS_%02d", 2:8),  Group = "SW_POS"),
    tibble(SAMPLE_ID = sprintf("ZA_NEG_%02d", 2:7), Group = "ZA_NEG"),
    tibble(SAMPLE_ID = sprintf("ZA_SUB_%02d", 1:8), Group = "ZA_SUB")
  ) %>%
    mutate(Country = if_else(startsWith(SAMPLE_ID, "SW_"),
                             "Switzerland", "South Africa"))

  metadata_raw <- readr::read_tsv(metadata_file, show_col_types = FALSE)
  if ("Sample_ID_Neu" %in% names(metadata_raw)) {
    metadata_ids <- trimws(as.character(metadata_raw$Sample_ID_Neu))
  } else if ("SAMPLE_ID" %in% names(metadata_raw)) {
    metadata_ids <- trimws(as.character(metadata_raw$SAMPLE_ID))
  } else {
    stop("Metadata must contain 'Sample_ID_Neu' or 'SAMPLE_ID'.")
  }
  missing_metadata <- setdiff(sample_metadata$SAMPLE_ID, metadata_ids)
  if (length(missing_metadata)) {
    stop("Expected samples not found in metadata: ",
         paste(missing_metadata, collapse = ", "))
  }
  stopifnot(!anyDuplicated(sample_metadata$SAMPLE_ID))

  expected_counts <- c(SW_NEG = 16, SW_SUB = 13, SW_POS = 7,
                       ZA_NEG = 6, ZA_SUB = 8)
  actual_counts <- table(factor(sample_metadata$Group,
                                levels = names(expected_counts)))
  stopifnot(identical(as.integer(actual_counts),
                      as.integer(expected_counts)))

  cat("\nCOHORT CHECK (expected 50 samples):\n")
  print(actual_counts)
  write_csv(sample_metadata, file.path(output_dir, "Included_50_samples.csv"))

  read_rank <- function(rank) {
    filename <- taxonomy_file("MetaPhlAn_GTDB", rank, "abundance")
    if (!file.exists(filename)) stop("Missing input: ", filename)

    raw <- readr::read_tsv(
      filename, show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    )
    stopifnot("clade_name" %in% names(raw))
    if (anyDuplicated(raw$clade_name)) stop("Duplicate taxonomic lineages")

    missing_samples <- setdiff(sample_metadata$SAMPLE_ID, names(raw))
    if (length(missing_samples)) {
      stop(rank, " input is missing: ",
           paste(missing_samples, collapse = ", "))
    }

    original_lineage <- raw$clade_name
    genus <- str_match(original_lineage, "(?:^|;)g__([^;]*)")[, 2]
    species <- str_match(original_lineage, "(?:^|;)s__([^;]*)")[, 2]
    taxon <- if (rank == "Genus") genus else species

    valid <- !is.na(taxon) & nzchar(trimws(taxon))
    taxonomy_qc <- tibble(
      Rank = rank,
      Input_lineages = length(valid),
      Unresolved_rank_excluded = sum(!valid),
      Named_lineages = sum(valid)
    )

    numeric_columns <- lapply(
      raw[sample_metadata$SAMPLE_ID],
      function(v) suppressWarnings(as.numeric(v))
    )
    x <- as.matrix(as.data.frame(numeric_columns, check.names = FALSE))
    rownames(x) <- original_lineage
    if (anyNA(x) || any(!is.finite(x)) || any(x < 0)) {
      stop(rank, " table contains NA, non-finite or negative abundances.")
    }
    x <- x[valid, , drop = FALSE]
    x <- t(x)  # sample x taxon
    stopifnot(identical(rownames(x), sample_metadata$SAMPLE_ID))

    lookup <- tibble(
      taxon_id = original_lineage[valid],
      taxon = trimws(taxon[valid]),
      genus = genus[valid],
      species_GTDB = species[valid],
      Rank = rank
    )
    list(matrix = x, lookup = lookup, qc = taxonomy_qc)
  }

  prepare_cohort <- function(taxon_data, country) {
    meta <- sample_metadata %>% filter(Country == country)
    x <- taxon_data$matrix[meta$SAMPLE_ID, , drop = FALSE]

    detected_n <- colSums(x > 0)
    keep <- detected_n >= minimum_prevalence_samples &
            apply(x, 2, function(v) length(unique(v)) > 1)
    x <- x[, keep, drop = FALSE]
    if (ncol(x) == 0) stop("No taxa after prevalence filter: ", country)
    stopifnot(identical(rownames(x), meta$SAMPLE_ID))

    list(
      x = x,
      meta = meta,
      lookup = taxon_data$lookup %>% filter(taxon_id %in% colnames(x)),
      qc = tibble(
        Rank = unique(taxon_data$lookup$Rank),
        Country = country,
        Samples = nrow(x),
        Taxa_after_prevalence_filter = ncol(x),
        Minimum_detected_samples = minimum_prevalence_samples
      )
    )
  }

  make_descriptives <- function(x, groups) {
    result <- lapply(levels(groups), function(g) {
      sub <- x[groups == g, , drop = FALSE]
      tibble(
        taxon_id = colnames(x),
        Group = g,
        n_group = nrow(sub),
        prevalence_n = colSums(sub > 0),
        prevalence_percent = 100 * colMeans(sub > 0),
        mean_abundance = colMeans(sub),
        median_abundance = apply(sub, 2, median)
      )
    })
    bind_rows(result) %>%
      tidyr::pivot_wider(
        names_from = Group,
        values_from = c(n_group, prevalence_n, prevalence_percent,
                        mean_abundance, median_abundance),
        names_glue = "{.value}_{Group}"
      )
  }

  run_indval <- function(prepped, rank, analysis, group_vector,
                         combinations = FALSE) {
    groups <- droplevels(factor(group_vector))
    stopifnot(length(groups) == nrow(prepped$x), !anyNA(groups))
    if (nlevels(groups) < 2) stop("At least two groups required.")

    cat("\n--- ", rank, " | ", analysis, " ---\n", sep = "")
    print(table(groups))

    set.seed(123)
    fit <- indicspecies::multipatt(
      prepped$x,
      groups,
      func = "IndVal.g",
      duleg = !combinations,
      min.order = 1,
      max.order = if (combinations) 2 else 1,
      control = permute::how(nperm = n_permutations)
    )

    sign <- as.data.frame(fit$sign, check.names = FALSE)
    flags <- grep("^s\\.", names(sign), value = TRUE)
    if (length(flags) != nlevels(groups)) {
      stop("Unexpected group columns in multipatt result: ",
           paste(names(sign), collapse = ", "))
    }
    group_names <- levels(groups)
    assignments <- apply(sign[, flags, drop = FALSE], 1, function(v) {
      ix <- which(as.numeric(v) == 1)
      if (!length(ix)) return(NA_character_)
      paste(group_names[ix], collapse = " + ")
    })

    results <- tibble(
      taxon_id = rownames(sign),
      indicator_group = assignments,
      IndVal_g_sqrt = as.numeric(sign$stat),
      p_value = as.numeric(sign$p.value)
    ) %>%
      mutate(p_adj_BH = p.adjust(p_value, method = "BH"),
             significant_BH = !is.na(p_adj_BH) & p_adj_BH < 0.05,
             Analysis = analysis, Rank = rank) %>%
      left_join(prepped$lookup, by = c("taxon_id", "Rank")) %>%
      left_join(make_descriptives(prepped$x, groups), by = "taxon_id") %>%
      select(Analysis, Rank, taxon, genus, species_GTDB, taxon_id,
             indicator_group, IndVal_g_sqrt, p_value, p_adj_BH,
             significant_BH, everything()) %>%
      arrange(p_adj_BH, p_value)

    prefix <- paste(rank, analysis, sep = "_")
    write_csv(results, file.path(output_dir, paste0(prefix, "_ALL.csv")))
    write_csv(results %>% filter(p_value < 0.05),
              file.path(output_dir, paste0(prefix, "_raw_p05.csv")))
    write_csv(results %>% filter(significant_BH),
              file.path(output_dir, paste0(prefix, "_BH_FDR05.csv")))

    cat("Taxa tested: ", nrow(results),
        " | raw p<0.05: ", sum(results$p_value < 0.05, na.rm = TRUE),
        " | BH-FDR<0.05: ", sum(results$significant_BH), "\n", sep = "")
    if (any(results$significant_BH)) {
      print(results %>% filter(significant_BH) %>%
              select(taxon, indicator_group, IndVal_g_sqrt,
                     p_value, p_adj_BH) %>% head(20), n = 20)
    }
    list(
      results = results,
      summary = tibble(
        Analysis = analysis, Rank = rank,
        n_samples = nrow(prepped$x),
        group_sizes = paste(names(table(groups)), as.integer(table(groups)),
                            collapse = "; "),
        n_taxa_tested = nrow(results),
        raw_p05 = sum(results$p_value < 0.05, na.rm = TRUE),
        BH_FDR05 = sum(results$significant_BH),
        permutations = n_permutations
      )
    )
  }

  all_results <- list()
  summary_list <- list()
  input_qc_list <- list()
  cohort_qc_list <- list()

  for (rank in c("Genus", "Species")) {
    td <- read_rank(rank)
    input_qc_list[[rank]] <- td$qc

    sw <- prepare_cohort(td, "Switzerland")
    za <- prepare_cohort(td, "South Africa")
    cohort_qc_list[[paste0(rank, "_SW")]] <- sw$qc
    cohort_qc_list[[paste0(rank, "_ZA")]] <- za$qc

    sw3 <- factor(sw$meta$Group, levels = c("SW_NEG", "SW_SUB", "SW_POS"))
    a <- run_indval(sw, rank, "SW_three_groups", sw3)
    all_results[[paste0(rank, "_SW3")]] <- a$results
    summary_list[[paste0(rank, "_SW3")]] <- a$summary

    sw2 <- factor(ifelse(sw$meta$Group == "SW_NEG", "MAP_negative", "MAP_positive"),
                  levels = c("MAP_negative", "MAP_positive"))
    b <- run_indval(sw, rank, "SW_binary_MAP", sw2)
    all_results[[paste0(rank, "_SWbinary")]] <- b$results
    summary_list[[paste0(rank, "_SWbinary")]] <- b$summary

    if (run_group_combinations) {
      c <- run_indval(sw, rank, "SW_group_combinations", sw3,
                      combinations = TRUE)
      all_results[[paste0(rank, "_SWpairs")]] <- c$results
      summary_list[[paste0(rank, "_SWpairs")]] <- c$summary
    }

    if (run_south_africa) {
      za2 <- factor(za$meta$Group, levels = c("ZA_NEG", "ZA_SUB"))
      d <- run_indval(za, rank, "ZA_binary_MAP", za2)
      all_results[[paste0(rank, "_ZA")]] <- d$results
      summary_list[[paste0(rank, "_ZA")]] <- d$summary
    }
  }

  summary_table <- bind_rows(summary_list)
  input_qc <- bind_rows(input_qc_list)
  cohort_qc <- bind_rows(cohort_qc_list)

  write_csv(summary_table, file.path(output_dir, "Analysis_overview.csv"))
  write_csv(input_qc, file.path(output_dir, "Taxonomy_input_QC.csv"))
  write_csv(cohort_qc, file.path(output_dir, "Cohort_prevalence_QC.csv"))

  notes <- tibble(
    Topic = c("Cohort", "Swiss primary", "Swiss secondary", "Group combinations",
              "South Africa", "Taxonomic levels", "Abundances", "Prevalence",
              "IndVal statistic", "Permutations", "BH correction",
              "Study-design limitation"),
    Explanation = c(
      "Exact alpha/beta microbiome cohort: Swiss n=36; South African n=14. Additional diagnostic MAP-detection samples excluded.",
      "SW_NEG (16), SW_SUB (13), SW_POS (7); single-group indicators.",
      "SW_NEG (16) vs SW_SUB+SW_POS (20); single-group indicators.",
      "Up to TWO of the Swiss groups may form an indicator combination; exploratory.",
      "ZA_NEG (6) versus ZA_SUB (8); separate exploratory analysis.",
      "MetaPhlAn GTDB standardized Genus_abundance.tsv and Species_abundance.tsv.",
      "Uses original standardized relative abundance values, not renormalized or resampled.",
      paste("A taxon must occur (>0) in at least", minimum_prevalence_samples,
            "samples within the corresponding cohort."),
      "indicspecies::multipatt(func='IndVal.g') reports the SQUARE ROOT of IndVal.g; higher is stronger.",
      paste(n_permutations, "permutations; random seed 123."),
      "Benjamini-Hochberg correction separately within each rank x analysis across all tested taxa; significance at adjusted p < 0.05.",
      "MAP status and farm origin are confounded in the Swiss cohort. Indicator associations are descriptive, not causal or independently MAP-specific."
    )
  )

  sheets <- c(
    list(Overview = summary_table,
         Sample_groups = sample_metadata,
         Taxonomy_QC = input_qc,
         Cohort_QC = cohort_qc,
         Read_me = notes),
    all_results
  )

  writexl::write_xlsx(
    sheets, path = file.path(output_dir, "Indicator_Taxa_Supplementary_Master.xlsx")
  )

  plot_top_indicators <- function(results, rank) {
    top <- results %>%
      filter(!is.na(indicator_group), !is.na(IndVal_g_sqrt)) %>%
      arrange(p_adj_BH, desc(IndVal_g_sqrt)) %>%
      slice_head(n = 20) %>%
      mutate(
        plot_label = make.unique(as.character(taxon)),
        plot_label = factor(plot_label, levels = rev(plot_label)),
        Status = if_else(significant_BH, "BH-FDR < 0.05", "Not BH-significant")
      )

    if (!nrow(top)) return(invisible(NULL))

    p <- ggplot(top, aes(x = plot_label, y = IndVal_g_sqrt,
                         fill = indicator_group)) +
      geom_col(width = 0.75, color = "grey35", linewidth = 0.2) +
      geom_point(data = top %>% filter(significant_BH),
                 aes(y = IndVal_g_sqrt), inherit.aes = TRUE,
                 shape = 21, size = 1.9, fill = "white", color = "black") +
      coord_flip() +
      scale_fill_manual(values = c(SW_NEG = "#49AFC7",
                                   SW_SUB = "#9479BA",
                                   SW_POS = "#EB907C")) +
      labs(x = NULL, y = "sqrt(IndVal.g)",
           fill = "Associated group",
           title = paste(rank, "— Swiss three-group indicator taxa"),
           subtitle = "Top 20 by BH-adjusted p-value; white dots indicate BH-FDR < 0.05") +
      theme_classic(base_size = 11) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold"))

    ggsave(file.path(output_dir, paste0("SW_three_groups_", rank, "_top20.png")),
           p, width = 9, height = 7, dpi = 400, bg = "white")
    ggsave(file.path(output_dir, paste0("SW_three_groups_", rank, "_top20.pdf")),
           p, width = 9, height = 7, bg = "white")
  }

  for (rank in c("Genus", "Species")) {
    plot_top_indicators(all_results[[paste0(rank, "_SW3")]], rank)
  }

  cat("\n================================================\n")
  cat("INDICATOR TAXA ANALYSIS COMPLETE\n")
  cat("================================================\n")
  print(summary_table, n = nrow(summary_table))
  cat("\nResults saved to:\n", normalizePath(output_dir), "\n")
  cat("\nIMPORTANT: Report BH-adjusted p-values. Associations\n")
  cat("are not independent of farm origin and do not establish causality.\n")
})
