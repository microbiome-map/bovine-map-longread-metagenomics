# 02_Taxonomic_Analysis.R: taxonomic composition, differential abundance and classifier profiles

# Map analysis methods to the original filenames in the Zenodo archive.
taxonomy_file <- function(method, rank, kind) {
  prefix <- ifelse(method == "MetaPhlAn_default", "MetaPhlAn", method)
  file.path("data", "taxonomy", method,
            paste0(prefix, "_", rank, "_", kind, ".tsv"))
}
# See README.md for inputs, provenance and validation status.


local({
  library(tidyverse)
  library(writexl)

  input_dir <- file.path("data", "taxonomy")
  meta_file <- file.path("data", "metadata.txt")

  output_dir <- "DVM_Taxonomic_Composition_and_Differential_Abundance"

  methods <- c("MetaPhlAn_GTDB", "ganon2")
  ranks <- c("Phylum", "Genus", "Species")

  group_order <- c(
    "SW_NEG", "SW_SUB", "SW_POS",
    "ZA_NEG", "ZA_SUB"
  )

  excluded_ids <- c("ZA_NEG_01", "ZA_NEG_08")
  alpha <- 0.05

  original_units <- c(
    MetaPhlAn_GTDB = "percent",
    ganon2 = "fraction"
  )

  top_n <- c(Phylum = 10L, Genus = 20L)

  cohorts <- list(
    Switzerland = c("SW_NEG", "SW_SUB", "SW_POS"),
    South_Africa = c("ZA_NEG", "ZA_SUB"),
    Combined = group_order
  )

  secondary_comparisons <- list(
    South_Africa = c("ZA_NEG", "ZA_SUB"),
    Country_NEG = c("SW_NEG", "ZA_NEG"),
    Country_SUB = c("SW_SUB", "ZA_SUB")
  )

  folders <- c(
    "Tables",
    "Statistics",
    "QC",
    "Data_used",
    "Figures",
    "Reproducibility"
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (folder in folders) {
    dir.create(
      file.path(output_dir, folder),
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  meta_raw <- read_tsv(
    meta_file,
    show_col_types = FALSE,
    name_repair = "check_unique"
  )

  if (!all(c("SAMPLE_ID", "GROUP") %in% names(meta_raw))) {
    stop("Metadata must contain SAMPLE_ID and GROUP.")
  }

  meta_all <- meta_raw %>%
    mutate(
      SAMPLE_ID = trimws(as.character(SAMPLE_ID)),
      GROUP = trimws(as.character(GROUP)),

      Cohort_group = case_when(
        str_starts(SAMPLE_ID, "SW_NEG_") &
          GROUP == "Healthy" ~ "SW_NEG",

        str_starts(SAMPLE_ID, "SW_SUB_") &
          GROUP == "Positive, Subclinical" ~ "SW_SUB",

        str_starts(SAMPLE_ID, "SW_POS_") &
          GROUP == "Positive, Clinical" ~ "SW_POS",

        str_starts(SAMPLE_ID, "ZA_NEG_") &
          GROUP == "ZA Healthy" ~ "ZA_NEG",

        str_starts(SAMPLE_ID, "ZA_SUB_") &
          GROUP == "ZA Positive, Subclinical" ~ "ZA_SUB",

        TRUE ~ NA_character_
      ),

      Exclusion_reason = case_when(
        SAMPLE_ID %in% excluded_ids ~
          "Insufficient sequencing depth",

        GROUP == "MAP detection" ~
          "MAP detection only",

        GROUP == "EX" ~
          "Excluded in metadata",

        is.na(Cohort_group) ~
          "Not an eligible microbiome group",

        TRUE ~ NA_character_
      )
    )

  if (anyNA(meta_all$SAMPLE_ID) ||
      any(!nzchar(meta_all$SAMPLE_ID))) {
    stop("Missing or blank sample IDs in metadata.")
  }

  if (anyDuplicated(meta_all$SAMPLE_ID)) {
    stop("Duplicate sample IDs in metadata.")
  }

  files <- expand_grid(
    Method = methods,
    Rank = ranks,
    Kind = c("normalized", "abundance")
  ) %>%
    mutate(
      Path = taxonomy_file(Method, Rank, Kind)
    )

  if (any(!file.exists(files$Path))) {
    stop(
      "Missing input files:\n",
      paste(
        files$Path[!file.exists(files$Path)],
        collapse = "\n"
      )
    )
  }

  get_sample_ids <- function(path) {
    header <- names(
      read_tsv(
        path,
        n_max = 0,
        show_col_types = FALSE,
        name_repair = "check_unique"
      )
    )

    if (!"clade_name" %in% header) {
      stop("Missing clade_name: ", path)
    }

    setdiff(header, "clade_name")
  }

  shared_ids <- Reduce(
    intersect,
    lapply(files$Path, get_sample_ids)
  )

  meta_all <- meta_all %>%
    mutate(
      Exclusion_reason = case_when(
        !is.na(Exclusion_reason) ~ Exclusion_reason,
        !SAMPLE_ID %in% shared_ids ~
          "Missing from at least one taxonomic profile",
        TRUE ~ NA_character_
      )
    )

  meta <- meta_all %>%
    filter(is.na(Exclusion_reason)) %>%
    mutate(
      Cohort_group = factor(
        Cohort_group,
        levels = group_order
      ),
      Country = if_else(
        str_starts(SAMPLE_ID, "SW_"),
        "Switzerland",
        "South Africa"
      )
    ) %>%
    arrange(Cohort_group, SAMPLE_ID)

  if (!nrow(meta)) {
    stop("No eligible samples remain.")
  }

  sample_counts <- meta %>%
    count(Cohort_group, .drop = FALSE)

  if (any(sample_counts$n < 2)) {
    print(sample_counts)
    stop("Every analysis group must contain at least two samples.")
  }

  cat("\nFINAL SAMPLE COUNTS:\n")
  print(sample_counts)

  write_csv(
    meta_all %>%
      select(SAMPLE_ID, GROUP, Cohort_group, Exclusion_reason),
    file.path(output_dir, "QC", "Sample_inclusion_exclusion.csv")
  )

  write_csv(
    meta,
    file.path(output_dir, "Data_used", "Metadata_used.csv")
  )

  clean_taxon <- function(x) {
    x <- sub("^.*[;|]", "", as.character(x))
    x <- sub("^[dkpcofgst]__", "", x)
    trimws(x)
  }

  harmonize_phylum <- function(x) {
    x <- if_else(
      str_detect(x, "^Bacillota(_[A-Z0-9]+)?$|^Firmicutes$"),
      "Bacillota",
      x
    )

    recode(
      x,
      "Proteobacteria" = "Pseudomonadota",
      "Bacteroidetes" = "Bacteroidota",
      "Actinobacteria" = "Actinomycetota",
      .default = x
    )
  }

  read_profile <- function(method, rank, kind) {

    path <- taxonomy_file(method, rank, kind)

    dat <- read_tsv(
      path,
      show_col_types = FALSE,
      name_repair = "check_unique"
    )

    if (!"clade_name" %in% names(dat)) {
      stop("Missing clade_name: ", path)
    }

    dat <- dat %>%
      select(clade_name, all_of(meta$SAMPLE_ID))

    values <- suppressWarnings(
      as.numeric(
        as.matrix(dat %>% select(-clade_name))
      )
    )

    if (any(!is.finite(values) | values < 0)) {
      stop(
        "Missing, non-numeric, infinite or negative abundance ",
        "values in included samples: ", path
      )
    }

    dat %>%
      pivot_longer(
        -clade_name,
        names_to = "SAMPLE_ID",
        values_to = "Abundance",
        values_transform = list(Abundance = as.numeric)
      ) %>%
      mutate(
        Taxon = clean_taxon(clade_name),
        Taxon = if (rank == "Phylum") {
          harmonize_phylum(Taxon)
        } else {
          Taxon
        }
      ) %>%
      filter(
        !is.na(Taxon),
        nzchar(Taxon),
        !str_detect(
          str_to_lower(Taxon),
          "^(unclassified|unassigned|unknown|unidentified|na|none)$"
        ),
        !str_detect(clade_name, "(^|[;|])t__")
      ) %>%
      group_by(SAMPLE_ID, Taxon) %>%
      summarise(
        Abundance = sum(Abundance),
        .groups = "drop"
      )
  }

  profiles_list <- list()
  qc_list <- list()

  for (method in methods) {
    for (rank in ranks) {

      key <- paste(method, rank, sep = "__")

      norm <- read_profile(method, rank, "normalized") %>%
        rename(Relative_abundance_input = Abundance)

      raw <- read_profile(method, rank, "abundance") %>%
        rename(Original_abundance = Abundance)

      if (original_units[[method]] == "percent") {
        raw <- raw %>%
          mutate(Original_abundance = Original_abundance / 100)
      } else if (original_units[[method]] != "fraction") {
        stop("Unknown original input unit for ", method)
      }

      dat <- full_join(
        norm,
        raw,
        by = c("SAMPLE_ID", "Taxon")
      ) %>%
        mutate(
          Relative_abundance_input =
            replace_na(Relative_abundance_input, 0),
          Original_abundance =
            replace_na(Original_abundance, 0)
        )

      if (!nrow(dat)) {
        stop("No retained taxa: ", key)
      }

      if (any(
        (dat$Relative_abundance_input > 0) !=
          (dat$Original_abundance > 0)
      )) {
        stop(
          "Original and normalized files disagree on ",
          "taxon presence/absence: ", key
        )
      }

      dat <- dat %>%
        group_by(SAMPLE_ID) %>%
        mutate(
          Retained_fraction = sum(Relative_abundance_input),
          Relative_abundance = if (Retained_fraction[1] > 0) {
            Relative_abundance_input / Retained_fraction[1]
          } else {
            rep(0, n())
          }
        ) %>%
        ungroup() %>%
        left_join(
          meta %>%
            select(SAMPLE_ID, Cohort_group, Country),
          by = "SAMPLE_ID"
        ) %>%
        mutate(
          Method = method,
          Rank = rank
        )

      qc <- dat %>%
        group_by(SAMPLE_ID) %>%
        summarise(
          Total_before = sum(Relative_abundance_input),
          Total_after = sum(Relative_abundance),
          Total_original = sum(Original_abundance),
          .groups = "drop"
        )

      qc <- meta %>%
        select(SAMPLE_ID, Cohort_group, Country) %>%
        left_join(qc, by = "SAMPLE_ID") %>%
        mutate(
          across(
            c(Total_before, Total_after, Total_original),
            ~ replace_na(.x, 0)
          ),
          Method = method,
          Rank = rank,
          Valid = abs(Total_after - 1) < 0.0001 &
            Total_original > 0
        )

      if (any(!qc$Valid)) {
        print(filter(qc, !Valid), n = Inf)
        stop("Invalid retained abundance totals: ", key)
      }

      if (any(qc$Total_original > 1.0001)) {
        stop(
          "Original totals exceed 1 after unit conversion: ",
          key, ". Check original_units and input files."
        )
      }

      profiles_list[[key]] <- dat
      qc_list[[key]] <- qc

      cat("Validated:", method, "|", rank, "\n")
    }
  }

  profiles <- bind_rows(profiles_list)
  sample_qc <- bind_rows(qc_list)

  write_csv(
    sample_qc,
    file.path(output_dir, "QC", "Sample_abundance_QC.csv")
  )

  complete_profile <- function(method, rank, groups) {

    ids <- meta %>%
      filter(as.character(Cohort_group) %in% groups) %>%
      pull(SAMPLE_ID)

    d <- profiles %>%
      filter(
        Method == method,
        Rank == rank,
        SAMPLE_ID %in% ids
      )

    taxa <- d %>%
      group_by(Taxon) %>%
      summarise(Total = sum(Original_abundance), .groups = "drop") %>%
      filter(Total > 0) %>%
      pull(Taxon) %>%
      sort()

    expand_grid(SAMPLE_ID = ids, Taxon = taxa) %>%
      left_join(
        d %>%
          select(
            SAMPLE_ID, Taxon,
            Relative_abundance,
            Original_abundance,
            Relative_abundance_input
          ),
        by = c("SAMPLE_ID", "Taxon")
      ) %>%
      mutate(
        across(
          c(
            Relative_abundance,
            Original_abundance,
            Relative_abundance_input
          ),
          ~ replace_na(.x, 0)
        )
      ) %>%
      left_join(
        meta %>%
          select(SAMPLE_ID, Cohort_group, Country),
        by = "SAMPLE_ID"
      )
  }

  abundance_summary <- function(d, grouping) {

    d %>%
      group_by(across(all_of(grouping))) %>%
      summarise(
        N = n(),
        N_detected = sum(Original_abundance > 0),
        Prevalence_percent = 100 * mean(Original_abundance > 0),

        Mean_original_percent = 100 * mean(Original_abundance),
        Median_original_percent = 100 * median(Original_abundance),
        Q1_original_percent = 100 * quantile(
          Original_abundance, 0.25,
          names = FALSE, type = 7
        ),
        Q3_original_percent = 100 * quantile(
          Original_abundance, 0.75,
          names = FALSE, type = 7
        ),

        Mean_percent = 100 * mean(Relative_abundance),
        Median_percent = 100 * median(Relative_abundance),
        Q1_percent = 100 * quantile(
          Relative_abundance, 0.25,
          names = FALSE, type = 7
        ),
        Q3_percent = 100 * quantile(
          Relative_abundance, 0.75,
          names = FALSE, type = 7
        ),
        .groups = "drop"
      )
  }

  taxon_palette <- c(
    "#F8766D", "#E69F00", "#A3A500", "#39B600",
    "#00BF7D", "#00BFC4", "#00B0F6", "#619CFF",
    "#B983FF", "#E76BF3", "#F564A3", "#7F7F7F",
    "#8C6D31", "#17BECF", "#BCBD22", "#9467BD",
    "#D62728", "#2CA02C", "#1F77B4", "#FF7F0E"
  )

  make_plots <- function(d, rank, method, cohort) {

    folder <- file.path(
      output_dir, "Figures", method, cohort
    )

    dir.create(folder, recursive = TRUE, showWarnings = FALSE)

    n_top <- top_n[[rank]]

    top_taxa <- d %>%
      group_by(Taxon) %>%
      summarise(
        Mean = mean(Relative_abundance),
        .groups = "drop"
      ) %>%
      arrange(desc(Mean), Taxon) %>%
      slice_head(n = n_top) %>%
      pull(Taxon)

    ids <- meta %>%
      filter(as.character(Cohort_group) %in% cohorts[[cohort]]) %>%
      arrange(Cohort_group, SAMPLE_ID) %>%
      pull(SAMPLE_ID)

    values <- d %>%
      mutate(
        Taxon_plot = if_else(
          Taxon %in% top_taxa, Taxon, "Other"
        )
      ) %>%
      group_by(SAMPLE_ID, Taxon_plot) %>%
      summarise(
        Value = sum(Relative_abundance),
        .groups = "drop"
      )

    plotting <- expand_grid(
      SAMPLE_ID = ids,
      Taxon_plot = c(top_taxa, "Other")
    ) %>%
      left_join(values, by = c("SAMPLE_ID", "Taxon_plot")) %>%
      mutate(Value = replace_na(Value, 0)) %>%
      left_join(
        meta %>% select(SAMPLE_ID, Cohort_group),
        by = "SAMPLE_ID"
      ) %>%
      mutate(
        SAMPLE_ID = factor(SAMPLE_ID, levels = ids),
        Cohort_group = factor(
          as.character(Cohort_group),
          levels = cohorts[[cohort]]
        ),
        Taxon_plot = factor(
          Taxon_plot,
          levels = c(top_taxa, "Other")
        )
      )

    plot_qc <- plotting %>%
      group_by(SAMPLE_ID) %>%
      summarise(Total = sum(Value), .groups = "drop")

    if (any(abs(plot_qc$Total - 1) > 0.0001)) {
      stop("Individual plotting totals are not 100%.")
    }

    group_means <- plotting %>%
      group_by(Cohort_group, Taxon_plot) %>%
      summarise(Mean = mean(Value), .groups = "drop")

    mean_qc <- group_means %>%
      group_by(Cohort_group) %>%
      summarise(Total = sum(Mean), .groups = "drop")

    if (any(abs(mean_qc$Total - 1) > 0.0001)) {
      stop("Group plotting totals are not 100%.")
    }

    colors <- setNames(
      c(
        rep(taxon_palette, length.out = length(top_taxa)),
        "#BDBDBD"
      ),
      c(top_taxa, "Other")
    )

    common_theme <- theme_bw(base_size = 11) +
      theme(panel.grid.major.x = element_blank())

    p_individual <- ggplot(
      plotting,
      aes(SAMPLE_ID, Value * 100, fill = Taxon_plot)
    ) +
      geom_col(width = 0.95) +
      facet_grid(
        ~ Cohort_group,
        scales = "free_x",
        space = "free_x",
        drop = TRUE
      ) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(
        breaks = seq(0, 100, 25),
        expand = expansion(mult = c(0, 0))
      ) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = paste(cohort, rank, method, sep = " | "),
        subtitle = paste0("Top ", n_top, " taxa + Other"),
        x = "Sample",
        y = "Renormalized relative abundance (%)",
        fill = rank
      ) +
      common_theme +
      theme(
        axis.text.x = element_text(
          angle = 90, hjust = 1, vjust = 0.5, size = 7
        )
      )

    p_mean <- ggplot(
      group_means,
      aes(Cohort_group, Mean * 100, fill = Taxon_plot)
    ) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(
        breaks = seq(0, 100, 25),
        expand = expansion(mult = c(0, 0))
      ) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = paste(cohort, rank, method, sep = " | "),
        subtitle = "Mean across individual samples",
        x = "MAP group",
        y = "Mean renormalized relative abundance (%)",
        fill = rank
      ) +
      common_theme

    base <- paste(cohort, rank, sep = "_")

    ggsave(
      file.path(folder, paste0(base, "_individual.png")),
      p_individual,
      width = max(11, length(ids) * 0.22 + 5),
      height = 7.5,
      dpi = 300
    )

    ggsave(
      file.path(folder, paste0(base, "_mean.png")),
      p_mean,
      width = 10,
      height = 7,
      dpi = 300
    )

    list(
      Means = group_means %>%
        mutate(
          Method = method,
          Cohort = cohort,
          Rank = rank,
          Mean_percent = Mean * 100
        ),
      QC = plot_qc %>%
        mutate(
          Method = method,
          Cohort = cohort,
          Rank = rank
        )
    )
  }

  all_overall <- list()
  all_group <- list()
  all_sample_richness <- list()
  all_richness_summary <- list()
  all_counts <- list()
  all_plot_means <- list()
  all_plot_qc <- list()

  for (cohort in names(cohorts)) {
    for (method in methods) {
      for (rank in ranks) {

        key <- paste(cohort, method, rank, sep = "__")

        d <- complete_profile(
          method, rank, cohorts[[cohort]]
        )

        if (!nrow(d)) next

        overall <- abundance_summary(d, "Taxon") %>%
          mutate(Cohort = cohort, Method = method, Rank = rank) %>%
          arrange(desc(Mean_percent), Taxon)

        by_group <- abundance_summary(
          d, c("Cohort_group", "Taxon")
        ) %>%
          mutate(Cohort = cohort, Method = method, Rank = rank)

        sample_richness <- d %>%
          group_by(SAMPLE_ID, Cohort_group, Country) %>%
          summarise(
            N_taxa = sum(Original_abundance > 0),
            .groups = "drop"
          ) %>%
          mutate(Cohort = cohort, Method = method, Rank = rank)

        richness_summary <- sample_richness %>%
          group_by(Cohort_group) %>%
          summarise(
            N_samples = n(),
            Median_taxa = median(N_taxa),
            Min_taxa = min(N_taxa),
            Max_taxa = max(N_taxa),
            .groups = "drop"
          ) %>%
          mutate(Cohort = cohort, Method = method, Rank = rank)

        counts <- tibble(
          Cohort = cohort,
          Method = method,
          Rank = rank,
          N_samples = n_distinct(d$SAMPLE_ID),
          N_distinct_taxa = n_distinct(d$Taxon)
        )

        all_overall[[key]] <- overall
        all_group[[key]] <- by_group
        all_sample_richness[[key]] <- sample_richness
        all_richness_summary[[key]] <- richness_summary
        all_counts[[key]] <- counts

        write_csv(
          d %>%
            mutate(
              Original_percent = Original_abundance * 100,
              Renormalized_percent = Relative_abundance * 100
            ),
          file.path(
            output_dir, "Data_used",
            paste0(key, ".csv")
          )
        )

        if (rank %in% c("Phylum", "Genus")) {
          plot_results <- make_plots(d, rank, method, cohort)
          all_plot_means[[key]] <- plot_results$Means
          all_plot_qc[[key]] <- plot_results$QC
        }

        cat(
          cohort, "|", method, "|", rank, "|",
          counts$N_distinct_taxa, "taxa detected\n"
        )
      }
    }
  }

  overall_abundances <- bind_rows(all_overall)
  group_abundances <- bind_rows(all_group)
  taxon_counts <- bind_rows(all_counts)

  safe_kw <- function(d) {

    if (n_distinct(d$Original_abundance) < 2) {
      return(tibble(
        H = 0,
        df = n_distinct(d$Cohort_group) - 1,
        p_value = 1,
        Test_status = "constant_values"
      ))
    }

    test <- kruskal.test(
      Original_abundance ~ droplevels(Cohort_group),
      data = d
    )

    tibble(
      H = unname(test$statistic),
      df = unname(test$parameter),
      p_value = test$p.value,
      Test_status = "tested"
    )
  }

  safe_wilcox <- function(d, pair) {

    x <- d$Original_abundance[
      as.character(d$Cohort_group) == pair[1]
    ]

    y <- d$Original_abundance[
      as.character(d$Cohort_group) == pair[2]
    ]

    if (!length(x) || !length(y)) {
      stop("A requested comparison group is missing.")
    }

    result <- tibble(
      Group1 = pair[1],
      Group2 = pair[2],
      N1 = length(x),
      N2 = length(y),
      Mean1_original_percent = 100 * mean(x),
      Mean2_original_percent = 100 * mean(y),
      Mean_difference_percentage_points =
        100 * (mean(x) - mean(y)),
      W = NA_real_,
      p_value = NA_real_,
      Test_status = "not_testable_constant_values"
    )

    if (n_distinct(c(x, y)) < 2) {
      return(result)
    }

    test <- wilcox.test(
      x, y,
      alternative = "two.sided",
      paired = FALSE,
      exact = FALSE,
      correct = TRUE
    )

    result$W <- unname(test$statistic)
    result$p_value <- test$p.value
    result$Test_status <- "tested"

    result
  }

  adjust_bh <- function(p) {
    result <- rep(NA_real_, length(p))
    valid <- is.finite(p)
    result[valid] <- p.adjust(p[valid], method = "BH")
    result
  }

  empty_pairwise <- tibble(
    Taxon = character(),
    Group1 = character(),
    Group2 = character(),
    N1 = integer(),
    N2 = integer(),
    Mean1_original_percent = double(),
    Mean2_original_percent = double(),
    Mean_difference_percentage_points = double(),
    W = double(),
    p_value = double(),
    Test_status = character(),
    q_BH = double(),
    Significant_BH = logical(),
    Method = character(),
    Rank = character(),
    Comparison = character()
  )

  sw_kw_list <- list()
  sw_pairwise_list <- list()

  sw_pairs <- list(
    c("SW_NEG", "SW_SUB"),
    c("SW_NEG", "SW_POS"),
    c("SW_SUB", "SW_POS")
  )

  for (method in methods) {
    for (rank in ranks) {

      key <- paste(method, rank, sep = "__")

      d <- complete_profile(
        method, rank, cohorts$Switzerland
      )

      if (!nrow(d)) next

      kw <- d %>%
        group_by(Taxon) %>%
        group_modify(~ safe_kw(.x)) %>%
        ungroup() %>%
        mutate(
          q_BH = adjust_bh(p_value),
          Significant_BH = q_BH < alpha,
          Method = method,
          Rank = rank,
          Comparison = "Switzerland"
        ) %>%
        arrange(q_BH, p_value, Taxon)

      significant_taxa <- kw %>%
        filter(Significant_BH) %>%
        pull(Taxon)

      pairwise <- empty_pairwise

      if (length(significant_taxa)) {

        pairwise <- map_dfr(
          significant_taxa,
          function(taxon) {

            taxon_data <- d %>% filter(Taxon == taxon)

            map_dfr(
              sw_pairs,
              ~ safe_wilcox(taxon_data, .x)
            ) %>%
              mutate(Taxon = taxon, .before = 1)
          }
        ) %>%
          group_by(Taxon) %>%
          mutate(
            q_BH = adjust_bh(p_value),
            Significant_BH = q_BH < alpha
          ) %>%
          ungroup() %>%
          mutate(
            Method = method,
            Rank = rank,
            Comparison = "Switzerland_posthoc"
          )
      }

      sw_kw_list[[key]] <- kw
      sw_pairwise_list[[key]] <- pairwise
    }
  }

  sw_kw <- bind_rows(sw_kw_list)
  sw_pairwise <- bind_rows(
    empty_pairwise,
    bind_rows(sw_pairwise_list)
  )

  secondary_list <- list()

  for (comparison in names(secondary_comparisons)) {
    for (method in methods) {
      for (rank in ranks) {

        pair <- secondary_comparisons[[comparison]]

        d <- complete_profile(method, rank, pair)

        if (!nrow(d)) next

        result <- d %>%
          group_by(Taxon) %>%
          group_modify(~ safe_wilcox(.x, pair)) %>%
          ungroup() %>%
          mutate(
            q_BH = adjust_bh(p_value),
            Significant_BH = q_BH < alpha,
            Method = method,
            Rank = rank,
            Comparison = comparison
          ) %>%
          arrange(q_BH, p_value, Taxon)

        key <- paste(comparison, method, rank, sep = "__")
        secondary_list[[key]] <- result
      }
    }
  }

  secondary_tests <- bind_rows(
    empty_pairwise,
    bind_rows(secondary_list)
  )

  sw_overall <- overall_abundances %>%
    filter(Cohort == "Switzerland") %>%
    select(-Cohort)

  sw_group_wide <- group_abundances %>%
    filter(Cohort == "Switzerland") %>%
    select(-Cohort) %>%
    mutate(Cohort_group = as.character(Cohort_group)) %>%
    pivot_wider(
      id_cols = c(Method, Rank, Taxon),
      names_from = Cohort_group,
      values_from = c(
        N, N_detected, Prevalence_percent,
        Mean_original_percent,
        Median_original_percent,
        Q1_original_percent,
        Q3_original_percent,
        Mean_percent, Median_percent,
        Q1_percent, Q3_percent
      ),
      names_glue = "{Cohort_group}_{.value}"
    )

  supplementary_swiss <- sw_overall %>%
    left_join(
      sw_group_wide,
      by = c("Method", "Rank", "Taxon")
    ) %>%
    left_join(
      sw_kw %>%
        select(
          Method, Rank, Taxon,
          H, df, p_value, q_BH,
          Significant_BH, Test_status
        ),
      by = c("Method", "Rank", "Taxon")
    ) %>%
    arrange(
      match(Method, methods),
      match(Rank, ranks),
      desc(Mean_original_percent),
      Taxon
    )

  kw_overview <- sw_kw %>%
    group_by(Method, Rank, Comparison) %>%
    summarise(
      Test = "Kruskal-Wallis",
      N_tests = sum(is.finite(p_value)),
      N_significant_BH = sum(Significant_BH, na.rm = TRUE),
      .groups = "drop"
    )

  pairwise_overview <- sw_pairwise %>%
    group_by(Method, Rank, Comparison) %>%
    summarise(
      Test = "Wilcoxon; BH within taxon",
      N_tests = sum(is.finite(p_value)),
      N_significant_BH = sum(Significant_BH, na.rm = TRUE),
      .groups = "drop"
    )

  secondary_overview <- secondary_tests %>%
    group_by(Method, Rank, Comparison) %>%
    summarise(
      Test = "Wilcoxon",
      N_tests = sum(is.finite(p_value)),
      N_significant_BH = sum(Significant_BH, na.rm = TRUE),
      .groups = "drop"
    )

  statistical_overview <- bind_rows(
    kw_overview,
    pairwise_overview,
    secondary_overview
  )

  table_outputs <- list(
    Taxonomic_summary = taxon_counts,
    Overall_abundances = overall_abundances,
    Group_abundances = group_abundances,
    Supplementary_Swiss = supplementary_swiss,
    Sample_taxon_richness = bind_rows(all_sample_richness),
    Richness_summary = bind_rows(all_richness_summary)
  )

  for (name in names(table_outputs)) {
    write_csv(
      table_outputs[[name]],
      file.path(
        output_dir, "Tables",
        paste0(name, ".csv")
      )
    )
  }

  stat_outputs <- list(
    Swiss_Kruskal_Wallis = sw_kw,
    Swiss_pairwise_Wilcoxon = sw_pairwise,
    Secondary_Wilcoxon = secondary_tests,
    Statistical_overview = statistical_overview
  )

  for (name in names(stat_outputs)) {
    write_csv(
      stat_outputs[[name]],
      file.path(
        output_dir, "Statistics",
        paste0(name, ".csv")
      )
    )
  }

  analysis_notes <- tribble(
    ~Item, ~Details,

    "Primary method",
    "MetaPhlAn_GTDB",

    "Supportive method",
    "ganon2",

    "Sample inclusion",
    paste(
      "Metadata-eligible samples shared across all 12 input files;",
      "MAP detection, EX, ZA_NEG_01 and ZA_NEG_08 excluded."
    ),

    "Taxon inclusion",
    paste(
      "Named retained taxa; exact placeholder labels and strain",
      "rows excluded. Only taxa detected in the selected cohort",
      "or comparison are summarized and tested."
    ),

    "Phylum harmonization",
    paste(
      "Bacillota subdivisions and Firmicutes merged as Bacillota;",
      "Proteobacteria, Bacteroidetes and Actinobacteria renamed."
    ),

    "Original units",
    paste(
      "MetaPhlAn_GTDB input percentages divided by 100;",
      "ganon2 original input assumed to contain fractions."
    ),

    "Abundance reporting",
    paste(
      "Original and renormalized summaries reported as percent.",
      "Zeros included in all means, medians and quartiles.",
      "Quartiles use quantile(type = 7)."
    ),

    "Original abundance interpretation",
    paste(
      "Original profiler estimates retain their source denominator;",
      "MetaPhlAn estimates are not directly measured read/base fractions."
    ),

    "Figures",
    paste(
      "Normalized input rescaled to 100% across retained taxa",
      "within each sample and rank. Top taxa selected by cohort",
      "mean across individuals. Group means are not renormalized again."
    ),

    "Swiss omnibus tests",
    paste(
      "Kruskal-Wallis on original abundances;",
      "BH separately within each method and rank."
    ),

    "Swiss post hoc tests",
    paste(
      "Only omnibus q < 0.05 taxa;",
      "two-sided unpaired Wilcoxon, exact = FALSE, correct = TRUE;",
      "BH across testable pairs within each taxon."
    ),

    "Post hoc multiplicity scope",
    paste(
      "Within-taxon pairwise q values do not establish global",
      "FDR control across all selected pairwise comparisons."
    ),

    "Secondary comparisons",
    paste(
      "ZA_NEG vs ZA_SUB; SW_NEG vs ZA_NEG; SW_SUB vs ZA_SUB.",
      "Wilcoxon with BH separately within method/rank/comparison."
    ),

    "Constant values",
    paste(
      "Constant omnibus data assigned p = 1.",
      "Constant two-group comparisons marked not testable",
      "and excluded from pairwise BH adjustment."
    ),

    "Design limitations",
    paste(
      "Exploratory relative-abundance associations.",
      "Farm clustering, confounding and compositionality are",
      "not modeled. Country comparisons also reflect differences",
      "in environment and sample processing."
    ),

    "Mycobacterium interpretation",
    "Species-level M. avium assignments do not establish MAP identity.",

    "Scope",
    paste(
      "Composition and differential abundance only.",
      "Alpha/beta diversity and indicator-species testing are separate."
    )
  )

  input_manifest <- files %>%
    mutate(
      MD5 = unname(tools::md5sum(Path))
    )

  input_manifest <- bind_rows(
    tibble(
      Method = "Metadata",
      Rank = NA_character_,
      Kind = "metadata",
      Path = meta_file,
      MD5 = unname(tools::md5sum(meta_file))
    ),
    input_manifest
  )

  write_csv(
    input_manifest,
    file.path(
      output_dir, "Reproducibility",
      "Input_file_manifest.csv"
    )
  )

  write_csv(
    analysis_notes,
    file.path(
      output_dir, "Reproducibility",
      "Analysis_methods.csv"
    )
  )

  writeLines(
    capture.output(sessionInfo()),
    file.path(
      output_dir, "Reproducibility",
      "R_sessionInfo.txt"
    )
  )

  master <- list(
    Sample_inclusion = meta_all %>%
      select(SAMPLE_ID, GROUP, Cohort_group, Exclusion_reason),

    Supplementary_Swiss = supplementary_swiss,

    Statistical_overview = statistical_overview,

    Swiss_KW = sw_kw,

    Swiss_pairwise = sw_pairwise,

    Secondary_Wilcoxon = secondary_tests,

    Taxon_counts = taxon_counts,

    Overall_abundances = overall_abundances,

    Group_abundances = group_abundances,

    Sample_QC = sample_qc,

    Sample_taxon_richness = bind_rows(all_sample_richness),

    Richness_summary = bind_rows(all_richness_summary),

    Figure_group_means = bind_rows(all_plot_means),

    Figure_sample_QC = bind_rows(all_plot_qc),

    Methods = analysis_notes,

    Input_files = input_manifest
  )

  write_xlsx(
    master,
    file.path(
      output_dir,
      "DVM_Taxonomic_Composition_and_DA_Master.xlsx"
    )
  )

  cat(
    "\n============================================\n",
    "COMPOSITION + DIFFERENTIAL ABUNDANCE COMPLETE\n",
    "============================================\n\n",
    sep = ""
  )

  cat("Included samples:\n")
  print(sample_counts, n = Inf)

  cat("\nDetected taxa:\n")
  print(taxon_counts, n = Inf)

  cat("\nStatistical overview:\n")
  print(statistical_overview, n = Inf)

  cat(
    "\nAll outputs saved to:\n",
    normalizePath(output_dir),
    "\n"
  )
})

# Five-classifier taxonomic comparison
local({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(ggplot2)

  input_dir <- file.path("data", "taxonomy")
  output_dir <- "DVM_taxonomic_comparison_5methods_v4"

  dir.create(output_dir, recursive = TRUE,
             showWarnings = FALSE)

  methods <- c(
    "MetaPhlAn_GTDB",
    "MetaPhlAn_default",
    "Kraken2",
    "Kraken2_conf10",
    "ganon2"
  )

  ranking_methods <- c(
    "MetaPhlAn_GTDB",
    "MetaPhlAn_default",
    "Kraken2",
    "ganon2"
  )

  method_labels <- c(
    MetaPhlAn_GTDB = "MetaPhlAn GTDB",
    MetaPhlAn_default = "MetaPhlAn default",
    Kraken2 = "Kraken2 0.00",
    Kraken2_conf10 = "Kraken2 0.10",
    ganon2 = "ganon2"
  )

  ranks <- c("Phylum", "Genus", "Species")
  group_order <- c("MAP-negative", "Subclinical", "Clinical")

  # Sample presence is derived from the deposited profiles, not a separate QC workbook.
  meta <- read_tsv(
    file.path("data", "metadata.txt"),
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  ) %>% mutate(SAMPLE_ID = as.character(SAMPLE_ID),
               GROUP = as.character(GROUP))

  get_header_ids <- function(method, rank) {
    f <- taxonomy_file(method, rank, "normalized")
    if (!file.exists(f)) stop("Missing classifier input: ", f)
    setdiff(names(read_tsv(f, n_max = 0, show_col_types = FALSE,
                            name_repair = "check_unique")), "clade_name")
  }
  method_ids <- lapply(methods, function(m) {
    Reduce(intersect, lapply(ranks, function(r) get_header_ids(m, r)))
  })
  names(method_ids) <- methods
  common_samples <- Reduce(intersect, method_ids)

  # SW_NEG_09 has no Kraken2-conf0.10 genus/species assignments in the archive.
  # Exclude empty profiles from the matched five-classifier comparison.
  empty_conf10 <- unique(unlist(lapply(c("Genus", "Species"), function(r) {
    f <- taxonomy_file("Kraken2_conf10", r, "normalized")
    x <- read_tsv(f, show_col_types = FALSE,
                  col_types = cols(.default = col_character()))
    ids <- intersect(setdiff(names(x), "clade_name"), common_samples)
    ids[vapply(x[ids], function(v) {
      nums <- suppressWarnings(as.numeric(v))
      all(is.na(nums)) || sum(nums, na.rm = TRUE) <= 0
    }, logical(1))]
  })))
  common_samples <- setdiff(common_samples, empty_conf10)
  presence_qc <- tibble(SAMPLE_ID = meta$SAMPLE_ID)
  for (m in methods) presence_qc[[m]] <- presence_qc$SAMPLE_ID %in% method_ids[[m]]
  presence_qc$No_conf10_assignments <- presence_qc$SAMPLE_ID %in% empty_conf10
  presence_qc$Included_all_five <- presence_qc$SAMPLE_ID %in% common_samples
  write_csv(presence_qc, file.path(output_dir, "Matched_sample_presence_QC.csv"))

  swiss <- meta %>%
    filter(
      SAMPLE_ID %in% common_samples,
      str_starts(SAMPLE_ID, "SW_"),
      !GROUP %in% c("EX", "Excluded", "MAP detection"),
      SAMPLE_ID != "SW_NEG_01"
    ) %>%
    mutate(
      MAP_group = case_when(
        str_starts(SAMPLE_ID, "SW_NEG_") ~ "MAP-negative",
        str_starts(SAMPLE_ID, "SW_SUB_") ~ "Subclinical",
        str_starts(SAMPLE_ID, "SW_POS_") ~ "Clinical",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(MAP_group)) %>%
    mutate(
      MAP_group = factor(MAP_group, levels = group_order)
    )

  if (nrow(swiss) != 35 || any(table(swiss$MAP_group) != c(15, 13, 7))) {
    stop("Expected matched Swiss cohort: 35 samples (15/13/7). ",
         "Check classifier inputs and metadata.")
  }
  stopifnot(nrow(swiss) > 0)
  stopifnot(!anyDuplicated(swiss$SAMPLE_ID))

  message("Swiss samples: ", nrow(swiss))
  print(table(swiss$MAP_group))

  harmonize_phylum <- function(x) {

    case_when(
      str_detect(x, "^Bacillota(_[A-Z0-9]+)?$") |
        x == "Firmicutes" ~ "Bacillota",

      x %in% c("Bacteroidetes", "Bacteroidota") ~
        "Bacteroidota",

      x %in% c("Proteobacteria", "Pseudomonadota") ~
        "Pseudomonadota",

      x %in% c("Actinobacteria", "Actinomycetota") ~
        "Actinomycetota",

      TRUE ~ x
    )
  }

  read_profile <- function(method, rank, type) {

    file <- taxonomy_file(method, rank, type)

    if (!file.exists(file)) return(NULL)

    dat <- read_tsv(file, show_col_types = FALSE)

    if (!"clade_name" %in% names(dat)) {
      stop("Missing clade_name: ", file)
    }

    samples <- swiss$SAMPLE_ID

    missing_samples <- setdiff(samples, names(dat))

    if (length(missing_samples) > 0) {
      stop(
        "Missing samples in ", file, ": ",
        paste(missing_samples, collapse = ", ")
      )
    }

    if (
      rank == "Species" &&
      any(
        str_detect(dat$clade_name, "(^|[;|])t__"),
        na.rm = TRUE
      )
    ) {
      stop("Strain entries found in species table: ", file)
    }

    dat %>%
      mutate(
        Original_taxon = sub("^.*[;|]", "", clade_name),
        Original_taxon = sub(
          "^[dkpcofgst]__", "", Original_taxon
        ),
        Taxon = if (rank == "Phylum") {
          harmonize_phylum(Original_taxon)
        } else {
          Original_taxon
        }
      ) %>%
      filter(
        !is.na(Taxon),
        nzchar(trimws(Taxon)),
        !str_detect(
          str_to_lower(Taxon),
          "^(unclassified|unknown|unassigned|na)$"
        )
      ) %>%
      select(Taxon, all_of(samples)) %>%
      pivot_longer(
        cols = all_of(samples),
        names_to = "SAMPLE_ID",
        values_to = "Value"
      ) %>%
      mutate(Value = as.numeric(Value)) %>%
      group_by(Taxon, SAMPLE_ID) %>%
      summarise(
        Value = sum(Value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        Method = method,
        Rank = rank,
        Measure = type
      )
  }

  normalized <- map_dfr(methods, function(m) {
    map_dfr(ranks, function(r) {
      read_profile(m, r, "normalized")
    })
  })

  original <- map_dfr(methods, function(m) {
    map_dfr(ranks, function(r) {
      read_profile(m, r, "abundance")
    })
  })

  ganon_counts <- map_dfr(ranks, function(r) {
    read_profile("ganon2", r, "counts")
  })

  complete_profiles <- function(dat) {

    dat %>%
      group_by(Method, Rank, Measure) %>%
      group_modify(~ {

        expand_grid(
          Taxon = unique(.x$Taxon),
          SAMPLE_ID = swiss$SAMPLE_ID
        ) %>%
          left_join(
            .x %>%
              select(Taxon, SAMPLE_ID, Value),
            by = c("Taxon", "SAMPLE_ID")
          ) %>%
          mutate(Value = replace_na(Value, 0))
      }) %>%
      ungroup()
  }

  normalized <- complete_profiles(normalized)

  original_qc <- original %>%
    group_by(Method, Rank, SAMPLE_ID) %>%
    summarise(
      Original_sum = sum(Value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      swiss %>% select(SAMPLE_ID, MAP_group),
      by = "SAMPLE_ID"
    ) %>%
    mutate(
      Original_sum_percent = if_else(
        Method == "ganon2",
        Original_sum * 100,
        Original_sum
      ),
      Metric = case_when(
        str_starts(Method, "MetaPhlAn") ~
          "Original MetaPhlAn abundance sum; NOT read classification %",
        str_starts(Method, "Kraken2") ~
          "Original Kraken2 abundance sum; verify denominator",
        Method == "ganon2" ~
          "Original ganon2 abundance sum; verify denominator"
      )
    )

  original_summary <- original_qc %>%
    group_by(Method, Rank) %>%
    summarise(
      N_samples = n(),
      Mean_percent = mean(Original_sum_percent),
      Median_percent = median(Original_sum_percent),
      Min_percent = min(Original_sum_percent),
      Max_percent = max(Original_sum_percent),
      .groups = "drop"
    )

  normalization_qc <- normalized %>%
    group_by(Method, Rank, SAMPLE_ID) %>%
    summarise(
      Normalized_sum = sum(Value, na.rm = TRUE),
      .groups = "drop"
    )

  zero_sum_samples <- normalization_qc %>%
    filter(Normalized_sum < 0.999)

  if (nrow(zero_sum_samples) > 0) {
    warning(
      "Some normalized sample/rank profiles do not sum to 1. ",
      "Check Zero_sum_QC before interpreting results."
    )
  }

  if (nrow(ganon_counts) > 0) {

    ganon_count_qc <- ganon_counts %>%
      group_by(Rank, SAMPLE_ID) %>%
      summarise(
        Assigned_rank_counts = sum(Value, na.rm = TRUE),
        .groups = "drop"
      )

  } else {

    ganon_count_qc <- tibble(
      Rank = character(),
      SAMPLE_ID = character(),
      Assigned_rank_counts = numeric()
    )
  }

  summary_all <- normalized %>%
    group_by(Method, Rank, Taxon) %>%
    summarise(
      Mean_abundance = mean(Value),
      Median_abundance = median(Value),
      Prevalence = mean(Value > 0),
      .groups = "drop"
    ) %>%
    mutate(
      Mean_percent = Mean_abundance * 100,
      Median_percent = Median_abundance * 100,
      Prevalence_percent = Prevalence * 100
    )

  group_summary <- normalized %>%
    left_join(
      swiss %>% select(SAMPLE_ID, MAP_group),
      by = "SAMPLE_ID"
    ) %>%
    group_by(Method, Rank, MAP_group, Taxon) %>%
    summarise(
      Mean_percent = mean(Value) * 100,
      Prevalence_percent = mean(Value > 0) * 100,
      .groups = "drop"
    )

  independent_top5 <- summary_all %>%
    group_by(Method, Rank) %>%
    arrange(
      desc(Mean_abundance),
      Taxon,
      .by_group = TRUE
    ) %>%
    slice_head(n = 5) %>%
    mutate(Position = row_number()) %>%
    ungroup() %>%
    select(
      Method, Rank, Position, Taxon,
      Mean_percent, Median_percent,
      Prevalence_percent
    )

  top5_compact <- independent_top5 %>%
    mutate(
      Entry = paste0(
        Taxon,
        " (",
        round(Mean_percent, 2),
        "%)"
      )
    ) %>%
    select(Method, Rank, Position, Entry) %>%
    pivot_wider(
      names_from = Position,
      values_from = Entry,
      names_prefix = "Top_"
    ) %>%
    arrange(
      factor(Method, levels = methods),
      factor(Rank, levels = ranks)
    )

  shared_ranking <- summary_all %>%
    filter(
      Method %in% ranking_methods,
      Rank %in% c("Genus", "Species")
    ) %>%
    select(Method, Rank, Taxon, Mean_abundance) %>%
    group_by(Rank) %>%
    complete(
      Taxon,
      Method = ranking_methods,
      fill = list(Mean_abundance = 0)
    ) %>%
    ungroup() %>%
    group_by(Rank, Taxon) %>%
    summarise(
      Combined_mean = mean(Mean_abundance),
      Methods_detected = sum(Mean_abundance > 0),
      Max_method_mean = max(Mean_abundance),
      .groups = "drop"
    ) %>%
    group_by(Rank) %>%
    arrange(
      desc(Combined_mean),
      Taxon,
      .by_group = TRUE
    ) %>%
    mutate(
      Shared_rank = row_number()
    ) %>%
    ungroup()

  shared_top20 <- shared_ranking %>%
    filter(Shared_rank <= 20)

  shared_comparison <- shared_top20 %>%
    select(Rank, Taxon, Shared_rank) %>%
    crossing(Method = methods) %>%
    left_join(
      summary_all,
      by = c("Rank", "Taxon", "Method")
    ) %>%
    mutate(
      Mean_percent = replace_na(Mean_percent, 0),
      Prevalence_percent = replace_na(
        Prevalence_percent, 0
      )
    ) %>%
    arrange(Rank, Shared_rank, Method)

  make_heatmap <- function(rank) {

    selected <- shared_top20 %>%
      filter(Rank == rank) %>%
      arrange(Shared_rank) %>%
      pull(Taxon)

    dat <- shared_comparison %>%
      filter(Rank == rank) %>%
      mutate(
        Method = factor(
          Method,
          levels = methods,
          labels = method_labels[methods]
        ),
        Taxon = factor(
          Taxon,
          levels = rev(selected)
        )
      )

    p <- ggplot(
      dat,
      aes(
        x = Method,
        y = Taxon,
        fill = Mean_percent
      )
    ) +
      geom_tile(
        color = "white",
        linewidth = 0.3
      ) +
      scale_fill_gradient(
        low = "white",
        high = "#21618C",
        name = "Mean (%)"
      ) +
      labs(
        title = paste0(
          "Shared top 20 ",
          tolower(rank)
        ),
        subtitle = paste(
          "Selected using four classifiers;",
          "Kraken2 0.10 shown for comparison"
        ),
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(
          angle = 35, hjust = 1
        ),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold")
      )

    ggsave(
      file.path(
        output_dir,
        paste0(rank, "_revised_top20_heatmap.png")
      ),
      p,
      width = 11,
      height = 8,
      dpi = 300
    )

    p
  }

  genus_heatmap <- make_heatmap("Genus")
  species_heatmap <- make_heatmap("Species")

  make_stacked_plot <- function(rank) {

    selected <- shared_top20 %>%
      filter(Rank == rank) %>%
      arrange(Shared_rank) %>%
      pull(Taxon)

    dat <- normalized %>%
      filter(Rank == rank) %>%
      mutate(
        Plot_taxon = if_else(
          Taxon %in% selected,
          Taxon,
          "Other"
        )
      ) %>%
      group_by(Method, SAMPLE_ID, Plot_taxon) %>%
      summarise(
        Value = sum(Value),
        .groups = "drop"
      )

    dat <- expand_grid(
      Method = methods,
      SAMPLE_ID = swiss$SAMPLE_ID,
      Plot_taxon = c(selected, "Other")
    ) %>%
      left_join(
        dat,
        by = c("Method", "SAMPLE_ID", "Plot_taxon")
      ) %>%
      mutate(Value = replace_na(Value, 0)) %>%
      left_join(
        swiss %>% select(SAMPLE_ID, MAP_group),
        by = "SAMPLE_ID"
      )

    plot_data <- dat %>%
      group_by(Method, MAP_group, Plot_taxon) %>%
      summarise(
        Mean_percent = mean(Value) * 100,
        .groups = "drop"
      )

    colors <- grDevices::hcl.colors(
      length(selected),
      palette = "Dynamic"
    )

    names(colors) <- selected
    colors <- c(colors, Other = "#BDBDBD")

    plot_data <- plot_data %>%
      mutate(
        Method = factor(
          Method,
          levels = methods,
          labels = method_labels[methods]
        ),
        MAP_group = factor(
          MAP_group,
          levels = group_order
        ),
        Plot_taxon = factor(
          Plot_taxon,
          levels = c(selected, "Other")
        )
      )

    p <- ggplot(
      plot_data,
      aes(
        x = MAP_group,
        y = Mean_percent,
        fill = Plot_taxon
      )
    ) +
      geom_col(width = 0.75) +
      facet_wrap(~ Method, ncol = 1) +
      scale_fill_manual(
        values = colors,
        drop = FALSE
      ) +
      coord_cartesian(ylim = c(0, 100)) +
      labs(
        title = paste0(
          "Shared top 20 ",
          tolower(rank),
          " — revised selection"
        ),
        subtitle = paste(
          "Swiss cohort;",
          "within-rank normalized abundances"
        ),
        x = NULL,
        y = "Mean relative abundance (%)",
        fill = rank
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(
          angle = 20, hjust = 1
        ),
        plot.title = element_text(face = "bold")
      )

    ggsave(
      file.path(
        output_dir,
        paste0(rank, "_revised_top20_stacked.png")
      ),
      p,
      width = 14,
      height = 15,
      dpi = 300
    )

    write_tsv(
      plot_data %>%
        mutate(across(where(is.factor), as.character)),
      file.path(
        output_dir,
        paste0(rank, "_revised_top20_plot_data.tsv")
      )
    )

    p
  }

  genus_stacked <- make_stacked_plot("Genus")
  species_stacked <- make_stacked_plot("Species")

  make_top5_plot <- function(rank) {

    dat <- independent_top5 %>%
      filter(Rank == rank) %>%
      mutate(
        Method = factor(
          Method,
          levels = methods,
          labels = method_labels[methods]
        ),
        Taxon_label = paste0(
          Taxon,
          " (", round(Mean_percent, 1), "%)"
        ),
        Taxon_label = factor(
          Taxon_label,
          levels = rev(unique(Taxon_label))
        )
      )

    p <- ggplot(
      dat,
      aes(
        x = Mean_percent,
        y = reorder_within_helper,
        fill = Method
      )
    )

    dat <- dat %>%
      mutate(
        Panel_taxon = paste(
          Method,
          Taxon,
          sep = "___"
        )
      ) %>%
      group_by(Method) %>%
      arrange(Mean_percent, .by_group = TRUE) %>%
      mutate(
        Panel_taxon = factor(
          Panel_taxon,
          levels = unique(Panel_taxon)
        )
      ) %>%
      ungroup()

    p <- ggplot(
      dat,
      aes(
        x = Mean_percent,
        y = reorder(
          Panel_taxon,
          Mean_percent
        )
      )
    ) +
      geom_col(fill = "#35769B") +
      facet_wrap(
        ~ Method,
        scales = "free_y",
        ncol = 1
      ) +
      scale_y_discrete(
        labels = function(x) {
          sub("^.*___", "", x)
        }
      ) +
      labs(
        title = paste(
          "Independent top 5",
          tolower(rank),
          "per classifier"
        ),
        subtitle = "Mean normalized relative abundance",
        x = "Mean relative abundance (%)",
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold")
      )

    ggsave(
      file.path(
        output_dir,
        paste0(rank, "_independent_top5.png")
      ),
      p,
      width = 11,
      height = 15,
      dpi = 300
    )

    p
  }

  make_top5_plot <- function(rank) {

    dat <- independent_top5 %>%
      filter(Rank == rank) %>%
      mutate(
        Method = factor(
          Method,
          levels = methods,
          labels = method_labels[methods]
        ),
        Taxon_label = paste0(
          Taxon,
          " [", as.character(Method), "]"
        )
      )

    p <- ggplot(
      dat,
      aes(
        x = Mean_percent,
        y = reorder(Taxon_label, Mean_percent)
      )
    ) +
      geom_col(fill = "#35769B") +
      facet_wrap(
        ~ Method,
        scales = "free_y",
        ncol = 1
      ) +
      scale_y_discrete(
        labels = function(x) {
          sub(" \\[.*$", "", x)
        }
      ) +
      labs(
        title = paste0(
          "Independent top 5 ",
          tolower(rank),
          " per classifier"
        ),
        subtitle = "Swiss samples; mean normalized abundance",
        x = "Mean relative abundance (%)",
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold")
      )

    ggsave(
      file.path(
        output_dir,
        paste0(rank, "_independent_top5.png")
      ),
      p,
      width = 11,
      height = 15,
      dpi = 300
    )

    p
  }

  phylum_top5_plot <- make_top5_plot("Phylum")
  genus_top5_plot <- make_top5_plot("Genus")
  species_top5_plot <- make_top5_plot("Species")

  p_qc <- ggplot(
    original_qc,
    aes(
      x = Method,
      y = Original_sum_percent
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      fill = "#D6EAF8"
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.45,
      size = 1
    ) +
    facet_wrap(
      ~ Rank,
      scales = "free_y"
    ) +
    labs(
      title = "Original abundance sums by classifier",
      subtitle = paste(
        "QC proxy only — NOT directly comparable",
        "read-classification fractions"
      ),
      x = NULL,
      y = "Original abundance sum (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      )
    )

  ggsave(
    file.path(
      output_dir,
      "Original_abundance_QC.png"
    ),
    p_qc,
    width = 13,
    height = 7,
    dpi = 300
  )

  write_xlsx(
    list(
      Samples = swiss,
      Original_QC_summary = original_summary,
      Original_QC_samples = original_qc,
      Normalization_QC = normalization_qc,
      Zero_sum_QC = zero_sum_samples,
      ganon_count_QC = ganon_count_qc,
      Top5_compact = top5_compact,
      Top5_detailed = independent_top5,
      Shared_top20 = shared_top20,
      Shared_comparison = shared_comparison,
      MAP_group_summary = group_summary,
      All_taxa = summary_all
    ),
    file.path(
      output_dir,
      "Taxonomic_comparison_5methods_v4.xlsx"
    )
  )

  message("\n===================================")
  message("ANALYSIS COMPLETE")
  message("===================================")

  message("\nORIGINAL ABUNDANCE QC:")
  print(original_summary, n = Inf)

  message("\nINDEPENDENT TOP 5:")
  print(top5_compact, n = Inf, width = Inf)

  message("\nREVISED SHARED TOP 20:")
  print(shared_top20, n = Inf)

  message("\nOUTPUT:")
  message(output_dir)
})
