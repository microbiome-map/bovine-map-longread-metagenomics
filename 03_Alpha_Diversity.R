# Bovine MAP metagenomics: alpha diversity

# Map analysis methods to the original filenames in the Zenodo archive.
taxonomy_file <- function(method, rank, kind) {
  prefix <- ifelse(method == "MetaPhlAn_default", "MetaPhlAn", method)
  file.path("data", "taxonomy", method,
            paste0(prefix, "_", rank, "_", kind, ".tsv"))
}
# Primary: MetaPhlAn GTDB; supportive: five-method matched Swiss comparison.

packages <- c("tidyverse", "openxlsx")

missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}

library(tidyverse)
library(openxlsx)

input_dir <- file.path("data", "taxonomy")

metadata_file <- file.path("data", "metadata.txt")
if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

output_dir <- "DVM_Alpha_Diversity_FINAL"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

methods <- c(
  "MetaPhlAn_GTDB",
  "MetaPhlAn_default",
  "ganon2",
  "Kraken2",
  "Kraken2_conf10"
)

primary_method <- "MetaPhlAn_GTDB"

ranks <- c("Genus", "Species")

metrics <- c(
  "Richness",
  "Shannon",
  "Gini_Simpson"
)

group_order <- c(
  "SW_NEG",
  "SW_SUB",
  "SW_POS",
  "ZA_NEG",
  "ZA_SUB"
)

excluded_ids <- c(
  "SW_NEG_01",
  "ZA_NEG_01",
  "ZA_NEG_08"
)

cohorts <- list(
  Switzerland = c("SW_NEG", "SW_SUB", "SW_POS"),
  South_Africa = c("ZA_NEG", "ZA_SUB"),
  Combined = group_order
)

min_farm_n <- 2

# Exact sample counts are reported, not hard-coded.

group_colors <- c(
  SW_NEG = "#2793A8",
  SW_SUB = "#9C83C6",
  SW_POS = "#EE8D79",
  ZA_NEG = "#72C5D7",
  ZA_SUB = "#C7B6E0"
)

method_colors <- c(
  MetaPhlAn_GTDB = "#4477AA",
  MetaPhlAn_default = "#66CCEE",
  ganon2 = "#228833",
  Kraken2 = "#EE8866",
  Kraken2_conf10 = "#AA3377"
)

metric_labels <- c(
  Richness = "Observed richness",
  Shannon = "Shannon index",
  Gini_Simpson = "Gini-Simpson index"
)

group_labels <- c(
  SW_NEG = "SW MAP-negative",
  SW_SUB = "SW subclinical",
  SW_POS = "SW clinical",
  ZA_NEG = "ZA MAP-negative",
  ZA_SUB = "ZA subclinical"
)

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x <- str_remove(x, "_3\\.2g.*$")
  x <- str_remove(x, "_metaphlan.*$")
  x <- str_remove(x, "\\.fastq(\\.gz)?$")
  x
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, p))
}

safe_range <- function(x, which = "min") {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  if (which == "min") min(x) else max(x)
}

safe_wilcox <- function(x, y, paired = FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (paired) {
    valid <- is.finite(x) & is.finite(y)
    x <- x[valid]
    y <- y[valid]

    if (length(x) < 2) return(NA_real_)
    if (all(x == y)) return(1)

  } else {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    if (length(x) < 2 || length(y) < 2) {
      return(NA_real_)
    }

    if (length(unique(c(x, y))) == 1) {
      return(1)
    }
  }

  suppressWarnings(
    wilcox.test(
      x, y,
      paired = paired,
      exact = FALSE
    )$p.value
  )
}

safe_kw <- function(d, group_col) {
  d <- d %>%
    filter(is.finite(value)) %>%
    mutate(
      Test_group = as.character(.data[[group_col]])
    ) %>%
    filter(
      !is.na(Test_group),
      nzchar(Test_group)
    )

  counts <- table(d$Test_group)

  if (
    length(counts) < 2 ||
    any(counts < 2) ||
    n_distinct(d$value) < 2
  ) {
    return(tibble(
      H = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      N = nrow(d),
      N_groups = length(counts)
    ))
  }

  test <- kruskal.test(value ~ Test_group, data = d)

  tibble(
    H = unname(test$statistic),
    df = unname(test$parameter),
    p_value = test$p.value,
    N = nrow(d),
    N_groups = length(counts)
  )
}

pairwise_tests <- function(d, group_col) {
  d <- d %>%
    filter(is.finite(value)) %>%
    mutate(
      Test_group = as.character(.data[[group_col]])
    ) %>%
    filter(!is.na(Test_group))

  groups <- unique(d$Test_group)

  if (length(groups) < 2) {
    return(tibble())
  }

  pairs <- combn(groups, 2, simplify = FALSE)

  map_dfr(pairs, function(pair) {
    x <- d$value[d$Test_group == pair[1]]
    y <- d$value[d$Test_group == pair[2]]

    tibble(
      Group_1 = pair[1],
      Group_2 = pair[2],
      N_1 = length(x),
      N_2 = length(y),
      Median_1 = safe_median(x),
      Median_2 = safe_median(y),
      p_value = safe_wilcox(x, y)
    )
  }) %>%
    mutate(
      p_BH = p.adjust(p_value, method = "BH")
    )
}

meta_raw <- read_tsv(
  metadata_file,
  col_types = cols(.default = col_character()),
  show_col_types = FALSE
)

id_candidates <- c(
  "SAMPLE_ID",
  "Sample_ID_Neu",
  "sample_id"
)

id_col <- id_candidates[
  id_candidates %in% names(meta_raw)
][1]

if (is.na(id_col)) {
  stop("No recognized sample-ID column in metadata.")
}

if (!"GROUP" %in% names(meta_raw)) {
  stop("Metadata must contain GROUP.")
}

farm_col <- names(meta_raw)[
  tolower(trimws(names(meta_raw))) == "farm number"
][1]

if (is.na(farm_col)) {
  stop(
    "Metadata must contain the column 'Farm number'."
  )
}

meta <- meta_raw %>%
  mutate(
    SAMPLE_ID = clean_id(.data[[id_col]]),
    GROUP = trimws(GROUP),

    Farm_number = trimws(
      as.character(.data[[farm_col]])
    ),

    Cohort_group = case_when(
      str_detect(SAMPLE_ID, "^SW_NEG_") ~ "SW_NEG",
      str_detect(SAMPLE_ID, "^SW_SUB_") ~ "SW_SUB",
      str_detect(SAMPLE_ID, "^SW_POS_") ~ "SW_POS",
      str_detect(SAMPLE_ID, "^ZA_NEG_") ~ "ZA_NEG",
      str_detect(SAMPLE_ID, "^ZA_SUB_") ~ "ZA_SUB",
      TRUE ~ NA_character_
    ),

    Country = case_when(
      str_starts(SAMPLE_ID, "SW_") ~ "Switzerland",
      str_starts(SAMPLE_ID, "ZA_") ~ "South Africa",
      TRUE ~ NA_character_
    ),

    Exclusion_reason = case_when(
      SAMPLE_ID %in% excluded_ids ~
        "Insufficient sequencing depth",

      str_to_lower(GROUP) == "map detection" ~
        "MAP detection only",

      str_to_lower(GROUP) %in% c("ex", "excluded") ~
        "Excluded",

      is.na(Cohort_group) ~
        "Not a microbiome sample",

      TRUE ~ NA_character_
    ),

    Farm_label = if_else(
      !is.na(Farm_number) &
        nzchar(Farm_number),
      paste("Farm", Farm_number),
      NA_character_
    )
  )

if (anyDuplicated(meta$SAMPLE_ID)) {
  stop("Duplicate SAMPLE_ID in metadata.")
}

meta_included <- meta %>%
  filter(is.na(Exclusion_reason)) %>%
  mutate(
    Cohort_group = factor(
      Cohort_group,
      levels = group_order
    )
  ) %>%
  arrange(Cohort_group, SAMPLE_ID)

cat("\nMETADATA SAMPLE COUNTS:\n")

print(
  meta_included %>%
    count(Cohort_group, .drop = FALSE)
)

read_profile <- function(method, rank) {

  path <- taxonomy_file(method, rank, "normalized")

  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }

  raw <- read_tsv(
    path,
    show_col_types = FALSE
  )

  if (!"clade_name" %in% names(raw)) {
    stop("Missing clade_name: ", path)
  }

  sample_cols <- setdiff(
    names(raw),
    "clade_name"
  )

  # Use only eligible metadata sample IDs.

  matched <- intersect(
    sample_cols,
    meta_included$SAMPLE_ID
  )

  if (!length(matched)) {
    stop("No matching samples in: ", path)
  }

  raw %>%
    select(
      clade_name,
      all_of(matched)
    ) %>%
    pivot_longer(
      -clade_name,
      names_to = "SAMPLE_ID",
      values_to = "Abundance"
    ) %>%
    mutate(
      Abundance = as.numeric(Abundance),

      Taxon = str_remove(
        clade_name,
        "^.*[;|]"
      ),

      Taxon = str_remove(
        Taxon,
        "^[dkpcofgst]__"
      ),

      Taxon = trimws(Taxon)
    ) %>%
    filter(
      !is.na(Taxon),
      nzchar(Taxon),

      !str_detect(
        str_to_lower(Taxon),
        "^(unclassified|unassigned|unknown|unidentified|na|none)$"
      ),

      !str_detect(
        clade_name,
        "(^|[;|])t__"
      )
    ) %>%
    group_by(
      SAMPLE_ID,
      Taxon
    ) %>%
    summarise(
      Abundance = sum(
        Abundance,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    mutate(
      Method = method,
      Rank = rank
    )
}

abundance_list <- list()

for (method in methods) {
  for (rank in ranks) {

    key <- paste(method, rank, sep = "__")

    abundance_list[[key]] <- read_profile(
      method,
      rank
    )

    cat("Loaded:", key, "\n")
  }
}

abundance_all <- bind_rows(abundance_list)

if (any(
  !is.finite(abundance_all$Abundance) |
  abundance_all$Abundance < 0
)) {
  stop("Invalid abundance values detected.")
}

calculate_alpha <- function(x) {

  richness <- sum(x > 0)
  total <- sum(x)

  if (total <= 0) {
    return(tibble(
      Richness = 0,
      Shannon = NA_real_,
      Gini_Simpson = NA_real_
    ))
  }

  p <- x / total
  p <- p[p > 0]

  tibble(
    Richness = richness,
    Shannon = -sum(p * log(p)),
    Gini_Simpson = 1 - sum(p^2)
  )
}

# All eligible metadata samples are included in QC,

sample_grid <- expand_grid(
  Method = methods,
  Rank = ranks,
  SAMPLE_ID = meta_included$SAMPLE_ID
)

abundance_qc <- abundance_all %>%
  group_by(Method, Rank, SAMPLE_ID) %>%
  summarise(
    Total_before = sum(Abundance),
    N_taxa_detected = sum(Abundance > 0),
    .groups = "drop"
  )

abundance_qc <- sample_grid %>%
  left_join(
    abundance_qc,
    by = c("Method", "Rank", "SAMPLE_ID")
  ) %>%
  mutate(
    Profile_status = case_when(
      is.na(Total_before) ~ "Missing",
      Total_before <= 0 ~ "Zero abundance",
      TRUE ~ "Valid"
    )
  ) %>%
  left_join(
    meta_included %>%
      select(
        SAMPLE_ID,
        Cohort_group,
        Country,
        Farm_number,
        Farm_label
      ),
    by = "SAMPLE_ID"
  )

cat("\nPROFILE QC:\n")

print(
  abundance_qc %>%
    count(Method, Rank, Profile_status),
  n = Inf
)

# Missing and zero profiles are not silently assigned

alpha_all <- abundance_all %>%
  group_by(
    Method,
    Rank,
    SAMPLE_ID
  ) %>%
  group_modify(
    ~ calculate_alpha(.x$Abundance)
  ) %>%
  ungroup() %>%
  left_join(
    meta_included %>%
      select(
        SAMPLE_ID,
        Cohort_group,
        Country,
        Farm_number,
        Farm_label
      ),
    by = "SAMPLE_ID"
  ) %>%
  mutate(
    Cohort_group = factor(
      Cohort_group,
      levels = group_order
    )
  )

alpha_all <- alpha_all %>%
  left_join(
    abundance_qc %>%
      select(
        Method,
        Rank,
        SAMPLE_ID,
        Profile_status
      ),
    by = c("Method", "Rank", "SAMPLE_ID")
  ) %>%
  mutate(
    across(
      all_of(metrics),
      ~ if_else(
        Profile_status == "Valid",
        .x,
        NA_real_
      )
    )
  )

alpha_long <- alpha_all %>%
  pivot_longer(
    cols = all_of(metrics),
    names_to = "Metric",
    values_to = "value"
  )

primary <- alpha_long %>%
  filter(
    Method == primary_method,
    is.finite(value)
  )

# metadata samples at BOTH taxonomic ranks.

primary_qc <- abundance_qc %>%
  filter(Method == primary_method)

if (any(primary_qc$Profile_status != "Valid")) {
  print(
    primary_qc %>%
      filter(Profile_status != "Valid"),
    n = Inf
  )

  stop(
    "Some MetaPhlAn GTDB profiles are missing or zero. ",
    "Review sample QC before biological analysis."
  )
}

primary_descriptive <- list()
primary_overall <- list()
primary_pairwise <- list()

for (cohort in names(cohorts)) {

  d <- primary %>%
    filter(
      as.character(Cohort_group) %in%
        cohorts[[cohort]]
    )

  primary_descriptive[[cohort]] <- d %>%
    group_by(
      Rank,
      Metric,
      Cohort_group
    ) %>%
    summarise(
      N = n(),
      Mean = mean(value),
      SD = safe_sd(value),
      Median = median(value),
      Q1 = safe_quantile(value, 0.25),
      Q3 = safe_quantile(value, 0.75),
      Minimum = min(value),
      Maximum = max(value),
      .groups = "drop"
    ) %>%
    mutate(Cohort = cohort)

  if (cohort == "South_Africa") {

    primary_overall[[cohort]] <- d %>%
      group_by(Rank, Metric) %>%
      group_modify(~ {
        x <- .x$value[
          .x$Cohort_group == "ZA_NEG"
        ]

        y <- .x$value[
          .x$Cohort_group == "ZA_SUB"
        ]

        tibble(
          Test = "Wilcoxon",
          H = NA_real_,
          df = NA_real_,
          p_value = safe_wilcox(x, y),
          N = length(x) + length(y)
        )
      }) %>%
      ungroup()

  } else {

    primary_overall[[cohort]] <- d %>%
      group_by(Rank, Metric) %>%
      group_modify(
        ~ safe_kw(.x, "Cohort_group")
      ) %>%
      ungroup() %>%
      mutate(Test = "Kruskal-Wallis")
  }

  primary_overall[[cohort]] <-
    primary_overall[[cohort]] %>%
    mutate(
      Cohort = cohort,
      p_BH = p.adjust(
        p_value,
        method = "BH"
      )
    )

  primary_pairwise[[cohort]] <- d %>%
    group_by(Rank, Metric) %>%
    group_modify(
      ~ pairwise_tests(.x, "Cohort_group")
    ) %>%
    ungroup() %>%
    mutate(Cohort = cohort)
}

primary_descriptive <- bind_rows(
  primary_descriptive
)

primary_overall <- bind_rows(
  primary_overall
)

primary_pairwise <- bind_rows(
  primary_pairwise
)

swiss_primary <- primary %>%
  filter(Country == "Switzerland")

# Determine eligible farms using UNIQUE animals,

farm_counts <- meta_included %>%
  filter(Country == "Switzerland") %>%
  count(
    Cohort_group,
    Farm_number,
    Farm_label,
    name = "N_animals"
  ) %>%
  mutate(
    Eligible = !is.na(Farm_number) &
      nzchar(Farm_number) &
      N_animals >= min_farm_n
  )

cat("\nSWISS FARM COUNTS:\n")

print(farm_counts, n = Inf)

eligible_farms <- farm_counts %>%
  filter(Eligible) %>%
  pull(Farm_number) %>%
  unique()

farm_data <- swiss_primary %>%
  filter(Farm_number %in% eligible_farms)

# Clinical singleton farms are excluded automatically.

farm_overall <- farm_data %>%
  group_by(Rank, Metric) %>%
  group_modify(
    ~ safe_kw(.x, "Farm_number")
  ) %>%
  ungroup() %>%
  mutate(
    Analysis = "Swiss farms overall",
    p_BH = p.adjust(
      p_value,
      method = "BH"
    )
  )

farm_pairwise <- farm_data %>%
  group_by(Rank, Metric) %>%
  group_modify(
    ~ pairwise_tests(.x, "Farm_number")
  ) %>%
  ungroup() %>%
  mutate(
    Analysis = "Swiss farms overall"
  )

farm_within_overall <- list()
farm_within_pairwise <- list()

for (grp in c("SW_NEG", "SW_SUB")) {

  d <- farm_data %>%
    filter(Cohort_group == grp)

  farm_within_overall[[grp]] <- d %>%
    group_by(Rank, Metric) %>%
    group_modify(
      ~ safe_kw(.x, "Farm_number")
    ) %>%
    ungroup() %>%
    mutate(
      MAP_group = grp,
      p_BH = p.adjust(
        p_value,
        method = "BH"
      )
    )

  farm_within_pairwise[[grp]] <- d %>%
    group_by(Rank, Metric) %>%
    group_modify(
      ~ pairwise_tests(.x, "Farm_number")
    ) %>%
    ungroup() %>%
    mutate(MAP_group = grp)
}

farm_within_overall <- bind_rows(
  farm_within_overall
)

farm_within_pairwise <- bind_rows(
  farm_within_pairwise
)

# Use the same Swiss animals across all five methods

swiss_ids <- meta_included %>%
  filter(Country == "Switzerland") %>%
  pull(SAMPLE_ID)

valid_method_ids <- abundance_qc %>%
  filter(
    SAMPLE_ID %in% swiss_ids,
    Profile_status == "Valid"
  ) %>%
  count(SAMPLE_ID, name = "N_valid") %>%
  filter(
    N_valid == length(methods) * length(ranks)
  ) %>%
  pull(SAMPLE_ID)

cat(
  "\nFIVE-METHOD SHARED SWISS SAMPLES:",
  length(valid_method_ids),
  "\n"
)

method_data <- alpha_long %>%
  filter(
    SAMPLE_ID %in% valid_method_ids,
    is.finite(value)
  )

method_descriptive <- method_data %>%
  group_by(
    Method,
    Rank,
    Metric
  ) %>%
  summarise(
    N = n(),
    Mean = mean(value),
    SD = safe_sd(value),
    Median = median(value),
    Q1 = safe_quantile(value, 0.25),
    Q3 = safe_quantile(value, 0.75),
    Minimum = min(value),
    Maximum = max(value),
    .groups = "drop"
  )

method_MAP_KW <- method_data %>%
  group_by(
    Method,
    Rank,
    Metric
  ) %>%
  group_modify(
    ~ safe_kw(.x, "Cohort_group")
  ) %>%
  ungroup() %>%
  group_by(Rank, Metric) %>%
  mutate(
    p_BH_methods = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  mutate(
    p_BH_all30 = p.adjust(
      p_value,
      method = "BH"
    )
  )

method_MAP_pairwise <- method_data %>%
  group_by(
    Method,
    Rank,
    Metric
  ) %>%
  group_modify(
    ~ pairwise_tests(.x, "Cohort_group")
  ) %>%
  ungroup()

method_pairs <- combn(
  methods,
  2,
  simplify = FALSE
)

paired_method_results <- list()

for (rank in ranks) {
  for (metric in metrics) {

    d <- method_data %>%
      filter(
        Rank == rank,
        Metric == metric
      )

    wide <- d %>%
      select(
        SAMPLE_ID,
        Method,
        value
      ) %>%
      pivot_wider(
        names_from = Method,
        values_from = value
      )

    results <- map_dfr(
      method_pairs,
      function(pair) {

        x <- wide[[pair[1]]]
        y <- wide[[pair[2]]]

        valid <- is.finite(x) & is.finite(y)

        tibble(
          Method_1 = pair[1],
          Method_2 = pair[2],
          N_pairs = sum(valid),
          Median_difference =
            safe_median(x[valid] - y[valid]),
          p_value = safe_wilcox(
            x, y,
            paired = TRUE
          )
        )
      }
    ) %>%
      mutate(
        Rank = rank,
        Metric = metric,
        p_BH = p.adjust(
          p_value,
          method = "BH"
        )
      )

    paired_method_results[[
      paste(rank, metric, sep = "__")
    ]] <- results
  }
}

paired_method_results <- bind_rows(
  paired_method_results
)

alpha_theme <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1,
      size = 10
    ),
    strip.text = element_text(
      face = "bold"
    ),
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 13
    )
  )

make_alpha_plot <- function(
    data,
    rank,
    xvar,
    colors,
    title,
    filename,
    xlabels = NULL
) {

  d <- data %>%
    filter(
      Rank == rank,
      is.finite(value)
    ) %>%
    mutate(
      Metric = factor(
        Metric,
        levels = metrics
      )
    )

  if (!nrow(d)) return(invisible(NULL))

  p <- ggplot(
    d,
    aes(
      x = .data[[xvar]],
      y = value,
      fill = .data[[xvar]]
    )
  ) +
    geom_boxplot(
      width = 0.65,
      alpha = 0.65,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      alpha = 0.80,
      size = 1.8,
      shape = 21,
      stroke = 0.25
    ) +
    facet_wrap(
      ~ Metric,
      nrow = 1,
      scales = "free_y",
      labeller = labeller(
        Metric = metric_labels
      )
    ) +
    scale_fill_manual(
      values = colors
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    alpha_theme

  if (!is.null(xlabels)) {
    p <- p +
      scale_x_discrete(
        labels = xlabels
      )
  }

  dir.create(
    dirname(filename),
    recursive = TRUE,
    showWarnings = FALSE
  )

  ggsave(
    filename,
    p,
    width = 12,
    height = 4.5,
    dpi = 300,
    bg = "white"
  )

  for (metric in metrics) {

    single <- d %>%
      filter(Metric == metric)

    p_single <- ggplot(
      single,
      aes(
        x = .data[[xvar]],
        y = value,
        fill = .data[[xvar]]
      )
    ) +
      geom_boxplot(
        width = 0.65,
        alpha = 0.65,
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.10,
        alpha = 0.8,
        size = 1.8,
        shape = 21,
        stroke = 0.25
      ) +
      scale_fill_manual(
        values = colors
      ) +
      labs(
        title = paste(title, metric),
        x = NULL,
        y = metric_labels[[metric]]
      ) +
      alpha_theme

    if (!is.null(xlabels)) {
      p_single <- p_single +
        scale_x_discrete(labels = xlabels)
    }

    ggsave(
      str_replace(
        filename,
        "\\.png$",
        paste0("_", metric, ".png")
      ),
      p_single,
      width = 5.5,
      height = 4.5,
      dpi = 300,
      bg = "white"
    )
  }

  invisible(p)
}

for (cohort in names(cohorts)) {

  d <- primary %>%
    filter(
      as.character(Cohort_group) %in%
        cohorts[[cohort]]
    ) %>%
    mutate(
      Cohort_group = factor(
        as.character(Cohort_group),
        levels = cohorts[[cohort]]
      )
    )

  for (rank in ranks) {

    make_alpha_plot(
      data = d,
      rank = rank,
      xvar = "Cohort_group",
      colors = group_colors,
      title = paste(
        cohort,
        rank,
        "MetaPhlAn GTDB",
        sep = " | "
      ),
      filename = file.path(
        output_dir,
        "Figures_MetaPhlAn_GTDB",
        cohort,
        paste0(
          cohort,
          "_",
          rank,
          "_Alpha.png"
        )
      ),
      xlabels = group_labels
    )
  }
}

farm_levels <- farm_counts %>%
  filter(Eligible) %>%
  distinct(Farm_number, Farm_label) %>%
  mutate(
    Farm_numeric = suppressWarnings(
      as.numeric(Farm_number)
    )
  ) %>%
  arrange(Farm_numeric, Farm_number) %>%
  pull(Farm_label)

farm_plot_data <- farm_data %>%
  mutate(
    Farm_label = factor(
      Farm_label,
      levels = farm_levels
    )
  )

farm_palette <- setNames(
  grDevices::hcl.colors(
    length(farm_levels),
    palette = "Dark 3"
  ),
  farm_levels
)

for (rank in ranks) {

  make_alpha_plot(
    data = farm_plot_data,
    rank = rank,
    xvar = "Farm_label",
    colors = farm_palette,
    title = paste(
      "Swiss farm comparison |",
      rank
    ),
    filename = file.path(
      output_dir,
      "Figures_MetaPhlAn_GTDB",
      "Swiss_Farms",
      paste0(
        "Swiss_Farms_",
        rank,
        ".png"
      )
    )
  )
}

for (grp in c("SW_NEG", "SW_SUB")) {

  d <- farm_plot_data %>%
    filter(Cohort_group == grp) %>%
    mutate(
      Farm_label = droplevels(Farm_label)
    )

  for (rank in ranks) {

    make_alpha_plot(
      data = d,
      rank = rank,
      xvar = "Farm_label",
      colors = farm_palette,
      title = paste(
        grp,
        "farm comparison |",
        rank
      ),
      filename = file.path(
        output_dir,
        "Figures_MetaPhlAn_GTDB",
        "Swiss_Farms",
        paste0(
          grp,
          "_Farms_",
          rank,
          ".png"
        )
      )
    )
  }
}

method_plot_data <- method_data %>%
  mutate(
    Method = factor(
      Method,
      levels = methods
    )
  )

for (rank in ranks) {

  make_alpha_plot(
    data = method_plot_data,
    rank = rank,
    xvar = "Method",
    colors = method_colors,
    title = paste(
      "Swiss five-method comparison |",
      rank
    ),
    filename = file.path(
      output_dir,
      "Figures_Method_Comparison",
      paste0(
        "Methods_",
        rank,
        ".png"
      )
    )
  )
}

for (rank in ranks) {

  d <- method_data %>%
    filter(
      Rank == rank,
      is.finite(value)
    ) %>%
    mutate(
      Method = factor(
        Method,
        levels = methods
      ),

      Cohort_group = factor(
        Cohort_group,
        levels = group_order[1:3]
      ),

      Metric = factor(
        Metric,
        levels = metrics
      )
    )

  p <- ggplot(
    d,
    aes(
      x = Cohort_group,
      y = value,
      fill = Cohort_group
    )
  ) +
    geom_boxplot(
      width = 0.65,
      alpha = 0.65,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.10,
      alpha = 0.7,
      size = 1.1
    ) +
    facet_grid(
      Metric ~ Method,
      scales = "free_y",
      labeller = labeller(
        Metric = metric_labels
      )
    ) +
    scale_fill_manual(
      values = group_colors
    ) +
    scale_x_discrete(
      labels = c(
        SW_NEG = "Negative",
        SW_SUB = "Subclinical",
        SW_POS = "Clinical"
      )
    ) +
    labs(
      title = paste(
        "Alpha diversity across methods |",
        rank
      ),
      x = NULL,
      y = NULL
    ) +
    alpha_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        size = 8
      )
    )

  ggsave(
    file.path(
      output_dir,
      "Figures_Method_Comparison",
      paste0(
        "Methods_MAP_",
        rank,
        ".png"
      )
    ),
    p,
    width = 17,
    height = 10,
    dpi = 300,
    bg = "white"
  )
}

sample_inclusion <- meta %>%
  select(
    SAMPLE_ID,
    GROUP,
    Cohort_group,
    Country,
    Farm_number,
    Farm_label,
    Exclusion_reason
  ) %>%
  mutate(
    Included = is.na(Exclusion_reason)
  )

sample_counts <- meta_included %>%
  count(
    Country,
    Cohort_group,
    name = "N_samples"
  )

outputs <- list(
  "Sample inclusion" = sample_inclusion,
  "Sample counts" = sample_counts,
  "Sample QC" = abundance_qc,
  "Alpha values" = alpha_all,

  "Primary descriptive" = primary_descriptive,
  "Primary overall tests" = primary_overall,
  "Primary pairwise" = primary_pairwise,

  "Farm counts" = farm_counts,
  "Farm overall" = farm_overall,
  "Farm pairwise" = farm_pairwise,
  "Farm within overall" = farm_within_overall,
  "Farm within pairwise" = farm_within_pairwise,

  "Method descriptive" = method_descriptive,
  "Method MAP KW" = method_MAP_KW,
  "Method MAP pairwise" = method_MAP_pairwise,
  "Paired methods" = paired_method_results
)

wb <- createWorkbook()

header_style <- createStyle(
  fontColour = "#FFFFFF",
  fgFill = "#284B63",
  textDecoration = "bold"
)

for (sheet in names(outputs)) {

  addWorksheet(wb, sheet)

  dat <- as.data.frame(outputs[[sheet]])

  writeData(
    wb,
    sheet,
    dat,
    headerStyle = header_style
  )

  freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )

  if (ncol(dat) > 0) {
    setColWidths(
      wb,
      sheet,
      cols = seq_len(ncol(dat)),
      widths = "auto"
    )
  }

  if (nrow(dat) > 0 && ncol(dat) > 0) {
    addFilter(
      wb,
      sheet,
      rows = 1,
      cols = seq_len(ncol(dat))
    )
  }
}

saveWorkbook(
  wb,
  file.path(
    output_dir,
    "DVM_Alpha_Diversity_Master.xlsx"
  ),
  overwrite = TRUE
)

write_csv(
  sample_inclusion,
  file.path(
    output_dir,
    "Sample_inclusion_QC.csv"
  )
)

cat("\n========================================\n")
cat("DVM ALPHA DIVERSITY COMPLETE\n")
cat("========================================\n")

cat("\nSAMPLE COUNTS:\n")
print(sample_counts, n = Inf)

cat("\nPRIMARY OVERALL TESTS:\n")
print(primary_overall, n = Inf)

cat("\nSWISS FARM TESTS:\n")
print(farm_overall, n = Inf)

cat("\nFIVE-METHOD MAP TESTS:\n")
print(method_MAP_KW, n = Inf)

cat("\nSHARED SWISS SAMPLES FOR METHODS:\n")
cat(length(valid_method_ids), "\n")

cat("\nOUTPUT DIRECTORY:\n")
cat(normalizePath(output_dir), "\n")
