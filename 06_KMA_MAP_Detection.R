# 06_KMA_MAP_Detection.R: IS900/F57 marker detection and statistical analyses
# See README.md for inputs, provenance and validation status.

# Marker evidence and detection predictors
local({
  packages <- c(
    "readxl",
    "dplyr",
    "tidyr",
    "ggplot2",
    "readr",
    "broom",
    "patchwork"
  )

  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1),
            quietly = TRUE)
  ]

  if (length(missing_packages) > 0) stop("Required packages missing: ", paste(missing_packages, collapse = ", "))

  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(broom)
    library(patchwork)
  })

  metadata_file <- file.path("data", "metadata.txt")
  kma_file <- file.path("data", "map_detection", "kma_IS900_F57_neu.xlsx")
  s1_file <- file.path("data", "map_detection", "S1_Sample_Master.xlsx")
  output_dir <- "MAP_detection_analysis"

  for (input in c(metadata_file, kma_file, s1_file)) {
    if (!file.exists(input)) {
      stop("Missing required input: ", input,
           "\nS1_Sample_Master.xlsx is Supplementary Table S1 from the publication; ",
           "place a copy under data/map_detection/ after downloading it from the journal.")
    }
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  to_numeric <- function(x) {
    suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = TRUE)))
  }

  is_positive <- function(x) {
    tolower(trimws(as.character(x))) %in% c("yes", "positive", "true", "1")
  }

  metadata <- readr::read_tsv(metadata_file,
                              col_types = readr::cols(.default = readr::col_character()),
                              show_col_types = FALSE)
  if (!all(c("SAMPLE_ID", "Animal_ID") %in% names(metadata))) {
    stop("Metadata must contain SAMPLE_ID and Animal_ID.")
  }
  id_lookup <- dplyr::bind_rows(
    metadata %>% transmute(raw_id = trimws(Animal_ID), id = trimws(SAMPLE_ID)),
    metadata %>% transmute(raw_id = trimws(SAMPLE_ID), id = trimws(SAMPLE_ID))
  ) %>%
    filter(!is.na(raw_id), nzchar(raw_id)) %>%
    distinct(raw_id, .keep_all = TRUE)
  if (anyDuplicated(id_lookup$raw_id)) stop("Duplicate raw KMA identifiers in metadata.")

  # S1 is published separately with the article, not part of the Zenodo archive.
  s1 <- readxl::read_excel(s1_file, sheet = 1, skip = 4)
  s1_fields <- c("Sample ID", "Country of origin", "Study group",
                 "qPCR result", "qPCR Ct", "Culture result",
                 "Host-depleted yield (Gbp)")
  if (!all(s1_fields %in% names(s1))) {
    stop("Published S1 has unexpected column names. Expected: ",
         paste(s1_fields, collapse = ", "))
  }
  samples <- s1 %>%
    transmute(
      id = trimws(as.character(`Sample ID`)),
      country = as.character(`Country of origin`),
      group = as.character(`Study group`),
      qpcr = as.character(`qPCR result`),
      ct = to_numeric(`qPCR Ct`),
      culture = as.character(`Culture result`),
      yield_gbp = to_numeric(`Host-depleted yield (Gbp)`)
    ) %>%
    filter(grepl("^(SW|ZA)_(NEG|SUB|POS)_[0-9]+$", id),
           !id %in% c("SW_POS_00", "SW_POS_01"))
  if (anyDuplicated(samples$id)) stop("Duplicate sample IDs in S1.")

  # Raw KMA spreadsheet consists of repeated per-sample reports.
  kma_columns <- c("result_file", "template", "score", "expected",
                   "template_length", "template_identity", "template_coverage",
                   "query_identity", "query_coverage", "depth", "q_value", "p_value")
  raw_kma <- readxl::read_excel(kma_file, sheet = 1, col_names = FALSE,
                                range = readxl::cell_cols("A:L"))
  names(raw_kma) <- kma_columns
  kma_hits <- raw_kma %>%
    filter(grepl("IS900|F57", as.character(template))) %>%
    mutate(
      raw_id = sub("_kma_IS900_F57.*$", "", as.character(result_file)),
      marker = case_when(grepl("IS900", template) ~ "is900",
                         grepl("F57", template) ~ "f57",
                         TRUE ~ NA_character_),
      breadth = to_numeric(template_coverage),
      depth = to_numeric(depth)
    ) %>%
    left_join(id_lookup, by = "raw_id")
  if (anyNA(kma_hits$id)) {
    stop("Raw KMA sample identifiers not matched in metadata: ",
         paste(unique(kma_hits$raw_id[is.na(kma_hits$id)]), collapse = ", "))
  }
  if (anyDuplicated(paste(kma_hits$id, kma_hits$marker, sep = "|"))) {
    stop("Duplicate sample-marker rows found in raw KMA workbook.")
  }
  if (any(!is.finite(kma_hits$breadth)) ||
      any(kma_hits$breadth < 0 | kma_hits$breadth > 100)) {
    stop("Invalid marker-reference coverage breadth values in raw KMA file.")
  }
  marker_wide <- kma_hits %>%
    select(id, marker, breadth, depth) %>%
    pivot_wider(names_from = marker, values_from = c(breadth, depth),
                names_glue = "{marker}_{.value}")

  if (any(!marker_wide$id %in% samples$id)) {
    warning("KMA evidence for samples not in S1 will be excluded: ",
            paste(setdiff(marker_wide$id, samples$id), collapse = ", "))
  }
  dat <- samples %>%
    left_join(marker_wide, by = "id") %>%
    mutate(
      is900_breadth = coalesce(is900_breadth, 0),
      f57_breadth = coalesce(f57_breadth, 0),
      is900_depth = coalesce(is900_depth, 0),
      f57_depth = coalesce(f57_depth, 0),
      is900 = is900_breadth >= 90,
      f57 = f57_breadth >= 90,
      either = is900 | f57,
      both = is900 & f57,
      conventional_positive = is_positive(qpcr) | is_positive(culture)
    )

  positive_s1 <- dat %>% filter(conventional_positive)
  cat("KMA threshold: >=90% of marker reference sequence covered.\n")
  cat("Conventionally positive (all available S1 samples): ",
      nrow(positive_s1), "; either marker: ", sum(positive_s1$either),
      "; both: ", sum(positive_s1$both), "\n", sep = "")
  if (nrow(positive_s1) != 34L || sum(positive_s1$either) != 21L ||
      sum(positive_s1$both) != 15L) {
    warning("KMA counts differ from reported 34 / 21 / 15. Check input versions and IDs.")
  }

  write_csv(
    dat,
    file.path(output_dir, "analysis_dataset.csv")
  )

  positive <- dat %>%
    filter(conventional_positive) %>%
    mutate(
      log10_yield = ifelse(
        yield_gbp > 0,
        log10(yield_gbp),
        NA_real_
      )
    )

  ctdata <- dat %>%
    filter(
      is_positive(qpcr),
      !is.na(ct)
    ) %>%
    mutate(
      log10_yield = ifelse(
        yield_gbp > 0,
        log10(yield_gbp),
        NA_real_
      ),
      detection_group = factor(
        either,
        levels = c(FALSE, TRUE),
        labels = c("Not detected", "Detected")
      )
    )

  cat("\nEligible samples:", nrow(dat))
  cat("\nConventionally MAP-positive:", nrow(positive))
  cat("\nSamples with available Ct:", nrow(ctdata), "\n")

  wilcox_comparison <- function(
      data, variable, detection, analysis_name
  ) {

    d <- data %>%
      filter(
        !is.na(.data[[variable]]),
        !is.na(.data[[detection]])
      )

    detected <- d[[variable]][d[[detection]]]
    undetected <- d[[variable]][!d[[detection]]]

    p <- NA_real_

    if (
      length(detected) >= 2 &&
      length(undetected) >= 2
    ) {
      p <- wilcox.test(
        detected,
        undetected,
        exact = FALSE
      )$p.value
    }

    tibble(
      analysis = analysis_name,
      n_detected = length(detected),
      n_undetected = length(undetected),
      median_detected = if (length(detected) > 0)
        median(detected) else NA_real_,
      median_undetected = if (length(undetected) > 0)
        median(undetected) else NA_real_,
      p_value = p
    )
  }

  wilcoxon_results <- bind_rows(
    wilcox_comparison(
      positive, "yield_gbp", "either",
      "Sequencing yield vs either marker"
    ),
    wilcox_comparison(
      positive, "yield_gbp", "both",
      "Sequencing yield vs both markers"
    ),
    wilcox_comparison(
      ctdata, "ct", "either",
      "Ct vs either marker"
    )
  )

  write_csv(
    wilcoxon_results,
    file.path(output_dir, "wilcoxon_results.csv")
  )

  spearman_analysis <- function(
      data, x, y, analysis_name
  ) {

    d <- data %>%
      filter(
        !is.na(.data[[x]]),
        !is.na(.data[[y]])
      )

    if (
      nrow(d) < 4 ||
      sd(d[[x]]) == 0 ||
      sd(d[[y]]) == 0
    ) {
      return(tibble(
        analysis = analysis_name,
        n = nrow(d),
        rho = NA_real_,
        p_value = NA_real_
      ))
    }

    test <- suppressWarnings(
      cor.test(
        d[[x]],
        d[[y]],
        method = "spearman",
        exact = FALSE
      )
    )

    tibble(
      analysis = analysis_name,
      n = nrow(d),
      rho = unname(test$estimate),
      p_value = test$p.value
    )
  }

  markers <- c(
    "is900_depth",
    "f57_depth",
    "is900_breadth",
    "f57_breadth"
  )

  spearman_results <- bind_rows(
    lapply(markers, function(marker) {
      spearman_analysis(
        positive,
        "yield_gbp",
        marker,
        paste("Sequencing yield vs", marker)
      )
    })
  )

  write_csv(
    spearman_results,
    file.path(output_dir, "spearman_results.csv")
  )

  fit_logistic <- function(
      data, formula, model_name
  ) {

    variables <- all.vars(formula)

    d <- data[
      complete.cases(data[, variables, drop = FALSE]),
      ,
      drop = FALSE
    ]

    if (
      nrow(d) < 10 ||
      length(unique(d$either)) < 2
    ) {
      return(tibble(
        model = model_name,
        term = "Insufficient data",
        n = nrow(d),
        OR = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_value = NA_real_
      ))
    }

    model <- glm(
      formula,
      data = d,
      family = binomial()
    )

    broom::tidy(model) %>%
      mutate(
        model = model_name,
        n = nrow(d),
        OR = exp(estimate),
        CI_low = exp(
          estimate - 1.96 * std.error
        ),
        CI_high = exp(
          estimate + 1.96 * std.error
        ),
        p_value = p.value
      ) %>%
      select(
        model, term, n,
        OR, CI_low, CI_high, p_value
      )
  }

  logistic_results <- bind_rows(
    fit_logistic(
      positive,
      either ~ log10_yield,
      "Sequencing yield only"
    ),
    fit_logistic(
      ctdata,
      either ~ ct,
      "Ct only"
    ),
    fit_logistic(
      ctdata,
      either ~ log10_yield + ct,
      "Sequencing yield + Ct"
    )
  )

  write_csv(
    logistic_results,
    file.path(
      output_dir,
      "logistic_regression_results.csv"
    )
  )

  theme_publication <- function() {
    theme_classic(base_size = 12) +
      theme(
        axis.title = element_text(
          size = 12,
          colour = "black"
        ),
        axis.text = element_text(
          size = 11,
          colour = "black"
        ),
        plot.tag = element_text(
          size = 15,
          face = "bold"
        ),
        plot.margin = margin(10, 12, 10, 10)
      )
  }

  p_yield <- ggplot(
    positive,
    aes(
      x = factor(
        either,
        levels = c(FALSE, TRUE),
        labels = c("Not detected", "Detected")
      ),
      y = yield_gbp
    )
  ) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      fill = "grey88"
    ) +
    geom_jitter(
      width = 0.12,
      height = 0,
      size = 2,
      alpha = 0.75
    ) +
    labs(
      x = "Metagenomic MAP detection",
      y = "Sequencing yield (Gbp)"
    ) +
    theme_publication()

  ggsave(
    file.path(output_dir, "yield_by_detection.png"),
    p_yield,
    width = 6,
    height = 4.5,
    dpi = 600
  )

  if (
    nrow(ctdata) >= 10 &&
    length(unique(ctdata$either)) == 2
  ) {

    ct_counts <- ctdata %>%
      count(detection_group)

    ct_p <- wilcox.test(
      ct ~ detection_group,
      data = ctdata,
      exact = FALSE
    )$p.value

    ct_p_label <- paste0(
      "Wilcoxon p = ",
      format.pval(
        ct_p,
        digits = 2,
        eps = 0.001
      )
    )

    ct_upper <- max(ctdata$ct, na.rm = TRUE)

    p_ct_box <- ggplot(
      ctdata,
      aes(
        x = detection_group,
        y = ct
      )
    ) +
      geom_boxplot(
        width = 0.48,
        outlier.shape = NA,
        fill = "grey88",
        colour = "black",
        linewidth = 0.5
      ) +
      geom_jitter(
        width = 0.12,
        height = 0,
        size = 2,
        alpha = 0.75
      ) +
      geom_text(
        data = ct_counts,
        aes(
          x = detection_group,
          y = ct_upper + 2,
          label = paste0("n = ", n)
        ),
        inherit.aes = FALSE,
        size = 3.5
      ) +
      annotate(
        "text",
        x = 1.5,
        y = ct_upper + 4,
        label = ct_p_label,
        size = 3.5
      ) +
      scale_y_continuous(
        expand = expansion(
          mult = c(0.04, 0.18)
        )
      ) +
      labs(
        x = "Metagenomic MAP detection",
        y = "qPCR Ct"
      ) +
      theme_publication()

    ct_model <- glm(
      either ~ ct,
      data = ctdata,
      family = binomial()
    )

    prediction_grid <- data.frame(
      ct = seq(
        min(ctdata$ct),
        max(ctdata$ct),
        length.out = 300
      )
    )

    predictions <- predict(
      ct_model,
      newdata = prediction_grid,
      type = "link",
      se.fit = TRUE
    )

    prediction_grid <- prediction_grid %>%
      mutate(
        probability = plogis(
          predictions$fit
        ),
        lower_95 = plogis(
          predictions$fit -
            1.96 * predictions$se.fit
        ),
        upper_95 = plogis(
          predictions$fit +
            1.96 * predictions$se.fit
        )
      )

    set.seed(123)

    p_ct_logistic <- ggplot(
      ctdata,
      aes(
        x = ct,
        y = as.integer(either)
      )
    ) +
      geom_ribbon(
        data = prediction_grid,
        aes(
          x = ct,
          ymin = lower_95,
          ymax = upper_95
        ),
        inherit.aes = FALSE,
        fill = "grey70",
        alpha = 0.45
      ) +
      geom_line(
        data = prediction_grid,
        aes(
          x = ct,
          y = probability
        ),
        inherit.aes = FALSE,
        linewidth = 0.9,
        colour = "black"
      ) +
      geom_jitter(
        width = 0,
        height = 0.025,
        size = 2,
        alpha = 0.75
      ) +
      scale_y_continuous(
        limits = c(-0.06, 1.06),
        breaks = seq(0, 1, 0.2),
        labels = function(x) {
          paste0(round(x * 100), "%")
        }
      ) +
      coord_cartesian(
        ylim = c(-0.06, 1.06)
      ) +
      labs(
        x = "qPCR Ct",
        y = "Probability of MAP marker detection"
      ) +
      theme_publication()

    combined_ct <- (
      p_ct_box + p_ct_logistic
    ) +
      plot_layout(
        widths = c(1, 1.15)
      ) +
      plot_annotation(
        tag_levels = "A"
      )

    ggsave(
      file.path(
        output_dir,
        "Figure_Ct_MAP_detection_combined.png"
      ),
      combined_ct,
      width = 11,
      height = 4.8,
      dpi = 600,
      bg = "white"
    )

    ggsave(
      file.path(
        output_dir,
        "Figure_Ct_MAP_detection_combined.pdf"
      ),
      combined_ct,
      width = 11,
      height = 4.8,
      bg = "white"
    )

    ggsave(
      file.path(output_dir, "ct_by_detection.png"),
      p_ct_box,
      width = 6,
      height = 5,
      dpi = 600
    )

    ggsave(
      file.path(
        output_dir,
        "ct_detection_probability.png"
      ),
      p_ct_logistic,
      width = 6,
      height = 5,
      dpi = 600
    )

    write_csv(
      prediction_grid,
      file.path(
        output_dir,
        "ct_predicted_probabilities.csv"
      )
    )

    print(combined_ct)
  }

  for (marker in markers) {

    p <- ggplot(
      positive,
      aes(
        x = yield_gbp,
        y = .data[[marker]]
      )
    ) +
      geom_point(
        size = 2,
        alpha = 0.75
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        linetype = "dashed",
        colour = "grey40"
      ) +
      labs(
        x = "Sequencing yield (Gbp)",
        y = marker
      ) +
      theme_publication()

    ggsave(
      file.path(
        output_dir,
        paste0("yield_vs_", marker, ".png")
      ),
      p,
      width = 6,
      height = 4.5,
      dpi = 600
    )
  }

  subgroup_counts <- dat %>%
    group_by(group, country) %>%
    summarise(
      n = n(),
      conventional_positive = sum(
        conventional_positive
      ),
      detected = sum(either, na.rm = TRUE),
      ct_available = sum(!is.na(ct)),
      .groups = "drop"
    )

  write_csv(
    subgroup_counts,
    file.path(output_dir, "subgroup_counts.csv")
  )

  cat("\n====================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("====================================\n")

  cat("\nWILCOXON RESULTS\n")
  print(wilcoxon_results)

  cat("\nSPEARMAN RESULTS\n")
  print(spearman_results)

  cat("\nLOGISTIC REGRESSION RESULTS\n")
  print(logistic_results)

  cat("\nSUBGROUP COUNTS\n")
  print(subgroup_counts)

  cat(
    "\nResults saved to:",
    normalizePath(output_dir),
    "\n"
  )
})
