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

  input_file <- file.path("data", "map_detection", "map_detection_master.xlsx")
  output_dir <- "MAP_detection_analysis"

  if (!file.exists(input_file)) {
    stop(
      "Required analysis workbook not found: ", input_file,
      "\nThe deposited kma_IS900_F57_neu.xlsx contains raw KMA hits only; ",
      "script 06 also needs the matched qPCR/culture and sequencing-QC ",
      "tables in sheets Samples, MAP Evidence (KMA), and Sequencing QC."
    )
  }

  dir.create(output_dir, showWarnings = FALSE)

  to_numeric <- function(x) {
    suppressWarnings(
      as.numeric(
        gsub(",", ".", as.character(x), fixed = TRUE)
      )
    )
  }

  is_positive <- function(x) {
    tolower(trimws(as.character(x))) %in%
      c("yes", "positive", "true", "1")
  }

  samples <- read_excel(
    input_file,
    sheet = "Samples"
  )

  sheet_names <- excel_sheets(input_file)

  kma_sheet <- sheet_names[
    trimws(sheet_names) == "MAP Evidence (KMA)"
  ]

  if (length(kma_sheet) != 1) {
    stop("Could not uniquely identify MAP Evidence (KMA) sheet.")
  }

  kma <- read_excel(
    input_file,
    sheet = kma_sheet
  )

  qc <- read_excel(
    input_file,
    sheet = "Sequencing QC"
  )

  dat <- samples %>%
    transmute(
      id = `Sample ID`,
      country = Country,
      group = `Clinical group`,
      qpcr = as.character(`qPCR result`),
      ct = to_numeric(`qPCR Ct`),
      culture = as.character(`Culture result`)
    ) %>%
    left_join(
      kma %>%
        transmute(
          id = `Sample ID`,

          is900 = is_positive(`IS900 hit reported`),
          f57 = is_positive(`F57 hit reported`),

          is900_depth = to_numeric(`IS900 depth (×)`),
          f57_depth = to_numeric(`F57 depth (×)`),

          is900_breadth = to_numeric(`IS900 breadth (%)`),
          f57_breadth = to_numeric(`F57 breadth (%)`)
        ),
      by = "id"
    ) %>%
    left_join(
      qc %>%
        transmute(
          id = `Sample ID`,
          yield_gbp = to_numeric(
            `Host-removed Yield (Gbp)`
          )
        ),
      by = "id"
    ) %>%
    filter(
      !id %in% c("SW_POS_00", "SW_POS_01")
    ) %>%
    mutate(
      conventional_positive =
        is_positive(qpcr) | is_positive(culture),

      either = is900 | f57,
      both = is900 & f57,

      is900_depth = ifelse(
        !is900 & is.na(is900_depth),
        0, is900_depth
      ),

      f57_depth = ifelse(
        !f57 & is.na(f57_depth),
        0, f57_depth
      ),

      is900_breadth = ifelse(
        !is900 & is.na(is900_breadth),
        0, is900_breadth
      ),

      f57_breadth = ifelse(
        !f57 & is.na(f57_breadth),
        0, f57_breadth
      )
    )

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
