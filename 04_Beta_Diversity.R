# 04_Beta_Diversity.R: Bray–Curtis diversity, PERMANOVA and figures

# Map analysis methods to the original filenames in the Zenodo archive.
taxonomy_file <- function(method, rank, kind) {
  prefix <- ifelse(method == "MetaPhlAn_default", "MetaPhlAn", method)
  file.path("data", "taxonomy", method,
            paste0(prefix, "_", rank, "_", kind, ".tsv"))
}
# See README.md for inputs, provenance and validation status.


local({
  required_packages <- c('tidyverse', 'vegan', 'patchwork', 'writexl')
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop('Install these packages first: ', paste(missing_packages, collapse = ', '))
  }
  suppressPackageStartupMessages({
    library(tidyverse)
    library(vegan)
    library(patchwork)
    library(writexl)
  })

  set.seed(12345)
  nperm <- 9999
  input_dir <- file.path('data', 'taxonomy')
  output_dir <- 'DVM_Beta_Diversity_FINAL'
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  methods <- c('MetaPhlAn_GTDB', 'MetaPhlAn_default', 'ganon2',
               'Kraken2', 'Kraken2_conf10')
  ranks <- c('Genus', 'Species')
  map_levels <- c('MAP-negative', 'Subclinical', 'Clinical')
  cohort_levels <- c('SW MAP-negative', 'SW subclinical', 'SW clinical',
                     'ZA MAP-negative', 'ZA subclinical')

  cohort_colours <- c('SW MAP-negative' = '#69B4C4',
                      'SW subclinical'  = '#AB99C8',
                      'SW clinical'     = '#F0AA9A',
                      'ZA MAP-negative' = '#8FD0DF',
                      'ZA subclinical'  = '#D1C2E3')
  map_colours <- c('MAP-negative' = '#69B4C4',
                   'Subclinical' = '#AB99C8',
                   'Clinical' = '#F0AA9A')
  country_colours <- c('Switzerland' = '#69B4C4',
                       'South Africa' = '#F0AA9A')
  farm_colours <- c('Farm 1' = '#E49CAD',
                    'Farm 2' = '#D8A772',
                    'Farm 3' = '#C2B84E',
                    'Farm 4' = '#87C36E',
                    'Farm 5' = '#55C4A3',
                    'Farm 6' = '#58BAC8',
                    'Farm 7' = '#51A8D4',
                    'Farm 8' = '#B39ADB',
                    'Unreplicated/unknown farm' = '#B8BBC2')
  method_labels <- c('MetaPhlAn_GTDB' = 'MetaPhlAn GTDB',
                     'MetaPhlAn_default' = 'MetaPhlAn default',
                     'ganon2' = 'ganon2',
                     'Kraken2' = 'Kraken2',
                     'Kraken2_conf10' = 'Kraken2 conf. 0.10')

  for (m in methods) for (rk in ranks) {
    f <- taxonomy_file(m, rk, 'normalized')
    if (!file.exists(f)) stop('Missing normalized input table: ', f)
  }

  metadata_file <- file.path('data', 'metadata.txt')
  if (!file.exists(metadata_file)) stop('Metadata file missing: ', metadata_file)
  message('Metadata source: ', metadata_file)

  metadata_raw <- readr::read_tsv(
    metadata_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )

  pick_col <- function(x, choices) {
    hit <- intersect(choices, names(x))
    if (!length(hit)) stop('Metadata must contain one of: ', paste(choices, collapse = ', '))
    str_trim(as.character(x[[hit[1]]]))
  }

  animal_id_column <- intersect(c('Animal_ID', 'animal_id', 'Animal ID'), names(metadata_raw))
  meta <- tibble(
    SAMPLE_ID = pick_col(metadata_raw, c('SAMPLE_ID', 'Sample_ID_Neu')),
    Animal_ID = if (length(animal_id_column))
      str_trim(as.character(metadata_raw[[animal_id_column[1]]]))
    else rep(NA_character_, nrow(metadata_raw)),
    Farm_number = pick_col(metadata_raw, c('Farm number')),
    GROUP = pick_col(metadata_raw, c('GROUP', 'group'))
  ) %>%
    filter(str_detect(SAMPLE_ID, '^(SW|ZA)_(NEG|SUB|POS)_[0-9]+$')) %>%
    mutate(
      Country = if_else(str_starts(SAMPLE_ID, 'SW_'), 'Switzerland', 'South Africa'),
      MAP_group = case_when(
        str_detect(SAMPLE_ID, '_NEG_') ~ 'MAP-negative',
        str_detect(SAMPLE_ID, '_SUB_') ~ 'Subclinical',
        str_detect(SAMPLE_ID, '_POS_') ~ 'Clinical'
      ),
      Cohort_group = case_when(
        str_starts(SAMPLE_ID, 'SW_NEG_') ~ 'SW MAP-negative',
        str_starts(SAMPLE_ID, 'SW_SUB_') ~ 'SW subclinical',
        str_starts(SAMPLE_ID, 'SW_POS_') ~ 'SW clinical',
        str_starts(SAMPLE_ID, 'ZA_NEG_') ~ 'ZA MAP-negative',
        str_starts(SAMPLE_ID, 'ZA_SUB_') ~ 'ZA subclinical'
      )
    )

  if (anyDuplicated(meta$SAMPLE_ID)) stop('Duplicate IDs in metadata.')

  # Farm assignments use the deposited Farm number column.
  # Check known sample-to-animal links when Animal_ID is provided.
  if (length(animal_id_column)) {
    expected_animal_to_new <- c(E10 = 'SW_NEG_16', E2 = 'SW_NEG_17',
                                G6 = 'SW_POS_06', S1 = 'SW_POS_08')
    for (animal in names(expected_animal_to_new)) {
      actual <- meta$SAMPLE_ID[meta$Animal_ID == animal]
      if (length(actual) != 1 || actual != expected_animal_to_new[[animal]]) {
        stop('Incorrect animal-to-sample mapping for ', animal)
      }
    }
  } else {
    warning('Animal_ID is absent from metadata; animal-to-sample cross-check skipped. ',
            'Re-upload the corrected metadata before publication.')
  }
  expected_farms <- c('SW_POS_06' = '6', 'SW_POS_08' = '7',
                      'SW_POS_07' = '34')
  for (id in names(expected_farms)) {
    observed_farm <- meta$Farm_number[meta$SAMPLE_ID == id]
    if (length(observed_farm) != 1 || observed_farm != expected_farms[[id]]) {
      stop('Unexpected farm assignment for ', id)
    }
  }

  qc_candidates <- c(
    file.path('DVM_Alpha_Diversity_FINAL', 'Sample_inclusion_QC.csv'),
    'Sample_inclusion_QC.csv'
  )
  qc_file <- qc_candidates[file.exists(qc_candidates)][1]
  if (!is.na(qc_file)) {
    inclusion <- readr::read_csv(qc_file, show_col_types = FALSE) %>%
      mutate(SAMPLE_ID = as.character(SAMPLE_ID),
             Included = as.logical(Included)) %>%
      filter(Included) %>% pull(SAMPLE_ID)
    if (!length(inclusion)) stop('Alpha sample-inclusion QC contains no included samples.')
    meta <- meta %>% filter(SAMPLE_ID %in% inclusion)
    if (!setequal(meta$SAMPLE_ID, inclusion)) {
      stop('Metadata and alpha inclusion file disagree on sample IDs. Check ', qc_file)
    }
    message('Reusing sample inclusion from: ', qc_file)
  } else {
    meta <- meta %>%
      filter(!GROUP %in% c('EX', 'Excluded', 'MAP detection'),
             !SAMPLE_ID %in% c('SW_NEG_01', 'ZA_NEG_01', 'ZA_NEG_08'))
    message('No alpha sample-inclusion QC found: applying the alpha script exclusions.')
  }

  meta <- meta %>% mutate(
    Farm_code = if_else(Country == 'Switzerland', Farm_number, NA_character_),
    Farm = if_else(Farm_code %in% as.character(1:8),
                   paste('Farm', Farm_code), NA_character_),
    MAP_group = factor(MAP_group, levels = map_levels),
    Country = factor(Country, levels = c('Switzerland', 'South Africa')),
    Cohort_group = factor(Cohort_group, levels = cohort_levels)
  ) %>% arrange(SAMPLE_ID)

  expected <- c('SW MAP-negative' = 16L, 'SW subclinical' = 13L,
                'SW clinical' = 7L, 'ZA MAP-negative' = 6L,
                'ZA subclinical' = 8L)
  observed <- table(meta$Cohort_group)
  if (!identical(unname(as.integer(observed[names(expected)])), unname(expected))) {
    stop('Sample counts differ from final alpha analysis. Found: ',
         paste(names(observed), observed, collapse = '; '))
  }
  message('Samples match alpha: SW 36 (16/13/7), ZA 14 (6/8).')

  replicated_farms <- meta %>%
    filter(Country == 'Switzerland', !is.na(Farm)) %>%
    count(Farm, name = 'n') %>%
    filter(n >= 2) %>% pull(Farm)
  if (length(replicated_farms) != 8) {
    warning('Expected eight replicated Swiss farms, got ', length(replicated_farms),
            ': ', paste(replicated_farms, collapse = ', '))
  }
  meta <- meta %>% mutate(
    Farm_test = if_else(as.character(Country) == 'Switzerland' &
                          Farm %in% replicated_farms, Farm, NA_character_),
    Farm_plot = if_else(is.na(Farm_test), 'Unreplicated/unknown farm', Farm_test)
  )

  farm_qc <- meta %>% filter(Country == 'Switzerland') %>%
    count(Farm, MAP_group, name = 'n', .drop = FALSE) %>%
    filter(n > 0)
  readr::write_csv(meta, file.path(output_dir, 'Sample_inclusion_and_real_farms.csv'))
  readr::write_csv(farm_qc, file.path(output_dir, 'Swiss_actual_farm_QC.csv'))

  read_profile <- function(method, rank, sample_ids) {
    f <- taxonomy_file(method, rank, 'normalized')
    raw <- readr::read_tsv(f, show_col_types = FALSE, name_repair = 'minimal')
    if (!'clade_name' %in% names(raw)) stop('No clade_name in ', f)
    missing <- setdiff(sample_ids, names(raw))
    if (length(missing)) stop(method, ' ', rank, ' missing samples: ',
                              paste(missing, collapse = ', '))

    terminal <- str_remove(basename(gsub('[|;]', '/', raw$clade_name)), '^[a-z]__')
    keep_taxa <- !is.na(raw$clade_name) & nzchar(raw$clade_name) &
      !is.na(terminal) & nzchar(trimws(terminal)) &
      !tolower(trimws(terminal)) %in% c('unclassified', 'unassigned', 'unknown')
    raw <- raw[keep_taxa, , drop = FALSE]
    taxon <- as.character(raw$clade_name)

    values <- lapply(raw[sample_ids], function(x) suppressWarnings(as.numeric(x)))
    na_counts <- vapply(values, function(x) sum(is.na(x)), integer(1))
    for (id in names(values)) {
      x <- values[[id]]
      if (anyNA(x) && !all(is.na(x))) {
        stop('Partial NA profile in ', method, ' / ', rank, ' / ', id,
             '; inspect input rather than replacing unknown abundance with zero.')
      }
      if (all(is.na(x))) values[[id]] <- rep(0, length(x))
    }
    mat <- t(as.matrix(as.data.frame(values, check.names = FALSE)))
    rownames(mat) <- sample_ids
    colnames(mat) <- make.unique(taxon, sep = '__duplicate_')
    storage.mode(mat) <- 'double'
    if (any(!is.finite(mat)) || any(mat < 0)) stop('Invalid abundance in ', f)
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    sums <- rowSums(mat)
    valid <- sums > 0
    if (any(valid)) mat[valid, ] <- sweep(mat[valid, , drop = FALSE], 1, sums[valid], '/')
    qc <- tibble(Method = method, Rank = rank, SAMPLE_ID = sample_ids,
                 Terminal_taxa_removed = sum(!keep_taxa),
                 Input_NA_cells = as.integer(na_counts[sample_ids]),
                 Valid_profile = valid,
                 Taxa_detected = rowSums(mat > 0),
                 Sum_after_normalization = rowSums(mat))
    list(matrix = mat, qc = qc)
  }

  all_ids <- meta$SAMPLE_ID
  profile_list <- list()
  qc_list <- list()
  for (method in methods) for (rank in ranks) {
    key <- paste(method, rank, sep = '__')
    res <- read_profile(method, rank, all_ids)
    profile_list[[key]] <- res$matrix
    qc_list[[key]] <- res$qc
  }
  profile_qc <- bind_rows(qc_list)
  readr::write_csv(profile_qc, file.path(output_dir, 'Input_profiles_QC.csv'))

  for (rank in ranks) {
    q <- profile_qc %>% filter(Method == 'MetaPhlAn_GTDB', Rank == rank)
    if (!all(q$Valid_profile)) stop('Invalid primary GTDB profile: ', rank)
  }

  matched_swiss <- profile_qc %>%
    filter(SAMPLE_ID %in% meta$SAMPLE_ID[meta$Country == 'Switzerland']) %>%
    group_by(SAMPLE_ID) %>%
    summarise(N_profiles = n(), Valid_all = all(Valid_profile), .groups = 'drop') %>%
    filter(N_profiles == length(methods) * length(ranks), Valid_all) %>%
    pull(SAMPLE_ID)
  matched_swiss <- meta$SAMPLE_ID[meta$SAMPLE_ID %in% matched_swiss]
  if (length(matched_swiss) != 35) {
    stop('Expected 35 matched Swiss samples, found ', length(matched_swiss),
         '. Inspect Input_profiles_QC.csv.')
  }
  message('Method comparison: ', length(matched_swiss), ' identical Swiss samples.')

  get_matrix <- function(method, rank, ids) {
    mat <- profile_list[[paste(method, rank, sep = '__')]][ids, , drop = FALSE]
    if (any(rowSums(mat) <= 0)) {
      stop('Zero-profile sample in tested subset: ', method, ' ', rank)
    }
    mat[, colSums(mat) > 0, drop = FALSE]
  }

  subset_matmeta <- function(mat, dat, ids) {
    ids <- as.character(ids)
    list(mat = mat[ids, , drop = FALSE],
         dat = dat[match(ids, dat$SAMPLE_ID), , drop = FALSE])
  }

  run_beta_test <- function(mat, dat, grouping, analysis, method, rank) {
    stopifnot(identical(rownames(mat), as.character(dat$SAMPLE_ID)))
    g <- droplevels(factor(dat[[grouping]]))
    group_counts <- table(g)
    blank <- tibble(Analysis = analysis, Method = method, Rank = rank,
                    Group_variable = grouping, N = nrow(dat),
                    N_groups = nlevels(g), Group_counts = paste(names(group_counts),
                       group_counts, sep = ':', collapse = '; '),
                    PERMANOVA_df = NA_real_, PERMANOVA_F = NA_real_,
                    PERMANOVA_R2 = NA_real_, PERMANOVA_p = NA_real_,
                    PERMDISP_F = NA_real_, PERMDISP_p = NA_real_,
                    Note = 'Not estimable')
    if (anyNA(g) || nlevels(g) < 2 || any(group_counts < 2) ||
        nrow(dat) <= nlevels(g) + 1) return(blank)

    bray <- vegan::vegdist(mat, method = 'bray')
    set.seed(12345)
    fit <- vegan::adonis2(bray ~ group, data = data.frame(group = g),
                          permutations = nperm)
    bd <- tryCatch({
      dis <- vegan::betadisper(bray, group = g, type = 'median',
                               bias.adjust = TRUE)
      set.seed(12345)
      perm <- vegan::permutest(dis, permutations = nperm)
      c(F = as.numeric(perm$tab[1, 'F']),
        P = as.numeric(perm$tab[1, 'Pr(>F)']))
    }, error = function(e) {
      warning('PERMDISP failed in ', analysis, ' ', method, ' ', rank, ': ',
              conditionMessage(e))
      c(F = NA_real_, P = NA_real_)
    })
    blank %>% mutate(
      PERMANOVA_df = as.numeric(fit$Df[1]),
      PERMANOVA_F = as.numeric(fit$F[1]),
      PERMANOVA_R2 = as.numeric(fit$R2[1]),
      PERMANOVA_p = as.numeric(fit$`Pr(>F)`[1]),
      PERMDISP_F = unname(bd['F']),
      PERMDISP_p = unname(bd['P']),
      Note = if_else(any(group_counts == 2),
                     'Contains a group with n=2: dispersion inference is fragile',
                     'Bray-Curtis; unrestricted permutations')
    )
  }

  run_pairwise <- function(mat, dat, grouping, analysis, method, rank) {
    groups <- names(table(droplevels(factor(dat[[grouping]]))))
    groups <- groups[!is.na(groups) & groups != '']
    if (length(groups) < 3) return(tibble())
    pairs <- combn(groups, 2, simplify = FALSE)
    res <- purrr::map_dfr(pairs, function(pair) {
      ids <- as.character(dat$SAMPLE_ID[as.character(dat[[grouping]]) %in% pair])
      sub <- subset_matmeta(mat, dat, ids)
      test <- run_beta_test(sub$mat, sub$dat, grouping, analysis,
                            method, rank)
      test %>% mutate(Group1 = pair[1], Group2 = pair[2], .after = Group_variable)
    })
    if (nrow(res)) {
      res <- res %>% mutate(
        PERMANOVA_p_BH_within_analysis = p.adjust(PERMANOVA_p, method = 'BH'),
        PERMDISP_p_BH_within_analysis = p.adjust(PERMDISP_p, method = 'BH')
      )
    }
    res
  }

  get_pcoa <- function(mat, dat, dataset, method, rank) {
    if (nrow(mat) < 4) stop('Insufficient PCoA samples: ', dataset)
    d <- vegan::vegdist(mat, method = 'bray')
    pco <- vegan::wcmdscale(d, k = 2, eig = TRUE, add = FALSE)
    pts <- as.data.frame(pco$points)
    if (ncol(pts) < 2) stop('Only one PCoA axis in: ', dataset)
    names(pts)[1:2] <- c('PCoA1', 'PCoA2')
    pts$SAMPLE_ID <- rownames(pco$points)
    eig <- as.numeric(pco$eig)
    positive <- eig[eig > 0]
    var1 <- 100 * eig[1] / sum(positive)
    var2 <- 100 * eig[2] / sum(positive)
    points <- left_join(dat, pts, by = 'SAMPLE_ID')
    if (anyNA(points$PCoA1) || anyNA(points$PCoA2)) {
      stop('PCoA and metadata did not join correctly: ', dataset)
    }
    list(points = points %>% mutate(Dataset = dataset, Method = method, Rank = rank),
         var1 = var1, var2 = var2,
         eigen_qc = tibble(Dataset = dataset, Method = method, Rank = rank,
                           N = nrow(mat), Number_negative_eigenvalues = sum(eig < -1e-9),
                           Sum_positive_eigenvalues = sum(positive),
                           PCoA1_positive_eigen_pct = var1,
                           PCoA2_positive_eigen_pct = var2))
  }

  fmt_p <- function(p) {
    if (length(p) == 0 || is.na(p)) return('NA')
    if (p < 0.001) return('< 0.001')
    sprintf('%.3f', p)
  }

  stat_subtitle <- function(test) {
    if (!nrow(test) || is.na(test$PERMANOVA_p[1])) return('PERMANOVA not estimable')
    paste0('n = ', test$N[1], ' | PERMANOVA R² = ',
           sprintf('%.3f', test$PERMANOVA_R2[1]),
           ', p = ', fmt_p(test$PERMANOVA_p[1]),
           ' | PERMDISP p = ', fmt_p(test$PERMDISP_p[1]))
  }

  plot_pcoa <- function(ord, aesthetic_col, palette, plot_title,
                        subtitle = NULL, restrict_ids = NULL) {
    dd <- ord$points
    if (!is.null(restrict_ids)) dd <- dd %>% filter(SAMPLE_ID %in% restrict_ids)
    dd$Plot_group <- factor(as.character(dd[[aesthetic_col]]), levels = names(palette))
    ggplot(dd, aes(x = PCoA1, y = PCoA2, fill = Plot_group)) +
      geom_hline(yintercept = 0, colour = '#DEE3E6', linewidth = 0.3) +
      geom_vline(xintercept = 0, colour = '#DEE3E6', linewidth = 0.3) +
      geom_point(shape = 21, size = 3.05, stroke = 0.40,
                 colour = '#34434B', alpha = 0.96) +
      scale_fill_manual(values = palette, drop = TRUE, name = NULL) +
      coord_equal() +
      labs(title = plot_title, subtitle = subtitle,
           x = sprintf('PCoA1 (%.1f%%)', ord$var1),
           y = sprintf('PCoA2 (%.1f%%)', ord$var2)) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_text(face = 'bold', size = 12),
        plot.subtitle = element_text(size = 9),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        legend.position = 'bottom',
        legend.text = element_text(size = 9),
        legend.key.height = grid::unit(0.35, 'cm'),
        panel.border = element_rect(colour = '#D6DCE1', fill = NA, linewidth = 0.5),
        plot.margin = margin(7, 8, 7, 8)
      ) +
      guides(fill = guide_legend(override.aes = list(size = 3.1), nrow = 2))
  }

  save_plot <- function(fig, filename, width = 11, height = 7) {
    ggplot2::ggsave(file.path(output_dir, filename), fig, width = width,
                    height = height, units = 'in', dpi = 300, bg = 'white')
    ggplot2::ggsave(file.path(output_dir, sub('\\.png$', '.pdf', filename)),
                    fig, width = width, height = height, units = 'in', bg = 'white')
  }

  primary_tests <- list()
  primary_pairs <- list()
  primary_ordinations <- list()
  eigen_qc <- list()
  group_sizes <- list()

  for (rank in ranks) {
    method <- 'MetaPhlAn_GTDB'
    mat <- get_matrix(method, rank, meta$SAMPLE_ID)
    swiss_ids <- meta$SAMPLE_ID[meta$Country == 'Switzerland']
    za_ids <- meta$SAMPLE_ID[meta$Country == 'South Africa']
    sw <- subset_matmeta(mat, meta, swiss_ids)
    za <- subset_matmeta(mat, meta, za_ids)

    datasets <- list()
    datasets[['Swiss_MAP']] <- c(sw, list(grouping = 'MAP_group'))
    datasets[['SouthAfrica_MAP']] <- c(za, list(grouping = 'MAP_group'))
    datasets[['Country_SW_vs_ZA']] <- list(mat = mat, dat = meta,
                                          grouping = 'Country')
    datasets[['Combined_5_groups']] <- list(mat = mat, dat = meta,
                                            grouping = 'Cohort_group')

    farm_ids <- swiss_ids[swiss_ids %in% meta$SAMPLE_ID[!is.na(meta$Farm_test)]]
    farm <- subset_matmeta(mat, meta, farm_ids)
    datasets[['Swiss_farms']] <- c(farm, list(grouping = 'Farm_test'))

    for (status in c('MAP-negative', 'Subclinical')) {
      dd <- sw$dat %>% filter(MAP_group == status, !is.na(Farm_test)) %>%
        add_count(Farm_test, name = 'n_within') %>% filter(n_within >= 2)
      label <- if (status == 'MAP-negative') 'Swiss_farms_within_NEG' else
        'Swiss_farms_within_SUB'
      ss <- subset_matmeta(mat, meta, dd$SAMPLE_ID)
      datasets[[label]] <- c(ss, list(grouping = 'Farm_test'))
    }

    for (status in c('MAP-negative', 'Subclinical')) {
      ids <- meta$SAMPLE_ID[as.character(meta$MAP_group) == status]
      ss <- subset_matmeta(mat, meta, ids)
      label <- if (status == 'MAP-negative') 'SW_vs_ZA_within_NEG' else
        'SW_vs_ZA_within_SUB'
      datasets[[label]] <- c(ss, list(grouping = 'Country'))
    }

    for (analysis in names(datasets)) {
      ds <- datasets[[analysis]]
      test <- run_beta_test(ds$mat, ds$dat, ds$grouping,
                            analysis, method, rank)
      primary_tests[[paste(rank, analysis, sep = '__')]] <- test
      pair <- run_pairwise(ds$mat, ds$dat, ds$grouping,
                           analysis, method, rank)
      if (nrow(pair)) primary_pairs[[paste(rank, analysis, sep = '__')]] <- pair
      group_sizes[[paste(rank, analysis, sep = '__')]] <-
        ds$dat %>%
        group_by(Group = as.character(.data[[ds$grouping]])) %>%
        summarise(N = n(), .groups = 'drop') %>%
        mutate(Analysis = analysis, Rank = rank, .before = 1)
    }

    swiss_ord <- get_pcoa(sw$mat, sw$dat, 'Switzerland', method, rank)
    za_ord <- get_pcoa(za$mat, za$dat, 'South Africa', method, rank)
    combined_ord <- get_pcoa(mat, meta, 'Swiss + South Africa', method, rank)
    primary_ordinations[[paste(rank, 'Swiss', sep='__')]] <- swiss_ord$points
    primary_ordinations[[paste(rank, 'ZA', sep='__')]] <- za_ord$points
    primary_ordinations[[paste(rank, 'Combined', sep='__')]] <- combined_ord$points
    eigen_qc[[paste(rank, 'Swiss', sep='__')]] <- swiss_ord$eigen_qc
    eigen_qc[[paste(rank, 'ZA', sep='__')]] <- za_ord$eigen_qc
    eigen_qc[[paste(rank, 'Combined', sep='__')]] <- combined_ord$eigen_qc

    lookup_test <- function(name) primary_tests[[paste(rank, name, sep = '__')]]
    p_sw <- plot_pcoa(swiss_ord, 'MAP_group', map_colours,
                      paste('A | Swiss MAP groups |', rank),
                      stat_subtitle(lookup_test('Swiss_MAP')))
    p_farm <- plot_pcoa(swiss_ord, 'Farm_plot', farm_colours,
                        paste('B | Swiss farms |', rank),
                        stat_subtitle(lookup_test('Swiss_farms'))) +
      labs(caption = 'Grey points: independent/unreplicated farms; excluded from farm test') +
      theme(plot.caption = element_text(size = 8, colour = '#59646A'))
    p_za <- plot_pcoa(za_ord, 'MAP_group', map_colours[1:2],
                      paste('South African MAP groups |', rank),
                      stat_subtitle(lookup_test('SouthAfrica_MAP')))
    p_country <- plot_pcoa(combined_ord, 'Country', country_colours,
                           paste('C | Switzerland vs South Africa |', rank),
                           stat_subtitle(lookup_test('Country_SW_vs_ZA')))
    p_combined <- plot_pcoa(combined_ord, 'Cohort_group', cohort_colours,
                            paste('Five MAP-status/cohort groups |', rank),
                            stat_subtitle(lookup_test('Combined_5_groups')))

    save_plot(p_sw, paste0('PCoA_', rank, '_Swiss_MAP.png'), 8, 6)
    save_plot(p_farm, paste0('PCoA_', rank, '_Swiss_Farms.png'), 9, 7)
    save_plot(p_za, paste0('PCoA_', rank, '_ZA_MAP.png'), 8, 6)
    save_plot(p_country, paste0('PCoA_', rank, '_Country.png'), 8, 6)
    save_plot(p_combined, paste0('PCoA_', rank, '_Combined_five_groups.png'), 9, 6)

    p_abc <- (p_sw + p_farm + p_country +
                plot_layout(ncol = 3, guides = 'collect')) &
      theme(legend.position = 'bottom')
    save_plot(p_abc, paste0('PCoA_', rank, '_ABC_combined.png'), 18, 7)
    p_geography <- (p_za + p_country + p_combined +
                      plot_layout(ncol = 3, guides = 'collect')) &
      theme(legend.position = 'bottom')
    save_plot(p_geography, paste0('PCoA_', rank, '_ZA_Country_Combined.png'), 18, 7)
  }

  primary_overall <- bind_rows(primary_tests) %>%
    group_by(Analysis) %>%
    mutate(PERMANOVA_p_BH_two_ranks = p.adjust(PERMANOVA_p, 'BH'),
           PERMDISP_p_BH_two_ranks = p.adjust(PERMDISP_p, 'BH')) %>%
    ungroup() %>%
    mutate(PERMANOVA_p_BH_all_primary = p.adjust(PERMANOVA_p, 'BH'),
           PERMDISP_p_BH_all_primary = p.adjust(PERMDISP_p, 'BH'))

  primary_pairwise <- bind_rows(primary_pairs)
  if (nrow(primary_pairwise)) {
    primary_pairwise <- primary_pairwise %>%
      mutate(PERMANOVA_p_BH_all_primary_pairwise = p.adjust(PERMANOVA_p, 'BH'),
             PERMDISP_p_BH_all_primary_pairwise = p.adjust(PERMDISP_p, 'BH'))
  }

  method_overall_list <- list()
  method_pair_list <- list()
  method_ord_list <- list()
  method_distance <- list()
  for (rank in ranks) {
    panels <- list()
    for (method in methods) {
      mat <- get_matrix(method, rank, matched_swiss)
      dat <- meta[match(matched_swiss, meta$SAMPLE_ID), , drop = FALSE]
      test <- run_beta_test(mat, dat, 'MAP_group', 'Swiss_MAP_matched35',
                            method, rank)
      pair <- run_pairwise(mat, dat, 'MAP_group', 'Swiss_MAP_matched35',
                           method, rank)
      method_overall_list[[paste(rank, method, sep = '__')]] <- test
      if (nrow(pair)) method_pair_list[[paste(rank, method, sep = '__')]] <- pair
      method_distance[[paste(rank, method, sep = '__')]] <-
        vegan::vegdist(mat, method = 'bray')
      ord <- get_pcoa(mat, dat, 'Swiss matched (35)', method, rank)
      method_ord_list[[paste(rank, method, sep = '__')]] <- ord$points
      eigen_qc[[paste(rank, method, sep = '__')]] <- ord$eigen_qc
      panels[[method]] <- plot_pcoa(ord, 'MAP_group', map_colours,
                                   method_labels[[method]], stat_subtitle(test))
      save_plot(panels[[method]], paste0('PCoA_', rank, '_Method_', method, '.png'),
                8, 6)
    }
    comparison_fig <- wrap_plots(panels, ncol = 3, guides = 'collect') &
      theme(legend.position = 'bottom')
    save_plot(comparison_fig, paste0('PCoA_', rank, '_Five_methods.png'),
              17, 11)
  }

  method_overall <- bind_rows(method_overall_list) %>%
    group_by(Rank) %>%
    mutate(PERMANOVA_p_BH_five_methods = p.adjust(PERMANOVA_p, 'BH'),
           PERMDISP_p_BH_five_methods = p.adjust(PERMDISP_p, 'BH')) %>%
    ungroup() %>%
    mutate(PERMANOVA_p_BH_all10 = p.adjust(PERMANOVA_p, 'BH'),
           PERMDISP_p_BH_all10 = p.adjust(PERMDISP_p, 'BH'))
  method_pairwise <- bind_rows(method_pair_list)
  if (nrow(method_pairwise)) {
    method_pairwise <- method_pairwise %>%
      mutate(PERMANOVA_p_BH_all_method_pairs = p.adjust(PERMANOVA_p, 'BH'),
             PERMDISP_p_BH_all_method_pairs = p.adjust(PERMDISP_p, 'BH'))
  }

  method_agreement <- list()
  for (rank in ranks) {
    pairs <- combn(methods, 2, simplify = FALSE)
    for (pair in pairs) {
      method_a <- pair[1]; method_b <- pair[2]
      ka <- paste(rank, method_a, sep = '__')
      kb <- paste(rank, method_b, sep = '__')
      d1 <- method_distance[[ka]]
      d2 <- method_distance[[kb]]
      if (!identical(attr(d1, 'Labels'), attr(d2, 'Labels'))) {
        stop('Sample ordering differs in distance matrices: ', ka, ' and ', kb)
      }
      set.seed(12345)
      mantel_test <- vegan::mantel(d1, d2, method = 'spearman',
                                  permutations = nperm)
      o1 <- method_ord_list[[ka]] %>%
        arrange(match(SAMPLE_ID, matched_swiss)) %>% select(PCoA1, PCoA2) %>% as.matrix()
      o2 <- method_ord_list[[kb]] %>%
        arrange(match(SAMPLE_ID, matched_swiss)) %>% select(PCoA1, PCoA2) %>% as.matrix()
      set.seed(12345)
      protest_test <- vegan::protest(o1, o2, permutations = nperm,
                                    symmetric = TRUE)
      method_agreement[[paste(rank, method_a, method_b, sep = '__')]] <-
        tibble(Rank = rank, Method1 = method_a, Method2 = method_b,
               N = length(matched_swiss),
               Mantel_Spearman_r = unname(mantel_test$statistic),
               Mantel_p = mantel_test$signif,
               Procrustes_correlation = unname(protest_test$t0),
               PROTEST_p = protest_test$signif)
    }
  }
  method_agreement <- bind_rows(method_agreement) %>%
    mutate(Mantel_p_BH_all20 = p.adjust(Mantel_p, 'BH'),
           PROTEST_p_BH_all20 = p.adjust(PROTEST_p, 'BH'))

  pcoa_coordinates <- bind_rows(c(primary_ordinations, method_ord_list))
  eigen_diagnostics <- bind_rows(eigen_qc)
  group_sizes <- bind_rows(group_sizes)

  result_tables <- list(
    'Primary_overall' = primary_overall,
    'Primary_pairwise' = primary_pairwise,
    'Primary_group_sizes' = group_sizes,
    'Method_MAP_overall' = method_overall,
    'Method_MAP_pairwise' = method_pairwise,
    'Method_concordance' = method_agreement,
    'PCoA_coordinates' = pcoa_coordinates,
    'PCoA_eigen_QC' = eigen_diagnostics,
    'Sample_QC' = profile_qc,
    'Swiss_real_farms' = farm_qc
  )

  for (tab in names(result_tables)) {
    readr::write_csv(result_tables[[tab]],
                     file.path(output_dir, paste0(tab, '.csv')))
  }
  writexl::write_xlsx(result_tables,
                      path = file.path(output_dir, 'DVM_Beta_Diversity_Master.xlsx'))

  readme <- c(
    'DVM BETA DIVERSITY - NOTES',
    paste0('Run date: ', Sys.Date()),
    'Distance: Bray-Curtis from rank-normalized classified-taxon abundances.',
    'Ordination: wcmdscale on ORIGINAL Bray-Curtis (no Lingoes correction).',
    'PCoA axis percentages use the sum of positive eigenvalues as denominator.',
    'PERMANOVA: vegan::adonis2, 9999 unrestricted permutations, seed 12345.',
    'PERMDISP: betadisper(type=median, bias.adjust=TRUE), permutest 9999.',
    'Inspect PERMDISP alongside PERMANOVA: heterogeneous dispersions may affect inference.',
    'Primary GTDB: Swiss 36, South African 14, combined 50.',
    'Classifier analyses: 35 Swiss animals valid in all 10 rank x method profiles.',
    'Farm: use the deposited metadata Farm number; farm >=2 animals overall.',
    'Clinical Dia* samples and other independent farms NEVER pooled as a real farm.',
    'Farm tests may be confounded by MAP status. Interpret as exploratory.',
    'Within-group farm tests require >=2 animals per farm in the corresponding MAP group.',
    'Pairwise PERMANOVA is exploratory; p_BH_within_analysis adjusts pairs within each test.',
    'Primary overall BH_all_primary adjusts all primary analyses and both ranks together.',
    'Method overall BH_all10 adjusts the five profilers x two ranks.',
    'Country/cohort comparisons may be affected by extraction and other sampling differences.',
    'A/B Swiss PCoAs reuse the exact same ordination coordinates.',
    'C Country panel uses an independent 50-sample ordination; compare axes within panel.',
    'Method PCoAs are independently computed and need not share orientation.',
    'Procrustes/Mantel measure methodological concordance, not equivalence of taxon names.',
    'NA zero-profile Kraken2 conf10 sample SW_NEG_09 is excluded from matched comparison.',
    'Swiss farms with n=2 can have unstable pairwise PERMDISP results.',
    'PNG at 300 dpi and vector PDF figures are both generated.'
  )
  writeLines(readme, file.path(output_dir, 'README_Analysis_Notes.txt'))

  cat('\nBETA-DIVERSITY ANALYSIS COMPLETE\n')
  cat('Output folder: ', normalizePath(output_dir), '\n', sep = '')
  cat('GTDB overall tests: ', nrow(primary_overall), '\n', sep = '')
  cat('GTDB pairwise tests: ', nrow(primary_pairwise), '\n', sep = '')
  cat('Method overall tests: ', nrow(method_overall), '\n', sep = '')
  cat('Method pairwise tests: ', nrow(method_pairwise), '\n', sep = '')
  cat('Method concordance tests: ', nrow(method_agreement), '\n\n', sep = '')
  print(primary_overall %>% select(Analysis, Rank, N, PERMANOVA_R2,
                                   PERMANOVA_p, PERMDISP_p,
                                   PERMANOVA_p_BH_all_primary), n = Inf)
})

# Manuscript beta-diversity figures
local({
  packages <- c(
    "readxl",
    "dplyr",
    "ggplot2",
    "patchwork"
  )

  missing_packages <- packages[
    !vapply(
      packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_packages) > 0) stop("Required packages missing: ", paste(missing_packages, collapse = ", "))

  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)

  possible_workbooks <- c(
    file.path(
      "DVM_Beta_Diversity_FINAL",
      "DVM_Beta_Diversity_Master.xlsx"
    ),
    "DVM_Beta_Diversity_Master.xlsx"
  )

  existing_workbooks <- possible_workbooks[
    file.exists(possible_workbooks)
  ]

  if (length(existing_workbooks) == 0) {
    stop(
      "DVM_Beta_Diversity_Master.xlsx not found. ",
      "Check getwd() and the input directory."
    )
  }

  workbook <- existing_workbooks[1]

  output_dir <- file.path(
    "DVM_Beta_Diversity_FINAL",
    "Publication_Plots"
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  x_limits <- c(-0.75, 0.75)
  y_limits <- c(-0.65, 0.65)

  tick_step <- 0.2

  x_breaks <- seq(-0.6, 0.6, by = tick_step)
  y_breaks <- seq(-0.6, 0.6, by = tick_step)

  map_colors <- c(
    "MAP-negative" = "#49AFC7",
    "Subclinical"  = "#9479BA",
    "Clinical"     = "#EB907C"
  )

  country_colors <- c(
    "Switzerland"  = "#238AA8",
    "South Africa" = "#DD775C"
  )

  cohort_colors <- c(
    "SW MAP-negative" = "#49AFC7",
    "SW subclinical"  = "#9479BA",
    "SW clinical"     = "#EB907C",
    "ZA MAP-negative" = "#087A86",
    "ZA subclinical"  = "#63418B"
  )

  farm_colors <- c(
    "Farm 1" = "#D55E00",
    "Farm 2" = "#E69F00",
    "Farm 3" = "#9A9C31",
    "Farm 4" = "#009E73",
    "Farm 5" = "#39AAB3",
    "Farm 6" = "#0072B2",
    "Farm 7" = "#7257AF",
    "Farm 8" = "#CC79A7",
    "Unreplicated/unknown farm" = "#B9BFC5"
  )

  country_shapes <- c(
    "Switzerland"  = 21,   # Circle
    "South Africa" = 24    # Triangle
  )

  method_labels <- c(
    "MetaPhlAn_GTDB"    = "MetaPhlAn GTDB",
    "MetaPhlAn_default" = "MetaPhlAn default",
    "ganon2"            = "ganon2",
    "Kraken2"           = "Kraken2",
    "Kraken2_conf10"    = "Kraken2 conf. 0.10"
  )

  methods <- c(
    "MetaPhlAn_GTDB",
    "MetaPhlAn_default",
    "ganon2",
    "Kraken2",
    "Kraken2_conf10"
  )

  ranks <- c("Genus", "Species")

  coords <- read_excel(
    workbook,
    sheet = "PCoA_coordinates"
  )

  eigen <- read_excel(
    workbook,
    sheet = "PCoA_eigen_QC"
  )

  primary <- read_excel(
    workbook,
    sheet = "Primary_overall"
  )

  method_tests <- read_excel(
    workbook,
    sheet = "Method_MAP_overall"
  )

  required_columns <- c(
    "SAMPLE_ID",
    "PCoA1",
    "PCoA2",
    "Dataset",
    "Method",
    "Rank",
    "MAP_group",
    "Country",
    "Cohort_group",
    "Farm_plot"
  )

  missing_columns <- setdiff(
    required_columns,
    names(coords)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing coordinate columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  stopifnot(
    all(is.finite(coords$PCoA1)),
    all(is.finite(coords$PCoA2))
  )

  if (
    any(coords$PCoA1 < x_limits[1]) ||
    any(coords$PCoA1 > x_limits[2]) ||
    any(coords$PCoA2 < y_limits[1]) ||
    any(coords$PCoA2 > y_limits[2])
  ) {
    stop(
      "Some points are outside the chosen axis limits. ",
      "Increase x_limits or y_limits."
    )
  }

  get_coordinates <- function(
      dataset,
      method,
      rank,
      expected_n
  ) {

    dat <- coords %>%
      filter(
        Dataset == dataset,
        Method == method,
        Rank == rank
      )

    if (nrow(dat) != expected_n) {
      stop(
        "Unexpected sample count for ",
        dataset, " / ", method, " / ", rank,
        ". Expected ", expected_n,
        ", found ", nrow(dat)
      )
    }

    dat
  }

  get_eigen <- function(dataset, method, rank) {

    dat <- eigen %>%
      filter(
        Dataset == dataset,
        Method == method,
        Rank == rank
      )

    if (nrow(dat) != 1) {
      stop(
        "Missing or duplicated eigen summary: ",
        dataset, " / ", method, " / ", rank
      )
    }

    dat
  }

  format_p <- function(p) {

    if (length(p) != 1 || is.na(p)) {
      return("NA")
    }

    if (p < 0.001) {
      return("< 0.001")
    }

    sprintf("= %.3f", p)
  }

  get_statistics <- function(
      table,
      rank,
      analysis,
      method = NULL
  ) {

    dat <- table %>%
      filter(Rank == rank)

    if ("Analysis" %in% names(dat)) {
      dat <- dat %>%
        filter(Analysis == analysis)
    }

    if (!is.null(method)) {
      dat <- dat %>%
        filter(Method == method)
    }

    if (nrow(dat) != 1) {
      stop(
        "Missing statistical result: ",
        analysis, " / ", rank
      )
    }

    sprintf(
      "n = %d | PERMANOVA R² = %.3f, p %s | PERMDISP p %s",
      dat$N,
      dat$PERMANOVA_R2,
      format_p(dat$PERMANOVA_p),
      format_p(dat$PERMDISP_p)
    )
  }

  make_pcoa <- function(
      dataset,
      method,
      rank,
      n,
      group_column,
      colors,
      title,
      subtitle,
      legend_columns = 3,
      show_country_shapes = FALSE,
      caption = NULL
  ) {

    dat <- get_coordinates(
      dataset,
      method,
      rank,
      n
    )

    eig <- get_eigen(
      dataset,
      method,
      rank
    )

    dat$Plot_group <- factor(
      as.character(dat[[group_column]]),
      levels = names(colors)
    )

    if (anyNA(dat$Plot_group)) {
      stop("Unmapped group in plot: ", title)
    }

    p <- ggplot(
      dat,
      aes(
        x = PCoA1,
        y = PCoA2
      )
    ) +
      geom_hline(
        yintercept = 0,
        color = "#DEE4E9",
        linewidth = 0.25
      ) +
      geom_vline(
        xintercept = 0,
        color = "#DEE4E9",
        linewidth = 0.25
      )

    if (show_country_shapes) {

      dat$Country_plot <- factor(
        as.character(dat$Country),
        levels = names(country_shapes)
      )

      if (anyNA(dat$Country_plot)) {
        stop(
          "Unmapped country in plot: ",
          title
        )
      }

      p <- p +
        geom_point(
          data = dat,
          mapping = aes(
            fill = Plot_group,
            shape = Country_plot
          ),
          size = 2.8,
          color = "#33414A",
          stroke = 0.35,
          alpha = 0.96
        ) +
        scale_shape_manual(
          name = "Country",
          values = country_shapes,
          drop = FALSE
        )

    } else {

      p <- p +
        geom_point(
          data = dat,
          mapping = aes(
            fill = Plot_group
          ),
          shape = 21,
          size = 2.8,
          color = "#33414A",
          stroke = 0.35,
          alpha = 0.96
        )

    }

    fill_title <- if (show_country_shapes) {
      "MAP status / cohort"
    } else {
      NULL
    }

    p <- p +
      scale_fill_manual(
        name = fill_title,
        values = colors,
        drop = TRUE
      )

    if (show_country_shapes) {

      p <- p +
        guides(
          fill = guide_legend(
            title = "MAP status / cohort",
            order = 1,
            ncol = legend_columns,
            byrow = TRUE,
            override.aes = list(
              shape = 21,
              size = 3
            )
          ),
          shape = guide_legend(
            title = "Country",
            order = 2,
            nrow = 1,
            override.aes = list(
              fill = "#B0B0B0",
              size = 3
            )
          )
        )

    } else {

      p <- p +
        guides(
          fill = guide_legend(
            title = NULL,
            ncol = legend_columns,
            byrow = TRUE,
            override.aes = list(
              shape = 21,
              size = 3
            )
          )
        )

    }

    p <- p +
      scale_x_continuous(
        breaks = x_breaks,
        labels = function(x) sprintf("%.1f", x)
      ) +
      scale_y_continuous(
        breaks = y_breaks,
        labels = function(x) sprintf("%.1f", x)
      ) +
      coord_fixed(
        ratio = 1,
        xlim = x_limits,
        ylim = y_limits,
        expand = FALSE
      ) +
      labs(
        title = title,
        subtitle = subtitle,
        x = sprintf(
          "PCoA1 (%.1f%%)",
          eig$PCoA1_positive_eigen_pct
        ),
        y = sprintf(
          "PCoA2 (%.1f%%)",
          eig$PCoA2_positive_eigen_pct
        ),
        caption = caption
      ) +
      theme_classic(
        base_size = 10
      ) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 11,
          margin = margin(b = 3)
        ),
        plot.subtitle = element_text(
          size = 8,
          color = "#4F5962",
          margin = margin(b = 7)
        ),
        axis.title = element_text(
          size = 9
        ),
        axis.text = element_text(
          size = 8
        ),
        panel.border = element_rect(
          color = "#D2DCE3",
          fill = NA,
          linewidth = 0.45
        ),
        plot.caption = element_text(
          size = 7.3,
          color = "#59656D",
          hjust = 0
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.title = element_text(
          face = "bold",
          size = 8,
          margin = margin(b = 3)
        ),
        legend.text = element_text(
          size = 8
        ),
        legend.spacing.x = grid::unit(
          0.12, "cm"
        ),
        legend.spacing.y = grid::unit(
          0.13, "cm"
        ),
        legend.key.width = grid::unit(
          0.32, "cm"
        ),
        legend.key.height = grid::unit(
          0.30, "cm"
        ),
        plot.margin = margin(
          11, 7, 7, 7
        )
      )

    return(p)
  }

  save_figure <- function(
      figure,
      filename,
      width,
      height
  ) {

    ggsave(
      filename = file.path(
        output_dir,
        paste0(filename, ".png")
      ),
      plot = figure,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )

    ggsave(
      filename = file.path(
        output_dir,
        paste0(filename, ".pdf")
      ),
      plot = figure,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )

  }

  combine_three <- function(p1, p2, p3) {

    (p1 + p2 + p3) +
      plot_layout(
        ncol = 3,
        guides = "keep"
      ) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(
          plot.tag = element_text(
            face = "bold",
            size = 14,
            color = "#182934"
          ),
          plot.tag.position = c(0.01, 1)
        )
      )

  }

  for (rank in ranks) {

    cat(
      "\nGenerating:", rank, "\n"
    )

    method <- "MetaPhlAn_GTDB"

    p_swiss <- make_pcoa(
      dataset = "Switzerland",
      method = method,
      rank = rank,
      n = 36,
      group_column = "MAP_group",
      colors = map_colors,
      title = "Swiss MAP status",
      subtitle = get_statistics(
        primary,
        rank,
        "Swiss_MAP"
      ),
      legend_columns = 3
    )

    p_farms <- make_pcoa(
      dataset = "Switzerland",
      method = method,
      rank = rank,
      n = 36,
      group_column = "Farm_plot",
      colors = farm_colors,
      title = "Swiss farms",
      subtitle = get_statistics(
        primary,
        rank,
        "Swiss_farms"
      ),
      legend_columns = 3,
      caption = paste(
        "Grey: independent/unreplicated farms;",
        "excluded from farm-specific tests."
      )
    )

    p_country <- make_pcoa(
      dataset = "Swiss + South Africa",
      method = method,
      rank = rank,
      n = 50,
      group_column = "Country",
      colors = country_colors,
      title = "Switzerland vs South Africa",
      subtitle = get_statistics(
        primary,
        rank,
        "Country_SW_vs_ZA"
      ),
      legend_columns = 2
    )

    p_za <- make_pcoa(
      dataset = "South Africa",
      method = method,
      rank = rank,
      n = 14,
      group_column = "MAP_group",
      colors = map_colors[1:2],
      title = "South African MAP status",
      subtitle = get_statistics(
        primary,
        rank,
        "SouthAfrica_MAP"
      ),
      legend_columns = 2
    )

    p_five <- make_pcoa(
      dataset = "Swiss + South Africa",
      method = method,
      rank = rank,
      n = 50,
      group_column = "Cohort_group",
      colors = cohort_colors,
      title = "Five MAP-status/cohort groups",
      subtitle = get_statistics(
        primary,
        rank,
        "Combined_5_groups"
      ),
      legend_columns = 3,
      show_country_shapes = TRUE
    )

    p_abc_swiss <- combine_three(
      p_swiss,
      p_farms,
      p_country
    )

    save_figure(
      p_abc_swiss,
      paste0(
        "PCoA_", rank,
        "_ABC_Swiss_MAP_Farms_Country"
      ),
      width = 20,
      height = 8.5
    )

    p_abc_geography <- combine_three(
      p_za,
      p_country,
      p_five
    )

    save_figure(
      p_abc_geography,
      paste0(
        "PCoA_", rank,
        "_ABC_ZA_Country_Combined"
      ),
      width = 20,
      height = 9.0
    )

    save_figure(
      p_swiss,
      paste0(
        "PCoA_", rank,
        "_Swiss_MAP"
      ),
      7.3, 6.5
    )

    save_figure(
      p_farms,
      paste0(
        "PCoA_", rank,
        "_Swiss_Farms"
      ),
      7.3, 7.2
    )

    save_figure(
      p_country,
      paste0(
        "PCoA_", rank,
        "_Country"
      ),
      7.3, 6.5
    )

    save_figure(
      p_za,
      paste0(
        "PCoA_", rank,
        "_ZA_MAP"
      ),
      7.3, 6.5
    )

    save_figure(
      p_five,
      paste0(
        "PCoA_", rank,
        "_Combined_Five_Groups"
      ),
      7.8, 7.3
    )

    method_plots <- lapply(
      methods,
      function(m) {

        make_pcoa(
          dataset = "Swiss matched (35)",
          method = m,
          rank = rank,
          n = 35,
          group_column = "MAP_group",
          colors = map_colors,
          title = method_labels[[m]],
          subtitle = get_statistics(
            method_tests,
            rank,
            "Swiss_MAP_matched35",
            method = m
          ),
          legend_columns = 3
        )

      }
    )

    p_methods <- wrap_plots(
      method_plots,
      ncol = 3,
      guides = "collect"
    ) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(
          plot.tag = element_text(
            face = "bold",
            size = 14,
            color = "#182934"
          ),
          plot.tag.position = c(0.01, 1)
        )
      ) &
      theme(
        legend.position = "bottom"
      )

    save_figure(
      p_methods,
      paste0(
        "PCoA_", rank,
        "_ABCDE_Five_Classifiers"
      ),
      width = 19,
      height = 12
    )

  }

  cat(
    "\n============================================\n",
    "DVM BETA-DIVERSITY FIGURES COMPLETE\n",
    "============================================\n",
    "Output directory:\n",
    normalizePath(output_dir),
    "\n\n",
    "Genus and species PCoAs created.\n",
    "PNG and PDF formats exported.\n",
    "Three-panel figures: A-C.\n",
    "Five-classifier figures: A-E.\n",
    "Common axis limits and ticks applied.\n",
    "Separate country shape legend included.\n",
    sep = ""
  )
})
