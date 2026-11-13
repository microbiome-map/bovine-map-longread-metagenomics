# 01_Pilot_Analysis.R for methodological validation: paired taxonomy and diversity
# See README.md for inputs, provenance and validation status.

# Pilot taxonomic composition (original and renormalised)
local({
  library(tidyverse)

  abundance_file <- file.path("data", "pilot", "Pilot_metaphlan_gtdb.txt")
  metadata_file <- file.path("data", "pilot", "Pilot_meta.txt")
  if (!file.exists(abundance_file)) stop("Missing pilot MetaPhlAn abundance table: ", abundance_file)
  if (!file.exists(metadata_file)) stop("Missing pilot metadata: ", metadata_file)

  abund <- read.delim(
    abundance_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  metadata <- read.delim(
    metadata_file,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Preserve the original pilot IDs while matching the published profile columns.
  metadata$original_sample_id <- metadata$sample_id
  metadata$sample_id <- paste0(metadata$animal_id_new, "_d", metadata$day)

  pilot_samples <- metadata$sample_id

  missing_samples <- setdiff(
    pilot_samples,
    colnames(abund)
  )

  if (length(missing_samples) > 0) {
    stop(
      paste(
        "Missing samples:",
        paste(missing_samples, collapse = ", ")
      )
    )
  }

  abund_long <- abund %>%
    select(
      clade_name,
      level,
      all_of(pilot_samples)
    ) %>%
    pivot_longer(
      cols = all_of(pilot_samples),
      names_to = "sample_id",
      values_to = "abundance_total_reads"
    ) %>%
    mutate(
      abundance_total_reads =
        as.numeric(abundance_total_reads)
    ) %>%
    left_join(
      metadata,
      by = "sample_id"
    )

  phylum <- abund_long %>%
    filter(level == 2) %>%
    mutate(
      Phylum = sub(
        ".*p__([^;]+).*",
        "\\1",
        clade_name
      ),

      Phylum = case_when(
        str_detect(Phylum, "^Bacillota") ~ "Bacillota",
        TRUE ~ Phylum
      )
    ) %>%

    group_by(
      sample_id,
      animal_id,
      day,
      Phylum
    ) %>%
    summarise(
      abundance_total_reads =
        sum(abundance_total_reads, na.rm = TRUE),
      .groups = "drop"
    )

  phylum <- phylum %>%
    group_by(sample_id) %>%
    mutate(
      classified_phylum_percent =
        sum(abundance_total_reads, na.rm = TRUE),

      abundance_normalized =
        abundance_total_reads /
        classified_phylum_percent * 100
    ) %>%
    ungroup()

  genus <- abund_long %>%
    filter(level == 6) %>%
    mutate(
      Genus = sub(
        ".*g__([^;]+).*",
        "\\1",
        clade_name
      )
    ) %>%
    filter(
      !is.na(Genus),
      Genus != ""
    )

  genus <- genus %>%
    group_by(sample_id) %>%
    mutate(
      classified_genus_percent =
        sum(abundance_total_reads, na.rm = TRUE),

      abundance_normalized =
        abundance_total_reads /
        classified_genus_percent * 100
    ) %>%
    ungroup()

  phylum_overall <- phylum %>%
    group_by(Phylum) %>%
    summarise(

      mean_total_reads_percent =
        mean(abundance_total_reads, na.rm = TRUE),

      median_total_reads_percent =
        median(abundance_total_reads, na.rm = TRUE),

      IQR_total_reads_percent =
        IQR(abundance_total_reads, na.rm = TRUE),

      min_total_reads_percent =
        min(abundance_total_reads, na.rm = TRUE),

      max_total_reads_percent =
        max(abundance_total_reads, na.rm = TRUE),

      mean_normalized_percent =
        mean(abundance_normalized, na.rm = TRUE),

      median_normalized_percent =
        median(abundance_normalized, na.rm = TRUE),

      IQR_normalized_percent =
        IQR(abundance_normalized, na.rm = TRUE),

      min_normalized_percent =
        min(abundance_normalized, na.rm = TRUE),

      max_normalized_percent =
        max(abundance_normalized, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(desc(mean_normalized_percent))

  genus_overall <- genus %>%
    group_by(Genus) %>%
    summarise(

      mean_total_reads_percent =
        mean(abundance_total_reads, na.rm = TRUE),

      median_total_reads_percent =
        median(abundance_total_reads, na.rm = TRUE),

      IQR_total_reads_percent =
        IQR(abundance_total_reads, na.rm = TRUE),

      min_total_reads_percent =
        min(abundance_total_reads, na.rm = TRUE),

      max_total_reads_percent =
        max(abundance_total_reads, na.rm = TRUE),

      mean_normalized_percent =
        mean(abundance_normalized, na.rm = TRUE),

      median_normalized_percent =
        median(abundance_normalized, na.rm = TRUE),

      IQR_normalized_percent =
        IQR(abundance_normalized, na.rm = TRUE),

      min_normalized_percent =
        min(abundance_normalized, na.rm = TRUE),

      max_normalized_percent =
        max(abundance_normalized, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(desc(mean_normalized_percent))

  phylum_day <- phylum %>%
    group_by(
      Phylum,
      day
    ) %>%
    summarise(

      mean_total =
        mean(abundance_total_reads, na.rm = TRUE),

      median_total =
        median(abundance_total_reads, na.rm = TRUE),

      IQR_total =
        IQR(abundance_total_reads, na.rm = TRUE),

      mean_normalized =
        mean(abundance_normalized, na.rm = TRUE),

      median_normalized =
        median(abundance_normalized, na.rm = TRUE),

      IQR_normalized =
        IQR(abundance_normalized, na.rm = TRUE),

      .groups = "drop"
    ) %>%

    mutate(
      day = paste0("Day", day)
    ) %>%

    pivot_wider(
      names_from = day,
      values_from = c(
        mean_total,
        median_total,
        IQR_total,
        mean_normalized,
        median_normalized,
        IQR_normalized
      ),
      names_glue = "{.value}_{day}"
    ) %>%

    arrange(
      desc(mean_normalized_Day0)
    )

  genus_day <- genus %>%
    group_by(
      Genus,
      day
    ) %>%
    summarise(

      mean_total =
        mean(abundance_total_reads, na.rm = TRUE),

      median_total =
        median(abundance_total_reads, na.rm = TRUE),

      IQR_total =
        IQR(abundance_total_reads, na.rm = TRUE),

      mean_normalized =
        mean(abundance_normalized, na.rm = TRUE),

      median_normalized =
        median(abundance_normalized, na.rm = TRUE),

      IQR_normalized =
        IQR(abundance_normalized, na.rm = TRUE),

      .groups = "drop"
    ) %>%

    mutate(
      day = paste0("Day", day)
    ) %>%

    pivot_wider(
      names_from = day,
      values_from = c(
        mean_total,
        median_total,
        IQR_total,
        mean_normalized,
        median_normalized,
        IQR_normalized
      ),
      names_glue = "{.value}_{day}"
    ) %>%

    arrange(
      desc(mean_normalized_Day0)
    )

  write.csv(
    phylum_overall,
    "PILOT_phyla_overall.csv",
    row.names = FALSE
  )

  write.csv(
    genus_overall,
    "PILOT_genera_overall.csv",
    row.names = FALSE
  )

  write.csv(
    phylum_day,
    "PILOT_phyla_day0_vs_day1.csv",
    row.names = FALSE
  )

  write.csv(
    genus_day,
    "PILOT_genera_day0_vs_day1.csv",
    row.names = FALSE
  )

  cat("\n========================================\n")
  cat("TOP PHYLA OVERALL\n")
  cat("========================================\n")

  print(
    phylum_overall %>%
      slice_head(n = 10),
    n = Inf,
    width = Inf
  )

  cat("\n========================================\n")
  cat("TOP GENERA OVERALL\n")
  cat("========================================\n")

  print(
    genus_overall %>%
      slice_head(n = 15),
    n = Inf,
    width = Inf
  )

  cat("\n========================================\n")
  cat("DAY 0 VS DAY 1 – TOP PHYLA\n")
  cat("========================================\n")

  print(
    phylum_day %>%
      slice_head(n = 10),
    n = Inf,
    width = Inf
  )

  cat("\n========================================\n")
  cat("NORMALIZATION CHECK\n")
  cat("========================================\n")

  phylum %>%
    group_by(sample_id) %>%
    summarise(
      total_original =
        sum(abundance_total_reads),
      total_normalized =
        sum(abundance_normalized),
      .groups = "drop"
    ) %>%
    print(n = Inf)

  cat("\nDone. Four CSV files generated.\n")
})

# Pilot paired alpha and beta diversity
local({
  input_dir  <- file.path("data", "pilot")
  output_dir <- "DVM_Pilot_Manuscript_Style"
  draw_pair_lines <- FALSE        # FALSE matches main-manuscript PCoAs
  set.seed(123)

  packages <- c("readr", "dplyr", "tidyr", "ggplot2", "vegan", "ape")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages first: install.packages(c('",
                            paste(missing, collapse = "', '"), "'))")
  suppressPackageStartupMessages({
    library(readr); library(dplyr); library(tidyr)
    library(ggplot2); library(vegan); library(ape)
  })
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  find_file <- function(choices) {
    paths <- file.path(input_dir, choices)
    found <- paths[file.exists(paths)]
    if (!length(found)) stop("Input missing. Checked: ", paste(paths, collapse = ", "))
    found[1]
  }
  abundance_file <- find_file("Pilot_metaphlan_gtdb.txt")
  metadata_file <- find_file("Pilot_meta.txt")

  animal_cols <- c(SW_NEG_01="#E76F51", SW_NEG_07="#2A9D8F",
                   SW_NEG_09="#457B9D", SW_NEG_11="#9D4EDD")
  day_cols <- c("Day 0"="#91BBD2", "Day 1"="#E9A18B")
  day_shapes <- c("Day 0"=16, "Day 1"=17)

  meta <- read_tsv(metadata_file, show_col_types=FALSE)
  required <- c("sample_id", "animal_id", "animal_id_new", "day")
  if (!all(required %in% names(meta))) stop("Metadata columns missing.")
  meta <- meta %>% mutate(
    original_sample_id=trimws(as.character(sample_id)),
    sample_id=paste0(trimws(as.character(animal_id_new)), "_d", day),
    animal_id_new=trimws(as.character(animal_id_new)),
    day=as.integer(day),
    Animal=factor(animal_id_new, levels=names(animal_cols)),
    Sampling_day=factor(day, levels=c(0,1), labels=c("Day 0","Day 1"))
  )
  if (nrow(meta)!=8 || anyDuplicated(meta$sample_id) ||
      anyNA(meta$Animal) || anyNA(meta$Sampling_day)) {
    stop("Expected eight unique samples from the four known animals.")
  }
  if (any(table(meta$Animal, meta$Sampling_day)!=1)) {
    stop("Each animal must have exactly one Day 0 and one Day 1 sample.")
  }
  raw <- read_tsv(abundance_file, show_col_types=FALSE)
  if (!all(c("clade_name","level") %in% names(raw))) stop("Invalid abundance table.")
  if (!all(meta$sample_id %in% names(raw))) stop("Some pilot samples are missing.")

  make_matrix <- function(level) {
    dat <- raw[raw$level==level, c("clade_name", meta$sample_id)]
    if (!nrow(dat)) stop("Taxonomic level missing: ", level)
    x <- as.matrix(dat[, meta$sample_id, drop=FALSE])
    storage.mode(x) <- "numeric"
    if (anyNA(x) || any(!is.finite(x)) || any(x<0)) stop("Invalid abundances.")
    rownames(x) <- make.unique(dat$clade_name)
    x <- t(x)   # samples x taxonomic lineages
    if (any(rowSums(x)==0)) stop("At least one sample lacks classified abundance.")
    x
  }
  matrices <- list(Genus=make_matrix(6), Species=make_matrix(7))
  write_tsv(meta, file.path(output_dir, "included_pilot_samples.tsv"))

  alpha <- bind_rows(lapply(names(matrices), function(rank) {
    x <- matrices[[rank]]
    tibble(sample_id=rownames(x), Rank=rank,
           `Observed richness`=vegan::specnumber(x),
           `Shannon index`=vegan::diversity(x, index="shannon"),
           `Gini-Simpson index`=vegan::diversity(x, index="simpson"))
  })) %>%
    left_join(meta %>% select(sample_id, Animal, Sampling_day), by="sample_id") %>%
    mutate(Rank=factor(Rank, levels=c("Genus","Species")))

  alpha_long <- alpha %>% pivot_longer(
    c(`Observed richness`, `Shannon index`, `Gini-Simpson index`),
    names_to="Metric", values_to="Value") %>%
    mutate(Metric=factor(Metric,
      levels=c("Observed richness", "Shannon index", "Gini-Simpson index")))

  paired_values <- alpha_long %>%
    select(Rank, Metric, Animal, Sampling_day, Value) %>%
    pivot_wider(names_from=Sampling_day, values_from=Value)

  alpha_tests <- paired_values %>% group_by(Rank, Metric) %>% summarise(
    n_pairs=n(), median_day0=median(`Day 0`), median_day1=median(`Day 1`),
    p_value=wilcox.test(`Day 1`, `Day 0`, paired=TRUE,
                        exact=FALSE, correct=TRUE)$p.value,
    .groups="drop")

  write_tsv(alpha, file.path(output_dir, "pilot_alpha_per_sample.tsv"))
  write_tsv(alpha_tests, file.path(output_dir, "pilot_alpha_paired_tests.tsv"))
  print(alpha_tests, n=nrow(alpha_tests))

  annotation <- alpha_long %>% group_by(Rank, Metric) %>% summarise(
    max_y=max(Value), min_y=min(Value), .groups="drop") %>%
    left_join(alpha_tests, by=c("Rank","Metric")) %>%
    mutate(x=1.5,
           y=max_y + 0.12*pmax(max_y-min_y, abs(max_y)*0.05),
           label=sprintf("Paired Wilcoxon p = %.3f", p_value))

  plot_alpha <- function(dat, ann, stem, width, height) {
    p <- ggplot(dat, aes(Sampling_day, Value)) +
      geom_boxplot(aes(fill=Sampling_day), width=0.54, alpha=0.50,
                   outlier.shape=NA, linewidth=0.5, colour="grey30") +
      geom_line(aes(group=Animal), colour="grey65", linewidth=0.45) +
      geom_point(aes(colour=Animal), size=2.2) +
      geom_text(data=ann, aes(x=x,y=y,label=label), inherit.aes=FALSE,
                colour="grey35", size=3) +
      facet_grid(Metric ~ Rank, scales="free_y", switch="y") +
      scale_fill_manual(values=day_cols) +
      scale_colour_manual(values=animal_cols) +
      scale_y_continuous(expand=expansion(mult=c(0.05,0.12))) +
      labs(x=NULL, y=NULL, fill="Sampling day", colour="Animal ID") +
      theme_bw(base_size=11) +
      theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_line(colour="grey93", linewidth=0.25),
            strip.background=element_rect(fill="grey97", colour=NA),
            strip.placement="outside", strip.text=element_text(face="bold"),
            legend.position="right")
    print(p)
    ggsave(file.path(output_dir,paste0(stem,".png")),p,
           width=width,height=height,dpi=400,bg="white")
    ggsave(file.path(output_dir,paste0(stem,".pdf")),p,
           width=width,height=height,bg="white")
    p
  }

  alpha_all <- plot_alpha(alpha_long, annotation,
                          "Pilot_Alpha_Boxplots_Genus_Species", 10.2, 8.5)
  alpha_species <- plot_alpha(filter(alpha_long, Rank=="Species"),
                              filter(annotation, Rank=="Species"),
                              "Pilot_Alpha_Boxplots_Species", 5.8, 7.6)

  pair_rows <- split(seq_len(nrow(meta)), as.character(meta$animal_id_new))
  if (length(pair_rows)!=4 || any(lengths(pair_rows)!=2)) stop("Pairing error.")
  flips <- as.matrix(expand.grid(rep(list(c(FALSE,TRUE)), 4)))
  flips <- flips[rowSums(flips)>0,,drop=FALSE]
  permutation_matrix <- matrix(rep(seq_len(nrow(meta)), nrow(flips)),
                               nrow=nrow(flips), byrow=TRUE)
  for (j in seq_along(pair_rows)) {
    r <- which(flips[,j])
    permutation_matrix[r, pair_rows[[j]][1]] <- pair_rows[[j]][2]
    permutation_matrix[r, pair_rows[[j]][2]] <- pair_rows[[j]][1]
  }
  stopifnot(nrow(permutation_matrix)==15)

  theme_pcoa <- function() {
    theme_classic(base_size=10) +
      theme(
        axis.line=element_line(colour="grey45",linewidth=0.38),
        axis.ticks=element_line(colour="grey45",linewidth=0.30),
        axis.text=element_text(colour="grey30",size=8.7),
        axis.title=element_text(colour="black",size=10),
        plot.title=element_text(face="bold",size=11,hjust=0),
        plot.subtitle=element_text(size=8.4,colour="grey40",hjust=0,
                                   margin=margin(b=6)),
        legend.position="bottom",
        legend.box="vertical",
        legend.title=element_text(face="plain",size=9),
        legend.text=element_text(size=8.5),
        legend.key.size=grid::unit(0.36,"cm"),
        legend.spacing.y=grid::unit(0.04,"cm"),
        plot.margin=margin(8,10,8,8)
      )
  }

  run_beta <- function(rank) {
    x <- matrices[[rank]]
    bray <- vegan::vegdist(x, method="bray")
    ord <- ape::pcoa(bray)  # same ordination as the earlier pilot
    explained <- 100*ord$values$Relative_eig[1:2]

    scores <- as.data.frame(ord$vectors[,1:2,drop=FALSE])
    names(scores) <- c("PCoA1","PCoA2")
    scores$sample_id <- rownames(scores)
    scores <- scores %>% left_join(
      meta %>% select(sample_id,Animal,Sampling_day), by="sample_id") %>%
      arrange(Animal,Sampling_day)

    paired <- vegan::adonis2(bray ~ Animal + Sampling_day,
                             data=meta, permutations=permutation_matrix,
                             by="margin")
    set.seed(123)
    unblocked <- vegan::adonis2(bray ~ Sampling_day,
                                data=meta,permutations=999,by="terms")

    extract <- function(model, design) {
      i <- match("Sampling_day",rownames(model))
      if (is.na(i)) stop("Cannot find Sampling_day in adonis2 results")
      tibble(Rank=rank,Analysis=design,
             R2=as.numeric(model$R2[i]),
             F=as.numeric(model$F[i]),
             p_value=as.numeric(model$`Pr(>F)`[i]))
    }
    stats <- bind_rows(
      extract(paired,"Animal-adjusted paired PERMANOVA"),
      extract(unblocked,"Unblocked PERMANOVA (legacy only)"))
    print(stats)

    axis_limit <- max(abs(c(scores$PCoA1,scores$PCoA2)))*1.16
    axis_limit <- max(axis_limit,0.10)

    p <- ggplot(scores,aes(x=PCoA1,y=PCoA2)) +
      geom_hline(yintercept=0,linewidth=0.4,colour="grey90") +
      geom_vline(xintercept=0,linewidth=0.4,colour="grey90")
    if (draw_pair_lines) {
      p <- p + geom_line(aes(group=Animal), colour="grey75",
                         linewidth=0.35, alpha=0.70)
    }
    p <- p +
      geom_point(aes(colour=Animal,shape=Sampling_day),
                 size=2.3,stroke=0.35) +
      scale_colour_manual(values=animal_cols,drop=FALSE) +
      scale_shape_manual(values=day_shapes,drop=FALSE) +
      coord_fixed(xlim=c(-axis_limit,axis_limit),
                  ylim=c(-axis_limit,axis_limit),expand=FALSE) +
      labs(title=paste0("Pilot study — ",rank),
           subtitle=sprintf("n = 4 animals (8 samples) | Paired PERMANOVA: R² = %.3f; p = %.3f",
                            stats$R2[1],stats$p_value[1]),
           x=sprintf("PCoA1 (%.1f%%)",explained[1]),
           y=sprintf("PCoA2 (%.1f%%)",explained[2]),
           colour="Animal ID",shape="Sampling day") +
      guides(shape=guide_legend(order=1,nrow=1),
             colour=guide_legend(order=2,nrow=2,byrow=TRUE)) +
      theme_pcoa()

    print(p)
    ggsave(file.path(output_dir,paste0("Pilot_PCoA_",rank,"_Manuscript.png")),
           p,width=5.7,height=5.7,dpi=400,bg="white")
    ggsave(file.path(output_dir,paste0("Pilot_PCoA_",rank,"_Manuscript.pdf")),
           p,width=5.7,height=5.7,bg="white")

    write_tsv(scores,file.path(output_dir,paste0("Pilot_PCoA_Scores_",rank,".tsv")))
    write_tsv(stats,file.path(output_dir,paste0("Pilot_PERMANOVA_",rank,".tsv")))
    pair_dist <- meta %>% group_by(animal_id_new) %>% summarise(
      Day0=sample_id[day==0],Day1=sample_id[day==1],.groups="drop") %>%
      mutate(Rank=rank,Bray_Curtis=mapply(function(a,b) as.matrix(bray)[a,b],Day0,Day1))
    write_tsv(pair_dist,file.path(output_dir,paste0("Pilot_Paired_Distances_",rank,".tsv")))
    list(plot=p,stats=stats)
  }

  beta_genus <- run_beta("Genus")
  beta_species <- run_beta("Species")
  write_tsv(bind_rows(beta_genus$stats,beta_species$stats),
            file.path(output_dir,"Pilot_PERMANOVA_Combined.tsv"))

  if (requireNamespace("patchwork",quietly=TRUE)) {
    combined <- (beta_genus$plot + beta_species$plot +
                   patchwork::plot_layout(ncol=2,guides="collect") +
                   patchwork::plot_annotation(tag_levels="A")) &
                theme(legend.position="bottom",
                      plot.tag=element_text(size=11,face="bold"))
    print(combined)
    ggsave(file.path(output_dir,"Pilot_PCoA_Genus_Species_AB.png"),
           combined,width=11.5,height=5.8,dpi=400,bg="white")
    ggsave(file.path(output_dir,"Pilot_PCoA_Genus_Species_AB.pdf"),
           combined,width=11.5,height=5.8,bg="white")
  } else {
    message("To also export the two-panel A/B PCoA, install.packages('patchwork')")
  }

  cat("\nPILOT REPLOT FINISHED. Files saved in: ",normalizePath(output_dir),"\n")
  cat("PCoA statistics and coordinates are unchanged by this visual redesign.\n")
})
