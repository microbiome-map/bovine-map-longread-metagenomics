#!/usr/bin/env Rscript
# Histogram of per-read divergence (minimap2 de tag) of reads mapped to MAP CP033688.1.
# Input: read_divergence.tsv from extract_read_divergence.py (sample, read, primary, mapq, aln_len, de)
# Usage: Rscript plot_read_divergence.R read_divergence.tsv read_divergence_histogram ["sample=Panel title;..."]
# Panels follow the order of samples in the input file. Writes <out>.pdf, <out>.png and <out>_summary.tsv.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else "read_divergence.tsv"
out    <- if (length(args) >= 2) args[2] else "read_divergence_histogram"

d <- read.delim(infile, colClasses = c("character", "character", "integer", "integer", "integer", "numeric"))
d <- d[d$primary == 1, ]                       # primary alignments only
d$de_pct <- 100 * d$de
# panel titles: optional third argument, e.g. "26-127=Isolate 26-127;Dia5=Metagenome Dia5"
samples <- unique(d$sample)
labels <- setNames(samples, samples)
if (length(args) >= 3) {
  for (kv in strsplit(args[3], ";")[[1]]) {
    kv <- strsplit(kv, "=", fixed = TRUE)[[1]]
    if (length(kv) == 2 && kv[1] %in% samples) labels[kv[1]] <- kv[2]
  }
}
d$panel <- factor(labels[d$sample], levels = labels[samples])

# summary table: share of reads below each threshold
s <- do.call(rbind, lapply(split(d, d$sample), function(x) data.frame(
  sample = x$sample[1], reads = nrow(x),
  median_de_pct = round(median(x$de_pct), 2),
  pct_below_3 = round(100 * mean(x$de_pct < 3), 1),
  pct_below_10 = round(100 * mean(x$de_pct < 10), 1),
  pct_10_or_more = round(100 * mean(x$de_pct >= 10), 1))))
write.table(s, paste0(out, "_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(s)

p <- ggplot(d, aes(x = pmin(de_pct, 40))) +
  geom_histogram(binwidth = 0.5, boundary = 0, fill = "grey35", colour = NA) +
  geom_vline(xintercept = c(3, 10), linetype = "dashed", linewidth = 0.4, colour = "grey10") +
  geom_text(data = data.frame(x = c(3, 10), lab = c("3%", "10%")), aes(x = x, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.5, hjust = -0.15, size = 3) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_continuous(breaks = seq(0, 40, 5), limits = c(0, 40.5), expand = c(0, 0)) +
  labs(x = "Per-read divergence from CP033688.1 (minimap2 de, %; values > 40% shown at 40%)",
       y = "Reads (log scale)") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3))
ggsave(paste0(out, ".pdf"), p, width = 170, height = 150, units = "mm")
ggsave(paste0(out, ".png"), p, width = 170, height = 150, units = "mm", dpi = 300)
