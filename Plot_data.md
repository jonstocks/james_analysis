In R
Upload data at genes:

```R
comb <- read.table("combined_gene_sums_with_coords.fixed.tsv",
                   header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")
```

<br>

Add length columns and edit column names:

```R
gene_len <- with(comb, end - start)
gene_len[gene_len <= 0] <- NA_real_  # guard against odd entries

# Identify columns
meta_cols   <- c("chrom","start","end","gene","strand")
all_cols    <- colnames(comb)
sample_cols <- setdiff(all_cols, meta_cols)

# Helper: extract bin size (50 or 500) from a column name
get_bin <- function(nm) {
  # look for "... norm_50 ..." or "... norm_500 ..."
  m <- regexec("(?:^|_)norm_(\\d+)(?:\\b|_)?", nm)
  hit <- regmatches(nm, m)
  sapply(hit, function(x) if (length(x) >= 2) as.numeric(x[2]) else NA_real_)
}

bins <- get_bin(sample_cols)

# Warn/default if any columns don't match the pattern
if (any(is.na(bins))) {
  warning("No bin size found for columns:\n  ",
          paste(sample_cols[is.na(bins)], collapse = ", "),
          "\nDefaulting those to 50.")
  bins[is.na(bins)] <- 50
}

# Split single-sample vs ratio tracks
ratio_idx <- grepl("_vs_", sample_cols)
rpkm_idx  <- !ratio_idx

out <- comb[meta_cols]
```

<br>

Work out RPKMs of transc/ChIP data

```R
# Function to scale a matrix of sums to per-gene means with per-column bin sizes
scale_to_mean <- function(mat, bin_vec) {
  # mat: nGenes x nCols; bin_vec: length nCols
  # scale_ij = bin_j / gene_len_i
  scale_mat <- outer(gene_len, bin_vec, function(gl, bs) bs / gl)
  mat * scale_mat
}

# RPKM-like tracks (single-sample): convert sum -> mean per gene
if (any(rpkm_idx)) {
  m <- as.matrix(comb[ sample_cols[rpkm_idx] ])
  m[!is.finite(m)] <- NA_real_
  m_scaled <- scale_to_mean(m, bins[rpkm_idx])
  out[ sample_cols[rpkm_idx] ] <- m_scaled
}

# Ratio tracks: convert sum -> mean log2FC per gene (still log2 units, not RPKM)
if (any(ratio_idx)) {
  m <- as.matrix(comb[ sample_cols[ratio_idx] ])
  m[!is.finite(m)] <- NA_real_
  m_scaled <- scale_to_mean(m, bins[ratio_idx])
  out[ paste0(sample_cols[ratio_idx], "_meanLog2FC") ] <- m_scaled
}
```

<br>

Export data:

```R
write.table(out, file = "combined_gene_RPKM_and_meanLog2FC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

<br>

Produce averages:

```R
avg <- out[,1:5]

avg$ES_2_50 <- (out$ESA_2_merged.mapq30_norm_50 + out$ESB_2_merged.mapq30_norm_50)/2
avg$ES_2_500 <- (out$ESA_2_merged.mapq30_norm_500 + out$ESB_2_merged.mapq30_norm_500)/2
avg$ES_4_50 <- (out$ESA_4_merged.mapq30_norm_50 + out$ESB_4_merged.mapq30_norm_50)/2
avg$ES_4_500 <- (out$ESA_4_merged.mapq30_norm_500 + out$ESB_4_merged.mapq30_norm_500)/2
avg$ES_4v2_50 <- (out$ESA_4_vs_2_merged.mapq30_norm_50_meanLog2FC + out$ESB_4_vs_2_merged.mapq30_norm_50_meanLog2FC)/2
avg$ES_4v2_500 <- (out$ESA_4_vs_2_merged.mapq30_norm_500_meanLog2FC + out$ESB_4_vs_2_merged.mapq30_norm_500_meanLog2FC)/2
avg$NP_2_50 <- (out$NPA_2_merged.mapq30_norm_50 + out$NPB_2_merged.mapq30_norm_50)/2
avg$NP_2_500 <- (out$NPA_2_merged.mapq30_norm_500 + out$NPB_2_merged.mapq30_norm_500)/2
avg$NP_4_50 <- (out$NPA_4_merged.mapq30_norm_50 + out$NPB_4_merged.mapq30_norm_50)/2
avg$NP_4_500 <- (out$NPA_4_merged.mapq30_norm_500 + out$NPB_4_merged.mapq30_norm_500)/2
avg$NP_4v2_50 <- (out$NPA_4_vs_2_merged.mapq30_norm_50_meanLog2FC + out$NPB_4_vs_2_merged.mapq30_norm_50_meanLog2FC)/2
avg$NP_4v2_500 <- (out$NPA_4_vs_2_merged.mapq30_norm_500_meanLog2FC + out$NPB_4_vs_2_merged.mapq30_norm_500_meanLog2FC)/2
avg$TTseq <- (out$TTseq_A.mapq30_norm_50 + out$TTseq_B.mapq30_norm_50 + out$TTseq_C.mapq30_norm_50 + out$TTseq_D.mapq30_norm_50)/4
avg$H3K27Ac <- (out$H3K27Ac_A.mapq30_norm_50 + out$H3K27Ac_B.mapq30_norm_50)/2
avg$NP_RNA <- (out$SRR28255210.mapq30_norm_50 + out$SRR28255211.mapq30_norm_50 + out$SRR28255212.mapq30_norm_50 + out$SRR28255213.mapq30_norm_50)/4
```

<br>

Make some plots:

```R
ggplot(avg) + geom_point(aes(x = log10(TTseq), y = ES_2_50, alpha=.5)) +
  theme_classic()
```

<br>

Filter the data:

```
avg_filter <- filter(avg, TTseq > 0)
avg_filter <- filter(avg_filter, NP_RNA > 0)
avg_filter$size <- avg_filter$end - avg_filter$start
avg_filter <- filter(avg_filter, size > 500)
avg_filter <- avg_filter %>%
  filter(!grepl("^Gm", gene))
avg_filter <- avg_filter %>%
  filter(!grepl("chrM", chrom))
avg_filter <- avg_filter %>%
  filter(!grepl("CT010467*", gene))
```

<br>

Work out openness and transc changes:

```R
avg_filter$diff_50 <- avg_filter$ES_4v2_50 - avg_filter$NP_4v2_50
avg_filter$diff_500 <- avg_filter$ES_4v2_500 - avg_filter$NP_4v2_500
avg_filter$transc_diff <- log10(avg_filter$NP_RNA) - log10(avg_filter$TTseq)
avg_filter$open_diff <- avg_filter$NP_4v2_500 - avg_filter$ES_4v2_500
```

<br> 

Make some more plots:

```R
ggplot(avg_filter) + geom_point(aes(x = log10(TTseq), y = ES_4v2_500, alpha=.5)) +
  theme_classic()
```

<br>

Upload gene +/- 10kb data, add column info and RPKMs:

```R
comb <- read.delim("combined_gene_sums_with_coords_slop10000.tsv",
                   sep = "\t", header = TRUE,
                   check.names = FALSE, quote = "", comment.char = "")

gene_len <- with(comb, end - start)
gene_len[gene_len <= 0] <- NA_real_  # guard against odd entries

meta_cols   <- c("chrom","start","end","gene","strand")
all_cols    <- colnames(comb)
sample_cols <- setdiff(all_cols, meta_cols)

get_bin <- function(nm) {
  # look for "... norm_50 ..." or "... norm_500 ..."
  m <- regexec("(?:^|_)norm_(\\d+)(?:\\b|_)?", nm)
  hit <- regmatches(nm, m)
  sapply(hit, function(x) if (length(x) >= 2) as.numeric(x[2]) else NA_real_)
}

bins <- get_bin(sample_cols)

if (any(is.na(bins))) {
  warning("No bin size found for columns:\n  ",
          paste(sample_cols[is.na(bins)], collapse = ", "),
          "\nDefaulting those to 50.")
  bins[is.na(bins)] <- 50
}

ratio_idx <- grepl("_vs_", sample_cols)
rpkm_idx  <- !ratio_idx

out <- comb[meta_cols]

# Function to scale a matrix of sums to per-gene means with per-column bin sizes
scale_to_mean <- function(mat, bin_vec) {
  # mat: nGenes x nCols; bin_vec: length nCols
  # scale_ij = bin_j / gene_len_i
  scale_mat <- outer(gene_len, bin_vec, function(gl, bs) bs / gl)
  mat * scale_mat
}

# RPKM-like tracks (single-sample): convert sum -> mean per gene
if (any(rpkm_idx)) {
  m <- as.matrix(comb[ sample_cols[rpkm_idx] ])
  m[!is.finite(m)] <- NA_real_
  m_scaled <- scale_to_mean(m, bins[rpkm_idx])
  out[ sample_cols[rpkm_idx] ] <- m_scaled
}

# Ratio tracks: convert sum -> mean log2FC per gene (still log2 units, not RPKM)
if (any(ratio_idx)) {
  m <- as.matrix(comb[ sample_cols[ratio_idx] ])
  m[!is.finite(m)] <- NA_real_
  m_scaled <- scale_to_mean(m, bins[ratio_idx])
  out[ paste0(sample_cols[ratio_idx], "_meanLog2FC") ] <- m_scaled
}
```

<br>

Save data:

```R
write.table(out, file = "combined_gene_RPKM_and_meanLog2FC.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
```

<br>

Work out averages:

```R
avg_10 <- out[,1:5]

avg_10$ES_2_50 <- (out$ESA_2_merged.mapq30_norm_50 + out$ESB_2_merged.mapq30_norm_50)/2
avg_10$ES_2_500 <- (out$ESA_2_merged.mapq30_norm_500 + out$ESB_2_merged.mapq30_norm_500)/2
avg_10$ES_4_50 <- (out$ESA_4_merged.mapq30_norm_50 + out$ESB_4_merged.mapq30_norm_50)/2
avg_10$ES_4_500 <- (out$ESA_4_merged.mapq30_norm_500 + out$ESB_4_merged.mapq30_norm_500)/2
avg_10$ES_4v2_50 <- (out$ESA_4_vs_2_merged.mapq30_norm_50_meanLog2FC + out$ESB_4_vs_2_merged.mapq30_norm_50_meanLog2FC)/2
avg_10$ES_4v2_500 <- (out$ESA_4_vs_2_merged.mapq30_norm_500_meanLog2FC + out$ESB_4_vs_2_merged.mapq30_norm_500_meanLog2FC)/2
avg_10$NP_2_50 <- (out$NPA_2_merged.mapq30_norm_50 + out$NPB_2_merged.mapq30_norm_50)/2
avg_10$NP_2_500 <- (out$NPA_2_merged.mapq30_norm_500 + out$NPB_2_merged.mapq30_norm_500)/2
avg_10$NP_4_50 <- (out$NPA_4_merged.mapq30_norm_50 + out$NPB_4_merged.mapq30_norm_50)/2
avg_10$NP_4_500 <- (out$NPA_4_merged.mapq30_norm_500 + out$NPB_4_merged.mapq30_norm_500)/2
avg_10$NP_4v2_50 <- (out$NPA_4_vs_2_merged.mapq30_norm_50_meanLog2FC + out$NPB_4_vs_2_merged.mapq30_norm_50_meanLog2FC)/2
avg_10$NP_4v2_500 <- (out$NPA_4_vs_2_merged.mapq30_norm_500_meanLog2FC + out$NPB_4_vs_2_merged.mapq30_norm_500_meanLog2FC)/2
avg_10$TTseq <- (out$TTseq_A.mapq30_norm_50 + out$TTseq_B.mapq30_norm_50 + out$TTseq_C.mapq30_norm_50 + out$TTseq_D.mapq30_norm_50)/4
avg_10$H3K27Ac <- (out$H3K27Ac_A.mapq30_norm_50 + out$H3K27Ac_B.mapq30_norm_50)/2
avg_10$NP_RNA <- (out$SRR28255210.mapq30_norm_50 + out$SRR28255211.mapq30_norm_50 + out$SRR28255212.mapq30_norm_50 + out$SRR28255213.mapq30_norm_50)/4
```

<br>

Merge with gene data:
```
avg_filter$gene -> genes
avg_10_filter <- filter(avg_10, avg_10$gene %in% genes)

# Merge the two dataframes by "gene"
combined <- merge(
  avg_10_filter,
  avg_filter,
  by = "gene",
  suffixes = c("_10", "_full")
)
```

<br>

Work out some differences:

```R
combined$diff_10 <- combined$ES_4v2_500_10 - combined$NP_4v2_500_10
combined$transcdiff <- log10(combined$TTseq_10) - log10(combined$NP_RNA)
```

<br>

Plot some stuff

```R
ggplot(comb_some, aes(x = log10(TTseq_full), y = ES_4v2_500_10)) +
  geom_point(alpha = 0.5 ) +
  geom_smooth() +
  theme_classic() +
  theme(legend.position = "none")
```

<br>

Bin data based on TTseq and add averages at bins to plot:
  
```R
comb_binned <- comb_some %>%
  mutate(logTTseq = log10(TTseq_full)) %>%
  mutate(bin = cut(logTTseq, breaks = 100)) %>%
  group_by(bin) %>%
  summarise(
    mean_logTTseq = mean(logTTseq, na.rm = TRUE),
    mean_ES = mean(ES_4v2_500_10, na.rm = TRUE),
    count = n()
  )

ggplot() +
  geom_point(data = comb_some, aes(x = log10(H3K27Ac_full), y = ES_4v2_500_10), alpha = 0.1) +
  geom_point(data = comb_binned, aes(x = mean_logTTseq, y = mean_ES, colour = count), size = 2) +
#  geom_smooth(data = comb_binned, aes(x = mean_logTTseq, y = mean_ES), se = FALSE, method = "loess", colour = "red") +
  scale_colour_viridis_c(option = "magma") +
  theme_classic() +
  labs(colour = "Genes in bin",
       x = "log10(H3K27Ac at gene)",
       y = "ES_4v2_500_10")
```

<br>

Determine threshold in which fit changes:

```R
library(segmented)

comb_some <- comb_some %>%
  mutate(logTTseq = log10(TTseq_full))

fit <- lm(ES_4v2_500_10 ~ logTTseq, data = comb_some)

segfit <- segmented(fit, seg.Z = ~ logTTseq)

summary(segfit)

plot(segfit)

fit <- lm(ES_4v2_500_10 ~ logTTseq, data = comb_some)
segfit <- segmented(fit, seg.Z = ~ logTTseq)
plot(segfit)
```

<br>

Work out correlations at each part:

```R
# Extract the estimated breakpoint
bp <- segfit$psi[2]  # psi contains the breakpoint estimates


#Segment 1: below breakpoint
seg1 <- comb_some[comb_some$logTTseq < bp, ]

# Segment 2: above breakpoint
seg2 <- comb_some[comb_some$logTTseq >= bp, ]

cor_seg1_pearson <- cor(seg1$logTTseq, seg1$ES_4v2_500_10, method = "pearson")
cor_seg2_pearson <- cor(seg2$logTTseq, seg2$ES_4v2_500_10, method = "pearson")
cor_seg1_spearman <- cor(seg1$logTTseq, seg1$ES_4v2_500_10, method = "spearman")
cor_seg2_spearman <- cor(seg2$logTTseq, seg2$ES_4v2_500_10, method = "spearman")
cat("Segment 1 (low TTseq):\n",
    " Pearson r =", cor_seg1_pearson,
    " Spearman rho =", cor_seg1_spearman, "\n\n")

cat("Segment 2 (high TTseq):\n",
    " Pearson r =", cor_seg2_pearson,
    " Spearman rho =", cor_seg2_spearman, "\n")

```

<br>

Work out correlations as TTSeq RPKM is filtered:

```R
comb_some <- comb_some %>%
  mutate(logTTseq = log10(TTseq_full))

thresholds <- seq(min(comb_some$logTTseq, na.rm = TRUE),
                  max(comb_some$logTTseq, na.rm = TRUE),
                  length.out = 100)

# For each threshold, compute correlation using only genes with TTseq above it
corr_df <- sapply(thresholds, function(t) {
  df_sub <- comb_some %>% filter(logTTseq >= t)
  if (nrow(df_sub) > 5) {  # only compute if enough points
    cor(df_sub$logTTseq, df_sub$ES_4v2_500_10, method = "pearson", use = "pairwise.complete.obs")
  } else {
    NA
  }
})
```

<br>

Work out and plot correlations across data:

```R

library(pheatmap)
library(ggcorrplot)

# Numeric columns only
num_windows <- windows_heatmap[sapply(windows_heatmap, is.numeric)]

# Correlation matrices
cor_pearson_windows <- cor(num_windows, method = "pearson", use = "pairwise.complete.obs")
cor_spearman_windows <- cor(num_windows, method = "spearman", use = "pairwise.complete.obs")



# Pearson
pheatmap(cor_pearson_windows,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Pearson correlation (1kb windows)",
         display_numbers = TRUE,
         number_color = "black",
         fontsize_number = 8)

# Spearman
pheatmap(cor_spearman_windows,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman correlation (1kb windows)",
         display_numbers = TRUE,
         number_color = "black",
         fontsize_number = 8)


num_comb <- comb_some_heatmap[sapply(comb_some_heatmap, is.numeric)]

cor_pearson_comb <- cor(num_comb, method = "pearson", use = "pairwise.complete.obs")
cor_spearman_comb <- cor(num_comb, method = "spearman", use = "pairwise.complete.obs")

pheatmap(cor_pearson_comb,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Pearson correlation (genes/gene+10kbs)",
         display_numbers = TRUE)

pheatmap(cor_spearman_comb,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Spearman correlation (genes/gene+10kbs)",
         display_numbers = TRUE)
```

<br>

Upload data of counts at 1kb bins across genome:

```R
library(dplyr)
library(GenomicRanges)
library(ggplot2)

comb <- read.table("merged_1000bp_windows.tsv",
                   header = TRUE, sep = "\t", check.names = FALSE,
                   quote = "", comment.char = "")

colnames(comb) <- c(
  colnames(comb)[1:3],  # keep first 3 columns as-is
  gsub("(_merged|\\.mapq30).*", "", colnames(comb)[4:length(colnames(comb))])
)

windows <- comb[,1:3]
windows$ES_2_50 <- (comb$ESA_2_merged.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$ESB_2_merged.mapq30_norm_50_sum_over_1000bp_windows.bed)/2
windows$ES_2_500 <- (comb$ESA_2_merged.mapq30_norm_500_sum_over_1000bp_windows.bed + comb$ESB_2_merged.mapq30_norm_500_sum_over_1000bp_windows.bed)/2
windows$ES_4_50 <- (comb$ESA_4_merged.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$ESB_4_merged.mapq30_norm_50_sum_over_1000bp_windows.bed)/2
windows$ES_4_500 <- (comb$ESA_4_merged.mapq30_norm_500_sum_over_1000bp_windows.bed + comb$ESB_4_merged.mapq30_norm_500_sum_over_1000bp_windows.bed)/2
windows$ES_4_v_2_50 <- (comb$ESA_4_vs_2_merged.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$ESB_4_vs_2_merged.mapq30_norm_50_sum_over_1000bp_windows.bed)/2
windows$ES_4_v_2_500 <- (comb$ESA_4_vs_2_merged.mapq30_norm_500_sum_over_1000bp_windows.bed + comb$ESB_4_vs_2_merged.mapq30_norm_500_sum_over_1000bp_windows.bed)/2
windows$H3K27Ac <- (comb$H3K27Ac_A.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$H3K27Ac_B.mapq30_norm_50_sum_over_1000bp_windows.bed)/2
windows$TTseq <- (comb$TTseq_A.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$TTseq_B.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$TTseq_C.mapq30_norm_50_sum_over_1000bp_windows.bed + comb$TTseq_D.mapq30_norm_50_sum_over_1000bp_windows.bed)/4

standard_chroms <- paste0("chr", c(1:19, "X", "Y"))

windows <- windows %>%
  filter(chrom %in% standard_chroms)
```

<br>

Upload ChromHMM data:

```R
chromhmm <- read.table(
  "mESC_E14_12_dense.annotated.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  comment.char = "",
  skip = 1
)
```

<br>

Add chromatin state to binned data:

```
# convert both to GRanges
chromhmm_gr <- GRanges(
  seqnames = chromhmm$V1,
  ranges = IRanges(start = chromhmm$V2 + 1, end = chromhmm$V3),
  state = chromhmm$V4
)

windows_gr <- GRanges(
  seqnames = windows$chrom,
  ranges = IRanges(start = windows$start + 1, end = windows$end)
)

# find overlaps
hits <- findOverlaps(windows_gr, chromhmm_gr)

# extract the state for each window
windows$state <- NA
windows$state[queryHits(hits)] <- chromhmm$V4[subjectHits(hits)]

windows <- filter(windows, TTseq > 0)
```

<br>

Plot:

```R
ggplot(windows, aes(x = state, y = ES_4_v_2_500)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,      # tilt labels 45 degrees
      hjust = 1,       # right-align so they don’t overlap
      size = 8         # smaller font
    )
  )
```
