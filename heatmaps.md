Generate bed files in R:

```R
df <- avg_filter

# Use the correct columns present in avg_filter
chrom_col  <- "chrom"
start_col  <- "start"      # 0-based, BED style
end_col    <- "end"
strand_col <- "strand"
name_col   <- "gene"
expr_col   <- "TTseq"      # pick your expression measure

# Keep rows with finite expression and valid coords
keep <- is.finite(df[[expr_col]]) &
  !is.na(df[[chrom_col]]) & !is.na(df[[start_col]]) &
  !is.na(df[[end_col]])   & !is.na(df[[strand_col]])

x <- df[keep, c(chrom_col, start_col, end_col, name_col, strand_col, expr_col)]
colnames(x) <- c("chrom","start","end","name","strand","expr")

# Compute tertiles (33% and 66%)
cuts <- quantile(x$expr, probs = c(1/3, 2/3), na.rm = TRUE, names = FALSE)
low  <- subset(x, expr <= cuts[1])
mid  <- subset(x, expr >  cuts[1] & expr <= cuts[2])
high <- subset(x, expr >  cuts[2])

# Helper to write BED6 (chrom start end name score strand)
write_bed6 <- function(d, path) {
  bed <- data.frame(d$chrom, d$start, d$end, d$name, 0, d$strand, check.names = FALSE)
  write.table(bed, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

write_bed6(high, "genes_high.bed")
write_bed6(mid,  "genes_medium.bed")
write_bed6(low,  "genes_low.bed")

cat("Wrote: genes_high.bed, genes_medium.bed, genes_low.bed\n")
```

Generate heatmaps on eddie:
```bash
#!/bin/bash
# SGE directives
#$ -S /bin/bash
#$ -cwd
#$ -N heatmap_TSS_CGI_exp
#$ -pe sharedmem 4
#$ -l h_vmem=5G
#$ -l h_rt=02:00:00
#$ -l rl9=TRUE
#$ -j y

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: qsub $0 <bw1> <bw2> [BED_DIR]" >&2
  echo "Example: qsub $0 ../bw/A.bw ../bw/B.bw ./beds" >&2
  exit 1
fi

BW1="$1"
BW2="$2"
BED_DIR="${3:-.}"

# Load environment modules
. /etc/profile.d/modules.sh

# Use the Rocky Linux 9 module tree (matches -l rl9=TRUE)
module use /exports/igmm/software/etc/el9/modules
module purge
module load igmm/apps/miniconda3/24.9.2

# Initialize conda for non-interactive shells
source "$(conda info --base)/etc/profile.d/conda.sh"

# Create deeptools env if needed (persists in your home ~/.conda)
if ! conda env list | awk '{print $1}' | grep -qx deeptools; then
  echo "[INFO] Creating conda env 'deeptools'..."
  conda create -y -n deeptools -c conda-forge -c bioconda deeptools
fi

# Activate env
conda activate deeptools

# Threads from SGE
THREADS="${NSLOTS:-1}"

# Validate inputs
for f in "$BW1" "$BW2" \
         "$BED_DIR/genes_high.bed" "$BED_DIR/genes_medium.bed" "$BED_DIR/genes_low.bed"; do
  [[ -e "$f" ]] || { echo "[ERROR] Missing input: $f" >&2; exit 2; }
done

# Output names based on input bigWig basenames
B1=$(basename "${BW1%.*}")
B2=$(basename "${BW2%.*}")
OUT_PREFIX="${B1}__${B2}.CGI_TSS_exp"

MAT="${OUT_PREFIX}.mat.gz"
PNG="${OUT_PREFIX}.png"

echo "[INFO] computeMatrix on:"
echo "  -S $BW1 $BW2"
echo "  -R $BED_DIR/genes_high.bed $BED_DIR/genes_medium.bed $BED_DIR/genes_low.bed"
echo "[INFO] Threads: $THREADS"

computeMatrix reference-point \
  -S "$BW1" "$BW2" \
  -R "$BED_DIR/genes_high.bed" "$BED_DIR/genes_medium.bed" "$BED_DIR/genes_low.bed" \
  -a 20000 -b 20000 --skipZeros --missingDataAsZero \
  --referencePoint TSS \
  -p "$THREADS" \
  -o "$MAT"

plotHeatmap -m "$MAT" -out "$PNG" \
  --regionsLabel "" "" "" \

echo "[DONE] Wrote:"
ls -lh "$MAT" "$PNG"
```
