Align data to genome:

```bash
#!/bin/bash -f
#$ -S /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00
#$ -l h_vmem=32g
#$ -l rl9=TRUE

. /etc/profile.d/modules.sh
MODULEPATH=$MODULEPATH:/exports/igmm/software/etc/el7/modules

module load igmm/apps/bowtie/2.2.6
module load igmm/apps/samtools/1.6

# input name (prefix, e.g. SRR29706795)
base=$1

# file names
file1="${base}_1.fastq"
file2="${base}_2.fastq"

# bowtie2 index (mm10)
bowtieindex="/exports/igmm/eddie/gilbert-lab/jstocks/james/genome/mm10"

echo ">>> Aligning paired-end reads for ${base}"
echo ">>> R1: $file1"
echo ">>> R2: $file2"
echo ">>> Index: $bowtieindex"

# run alignment
bowtie2 -p 16 -x $bowtieindex -1 $file1 -2 $file2 -S ${base}.sam

# convert, sort, and index
samtools view -bS ${base}.sam > ${base}.bam
rm ${base}.sam

samtools flagstat ${base}.bam | tee ${base}.flagstat.txt

samtools sort -o ${base}.sorted.bam ${base}.bam
samtools index ${base}.sorted.bam

echo ">>> Finished alignment for ${base}"
```
<br>

Generate BWs:

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe sharedmem 16
#$ -l h_rt=24:00:00
#$ -l rl9=TRUE

set -euo pipefail
. /etc/profile.d/modules.sh
module load igmm/apps/apptainer/1.3.4

# Path to deepTools container SIF
SIF="/exports/eddie/scratch/$USER/containers/deeptools-3.5.5.sif"

# Optional samtools container (leave empty if using system samtools)
SAMTOOLS_SIF=""

# Path to mm10 blacklist (download it first)
BLACKLIST="mm10-blacklist.v2.bed"

# Temporary/cache directories
export APPTAINER_CACHEDIR="/exports/eddie/scratch/$USER/apptainer-cache"
export APPTAINER_TMPDIR="/exports/eddie/scratch/$USER/apptainer-tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

THREADS=${NSLOTS:-16}
BIN=500

# Check deepTools container exists
if [[ ! -f "$SIF" ]]; then
  echo "deepTools container not found at: $SIF"
  exit 2
fi

# Check blacklist exists
if [[ ! -f "$BLACKLIST" ]]; then
  echo "Blacklist file not found at: $BLACKLIST"
  echo "Please download it using:"
  echo "wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz && gunzip mm10-blacklist.v2.bed.gz"
  exit 2
fi

# Helper: index BAMs if needed
index_bam() {
  local bam="$1"
  if [[ -f "${bam}.bai" ]]; then
    return 0
  fi
  echo "Index missing for $bam — attempting to create..."
  if command -v samtools >/dev/null 2>&1; then
    samtools index -@ "$THREADS" "$bam"
  elif [[ -n "$SAMTOOLS_SIF" && -f "$SAMTOOLS_SIF" ]]; then
    apptainer exec -B "$PWD" -W "$PWD" "$SAMTOOLS_SIF" samtools index -@ "$THREADS" "$bam"
  else
    echo "No samtools available to create index for $bam. Please index your BAMs first."
    exit 1
  fi
}

shopt -s nullglob
found_any=false

for bam in *.sorted.bam; do
  found_any=true
  prefix=${bam%.sorted.bam}
  out="${prefix}_norm_${BIN}.bw"

  echo ">>> Generating BigWig for: $bam"
  echo ">>> Output: $out"

  # Ensure index exists
  index_bam "$bam"

  # Run bamCoverage inside Apptainer (with blacklist filter)
  apptainer exec -B "$PWD" -W "$PWD" "$SIF" \
    bamCoverage \
      -b "$bam" \
      -o "$out" \
      -p "$THREADS" \
      --binSize "$BIN" \
      --normalizeUsing RPKM \
      --blackListFileName "$BLACKLIST"

done

if ! $found_any; then
  echo "No *.sorted.bam files found in $(pwd)"
fi
```

<br>

or generate log2 difference BWs between data:
```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe sharedmem 16
#$ -l h_rt=24:00:00
#$ -l rl9=TRUE

set -euo pipefail
. /etc/profile.d/modules.sh
module load igmm/apps/apptainer/1.3.4

# Path to deepTools container
SIF="/exports/eddie/scratch/$USER/containers/deeptools-3.5.5.sif"

# Optional samtools container (leave empty if not used)
SAMTOOLS_SIF=""

# mm10 blacklist BED file
BLACKLIST="mm10-blacklist.v2.bed"

# Apptainer temporary/cache locations
export APPTAINER_CACHEDIR="/exports/eddie/scratch/$USER/apptainer-cache"
export APPTAINER_TMPDIR="/exports/eddie/scratch/$USER/apptainer-tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

THREADS=${NSLOTS:-16}
BIN=500

# Check for deepTools container
if [[ ! -f "$SIF" ]]; then
  echo "deepTools container not found at: $SIF"
  exit 2
fi

# Helper function: index BAMs if needed
index_bam() {
  local bam="$1"
  if [[ -f "${bam}.bai" ]]; then
    return 0
  fi
  echo "Index missing for $bam — attempting to create..."
  if command -v samtools >/dev/null 2>&1; then
    samtools index -@ "$THREADS" "$bam"
  elif [[ -n "$SAMTOOLS_SIF" && -f "$SAMTOOLS_SIF" ]]; then
    apptainer exec -B "$PWD" -W "$PWD" "$SAMTOOLS_SIF" samtools index -@ "$THREADS" "$bam"
  else
    echo "No samtools available to create index for $bam. Please index your BAMs first."
    exit 1
  fi
}

shopt -s nullglob
found_any=false

for bam4 in *_4_merged.mapq30.sorted.bam; do
  found_any=true
  prefix=${bam4%_4_merged.mapq30.sorted.bam}
  bam2="${prefix}_2_merged.mapq30.sorted.bam"
  out="${prefix}_4_vs_2_merged.mapq30_norm_${BIN}.bw"

  if [[ ! -f "$bam2" ]]; then
    echo "Skipping $prefix — matching file $bam2 not found"
    continue
  fi

  # Ensure indexes exist
  index_bam "$bam4"
  index_bam "$bam2"

  echo ">>> Comparing: $bam4 vs $bam2 (threads: $THREADS, bin: $BIN)"
  echo ">>> Output:   $out"
  echo ">>> Using blacklist: $BLACKLIST"

  # Run bamCompare inside Apptainer
  apptainer exec -B "$PWD" -W "$PWD" "$SIF" \
    bamCompare \
      -p "$THREADS" \
      --binSize "$BIN" \
      --scaleFactorsMethod None \
      --normalizeUsing RPKM \
      --operation log2 \
      --blackListFileName "$BLACKLIST" \
      -b1 "$bam4" \
      -b2 "$bam2" \
      -o "$out"
done

if ! $found_any; then
  echo "No *_4_merged.mapq30.sorted.bam files found in $(pwd)"
fi
```

<br>

Can change bin to whatever

<br>

Convert BWs to bedgraphs:

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe sharedmem 10
#$ -l h_rt=24:00:00
#$ -l rl9=TRUE

set -euo pipefail
. /etc/profile.d/modules.sh
module purge

# Load UCSC tools (for bigWigToBedGraph). Keep libpng if your UCSC build needs it.
module load igmm/apps/ucsc/v456 || module load ucsc || true
module load igmm/libs/libpng/1.6.44 || true

# Sanity check
if ! command -v bigWigToBedGraph >/dev/null 2>&1; then
  echo "bigWigToBedGraph not found. Try a different UCSC module or use a Conda env (bioconda: ucsc-bigwigtobedgraph)." >&2
  exit 2
fi

THREADS=${NSLOTS:-1}

# Pick up all .bw and .bigWig files in the current directory
shopt -s nullglob
files=( *.bw *.bigWig )
if (( ${#files[@]} == 0 )); then
  echo "No .bw or .bigWig files found in $(pwd)"
  exit 0
fi

echo "Found ${#files[@]} bigWigs. Converting with up to $THREADS concurrent jobs..."

convert_one() {
  local in="$1"
  local base="${in##*/}"
  base="${base%.bw}"
  base="${base%.bigWig}"
  local out="${base}.bedgraph"

  if [[ -s "$out" ]]; then
    echo "Skipping existing: $out"
    return 0
  fi

  echo "Converting: $in -> $out"
  bigWigToBedGraph "$in" "$out"
}

# Run up to $THREADS conversions concurrently (portable throttling)
for f in "${files[@]}"; do
  convert_one "$f" &
  # throttle
  while (( $(jobs -pr | wc -l) >= THREADS )); do
    sleep 0.5
  done
done
wait

# Count outputs robustly
shopt -s nullglob
bg=( *.bedgraph )
count=${#bg[@]}
echo "Done. Wrote $count bedGraph files in $(pwd)"
```

<br>

Count reads at genes:

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe sharedmem 20
#$ -l h_vmem=10G
#$ -l h_rt=48:00:00
#$ -l rl9=TRUE

set -euo pipefail

# Load environment
. /etc/profile.d/modules.sh
module purge
module load igmm/apps/BEDTools/2.31.1 || module load bedtools || true

# Inputs
GENE_BED="mm10.gencode.vM25.genes.standardChroms.bed"  # BED6: chr start end gene 0 strand
GENOME_FASTA="mm10.fa"
GENOME_SIZES="mm10.chrom.sizes"                         # two columns: chr size

# Output
output_dir="new_output_geneBodies"
mkdir -p "$output_dir"

THREADS=${NSLOTS:-1}

# Checks
bedtools --version || { echo "bedtools not found"; exit 2; }
[[ -f "$GENE_BED" ]]      || { echo "Gene BED not found: $GENE_BED" >&2; exit 1; }
[[ -f "$GENOME_FASTA" ]]  || { echo "Genome FASTA not found: $GENOME_FASTA" >&2; exit 1; }
[[ -f "${GENOME_FASTA}.fai" ]] || { echo "FASTA index missing: run samtools faidx $GENOME_FASTA" >&2; exit 1; }
[[ -f "$GENOME_SIZES" ]]  || { echo "Chrom sizes not found: $GENOME_SIZES" >&2; exit 1; }

# Step 1: filter gene BED to supported contigs and sort using genome order
FILTERED_GENE_BED="${GENE_BED%.bed}_filtered.bed"
SORTED_GENE_BED="${GENE_BED%.bed}_sorted.bed"

echo "Filtering gene BED against ${GENOME_FASTA}.fai ..."
awk 'BEGIN{while((getline<("'"$GENOME_FASTA"'.fai"))>0) ok[$1]=1} ok[$1]' "$GENE_BED" > "$FILTERED_GENE_BED"

echo "Sorting gene BED with genome order: $GENOME_SIZES ..."
LC_ALL=C bedtools sort -g "$GENOME_SIZES" -i "$FILTERED_GENE_BED" > "$SORTED_GENE_BED"
echo "Sorted gene BED: $SORTED_GENE_BED"

# Gather bedGraphs
shopt -s nullglob
files=( *.bedgraph )
if (( ${#files[@]} == 0 )); then
  echo "No .bedgraph files found in $(pwd)"
  exit 0
fi
echo "Found ${#files[@]} bedGraphs. Processing with up to $THREADS concurrent jobs ..."
echo "Output directory: $output_dir"

process_one() {
  local bedgraph_file="$1"
  local base="${bedgraph_file%.bedgraph}"
  local sorted_bedgraph_file="${output_dir}/${base}_sorted.bedgraph"
  local output_sum="${output_dir}/${base}_sum_over_genes.bed"

  if [[ -s "$output_sum" ]]; then
    echo "Skipping existing: $output_sum"
    return 0
  fi

  # Strip UCSC headers and sort consistently with genome order
  echo "Sorting (and stripping headers): $bedgraph_file -> $sorted_bedgraph_file"
  awk 'BEGIN{FS=OFS="\t"} !/^track/ && !/^browser/ && NF>=4 {print $1,$2,$3,$4}' "$bedgraph_file" \
    | LC_ALL=C bedtools sort -g "$GENOME_SIZES" -i - > "$sorted_bedgraph_file"

  # Sanity: ensure 4 columns
  if ! awk 'BEGIN{ok=1} NF<4 {ok=0; exit} END{exit !ok}' "$sorted_bedgraph_file"; then
    echo "ERROR: $sorted_bedgraph_file does not have >=4 columns after sorting" >&2
    return 1
  fi

  echo "Mapping sum over genes: $sorted_bedgraph_file"
  bedtools map -a "$SORTED_GENE_BED" -b "$sorted_bedgraph_file" -c 4 -o sum > "$output_sum"
  echo "Generated: $output_sum"
}

# Run with throttling
for f in "${files[@]}"; do
  process_one "$f" &
  while (( $(jobs -pr | wc -l) >= THREADS )); do
    sleep 0.5
  done
done
wait

# Count outputs
shopt -s nullglob
done_files=( "$output_dir"/*_sum_over_genes.bed )
echo "All done. Wrote ${#done_files[@]} files to $output_dir"
```

<br>

Count reads at genes +/- bin either side:

```bash
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe sharedmem 10
#$ -l h_vmem=10G
#$ -l h_rt=48:00:00
#$ -l rl9=TRUE

set -euo pipefail

# Load environment
. /etc/profile.d/modules.sh
module purge
module load igmm/apps/BEDTools/2.31.1 || module load bedtools || true

# Inputs
GENE_BED="mm10.gencode.vM25.genes.standardChroms.bed"  # BED6: chr start end gene 0 strand
GENOME_FASTA="mm10.fa"
GENOME_SIZES="mm10.chrom.sizes"                         # two columns: chr size
FLANK_BP=10000                                          # extend each gene by this many bp on both sides

# Output
output_dir="new_output_geneBodies_pm${FLANK_BP}"
mkdir -p "$output_dir"

THREADS=${NSLOTS:-1}

# Checks
bedtools --version || { echo "bedtools not found"; exit 2; }
[[ -f "$GENE_BED" ]]      || { echo "Gene BED not found: $GENE_BED" >&2; exit 1; }
[[ -f "$GENOME_FASTA" ]]  || { echo "Genome FASTA not found: $GENOME_FASTA" >&2; exit 1; }
[[ -f "${GENOME_FASTA}.fai" ]] || { echo "FASTA index missing: run samtools faidx $GENOME_FASTA" >&2; exit 1; }
[[ -f "$GENOME_SIZES" ]]  || { echo "Chrom sizes not found: $GENOME_SIZES" >&2; exit 1; }

# Step 1: filter gene BED to supported contigs
FILTERED_GENE_BED="${GENE_BED%.bed}_filtered.bed"
awk 'BEGIN{while((getline<("'"$GENOME_FASTA"'.fai"))>0) ok[$1]=1} ok[$1]' "$GENE_BED" > "$FILTERED_GENE_BED"

# Step 2: extend genes by ±FLANK_BP and sort using genome order
SLOPPED_GENE_BED="${GENE_BED%.bed}_slop${FLANK_BP}.bed"
SORTED_GENE_BED="${GENE_BED%.bed}_slop${FLANK_BP}_sorted.bed"

# bedtools slop clamps to [0, chromSize]; -s is not needed for symmetric ± extension
bedtools slop -i "$FILTERED_GENE_BED" -g "$GENOME_SIZES" -b "$FLANK_BP" > "$SLOPPED_GENE_BED"
LC_ALL=C bedtools sort -g "$GENOME_SIZES" -i "$SLOPPED_GENE_BED" > "$SORTED_GENE_BED"
echo "Using extended gene regions: $SORTED_GENE_BED"

# Gather bedGraphs
shopt -s nullglob
files=( *.bedgraph )
if (( ${#files[@]} == 0 )); then
  echo "No .bedgraph files found in $(pwd)"
  exit 0
fi
echo "Found ${#files[@]} bedGraphs. Processing with up to $THREADS concurrent jobs ..."
echo "Output directory: $output_dir"

process_one() {
  local bedgraph_file="$1"
  local base="${bedgraph_file%.bedgraph}"
  local sorted_bedgraph_file="${output_dir}/${base}_sorted.bedgraph"
  local output_sum="${output_dir}/${base}_sum_over_genes_slop${FLANK_BP}.bed"

  if [[ -s "$output_sum" ]]; then
    echo "Skipping existing: $output_sum"
    return 0
  fi

  # Strip headers (track/browser) and sort consistently with genome order
  awk 'BEGIN{FS=OFS="\t"} !/^track/ && !/^browser/ && NF>=4 {print $1,$2,$3,$4}' "$bedgraph_file" \
    | LC_ALL=C bedtools sort -g "$GENOME_SIZES" -i - > "$sorted_bedgraph_file"

  # Map: sum or mean depending on what you want
  # Sum (as before):
  bedtools map -a "$SORTED_GENE_BED" -b "$sorted_bedgraph_file" -c 4 -o sum > "$output_sum"
  # If you prefer mean across the extended region, use: -o mean

  echo "Generated: $output_sum"
}

# Run with throttling
for f in "${files[@]}"; do
  process_one "$f" &
  while (( $(jobs -pr | wc -l) >= THREADS )); do
    sleep 0.5
  done
done
wait

done_files=( "$output_dir"/*_sum_over_genes_slop${FLANK_BP}.bed )
echo "All done. Wrote ${#done_files[@]} files to $output_dir"
```

<br>

Combine files into one tsv:

```bash
# cd new_output_geneBodies_pm10000

files=( *_sum_over_genes_slop10000.bed )
[[ ${#files[@]} -gt 0 ]] || { echo "No *_sum_over_genes_slop10000.bed files found"; exit 1; }

out=combined_gene_sums_with_coords_slop10000.tsv

# Header
{
  printf "chrom\tstart\tend\tgene\tstrand"
  for f in "${files[@]}"; do
    printf "\t%s" "${f%_sum_over_genes_slop10000.bed}"
  done
  printf "\n"
} > "$out"

# Body: coords + gene + strand from first file; col7 values from each file
awk 'BEGIN{OFS="\t"}
     FNR==1 { f++ }
     {
       if (f==1) coord[FNR]=$1"\t"$2"\t"$3"\t"$4"\t"$6
       v=$7; if (v=="" || v=="." || v=="NA") v=0
       vals[FNR,f]=v
       if (FNR>max) max=FNR
     }
     END{
       for (i=1;i<=max;i++){
         printf "%s", coord[i]
         for (j=1;j<=f;j++){
           key=i SUBSEP j
           val=(key in vals)? vals[key] : 0
           printf "\t%s", val
         }
         printf "\n"
       }
     }' "${files[@]}" >> "$out"

echo "Wrote: $out"
```



