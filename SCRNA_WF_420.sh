#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=budankm@ccf.org
#SBATCH --job-name=scrna_seq_workflow
#SBATCH --time=10:00:00
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=100G
#SBATCH --output=/home/budankm/isilon/Meghana/slurm/slurm_out_%j.out
#SBATCH --error=/home/budankm/isilon/Meghana/slurm/slurm_err_%j.err
#####################################################################
# scRNA-seq Workflow Script
# === Activate Conda Environment ===
source /home/budankm/miniforge3/bin/activate
conda activate bioinformatics

# === Directories ===
WORKDIR="/home/budankm/isilon/XieLab_NGSData/PIP_9_12_2024/Demultiplexed"
cd "$WORKDIR"
# === Create Output Directories ===
mkdir -p trimmed_reads fastp_qc_output raw_qc_output multiqc_report_trimmed

# === Step 1: FASTQC + FastP Trimming ===
for sample in $(ls ATAC_Meghana_*_L008_R1_001.fastq.gz); do
    sample_name=$(echo $sample | sed 's/_L008_R1_001.fastq.gz//')

    echo "Processing sample: $sample_name"

    # Run FastQC on raw reads
    fastqc -t 4 --outdir=raw_qc_output ${sample} ${sample_name}_L008_R2_001.fastq.gz

    # Decompress FASTQ reads
    gunzip -c ${sample} > trimmed_reads/${sample_name}_R1.fastq
    gunzip -c ${sample_name}_L008_R2_001.fastq.gz > trimmed_reads/${sample_name}_R2.fastq

    # Trim with FastP
    fastp -i trimmed_reads/${sample_name}_R1.fastq -o trimmed_reads/${sample_name}_R1_trimmed.fastq \
          -I trimmed_reads/${sample_name}_R2.fastq -O trimmed_reads/${sample_name}_R2_trimmed.fastq \
          --cut_front --cut_tail --qualified_quality_phred 20 --detect_adapter_for_pe \
          --length_required 30 --thread 8 \
          --html fastp_qc_output/fastp_report_${sample_name}.html \
          --json fastp_qc_output/fastp_report_${sample_name}.json

    # Run FastQC on trimmed reads
    fastqc -t 4 --outdir=fastp_qc_output trimmed_reads/${sample_name}_R1_trimmed.fastq trimmed_reads/${sample_name}_R2_trimmed.fastq

done

# === Step 2: STAR Alignment ===

# Define directories
STAR_INDEX_DIR="/home/budankm/isilon/Meghana/starindexGRCh38"
FASTQ_DIR="/home/budankm/isilon/XieLab_NGSData/250318_lh00134_0671_A22MKLCLT4/DemultiplexPIP_ATAC_MB/8fastqfiles/trimmed_reads"
OUTPUT_DIR="/home/budankm/isilon/XieLab_NGSData/250318_lh00134_0671_A22MKLCLT4/DemultiplexPIP_ATAC_MB/8fastqfiles/trimmed_reads/star_output"
mkdir -p "$OUTPUT_DIR"

# Confirm STAR is installed
echo "Checking if STAR is installed..."
if which STAR > /dev/null; then
  echo "STAR version: $(STAR --version)"
else
  echo "ERROR: STAR is not installed in your environment."
  exit 1
fi

# Run STAR alignment
for R1 in ${FASTQ_DIR}/*_R1_trimmed.fastq; do
    R2="${R1/_R1_trimmed.fastq/_R2_trimmed.fastq}"
    SAMPLE_NAME=$(basename "$R1" _R1_trimmed.fastq)
    OUTPUT_PREFIX="${OUTPUT_DIR}/${SAMPLE_NAME}_aligned"

    echo "Aligning: $SAMPLE_NAME"

    STAR --runThreadN 12 \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$OUTPUT_PREFIX" \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts \
         --outSAMattributes NH HI AS nM MD \
         --sjdbOverhang 150
done
# === Step 3: MultiQC ===

# Run MultiQC on all QC and alignment outputs
multiqc raw_qc_output fastp_qc_output "$OUTPUT_DIR" -o multiqc_report_combined


