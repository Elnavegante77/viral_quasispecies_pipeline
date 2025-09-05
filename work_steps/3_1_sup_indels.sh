#!/usr/bin/env bash

#SBATCH -J 3_1_varscan
#SBATCH --cpus-per-task=16
#SBATCH --mem=400gb
#SBATCH --time=24:00:00
#SBATCH --error=varscan.%J.err
#SBATCH --output=varscan.%J.out

# Cargar módulo de VarScan
module load varscan/2.4.6

# Variables de entrada
export DATA=Sample
FSCRATCH=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94
MPIL_FILE=${FSCRATCH}/temp_dir_SNV/alignment/${DATA}.mpileup
OUT_DIR=${FSCRATCH}/temp_dir_SNV/varscan_output2

mkdir -p $OUT_DIR

# Llamado de SNPs
varscan mpileup2snp $MPIL_FILE \
  --min-coverage 10 \
  --min-reads2 4 \
  --min-avg-qual 20 \
  --min-var-freq 0.01 \
  --output-vcf 1 \
  > $OUT_DIR/${DATA}_snps.vcf

# Llamado de INDELs
varscan mpileup2indel $MPIL_FILE \
  --min-coverage 10 \
  --min-reads2 4 \
  --min-avg-qual 20 \
  --min-var-freq 0.01 \
  --output-vcf 1 \
  > $OUT_DIR/${DATA}_indels.vcf

# Comprobación
if [[ -s "$OUT_DIR/${DATA}_snps.vcf" && -s "$OUT_DIR/${DATA}_indels.vcf" ]]; then
    echo "VarScan completado. Archivos VCF generados:"
    echo "- SNPs:   $OUT_DIR/${Sample7_S7}_snps.vcf"
    echo "- Indels: $OUT_DIR/${Sample7_S7}_indels.vcf"
else
    echo "ERROR: Uno o más archivos VCF están vacíos o no se generaron."
  #  exit 1
fi
