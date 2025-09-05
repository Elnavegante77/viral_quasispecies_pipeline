#!/usr/bin/env bash

#SBATCH -J 3_1_varscan
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=04:00:00
#SBATCH --error=varscan.%J.err
#SBATCH --output=varscan.%J.out

module load java/1.8.0  # VarScan necesita Java
module load varscan/2.4.6

export DATA=SAMPLE  
FSCRATCH=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94
MPIL_FILE=${FSCRATCH}/temp_dir_SNV/alignment/${DATA}.mpileup
OUT_DIR=${FSCRATCH}/temp_dir_SNV/varscan_output

mkdir -p "$OUT_DIR"

if [[ ! -s "$MPIL_FILE" ]]; then
    echo "ERROR: Archivo mpileup $MPIL_FILE no existe o está vacío."
    exit 1
fi

# Llamado de SNPs
java -jar $(which varscan.jar) mpileup2snp "$MPIL_FILE" \
  --min-coverage 10 --min-reads2 4 --min-avg-qual 20 --min-var-freq 0.01 --p-value 0.01 --output-vcf 1 > "$OUT_DIR/${DATA}_snps.vcf"

# Llamado de INDELs
java -jar $(which varscan.jar) mpileup2indel "$MPIL_FILE" \
  --min-coverage 10 --min-reads2 4 --min-avg-qual 20 --min-var-freq 0.01 --p-value 0.01 --output-vcf 1 > "$OUT_DIR/${DATA}_indels.vcf"

# Comprobación de salida
if [[ -s "$OUT_DIR/${DATA}_snps.vcf" && -s "$OUT_DIR/${DATA}_indels.vcf" ]]; then
    echo "VarScan completado. Archivos VCF generados:"
    echo "- SNPs:   $OUT_DIR/${DATA}_snps.vcf"
    echo "- Indels: $OUT_DIR/${DATA}_indels.vcf"
else
    echo "ERROR: Uno o más archivos VCF están vacíos o no se generaron."
    exit 1
fi
