#!/usr/bin/env bash

#SBATCH -J 3_3_quasitools.sh
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --error=complexity_%x_%J.err
#SBATCH --output=complexity_%x_%J.out

# Cargar módulo si es necesario
# module load quasitools

export SAMPLE=SAMPLE   # Variable general, será definida al lanzar el job
FSCRATCH=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94
INPUT_BAM=${FSCRATCH}/temp_dir_SNV/alignment/${SAMPLE}.sorted.bam
OUTPUT_DIR=${FSCRATCH}/temp_dir_SNV/complexity

mkdir -p "$OUTPUT_DIR"

if [[ ! -f "$INPUT_BAM" ]]; then
  echo "ERROR: No se encontró el archivo BAM: $INPUT_BAM"
  exit 1
fi

quasitools complexity -v -i "$INPUT_BAM" -o "${OUTPUT_DIR}/${SAMPLE}_complexity.txt"

if [[ $? -eq 0 ]]; then
  echo "Complejidad calculada: ${OUTPUT_DIR}/${SAMPLE}_complexity.txt"
else
  echo "Error al calcular la complejidad para $SAMPLE"
  exit 1
fi
