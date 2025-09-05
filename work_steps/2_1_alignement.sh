#!/usr/bin/env bash

####################################################################
## Created by the Bioinformatics Core of the University of Malaga ##
## email: bio@scbi.uma.es                                         ##
## Date: 08-11-2022                                               ##
## Version 1.1                                                    ##
####################################################################

#SBATCH -J 2_1_alignment
#SBATCH --cpus-per-task=16
#SBATCH --mem=52gb
#SBATCH --time=24:00:00
#SBATCH --constraint=cal
#SBATCH --error=align.%J.err
#SBATCH --output=align.%J.out

unset LANG

export EXPERIMENTO=SNV
export DATA=SAMPLE
export TEMPORAL_DIR=${FSCRATCH}/temp_dir_${EXPERIMENTO}
RESULTS_DIR=${FSCRATCH}/temp_dir_${EXPERIMENTO}/alignment
mkdir -p $RESULTS_DIR

hostname

module load bwa/0.7.5a
module load samtools/1.12

# Variables importantes
INDEX_DIR=${TEMPORAL_DIR}/index
INDEX_BASE=${INDEX_DIR}/nodavirus_index
INPUTREF_GENOME=$GENOME_DIR/nodavirus_rna1.fasta

READ1=${TEMPORAL_DIR}/${DATA}_fastp/${DATA}_paired_1.fastq.gz
READ2=${TEMPORAL_DIR}/${DATA}_fastp/${DATA}_paired_2.fastq.gz

OUTPUT_BAM=${RESULTS_DIR}/${DATA}.bam
OUTPUT_BAM_SORTED=${RESULTS_DIR}/${DATA}.sorted.bam

# Comprobaciones de archivos
if [[ ! -f "$INPUTREF_GENOME" ]]; then
  echo "ERROR: El archivo de referencia $INPUTREF_GENOME no existe."
  exit 1
fi

if [[ ! -f "$READ1" ]]; then
  echo "ERROR: El archivo de lectura $READ1 no existe."
  exit 1
fi

if [[ ! -f "$READ2" ]]; then
  echo "ERROR: El archivo de lectura $READ2 no existe."
  exit 1
fi

# Ejecutar bwa mem con los parámetros dados y piping a samtools para filtrar y ordenar
time bwa mem -B 20 -A 3 -O 30 -E 3 -t 16 $INDEX_BASE $READ1 $READ2 | \
  samtools view -bS -q 30 -F 4 - | \
  samtools sort -o $OUTPUT_BAM_SORTED -

if [ $? -ne 0 ]; then
  echo "ERROR: bwa mem o samtools falló."
  exit 1
fi

# Crear mpileup
time samtools mpileup -BQ 20 -q 30 -d 10000000 -f $INPUTREF_GENOME $OUTPUT_BAM_SORTED > ${RESULTS_DIR}/${DATA}.mpileup

if [ $? -ne 0 ]; then
  echo "ERROR: samtools mpileup falló."
  exit 1
fi

# Indexar BAM ordenado
samtools index $OUTPUT_BAM_SORTED

# Contar nucleótidos
samtools view $OUTPUT_BAM_SORTED | awk -F '\t' '{print $10}' | fold -w 1 | sort | uniq -c > ${RESULTS_DIR}/${DATA}_nucleotide_count.txt

# Limpiar BAM no ordenado
if [ -f "$OUTPUT_BAM" ]; then
  rm $OUTPUT_BAM
fi

echo "Alignment and mpileup completed."
echo "Sorted BAM file: $OUTPUT_BAM_SORTED"
echo "Mpileup file: ${RESULTS_DIR}/${DATA}.mpileup"
echo "Nucleotide count file: ${RESULTS_DIR}/${DATA}_nucleotide_count.txt"
