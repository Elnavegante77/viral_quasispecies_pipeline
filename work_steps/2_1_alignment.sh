#!/usr/bin/env bash

####################################################################
## Created by the Bioinformatics Core of the University of Malaga ##
## email: bio@scbi.uma.es                                         ##
## Date: 08-11-2022                                               ##
## Version 1.0                                                    ##
####################################################################

#SBATCH -J 2_1_alignment
#SBATCH --cpus-per-task=16
#SBATCH --mem=52gb
#SBATCH --time=24:00:00
#SBATCH --constraint=cal
#SBATCH --error=align.%J.err
#SBATCH --output=align.%J.out

# Desactivar la variable LANG para evitar problemas de localización
unset LANG

# Definir el directorio de resultados
#RESULTS_DIR=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94/temp_dir_SNV/alignment
export EXPERIMENTO=SNV
export DATA=SAMPLE
RESULTS_DIR=${FSCRATCH}/temp_dir_${EXPERIMENTO}/alignment
mkdir -p $RESULTS_DIR

# Imprimir el nombre del host
hostname

# Cargar los módulos necesarios
module load bwa/0.7.5a
module load samtools/1.12

# Directorio de datos de lectura (lecturas FASTQ)
#READS_DIR=~/reads

# Archivos de lectura
#READ1=${READS_DIR}/sample_R1_001.fastq.gz
#READ2=${READS_DIR}/sample_R2_001.fastq.gz

# 
export INDEXADO =TEMPORAL_DIR/index/nodavirus_index

#!/usr/bin/env bash




# Definir las variables de entorno que voy a usar referencia 
#export GENOME_DIR=~/references
#export INPUTREF_GENOME=$GENOME_DIR/nodavirus_rna1.fasta

#export SAMPLE=SAMPLE  # Reemplazar con el nombre de la muestra
export READ1=TEMPORAL_DIR/${DATA}_fastp/${DATA}_paired_1.fastq.gz
export READ2=TEMPORAL_DIR/${DATA}_fastp/${DATA}_paired_2.fastq.gz
#export OUTPUT_DIR=~/alignments
#export OUTPUT_BAM=$OUTPUT_DIR/${SAMPLE}.bam
export OUTPUT_DIR=$RESULTS_DIR
export OUTPUT_BAM=$OUTPUT_DIR/${DATA}.bam

# Verificar si los archivos de entrada existen
if [[ ! -f "$INPUTREF_GENOME" ]]; then
  echo "ERROR: El archivo de referencia $INPUTREF_GENOME no existe."
#  exit 1
fi

if [[ ! -f "$READ1" ]]; then
  echo "ERROR: El archivo de lectura $READ1 no existe."
#  exit 1
fi

if [[ ! -f "$READ2" ]]; then
  echo "ERROR: El archivo de lectura $READ2 no existe."
#  exit 1
fi

# Crear el directorio de salida si no existe
mkdir -p $RESULTS_DIR

# Indexar el archivo de referencia si no está indexado
#if [[ ! -f "${INPUTREF_GENOME}.bwt" ]]; then
#  bwa index $INPUTREF_GENOME
#fi
time bwa mem -B 20 -A 3 -O 30 -E 3 -t 8 $INPUTREF_GENOME $READ1 $READ2 | samtools view -bS -q 30 -F 4 | samtools sort -o ${OUTPUT_BAM%.bam}.sorted.bam -  #el read.bam se crea en el directo actual 
#samtools mpileup -BQ 20  -q 30 -d10000000 -f /mnt/scratch/users/pab_001_uma/luisdiaz/quasiflow/nodavirus/RNA2_new/Sparus_5_3/bwa_0000/virus_ref.fasta reads.bam > Data.mpileup
#samtools mpileup -BQ 20  -q 30 -d10000000 -f $INPUTREF_GENOME reads.bam > Data.mpileup // referencia
samtools mpileup -BQ 20 -q 30 -d 10000000 -f $INPUTREF_GENOME ${OUTPUT_BAM%.bam}.sorted.bam > $OUTPUT_DIR/${DATA}.mpileup
if [ ! -f reads.bam ] || [ ! -f Data.mpileup ]; then  #reciclado
echo "ERROR"
# exit 1
fi


# Indexar el archivo BAM ordenado
samtools index ${OUTPUT_BAM%.bam}.sorted.bam

# Limpiar el archivo BAM no ordenado

#samtools view reads.bam | awk -F '\t' '{print $10}' | fold -w 1 | sort | uniq -c > nucleotide_count //referencia 
samtools view ${OUTPUT_BAM%.bam}.sorted.bam | awk -F '\t' '{print $10}' | fold -w 1 | sort | uniq -c > $OUTPUT_DIR/${DATA}_nucleotide_count.txt

# Limpiar el archivo BAM no ordenado
rm $OUTPUT_BAM


echo "Alignment and mpileup completed. Sorted BAM file: ${OUTPUT_BAM%.bam}.sorted.bam"
echo "Mpileup file: $OUTPUT_DIR/${DATA}.mpileup"
echo "Nucleotide count file: $OUTPUT_DIR/${DATA}_nucleotide_count.txt"