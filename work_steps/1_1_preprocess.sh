#!/usr/bin/env bash

####################################################################
## Created by the Bioinformatics Core of the University of Malaga ##
## email: bio@scbi.uma.es                                         ##
## Date: 08-11-2022                                               ##
## Version 1.1                                                    ##
####################################################################

#SBATCH -J 1_1_preprocess.sh
#SBATCH --cpus-per-task=50
#SBATCH --mem=400gb
#SBATCH --time=14:00:00
#SBATCH --error=preprocess.%J.err
#SBATCH --output=preprocess.%J.out
#SBATCH --constraint=cal

hostname

module load fastp/0.23.4
module load ribodetector/0.3.0

# Variables que deben ser definidas desde el script padre
 export DATA=SAMPLE
 export ORIG=ORIGEN_DATOS
 export DIRTMP=TEMPORAL_DIR
# export LEC_PAR1 y LEC_PAR2 ya definidos en el entorno principal

if [ "$DATA" == "NA" ]; then
    echo "No hay muestra (DATA=NA), saliendo sin hacer nada"
    echo 'sequences in total' > ribo.log
    exit 0
fi

echo "Procesando muestra: $DATA"

# Crear directorios de salida en temporal
mkdir -p $DIRTMP/${DATA}_fastp
cd $DIRTMP/${DATA}_fastp

# Definir rutas de los archivos de entrada
file1="${ORIG}/${DATA}${LEC_PAR1}"
file2="${ORIG}/${DATA}${LEC_PAR2}"

echo "Archivo R1: $file1"
ls -l $file1
echo "Archivo R2: $file2"
ls -l $file2

if [ ! -e "$file1" ] || [ ! -e "$file2" ]; then
    echo "ERROR: No existen uno o ambos ficheros de entrada"
    echo "Fichero R1: $file1"
    echo "Fichero R2: $file2"
    echo "ERROR: Ficheros de entrada no encontrados" >> $LOG
    exit 1
fi

# Ejecutar fastp
echo "Ejecutando fastp..."
time fastp \
    -i $file1 \
    -I $file2 \
    -o ${DATA}_paired_1.fastq.gz \
    -O ${DATA}_paired_2.fastq.gz \
    -q 20 \
    -l 40 \
    -g \
    --thread 50 \
    --html ${DATA}_fastp.html \
    --json ${DATA}_fastp.json

# Ejecutar ribodetector para eliminar rRNA
echo "Ejecutando ribodetector..."
time ribodetector_cpu \
    -t 50 \
    -l 75 \
    -i ${DATA}_paired_1.fastq.gz ${DATA}_paired_2.fastq.gz \
    -e rrna \
    --log ribo.log \
    -o ${DATA}_no_rRNA_1.fastq.gz ${DATA}_no_rRNA_2.fastq.gz

echo "Preprocesamiento de $DATA finalizado."

echo "General time"
times
