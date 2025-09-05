#!/usr/bin/env bash

####################################################################
## Created by the Bioinformatics Core of the University of Malaga ##
## email: bio@scbi.uma.es                                         ##
## Date: 08-11-2022                                               ##
## Version 1.1                                                    ##
####################################################################

#SBATCH -J 1_2_unify_index_${EXPERIMENTO}
#SBATCH --cpus-per-task=16
#SBATCH --mem=52gb
#SBATCH --time=24:00:00
#SBATCH --constraint=cal
#SBATCH --error=index.%J.err
#SBATCH --output=index.%J.out

module load bowtie/v1_0.12.8
unset LANG

if [ "$SAMPLE" == "NA" ] || [ -z "$EXPERIMENTO" ]; then
    echo 'Finished succesfully' > index_NA.out
    exit 0
fi

# Comprobar que las variables necesarias están definidas
if [ -z "$TEMPORAL_DIR" ] || [ -z "$REF_GENOME" ]; then
    echo "ERROR: Las variables TEMPORAL_DIR o REF_GENOME no están definidas."
    exit 1
fi

# Crear directorio para índice si no existe
mkdir -p ${TEMPORAL_DIR}/index

cd ${TEMPORAL_DIR}/index || { echo "ERROR: No se pudo acceder al directorio de índice"; exit 1; }

hostname

# Nombre base para el índice
INDEX_BASE="nodavirus_index"

echo "Indexando el genoma de referencia con bowtie-build..."
bowtie-build $REF_GENOME ${INDEX_BASE} > bowtie_build.log 2>&1

if [ $? -eq 0 ]; then
    echo "Índices del genoma de referencia creados correctamente en ${TEMPORAL_DIR}/index/${INDEX_BASE}."
else
    echo "Error al crear los índices del genoma de referencia. Revisa bowtie_build.log"
    exit 1
fi
