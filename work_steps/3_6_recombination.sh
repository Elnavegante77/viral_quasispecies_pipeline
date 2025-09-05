#!/bin/bash
#SBATCH --job-name=ViReMa_S5
#SBATCH --output=slurm-ViReMa_S5_%j.out
#SBATCH --error=slurm-ViReMa_S5_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00

# Cargar Python si es necesario
module load python/3.8

# Variables de referencia y archivos de entrada
export INPUTREF_GENOME=~/references/nodavirus_rna2.fasta
export READS1=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94/temp_dir_SNV/Sample5_S5_fastp/Sample5_S5_paired_1.fastq.gz
export READS2=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94/temp_dir_SNV/Sample5_S5_fastp/Sample5_S5_paired_2.fastq.gz
export VIREMA_OUT=~/quasiflow/virema/Sample5

# Crear carpeta de salida
mkdir -p $VIREMA_OUT

# Ejecutar ViReMa
python3 /mnt/home/soft/virema/programs/x86_64/0.29/ViReMa.py \
    "$INPUTREF_GENOME" \
    "$READS1,$READS2" \
    "$VIREMA_OUT" \
    --Defuzz 0.01 \
    --N 2 \
    --X 10 \
    --p 4 \
    --Output_Tag Sample5_S5 \
    --Aligner_Directory /mnt/home/soft/bowtie/programs/x86_64/1.3.1/

echo "? ViReMa finalizado para Sample5"
