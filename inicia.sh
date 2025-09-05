module load aps
export EXPERIMENTO=SNV

export INPUTDIR=~/quasiflow
export INPUTFILE=${INPUTDIR}/entradas.lis
export ORIGEN_DATOS=~/reads
export TEMPORAL_DIR=${FSCRATCH}/temp_dir_${EXPERIMENTO}

export LEC_PAR1=_R1_001.fastq.gz
export LEC_PAR2=_R2_001.fastq.gz

export LOG=$INPUTDIR/envios.log

export JOB_SOURCE=${INPUTDIR}/work_steps


export GENOME_DIR=~/references
export INPUTREF_GENOME=$GENOME_DIR/nodavirus_rna1.fasta
export ANNOT_GENOME=$GENOME_DIR/RNA1.gff3
export GTF=   # $GENOME_DIR/gtf_file

