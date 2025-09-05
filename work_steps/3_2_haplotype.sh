#!/usr/bin/env bash
#SBATCH -J 3_2_haplotype_combined
#SBATCH --cpus-per-task=8
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=3_2_haplotype.%J.out
#SBATCH --error=3_2_haplotype.%J.err

unset LANG

# Cargar módulos
module load cliquesnv/2.0.3
module load samtools/1.21

# Asegurarse de que clique-snv está en PATH
export PATH=/mnt/home/soft/cliquesnv/programs/x86_64/2.0.3/bin:$PATH

# Ajustar memoria máxima de Java
export JAVA_OPTS="-Xmx400G"

# Variables
DATA=Sample
FSCRATCH=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94
BAM=${FSCRATCH}/temp_dir_SNV/alignment/${DATA}.sorted.bam
OUTDIR=${FSCRATCH}/temp_dir_SNV/cliquesnv_output
mkdir -p "$OUTDIR"

FINAL_HAPLO="${OUTDIR}/${DATA}_haplotypes_combined.txt"
> "$FINAL_HAPLO"

# Comprobar BAM
[[ ! -f "$BAM" ]] && { echo "ERROR: BAM $BAM no existe"; exit 1; }

# Preparar header para fragmentos
HEADER="$OUTDIR/header.sam"
samtools view -H "$BAM" > "$HEADER"

# Dividir BAM en fragmentos de 2 millones de lecturas
samtools view "$BAM" | split -l 1000000 - "$OUTDIR/${DATA}_part_"

# Procesar fragmentos
for part in "$OUTDIR"/${DATA}_part_*; do
    SAM_FILE="${part}.sam"
    cat "$HEADER" "$part" > "$SAM_FILE"
    FRAG_OUT="${SAM_FILE}_out"
    mkdir -p "$FRAG_OUT"
    
    echo "Ejecutando CliqueSNV en $SAM_FILE..."
    clique-snv -m snv-illumina -in "$SAM_FILE" -outDir "$FRAG_OUT" -threads 4
    if [[ $? -ne 0 ]]; then
        echo "ERROR: CliqueSNV falló en $SAM_FILE"
        exit 1
    fi

    # Combinar haplotipos: sumar frecuencias de haplotipos idénticos
    HAPLO_FILE="$FRAG_OUT/${DATA}_cliquesnv.haplo.freq"
    if [[ -f "$HAPLO_FILE" ]]; then
        awk '{freq[$1]+=$2} END {for(h in freq) print h"\t"freq[h]}' "$HAPLO_FILE" >> "${FINAL_HAPLO}.tmp"
    else
        echo "WARNING: No se encontró archivo de frecuencias en $FRAG_OUT"
    fi
done

# Sumarizar haplotipos idénticos de todos los fragmentos
awk '{freq[$1]+=$2} END {for(h in freq) print h"\t"freq[h]}' "${FINAL_HAPLO}.tmp" > "$FINAL_HAPLO"
rm -f "${FINAL_HAPLO}.tmp"

echo "Todos los fragmentos procesados correctamente."
echo "Archivo final de haplotipos combinado: $FINAL_HAPLO"
