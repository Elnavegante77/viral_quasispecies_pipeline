
#!/usr/bin/env bash

#SBATCH -J 4_resultados.sh
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
#SBATCH --constraint=cal
#SBATCH --error=get_results.%J.err
#SBATCH --output=get_results.%J.out

# Parámetros que se exportan al lanzar: SAMPLE (nombre de la muestra)
export SAMPLE=${SAMPLE:-SAMPLE}  # valor por defecto SAMPLE si no se pasa
FSCRATCH=/mnt2/fscratch/users/colabscbi_bio_uma/nzdb94

# Directorios
VARSCAN_DIR=${FSCRATCH}/temp_dir_SNV/varscan_output
CLIQUESNV_DIR=${FSCRATCH}/temp_dir_SNV/cliquesnv_output
PEHAPLO_DIR=${FSCRATCH}/temp_dir_SNV/pehaplo_output
QUASITOOLS_DIR=${FSCRATCH}/temp_dir_SNV/complexity
BWA_DIR=${FSCRATCH}/temp_dir_SNV/alignment

FINAL_DIR=${FSCRATCH}/temp_dir_SNV/final_results/${SAMPLE}
mkdir -p "$FINAL_DIR"

copy_if_exists() {
    local src=$1
    local dst=$2
    if [ -e "$src" ]; then
        cp -r "$src" "$dst"
    else
        echo "Warning: $src no encontrado, se omite."
    fi
}

echo "Recopilando resultados para $SAMPLE..."

# VarScan SNPs e Indels
copy_if_exists ${VARSCAN_DIR}/${SAMPLE}_snps.vcf $FINAL_DIR/
copy_if_exists ${VARSCAN_DIR}/${SAMPLE}_indels.vcf $FINAL_DIR/

# CliqueSNV resultados (carpeta o archivos con prefijo)
copy_if_exists ${CLIQUESNV_DIR}/${SAMPLE}_cliquesnv* $FINAL_DIR/

# PEHaplo (ajusta si tienes resultados específicos)
copy_if_exists ${PEHAPLO_DIR}/${SAMPLE}_pehaplo* $FINAL_DIR/

# Quasitools complexity
copy_if_exists ${QUASITOOLS_DIR}/${SAMPLE}_complexity.txt $FINAL_DIR/

# BAM alineado
copy_if_exists ${BWA_DIR}/${SAMPLE}.sorted.bam $FINAL_DIR/

echo "Generando reporte resumen..."

REPORT=$FINAL_DIR/${SAMPLE}_summary_report.txt
echo "Informe combinado para $SAMPLE" > $REPORT
echo "===============================" >> $REPORT
echo "" >> $REPORT

if [ -f "${FINAL_DIR}/${SAMPLE}_snps.vcf" ]; then
  echo "VarScan SNPs (primeras 20 líneas):" >> $REPORT
  head -n 20 "${FINAL_DIR}/${SAMPLE}_snps.vcf" >> $REPORT
  echo "" >> $REPORT
fi

if [ -f "${FINAL_DIR}/${SAMPLE}_indels.vcf" ]; then
  echo "VarScan Indels (primeras 20 líneas):" >> $REPORT
  head -n 20 "${FINAL_DIR}/${SAMPLE}_indels.vcf" >> $REPORT
  echo "" >> $REPORT
fi

if ls ${FINAL_DIR}/${SAMPLE}_cliquesnv* 1> /dev/null 2>&1; then
  echo "Archivos CliqueSNV:" >> $REPORT
  ls ${FINAL_DIR}/${SAMPLE}_cliquesnv* >> $REPORT
  echo "" >> $REPORT
fi

if [ -f "${FINAL_DIR}/${SAMPLE}_complexity.txt" ]; then
  echo "Quasitools Complexity:" >> $REPORT
  head -n 20 "${FINAL_DIR}/${SAMPLE}_complexity.txt" >> $REPORT
  echo "" >> $REPORT
fi

echo "Reporte generado en $REPORT"
echo "Finalizado."
