#!/bin/bash

mkdir -p Descargado/fastq 01_FastQC 02_Trim/fastp 03_Alignment 04_BAMs 04_BAMs/stats  04_BAMs/depth 04_BAMs/bamstats 05_Unaligned/kraken2 06_VCF 07_Consensus 08_Linajes/MSA 09_Arboles 10_Reports
unzip zip -d Descargado/fastq/
rm    zip
ls -A Descargado/fastq > samples

tarea (){
    local i=$1
    fastqc -o 01_FastQC Descargado/fastq/${i} --threads 2
    fastp -i Descargado/fastq/${i} -o 02_Trim/${i} --dont_eval_duplication --cut_tail 25 --cut_window_size 5 --average_qual 20 --length_required 100 --html 02_Trim/fastp/${i%%.*}.html --thread 2
    bwa mem -t 2 -P -S -Y -M -o 03_Alignment/${i%.*}.sam ../reference/NC_045512.2.fasta  02_Trim/${i}
    samtools view -Sb 03_Alignment/${i%.*}.sam -@ 2 -o 04_BAMs/${i%.*}.bam
    samtools sort 04_BAMs/${i%.*}.bam -o 04_BAMs/${i%.*}.sort.bam --threads 2
    samtools index 04_BAMs/${i%.*}.sort.bam -@ 2 && rm 04_BAMs/${i%.*}.bam
    samtools flagstat --threads 2 04_BAMs/${i%.*}.sort.bam > 04_BAMs/bamstats/${i%.*}.bamstats
    /home/usuario/Programs/samtools-1.15/samtools coverage 04_BAMs/${i%.*}.sort.bam -o 04_BAMs/stats/${i%.*} -w 29903
    samtools view -b -f 4 04_BAMs/${i%.*}.sort.bam -o 05_Unaligned/${i%.*}.unmap.bam
    samtools fastq  05_Unaligned/${i%.*}.unmap.bam > 05_Unaligned/${i}.unmap
    samtools depth 04_BAMs/${i%.*}.sort.bam > 04_BAMs/depth/${i%.*}.tsv
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference ../reference/NC_045512.2.fasta 04_BAMs/${i%.*}.sort.bam | ivar variants -p 06_VCF/${i%.*} -q 25 -t 0.3 -m 10 -r ../reference/NC_045512.2.fasta -g ../reference/NC_045512.2.gff3
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference ../reference/NC_045512.2.fasta 04_BAMs/${i%.*}.sort.bam | ivar consensus -p 07_Consensus/${i%.*} -q 25 -t 0.8 -m 10 -n N -i ${i%.*}
    pangolin 07_Consensus/${i%.*}.fa -o 08_Linajes --outfile Linaje_${i%.*}.csv --threads 2

}

for i in `cat samples`; do tarea "$i" & done

for i in `cat samples`; do kraken2 --db /home/usuario/miniconda3/envs/tormes-1.3.0/db/standard --threads 30 --output 05_Unaligned/kraken2/${i%.*}.unmap.out --use-names  --report 05_Unaligned/kraken2/${i%.*}.unmap.report 05_Unaligned/${i}.unmap; done
for i in `cat samples`; do cat 04_BAMs/stats/${i%.*}| grep 'Mean coverage'  ; done | awk '{print $NF}' > 04_BAMs/stats/avg_cov
for i in `cat samples`; do echo ${i} >> 04_BAMs/depth/resumen && awk  '{if($3>30)print$3}' < 04_BAMs/depth/${i%.*}.tsv | wc -l >> 04_BAMs/depth/bases_cov30x && cat 04_BAMs/depth/${i%.*}.tsv | wc -l >> 04_BAMs/depth/total_bases; done
paste 04_BAMs/depth/resumen 04_BAMs/stats/avg_cov 04_BAMs/depth/bases_cov30x 04_BAMs/depth/total_bases | awk '{$5=100*$3/$4}1' > 04_BAMs/depth/stats.tsv
sed -i  '1i Muestra Cov_promedio Bases>30x Total_bases %COV30x' 04_BAMs/depth/stats.tsv 
echo N_totales > 07_Consensus/N_totales && for i in `cat samples`; do cat 07_Consensus/${i%.*}.fa | grep -o N |wc -l; done >> 07_Consensus/N_totales
paste 04_BAMs/depth/stats.tsv 07_Consensus/N_totales > 07_Consensus/report
for i in `cat samples`; do cat 08_Linajes/Linaje_${i%.*}.csv | awk 'FNR==2' 08_Linajes/*.csv > 08_Linajes/resumen.tsv && sed -i '1i Muestra,	lineage,	conflict,	ambiguity,	scorpio_call,	scorpio_supp,	scorpio_conflict,	scorpio_notes,	version,	pangolin_vers,	scorpio_vers,	constellation_vers,	is_designated,	qc_status,	qc_notes,	note' 08_Linajes/resumen.tsv; done
cut -d ',' -f 2,5 < 08_Linajes/resumen.tsv | sed 's/\,/ /g'  > 08_Linajes/linajes
paste 07_Consensus/report 08_Linajes/linajes > 08_Linajes/report

for i in 07_Consensus/*.fa; do (cat "${i}"; echo) ; done > 09_Arboles/concat.fa
mafft --thread 16 --maxiterate 1000 --globalpair 09_Arboles/concat.fa > 09_Arboles/mafft
fasttree -nt 09_Arboles/mafft > 09_Arboles/local_tree.nwk
plottree 09_Arboles/local_tree.nwk -s 8 -w 6 -l 8 -o 09_Arboles/local_tree.png -c

multiqc . --outdir 10_Reports
rm 02_Trim/*.fastq && rm 03_Alignment/*.sam 

python3  ../parser.py -i .


# TODO
# El multiQC tiene tres tipos de input: 1. 602200293409_IonCode_0130, 2. Placa_125_602200293409_IonCode_0130 y 3. Placa_125_602200296498_IonCode_0129.sort (renombrar antes de lanzar al menos con el prefijo Placa_125_)
# Generar un archivo para subir luego los fasta concatenados a GISAID