#!/bin/sh

for i in $(cat prj_list)                                                              #Aquí la lista de PRJNA/PRJEB que quieras desargar
do 
        conda activate eutils                                                          #Para que deje activar entornos dentro del bash hay que lanzarlo con 'bash -i bucle.sh'
	mkdir -p ./fastqs/$i
        samples=`esearch -db sra -query $i | efetch --format runinfo | grep PAIRED | cut -d ',' -f 1 | grep -E "*RR"`
        #echo $samples > lista_SRR.txt
        arr=($samples)
        echo ${arr[@]} | xargs -n 200 | while read line; do
                conda activate eutils
                echo $line | xargs prefetch -p                                         #Con esto descarga los ficheros SRA
                echo $line | xargs fasterq-dump -p -e 20 --split-files -O ./fastqs/$i  #Y con esto los transforma, hay que borrar ambos para que no pete el disco
                echo $line" descargados, activando autosnippy:"
                conda activate autosnippy 
                echo 'autosnippy activado, arrancando pipeline:'
                inic=`echo $line |head -n1| cut -d " " -f1`
                mkdir -p ./$i/$inic
		autosnippy.py -i ./fastqs/$i \
-r /home/laura/DATABASES/REFERENCES/ancestorII/MTB_ancestorII_reference.fasta \
--snpeff_database Mycobacterium_tuberculosis_h37rv \
-V /home/laura/DATABASES/Anotacion/H37rv/resistance_h37rv.vcf \
-V /home/laura/DATABASES/Anotacion/H37rv/Lineages_col.vcf \
-o ./$i/$inic
                echo 'autosnippy terminado'
                rm -r ./fastqs/$i/*.fastq
                rm -r ./$i/$inic/Stats/Coverage/*.cov
                mv    ./$i/$inic/Variants/*/*.bam /media/NASII/Datos/BAMs_TBDB
                rm -r /media/laura/4T_11/TBDB/SRA/sra/*.sra
                echo 'Carpeta aligerada, siguiente chunk:' 
done
done