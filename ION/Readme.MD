# Intro
This is a pipeline for SARS-CoV2 sequencing data analysis.
The main script, ion.sh, takes in the output of our IonTorrent sequencer (a zip file containing all fastq files from the run)
and outputs a bunch of report files to assess if sequencing went right and which lineages/variants of concern/mutations in spike protein do we have.

Input required:
-zip file  (named as zip by default by the IonTorrent server) containing one fastq file per sample.
-reference (a NC_045512.2.fasta file containing the consensus SARS-Cov2 nucleotide sequence)
-VOC.csv   (a comma-separated value table containing the expected mutations in the spike for every Pango lineage, you may freely download my example sheet)

Outputs:
-local_tree.png      (a phylogenetic tree in PNG format showing distances among samples, generated with MAFFT-fasttree-plottree)
-multiqc_report.html (a multiQC summarizing quality parameters, contaminant reads presence and lineage prediction)
-short/report.xlsx   (a spreadsheet detailing quality parameters, lineages, and unexpected mutations on spike protein for every sample)

Files should be organized like this:
```
Main_folder
├──ion.sh
├──parser.py
├──reference
├  ├──NC_045512.2.fasta
├  ├──VOC.csv
├──Run_xxxx_folder
   ├──zip
├──Run_xxxy_folder
   ├──zip    
├──Run_xxxz_folder
   ├──zip   
 ```
And ion.sh should be launch with bash ../ion.sh inside every Run_xxx_folder.

# Usage

``` 
ion.sh [-u] [-p] <zip>
   -u    Extracts unmapped reads and assign them taxonomic information (requires kraken2 and kraken2 database)
   -p    Plots a phylogenetic tree based on FastTree output 
   <zip> Zip file containing one unpaired .fastq file per sample
```   

ion.sh uses the following programs, check them out and install them in a virtual environment (you can use the pangolin.yaml in this repository):

-[FastQC==0.11](https://github.com/s-andrews/FastQC)
-[FastP==0.23](https://github.com/OpenGene/fastp)
-[BWA-mem2==2.2.1](https://github.com/bwa-mem2/bwa-mem2)
-[SAMtools==1.6](http://www.htslib.org/)
-[Pangolin==4.1.2](https://cov-lineages.org/resources/pangolin/installation.html)
-[Kraken2==2.1](https://github.com/DerrickWood/kraken2)
-[Mafft==7.505](https://mafft.cbrc.jp/alignment/software/linux.html)
-[FastTree==2.1.11](http://www.microbesonline.org/fasttree/)
-[PlotTree](https://github.com/katholt/plotTree)

In order to copy the virtual environment I use to run this script, download the yml file and then:
```
conda create --name your_new_conda_env --file pangolin.yml
```

Only tested in Ubuntu20.04
