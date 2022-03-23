#!/usr/bin/env python

""""
Este módulo va a contener las funciones externas
que serán llamadas mediante import desde la f(x)
principal
"""
from Bio import SeqIO
import os
from os.path import isfile, join
import colorama as col
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


# Function to get average of all elements in a list

def Average(lst):
    return sum(lst) / len(lst)

# Function to check if directory exists, if it doesn't, it creates it

def check_create_dir(path):                                               # Takes in a path
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)                                                    # If it doesn't exist, makes it

# Funcion to get rid of all bad quality reads

def filter_qual(reads_iniciales, q, n):
    
    print("Checking reads passing the quality filter")
    reads_buenos= []
    reads_malos = []
    for record in reads_iniciales:

        promedio = Average(record.letter_annotations["phred_quality"])                                                   # Gets average quality for every read
        enes     = record.seq.count('N')                                                                                 # Gets number of N in read sequence

        if ((promedio>=q) and (enes<=n)):                                                                                # If the read pass the quality/N filters
            reads_buenos.append(record)                                                                                  # It is appended for the Sample.filtered.fastq
        else:
            reads_malos.append(record)                                                                                   # Some QC read filter software include an option to keep the discarded reads

    print("Found %d good quality reads and %d bad quality reads" % (len(reads_buenos), len(reads_malos)))

    return reads_buenos, reads_malos                                                                                     # So, we also return them, written under control of a --keep option



# Function to put the potentially contaminated reads apart

def filter_contam(lista_buenos, contam_dir, outdir):                                                                                 # Ideally, infile should be the output of qual_filter to reduce number of comparisons (here we've got ~NxM with N reads and M contaminant fastas)

    contam_files     = [join(contam_dir, f) for f in os.listdir(contam_dir) if (isfile(join(contam_dir, f)) and (f.endswith(('fasta','faa','fa','fna'))))]

    df=pd.DataFrame(columns=['Organismo','Nº reads'])
    lista_malos      = []                                                                                                      # In this list I will store every read present in any of the contam fasta files

    outdir=str(Path(outdir).absolute())

    for i in contam_files:

        fasta=open(i)
        
        id              = i.split(sep='/')[-1].split(sep='.')[0]
        nombre_fastq    = "contam_" + id + ".fastq"
        path_contam_out = outdir + '/discarded_reads/' + nombre_fastq                                                                          # This part will take the names of the potential contaminants to name the contaminant.fastq files 
        path_filt_out   = outdir + '/final.fastq'
        path_excel_out  = outdir + '/summary.xlsx'
        path_png_out    = outdir + '/Reads_per_organism.png'
        
        for contam in SeqIO.parse(fasta, "fasta"):                                                                             # In the example, there's only one record per fasta, but this program also accepts fastas with 2+ contigs/chromosomes 

            lista_contam = []                                                                                                  # In this list we are going to store the reads that are present in the contaminant species fasta file
            
            print("Analizing reads potentially contaminated by %s" % contam.id)
            
            for read in lista_buenos:                                                                                          # Takes in the list of good quality reads
                if read.seq in contam.seq:                                                                                     # If it finds any of them in the contaminant list, it is put apart to a)write a contaminant.fastq file and b)make up a list of contaminated reads.
                    lista_contam.append(read)
                    lista_malos.append(read)
 
            print("Found %d reads potentially contaminated by %s, writing them to discarded_reads/%s" % (len(lista_contam), contam.id, nombre_fastq))
            
            df.loc[df.shape[0]]=[contam.id, len(lista_contam)]                                                                 # Add another row to the pandas dataframe woth organism and number of reads mapped to that organism

            with open(path_contam_out, 'w') as contam_handle:

                SeqIO.write(lista_contam, contam_handle, "fastq")                                                              # Write a fastq file containing the "contaminated reads" 
            contam_handle.close()

    for reads_buenos in lista_buenos:                                                                                           # Now I will filter out every read stored in lista_malos and write the filtered fastq file
        for reads_malos in lista_malos:
            if reads_buenos.seq in reads_malos.seq:
                lista_buenos.remove(reads_buenos)

    print("Final fastq file has %d clean reads, and is written to %s" % (len(lista_buenos), path_filt_out))
    print("Writing excel file with number of reads per organism")
    df.loc[df.shape[0]]=["Sample", len(lista_buenos)]                                                                           # I let the name "Sample" instead of human so it may be used for filtering other fastq files
    df.to_excel(path_excel_out)

    print("Plotting number of reads per organism")
    ax=df.plot.bar(x='Organismo', y='Nº reads', rot=0, logy=True)                                                                  # I plot the pandas DataFrame y axis in logarithmic scale because human reads are by far the vast majority
    ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha="right", fontsize=7)

    plt.show()
    plt.savefig(path_png_out)


    with open(path_filt_out, 'w') as final_handle:
        SeqIO.write(lista_buenos, final_handle, "fastq")
    final_handle.close()




# This adds a comparison method to class SeqRecord, this will allow us to compare reads in different fastq files (not implemented in biopython)
# taken from https://stackoverflow.com/questions/62065729/indexing-for-seqrecords

def equal_seqs(self, other):
    if not isinstance(other, SeqRecord):
        raise NotImplementedError('Comparsion on wrong types!')
    else:
        return self.seq == other.seq # You can change it to whatever you want.

SeqRecord.__eq__ = equal_seqs

__version__ = '1.0'