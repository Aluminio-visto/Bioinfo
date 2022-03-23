#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
import argparse
import os
from utils import filter_qual, filter_contam, check_create_dir
from pathlib import Path
from os import listdir
from os.path import isfile, join
import colorama as col


"""
Este módulo contiene la parte principal del programa filtrador.py
Acepta como inputs:
-Fichero .fastq                                     (str)
-Directorio con fastas de potenciales contaminantes (str)
-Número máximo de Ns aceptados por read en el fastq (int)
-Valor mínimo de calidad promedio aceptado por read (int)
-Directorio de salida                               (str)
-Importa módulos propios                             (x)

Filtra según:
-Número de N                                         (x)
-Calidad media                                       (x)
-Presencia de contaminantes                          (x)

Saca:
-Pandas.df con nº reads/organismo contaminante       (x)
-Excel     con nº reads/organismo contaminante       (x)
-Barplot   con nº reads/organismo contaminante       (x)
-Fastq filtrado                                      (x)

Avisa de:
-Un mensaje verde si un parámetro falla.             (x)
-Los pasos que va dando.                             (x)
-Los reads iniciales, intermedios y finales          (x)
-Si el nº reads =/= entre fastq, excel, plot...      (x)

En adelante seguiré en inglés porque quiero subirlo luego a 
mi github por temas de curriculum, espero que no sea molestia.
"""


def get_arguments():
    """
    Esta función permite introducir los parámetros al llamar al programa en la terminal
    """

    parser = argparse.ArgumentParser(prog = 'filtrador.py', description = 'Programa para filtrar archivos fastq de lecturas cortas o largas. Permite filtrar por cantidad de Ns, calidad promedio y contaminantes')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input_file', dest="input_file", type=str, required=True, help="Required. Input fastq file to be filtered")

    input_group.add_argument('-c', '--contaminant_dir', dest="contaminant_dir", type=str, required=True, help="Subfolder containing a collection of .fasta files with potential contaminant sequences that should be removed from fastq file")

    input_group.add_argument('-n', '--max_n', dest="max_n", type=int, required=True, help="Required. Maximum number of Ns per read accepted")

    input_group.add_argument('-q', '--min_qual', dest="min_qual", type=int, required=True, help="Required. Minimum average quality per read")

    output_group = parser.add_argument_group('Output', 'Output parameters')

    output_group.add_argument('-o', '--out_dir', dest='out_dir', type= str, required=False, default='.', help='Final output folder for filtered fastqs, plots and tables. Inside this folder, filtrador.py will create a subfolder for discarded reads')

    output_group.add_argument('-k', '--keep', dest="keep", required=False, action='store_true', help="Make an extra fastq files with the discarded, bad quality reads")

    arguments = parser.parse_args()

    return arguments

args = get_arguments()

col.init()    # Initializing colorama for warning messages


# Now, let's parse the input arguments for ensuring A) paths are absolute to avoid potential conflicts and B) datatypes are consistent, even when already specified in get_arguments.

input_file       = Path(args.input_file).absolute()

contam_dir       = Path(args.contaminant_dir).absolute()

outdir           = Path(args.out_dir)

max_n            = int(args.max_n)

min_q            = int(args.min_qual)


# Now, we are going to open the input file and parse the reads:

def filter_reads(infile, outdir, n, q, contam_dir):

    infile=str(infile)

    if os.path.isfile(infile):
        print("Parsing and filtering %s fastq file" % infile)
    else:
        print(col.Fore.GREEN + "Input Fastq file not found")                               # Checks input_file exists:

    if not infile.endswith(('fq','fastq','fq.gz','fastq.gz')):
        print(col.Fore.GREEN +  "Input file must be a fastq file, try again")              # Checks input_file is, indeed, a fastq file:

    if not os.path.isdir(outdir):
        print(col.Fore.GREEN +  "Output directory does not seem to be an actual directory, I'll try to create it:")              # Checks output_dir is, indeed, an existing dir:

        try:
            os.mkdir(outdir)
        except OSError:
            print(col.Fore.GREEN + "Output directory could not be created, please try setting another direction or creating it manually")
        else:
            print(col.Fore.WHITE + "Output directory created successfully")



    nombre          = infile.split(sep='/')[-1]
    nombre          = nombre.split(sep='.')[0]
    discarded       = str(Path(outdir).absolute()) + '/discarded_reads/'
    out_discard     = str(Path(outdir).absolute()) + '/discarded_reads/' + nombre + '.low_q.fastq'
    input_handle    = open(infile, 'r')

    check_create_dir(discarded)

    reads_iniciales = [record for record in SeqIO.parse(input_handle, "fastq")]
    print("%d initial reads were found in %s fastq file." % (len(reads_iniciales), nombre))                                   # Check number of initial reads in fastq file

    reads_buenos, reads_malos    = filter_qual(reads_iniciales, q, n)                                                         # Filter fastq file by quality and 

    if args.keep == True:                                                                                                     # We may want to keep them for further applications
        with open(out_discard, 'w') as discard:

            SeqIO.write(reads_malos, discard, "fastq")
            print("%d low quality reads were stored in %s" % (len(reads_malos), out_discard))
        discard.close()

    print("%d reads passed the quality filter, a %s from the initial reads count" % (len(reads_buenos), "{:.2%}".format(len(reads_buenos)/len(reads_iniciales))))
    print("Now we check contaminant sequences presence")

    filter_contam(reads_buenos, contam_dir, outdir)


if __name__=="__main__":
    filter_reads(input_file, outdir, max_n, min_q, contam_dir)


#print('Versión:', utils.__version__)
#print('Nombre:', utils.__name__)