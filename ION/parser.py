#!/usr/bin/env python

#############################################################
# Jorge R Grande - HUMV - Santander
# parser.py takes in a folder containing the output of ion.sh
# and outputs a short report summarizing stats from the IonTorrent run
# including, if any, unexpected mutations on spike protein 
# with their respective relative abundance.

from logging import RootLogger
import os
import pandas as pd
import numpy as np
import math
import argparse
import colorama as col
from pathlib import Path
from os import listdir
from os.path import isfile, join

def get_arguments():

    """
    Esta función permite introducir los parámetros al llamar al programa en la terminal
    """

    parser = argparse.ArgumentParser(prog = 'parser.py', description = 'Parser.py is part of the ION analysis pipeline. It takes in the resulting csv files from ion.sh and outputs the final report in EXCEL format')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input_folder', dest="input_folder", required=True, help="Required. Input folder to be summarized", type=os.path.abspath)

    input_group.add_argument('-v', '--voc_file', dest="voc_file", required=False, default="/home/usuario/Seqs/Ion/reference/VOC.csv", help="VOCs file for expected lineage spike mutations", type=os.path.abspath)

    output_group = parser.add_argument_group('Output', 'Output parameters')

    output_group.add_argument('-o', '--out_dir', dest='out_dir', required=False, help='Final output folder for reports, same as input_folder by default', type=os.path.abspath)

    arguments = parser.parse_args()

    return arguments

args = get_arguments()

col.init()    # Initializing colorama for warning messages

#################
# Taking in args:

input_folder        = args.input_folder

if args.out_dir:
    out_folder   = args.out_folder
else:
    out_folder   = input_folder

linajes_txt         = args.voc_file

Placa               = input_folder.split(sep='/')[-1]

Root                = input_folder.strip(Placa)

muestras_txt = input_folder + "/samples"           # Creo que no hace falta ninguna de estas "/"

report_txt   = input_folder + "/08_Linajes/report"

vcf_folder   = input_folder + "/06_VCF/"

out_report   = out_folder + '/report.xlsx'

out_short    = out_folder + '/short.xlsx'

print("Analizando el run ", Placa, "\n")


# ########################
# Sacar lista de muestras


muestras=[]
with open(muestras_txt) as f:
    for line in f:
        muestras.append(line.split(sep='.')[0])



# Fabricar diccionario linaje:mutaciones del spike

linajes     = pd.read_csv(linajes_txt, sep=',')
linajes     = linajes.fillna('')
diccionario = linajes.to_dict('list')


# Extraer info del report de ion.sh

d_muts      = {}
d_htzs      = {}
report = pd.read_csv(report_txt, sep='\s+|\t', on_bad_lines='skip', engine='python', skiprows=1, index_col=False, decimal=',', names=['Muestra','Cov_promedio','Bases>30x','Total_bases','%COV30x','N_total','lineage','scorpio_call','scorpio','Mutaciones'])

for n, i in enumerate(muestras):
    
    # Sacar lista de mutaciones en el spike de cada muestra
    tsv             = vcf_folder+i+'.tsv'
    df              = pd.read_csv(tsv, sep='\t')                                               # Abrir output de ivar variants
    df              = df[df['POS'].between(21563,25384)]                                       # Tomar las variantes del spike
    df              = df.query('GFF_FEATURE == "surface glycoprotein"')                        # Esto evita los INDELs, cambiar más adelante
    df              = df[df.REF_AA != df.ALT_AA]
    df['POS_AA']    = ((df['POS']-21560)/3).astype(int)                                        # Traducir posición_nt -> pos_aa en columna nueva
    df['Mutación']  = df['REF_AA'] + df['POS_AA'].astype(str) + df['ALT_AA']                   # columna nueva con formato E484A
    df['HTZ']       = df['ALT_FREQ'].map("{:,.2%}".format)                                             # columa con porcentaje de fijación de la mutación
    mutaciones      = df['Mutación'].tolist()                                                  # Pasar columna a lista con todas las mutaciones (en Spike) de la muestra
    htz             = df['HTZ'].tolist()                                                       # lo mismo
    if all(x in mutaciones for x in ['S371F','S371P']):
        for mut in range(len(mutaciones)):
            if mutaciones[mut] in ['S371F','S371P']:
                mutaciones[mut] = 'S371L'
    if all(x in mutaciones for x in ['K417T', 'K417N']):
        for mut in range(len(mutaciones)):
            if mutaciones[mut] in ['K417T', 'K417N']:
                mutaciones[mut] = 'K417T'

    # Abrir informe que sale del pipeline para extraer el linaje de cada muestra (completar para linajes raros y Unassigned)
    
    # report['scorpio']          = report['scorpio'].str.replace(r"\(|\)", "", regex=True)
    # report['linaje']           = report['scorpio'].str.rstrip('-like')                       # Genero una columna (lin) desde scorpio: (Ba.5-like)->BA.5
    
    report['linaje'] = report['lineage'].str.split('.', n=4)
    report           = report.fillna('')
    report['linaje'] = ['.'.join(map(str, l[:4])) for l in report['linaje']]
        
    # Sacar columna con mutaciones propias de cada linaje
    report['linaje']           = report['linaje'].map(diccionario)                             # Sustituyo el valor BA.5 (key) por sus mutaciones (values) en diccionario

    # Actualizar el diccionario Muestra:mutaciones
    fastq            = i+'.fastq'
    d_muts.setdefault(fastq, mutaciones)
    d_htzs.setdefault(fastq, htz)

    


# Asignar a cada muestra sus mutaciones (con el porcentae de reads de cada mutación)
report['Mutaciones'] = report['Muestra'].map(d_muts)                                                # Meter columna mutaciones:muestra
report['Percent']    = report['Muestra'].map(d_htzs)                                                # Meter columna mutaciones:muestra
report = report.fillna('')
report['Diccionario']= report[['Mutaciones','Percent']].apply(lambda data: {k:[y for x, y in zip(data[0], data[1]) if x == k] for k in data[0]}, axis=1)

# Sacar diferencia entre mutaciones propias del linaje y mutaciones observadas
report = report.drop(columns=['Bases>30x','Total_bases'])   
report['Mutaciones inesperadas'] = report[['Mutaciones','linaje']].apply(lambda x: [i for i in x[0] if i not in x[1]], axis=1)
report['Diccionario inesperada'] = report[['Diccionario','linaje']].apply(lambda x: (dict((k, x[0][k]) for k in x[0] if k not in x[1])), axis=1)

# Quitar morralla
report['Muestra']   = report['Muestra'].str.replace(r'.fastq', '', regex=True)                                                            # Quitar el .fastq en los nombres de muestra
report_short=report.drop(columns=['Mutaciones','linaje','Diccionario','Percent','Mutaciones inesperadas'])                                # Este solo valdrá cuando me crea el resultado final
report_short.columns=['Muestra','Cov_promedio','% Bases >30x COV','Bases no resueltas','Linaje','Variante','Scorpio','Mutaciones inesperadas [% de reads mutados]']
report_short=report_short.query('Linaje!="Unassigned"')
report_short=report_short.loc[report_short['% Bases >30x COV'].astype(float)>90]

report.to_excel(out_report)
report_short.to_excel(out_short)

print(report_short)

