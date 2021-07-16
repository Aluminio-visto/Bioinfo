# Pairwise Distances and Dendrograms from VCF files

First we are going to pass a function to unzip our files, creating a folder for this purpose, the first function will parse the information about called variants in one VCF file, converting that table into a Pandas DataFrame, ignoring the ##INFO lines and headers, once we check this function works we'll parse all the 26 samples contained in our folder, extracting only the essential info (position, reference nucleotide and variant, with option to include INDELs), concatenating them in a single DataFrame. This DataFrame may be converted to a presence/absence table -very useful and common in ecology studies and metagenomics- thanks to the crosstab Pandas method. This table tells us which samples share the called variants in all VCFs, hence allowing to calculate the distances among samples. Finally, we may use these pairwise distances (and scipy has a great library for this kinds of calculations) to elaborate the dendrograms and heatmaps that visually explicit the differences or the phylogenetic relationships among our samples.

In [1]:

```python
import pandas as pd
import numpy as np
import os
import re
import zipfile
```

#### Iteration 1 | Create a function to read a single VCF 

In [2]:

```python
#A general function that creates a folder and unzips the content of <input_path> in that folder, 
#returning a list of files with their corresponding path

def unzip(input_path):
    output_path=input_path.rstrip(".zip")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print("Folder" , output_path ,  "created")
    else:    
        print("Folder" , output_path ,  "already exists") 
    zip_ref = zipfile.ZipFile(input_path,'r')
    zip_ref.extractall(output_path)
    zip_ref.close()
    file_list = [output_path +'\\' + x for x in os.listdir(output_path)]
    return file_list
```

In [3]:

```python
#Then, a function that takes one of these files and parses it, with options for including:
#possible sequencing artifacts (ambiguous), labelled as 0/1 when no heterozigosis is possible in M.tuberculosis
#complete record or brief reference to position, reference and alternative allele.
#INDELs
def parse_one(file, include_ambiguous=True, completo=True, include_INDEL=True):    
                                                     #include_ambiguos allows us to avoid those probably flawed alleles marked as 0/1
    tabla=[]  
    encabezado=[]
    
    with open(file) as file:
        f=file.readlines()                           #Read the file line by line
   
        for line in f:                               
            
            if re.match("^#CHROM", line):            #if line.startswith #CHROM it's the header
                encabezado=line.split("\t")
                muestra=encabezado[-1].rstrip("\n")  #the last column of these headers is the name of the sample


            if re.match("[^#]", line):                #every line that not startswith # is a variant call
                linea=line.split("\t")                #splitting the variant call into a list (only some info will be useful)
                alt=linea[4]                          #5th column describes alternative allele(s)
                alt=alt.split(',')                    #splits alternative alleles (if any)
                if incluir_ambiguos==False:           #if include_ambiguous is un-selected
                    if re.match("[^0]", linea[-1]):   #make sure no 0/1 genotype is indicated
                        if (incluir_INDEL==True):     #if include INDEL is selected
                            tabla.append(linea)       #append the whole line
                        else:                         #if include INDEL is not selected
                            for i in range(len(alt)): #for every alt allele 
                                if (len(alt[i])==1) and (len(linea[3])==1): #make sure alt and ref are 1 bp long (if not, we have an INDEL)
                                    tabla.append(linea) #append the whole line
                else:                                  #if include ambiguous is allowed, no 0/1 checking is necessary
                    if (incluir_INDEL==True):          #if include INDEL is selected
                            tabla.append(linea)        #append the whole line
                    else:                              #if include INDEL is not selected
                        for i in range(len(alt)):      #for every alt allele 
                            if (len(alt[i])==1) and (len(linea[3])==1): #make sure alt and ref are 1 bp long (if not, we have an INDEL)
                                tabla.append(linea)    #append the whole line
          

        df= pd.DataFrame(tabla, columns=encabezado)    #make the DataFrame
        df.index=[muestra]*len(df)                     #sample name will be useful when building the whole dataset
        
        print("Código de muestra:", muestra)           #this line is disposable
    if completo==True:
        return df                                 #show the whole DataFrame (by default)
    
    else:
        return df[["POS","REF","ALT"]]            #for most of the applications, only position and reference->variant is needed
            
        
```

In [4]:

```python
#let's check it out:
parse_one(r'C:\\Users\\jorge\\Downloads\\pedroscampoy iisgm-bioinfo-test main data\\ALM23685B3COL6.combined.hf.SNP.final.vcf', incluir_ambiguos=True, completo=True).head()
Código de muestra: ALM23685B3COL6
```

Out[4]:

|                |  #CHROM |   POS |   ID |  REF |  ALT |      QUAL | FILTER |                                              INFO |         FORMAT |            ALM23685B3COL6\n |
| -------------: | ------: | ----: | ---: | ---: | ---: | --------: | -----: | ------------------------------------------------: | -------------: | --------------------------: |
| ALM23685B3COL6 | MTB_anc |  1701 |    . |    T |    C | 140323.33 |   PASS | AC=2;AF=1.00;AN=2;BaseQRankSum=0.814;DP=95;Exc... | GT:AD:DP:GQ:PL | 1/1:0,95:95:99:3679,283,0\n |
| ALM23685B3COL6 | MTB_anc |  2532 |    . |    C |    T |  47660.33 |   PASS | AC=2;AF=1.00;AN=2;DP=38;ExcessHet=3.0103;FS=0.... | GT:AD:DP:GQ:PL | 1/1:0,38:38:99:1460,114,0\n |
| ALM23685B3COL6 | MTB_anc |  8040 |    . |    G |    A |  88799.33 |   PASS | AC=2;AF=1.00;AN=2;DP=73;ExcessHet=3.0103;FS=0.... | GT:AD:DP:GQ:PL | 1/1:0,73:73:99:2791,220,0\n |
| ALM23685B3COL6 | MTB_anc |  9143 |    . |    C |    T |  38105.33 |   PASS | AC=2;AF=1.00;AN=2;DP=18;ExcessHet=3.0103;FS=0.... | GT:AD:DP:GQ:PL |   1/1:0,18:18:54:675,54,0\n |
| ALM23685B3COL6 | MTB_anc | 13460 |    . |    G |    A |  33355.32 |   PASS | AC=2;AF=1.00;AN=2;DP=22;ExcessHet=3.0103;FS=0.... | GT:AD:DP:GQ:PL |   1/1:0,22:22:66:851,66,0\n |

#### **Iteration 2 | Combine present SNP into a presence matrix**

In [5]:

```python
#this simpler version of above function outputs a Pandas Series with Position(Reference->Alternative)
#in sake of simplicity, I only left include_INDELs option

def lista_muts(file, incluir_INDEL=True):  
    
    tabla=[]                                           
    
    with open(file) as file:
        f=file.readlines()
   
        for line in f:
            if re.match("^#CHROM", line):             
                encabezado=line.split("\t")           
                muestra=encabezado[-1].rstrip("\n")   #take the sample's name, and keep going the same way as before                
            if re.match("[^#]", line):                
                linea=line.split("\t")                
                alt=linea[4]                          
                alt=alt.split(',')                    
                if re.match("[^0]", linea[-1]):       
                        if (incluir_INDEL==True):      
                            for i in range(len(alt)):
                                posic=linea[1] + '(' + linea[3] + '>' + alt[i] + ')' #takes position, reference and variant
                                tabla.append(posic)
                        else:
                            for i in range(len(alt)):
                                if (len(alt[i])==1) and (len(linea[3])==1):              #again, taking only 1bp long vars/refs
                                    posic=linea[1] + '(' + linea[3] + '>' + alt[i] + ')' #if include_INDEL is not selected
                                    tabla.append(posic)
                                
                                
            
    df= pd.DataFrame(tabla, columns=["Posición (Ref>Alt)"])  
    df.index=[muestra]*len(df)                   
    
    return df
```

In [6]:

```python
#let's check it out:
lista_alt=lista_muts(r'C:\\Users\\jorge\\Downloads\\pedroscampoy iisgm-bioinfo-test main data\\ALM23685B3COL6.combined.hf.SNP.final.vcf', incluir_INDEL=True)
lista_alt.head(20)
```

Out[6]:

|                | Posición (Ref>Alt) |
| -------------: | -----------------: |
| ALM23685B3COL6 |          1701(T>C) |
| ALM23685B3COL6 |          2532(C>T) |
| ALM23685B3COL6 |          8040(G>A) |
| ALM23685B3COL6 |          9143(C>T) |
| ALM23685B3COL6 |         13460(G>A) |
| ALM23685B3COL6 |         14251(G>A) |
| ALM23685B3COL6 |         14401(G>A) |
| ALM23685B3COL6 |         15117(G>C) |
| ALM23685B3COL6 |         16055(C>A) |
| ALM23685B3COL6 |         17608(G>C) |
| ALM23685B3COL6 |         23174(A>C) |
| ALM23685B3COL6 |         27463(C>G) |
| ALM23685B3COL6 |         27469(G>A) |
| ALM23685B3COL6 |         30496(A>G) |
| ALM23685B3COL6 |         33457(C>T) |
| ALM23685B3COL6 |         33551(T>G) |
| ALM23685B3COL6 |         36538(C>T) |
| ALM23685B3COL6 |         37763(G>C) |
| ALM23685B3COL6 |         39758(G>A) |
| ALM23685B3COL6 |         42281(A>C) |

#### Iteration 3 | Parse a collection of VCF files in .zip and create a presence/absence table of  SNPs/INDELs

In [7]:

```python
def sacar_multidf(path_in):       #this function unzips the files in <path_to_zipfile> and parses all its files
                                  #creating only one, long list of sample | position (ref->alt)
    
    diccionario={}
    file_list=unzip(path_in)
    
    for file in file_list:
        
        diccionario[file]=lista_muts(file, incluir_INDEL=True) 
        #including INDELs or not expands/contracts the list in 48rows (less than 2 INDELs per sample, on average)
    
    multi_df=pd.concat(list(diccionario.values()))
    pd.set_option('display.max_rows', 50)          #Here we may ask Pandas to show more rows tan usual
    return(multi_df)
        
    
```

In [8]:

```python
#Checking the function
sacar_multidf(r'C:\Users\jorge\Downloads\pedroscampoy iisgm-bioinfo-test main data.zip') 
Dir  C:\Users\jorge\Downloads\pedroscampoy iisgm-bioinfo-test main data  ya existe
```

Out[8]:

|                  | Posición (Ref>Alt) |
| ---------------: | -----------------: |
| 10082989-0-COL14 |          1701(T>C) |
| 10082989-0-COL14 |          2532(C>T) |
| 10082989-0-COL14 |          8040(G>A) |
| 10082989-0-COL14 |          9143(C>T) |
| 10082989-0-COL14 |         13460(G>A) |
|              ... |                ... |
|   ALM93896B3COL6 |       4397324(G>C) |
|   ALM93896B3COL6 |       4407588(C>T) |
|   ALM93896B3COL6 |       4408156(A>C) |
|   ALM93896B3COL6 |       4408920(G>A) |
|   ALM93896B3COL6 |       4408923(T>C) |

21306 rows × 1 columns

In [9]:

```python
#With this function we may take the whole DataFrame created above and make a presence/absence table,
#very useful for further distance calculation
def tabla_io(input_path):                                     
    df=sacar_multidf(input_path)
    tabla_final=pd.crosstab(df["Posición (Ref>Alt)"], df.index)
    pd.set_option('display.max_rows', 1000)
    pd.set_option('display.max_columns', 26)
    return tabla_final
```

In [10]:

```python
io=tabla_io(r'C:\Users\jorge\Downloads\pedroscampoy iisgm-bioinfo-test main data.zip') 

#io will be the name of the presence/absence table with which we'll make the pairwise distance calculations
Dir  C:\Users\jorge\Downloads\pedroscampoy iisgm-bioinfo-test main data  ya existe
```

In [11]:

```py
#checking the function works as supposed to
tabla_io(r'C:\Users\jorge\Downloads\pedroscampoy iisgm-bioinfo-test main data.zip').head()
Dir  C:\Users\jorge\Downloads\pedroscampoy iisgm-bioinfo-test main data  ya existe
```

Out[11]:

|              col_0 | 10082989-0-COL14 | 10082989-0-COL2 | 10082989-0-COL3 | 10082989-0-COL7 | 10105494-0-COL1 | 10105494-0-COL2 | 10105494-0-COL40 | AL10082989COL3 | AL10105494COL0 | ALM10082989B3CO11R18 | ALM10082989B3COL2018 | ALM10105494B3COL3018 | ALM23685B3COL2 | ALM23685B3COL5 | ALM23685B3COL6 | ALM23685B3COL7 | ALM23685B3COL8 | ALM23685B3COL9 | ALM93896B2COL31 | ALM93896B3COL1 | ALM93896B3COL11 | ALM93896B3COL12 | ALM93896B3COL2 | ALM93896B3COL3 | ALM93896B3COL5 | ALM93896B3COL6 |
| -----------------: | ---------------: | --------------: | --------------: | --------------: | --------------: | --------------: | ---------------: | -------------: | -------------: | -------------------: | -------------------: | -------------------: | -------------: | -------------: | -------------: | -------------: | -------------: | -------------: | --------------: | -------------: | --------------: | --------------: | -------------: | -------------: | -------------: | -------------: |
| Posición (Ref>Alt) |                  |                 |                 |                 |                 |                 |                  |                |                |                      |                      |                      |                |                |                |                |                |                |                 |                |                 |                 |                |                |                |                |
|       1012322(G>A) |                1 |               1 |               1 |               1 |               1 |               1 |                1 |              1 |              1 |                    1 |                    1 |                    1 |              1 |              1 |              1 |              1 |              1 |              1 |               1 |              1 |               1 |               1 |              1 |              1 |              1 |              1 |
|       1014815(G>T) |                1 |               1 |               1 |               1 |               1 |               1 |                1 |              1 |              1 |                    1 |                    1 |                    1 |              1 |              1 |              1 |              1 |              1 |              1 |               1 |              1 |               1 |               1 |              1 |              1 |              1 |              1 |
|       1024346(G>A) |                1 |               1 |               1 |               1 |               1 |               1 |                1 |              1 |              1 |                    1 |                    1 |                    1 |              1 |              1 |              1 |              1 |              1 |              1 |               1 |              1 |               1 |               1 |              1 |              1 |              1 |              1 |
|        103600(G>A) |                1 |               1 |               1 |               1 |               1 |               1 |                1 |              1 |              1 |                    1 |                    1 |                    1 |              0 |              1 |              1 |              0 |              1 |              1 |               1 |              1 |               1 |               1 |              1 |              1 |              1 |              1 |
|       1040050(C>T) |                1 |               1 |               1 |               1 |               1 |               1 |                1 |              1 |              1 |                    1 |                    1 |                    1 |              1 |              1 |              1 |              1 |              1 |              1 |               1 |              1 |               1 |               1 |              1 |              1 |              1 |              1 |

#### Iteration 4 | Calculate the SNP distance between all samples

In [12]:

```python
#There's a handful of methods to calculate pairwise distances among samples, the simplest is to compare row by row
#two samples and assign 0 if both row values are the same or 1 if different. Simple but computationally economic as stated in
#Jandrasits, 2019 https://doi.org/10.1371/journal.pcbi.1007527

import numpy as np

def distancias_matches(tabla_io):
    n=len(tabla_io.columns)
    tabla_distancias=np.zeros((n,n))
    for i in range(len(tabla_io.columns)):
        for j in range(len(tabla_io.columns)):
            resta=(tabla_io.iloc[:,j]-tabla_io.iloc[:,i]).abs().sum()
            tabla_distancias[i,j]=resta
            resta=0
    return(tabla_distancias)
    
dm=pd.DataFrame(distancias_matches(io), index=io.columns, columns=io.columns)   
dm   #We may observe at first sight that sample #13 ALM23685B3COL2 is unusually deviated from the rest
```

Out[12]:

|                      |     0 |     1 |     2 |     3 |     4 |     5 |     6 |     7 |     8 |     9 |    10 |    11 |    12 |    13 |    14 |    15 |    16 |    17 |    18 |    19 |    20 |    21 |    22 |    23 |    24 |    25 |
| -------------------: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: | ----: |
|                col_0 |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |       |
|     10082989-0-COL14 |   0.0 |   8.0 |  30.0 |  14.0 |  15.0 |  10.0 |  39.0 |  10.0 |  20.0 |  73.0 |  19.0 |   8.0 | 201.0 |  14.0 |  46.0 |  33.0 |  90.0 |  73.0 |  16.0 |  20.0 |  15.0 |  46.0 |  14.0 |  14.0 |  18.0 |  16.0 |
|      10082989-0-COL2 |   8.0 |   0.0 |  30.0 |   6.0 |   9.0 |   4.0 |  43.0 |   8.0 |  18.0 |  75.0 |  23.0 |   8.0 | 207.0 |  12.0 |  50.0 |  37.0 |  96.0 |  75.0 |  18.0 |  20.0 |  15.0 |  46.0 |  14.0 |  12.0 |  16.0 |  16.0 |
|      10082989-0-COL3 |  30.0 |  30.0 |   0.0 |  30.0 |  31.0 |  30.0 |  45.0 |  32.0 |  32.0 |  79.0 |  41.0 |  32.0 | 193.0 |  34.0 |  58.0 |  51.0 |  92.0 |  81.0 |  30.0 |  36.0 |  31.0 |  52.0 |  30.0 |  26.0 |  34.0 |  34.0 |
|      10082989-0-COL7 |  14.0 |   6.0 |  30.0 |   0.0 |   7.0 |   4.0 |  45.0 |   8.0 |  16.0 |  77.0 |  25.0 |  10.0 | 211.0 |  14.0 |  52.0 |  41.0 | 102.0 |  75.0 |  20.0 |  20.0 |  15.0 |  46.0 |  16.0 |  12.0 |  18.0 |  16.0 |
|      10105494-0-COL1 |  15.0 |   9.0 |  31.0 |   7.0 |   0.0 |   5.0 |  44.0 |   9.0 |  19.0 |  74.0 |  28.0 |   7.0 | 210.0 |  15.0 |  53.0 |  38.0 | 101.0 |  76.0 |  19.0 |  23.0 |  12.0 |  49.0 |  19.0 |  11.0 |  19.0 |  19.0 |
|      10105494-0-COL2 |  10.0 |   4.0 |  30.0 |   4.0 |   5.0 |   0.0 |  43.0 |   6.0 |  16.0 |  77.0 |  23.0 |   8.0 | 211.0 |  10.0 |  52.0 |  37.0 | 100.0 |  77.0 |  18.0 |  20.0 |  11.0 |  44.0 |  14.0 |  10.0 |  14.0 |  14.0 |
|     10105494-0-COL40 |  39.0 |  43.0 |  45.0 |  45.0 |  44.0 |  43.0 |   0.0 |  43.0 |  43.0 |  74.0 |  42.0 |  43.0 | 190.0 |  43.0 |  55.0 |  54.0 | 103.0 |  78.0 |  33.0 |  39.0 |  36.0 |  59.0 |  35.0 |  37.0 |  37.0 |  37.0 |
|       AL10082989COL3 |  10.0 |   8.0 |  32.0 |   8.0 |   9.0 |   6.0 |  43.0 |   0.0 |  18.0 |  77.0 |  21.0 |   6.0 | 209.0 |  10.0 |  52.0 |  37.0 | 100.0 |  75.0 |  18.0 |  18.0 |  13.0 |  44.0 |  14.0 |  12.0 |  14.0 |  14.0 |
|       AL10105494COL0 |  20.0 |  18.0 |  32.0 |  16.0 |  19.0 |  16.0 |  43.0 |  18.0 |   0.0 |  83.0 |  33.0 |  18.0 | 205.0 |  24.0 |  60.0 |  47.0 |  96.0 |  79.0 |  24.0 |  24.0 |  23.0 |  48.0 |  14.0 |  20.0 |  24.0 |  26.0 |
| ALM10082989B3CO11R18 |  73.0 |  75.0 |  79.0 |  77.0 |  74.0 |  77.0 |  74.0 |  77.0 |  83.0 |   0.0 |  70.0 |  71.0 | 180.0 |  77.0 |  81.0 |  80.0 | 117.0 |  94.0 |  65.0 |  75.0 |  70.0 |  91.0 |  73.0 |  73.0 |  75.0 |  69.0 |
| ALM10082989B3COL2018 |  19.0 |  23.0 |  41.0 |  25.0 |  28.0 |  23.0 |  42.0 |  21.0 |  33.0 |  70.0 |   0.0 |  23.0 | 194.0 |  25.0 |  47.0 |  40.0 |  91.0 |  68.0 |  23.0 |  29.0 |  22.0 |  45.0 |  25.0 |  23.0 |  21.0 |  17.0 |
| ALM10105494B3COL3018 |   8.0 |   8.0 |  32.0 |  10.0 |   7.0 |   8.0 |  43.0 |   6.0 |  18.0 |  71.0 |  23.0 |   0.0 | 207.0 |  12.0 |  50.0 |  35.0 |  98.0 |  73.0 |  18.0 |  18.0 |  15.0 |  44.0 |  14.0 |  14.0 |  16.0 |  16.0 |
|       ALM23685B3COL2 | 201.0 | 207.0 | 193.0 | 211.0 | 210.0 | 211.0 | 190.0 | 209.0 | 205.0 | 180.0 | 194.0 | 207.0 |   0.0 | 203.0 | 187.0 | 192.0 | 161.0 | 198.0 | 195.0 | 199.0 | 200.0 | 191.0 | 197.0 | 201.0 | 199.0 | 201.0 |
|       ALM23685B3COL5 |  14.0 |  12.0 |  34.0 |  14.0 |  15.0 |  10.0 |  43.0 |  10.0 |  24.0 |  77.0 |  25.0 |  12.0 | 203.0 |   0.0 |  48.0 |  33.0 |  92.0 |  75.0 |  22.0 |  24.0 |  17.0 |  44.0 |  18.0 |  14.0 |  16.0 |  16.0 |
|       ALM23685B3COL6 |  46.0 |  50.0 |  58.0 |  52.0 |  53.0 |  52.0 |  55.0 |  52.0 |  60.0 |  81.0 |  47.0 |  50.0 | 187.0 |  48.0 |   0.0 |  53.0 |  92.0 |  79.0 |  46.0 |  52.0 |  45.0 |  62.0 |  50.0 |  44.0 |  48.0 |  50.0 |
|       ALM23685B3COL7 |  33.0 |  37.0 |  51.0 |  41.0 |  38.0 |  37.0 |  54.0 |  37.0 |  47.0 |  80.0 |  40.0 |  35.0 | 192.0 |  33.0 |  53.0 |   0.0 |  83.0 |  74.0 |  33.0 |  41.0 |  30.0 |  57.0 |  37.0 |  39.0 |  37.0 |  35.0 |
|       ALM23685B3COL8 |  90.0 |  96.0 |  92.0 | 102.0 | 101.0 | 100.0 | 103.0 | 100.0 |  96.0 | 117.0 |  91.0 |  98.0 | 161.0 |  92.0 |  92.0 |  83.0 |   0.0 | 105.0 |  90.0 |  92.0 |  93.0 |  98.0 |  88.0 |  92.0 |  92.0 |  96.0 |
|       ALM23685B3COL9 |  73.0 |  75.0 |  81.0 |  75.0 |  76.0 |  77.0 |  78.0 |  75.0 |  79.0 |  94.0 |  68.0 |  73.0 | 198.0 |  75.0 |  79.0 |  74.0 | 105.0 |   0.0 |  65.0 |  65.0 |  66.0 |  87.0 |  71.0 |  71.0 |  67.0 |  69.0 |
|      ALM93896B2COL31 |  16.0 |  18.0 |  30.0 |  20.0 |  19.0 |  18.0 |  33.0 |  18.0 |  24.0 |  65.0 |  23.0 |  18.0 | 195.0 |  22.0 |  46.0 |  33.0 |  90.0 |  65.0 |   0.0 |  18.0 |  15.0 |  48.0 |  14.0 |  14.0 |  20.0 |  18.0 |
|       ALM93896B3COL1 |  20.0 |  20.0 |  36.0 |  20.0 |  23.0 |  20.0 |  39.0 |  18.0 |  24.0 |  75.0 |  29.0 |  18.0 | 199.0 |  24.0 |  52.0 |  41.0 |  92.0 |  65.0 |  18.0 |   0.0 |  21.0 |  44.0 |  18.0 |  22.0 |  22.0 |  18.0 |
|      ALM93896B3COL11 |  15.0 |  15.0 |  31.0 |  15.0 |  12.0 |  11.0 |  36.0 |  13.0 |  23.0 |  70.0 |  22.0 |  15.0 | 200.0 |  17.0 |  45.0 |  30.0 |  93.0 |  66.0 |  15.0 |  21.0 |   0.0 |  47.0 |  15.0 |  13.0 |  15.0 |  11.0 |
|      ALM93896B3COL12 |  46.0 |  46.0 |  52.0 |  46.0 |  49.0 |  44.0 |  59.0 |  44.0 |  48.0 |  91.0 |  45.0 |  44.0 | 191.0 |  44.0 |  62.0 |  57.0 |  98.0 |  87.0 |  48.0 |  44.0 |  47.0 |   0.0 |  38.0 |  46.0 |  42.0 |  42.0 |
|       ALM93896B3COL2 |  14.0 |  14.0 |  30.0 |  16.0 |  19.0 |  14.0 |  35.0 |  14.0 |  14.0 |  73.0 |  25.0 |  14.0 | 197.0 |  18.0 |  50.0 |  37.0 |  88.0 |  71.0 |  14.0 |  18.0 |  15.0 |  38.0 |   0.0 |  18.0 |  16.0 |  16.0 |
|       ALM93896B3COL3 |  14.0 |  12.0 |  26.0 |  12.0 |  11.0 |  10.0 |  37.0 |  12.0 |  20.0 |  73.0 |  23.0 |  14.0 | 201.0 |  14.0 |  44.0 |  39.0 |  92.0 |  71.0 |  14.0 |  22.0 |  13.0 |  46.0 |  18.0 |   0.0 |  16.0 |  20.0 |
|       ALM93896B3COL5 |  18.0 |  16.0 |  34.0 |  18.0 |  19.0 |  14.0 |  37.0 |  14.0 |  24.0 |  75.0 |  21.0 |  16.0 | 199.0 |  16.0 |  48.0 |  37.0 |  92.0 |  67.0 |  20.0 |  22.0 |  15.0 |  42.0 |  16.0 |  16.0 |   0.0 |  16.0 |
|       ALM93896B3COL6 |  16.0 |  16.0 |  34.0 |  16.0 |  19.0 |  14.0 |  37.0 |  14.0 |  26.0 |  69.0 |  17.0 |  16.0 | 201.0 |  16.0 |  50.0 |  35.0 |  96.0 |  69.0 |  18.0 |  18.0 |  11.0 |  42.0 |  16.0 |  20.0 |  16.0 |   0.0 |

In [13]:

```python
#Jaccard Index is another, more common way to measure pairwise distances between groups, specially regarding those situations 
#in which absence is way more common than presence (e.g., when looking for species present or absent in some ecosystems)

from sklearn.metrics import jaccard_score

def distancias_Jaccard(tabla_io):
    n=len(tabla_io.columns)
    tabla_Jaccard=np.zeros((n,n))
    for i in range(len(tabla_io.columns)):
        for j in range(len(tabla_io.columns)):
            resta_jac=jaccard_score(tabla_io.iloc[:,j], tabla_io.iloc[:,i], average='macro')
            tabla_Jaccard[i,j]=resta_jac
    return tabla_Jaccard
    

jac=pd.DataFrame(data=distancias_Jaccard(io), index=io.columns, columns=io.columns)
jac
```

Out[13]:

|                col_0 | 10082989-0-COL14 | 10082989-0-COL2 | 10082989-0-COL3 | 10082989-0-COL7 | 10105494-0-COL1 | 10105494-0-COL2 | 10105494-0-COL40 | AL10082989COL3 | AL10105494COL0 | ALM10082989B3CO11R18 | ALM10082989B3COL2018 | ALM10105494B3COL3018 | ALM23685B3COL2 | ALM23685B3COL5 | ALM23685B3COL6 | ALM23685B3COL7 | ALM23685B3COL8 | ALM23685B3COL9 | ALM93896B2COL31 | ALM93896B3COL1 | ALM93896B3COL11 | ALM93896B3COL12 | ALM93896B3COL2 | ALM93896B3COL3 | ALM93896B3COL5 | ALM93896B3COL6 |
| -------------------: | ---------------: | --------------: | --------------: | --------------: | --------------: | --------------: | ---------------: | -------------: | -------------: | -------------------: | -------------------: | -------------------: | -------------: | -------------: | -------------: | -------------: | -------------: | -------------: | --------------: | -------------: | --------------: | --------------: | -------------: | -------------: | -------------: | -------------: |
|                col_0 |                  |                 |                 |                 |                 |                 |                  |                |                |                      |                      |                      |                |                |                |                |                |                |                 |                |                 |                 |                |                |                |                |
|     10082989-0-COL14 |         1.000000 |        0.813465 |        0.656183 |        0.700117 |        0.702725 |        0.756029 |         0.622405 |       0.766845 |       0.665668 |             0.555655 |             0.724869 |             0.828605 |       0.425890 |       0.722534 |       0.607702 |       0.643762 |       0.533775 |       0.555655 |        0.740544 |       0.710402 |        0.732535 |        0.595796 |       0.758392 |       0.741745 |       0.708137 |       0.732490 |
|      10082989-0-COL2 |         0.813465 |        1.000000 |        0.633495 |        0.820004 |        0.769712 |        0.864316 |         0.576587 |       0.773083 |       0.656116 |             0.534430 |             0.657899 |             0.804812 |       0.409622 |       0.720230 |       0.567362 |       0.592844 |       0.503096 |       0.534430 |        0.689412 |       0.685191 |        0.702725 |        0.576389 |       0.732496 |       0.742941 |       0.704874 |       0.704874 |
|      10082989-0-COL3 |         0.656183 |        0.633495 |        1.000000 |        0.625210 |        0.629471 |        0.616520 |         0.642300 |       0.609106 |       0.640684 |             0.573000 |             0.616119 |             0.633306 |       0.458344 |       0.610435 |       0.598429 |       0.594822 |       0.561511 |       0.565956 |        0.676084 |       0.645340 |        0.651891 |        0.617696 |       0.669727 |       0.689143 |       0.639905 |       0.639905 |
|      10082989-0-COL7 |         0.700117 |        0.820004 |        0.625210 |        1.000000 |        0.801443 |        0.854795 |         0.556956 |       0.760017 |       0.670621 |             0.522122 |             0.628186 |             0.756029 |       0.400803 |       0.673621 |       0.550129 |       0.557628 |       0.480752 |       0.529694 |        0.654930 |       0.675735 |        0.691197 |        0.569464 |       0.694303 |       0.732080 |       0.668008 |       0.694303 |
|      10105494-0-COL1 |         0.702725 |        0.769712 |        0.629471 |        0.801443 |        1.000000 |        0.840812 |         0.574148 |       0.757876 |       0.649577 |             0.540586 |             0.615147 |             0.829206 |       0.405473 |       0.678708 |       0.554798 |       0.589918 |       0.489765 |       0.532966 |        0.682372 |       0.657899 |        0.752933 |        0.562877 |       0.672170 |       0.764355 |       0.672170 |       0.672170 |
|      10105494-0-COL2 |         0.756029 |        0.864316 |        0.616520 |        0.854795 |        0.840812 |        1.000000 |         0.561333 |       0.796483 |       0.657299 |             0.517259 |             0.638017 |             0.784773 |       0.398798 |       0.730980 |       0.543290 |       0.576138 |       0.482668 |       0.517259 |        0.668008 |       0.665668 |        0.743537 |        0.574148 |       0.711774 |       0.756029 |       0.711774 |       0.711774 |
|     10105494-0-COL40 |         0.622405 |        0.576587 |        0.642300 |        0.556956 |        0.574148 |        0.561333 |         1.000000 |       0.569105 |       0.597453 |             0.606365 |             0.636350 |             0.590748 |       0.471385 |       0.583797 |       0.631700 |       0.603069 |       0.542170 |       0.592014 |        0.680357 |       0.651813 |        0.645340 |        0.605083 |       0.661034 |       0.635514 |       0.647697 |       0.647697 |
|       AL10082989COL3 |         0.766845 |        0.773083 |        0.609106 |        0.760017 |        0.757876 |        0.796483 |         0.569105 |       1.000000 |       0.643308 |             0.522122 |             0.669465 |             0.838572 |       0.404208 |       0.744131 |       0.550129 |       0.584669 |       0.486769 |       0.529694 |        0.679079 |       0.699077 |        0.721529 |        0.581261 |       0.722534 |       0.732080 |       0.722534 |       0.722534 |
|       AL10105494COL0 |         0.665668 |        0.656116 |        0.640684 |        0.670621 |        0.649577 |        0.657299 |         0.597453 |       0.643308 |       1.000000 |             0.518827 |             0.596913 |             0.679079 |       0.418913 |       0.598868 |       0.536217 |       0.552807 |       0.514691 |       0.533152 |        0.652549 |       0.670060 |        0.638017 |        0.584601 |       0.758392 |       0.665668 |       0.643042 |       0.623631 |
| ALM10082989B3CO11R18 |         0.555655 |        0.534430 |        0.573000 |        0.522122 |        0.540586 |        0.517259 |         0.606365 |       0.522122 |       0.518827 |             1.000000 |             0.589763 |             0.559062 |       0.518366 |       0.531522 |       0.595883 |       0.574795 |       0.553919 |       0.594605 |        0.600060 |       0.568759 |        0.569543 |        0.556685 |       0.564178 |       0.555655 |       0.556580 |       0.579808 |
| ALM10082989B3COL2018 |         0.724869 |        0.657899 |        0.616119 |        0.628186 |        0.615147 |        0.638017 |         0.636350 |       0.669465 |       0.596913 |             0.589763 |             1.000000 |             0.675628 |       0.448576 |       0.647456 |       0.631444 |       0.631476 |       0.550052 |       0.597628 |        0.705887 |       0.667623 |        0.697524 |        0.632337 |       0.680347 |       0.683791 |       0.718328 |       0.760175 |
| ALM10105494B3COL3018 |         0.828605 |        0.804812 |        0.633306 |        0.756029 |        0.829206 |        0.784773 |         0.590748 |       0.838572 |       0.679079 |             0.559062 |             0.675628 |             1.000000 |       0.413534 |       0.742941 |       0.579894 |       0.622220 |       0.504717 |       0.551249 |        0.708137 |       0.724656 |        0.723309 |        0.601145 |       0.750356 |       0.732496 |       0.723899 |       0.723899 |
|       ALM23685B3COL2 |         0.425890 |        0.409622 |        0.458344 |        0.400803 |        0.405473 |        0.398798 |         0.471385 |       0.404208 |       0.418913 |             0.518366 |             0.448576 |             0.413534 |       1.000000 |       0.418518 |       0.484129 |       0.462725 |       0.566744 |       0.486019 |        0.442197 |       0.438759 |        0.428605 |        0.473665 |       0.436757 |       0.425890 |       0.433208 |       0.429684 |
|       ALM23685B3COL5 |         0.722534 |        0.720230 |        0.610435 |        0.673621 |        0.678708 |        0.730980 |         0.583797 |       0.744131 |       0.598868 |             0.531522 |             0.647456 |             0.742941 |       0.418518 |       1.000000 |       0.584601 |       0.629479 |       0.519701 |       0.539059 |        0.653741 |       0.652549 |        0.686440 |        0.594746 |       0.689412 |       0.722534 |       0.714715 |       0.714715 |
|       ALM23685B3COL6 |         0.607702 |        0.567362 |        0.598429 |        0.550129 |        0.554798 |        0.543290 |         0.631700 |       0.550129 |       0.536217 |             0.595883 |             0.631444 |             0.579894 |       0.484129 |       0.584601 |       1.000000 |       0.628596 |       0.587653 |       0.602737 |        0.624199 |       0.602961 |        0.616198 |        0.610646 |       0.597280 |       0.619095 |       0.607928 |       0.597280 |
|       ALM23685B3COL7 |         0.643762 |        0.592844 |        0.594822 |        0.557628 |        0.589918 |        0.576138 |         0.603069 |       0.584669 |       0.552807 |             0.574795 |             0.631476 |             0.622220 |       0.462725 |       0.629479 |       0.628596 |       1.000000 |       0.595178 |       0.596359 |        0.663119 |       0.622291 |        0.669727 |        0.600727 |       0.629076 |       0.602032 |       0.629076 |       0.642751 |
|       ALM23685B3COL8 |         0.533775 |        0.503096 |        0.561511 |        0.480752 |        0.489765 |        0.482668 |         0.542170 |       0.486769 |       0.514691 |             0.553919 |             0.550052 |             0.504717 |       0.566744 |       0.519701 |       0.587653 |       0.595178 |       1.000000 |       0.585741 |        0.544643 |       0.545173 |        0.525986 |        0.564041 |       0.547619 |       0.527316 |       0.534654 |       0.522075 |
|       ALM23685B3COL9 |         0.555655 |        0.534430 |        0.565956 |        0.529694 |        0.532966 |        0.517259 |         0.592014 |       0.529694 |       0.533152 |             0.594605 |             0.597628 |             0.551249 |       0.486019 |       0.539059 |       0.602737 |       0.596359 |       0.585741 |       1.000000 |        0.600060 |       0.607817 |        0.585714 |        0.569581 |       0.571919 |       0.563444 |       0.587849 |       0.579808 |
|      ALM93896B2COL31 |         0.740544 |        0.689412 |        0.676084 |        0.654930 |        0.682372 |        0.668008 |         0.680357 |       0.679079 |       0.652549 |             0.600060 |             0.705887 |             0.708137 |       0.442197 |       0.653741 |       0.624199 |       0.663119 |       0.544643 |       0.600060 |        1.000000 |       0.752469 |        0.756749 |        0.602367 |       0.779575 |       0.765910 |       0.710402 |       0.732206 |
|       ALM93896B3COL1 |         0.710402 |        0.685191 |        0.645340 |        0.675735 |        0.657899 |        0.665668 |         0.651813 |       0.699077 |       0.670060 |             0.568759 |             0.667623 |             0.724656 |       0.438759 |       0.652549 |       0.602961 |       0.622291 |       0.545173 |       0.607817 |        0.752469 |       1.000000 |        0.703805 |        0.635379 |       0.746081 |       0.689716 |       0.704931 |       0.746081 |
|      ALM93896B3COL11 |         0.732535 |        0.702725 |        0.651891 |        0.691197 |        0.752933 |        0.743537 |         0.645340 |       0.721529 |       0.638017 |             0.569543 |             0.697524 |             0.723309 |       0.428605 |       0.686440 |       0.616198 |       0.669727 |       0.525986 |       0.585714 |        0.756749 |       0.703805 |        1.000000 |        0.593223 |       0.749199 |       0.760183 |       0.749199 |       0.803828 |
|      ALM93896B3COL12 |         0.595796 |        0.576389 |        0.617696 |        0.569464 |        0.562877 |        0.574148 |         0.605083 |       0.581261 |       0.584601 |             0.556685 |             0.632337 |             0.601145 |       0.473665 |       0.594746 |       0.610646 |       0.600727 |       0.564041 |       0.569581 |        0.602367 |       0.635379 |        0.593223 |        1.000000 |       0.655374 |       0.595796 |       0.630827 |       0.630827 |
|       ALM93896B3COL2 |         0.758392 |        0.732496 |        0.669727 |        0.694303 |        0.672170 |        0.711774 |         0.661034 |       0.722534 |       0.758392 |             0.564178 |             0.680347 |             0.750356 |       0.436757 |       0.689412 |       0.597280 |       0.629076 |       0.547619 |       0.571919 |        0.779575 |       0.746081 |        0.749199 |        0.655374 |       1.000000 |       0.708137 |       0.748108 |       0.748108 |
|       ALM93896B3COL3 |         0.741745 |        0.742941 |        0.689143 |        0.732080 |        0.764355 |        0.756029 |         0.635514 |       0.732080 |       0.665668 |             0.555655 |             0.683791 |             0.732496 |       0.425890 |       0.722534 |       0.619095 |       0.602032 |       0.527316 |       0.563444 |        0.765910 |       0.689716 |        0.760183 |        0.595796 |       0.708137 |       1.000000 |       0.732490 |       0.685191 |
|       ALM93896B3COL5 |         0.708137 |        0.704874 |        0.639905 |        0.668008 |        0.672170 |        0.711774 |         0.647697 |       0.722534 |       0.643042 |             0.556580 |             0.718328 |             0.723899 |       0.433208 |       0.714715 |       0.607928 |       0.629076 |       0.534654 |       0.587849 |        0.710402 |       0.704931 |        0.749199 |        0.630827 |       0.748108 |       0.732490 |       1.000000 |       0.748108 |
|       ALM93896B3COL6 |         0.732490 |        0.704874 |        0.639905 |        0.694303 |        0.672170 |        0.711774 |         0.647697 |       0.722534 |       0.623631 |             0.579808 |             0.760175 |             0.723899 |       0.429684 |       0.714715 |       0.597280 |       0.642751 |       0.522075 |       0.579808 |        0.732206 |       0.746081 |        0.803828 |        0.630827 |       0.748108 |       0.685191 |       0.748108 |       1.000000 |

In [14]:

```python
#Pearson correlation coefficient is, maybe the most convenient way to measure pairwise distance from a Pandas DataFrame, 
#as it is a Pandas DataFrame method itself:

io.corr()

#One more time, the #13 sample has almost no correlation with the rest 
```

Out[14]:

|                col_0 | 10082989-0-COL14 | 10082989-0-COL2 | 10082989-0-COL3 | 10082989-0-COL7 | 10105494-0-COL1 | 10105494-0-COL2 | 10105494-0-COL40 | AL10082989COL3 | AL10105494COL0 | ALM10082989B3CO11R18 | ALM10082989B3COL2018 | ALM10105494B3COL3018 | ALM23685B3COL2 | ALM23685B3COL5 | ALM23685B3COL6 | ALM23685B3COL7 | ALM23685B3COL8 | ALM23685B3COL9 | ALM93896B2COL31 | ALM93896B3COL1 | ALM93896B3COL11 | ALM93896B3COL12 | ALM93896B3COL2 | ALM93896B3COL3 | ALM93896B3COL5 | ALM93896B3COL6 |
| -------------------: | ---------------: | --------------: | --------------: | --------------: | --------------: | --------------: | ---------------: | -------------: | -------------: | -------------------: | -------------------: | -------------------: | -------------: | -------------: | -------------: | -------------: | -------------: | -------------: | --------------: | -------------: | --------------: | --------------: | -------------: | -------------: | -------------: | -------------: |
|                col_0 |                  |                 |                 |                 |                 |                 |                  |                |                |                      |                      |                      |                |                |                |                |                |                |                 |                |                 |                 |                |                |                |                |
|     10082989-0-COL14 |         1.000000 |        0.784610 |        0.530348 |        0.597857 |        0.591532 |        0.719482 |         0.475847 |       0.721339 |       0.511919 |             0.394081 |             0.645576 |             0.796305 |       0.252681 |       0.627010 |       0.463766 |       0.510388 |       0.374110 |       0.394081 |        0.662770 |       0.616083 |        0.642427 |        0.420567 |       0.690132 |       0.658343 |       0.600452 |       0.645292 |
|      10082989-0-COL2 |         0.784610 |        1.000000 |        0.512137 |        0.784265 |        0.704747 |        0.854334 |         0.384471 |       0.711471 |       0.496882 |             0.365366 |             0.536972 |             0.765569 |       0.206983 |       0.619193 |       0.385544 |       0.412550 |       0.299684 |       0.365366 |        0.587260 |       0.593712 |        0.597280 |        0.400854 |       0.664276 |       0.669519 |       0.611409 |       0.611409 |
|      10082989-0-COL3 |         0.530348 |        0.512137 |        1.000000 |        0.509008 |        0.494247 |        0.508758 |         0.480848 |       0.464286 |       0.494994 |             0.374710 |             0.417385 |             0.486185 |       0.294237 |       0.438708 |       0.392719 |       0.369179 |       0.378718 |       0.356708 |        0.554343 |       0.483746 |        0.516819 |        0.432284 |       0.546012 |       0.601056 |       0.481053 |       0.481053 |
|      10082989-0-COL7 |         0.597857 |        0.784265 |        0.509008 |        1.000000 |        0.758746 |        0.833980 |         0.335783 |       0.687596 |       0.536116 |             0.335545 |             0.479166 |             0.694578 |       0.167720 |       0.530118 |       0.343421 |       0.317295 |       0.213228 |       0.366984 |        0.524206 |       0.589553 |        0.583549 |        0.396102 |       0.602553 |       0.659598 |       0.545831 |       0.602553 |
|      10105494-0-COL1 |         0.591532 |        0.704747 |        0.494247 |        0.758746 |        1.000000 |        0.826717 |         0.370349 |       0.688222 |       0.480030 |             0.379464 |             0.427617 |             0.798949 |       0.175783 |       0.536849 |       0.337153 |       0.397693 |       0.234969 |       0.351075 |        0.567402 |       0.527377 |        0.686268 |        0.350943 |       0.539639 |       0.703035 |       0.539639 |       0.539639 |
|      10105494-0-COL2 |         0.719482 |        0.854334 |        0.508758 |        0.833980 |        0.826717 |        1.000000 |         0.369735 |       0.749175 |       0.518360 |             0.337179 |             0.524355 |             0.757301 |       0.170511 |       0.652889 |       0.337570 |       0.396286 |       0.239313 |       0.337179 |        0.572920 |       0.588642 |        0.702522 |        0.435499 |       0.657846 |       0.719482 |       0.657846 |       0.657846 |
|     10105494-0-COL40 |         0.475847 |        0.384471 |        0.480848 |        0.335783 |        0.370349 |        0.369735 |         1.000000 |       0.376504 |       0.411464 |             0.440873 |             0.476277 |             0.402184 |       0.310067 |       0.393127 |       0.463407 |       0.393443 |       0.312973 |       0.408089 |        0.582197 |       0.511807 |        0.526297 |        0.402258 |       0.548596 |       0.508038 |       0.519022 |       0.519022 |
|       AL10082989COL3 |         0.721339 |        0.711471 |        0.464286 |        0.687596 |        0.688222 |        0.749175 |         0.376504 |       1.000000 |       0.474375 |             0.335545 |             0.579859 |             0.824242 |       0.189552 |       0.667036 |       0.343421 |       0.403794 |       0.242096 |       0.366984 |        0.578852 |       0.640675 |        0.643906 |        0.435048 |       0.659274 |       0.659598 |       0.659274 |       0.659274 |
|       AL10105494COL0 |         0.511919 |        0.496882 |        0.494994 |        0.536116 |        0.480030 |        0.518360 |         0.411464 |       0.474375 |       1.000000 |             0.269813 |             0.366975 |             0.540046 |       0.218163 |       0.356417 |       0.255167 |       0.271059 |       0.305648 |       0.319520 |        0.489974 |       0.535258 |        0.451570 |        0.389779 |       0.690132 |       0.511919 |       0.465932 |       0.421091 |
| ALM10082989B3CO11R18 |         0.394081 |        0.365366 |        0.374710 |        0.335545 |        0.379464 |        0.337179 |         0.440873 |       0.335545 |       0.269813 |             1.000000 |             0.440458 |             0.419140 |       0.365081 |       0.337981 |       0.409357 |       0.374451 |       0.326268 |       0.405186 |        0.489301 |       0.387602 |        0.431362 |        0.319589 |       0.398385 |       0.394081 |       0.375552 |       0.444051 |
| ALM10082989B3COL2018 |         0.645576 |        0.536972 |        0.417385 |        0.479166 |        0.427617 |        0.524355 |         0.476277 |       0.579859 |       0.366975 |             0.440458 |             1.000000 |             0.555623 |       0.293201 |       0.501717 |       0.481380 |       0.456298 |       0.372702 |       0.460725 |        0.598767 |       0.522299 |        0.590756 |        0.475963 |       0.551099 |       0.565976 |       0.624228 |       0.697357 |
| ALM10105494B3COL3018 |         0.796305 |        0.765569 |        0.486185 |        0.694578 |        0.798949 |        0.757301 |         0.402184 |       0.824242 |       0.540046 |             0.419140 |             0.555623 |             1.000000 |       0.202090 |       0.660644 |       0.397393 |       0.467804 |       0.278973 |       0.393042 |        0.608071 |       0.649981 |        0.627046 |        0.445749 |       0.680331 |       0.642550 |       0.633246 |       0.633246 |
|       ALM23685B3COL2 |         0.252681 |        0.206983 |        0.294237 |        0.167720 |        0.175783 |        0.170511 |         0.310067 |       0.189552 |       0.218163 |             0.365081 |             0.293201 |             0.202090 |       1.000000 |       0.242299 |       0.326057 |       0.299472 |       0.449419 |       0.285984 |        0.291906 |       0.258489 |        0.259585 |        0.304433 |       0.279376 |       0.252681 |       0.263520 |       0.247664 |
|       ALM23685B3COL5 |         0.627010 |        0.619193 |        0.438708 |        0.530118 |        0.536849 |        0.652889 |         0.393127 |       0.667036 |       0.356417 |             0.337981 |             0.501717 |             0.660644 |       0.242299 |       1.000000 |       0.424145 |       0.497841 |       0.351467 |       0.365539 |        0.501387 |       0.510322 |        0.559010 |        0.440901 |       0.572114 |       0.627010 |       0.621833 |       0.621833 |
|       ALM23685B3COL6 |         0.463766 |        0.385544 |        0.392719 |        0.343421 |        0.337153 |        0.337570 |         0.463407 |       0.343421 |       0.255167 |             0.409357 |             0.481380 |             0.397393 |       0.326057 |       0.424145 |       1.000000 |       0.459002 |       0.413493 |       0.424531 |        0.478765 |       0.416455 |        0.480511 |        0.417977 |       0.418568 |       0.493566 |       0.445945 |       0.418568 |
|       ALM23685B3COL7 |         0.510388 |        0.412550 |        0.369179 |        0.317295 |        0.397693 |        0.396286 |         0.393443 |       0.403794 |       0.271059 |             0.374451 |             0.456298 |             0.467804 |       0.299472 |       0.497841 |       0.459002 |       1.000000 |       0.456212 |       0.426680 |        0.533249 |       0.436456 |        0.564091 |        0.392356 |       0.462481 |       0.407819 |       0.462481 |       0.493892 |
|       ALM23685B3COL8 |         0.374110 |        0.299684 |        0.378718 |        0.213228 |        0.234969 |        0.239313 |         0.312973 |       0.242096 |       0.305648 |             0.326268 |             0.372702 |             0.278973 |       0.449419 |       0.351467 |       0.413493 |       0.456212 |       1.000000 |       0.395992 |        0.376670 |       0.361888 |        0.340648 |        0.362373 |       0.396183 |       0.351290 |       0.354252 |       0.312321 |
|       ALM23685B3COL9 |         0.394081 |        0.365366 |        0.356708 |        0.366984 |        0.351075 |        0.337179 |         0.408089 |       0.366984 |       0.319520 |             0.405186 |             0.460725 |             0.393042 |       0.285984 |       0.365539 |       0.424531 |       0.426680 |       0.395992 |       1.000000 |        0.489301 |       0.490494 |        0.479955 |        0.350944 |       0.421218 |       0.418935 |       0.466884 |       0.444051 |
|      ALM93896B2COL31 |         0.662770 |        0.587260 |        0.554343 |        0.524206 |        0.567402 |        0.572920 |         0.582197 |       0.578852 |       0.489974 |             0.489301 |             0.598767 |             0.608071 |       0.291906 |       0.501387 |       0.478765 |       0.533249 |       0.376670 |       0.489301 |        1.000000 |       0.680639 |        0.688826 |        0.416303 |       0.722977 |       0.705969 |       0.603916 |       0.643603 |
|       ALM93896B3COL1 |         0.616083 |        0.593712 |        0.483746 |        0.589553 |        0.527377 |        0.588642 |         0.511807 |       0.640675 |       0.535258 |             0.387602 |             0.522299 |             0.649981 |       0.258489 |       0.510322 |       0.416455 |       0.436456 |       0.361888 |       0.490494 |        0.680639 |       1.000000 |        0.600856 |        0.484914 |       0.672038 |       0.575671 |       0.597784 |       0.672038 |
|      ALM93896B3COL11 |         0.642427 |        0.597280 |        0.516819 |        0.583549 |        0.686268 |        0.702522 |         0.526297 |       0.643906 |       0.451570 |             0.431362 |             0.590756 |             0.627046 |       0.259585 |       0.559010 |       0.480511 |       0.564091 |       0.340648 |       0.479955 |        0.688826 |       0.600856 |        1.000000 |        0.409222 |       0.673395 |       0.690142 |       0.673395 |       0.761065 |
|      ALM93896B3COL12 |         0.420567 |        0.400854 |        0.432284 |        0.396102 |        0.350943 |        0.435499 |         0.402258 |       0.435048 |       0.389779 |             0.319589 |             0.475963 |             0.445749 |       0.304433 |       0.440901 |       0.417977 |       0.392356 |       0.362373 |       0.350944 |        0.416303 |       0.484914 |        0.409222 |        1.000000 |       0.548867 |       0.420567 |       0.492297 |       0.492297 |
|       ALM93896B3COL2 |         0.690132 |        0.664276 |        0.546012 |        0.602553 |        0.539639 |        0.657846 |         0.548596 |       0.659274 |       0.690132 |             0.398385 |             0.551099 |             0.680331 |       0.279376 |       0.572114 |       0.418568 |       0.462481 |       0.396183 |       0.421218 |        0.722977 |       0.672038 |        0.673395 |        0.548867 |       1.000000 |       0.600452 |       0.670442 |       0.670442 |
|       ALM93896B3COL3 |         0.658343 |        0.669519 |        0.601056 |        0.659598 |        0.703035 |        0.719482 |         0.508038 |       0.659598 |       0.511919 |             0.394081 |             0.565976 |             0.642550 |       0.252681 |       0.627010 |       0.493566 |       0.407819 |       0.351290 |       0.418935 |        0.705969 |       0.575671 |        0.690142 |        0.420567 |       0.600452 |       1.000000 |       0.645292 |       0.555612 |
|       ALM93896B3COL5 |         0.600452 |        0.611409 |        0.481053 |        0.545831 |        0.539639 |        0.657846 |         0.519022 |       0.659274 |       0.465932 |             0.375552 |             0.624228 |             0.633246 |       0.263520 |       0.621833 |       0.445945 |       0.462481 |       0.354252 |       0.466884 |        0.603916 |       0.597784 |        0.673395 |        0.492297 |       0.670442 |       0.645292 |       1.000000 |       0.670442 |
|       ALM93896B3COL6 |         0.645292 |        0.611409 |        0.481053 |        0.602553 |        0.539639 |        0.657846 |         0.519022 |       0.659274 |       0.421091 |             0.444051 |             0.697357 |             0.633246 |       0.247664 |       0.621833 |       0.418568 |       0.493892 |       0.312321 |       0.444051 |        0.643603 |       0.672038 |        0.761065 |        0.492297 |       0.670442 |       0.555612 |       0.670442 |       1.000000 |

In [15]:

```python
#scipy has lots of different methods to measure pairwise distance inside the pdist function, let's try Hamming distance:
#don't forget to use squareform to get the distance matrix 
#and to use the transposed matrix of presence/absence or we'll obtain the distances between pairs of variants

from scipy.spatial.distance import pdist, squareform

y=squareform(scipy.spatial.distance.pdist(io.T, metric='hamming')) 

y=pd.DataFrame(data=y, index=io.columns, columns=io.columns)
```



#### Iteration 5 | Represent distance in heatmaps and phylogenetic trees

In [16]:

```python
#Maybe the most useful way to discriminate different samples from those which are similar amongst them is by plotting.
#Let's use pyplot, the ggplot equivalent in Python
#Here, I'm using tril to get only the trinagular form of these heatmaps (we might use triufor the other triangle)

import matplotlib.pyplot as plt
import numpy as np
from numpy import tril
from scipy.spatial.distance import pdist, squareform

n_clusters=len(io.columns)


for index, metric in enumerate(['braycurtis', 'hamming', 'jaccard']):           
    
    #Bray-Curtis is more usual in metagenomics experiments, when abundances are specified
    
    avg_dist = np.zeros((n_clusters, n_clusters))                               
    plt.figure(figsize=(14, 14))
    avg_dist=squareform(scipy.spatial.distance.pdist(io.T, metric=metric))      #Do not forget the transposed DataFrame.
    avg_dist=tril(avg_dist)                                                     #disposable
    avg_dist /= avg_dist.max()
    plt.imshow(avg_dist, cmap=plt.cm.gnuplot2, vmin=0)
    plt.xticks(range(n_clusters), io.columns, rotation=90, fontsize=8)
    plt.yticks(range(n_clusters), io.columns)
    plt.colorbar()
    plt.suptitle("%s distances" % metric, size=24)
    plt.tight_layout()
```



![Hamming](https://user-images.githubusercontent.com/77884314/124262521-743dff80-db32-11eb-826b-c31b26ebe103.jpg)


![Jaccard](https://user-images.githubusercontent.com/77884314/124262489-69836a80-db32-11eb-9a84-93bfe0c844e2.jpg)


In [17]:

```python
#For the dendrogram and heatmaps I'm going to use the matrix distance from Jaccard index method 
#and the WPGMA algorithm to conduct the tree formation. Finally cophenet function will estimate
#correlation between the "actual distances" and the calculated tree.

from scipy.cluster.hierarchy import dendrogram, linkage, cophenet

Z=linkage(io.T, method='weighted', metric='jaccard') 

#print(Z[:30]) #In less than 30 iterations the linkage function has clustered all samples

c, coph_dists = cophenet(Z, pdist(io.T)) 
c

  

0.9857297779347147
```

In [18]:

```python
#Finally, we've got a dendrogram like this:

plt.figure(figsize=(25, 10))
plt.title('Final dendrogram', size =40)
plt.xlabel('Samples', size=18)
plt.ylabel('Distances', size=18)
dendrogram(Z, labels=io.columns, leaf_rotation=90, leaf_font_size=18)
plt.show()

#We may see that, indeed, sample ALM23685B3 COL2 was way too different from the rest while there are 16 other samples
#really close to each other, indicating recent contagion events. It is important to note that M.tuberculosis has a 
#very low mutation rate (0.3-0.5 mutations per year and genome)
```

![Final Dendrogram](https://user-images.githubusercontent.com/77884314/124262443-5a9cb800-db32-11eb-8658-364c14c8dab0.jpg)


In [19]:

```python
#Depending on the objectives of the experiment we might be interested in:
#A) finding out the recent story of an outbreak and, so, repeat the tree with the most closely related samples.
#B) focus only in the phylogeny of quite different strains and,then, prune the tree to the ten most distant branches:

#A)
y=io[["10082989-0-COL14","10082989-0-COL2","10082989-0-COL7","10105494-0-COL1","10105494-0-COL2",
      "ALM93896B2COL31","ALM93896B3COL1","ALM93896B3COL6","ALM93896B3COL11","ALM93896B3COL5",
      "ALM93896B3COL2", "ALM93896B3COL5","ALM93896B3COL3", "ALM10105494B3COL3018", "AL10082989COL3", "AL10105494COL0"]]
Y=linkage(y.T, method='weighted', metric='jaccard') 

plt.subplot(1,2,2)
dendrogram(Y, labels=y.columns, leaf_rotation=90, leaf_font_size=10, show_contracted=True)
plt.rcParams["figure.figsize"] = (25,3)
plt.title('Recent outbreak')
plt.xlabel('Samples')
plt.ylabel('Distance')

#B)
plt.subplot(1,2,1)
dendrogram(Z, labels=io.columns, leaf_rotation=90, leaf_font_size=10, truncate_mode='lastp', p=11, show_leaf_counts=True)
fig.tight_layout()

plt.title('Pruned tree')
plt.xlabel('Samples')
plt.ylabel('Distance')

plt.show()


#It's important to pay attention to the different distance axes
```

![Pruned tree](https://user-images.githubusercontent.com/77884314/124262174-fe399880-db31-11eb-874a-758be07bd8e9.jpg)
 )
