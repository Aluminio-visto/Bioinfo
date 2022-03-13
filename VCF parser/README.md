## VCF parser

This script is part of a technical test I had to solve for a bioinformatics job. 

VCF files are the output of a typical NGS pipeline. In these files, SNPs and InDels in every chromosome position are listed for one specimen/colony/individual so that, comparing variants among individuals give us an idea of the relationship or the evolutive distance among them. After calculating these differences or "distances" among individuals we would be able to build their phylogenetic tree and estimate if they belong to a recent outbreak, having a close common ancestor or if they are more likely to be distant relatives.

The task was divided in 5 sub-tasks:

#### Iteration 0 | Download repository

#### Iteration 1 | Parse VCF files to table/dataframe

Create a function to read a single VCF

#### Iteration 2 | Extract relevant information from parsed VCF

M. tuberculosis is a haploid organism but the variants are usually called (variant calling step) as diploid, hence you will see the usual diploid genotyping (0/0, 0/1, 1/1).
With the correct information analyzed, filter the SNPs actually present on each sample, this can be a different function.

```python
def parsear_uno(file, incluir_ambiguos=True, completo=True, incluir_INDEL=True):    #incluir_ambiguos nos permite seleccionar -o no- las variantes marcadas 0/1
    
    tabla=[]  
    encabezado=[]
    
    with open(file) as file:
        f=file.readlines()                           #lee el archivo VCF línea por línea
   
        for line in f:                               #y en cada una de ellas
            
            if re.match("^#CHROM", line):            #si comienza por #CHROM es el encabezado
                encabezado=line.split("\t")
                muestra=encabezado[-1].rstrip("\n") #y la última columna del encabezado es el nombre de la muestra


            if re.match("[^#]", line):                #todo lo que no empieza por # es una variante
                linea=line.split("\t")                #generamos una sub-lista: solo nos interesan algunos atributos
                alt=linea[4]                          #la 5ª columna describe la(s) mutacion(es)
                alt=alt.split(',')                    #separa los alelos alternativos múltiples (si los hubiere)
                if incluir_ambiguos==False:           #Para que evite los marcados como heterocigotos por GATK
                    if re.match("[^0]", linea[-1]):   #si no está marcado como heterocigoto (0/1, artefacto de secuenciación)
                        if (incluir_INDEL==True):     #si podemos incluir INDEL (por defecto), meter en la lista ALT
                            tabla.append(linea)       #metemos la linea completa
                        else:
                            for i in range(len(alt)):
                                if (len(alt[i])==1) and (len(linea[3])==1):   #consideramos que los SNP miden 1 base en ambos y los INDEL no
                                    tabla.append(linea) #solo si los strings en Ref y el/los de Alt miden uno se incluye la línea
                else:
                    tabla.append(linea)               #si le pedimos en la función que nos acepte los ambiguos no busca re.match

        df= pd.DataFrame(tabla, columns=encabezado)
        df.index=[muestra]*len(df)                   #esto nos será de utilidad cuando pasemos todos los archivos de golpe
        
        print("Código de muestra:", muestra)         #esto en realidad es prescindible
    if completo==True:
        return df                                 #por defecto nos muestra el dataframe completo
    
    else:
        return df[["POS","REF","ALT"]]            #esto para cuando solo necesitemos mostrar posicion y alelo referencia/mutado
            
```

|#CHROM|	POS	|ID	|REF	|ALT	|QUAL	|FILTER	|INFO	|FORMAT	|
|------|------|---|-----|-----|-----|-------|-----|-------|
|ALM23685B3COL6	|MTB_anc	|1701	|.	|T	|C	|140323.33	|PASS	|AC=2;AF=1.00;AN=2;BaseQRankSum=0.814;DP=95;Exc...	|GT:AD:DP:GQ:PL	|1/1:0,95:95:99:3679,283,0\n|
|ALM23685B3COL6	|MTB_anc	|2532	|.	|C	|T	|47660.33	|PASS	|AC=2;AF=1.00;AN=2;DP=38;ExcessHet=3.0103;FS=0....	|GT:AD:DP:GQ:PL	|1/1:0,38:38:99:1460,114,0\n|
|ALM23685B3COL6	|MTB_anc	|8040	|.	|G	|A	|88799.33	|PASS	|AC=2;AF=1.00;AN=2;DP=73;ExcessHet=3.0103;FS=0....	|GT:AD:DP:GQ:PL	|1/1:0,73:73:99:2791,220,0\n|
|...|



#### Iteration 3 | Combine present SNP into a presence matrix

Merge extracted information into a matrix to keep track of relevant info such as:
-Sample name
-Position
-Mutation (Reference allele and Alternate allele).

```python
def lista_muts(file, incluir_INDEL=True):   #con esto sacamos una lista posicion-mutacion con indice=nombre de muestra
                                            #a partir de esta tabla poodemos sacar de fomra sencilla la tabla pres/ausencia
    tabla=[]                                           
    
    with open(file) as file:
        f=file.readlines()
   
        for line in f:
            if re.match("^#CHROM", line):             #de la linea de encabezado ya solo nos interesa el nombre de la muestra
                encabezado=line.split("\t")           #mejor usar este nombre que el de la ruta (facil de modificar por error)
                muestra=encabezado[-1].rstrip("\n")   #saca el código de muestra
                
            if re.match("[^#]", line):                #las lineas que no empiezan por #son variantes
                linea=line.split("\t")                #convertimos la linea en lista, solo nos interesan algunos atributos
                alt=linea[4]                          #la 5ª columna describe la(s) mutacion(es)
                alt=alt.split(',')                    #separa los alelos alternativos múltiples (si los hubiere)
                if re.match("[^0]", linea[-1]):       #si el genotipo no es ambiguo (0/1 en la última columna)
                        if (incluir_INDEL==True):      #si podemos incluir INDEL (por defecto), meter en la lista ALT
                            for i in range(len(alt)):
                                posic=linea[1] + '(' + linea[3] + '>' + alt[i] + ')' #saca la posicion y la mutacion
                                tabla.append(posic)
                        else:
                            for i in range(len(alt)):
                                if (len(alt[i])==1) and (len(linea[3])==1):   #consideramos que los SNP miden 1 base en ambos y los INDEL no
                                    posic=linea[1] + '(' + linea[3] + '>' + alt[i] + ')' #saca la posicion y la mutacion
                                    tabla.append(posic)
                                
                                
            
    df= pd.DataFrame(tabla, columns=["Posición (Ref>Alt)"])  #en realidad nos ha quedado un pd.Series (1D) 
    df.index=[muestra]*len(df)                   
    
    return df
    
   def sacar_multidf(path_in):       #con esta función unimos todas las tablas posicion-mutacion como paso previo a la tabla final
    
    diccionario={}
    file_list=unzip(path_in)
    
    for file in file_list:
        
        diccionario[file]=lista_muts(file, incluir_INDEL=True) 
        #si seleccionamos True or False veremos que hay 48 filas de diferencia, o 48 INDEL entre las 26 muestras
    
    multi_df=pd.concat(list(diccionario.values()))
    pd.set_option('display.max_rows', 50)
    return(multi_df)
        
```

|Sample| Position (Ref>Alt)|
|------|-------------------|
|10082989-0-COL14	|1701(T>C)|
|10082989-0-COL14	|2532(C>T)|
|10082989-0-COL14	|8040(G>A)|
|10082989-0-COL14	|9143(C>T)|
|10082989-0-COL14	|13460(G>A)|
...	...
|ALM93896B3COL6	|4397324(G>C)|
|ALM93896B3COL6	|4407588(C>T)|
|ALM93896B3COL6	|4408156(A>C)|
|ALM93896B3COL6	|4408920(G>A)|
|ALM93896B3COL6	|4408923(T>C)|

The preferred format is a binary Presence/Absence matrix

```python
def tabla_io(input_path):                                     #con la tabla de muestras-posiciones-mutaciones podemos sacar
                                                              #la tabla presencia/ausencia
    df=sacar_multidf(input_path)
    tabla_final=pd.crosstab(df["Posición (Ref>Alt)"], df.index)
    pd.set_option('display.max_rows', 1000)
    pd.set_option('display.max_columns', 26)
    return tabla_final
```
| Pos\Sample	|10082989-0-COL14	|10082989-0-COL2	|10082989-0-COL3	|10082989-0-COL7	|...|
|-------------|-----------------|----------------|-------------------|---------------|---|
|1012322(G>A)	|1	|1	|1	|1	|...|
|1014815(G>T)	|1	|1	|0	|0	|...|
|1024346(G>A)	|1	|1	|1	|0	|...|
|...|...|...|...|...|...|

#### Iteration 4 | Calculate the SNP distance between all samples

Determine the pairwise distance between each pair of samples
Determine simple pairwise distance (number of SNPs)
```python
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
```
Pairwise Distance by Jaccard Index
```python
from sklearn.metrics import jaccard_score

def distancias_Jaccard(tabla_io):
    n=len(tabla_io.columns)
    tabla_Jaccard=np.zeros((n,n))
    for i in range(len(tabla_io.columns)):
        for j in range(len(tabla_io.columns)):
            resta_jac=jaccard_score(tabla_io.iloc[:,j], tabla_io.iloc[:,i], average='macro')
            tabla_Jaccard[i,j]=resta_jac
    return tabla_Jaccard
    
#que ignora las posiciones que no comparte ninguna de las dos muestras (en realidad en nuestro ejemplo esto es menos 
#práctico, pero si empleásemos todas las posiciones genómicas (mutadas y no mutadas) sería de mayor aplicación.

jac=pd.DataFrame(data=distancias_Jaccard(io), index=io.columns, columns=io.columns)
jac
```

|Position\Sample	|10082989-0-COL14	|10082989-0-COL2	|10082989-0-COL3	|10082989-0-COL7|...|
|-------------|----------|-----------------|---------------|-----------|----|
|10082989-0-COL14	|1.000000	|0.784610	|0.530348	|0.597857	|...|
|10082989-0-COL2	|0.784610	|1.000000	|0.512137	|0.784265	|...|
|10082989-0-COL3	|0.530348	|0.512137	|1.000000	|0.509008 |...|
|10082989-0-COL7	|0.597857	|0.784265	|0.509008 |1.000000 |...|
|...|...|...|...|...|...|

With this we may already generate the distance graph:

```python

import matplotlib.pyplot as plt
import numpy as np
from numpy import tril
from scipy.spatial.distance import pdist, squareform

n_clusters=len(io.columns)


for index, metric in enumerate(['braycurtis', 'hamming', 'jaccard']):           #Bray-Curtis suele emplearse más en metagenómica (cuando tenemos medidas de abundancias)
    avg_dist = np.zeros((n_clusters, n_clusters))                               #la tabla de distancias
    plt.figure(figsize=(14, 14))
    avg_dist=squareform(scipy.spatial.distance.pdist(io.T, metric=metric))      #importante usar la matriz traspuesta para que no saque las distancias entre posic.
    avg_dist=tril(avg_dist)                                                     #es prescindible, solo una cuestión estética
    avg_dist /= avg_dist.max()
    plt.imshow(avg_dist, cmap=plt.cm.gnuplot2, vmin=0)
    plt.xticks(range(n_clusters), io.columns, rotation=90, fontsize=8)
    plt.yticks(range(n_clusters), io.columns)
    plt.colorbar()
    plt.suptitle("%s distances" % metric, size=24)
    plt.tight_layout()
```
![output](https://user-images.githubusercontent.com/77884314/158076163-4ff5609a-e160-4ee5-8fdb-9dee67a3fed7.png)


#### Iteration 5 | BONUS - Include INDELS

We have been using the term SNP distance but INDELS are also useful as phylogenetic marker
Add subtle changes to the functions to include INDELS in the distance calculation

#### Iteration 6 | BONUS - Represent distance in a phylogenetic tree

You can represent this distance in a dendrogram, using any method you find suitable.

```python
#para el cálculo del dendrograma voy a emplear la matriz de distancias que se obtiene sacando el índice de Jaccard 
#y el algoritmo WPGMA 'weighted' que arroja una correlación óptima cluster-distancias reales. Si las distancias fueran euclid
#podríamos usar otros algoritmos un poco más complejos como 'ward' uqe contemplan covarianzas entre muestras
#pero creo que por ahora es más que suficiente.

from scipy.cluster.hierarchy import dendrogram, linkage, cophenet

Z=linkage(io.T, method='weighted', metric='jaccard') #usamos el algoritmo WPGMA para ir buscando el clustering óptimo.

#print(Z[:30]) #podemos observar cómo la 3ªcolumna (distancia) va convergiendo hacia 1 
             #y la 4ªcolumna(nºmuestras) aumenta hasta incluir las 26 elementos.

c, coph_dists = cophenet(Z, pdist(io.T))  #esta f(x) nos indica la correlación entre nuestra matriz de distancias 
            #y las distancias que implicaría el clustering por nuestro método (98.6% en nuestro caso es casi perfecto)
plt.figure(figsize=(25, 10))
plt.title('Dendrograma final', size =40)
plt.xlabel('Muestras', size=18)
plt.ylabel('Distancias', size=18)
dendrogram(Z, labels=io.columns, leaf_rotation=90, leaf_font_size=18)
plt.show()

```

![image](https://user-images.githubusercontent.com/77884314/158075238-8e140924-30ff-4929-b482-fd0dc5f4378d.png)

