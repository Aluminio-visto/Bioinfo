## VCF parser

This script is part of a technical test I had to solve for a bioinformatics job. 

VCF files are the output of a typical NGS pipeline. In these files, SNPs and InDels in every chromosome position are listed for one specimen/colony/individual so that, comparing variants among individuals give us an idea of the relationship or the evolutive distance among them. After calculating these differences or "distances" among individuals we would be able to build their phylogenetic tree and estimate if they belong to a recent outbreak, having a close common ancestor or if they are more likely to be distant relatives.

The task was divided in 5 sub-tasks:

#### Iteration 0 | Download repository

#### Iteration 1 | Parse VCF files to table/dataframe

Create a function to read a single VCF

#### Iteration 2 | Extract relevant information from parsed VCF

M. tuberculosis is a haploid organism but the variants were called (variant calling step) as diploid, hence you will see the usual diploid genotyping (0/0, 0/1, 1/1).
With the correct information analyzed, filter the SNPs actually present on each sample, this can be a different function.

#### Iteration 3 | Combine present SNP into a presence matrix

Merge extracted information into a matrix to keep track of relevant info such as:
-Sample name
-Position
-Mutation (Reference allele and Alternate allele).

The preferred format is a binary Presence/Absence matrix

#### Iteration 4 | Calculate the SNP distance between all samples

Determine the pairwise distance between each pair of samples

#### Iteration 5 | BONUS - Include INDELS

We have been using the term SNP distance but INDELS are also useful as phylogenetic marker
Add subtle changes to the functions to include INDELS in the distance calculation

#### Iteration 6 | BONUS - Represent distance in a phylogenetic tree

You can represent this distance in a dendrogram, using any method you find suitable.

To follow along with the matrix format and the use of python, you are encouraged to use Scipy library, specifically linkage and dendrogram

![image](https://user-images.githubusercontent.com/77884314/158075238-8e140924-30ff-4929-b482-fd0dc5f4378d.png)

