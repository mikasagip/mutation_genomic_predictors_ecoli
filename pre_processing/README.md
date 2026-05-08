The two GFF files are the original annotation files from NBCI (for NC_000913.2 - <https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.2> for NC_000913.3 - <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/>). These were then exported as csv files. 

The make_clean.py file takes in these csv files and outputs 2 clean csv files called genes_2 and genes_3.

compare_refs.py takes in the same input as make_clean.py and simply compare the seq types of different lines. 

compare_gene_names.py takes genes_2 and genes_3 in as input and compares the genes names associated with each locus tag to ensure the two reference genomes agree that each loc tag is the same gene.

I manually moved the genes_2 and genes_3 into the mutation_density folder after they were produced