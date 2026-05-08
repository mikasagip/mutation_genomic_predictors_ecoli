The following are original data files and where they came from:

* The tjaden_transcripts_data came from Tjaden, 2023, Escherichia coli transcriptome assembly from a compendium of RNA-seq data sets. <https://dataverse.harvard.edu/api/access/datafile/6767650>

* The schmidt_2016_tables9 came from Schmidt et al, 2016, The quantitative and condition-dependent Escherichia coli proteome. <https://pmc.ncbi.nlm.nih.gov/articles/PMC4888949/#SM> . The lb_protein_data is the export of just table S9 from this original.

* The paxdb_protein_data came from <https://pax-db.org/dataset/511145/3616268379/> . The paxdb_protein_data_no_top_text had all the top comments from paxdb_protein_data removed manually.

* The msystems_essentiality data came from <https://journals.asm.org/doi/10.1128/msystems.00896-22> . This has both the PEC essentiality dataset and the Gerdes essentiality dataset associated with locus tags. The original version can be found at PEC: <https://shigen.nig.ac.jp/ecoli/pec/download.jsp> . Gerdes: <https://pubmed.ncbi.nlm.nih.gov/13129938/>

The collect_properties.py file takes in the following files in: 
* tjaden_transcripts_data
* lb_protein_data
* paxdb_protein_data_no_top_text
* msystems_essentiality
* mut_density

and outputs the all_properties.csv file. 

The skewness figure was moved to the figures folder manually

The all_properties.csv file was manually moved into the folder stats_and_graphs once it was created