Mutation data files are data_foster and data_wei. For ref genome NC_000913.2 annotation file is genes_2, and sequence file is sequence_2.For ref genome NC_000913.3 annotation file is genes_3, and sequence file is sequence_3.

The mutation_density.py file takes in all 6 of these as inputs and outputs  the mut_density.csv file. 

The selection_test.py file takes in all 6 of these as inputs and outputs all the rates_{dataset_name}_pos{1,2,3} files and the selection_test_results_py.csv. 

The selection_stats R file takes in the rates_{dataset_name}_pos{1,2,3} files created bt selection_test.py and outputs the selection_test_respots_r.csv file, the formatted_rates_table, as well as the qq_composite and the hist_composite figures. The table and two figures were manually moved to the figures folder. 

The mut_density.csv file was moved into the properties folder manually.