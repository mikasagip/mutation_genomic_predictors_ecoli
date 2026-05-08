import glob 
import shutil
import os
from pathlib import Path
import re

fname = "data_wei.csv"
out_file = open(fname, "w")
#creates a variable with file name 
#opens file under write mode and names it 

out_file.write("LineName,Position,RefBase,ObservedBase\n")
#writes the new headers in the outfile, these are the same as foster 2015 data headers

all_wei = "all_wei2.csv"
f_o = open(all_wei, "r") 
f_data = f_o.read()
f_o.close()
#opens the data file, reads the data and closes

data = f_data.splitlines()
#splits the lines in the wei original file 


#here we are going to filter the data to exclude inapproriate data
cnt = 0

for line in data:
    cnt = cnt +1 
    
    if cnt == 1:
        continue
    #excludes the header from being filtered out
    
    if "A" in line[0] or "B" in line[0]:
        continue 
    #filters out any data lines which are of MMR- types 
    
    split = line.split(',')
    #splits the string into strings in a list 
    
    mutation_raw = split[4]
    if mutation_raw == "-":
        continue
    #if the two codon mutation data is not available, skips that line
    #this filters out any intergenic areas
    
    codons = mutation_raw.split('->')
    ref_codon = codons[0]
    mut_codon = codons[1]
    #seperates the codon mutation into two variables, the original codon
    #and the mutated codon
    
    clone_name = split[0] 
    line_name = split[1]
    position = split[2] 

    good_data = [f"{clone_name}_{line_name}", f"{position}"]
    #writes the clone and line name together into new list
    #also writes the poistion into the new list
    
    base_cnt = 0
    for base in ref_codon:
        if base != mut_codon[base_cnt]:
            good_data.append(base)
            good_data.append(mut_codon[base_cnt])
            base_cnt = base_cnt + 1
        else:
            base_cnt = base_cnt + 1
    #starts another counter 
    #iterates for each base in ref_codon
    #Checks if it is the same base as the same position in mut_codon
    #if it is different, appends the ref base and mut base to good data 
    
    out_file.write(",".join(good_data) + "\n")
    #writes the data into the outfile in csv format
    
    
out_file.close()
    
    

    
    
    
    
    

