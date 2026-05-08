from collections import defaultdict
import re
from scipy.stats import skew
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

f_o = open("mut_density.csv", "r")
mut_density_raw = f_o.read()
f_o.close()

f_o1 = open("paxdb_protein_data_no_top_text.txt", "r")
protein_raw = f_o1.read()
f_o1.close()

f_o2 = open("msystems_essentiality.csv", "r")
essentiality_raw = f_o2.read()
f_o2.close()

f_o3 = open("tjaden_transcripts_data.csv", "r")
transcription_raw = f_o3.read()
f_o3.close()

f_o4 = open("lb_protein_data.csv", "r")
protein_lb_raw = f_o4.read()
f_o4.close()
#opens all the data files

fname = "all_properties.csv"
out_file = open(fname, "w")
#create and outfile and open it in w mode

out_file.write(f"locus_tag,mutation_density,mutation_count,start,middle,end,avg_length,strand,type,gc_content,cpg_excess,asymmetry,binary_asymmetry,distance_from_oric,relative_distance_from_oric,relative_distance_oscillating_triangle,relative_distance_oscillating_cosine,essentiality_pec,binary_pec,essentiality_gerdes,binary_gerdes,transcript,protein_abundance_pax,protein_abundance_lb,gene_names\n")
#writes headers to outfile

mut_density_lines = mut_density_raw.splitlines()
mut_density_lines = mut_density_lines[1:]
#split lines in density data and remove header

ess_lines = essentiality_raw.splitlines()
ess_lines = ess_lines[2:]
#split lines in essentiality data 
#this has two header lines

trans_lines = transcription_raw.splitlines()
trans_lines = trans_lines[1:]
#split lines in transcription fdata and remove header

pro_lines = protein_raw.splitlines()
pro_lines = pro_lines[1:]
#splits lines in protein data and removes header

pro_lb_lines = protein_lb_raw.splitlines()
pro_lb_lines = pro_lb_lines[3:]
#splits lines in lb protein data 
#this file has 3 header lines 

dict_lt2trans = defaultdict(list)
#makes a dictionary for locus tag to transcript abundance data

dict_lt2pro = {}
#creates a dictionary for locus tag to protein abundance data

dict_gn2pro = {}
#creates a dictionary for gene name to protein abundance data

dict_lt2ess = {}
#makes dictionary for locus tag to essen data

#not_num_trans = defaultdict(list)
#creates dictionary for checking if all trnacripts are number values

for ln in trans_lines:
    bits = ln.split(",")
    loc_tag_trans = bits[0]

    for transcript in bits[7:]:
        if transcript == "-":
            transcript = 0.0
            #in RNAseq - represents no transcripts found, therefore 0
        else: 
            transcript = float(transcript)
            #modify into a float
        
        dict_lt2trans[loc_tag_trans].append(transcript)
        #add the transcript value into the dict

        #try:
            #val = float(transcript)
            #dict_lt2trans[loc_tag].append(val)
            #if transcript number is a float adds to dictionary
        #except ValueError:
            #not_num_trans[loc_tag].append(transcript)
            #if transcript is not a float adds to bad dict

#print(f"there are {len(not_num_trans)} loc tags with bad transcripts")
#to check how many non-numerical trasncripts there may be
#there are none

sk_list = []
#make list for skewness coefficients 

dict_lt2transmedian = {}
#make dictionary for transcription medians

skew_big = []
skew_ok = []

for loc_tag1, values in dict_lt2trans.items():

    values = np.array(values)
    #make values into an array 

    trans_median = np.median(values)
    #calculates the median transcript abundance for the loc tag
    dict_lt2transmedian[loc_tag1] = trans_median
    #adds to the medians dictionary

    sk = skew(values)
    #computes the skewness of the distribution
    if sk > 1:
        skew_big.append(sk)
    else:
        skew_ok.append(sk)
    sk_list.append(sk)
    #adds the skewness coefficient into the list

print(f"There are {len(skew_big)} genes with skewness greater than 1")
print(f"There are {len(skew_ok)} genes with ok skewness")

sk_list = np.array(sk_list)
#make list into an array

plt.hist(sk_list, bins = 60, color = "#f7c767", edgecolor="black")
plt.axvline(1, linestyle="--", color = "blue")
plt.xlabel("Moment coefficient of skewness")
plt.ylabel("Count")
plt.show()
#makes a histogram of the skewness data

for ln1 in ess_lines:

    bits_ess = ln1.split(",")
    #splits the line 

    loc_tag = bits_ess[5]
    pec = bits_ess[8]
    gerdes = bits_ess[9]
    #captruring loc tag, as well as two sets of essentiality data (good pec KO data, and bad gerdes tn-seq data)

    dict_lt2ess[loc_tag] = {
        "pec": pec,
        "gerdes": gerdes
    }
    #add both to the dictionary

for ln2 in pro_lines:
    bits1 = ln2.split("\t")
    #splits the lines based on tab 

    middle = bits1[1]
    middle = middle.split(".")
    #split this item based on .

    loc_tag = middle[1]
    #find the loc tag

    pro = bits1[2]
    #saves the protein abundance 

    if pro == 1e-07:
        pro = 0.0000001

    dict_lt2pro[loc_tag] = pro
    #saves this info into dict

for ln3 in pro_lb_lines:
    bits3 = ln3.split(",")
    #split the line
    gene_name = bits3[2]
    gene_name = gene_name.lower()
    #make it lower case for later comparison
    pro_lb = bits3[16]
    #gets gene name and protein abundance

    dict_gn2pro[gene_name] = pro_lb
    #adds the abundance to the dictionary 

dict_lt2pro_final = {}
#makes dictionaries to later match locus tag to values currently only associated with gene names

bad_distance_index = []
no_pro = []
no_trans = []
no_essen_pec = []
no_essen_gerdes = []



oric = 3925860
L = 715792
genome_length = oric + L
half_genome = (oric + L)/2
ter = half_genome - L
#to calculate relative distance, take dist / half genome length

for gn in mut_density_lines:
    bits2 = gn.split(",")
    loc_tag = bits2[0]
    mut_den = bits2[1]
    mut_count = bits2[2]
    avg_len = float(bits2[3])
    strt = float(bits2[4])
    nd = bits2[5]
    strand = bits2[6]
    seq_type = bits2[7]
    gc_content = bits2[8]
    cpg_excess = bits2[9]
    gene_names = bits2[10]
    #split the line into components and name them 

    middle = strt + (avg_len/2)
    #finds the middle position of the gene

    if middle >= 0 and middle <= ter :
        replichore = "right"
    elif middle > ter and middle <= oric :
        replichore = "left"
    elif middle >  oric :
        replichore = "right"
    else:
        bad_distance_index.append(loc_tag)
    #indexes the distance index according to the middle position of the gene
    #assigns right or left replichore based on the middle position of the gene

    raw_dist = abs(middle - oric)
    circular_dist = min(raw_dist, genome_length - raw_dist)
    #calculates the shortest distance on a circle of genome length
    dist_index = circular_dist / half_genome
    #makes the index from 0 to 1 by considering within replichore length

    dist_osc_tri = 1 - abs(dist_index - 0.5) / 0.5
    #calculates relative distance from half of replichore in a linear fashion (traingle)
    dist_osc_cos = np.cos(np.pi * (dist_index - 0.5))
    #calculates relative distance from half of replichore in a curved fashion (using cosine curve)

    if replichore == "right":
        if strand == "+":
            asymmetry = "leading"
            bin_asym = 1
        else: 
            asymmetry = "lagging"
            bin_asym = 0
    elif replichore =="left":
        if strand == "-":
            asymmetry = "leading" 
            bin_asym = 1
        else: 
            asymmetry = "lagging"
            bin_asym = 0
    #associated leading or lagging based on replichore and strand, also numeric binary versions

    gene_names = gene_names.split(";")
    #splits the gene names based on ;
    
    pro_match = 0
    #start essen and pro match counters

    for gene in gene_names:
        gene_lower = gene.lower()
        #lowercase the gene name for the comparison

        if pro_match == 0:
            for dict_gn, value in dict_gn2pro.items():
                if gene_lower == dict_gn:
                    pro_lb = value
                    pro_pax = dict_lt2pro.get(loc_tag, "N/A")
                    dict_lt2pro_final[loc_tag] = {
                        "pro_pax": pro_pax,
                        "pro_lb": value
                    }
                    #add them to the protein final dictionary
                    
                    pro_match += 1
                    #get the matching value and increment counter
                    break 
                    #if match found, stop comparing names
    
    if pro_match == 0: 
        #if no match was found
        pro_pax = dict_lt2pro.get(loc_tag, "N/A")
        dict_lt2pro_final[loc_tag] = {
            "pro_pax": pro_pax,
            "pro_lb": 0
        }
    
    pro_entry = dict_lt2pro_final.get(loc_tag, {"pro_pax": "N/A", "pro_lb": "N/A"})
    pro_pax = pro_entry["pro_pax"]
    pro_lb = pro_entry["pro_lb"]
    if pro_pax == "N/A":
        no_pro.append(loc_tag)
    #gets the proptein data, if the pax DB dataset is missing it flags it and adds in N/A

    essen_entry = dict_lt2ess.get(loc_tag, {"pec": "N/A", "gerdes": "N/A"})
    essen_pec = essen_entry["pec"]
    essen_gerdes = essen_entry["gerdes"]
    if essen_pec == "N/A":
        no_essen_pec.append(loc_tag)
    if essen_gerdes == "N/A":
        no_essen_gerdes.append(loc_tag)
    #gets the essentiality data, if either is missing, flags them and adds in N/A

    trans = dict_lt2transmedian.get(loc_tag, "N/A")
    if trans == "N/A":
        no_trans.append(loc_tag)
    #gets the transcript info, if there is none, adds in N/A

    gene_names = ";".join(gene_names)
    #joins all list items in the gene_names list into one string so the output is organised

    if essen_pec == "NE":
        bin_pec = 0
    elif essen_pec == "E":
        bin_pec = 1

    if essen_gerdes == "NE":
        bin_gerdes = 0
    elif essen_gerdes == "E":
        bin_gerdes = 1

    #add in binary numeric versions of the essentiality data




    if (
    pro_pax != "N/A" and
    essen_pec != "N/A" and
    essen_gerdes != "N/A" and
    trans != "N/A"
    ):
        #filters out all ones that have N/A in any of the fields
        out_file.write(f"{loc_tag},{mut_den},{mut_count},{strt},{middle},{nd},{avg_len},{strand},{seq_type},{gc_content},{cpg_excess},{asymmetry},{bin_asym},{circular_dist},{dist_index},{dist_osc_tri},{dist_osc_cos},{essen_pec},{bin_pec},{essen_gerdes},{bin_gerdes},{trans},{pro_pax},{pro_lb},{gene_names}\n")
        #writes out all info to outfile


print(f"There are {len(bad_distance_index)} gene with a strange distance issue")
print(f"There are {len(no_trans)} gene with no transcript data")
print(f"There are {len(no_pro)} gene with no pax protein data")
print(f"There are {len(no_essen_pec)} gene with no essen pec data")
print(f"There are {len(no_essen_gerdes)} gene with no essen gerdes data")

out_file.close()