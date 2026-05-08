import re
from collections import defaultdict

f_o = open("genes_2.csv", "r")
genes_2 = f_o.read()
f_o.close()

f_o1 = open("genes_3.csv", "r")
genes_3 = f_o1.read()
f_o1.close()

f_o2 = open("data_foster.csv", "r")
data_foster = f_o2.read()
f_o2.close()

f_o3 = open("data_wei.csv", "r")
data_wei = f_o3.read()
f_o3.close()

seq_o = open("sequence_3.fna", "r")
seq_raw = seq_o.read()
seq_o.close()

#import in all the relevant files and store the information from them 

fname = "mut_density.csv"
out_file = open(fname, "w")
#create and outfile and open it in w mode

wei_lines_raw = data_wei.splitlines()
foster_lines_raw = data_foster.splitlines()
gn2_lines_raw = genes_2.splitlines()
gn3_lines_raw = genes_3.splitlines()
seq_lines_raw = seq_raw.splitlines()
#split all the data into lines 

wei_lines = wei_lines_raw[1:]
foster_lines = foster_lines_raw[1:]
gn2_lines = gn2_lines_raw[1:]
gn3_lines= gn3_lines_raw[1:]
#removes the headers from all data

whole_seq = "".join(seq_lines_raw[1:])
whole_seq = re.sub('[^A-Za-z]', "", whole_seq).upper()
#joins all sequence lines other than title line and ensures all is upper case

def reverse_complement(sequence):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    rev_seq = reversed(sequence)
    good_seq = []
    for b in rev_seq:
        good_b = comp[b]
        good_seq.append(good_b)
    good_seq = "".join(good_seq)
    return good_seq
#defines function for creating a reverse complement of the sequence, to be used for - strand genes 

def calc_gc_cpg(seq):

    length = len(seq)

    g = seq.count("G")
    c = seq.count("C")
    #counts all g and c

    gc_content = (g + c) / length
    #calculates gc content based on length and count

    obs_cpg = seq.count("CG")
    #counts cg pairs

    
    if g == 0 or c == 0:
        exp_cpg = 0
        #if not g or c then expected is 0
    else:
        exp_cpg = (c * g) / length

    if exp_cpg == 0:
        cpg_excess = 0
        #this would be because there were no g or c
    else:
        cpg_excess = obs_cpg / exp_cpg
        #calculates the excess

    return gc_content, cpg_excess

mut_list_foster = []
mut_list_wei = []
#creates lists to store mutations

for ln in foster_lines:
    bits = ln.split(",")
    position = int(bits[1])
    #extracts the positions from the data
    mut_list_foster.append(position)
    #adds the position to the mutations list

for ln in wei_lines:
    bits1 = ln.split(",")
    position1 = int(bits1[1])
    #extracts the positions from the data
    mut_list_wei.append(position1)
    #adds the position to the mutations list

dict_lt2mut = defaultdict(list)
#creates a dictionary for locus tags to mutatuions lists
dict_lt2mut2 = defaultdict(list)
#sets up dict for mutations in .2
dict_lt2lens = defaultdict(list)
#creates a dictionary for locus tags to gene length from both data sets
dict_lt2lens2 = {}
#sets up dictionary for lengths of gene in .2

dict_lt2annot2 = {}
dict_lt2annot3 = {}
#makes dictionaries for storing annotation info

for gn in gn2_lines:
    bits2 = gn.split(",")
    loc_tag2 = bits2[1]
    strt2 = int(bits2[2])
    nd2 = int(bits2[3])
    strand2 = bits2[4]
    seq_type2 = bits2[5]
    gene_names2 = bits2[6]

    gn_len2 = nd2 - strt2 + 1
    #calc gene length 

    dict_lt2annot2[loc_tag2] = {
        "strand": strand2,
        "type": seq_type2,
        "gene_names": gene_names2
    }
    #adds loc tag to list

    dict_lt2lens2[loc_tag2] = (gn_len2)
    #adds the gene length to the dictionary

    for mut in mut_list_foster:
        if mut >= strt2 and mut<= nd2:
            dict_lt2mut2[loc_tag2].append(mut)
            #if the mutation from foster is between the start and end of the gene, adds the mutation position into the dictionary


for gn in gn3_lines:
    bits3 = gn.split(",")
    loc_tag3 = bits3[1]
    strt3 = int(bits3[2])
    nd3 = int(bits3[3])
    strand3 = bits3[4]
    seq_type3 = bits3[5]
    gene_names3 = bits3[6]

    gn_len3 = nd3 - strt3 + 1
    #calc gene length 

    dict_lt2annot3[loc_tag3] = {
        "strand": strand3,
        "type": seq_type3,
        "gene_names": gene_names3,
        "start": strt3,
        "end": nd3
    }
    #adds loc tag to dict

    if loc_tag3 in dict_lt2lens2:
        #makes sure the loc tags are common
        dict_lt2lens[loc_tag3].append(gn_len3)
        gn_len2 = dict_lt2lens2[loc_tag3] 
        dict_lt2lens[loc_tag3].append(gn_len2)
        #adds the gene lengths to the dictionary
        for mut in mut_list_wei:
            if mut >= strt3 and mut<= nd3:
                dict_lt2mut[loc_tag3].append(mut)
                #if the mutation from wei is between the start and end of the gene, adds the mutation position into the dictionary
        
        dict_lt2mut[loc_tag3].extend(dict_lt2mut2.get(loc_tag3, []))
        #then adds all of the foster mutations for the same locus tag into the dictionary to combine the two
                
    
#so now we have a dictionary that maps locus tags with mutations using the appropriate reference annotations

dict_lt2annot_both = {}

for loc_tag in dict_lt2annot2:
    if loc_tag in dict_lt2annot3:
        strand_both = dict_lt2annot3[loc_tag]["strand"]
        type_both = dict_lt2annot3[loc_tag]["type"]
        gene_names_both = dict_lt2annot3[loc_tag]["gene_names"]
        start_both = dict_lt2annot3[loc_tag]["start"]
        end_both = dict_lt2annot3[loc_tag]["end"]

        dict_lt2annot_both[loc_tag] = {
        "strand": strand_both,
        "type": type_both,
        "gene_names": gene_names_both,
        "start": start_both,
        "end": end_both
    }
        
#compares the locus tags and only saves the relevant information for ones that are in both datasets
#this is so we dont get biases in the results that show less mutation in ones that are only present in one or the other


dict_lt2avglen = {}

for loc_tag, lengths in dict_lt2lens.items():
    avg_len = sum(lengths) / len(lengths)
    #calculates the average length

    dict_lt2avglen[loc_tag] = avg_len
    #adds the locus tag and average length into the dictionary


out_file.write("locus_tag,mutation_density,mutation_count,avg_length,start,end,strand,type,gc_content,cpg_excess,gene_names\n")
#writes headers to the out file


for loc_tag in dict_lt2avglen:

    seq_type =dict_lt2annot_both[loc_tag]["type"]
    #retrieves the sequence type

    if seq_type == "CDS":
        #filters out non-CDS

        strand = dict_lt2annot_both[loc_tag]["strand"]
        gene_names =dict_lt2annot_both[loc_tag]["gene_names"]
        start = dict_lt2annot_both[loc_tag]["start"]
        end = dict_lt2annot_both[loc_tag]["end"]
        #retrieves the annotationa from the dictionary

        mut_count = len(dict_lt2mut[loc_tag])
        avg_len = dict_lt2avglen[loc_tag]
        #gets the mutation count and gene length

        den = mut_count / avg_len
        #calculates the mutation density

        start0 = start - 1
        end0 = end 
        #convert start and end to python indexing

        gene_seq = whole_seq[start0:end0]
        #retrieves the gene sequence from the whole genome sequence

        if strand == "-":
            gene_seq = reverse_complement(gene_seq)
            #takes the reverse complement if the strand is -

        gc_content, cpg_excess = calc_gc_cpg(gene_seq)
        #uses the previously defined function to calculate gc content and cpg excess

        out_file.write(f"{loc_tag},{den},{mut_count},{avg_len},{start},{end},{strand},{seq_type},{gc_content},{cpg_excess},{gene_names}\n")


out_file.close()
