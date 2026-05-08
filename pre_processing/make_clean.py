import glob
import re

gff_list = glob.glob("ncbi*.csv")
#imports in the relevant data files 

for f in gff_list:
    f_o = open(f, "r")
    f_r = f_o.read()
    f_o.close()
    #opens file, reads, and closes

    all_lines = f_r.splitlines()
    #splits line of data 

    dict_lt2gn = {}
    dict_lt2start = {}
    dict_lt2end = {}
    dict_lt2strand = {}
    dict_lt2type = {}
    #create all necessery dictionaries

    for ln in all_lines:
        first_char = ln[0]
        #saves the first character of the line 

        ln = ln + "£"
        #adds a pound sign to the end of each line 
        if first_char != "#":
            #this is because some info about the data set at the top, all those lines start with #
            bits = ln.split(",")
            #splits the line according to ,
            genome = bits[0]
            seq_type = bits[2]
            strt = int(bits[3])
            nd = int(bits[4])
            strnd = bits[6]
            annot = bits[8]
            #names the parts of the data, some specifically as int
            

            pat_lt = re.compile("locus_tag=(b[0-9]+?)[^0-9]")
            #create the pattern for locus tags
            if re.search(pat_lt,annot):
                #searches for the pattern in annotation
                find_lt = re.search(pat_lt,annot)
                #returns a matching object
                loc_tag = find_lt.group(1)
                #saves the first captured object as loc_tag
                
                if seq_type == "gene":
                    #if it is specified as a gene
                    dict_lt2start[loc_tag] = strt
                    dict_lt2end[loc_tag] = nd
                    dict_lt2strand[loc_tag] = strnd
                    #saves the start, end, and seq type into relevant dictionaries with the loc_tag as the key
                else:
                    if loc_tag not in dict_lt2type.keys():
                        dict_lt2type[loc_tag] = seq_type
                        #if this specific loc_tag is not already in, add it with the sequence type

                pat_gn1 = re.compile("gene=(.*?);")
                #creates pattern foe gene find
                if re.search(pat_gn1,annot):
                    #searches for the pattern in annotation
                    find_gn = re.search(pat_gn1,annot)
                    #returns a matching object
                    gn1 = find_gn.group(1)
                    #saves the first captured object as loc_tag
                    try:
                        #print(loc_tag)
                        if gn1 not in dict_lt2gn[loc_tag]: 
                            dict_lt2gn[loc_tag].append(gn1)
                            #if gene name isnt already in the dictionary for this key, append it
                    except:
                        dict_lt2gn[loc_tag] = []
                        dict_lt2gn[loc_tag].append(gn1)
                        #if list doesnt exist, create it and append

                pat_gn1a = re.compile("Name=(.*?);")
                if re.search(pat_gn1a,annot):
                    find_gn = re.search(pat_gn1a,annot)
                    gn1a = find_gn.group(1)
                    try:
                        #print(loc_tag)
                        if gn1a not in dict_lt2gn[loc_tag]: 
                            dict_lt2gn[loc_tag].append(gn1a)
                    except:
                        dict_lt2gn[loc_tag] = []
                        dict_lt2gn[loc_tag].append(gn1a)
                    #Same as previous loop but for diff gene name pattern
                pat_gn2 = re.compile("gene_synonym=(.*?);")
                if re.search(pat_gn2,annot):
                    find_gn = re.search(pat_gn2,annot)
                    gn2 = find_gn.group(1)
                    syns = gn2.split(":")
                    for s in syns: 

                        try:
                            #print(loc_tag)
                            if s not in dict_lt2gn[loc_tag]: 
                                dict_lt2gn[loc_tag].append(s)
                        except:
                            dict_lt2gn[loc_tag] = []
                            dict_lt2gn[loc_tag].append(s)
                        #Same again but for diff name pattern
    print(f"there are {len(dict_lt2gn)} in {f}")
    #prints the number of keys in the each dictionary
    if re.search("3", f):
        f_out = open("genes_3.csv", "w")
    else: 
        f_out = open("genes_2.csv", "w")
    f_out.write("genome,locus_tag,start,end,strand,type,gene_names\n")
    #finds which ref genome it is in the file name and creates the out_file accordingly
    #the writes headers to these new files

    for key in dict_lt2start:
        if key in dict_lt2type:
            gene_list = ";".join(dict_lt2gn[key])
            f_out.write(f"{genome},{key},{dict_lt2start[key]},{dict_lt2end[key]},{dict_lt2strand[key]},{dict_lt2type[key]},{gene_list}\n")

    f_out.close()




