import glob 
import re

gff_list = glob.glob("ncbi*.csv")
print(gff_list)

for f in gff_list:
    print("processing", f)
    f_o = open(f, "r")
    f_r = f_o.read()
    f_o.close()

    all_lines = f_r.splitlines()

    dict_lt2type = {}

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
            match = re.search(pat_lt, annot)
            if not match:
                continue
            loc_tag = match.group(1)
            dict_lt2type.setdefault(loc_tag, [])
            if seq_type not in dict_lt2type[loc_tag]:
                dict_lt2type[loc_tag].append(seq_type)

    print(f"there are {len(dict_lt2type)} in {f}")

    count = sum(1 for v in dict_lt2type.values() if len(v) == 1)

    print(f"there are {count} keys with only one type associated")

    single_type_dict = {k: v for k, v in dict_lt2type.items() if len(v) == 1}

    print("Locus tags with exactly one type:")
    for k, v in single_type_dict.items():
        print(f"{k}: {v[0]}")  #v is a list with one element
    



