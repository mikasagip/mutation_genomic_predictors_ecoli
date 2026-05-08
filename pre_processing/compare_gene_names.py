import glob 
import re
from collections import defaultdict

flist = glob.glob("genes*.csv")

dict_lt2gn2 = defaultdict(list)
dict_lt2gn3 = defaultdict(list)
#create a dictionary for each dataset, makes them be key:list by default

cnt = 2
#start a counter for file number 

for f in flist: 
    f_o = open(f, "r")
    f_r = f_o.read()
    f_o.close()
    #opens, reads, closes

    all_lines = f_r.splitlines()
    #splits the lines of data

    for ln in all_lines[1:]:
        bits = ln.split(",")
        loc_tag = bits[1]
        gnames = bits[6].split(";")
        #names the variables

        if cnt == 2:
            dict_lt2gn2[loc_tag].extend(gnames)
        #if its data set 2, add to dict 2

        if cnt == 3:
            dict_lt2gn3[loc_tag].extend(gnames)
        #if its data set 3, add to dict 3
    if cnt == 2:
        print(f"there are {len(dict_lt2gn2)} in {f}")
    
    if cnt == 3:
        print(f"there are {len(dict_lt2gn3)} in {f}")
    
    cnt = cnt + 1

matches ={}
#Makes a dictionary for writing in the agreed upon name 
no_match_keys = []
#creates a list to store problematic keys

for key in dict_lt2gn2:
    if key in dict_lt2gn3:

        match_found = 0
        #starts counter of match numbers

        for name2 in dict_lt2gn2[key]:
            for name3 in dict_lt2gn3[key]:
                
                if name2.lower() == name3.lower():
                    #checks in a case-insensitive manner
                    matches[key] = name2
                    #appends to the matches dictionary
                    match_found = match_found + 1
            
            if match_found == 1:
                break
                #stop iterating through dict3 names if already found a match

            
        if match_found == 0:
            no_match_keys.append(key)
            #adds the key to the no matches list if doesnt find one
        
print(f"There are {len(no_match_keys)} keys with no matching names")



