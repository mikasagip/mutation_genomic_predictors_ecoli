import glob 
import shutil 
import os
from pathlib import Path
import re

fname = "all_wei.csv1"
out_file = open(fname, "w")
#creates a variable with file name 
#opens file under write mode and names it 

cnt = 0
#starts counter

f_list = glob.glob("*.csv")
#creates a list of files that all end in .csv

for f in f_list:
    cnt = cnt +1
    #increase counter 
    
    f_o = open(f, "r") 
    f_data = f_o.read()
    f_o.close()
    #opens file, reads data from open file and saves under f_data, then closes file
    
    all_lines = f_data.splitlines()
    header = all_lines[0]
    data_lines = all_lines[1:]
    #splits the lines, firet line is saved as header. data_lines is then the rest if the lines 
    
    if cnt == 1:
        out_file.write(f"{header}\n")
    #if it is the first file, write the header to the outfile 
    
    for ln in data_lines:
        out_file.write(f"{ln}\n")
    #write every line in the file other than header to the outfile. \n adds a new line at end. 

out_file.close()
#closes the out file 

path = os.getcwd()
#gets the current working directory and saves it as path
mydir = os.path.basename(path)
#mydir is then just the last part of path (aka folder name)
dest = re.sub(mydir, "Foster_Wei", path)
#creates a new variable called dest, this is path but take out mydir and swap it with string (aka different folder name)
newf = re.sub("1", "", fname)
#newf takes the original outfile and removes the 1 from the end of it
shutil.move(fname, os.path.join(dest, newf))
#os.path.join adds newf to the end of the dest path 
#then shutil moves the fname file to the new destination and renames it to newf variable, and deletes the original file. 