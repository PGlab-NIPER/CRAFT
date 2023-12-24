#main.py file
#This file calls other module of CRAFT to identify the cavities inside a given protein

import os
import shutil
import csv

from craft import ReadPdb
from craft import Base
from craft import Tetrahedral
from craft import Properties as prop
from craft.CraftVariables import craft

cwd = os.getcwd()
PDB_dir = os.path.join(cwd, 'PDB')

print("cwd",os.path.join(cwd,"Results"))
if os.path.exists(os.path.join(cwd,"Results")):
    shutil.rmtree(os.path.join(cwd,"Results"))
os.makedirs(os.path.join(cwd,"Results"))
os.makedirs(os.path.join(cwd,"Results","Cavity_files"))
os.makedirs(os.path.join(cwd,"Results","Property files"))

def rewrite(filename):
    f=open(os.path.join(PDB_dir,filename+".pdb"), 'r')
    lines = f.readlines()
    modified_lines = []
    for i in range(len(lines)):
        line = lines[i].strip()
        word=line[0:6]
        if word=='ANISOU':
            continue
        else:
            modified_lines.append(line)

    f=open(os.path.join(PDB_dir,filename+".pdb"), 'w')
    f.write('  \n'.join(modified_lines))

parameter={"MCR":1.4,
           "Peeling off radius":10,
           "susceptible radius":4,
           "susceptible volume":4,
           "Extend forbidden face": False,
           "Extend boundary":False,
           "cavity atoms":15}
           
def define_parameter():
    def validate_input(prompt, min_val, max_val):
        try:
            value = float(input(prompt)) if min_val != "True" and min_val != "False" else input(prompt)
            if min_val != "True" and min_val != "False":
                if min_val <= value <= max_val:
                    return value
                else:
                    print(f"Please enter a value between {min_val} and {max_val}")
            else:
                if value == "True" or value == "False":
                    return value == "True"
                else:
                    print("Please enter 'True' or 'False'")
        except ValueError:
            print("Please enter a valid input")

    parameter["MCR"] = validate_input("Enter MCR value between 1 to 2: ", 1, 2)
    parameter["Peeling off radius"] = int(validate_input("Enter integer Peeling off radius between 8 to 15: ", 8, 15))
    parameter["susceptible radius"] = int(validate_input("Enter integer susceptible radius between 2 to 5: ", 2, 5))
    parameter["susceptible volume"] = int(validate_input("Enter integer susceptible volume between 2 to 5: ", 2, 5))
    parameter["Extend boundary"] = validate_input("Set Extend boundary true by entering True or False: ", "True", "False")
    parameter["Extend forbidden face"] = validate_input("Set Extend forbidden face true by entering True or False: ", "True", "False")
    parameter["cavity atoms"] = int(validate_input("Enter integer cavity atoms between 10 to 50: ", 10, 50))
    return parameter
       

word = input("Would you like to process a single .pdb file using CRAFT? Please answer with Y or N.")
if word.casefold() == "Y".casefold():
    word = input("want to proceed with default parameter. Please answer with Y or N.")
    if word.casefold()=="N".casefold():
        parameter=define_parameter()
    print("parameter used", parameter)
    fn = (input("Enter the file name with .pdb extension ")).split(".")
    if(fn[1] == "pdb"):
        try:
            rf = open(os.path.join(PDB_dir,fn[0]+".pdb"), "r")
        except IOError:
            print("File not Found ")
        else:
            rewrite(fn[0])
            rPDB=ReadPdb.ReadPdb(rf)
            print("The given protein contains : ")
            rPDB.display_entries()
            rPDB.selection()
            print()
            print("Selected entries are :: ")
            print()
            rPDB.display_selected()
            status=rPDB.read_coords(rf)
            if status==1:
                rPDB.select_position()
                base=Base.Base(cwd)
                tetra=Tetrahedral.Tetrahedral(base)
                tetra.delimiter_tetrahedral(parameter["Peeling off radius"],parameter["Extend boundary"])    #peeling off radius=10
                tetra.susceptible_tetra(parameter["susceptible volume"],parameter["susceptible radius"])       #susceptible tetrahedral volume and radii 4,4
                tetra.cavity_scan(rPDB,parameter["MCR"],parameter["Extend forbidden face"],parameter["cavity atoms"])        #defined MCR==1.4
                prop.Properties(tetra,base,fn[0],parameter["MCR"],cwd)
            else:
                print("File is not in standard format. please provide standard pdb file")
else:
    print("You can process multiple .pdb files by providing pdb files separated by comma example abc1,xyz2,pqr3 ")
    print("Multiple files will be processed using default parameters, considering all chains within each file while excluding ligands.")
    word = input("want to proceed with default parameter. Please answer with Y or N.")
    if word.casefold()=="N".casefold():
        parameter=define_parameter()
    print("parameter used", parameter)
    input_files = list(input("Enter the file names without .pdb extension in a list ").split(','))
    non_std = open(os.path.join(cwd,"Results","Non-standard PDB.csv"), mode='w', newline='')  #store non standard pdbs which not excuted by program
    writer = csv.writer(non_std)
    writer.writerow(["PDB","Status"])
    
    for file in input_files:
        try:
            print(file)
            rf = open(os.path.join(PDB_dir,file+".pdb"), "r")
        except IOError:
            print("File not Found ",file,".pdb")
        else:
            print("file : ",file)
            craft.reset()
            rewrite(file)
            rPDB=ReadPdb.ReadPdb(rf)
            rPDB.multi_selection()
            status=rPDB.read_coords(rf)
            if status==1:
                base=Base.Base(cwd)
                tetra=Tetrahedral.Tetrahedral(base)
                tetra.delimiter_tetrahedral(parameter["Peeling off radius"],parameter["Extend boundary"])    #peeling off radius=10
                tetra.susceptible_tetra(parameter["susceptible volume"],parameter["susceptible radius"])       #susceptible tetrahedral volume and radii 4,4
                tetra.cavity_scan(rPDB,parameter["MCR"],parameter["Extend forbidden face"],parameter["cavity atoms"])        #defined MCR==1.4
                prop.Properties(tetra,base,file,parameter["MCR"],cwd)
            else:
                row=[file]
                row.append(status)
                writer.writerow(row)
                #pass
            
        
        
        
