#readpdb.py file
#File to read the standard .pdb format to estimate cavities inside it

from .CraftVariables import craft

class ReadPdb:
    def __init__(self, rf):   
        self.anisou=False
        
        word=rf.read(6).strip()
        at_count=0
        while(word!= "HETNAM" and word!= "ATOM"): 
            line=rf.readline()
            word=rf.read(6).strip()
        if(word == "HETNAM"):                   #Search for HETNAM in given pdb
            while(word!= "FORMUL"):
                line=rf.readline()
                line=" ".join(line.split())
                line=line.split(" ")
                if(len(line[0])==3):
                    craft.het_name[line[0]]=" ".join(line[1:])
                else:
                    if line[1] in craft.het_name.keys():
                        craft.het_name[line[1]]+=" ".join(line[2:])
                word=rf.read(6).strip()
            while(word!= "ATOM"):               #Search for chains present in given pdb
                word=rf.read(6).strip()
                if(word=='ATOM'):
                    rf.read(15)
                    craft.chain.append(rf.read(1))
                    at_count+=1                 #at_count=variable to count atoms present in each chain or ligand 
                line=rf.readline()
        else:
            if(word=='ATOM'):
                rf.read(15)
                craft.chain.append(rf.read(1))
                line=rf.readline()
                at_count+=1
        word=rf.read(6).strip()
        while(word!= "END"):
            if(word == "TER"):
                #craft.chain_at[craft.chain[len(craft.chain)-1]]=int(rf.read(5).split()[0]) #Picking atom no. present in the pDB file
                craft.chain_at[craft.chain[len(craft.chain)-1]]=at_count
                at_count=0
                line=rf.readline()
                word=rf.read(6).strip()
                if word == "END":
                    break
                if(word == "ATOM"):
                    rf.read(15)
                    craft.chain.append(rf.read(1))
            if(word == "HETATM"):                   #Search for hets present in chains of given pdb
                #at_no=int(rf.read(5).split()[0])
                rf.read(11)
                word1=rf.read(3).strip()
                if(word1 == "HOH"):
                    break
                else:
                    rf.read(1)
                    word2=rf.read(1)
                    if word2 in craft.het.keys():
                        if word1 in craft.het[word2]:
                            at_count+=1
                            craft.het_at[word2][lig_count]=[word1, at_count]
                            line=rf.readline()
                            word=rf.read(6).strip()
                            continue
                        else:
                            craft.het[word2].append(word1)
                            at_count=0
                            lig_count+=1
                            craft.het_at[word2].append([word1, at_count+1])
                    else:
                        craft.het[word2]=[word1]
                        craft.het_seq.append(word2)
                        craft.het_at[word2]=[[]]
                        at_count=0
                        lig_count=0
                        craft.het_at[word2][lig_count]=[word1, at_count+1]
            if(word=="ANISOU"):
                self.anisou=True
                line=rf.readline()
                word=rf.read(6).strip()
            else:
                at_count+=1
                line=rf.readline()
                word=rf.read(6).strip()
    
    '''Display the different entries present in given pdb'''
    def display_entries(self):              
        for chain in craft.chain:
            print(" Chain : ",chain)
            if chain in craft.het.keys():
                print("         having : ")
                for value in craft.het[chain]:
                    print("                ",value, end="")
                    if value in craft.het_name.keys():
                        print("      i.e. ",craft.het_name[value])
                    else:
                        print()
            else:
                print()
    
    '''Ask user to make selection according to available entries in the pdb'''
    def selection(self):                     
        craft.sel_chain=list(input("Please select chains in capital letters seprated by comma for pocket prediction : ").split(","))
        print("craft.sel_chain : ",craft.sel_chain)
        if all(x in craft.chain for x in craft.sel_chain):
            print()
            print("Make ligand selection for selected chain :: ")
            print("Please select hets in capital letters seprated by comma for pocket prediction : or avoid selection by entering \"NA\"  :: ")
            for key in craft.sel_chain:
                print()
                sel_het=[]
                if key in craft.het.keys():
                    print("Select hets for chain : ",key," : out of : ",craft.het[key])
                    sel_het=list(input(" Please select  :: ").split(","))
                if all(x in craft.het[key] for x in sel_het):
                    craft.sel_hets[key]=sel_het
                else:
                    if sel_het[0] == "NA":
                        pass
                    else:
                        print("You have made wrong selection. Please make selection as per hets present in the chain")
                        quit()
        else:
            print("Please provides correct chain in capital letters present in the listed chains ")
            quit()
    
    #Function used for multiple file at a time
    def multi_selection(self):
        craft.sel_hets={}
        craft.sel_chain=craft.chain
    
    
    '''Display the selected entries present in given pdb'''
    def display_selected(self):              
        for chain in craft.sel_chain:
            print(" Chain : ",chain)
            if chain in craft.sel_hets.keys():
                print("         having : ")
                for value in craft.sel_hets[chain]:
                    print("                ",value, end="")
                print()
            else:
                print()
    
    '''Read atom position, residues, 3D coordinates and atom type for selected entries'''
    def read_coords(self,rf):               
        rf.seek(0)
        word=rf.read(6).strip()
        while(word!= "ATOM"): 
            line=rf.readline()
            word=rf.read(6).strip()    
        '''Reading of selected chains'''         
        for chain in craft.chain:
            index=craft.chain.index(chain)           #Index use to find the location of chain in the dictionary
            if chain in craft.sel_chain:
                if index ==0:
                    at_no=1
                else:
                    index-=1
                    at_no=1
                    while index >= 0:
                        at_no+=craft.chain_at[craft.chain[index]]+1       #at_no variable find the atom no. using store data in craft.chain_at,craft.het_at
                        index-=1
                while word!="TER":
                    read_atno=int(rf.read(5).split()[0])
                    if (at_no == read_atno):                        #read_atno reads atom number from the pdb file
                        rf.read(1)
                        craft.atom_type.append(rf.read(4).strip())
                        craft.res.append(rf.read(4).strip())
                        rf.read(1)
                        craft.chain_id.append(rf.read(1))
                        craft.resn.append(rf.read(4).strip())
                        rf.read(4)#previously it was 3
                        craft.pts.append([float(rf.read(8).split()[0]),float(rf.read(8).split()[0]),float(rf.read(8).split()[0])])
                        rf.read(22)
                        craft.atom.append(rf.read(2).strip())
                        craft.charge.append(rf.read(2).strip())
                        at_no+=1
                        line=rf.readline()
                        word=rf.read(6).strip()
                    else:
                        if self.anisou== True:
                            rf.readline()
                            word=rf.read(6).strip()
                        else:
                            print("Warning: atom number ",at_no," missing from the given pdb")
                            print(at_no,"  :  ",read_atno)
                            return 0
                            #rf.close()
                            #quit()
                line=rf.readline()
                word=rf.read(6).strip()
            else:
                count=0
                while(count<=craft.chain_at[craft.chain[index]]):
                    rf.readline()
                    count+=1
                word=rf.read(6).strip()
        
        '''Reading of selected ligands'''
        at_no=1
        for chain in craft.chain:
            at_no+=craft.chain_at[chain]+1
        for chain in craft.het_seq:
            if chain in craft.sel_hets.keys():
                lig_count=0
                for lig in craft.het_at[chain]:
                    if lig[0] in craft.sel_hets[chain]:
                        for x in range(craft.het_at[chain][lig_count][1]):
                            read_atno=int(rf.read(5).strip())
                            if (at_no == read_atno):
                                rf.read(1)
                                craft.atom_type.append(rf.read(4).strip())
                                rf.read(1)
                                craft.res.append(rf.read(3).strip())
                                rf.read(1)
                                craft.chain_id.append(rf.read(1))
                                rf.read(1)
                                craft.resn.append(rf.read(4).strip())
                                rf.read(3)
                                a=[float(rf.read(8).split()[0]),float(rf.read(8).split()[0]),float(rf.read(8).split()[0])]
                                craft.pts.append(a)
                                rf.read(22)
                                craft.atom.append(rf.read(2).strip())
                                craft.charge.append(rf.read(2).strip())
                                at_no+=1
                            else:
                                print("Warning: Hetatom number ",at_no," missing from the given pdb")
                                return 0
                                #rf.close()
                                #quit()
                            line=rf.readline()
                            rf.read(6)
                    else:
                        at_no+=craft.het_at[chain][lig_count][1]
                        for x in range(craft.het_at[chain][lig_count][1]):
                            rf.readline()
                            rf.read(6)
                    lig_count+=1
            else:
                lig_count=len(craft.het_at[chain])
                for lig in range(lig_count):
                    at_no+=craft.het_at[chain][lig][1]
                    for x in range(craft.het_at[chain][lig][1]):
                        rf.readline()
                        rf.read(6)
        
        #rf.close()
        craft.resn = list(map(int, craft.resn))
        return 1
    
    '''Function to select the alternate position for the selected residues'''
    def select_position(self): 
        i=0
        alt_position=[]
        while i<len(craft.res):
            list=[]
            alt=[]
            if len(craft.res[i])>3:
                list.append(craft.resn[i])
                list.append(craft.chain_id[i])
                list.append(craft.res[i][1:4])
                
                while craft.resn[i]==int(list[0]):
                    if len(craft.res[i])>3:
                        if craft.res[i][0] not in alt:
                            alt.append(craft.res[i][0])
                    i+=1
                list.append(alt)
                alt_position.append(list)
            else:
                i+=1
        
        sel_position=[]
        print("Make a position selection for alternate postion of following residues: ")
        print("alter : ",alt_position)
        for x in alt_position:
            print(" Alternate position present for ",x, " having current occupancy ",1/len(x[3]))        
            sel_position.append(input(" Please select  :: "))
            
        i=0
        for x in alt_position:
            index=craft.resn.index(x[0])
            while craft.chain_id[index]!=x[1]:
                index=craft.resn.index(x[0],index+1)

            while craft.resn[index]==x[0]:
                if len(craft.res[index])>3:
                    if craft.res[index][0]!=sel_position[i]:
                        del craft.res[index]
                        del craft.resn[index]
                        del craft.chain_id[index]
                        del craft.atom_type[index]
                        del craft.atom[index]
                        del craft.charge[index]
                        del craft.pts[index]
                    else:
                        craft.res[index]=craft.res[index][1:4]
                        index+=1
                else:
                    index+=1
            i+=1
        
        #craft.chain_at=len(craft.res)
        
        
    
        
        
        
        
        
        
        
        
        
               