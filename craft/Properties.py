#properties.py file

# This file characterizes the cavities identified by the Flow Transfer algorithm.
# The properties defined here aid in identifying druggable cavities among them

import os
from operator import itemgetter
from .CraftVariables import craft
import numpy as np
from itertools import *
import freesasa
import csv


class Properties():
    def __init__(self,tetrahedral,base,file,Probe,cwd):
        self.node_radii=tetrahedral.node_radii
        self.vol=base.volumes
        self.rad=base.circumradius
        self.cells=base.cells
        self.tri = base.tri
        craft.res=np.array(craft.res)                      
        craft.resn=np.array(craft.resn)                    
        craft.chain_id=np.array(craft.chain_id)            
        craft.atom_type=np.array(craft.atom_type)          
        craft.atom=np.array(craft.atom)                    
        self.sequence=base.Sequence()
        
        ids = list(range(len(craft.cavities)))
        self.run_function(ids,base,cwd,Probe)
        #Rank cavity volume list in decending order
        indices=(np.argsort(craft.cavity_prop["Volume"])[::-1]).tolist()                 
     
        #Arrange all pockets in decreasing order by Rank score 
        craft.cavity_prop= {key: [craft.cavity_prop[key][idx] for idx in indices] for key in craft.cavity_prop}        
        craft.cavity_details= {key: [craft.cavity_details[key][idx] for idx in indices] for key in craft.cavity_details}     
        self.write_Cfile(cwd,file)
        self.write_Pfile(cwd,file)
        self.displayProp()
        '''for x in range(len(ids)):
            print(x+1," : ",craft.cavity_details["Catoms"][x])
            print()'''
            
    #Function that all property functions to characterize each identified cavity
    def run_function(self,ids,base,cwd,Probe):
        ARs = self.aromaticity(ids)
        ADs = self.acceptor_donor(ids, base,cwd)
        SASAs= self.cavity_SASA(ids,Probe)
        Vols=self.Cavity_volume(ids,base)
        kds = self.Cavity_kdscore(ids)
        ENTs =self.entrance(ids)
        self.cavity_info(ids)
        craft.cavity_prop["Aromaticity"]=ARs
        craft.cavity_prop["Donor"]=ADs[0]
        craft.cavity_prop["Acceptor"]=ADs[1]
        craft.cavity_prop["DA ratio"]=ADs[2]
        craft.cavity_prop["PSASA"]=SASAs[0]
        craft.cavity_prop["NSASA"]=SASAs[1]
        craft.cavity_prop["TSASA"]=SASAs[2]
        craft.cavity_prop["Volume"]=Vols
        craft.cavity_prop["HP Kd"]=kds[0]
        craft.cavity_prop["HB Kd"]=kds[1]
        craft.cavity_prop["Kd Ratio"]=kds[2]
        craft.cavity_prop["Percent Contribution"]=kds[3]
        craft.cavity_prop["Exposure"]=ENTs[0]
        craft.cavity_prop["Enclosure"]=ENTs[1]
        craft.cavity_prop["Door"]=ENTs[2]
    
    #Function provides complete cavity surface information
    def cavity_info(self,ids):     
        Catoms,Cres,Cresn,Cchain,Catom_type=[],[],[],[],[]
        for id in ids:
            cavity=np.array(craft.cavities[id])
            Catoms.append((cavity+1).tolist())
            Cres.append(craft.res[cavity].tolist())
            Cresn.append(craft.resn[cavity].tolist())
            Cchain.append(craft.chain_id[cavity].tolist())
            Catom_type.append(craft.atom_type[cavity].tolist())
        craft.cavity_details["Catoms"]=Catoms
        craft.cavity_details["Cres"]=Cres
        craft.cavity_details["Cresn"]=Cresn
        craft.cavity_details["Cchain"]=Cchain
        craft.cavity_details["Catom_type"]=Catom_type
        
    #Function to display cavity properties
    def displayProp(self):    
        print(len(craft.cavity_prop))
        for C_id in range(len(craft.cavity_prop["Aromaticity"])):
            print("Cavity number : ",C_id+1)
            print("Cavity atoms : ",craft.cavity_details["Catoms"][C_id])
            print("Cavity residues name : ",craft.cavity_details["Cres"][C_id])
            print("Cavity residue no  : ",craft.cavity_details["Cresn"][C_id])
            print("chain : ",craft.cavity_details["Cchain"][C_id])
            print("atom type : ",craft.cavity_details["Catom_type"][C_id])
            print("number of aromatic ring : ",craft.cavity_prop["Aromaticity"][C_id])
            print("number of donor : ",craft.cavity_prop["Donor"][C_id])
            print("number of acceptor : ",craft.cavity_prop["Acceptor"][C_id])
            print("DA ratio : ",craft.cavity_prop["DA ratio"][C_id])
            print("PSASA : ",craft.cavity_prop["PSASA"][C_id])
            print("NSASA : ",craft.cavity_prop["NSASA"][C_id])
            print("TSASA : ",craft.cavity_prop["TSASA"][C_id])
            print("volume : ",craft.cavity_prop["Volume"][C_id])
            print("HP Kd : ",craft.cavity_prop["HP Kd"][C_id])
            print("HB Kd : ",craft.cavity_prop["HB Kd"][C_id])
            print("Kd Ratio : ",craft.cavity_prop["Kd Ratio"][C_id])
            #print("Percent Contribution: ",craft.cavity_prop["Percent Contribution"][C_id])
            print("Exposure: ",craft.cavity_prop["Exposure"][C_id])
            print("Enclosure: ",craft.cavity_prop["Enclosure"][C_id])
            print("Door: ",craft.cavity_prop["Door"][C_id])
            print()
        
    # Function that calculate number of aromatic ring present inside a cavity
    def aromaticity(self,ids):              
        ring_atom={'HIS':['CG','ND1','CD2','CE1','NE2'],
                    'TYR':['CG', 'CD1', 'CD2','CE1', 'CE2','CZ' ],
                    'PHE':['CG', 'CD1', 'CD2','CE1', 'CE2','CZ'],
                    'TRP':[['CG', 'CD1','NE1','CE2','CD2'],[ 'CD2','CE2','CE3','CZ2','CZ3','CH2']]
                    }
        result=[]
        for id in ids:
            count=0
            for k, g in groupby( enumerate(craft.cavities[id]), lambda x: x[1]-x[0] ) :
                atoms=list(map(itemgetter(1), g))           #Give a list of sequencial atoms from the cavity atom
                if len(atoms)>=5:
                    res_atom={}
                    r_no=craft.resn[atoms[0]]
                    res=craft.res[atoms[0]]
                    res_atom[res]=[]
                    for atom in atoms:
                        if craft.resn[atom]==r_no:            #Collect amino acid atoms for a single residue present inside the cavity
                            res_atom[res].append(atom)
                        else:
                            if len(res_atom[res])<5:            #If collected atoms are less than 5 then no aromaticity check required
                                del res_atom[res]
                            r_no=craft.resn[atom]
                            res=craft.res[atom]
                            res_atom[res]=[atom]
                    for key in list(res_atom.keys()):           
                        if (key not in ring_atom) | ((len(res_atom[key]))<5) :  #Aromaticity check will be done only for aromatic residues and having sequence number >5
                            del res_atom[key]
                        else:
                            atom_type=[]
                            charge=[]
                            c_sum=0                              #Total charge on a residue
                            Pring=0                              #Pring represent number of rings present in reference residues like 2 in TRP
                            for value in res_atom[key]:
                                atom_type.append(craft.atom_type[value])
                                if key=='HIS':
                                    charge.append(craft.charge[value])
                            for x in charge:                      # His is aromatic above Ph7
                                if x=='':
                                    continue
                                else:
                                    if len(x)==2:
                                        x=x[1]+x[0]
                                    c_sum=c_sum+int(x)                #Calculate charge for His which show aromatic behaviour when imidazole ring remain netural(basic pH)
                            if key == 'TRP':                      #One aromatic ring can participate inside the cavity
                                for r in ring_atom[key]:
                                    if all(elem in atom_type  for elem in r):       #Check presence of all atom of reference residue inside the atom_type
                                        Pring+=1
                            else:
                                if all(elem in atom_type  for elem in ring_atom[key]):
                                    Pring+=1
                            if (Pring>0 and c_sum <= 0):
                                count+=1
            #print("aromatic ring : ",count)
            result.append(count)
        return result

    #Function estimate number of donors and acceptors present inside a cavity
    def acceptor_donor(self,ids,base,cwd):           
        donors=[]
        acceptors=[]
        da_ratio=[]
        valency={'O':2,'N':3,'F':1}
        aa=["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU","MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"] #"ACE","NME"
        cf = open(cwd+"/docs/amino_connectivity.txt", "r")
        for id in ids:
            donor,acceptor=0,0
            i=0
            cf.seek(0)
            patoms,atom_charge=base.polar_connectivity(craft.cavities[id])           #Provide information about polar atoms
            for atom,res in patoms:
                if res in aa:
                    connect=0                                   #Connect find number of bond one atom can form
                    word=cf.read(3).strip()
                    while word != res:
                        word=cf.read(3).strip()
                        line=cf.readline()
                    if word == "ALA":
                        line=cf.readline()
                    word1=cf.read(4).strip()
                    while word1!="END":
                        word2=cf.read(4).strip()
                        if (word1 == atom or word2 == atom):
                            a=cf.read(1)                        #a reads the number of bond form by a selected atom with other atom from the file
                            connect+=int(a)
                        line=cf.readline()
                        word1=cf.read(4).strip()
                    cf.seek(0)
                    if atom_charge[i][1] =='':
                        charge=0
                    else:
                        str=atom_charge[i][1]
                        if len(str)==1 and type(str)==int:
                            charge=str[0]
                        elif len(str)==2:
                            charge=str[1]+str[0]
                        else:
                            charge=0
                        #charge=str[1]+str[0]                   #last change
                    nH=valency[atom_charge[i][0]]-connect+int(charge)   #nH represent no. of hydrogen atoms attached to selected atom
                    if nH == 0:
                        acceptor+=1
                    elif int(charge) >0:
                        donor+=1
                    else:
                        acceptor+=1
                        donor+=1
                    i+=1
                else:
                    continue
            donors.append(donor)
            acceptors.append(acceptor)
            if acceptor !=0:
                da_ratio.append(round(donor/acceptor,2))
            else:
                da_ratio.append("NA")
        cf.close()
        return donors,acceptors,da_ratio       
        
    #Function eastimate the SASA for identified cavities using freesasa module
    def cavity_SASA(self,ids,Probe):
        psasa,nsasa,tsasa=[],[],[]
        default_options = {'hetatm' : True, 'hydrogen' : True,'join-models' : False,'skip-unknown' : False,'halt-at-unknown' : False}
        freesasa.setVerbosity(freesasa.nowarnings)
        eatoms=['O','N','F','S']                    #eatoms is electronegative atoms
        result=freesasa.calcCoord(craft.pts.flatten(),self.node_radii,freesasa.Parameters({'probe-radius' : Probe}))
        
        '''structure = freesasa.Structure(cwd+"/PDB/"+File,options =default_options)
        structure.setRadii(self.node_radii)
        result1 = freesasa.calc(structure,freesasa.Parameters({'probe-radius' : Probe}))'''
        
        for id in ids:
            sasa1,sasa2=0,0
            for x in craft.cavities[id]:
                if(craft.atom[x] in eatoms):
                    sasa1+=result.atomArea(x)
                else:
                    sasa2+=result.atomArea(x)
            psasa.append(round(sasa1,2))
            nsasa.append(round(sasa2,2))
            tsasa.append(round(sasa1+sasa2,2))
        return psasa,nsasa,tsasa
    
    #Function calculate volume of a cavity
    def Cavity_volume(self,ids,base):    
        volumes=[]
        for id in ids:
            Cavity_vol=0
            for tetra in craft.connect_tetra[id]:
                bt_vol=self.vol[tetra] 
                nodes=self.cells[tetra]
                onode=np.tile(np.array(nodes), (3,1))                                    #onode=occupied node
                rnodes_token=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])                 #rnode=remaining nodes
                rnodes=nodes[rnodes_token]
                onode_coords=craft.pts[onode.T]
                rnodes_coords=craft.pts[rnodes]
                edge_coords=onode_coords-rnodes_coords
                edge_dist=np.sqrt(np.einsum("ijk,ijk->ij",edge_coords,edge_coords))
                m=self.node_radii[onode.T]                                               #m represent the radius of occupied nodes
                n=edge_dist-m                                                            #n represent the radius of remaining atoms
    
                m_rnodes=np.einsum("ijk,ij->ijk",rnodes_coords,m)
                n_onode=np.einsum("ijk,ij->ijk",onode_coords,n)
                intersect_coords=(m_rnodes+n_onode)/edge_dist[:,:,np.newaxis]
                intersect_coords=np.transpose(intersect_coords, (1, 0, 2))
                edge_coords=intersect_coords-intersect_coords[[1,2,0]]
                intersect_dist=np.sqrt(np.einsum("ijk,ijk->ij", edge_coords,edge_coords)).T
                circumcenters=base.circumcircle_center_3d(intersect_coords,intersect_dist)
                onode_coords=craft.pts[nodes]                                                #dimensions are reduced
                onode_radii=self.node_radii[nodes]
                Scap_coords=onode_coords-circumcenters
                Scap_dist=np.sqrt(np.einsum("ij,ij->i", Scap_coords,Scap_coords)).T
                H=onode_radii-Scap_dist
                H=np.where(H < 0, 0, H)
                Scap_vol=((np.pi*H*H)*(3*onode_radii-H))/3                                  # Spherical cap volume is calculated for onodes
                st_coords=np.concatenate((intersect_coords, onode_coords[np.newaxis, :, :]), axis=0)
                st_coords=np.transpose(st_coords, (1, 0, 2))
                st_vol=base.st_volume(st_coords)                                            # st_vol represents small tetrahedral volume
                occ_vol=np.sum(Scap_vol)+np.sum(st_vol)
                unocc_vol=bt_vol-occ_vol
                if unocc_vol>=0:
                    Cavity_vol+=unocc_vol
            volumes.append(round(Cavity_vol,2))
        return volumes
    
    #Function to estimate the kd values of a cavity
    def Cavity_kdscore(self,ids):       
        nAPA={"ALA":5,"CYS":6,"ASP":8,"GLU":9,"PHE":11,"GLY":4,"HIS":10,"ILE":8,"LYS":9,"LEU":8,"MET":8,"ASN":8,"PRO":7,"GLN":9,
        "ARG":11,"SER":6,"THR":7,"VAL":7,"TRP":14,"TYR":12}                                                             #nAPA is number of non hydrogen atoms per amino acid
        
        kd_value={"ALA":1.8,"CYS":2.5,"ASP":-3.5,"GLU":-3.5,"PHE":2.8,"GLY":-0.4,"HIS":-3.2,"ILE":4.5,"LYS":-3.9,"LEU":3.8,
        "MET":1.9,"ASN":-3.5,"PRO":-1.6,"GLN":-3.5, "ARG":-4.5,"SER":-0.8,"THR":-0.7,"VAL":4.2,"TRP":-0.9,"TYR":-1.3}   #kd_value of each amino acid
        
        kd_hp,kd_hb,kd_ratio,percent=[],[],[],[]
        for id in ids:
            phobic,philic=0,0                               # kd value for hydrophobic/hydrophilic residues
            res_atom={}
            r_no=craft.resn[craft.cavities[id][0]]                
            aa=craft.res[craft.cavities[id][0]]
            chain=craft.chain_id[craft.cavities[id][0]]
            count=0
            for atom in craft.cavities[id]:
                if craft.atom[atom]!='H':
                    if r_no == craft.resn[atom]:
                        count+=1
                    else:
                        res_atom[r_no]=[count,aa,chain]
                        r_no=craft.resn[atom]
                        aa=craft.res[atom]
                        chain=craft.chain_id[atom]
                        count=1
                        res_atom[r_no]=[count,aa,chain]
            pCon={}                                             #Pecentage contribution of non hydrogen atoms per amino acid in a cavity
            for key in res_atom.keys():
                if res_atom[key][1] in nAPA.keys():
                    ratio=res_atom[key][0]/nAPA[res_atom[key][1]]
                    string=res_atom[key][2]+res_atom[key][1]+str(key)
                    pCon[string]=round((ratio*100),2)
                    if kd_value[res_atom[key][1]]>0:
                        phobic=phobic+(ratio*kd_value[res_atom[key][1]])
                    else:
                        philic=philic+(ratio*kd_value[res_atom[key][1]])
            kd_hp.append(round(-philic,2))
            kd_hb.append(round(phobic,2))
            percent.append(pCon)
            if philic!=0:
                kd_ratio.append(round(-phobic/philic,2))
            else:
                kd_ratio.append("NA")
        return kd_hp,kd_hb,kd_ratio,percent
    
    #Function to find neghbouring faces which shares common edges
    def neighbor(self,seed,cavity_face):                
        face_neigh=[] 
        for s in seed:
            face=cavity_face[s]
            for x in range(len(cavity_face)):
                catom=np.intersect1d(face,cavity_face[x])   # catom represents common atom between to given faces
                if len(catom)==2:
                    face_neigh.append(x)       
        return face_neigh
    
    #Function to estimate the cavity entrances
    def entrance(self,ids):
        exposures,enclosures,doors=[],[],[]
        for id in ids:
            cavity_face=[]
            connect_tetra=np.array(craft.connect_tetra[id])
            connect_face=np.array(craft.connect_face[id])
            
            for tetra in connect_tetra:
                neigh=self.tri.neighbors[tetra]
                neigh=np.delete(neigh, np.where(neigh == -1))
                outer_neigh=np.setdiff1d(neigh, connect_tetra)          #outer_neigh, keep record of neighbour tetrahedral which are not part of cavity
                tetrapts=self.cells[tetra]
                for tetra2 in outer_neigh:
                    ntetrapts=self.cells[tetra2]
                    shared_face=np.intersect1d(tetrapts, ntetrapts)     # Shared face with neighbour tetrahrdeal 
                    present = np.any(np.all(shared_face == connect_face, axis=1))      # Shared face is present in cavity connect face then that face not present at cavity entrance
                    if present== False:
                        cavity_face.append(shared_face.tolist())
            cavity_face=np.array(cavity_face)
            eatom=len(np.unique(cavity_face.flatten()))
            tatom=len(craft.cavities[id])
            exposure=round(eatom/tatom,2)
            enclosure=round((tatom-eatom)/tatom,2)
            
            idx=np.array(range(len(cavity_face)),dtype=np.int64)         
            door={}
            door_no=1
            while len(idx)!=0:
                seed=np.array([idx[0]])
                door[door_no]=seed
                idx=np.delete(idx, 0)
                covered_faces=np.array([],dtype=np.int64)
                while len(seed) !=0:
                    covered_faces=np.append(covered_faces,seed)
                    neigh=np.unique(self.neighbor(seed,cavity_face)) # neighour function is used to find all connected face for a entrance
                    permitted_neigh=np.setdiff1d(neigh,covered_faces)
                    idx=np.setdiff1d(idx,permitted_neigh)
                    door[door_no]=np.append(door[door_no],permitted_neigh)
                    seed=permitted_neigh
                door_no+=1  
            exposures.append(exposure)
            enclosures.append(enclosure)
            doors.append(len(door))
        return exposures,enclosures,doors

    def normalize_array(self,arr):
        if all(element == 0 for element in arr): 
            return arr
        else:
            minimum = np.min(arr)
            maximum = np.max(arr)
            normalized_arr = (arr - minimum) / (maximum - minimum)
            return normalized_arr

    #Function to write cavity file to cmd folder
    def write_Cfile(self,cwd,file):
        wf = open(os.path.join(cwd,"Results","Cavity_files",file.split(".")[0]+".Cav"), mode='w', newline='')
        a_no=1
        for cavity in range(len(craft.cavities)):
            for atom in range(len(craft.cavity_details["Catom_type"][cavity])):
                wf.write("CATOM".ljust(6))
                wf.write(str(a_no).rjust(5))
                wf.write(" ")
                wf.write(craft.cavity_details["Catom_type"][cavity][atom].ljust(4))
                wf.write(" ")
                wf.write(craft.cavity_details["Cres"][cavity][atom].rjust(3))
                wf.write(" ")
                wf.write(craft.cavity_details["Cchain"][cavity][atom].rjust(1))
                wf.write(str(craft.cavity_details["Cresn"][cavity][atom]).rjust(4))
                wf.write(" ".rjust(4))
                point=craft.pts[craft.cavity_details["Catoms"][cavity][atom]-1]
                wf.write(str(point[0]).rjust(8))
                wf.write(str(point[1]).rjust(8))
                wf.write(str(point[2]).rjust(8))
                wf.write(" ")
                wf.write("Cavity".ljust(6))
                wf.write(str(cavity).ljust(3))
                wf.write("\n")
                a_no+=1
    
    def write_Pfile(self,cwd,file):
        wf = open(os.path.join(cwd,"Results","Property files",file.split(".")[0]+".Pro"), mode='w', newline='')
        writer = csv.writer(wf)
        header=["Cavity","Aromaticity","Donor","Acceptor","DA ratio","PSASA","NSASA","TSASA","Volume","HP Kd","HB Kd","Kd Ratio",
        "Exposure","Enclosure","Door"] 
        writer.writerow(header)
        for cavity in range(len(craft.cavities)):
            data=[]
            for x in header:
                if x =="Cavity":
                    data.append(cavity+1)
                else:
                    data.append(craft.cavity_prop[x][cavity])
            writer.writerow(data)
        

            









        
            
            