#base.py file

#This module handles the building of Delaunay triangulation, estimate volume, circumradius, and
# circumcenter of the delaunay tetrahera. It provide fundamental functions for the tetrahedral, and properties module

from scipy.spatial import Delaunay
from .CraftVariables import craft
import numpy as np
import time

class Base():
    def __init__(self,cwd):
        craft.pts=np.array(craft.pts)
        self.build_triangulation()              #Build delaunay triangulation for given point set
        self.cells_volume()                     #Calling function to calculate volume of delaunay cells
        self.cells_circumcenter()               #Calling function to calculate circumcenter of Delaunay cells
        self.cells_circumradius()               #Calling function to calculate circumradius of delaunay cells
        
        tatom=list(set(craft.atom))            #Total type of atom present in pdb
        try:                                    #Reading vander Waal radius for all atom type present in the pdb             
            f = open(cwd+"/docs/atomic_radius_freesasa.txt", "r")
        except IOError:
            print("Unable to locate atomic radius file ")
        else:
            ratom=f.read(2).strip()
            while True:
                for x in tatom:
                    if ratom.casefold() == x.casefold():
                        f.read(1)
                        craft.atom_radii[x]=float(f.read(4).strip())
                line = f.readline()
                if not line:
                    break
                ratom=f.read(2).strip()
            if len(craft.atom_radii) == len(tatom):
                pass
            else:
                print("Warning: vander Waals radius not known for atoms : ",[i for i in tatom if i not in craft.atom_radii.keys()])
                quit()
    
    #Function to create Delaunay trianguations for selected points using scipy package
    def build_triangulation(self):
        #Create Delaunay trianguations for selected points using scipy package
        self.tri = Delaunay(craft.pts)                     
        self.cells=self.tri.simplices
        self.cells = np.sort(self.cells, axis=1)            #Sort the cells according to the rows
        present = np.zeros(len(craft.pts), dtype=bool) 
        present[self.cells] = True
        assert np.all(present), "There are {}  nodes which are not present in the mesh".format(np.sum(~present))
        
        #cells = {"nodes": self.cells}
        nodes = self.cells.T                                #nodes (4Xlen(self.cells)) represent four nodes of tetrahedral present inside the cells
        
        #token (2X3X4) arrange (edge->face->cells) of a tetraheral
        self.token = np.array(
            [
                [[2, 3], [3, 1], [1, 2]],
                [[3, 0], [0, 2], [2, 3]],
                [[0, 1], [1, 3], [3, 0]],
                [[1, 2], [2, 0], [0, 1]],
            ]
        ).T
        

        #token_hierarchy (2X3X4X cells) arrange (edge->face->tetra->cells) of given point set
        #Edge m is oppsite of node m in each face
        self.token_hierarchy = nodes[self.token]
        
        
        #edge_coords (3X4 X cells X 3) arrange (>face->tetra->cells->coordinates) of given point set
        self.edge_coords = (
            craft.pts[self.token_hierarchy[1]]
            - craft.pts[self.token_hierarchy[0]]
        )
        
        #ei_dot_ej represent Einstein's summation of edge_coords, using substring "ijkl, ijkl->ijk" 
        self.ei_dot_ej = np.einsum(
            "ijkl, ijkl->ijk",
            self.edge_coords[[1, 2, 0]],
            self.edge_coords[[2, 0, 1]]
        )
        
        # Zeta and alpha are calcualted for each cells of given point set
        eij = self.ei_dot_ej            
        self.zeta = (
            -eij[2, [1, 2, 3, 0]] * eij[1] * eij[2]
            - eij[1, [2, 3, 0, 1]] * eij[2] * eij[0]
            - eij[0, [3, 0, 1, 2]] * eij[0] * eij[1]
            + eij[0] * eij[1] * eij[2]
        )
        
        
    #Function to calculate cell volume   
    def cells_volume(self):                                 
        self.volumes = np.sqrt(np.sum(self.zeta, axis=0) / 72.0)
        
    #Function to calculate cell circumcenter    
    def cells_circumcenter(self):
        alpha = self.zeta / np.sum(self.zeta, axis=0)
        self.circumcenters = np.sum(
            alpha[None].T * craft.pts[self.cells], axis=1
        )
    
    #Function to calculate cell circumradius
    def cells_circumradius(self):        
        dist = craft.pts[self.token_hierarchy[0, 0, 0]] - self.circumcenters
        self.circumradius = np.sqrt(np.einsum("ij,ij->i", dist, dist))
    
    #Function calculate circumcenter of the given points
    def circumcircle_center_3d(self,points,dist):
        def coeff(dist):
            dist2=np.square(dist)
            c=(dist2[:,1]*(dist2[:,2]+dist2[:,0]-dist2[:,1]))/((np.square(dist[:,2]+dist[:,0])-dist2[:,1])*(dist2[:,1]-np.square(dist[:,2]-dist[:,0])))
            return c
        coeff_value=np.vstack((coeff(dist),coeff(dist[:, [1, 2, 0]]),coeff(dist[:, [2, 0,1]])))
        point_coeff=coeff_value[:, :, np.newaxis] * points
        circumcenters=point_coeff[0]+point_coeff[1]+point_coeff[2]
        return circumcenters
    
    #Function to calculate circumcenter of given tetrahedral points)
    def st_volume(self,node_coords):             
        alpha_mat =np.concatenate((node_coords, np.ones((node_coords.shape[0], node_coords.shape[1], 1))),axis=2)
        #alpha_mat = np.hstack((node_coords,np.ones((4,1))))
        self.alpha=np.linalg.det(alpha_mat)
        return abs(self.alpha/6)
          
    #Function to find sequence for given PDB     
    def Sequence(self):             
        code={"ALA" :'A',"CYS" :'C',"ASP" : 'D',"GLU" :'E',"PHE" :'F',"GLY" :'G',"HIS" :'H',"ILE" :'I',"LYS" :'K',"LEU" :'L',
        "MET" :'M',"ASN" :'N',"PRO" :'P',"GLN" :'Q',"ARG" : 'R' ,"SER" :'S',"THR" :'T',"VAL" :'V',"TRP" :'W',"TYR" :'Y'}
        seq={}
        if len(craft.chain_id)==len(craft.resn):
            i=0
            while i < len(craft.resn):
                if craft.chain_id[i] in seq:
                    if craft.resn[i-1]==craft.resn[i]:
                        i+=1
                        continue
                    elif (craft.resn[i-1]+1)==craft.resn[i]:
                        if craft.res[i] in code:
                            seq[craft.chain_id[i]].append(code[craft.res[i]])
                    else:
                        if craft.res[i] in code:
                            gap=craft.resn[i]-craft.resn[i-1]-1
                            for x in range(gap): 
                                seq[craft.chain_id[i]].append('-')
                            seq[craft.chain_id[i]].append(code[craft.res[i]])
                else:                                                           #last change 
                    if len(craft.res[i])==4:
                        word=craft.res[i]
                        seq[craft.chain_id[i]]=[code[word[1:]]]
                    elif len(craft.res[i])==3:
                        seq[craft.chain_id[i]]=[code[craft.res[i]]]
                i+=1
        else:
            print(" Error in file reading ")
        
        return seq
    
    #Function to calculate volume of same tetrahedral points which were given to previous circumradius function
    def volume(self,node_coords):               
        alpha_mat = np.hstack((node_coords,np.ones((4,1))))
        alpha=np.linalg.det(alpha_mat)
        return abs(alpha/6)
    
    #Function to find polar atoms
    def polar_connectivity(self,Catom):
        patoms=[]
        atom_charge=[]              #Contain atom  and charge
        eatoms=['O','N','F']        #Eatoms contains electronegative atoms
        for x in Catom:
            if(craft.atom[x] in eatoms):
                patoms.append([craft.atom_type[x],craft.res[x]])
                atom_charge.append([craft.atom[x],craft.charge[x]])
        return patoms,atom_charge
        
                
                
                
                
                
                
                
                
                
                
                



 
                
                
                
                
                
    