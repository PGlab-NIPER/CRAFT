# shared.py

# This file contains important shared variables used by multiple modules within the package.
# Define shared variables used across the package

class CraftVariables:
    def __init__(self):
        self.chain=[]           #Chain present in given pdb
        self.sel_chain=[]       #Selected chains by user
        self.chain_at={}        #Track no. of atoms present in each chain

        self.het= {}            #Track hets using chain as key 
        self.het_seq=[]         #Track hetro atom sequence in pdb file
        self.het_name={}        #Track hets name using het code as key 
        self.het_at={}          #Track no. of atoms present in hets of chains 
        self.sel_hets={}        #Selected hetro by user
        
        self.res=[]             #Three letter residue code
        self.chain_id=[]        #chain id of each atom present in pdb file
        self.resn=[]            #Residues number
        self.atom_type=[]       #Record atom position i.e. alpha,beta
        self.atom=[]            #Atom type i.e. C,N,O
        self.charge=[]          #Charge on each atom
        self.pts=[]             #3D coordinates of atoms

        self.atom_radii={}      #van der Waals radius of atom

        self.suscept_tetra=[]    #susceptible tetrahedral
        self.cavities=[]         #identified cavities
        self.connect_face=[]     #connected faces of cavities
        self.connect_tetra=[]    #tetrahera present in cavities

        self.cavity_prop={}      #estimated cavity properties
        self.cavity_details={}   #identified cavity surface details
        
    def reset(self):
        self.__init__()

craft=CraftVariables()    #create instance of share class


'''These variables are accessed by other modules within the src package
to perform operations related to chain selection and atom tracking.
Ensure proper handling and updating of these variables to maintain
consistency across the package.'''
