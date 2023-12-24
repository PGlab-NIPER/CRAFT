#tetrahedral.py file

#This file basically run through the developed Delaunay triangulation using the concept of Flow Transfer Algorithm
#It identify the diverse types of cavities present within the proteins

from .CraftVariables import craft
import numpy as np 

class Tetrahedral():
    def __init__(self, base):
        self.vol=base.volumes
        self.rad=base.circumradius
        self.edge_coords =base.edge_coords
        self.token_hierarchy=base.token_hierarchy
        self.node_radii=np.array([craft.atom_radii[x]  for x in craft.atom])
        self.face_token = np.array([[ 2, 3,1], [ 3, 0,2], [ 0, 1,3], [1, 2,0]])
        
        self.tri = base.tri   
        self.cells=base.cells
        print("points are : ",len(craft.pts))
        print("cells length is : ",len(self.cells))
    
    # Function identify the delimiter tetrahedra
    def delimiter_tetrahedral(self,Peel,Extend_boundary=False):
        surface_tetra=np.unique(np.where(self.tri.neighbors==-1)[0])                       #Search for surface tetrahedrals present on the mesh 
        radius=self.rad[surface_tetra]
        peeled=surface_tetra[np.where(radius>Peel)[0]]
        old_neigh=np.delete(np.unique(self.tri.neighbors[peeled]),0)
        while len(old_neigh) !=0:
            old_neigh=np.array([i for i in old_neigh if i not in surface_tetra])
            radius=self.rad[old_neigh.astype(int)]
            peeled=old_neigh[np.where(radius>Peel)[0]]
            if bool(int(Extend_boundary))== False:
                surface_tetra=np.append(surface_tetra,old_neigh)                    #To add small and big peeled off tetrahedrals
            else:
                surface_tetra=np.append(surface_tetra,peeled)                       #To add only big peeled off tetrahedrals
            new_neigh=np.unique(self.tri.neighbors[peeled.astype(int)])
            if len(new_neigh) !=0 and new_neigh[0]== -1:
                new_neigh=np.delete(new_neigh,0)
            old_neigh=new_neigh
        surface_batom=np.unique(self.cells[surface_tetra.astype(int)])
        present=np.any(np.isin(self.cells,surface_batom),axis=1)                    #When tetrahedral contain any one delimiter atoms
        #present=np.all(np.isin(self.cells,surface_batom),axis=1)                   #When tetrahedral contain delimiter atoms 
        self.boundary_tetra=np.where(present==True)[0]
        #print("boundary : ",self.boundary_tetra)
    
    #Function identify the susceptible tetrahedra in the Delaunay mesh
    def susceptible_tetra(self,Rad,Vol):
        for tetra in range (len(self.cells)):
            if self.vol[tetra]>=Vol :
                if self.rad[tetra]>=Rad:
                    if tetra in self.boundary_tetra:              #Eliminating boundary tetrahedrals being part of susceptible tetrahedrals
                        continue
                    else:
                        craft.suscept_tetra.append(tetra)
    
    #Function estimate the MCR that can pass through a face
    #It also identify tetrahedra face which are forbidden or nonforbidden 
    #which is further used to flow transfer from a tetrahedron to neighbouring tetrahedron
    '''
    The radius of MCR (r) can be easily determined by solving the subsequent set of Euclidean distance equations.
    x^2+y^2=〖(r_2+r)〗^2                      			          (1)
    〖(x-a)〗^2+y^2=〖(r_3+r)〗^2         			          (2)
    〖(x-X)〗^2+〖(x-Y)〗^2=〖(r_1+r)〗^2  			          (3)

    Equation 1 & 2: Calculate x by subtracting them
    x = ((r_2^2 - r_3^2 + 2r_2*r - 2r_3*r + a^2)) / 2a
    Simplify it by separating constant terms
    x = d_1 * r + c_1             (4)
    where c_1 = ((r_2^2 - r_3^2 + a^2)) / 2a, d_1 = ((r_2 - r_3)) / a
    
    Equation 1 & 3: Calculate y by subtracting them
    y = ((r_2^2 - r_1^2 + 2r_2*r - 2r_1*r + X^2 + Y^2 - 2X*x)) / 2Y
    Simplify it by separating constant terms
    y = d_2 * r + c_2 - P * x      (5)
    where c_2 = ((r_2^2 - r_1^2 + X^2 + Y^2)) / 2Y, d_2 = ((r_2 - r_1)) / Y, P = X / Y
    
    Substitute the value of x (equation 3) in equation 2.6 to eliminate x
    y = d_2 * r - R * r - Q        (6)
    where Q = c_2 - P * c_1, R = P * d_1
    
    Solve for the MCR(r) that can tangentially pass through vertexes with distinct radius circles in the XY plane
    m = d_1^2 + d_2^2 + R^2 - 2d_2 * R - 1
    n = 2 * (d_2 * Q + d_1 * c_1 - R * Q - r_2)
    s = c_1^2 + Q^2 - r_2^2
    Solve the quadratic equation mr^2 + nr + s = 0 to calculate MCR(r)
    D = n^2 - 4 * m * s
    r = ((-n - √D)) / (2 * m)'''
    
    def faceoverlap(self,tetra, Probe=1.4):
        forbid_faces={}
        forbid_radii={}
        
        #face_token = np.array([[ 2, 3,1], [ 3, 0,2], [ 0, 1,3], [1, 2,0]]).T            #[B,C,A]

        edges3d_sum=np.sum(np.square(self.edge_coords[:,:,tetra,:]),axis=3)
        dist3d=np.sqrt(edges3d_sum)                                                     #sides [a,b,c]

        X=(edges3d_sum[2]+edges3d_sum[0]-edges3d_sum[1])/(2*dist3d[0])
        Y=np.sqrt(edges3d_sum[2]-np.square(X))
       
        edge2d=np.array([np.zeros([4, len(tetra)], dtype = int),dist3d[0],X,Y])
        
        #3 X 4 X triangle X 2   #[B,C,A]
        point2d=np.array([
                np.stack((edge2d[0], edge2d[0]), axis=2),
                np.stack((edge2d[1], edge2d[0]), axis=2),
                np.stack((edge2d[2], edge2d[3]), axis=2)
                ])
        
        nodes=self.cells[tetra].T
        faces=nodes[self.face_token.T]
 
        face_radii=self.node_radii[faces]     #[B,C,A]          #Represent n fraction of section
        face_radii2=np.square(face_radii)
        
        c1=(face_radii2[0]-face_radii2[1]+np.square(point2d[1,:,:,0]))/(2*point2d[1,:,:,0])
        d1=(face_radii[0]-face_radii[1])/point2d[1,:,:,0]

        c2=(face_radii2[0]-face_radii2[2]+np.square(point2d[2,:,:,0])+np.square(point2d[2,:,:,1]))/(2*point2d[2,:,:,1])
        d2=(face_radii[0]-face_radii[2])/point2d[2,:,:,1]
        P=point2d[2,:,:,0]/point2d[2,:,:,1]
        
        Q=c2-(P*c1)
        R=P*d1
        
        m=np.square(d1)+np.square(d2)-(2*d2*R)+(R*R)-np.ones((4, 1))
        n=2*((d2*Q)+(d1*c1)-(R*Q)-face_radii[0])
        s=np.square(Q)+np.square(c1)-face_radii2[0]
        
        D=n*n-(4*m*s)
        up_radius=(-n-np.sqrt(D))/(2*m)

        forbid_token=np.where(up_radius <=Probe)    #<=
        
        for cell,token in zip(forbid_token[1],forbid_token[0]):
            if tetra[cell] in forbid_faces.keys():
                forbid_faces[tetra[cell]].append(faces[:,token,cell].tolist())
                forbid_radii[tetra[cell]].extend([up_radius[token,cell].tolist()])
            else:
                forbid_faces[tetra[cell]]=[faces[:,token,cell].tolist()]
                forbid_radii[tetra[cell]]=[up_radius[token,cell].tolist()]   
        return forbid_faces, forbid_radii
            
    #Function to extend the forbidden faces of a cavity
    def overlap_forbidden(self,extended_atom,face): 
        face_coords=craft.pts[np.array(face)]
        pt_coords=craft.pts[extended_atom]
        edge_coords=pt_coords-face_coords
        edge_dist=np.sqrt(np.einsum("ij,ij->i",edge_coords,edge_coords))
        
        face_radii=self.node_radii[face]
        edge_radii=self.node_radii[extended_atom]+face_radii
        overlapped_atom=np.sum(edge_dist <= edge_radii)
        
        if overlapped_atom>=2:
            return True
        else:
            return False 
    
    #Function to scan cavity for all the susceptible tetrahedrals present inside the protein
    def cavity_scan(self, readpdb, Probe, Extend_forbidden_face=False, Ntetra=15):
        while len(craft.suscept_tetra)!=0:   #Scan one by one all susceptible tetrahedrals
            entire_tetras=np.array([],dtype=np.int64)
            connect_tetra=np.array([],dtype=np.int64)
            connect_face=[]
            seed=np.array([craft.suscept_tetra[0]])
            while len(seed) !=0:
                entire_tetras=np.append(entire_tetras,seed)
                boundary=np.intersect1d(seed, self.boundary_tetra)
                seed=np.setdiff1d(seed, boundary) #set unique to make fast
                
                connect_tetra=np.append(connect_tetra,seed)
                forbid_faces, forbid_radii=self.faceoverlap(seed, Probe)
                
                permitted_neigh=np.array([],dtype=np.int64)
                for tetra in seed:
                    neigh=np.unique(self.tri.neighbors[tetra])
                    neigh=np.delete(neigh, np.where(neigh == -1))
                    neigh=np.setdiff1d(neigh,entire_tetras)
                    if tetra in forbid_faces.keys():
                        forbid_token=[]
                        neigh_pts=self.cells[neigh]
                        for i, face in enumerate(forbid_faces[tetra]):
                            for j, ntetra in enumerate(neigh_pts):
                                contain =  all(elem in ntetra  for elem in face)
                                if contain:
                                    forbid_token.append(j)
                                    if (bool(int(Extend_forbidden_face))==True) and forbid_radii[tetra][i]>(Probe*0.5) :
                                        extended_atom=list(set(ntetra) ^ set(face))[0]
                                        if (self.overlap_forbidden(extended_atom, face)==True ):
                                            connect_tetra=np.append(connect_tetra,neigh[j])
                                            extended_faces=ntetra[np.sort(self.face_token)]
                                            token = np.all(extended_faces == np.sort(face), axis=1)
                                            extended_faces=extended_faces[~token]
                                            connect_face.extend(extended_faces)
                                        else:
                                            connect_face.append(np.sort(face))
                                    else:
                                        connect_face.append(np.sort(face))
                                    break
                        non_forbid_token=list(set(range(len(neigh))) ^ set(forbid_token))
                        permitted_neigh=np.append(permitted_neigh,neigh[non_forbid_token])
                    else:
                        permitted_neigh=np.append(permitted_neigh,neigh)
                seed=np.unique(permitted_neigh)
            craft.suscept_tetra = np.setdiff1d(craft.suscept_tetra,entire_tetras)
            
            cavity_satom=np.unique(self.cells[connect_tetra])
            if len(cavity_satom)>Ntetra:
                craft.cavities.append(cavity_satom.tolist())
                craft.connect_tetra.append(connect_tetra.tolist())
                connect_face=np.unique(np.array(connect_face), axis=0).tolist()
                craft.connect_face.append(connect_face)
                


    




        
        
        
        