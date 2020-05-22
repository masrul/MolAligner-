import numpy as np 

class Aligner:
    
    def __init__(self,coordFile):
        self.coordFile=coordFile 
        self._natoms=0
        self._symbols=[]
        self._resids = [] 
        self._atomids = [] 
        self._resnames = [] 
        self._x=[]
        self._y=[]
        self._z=[] 
        self.read() 

    def read(self): 
        fileHandler=open(self.coordFile,"r")
        lines=fileHandler.readlines() 

        for line in lines: 
            if line.startswith('HETATM') or line.startswith('ATOM'): 
                self._atomids.append(int(line[6:11]))
                self._symbols.append(line[12:16])  
                self._resnames.append(line[17:20])
                self._resids.append(int(line[22:26]))
                self._x.append(float(line[30:38]))    
                self._y.append(float(line[38:46]))    
                self._z.append(float(line[46:54])) 

                self._natoms +=1 
        
        fileHandler.close() 


    
    def align(self,iatom= None,jatom = None,target_dir= None): 
        
        if not iatom or not jatom: 
            iatom,jatom = self.findMoleculeAxis() 

        if not target_dir: 
            target_dir = [0.0, 0.0, 1.0]

        coords=np.zeros((3,self._natoms),dtype=float)
        coords[0,:]=self._x
        coords[1,:]=self._y
        coords[2,:]=self._z  

        vec1=[
                coords[0,iatom-1] - coords[0,jatom-1],
                coords[1,iatom-1] - coords[1,jatom-1], 
                coords[2,iatom-1] - coords[2,jatom-1]
            ]
        coords=do_align(vec1,target_dir,coords)  
        
        self._x = coords[0,:]
        self._y = coords[1,:]
        self._z = coords[2,:] 
    


    
    def write(self,outfile = 'out.pdb'): 

        filehandler = open(outfile,'w') 
        

        for i in range(self._natoms): 
            filehandler.write('%-6s%5d %-4s %3s  %4d    %8.3f%8.3f%8.3f\n'%
              (
                'ATOM', 
                self._atomids[i], 
                self._symbols[i],  
                self._resnames[i],   
                self._resids[i],  
                self._x[i],self._y[i],self._z[i]
              )
            )
        
        filehandler.close() 
    
    def translate(self,transPos): 
        
        xcom = sum(self._x)/self._natoms 
        ycom = sum(self._y)/self._natoms 
        zcom = sum(self._z)/self._natoms 

        for i in range(self._natoms): 
            self._x[i] = self._x[i] - xcom + transPos[0]
            self._y[i] = self._y[i] - ycom + transPos[1]
            self._z[i] = self._z[i] - zcom + transPos[2]

    def findMoleculeAxis(self):  

        pairs = [] 
        dist = [] 

        for i in range(self._natoms-1): 
            for j in range(i+1,self._natoms): 
                dx = self._x[i] - self._x[j] 
                dy = self._y[i] - self._y[j] 
                dz = self._z[i] - self._z[j] 
                
                dist.append(dx**2 + dy**2 + dz**2) 
                pairs.append((i,j))

        argmax = np.argmax(dist) 

        return pairs[argmax][0] + 1, pairs[argmax][1] + 1 
    
    @property 
    def x(self): 
        return self._x 
    
    @property 
    def y(self): 
        return self._y 

    @property 
    def z(self): 
        return self._z  

    @property 
    def natoms(self): 
        return self._natoms   

    @property 
    def symbols(self): 
        return self._symbols  
    
    @property 
    def resnames(self): 
        return self._resnames   
    
    @property 
    def atomids(self): 
        return self._atomids 

def get_rot_matrix(v,c): 

    rotMat=np.zeros((3,3),dtype=float)
    iMat=np.zeros((3,3),dtype=float)
    vx=np.zeros((3,3),dtype=float)

    iMat[0,0]=1.0
    iMat[1,1]=1.0
    iMat[2,2]=1.0

    vx[0,0]=0.0 
    vx[0,1]=-v[2] 
    vx[0,2]=v[1]  

    vx[1,0]=v[2]
    vx[1,1]=0.0 
    vx[1,2]=-v[0]

    vx[2,0]=-v[1]
    vx[2,1]=v[0]
    vx[2,2]=0.0  

    factor=1.0/(1.0+c)

    rotMat=iMat + vx + np.matmul(vx,vx) * factor 
    
    return rotMat 

def do_align(u,v,coords): 
    
    u=u/np.linalg.norm(u)
    v=v/np.linalg.norm(v) 

    normal=np.cross(u,v)
    c=np.dot(u,v)
    
    rotMat=get_rot_matrix(normal,c)
    coords=np.matmul(rotMat,coords) 
    
    return(coords)

