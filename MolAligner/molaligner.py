import numpy as np 
import math 

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
        
        self._x = list(coords[0,:])
        self._y = list(coords[1,:])
        self._z = list(coords[2,:])  
    


    
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
    
    def moveTo(self,transPos,atomID = None): 

        '''
        tranPos: coordinate of target point, where COM or atomID will
         be moved  
        '''
        
        if  atomID: 
            xref = self.x[atomID-1] 
            yref = self.y[atomID-1] 
            zref = self.z[atomID-1] 
        else: 
            xcom = sum(self._x)/self._natoms 
            ycom = sum(self._y)/self._natoms 
            zcom = sum(self._z)/self._natoms 
            
            xref = xcom; yref = ycom; zref = zcom 



        for i in range(self._natoms): 
            self._x[i] = self._x[i] - xref + transPos[0]
            self._y[i] = self._y[i] - yref + transPos[1]
            self._z[i] = self._z[i] - zref + transPos[2] 


    def moveBy(self, transVector): 
        
        for i in range(self.natoms): 
            self.x[i] += transVector[0] 
            self.y[i] += transVector[1] 
            self.z[i] += transVector[2] 


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
    
    
    def merge(self,other): 
        
        self._natoms += other._natoms  
        self._symbols += other._symbols 
        self._resids  += [resid+self._resids[-1]  for resid  in other._resids]  
        self._atomids += [atomid+self._atomids[-1]  for atomid  in other._atomids] 
        self._resnames += other._resnames 
        self._x += other._x 
        self._y += other._y 
        self._z += other._z 
    
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

def getAlignMatrix(v,c): 

    '''
    https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d 
    '''

    alignMat=np.zeros((3,3),dtype=float)
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

    alignMat=iMat + vx + np.matmul(vx,vx) * factor 
    
    return alignMat 

def getRotMatrix(rotAxis,angle): 
    rotMat=np.zeros((3,3),dtype=float)
    
    angle = angle * 0.0174533  # degree to radian 

    cos = math.cosine(angle)
    sin = math.sine(anlge) 
    _cos_ = 1 - cos
    
    u = rotAxis / np.linalg.norm(u) 
    
    ux = u[0]; uy = u[1]; uz = u[2] 

    uxy = ux * uy 
    uyz = uy * uz 
    uzx = uz * ux 
    uxx = ux * ux 
    uyy = uy * uy
    uzz = uz * uz 

    rotMat[0,0] = cos + (uxx * _cos_) 
    rotMat[0,1] = (uxy * _cos_) - (uz * sin) 
    rotMat[0,2] = (uzx * _cos_) + (uy * sin)

    rotMat[1,0] = (uxy * _cos_) + (uz * sin) 
    rotMat[1,1] = cos + (uyy * _cos_) 
    rotMat[1,2] = (uyz * _cos_) - (ux * sin) 

    rotMat[2,0] = (uzx * _cos_) - (uy * sin) 
    rotMat[2,1] = (uyz * _cos_) + (ux * sin)  
    rotMat[2,2] = cos + (uzz*_cos_) 

    
    return rotMat 



def do_align(u,v,coords): 
    
    u=u/np.linalg.norm(u)
    v=v/np.linalg.norm(v) 

    normal=np.cross(u,v)
    c=np.dot(u,v) 

    if abs(abs(c)-1)<10e-10:  # if vectors are antiparallel 
        coords= coords * -1 
        return coords 

    
    alignMat=getAlignMatrix(normal,c)
    coords=np.matmul(alignMat,coords) 
    
    return(coords)


# m0 =Aligner('hemcel.pdb') 
# m0.align(5,294,[0,0,1]) 
# m0.moveTo([10,0,0]) 
# # molecule0.write('out0.pdb') 

# m1 = Aligner('cellulose_2.pdb') 
# m1.align(2,246,[0,0,1]) 
# m1.moveTo([0,0,0]) 

# m0.merge(m1)
# m0.write('out-merged.pdb') 

