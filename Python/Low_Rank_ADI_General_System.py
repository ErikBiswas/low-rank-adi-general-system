
import time
import math
import scipy as sp
import scipy.io as sio
import numpy as np
from numpy import zeros, diag, empty,ones
from scipy import sparse
from scipy.sparse.linalg import norm as snorm
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy.sparse import diags
from scipy.sparse import spdiags
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from numpy.lib.scimath import sqrt as csqrt



class lradi:
 
    def lradi_cg(self):
        
        myF2=sio.loadmat('C:\\Users\\Erik\\Documents\\MATLAB\\CDPlayer.mat')
        #myF2=sio.loadmat('D:\\STUDY MATERIAL\\LR_ADI General System\\NumPy\\beamm.mat')
        #myF2=sio.loadmat('D:\\STUDY MATERIAL\\LR_ADI General System\\NumPy\\eady.mat')
        #myF2=sio.loadmat('D:\\STUDY MATERIAL\\LR_ADI General System\\NumPy\\build.mat')
        
        A=myF2['A']
        B=myF2['B']
        C=myF2['C']
        
        A = sparse.csr_matrix(A)
        B = sparse.csr_matrix(B)
        C = sparse.csr_matrix(C)
        E = sparse.eye(A.shape[0])
        
        maxiter1 = 200
        restol = math.pow(10,-8)
        l=10
        
        Zc,res1=self.adaptiveshift_adi_c(E,A,B,l,maxiter1,restol)
        print("\nShape of Zc = ", Zc.shape)
        z = C * Zc
        zt = csc_matrix(z.T) * csc_matrix(z)   
        h2norm = math.sqrt(np.trace(csc_matrix(zt).toarray()))
        print()
        print("H2NORM :: ",h2norm)
        
        
               
        #plt.title('H2norm')
        plt.xlabel('iteration')
        plt.ylabel('residual norm')
        plt.loglog(res1,'o-')
        plt.savefig('build')
        plt.show()
        
        return Zc,res1
        
    def lp_s(self,p1,p2,sset):
        
        max_rrr = -1
        ind = 0
        
        for k in range(len(sset)):
            x = sset[k] 
            rr=1
            
            if p2 is 0: 
                p1 = np.array(p1)
                x = np.array(x)
                rr= rr * np.abs(p1-x)/np.abs(p1+x)
                
            else:
                for kk in range(len(p2)):
                    rr= rr * np.abs(p2[kk]-x)/np.abs(p2[kk]+x)
                
            if rr > max_rrr:
                max_rrr = rr
                ind =k
                
        return round(max_rrr,2),ind               
             
    def lp_minimax(self,rw2,l2):
        
        if len(rw2)<10:
            print("Length of (rw) must be at least 10.")
        
        max_rr = +float('Inf')
        conjP = np.array([])
        for i in range(len(rw2)):
            max_r,idn = self.lp_s(rw2[i],0,rw2)
            if max_r < max_rr:
                p0 = rw2[i]
                max_rr = max_r
            #endp0
        #end
        
        if p0.imag == 0 or p0.imag == 0.0:
            conjP = np.append(conjP,p0.real)
        else:
            conjP=np.append(p0,np.conjugate(p0))
                     
        #end
        max_r,idn = self.lp_s(0,conjP,rw2)
        
        while len(conjP) < 10:
            p0 = rw2[idn]
            
            if p0.imag == 0 or p0.imag == 0.0:
                conjP = np.append(conjP,p0.real)
                
            else:
                conjP = np.append(conjP,p0)
                conjP=np.append(conjP,np.conjugate(p0))
                
            #end
            max_r,idn = self.lp_s(0,conjP,rw2)
        #end
        return conjP  
        
    def adaptive_ADIshift(self,E,A,V,l):    

        V =csc_matrix(V)
        An = V.T * A * V
        En = V.T * E * V
              
        fullAn = sparse.csr_matrix(An).toarray()
        fullEn = sparse.csr_matrix(En).toarray()
        
        solve = np.linalg.solve(fullEn,fullAn)
        
        rw0 = np.linalg.eigvals(solve)
        
        rw = np.array([])
        for j in range(rw0.shape[0]):
                       
            if rw0[j].real <0:
                rw = np.append(rw,rw0[j])
        
        p = self.lp_minimax(rw,l)
        
        
        return p
    def adaptiveshift_adi_c(self,E,A,B,r,maxiter1,restol):
        
        i=0;
        ip=-1       
        m3=B.shape[1] 
        #res1=zeros((1,maxiter1+1))
        res1 = np.array([])
        Zc=np.zeros(B.shape)
        
        m1 = E.shape
        m1 = m1[0]
        W = B
                                      
        bnorm = sp.sparse.linalg.norm(W.T*W,'fro')
        
        Bn= sp.sparse.random(m1,100,.8)
        ful = sparse.csr_matrix(Bn).toarray()
        Q, R = np.linalg.qr(ful)
        Q = sparse.csc_matrix(Q)
        ps = self.adaptive_ADIshift(E,A,Q,r)
        l=len(ps)
        re = None
              
        while i < maxiter1:
            if ip < l-1:
                ip = ip+1
            else:
                s3=Zc.shape
                m2=s3[1]
                Zcn = Zc[:,m2-l*m3:]
                
                Q, R = np.linalg.qr(Zcn)
                ps = self.adaptive_ADIshift(E,A,Q,r)
                
                l=len(ps)
                ip=0   
 
            pc=ps[ip]
            if pc.imag == 0:
                pc=pc.real
                
            vir = E*pc
            Atil = vir + A      
            Btil = csc_matrix(W)           
            Xtil = sp.sparse.linalg.spsolve(csc_matrix(Atil),csc_matrix(Btil))        
            Xc = Xtil            
            Vc=Xc           
            con = np.isreal(pc)
            
            if con == True:                
                Vc = Vc.reshape(-1,1) #Is is used to convert a row vector to a column vector
                k = csqrt((-2*pc.real)) * Vc.real             
                k1 = csc_matrix(k).toarray()                
                Zc = np.column_stack((Zc,k1))
                W = W - 2 * pc.real * E * Vc
                W = csc_matrix(W)
                 
            else:               
                beta = pc.real / pc.imag
                gam = 2* np.sqrt(-pc.real)               
                vc = Vc.reshape(-1,1)#.toarray()               
                be = gam*(vc.real+beta*vc.imag)               
                Zc = np.column_stack((Zc,csc_matrix(be).toarray()))
                ga = gam * np.sqrt((beta**2)+1) * vc.imag
                Zc = np.column_stack((Zc,csc_matrix(ga).toarray()))
                # Vc = sparse.csc_matrix.multiply(Vc.real,Vc.imag)              
                #vc = Vc.toarray()
                vc = vc.real + beta * vc.imag
                # Vc = csc_matrix(Vc).toarray()+beta               
                pc = np.conjugate(pc) 
                Mul = E*vc
                fulM=csr_matrix(Mul).toarray()*pc.real
                #fulm2 = sparse.csc_matrix.multiply(fulMul,pc.real)
                W = csr_matrix(W).toarray() - 4 * fulM
                i = i+1
                ip=ip+1
                W = csc_matrix(W)
                
                re =sp.sparse.linalg.norm(W.T*W,'fro') / bnorm
            
                print("Step: ",i + 1," Normalized Residual: ",re)
            
                old = res1
                res1 = np.append(old, re)
            
                if re<restol:
                    break
           
            i=i+1
        
        return Zc, res1
        
       
        
"""MAIN RUN CODE"""   


Zm = lradi()
Zc,res1 = Zm.lradi_cg()


print('\n#####################################################################')
print("\n*******************Program Terminated Successfully*******************\n")  
print('#####################################################################') 
"""
elapsed_time = time.clock(code_to_test, number=1000)/1000
print( elapsed_time)
time.clock()     
Zm = lradi()
Zc,res1 = Zm.lradi_cg()
print('\n#####################################################################')
print("\n*******************Program Terminated Successfully*******************\n")  
print('#####################################################################')  
"""       
