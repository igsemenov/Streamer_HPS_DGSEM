###############################################################################
" pylib_dgs "
###############################################################################
import numpy as np
import pickle
###############################################################################
class dgs_box1d:
#------------------------------------------------------------------------------
    def __init__(self,nod):
#------------------------------------------------------------------------------

        n=nod.n
        m=nod.m

        self.n=n
        self.m=m

        self.h=nod.h
        self.g=1./nod.h

        self.f=[np.zeros_like(nod.x),np.zeros_like(nod.x)]

        self.ax=np.zeros_like(nod.x)
        self.bx=np.zeros_like(nod.x)

        self.fx=np.zeros_like(nod.x)
        self.gx=np.zeros_like(nod.x)

        self.df=np.zeros_like(nod.x)

        self.a=[np.zeros((1,m)),np.zeros((1,m))]
        self.b=[np.zeros((1,m)),np.zeros((1,m))]

        self.H=[np.zeros((1,m)),np.zeros((1,m))]
        self.F=[np.zeros((1,m)),np.zeros((1,m))]
        self.G=[np.zeros((1,m)),np.zeros((1,m))]

        self.B=[np.zeros_like(nod.x),np.zeros_like(nod.x)]

#------------------------------------------------------------------------------
    def delta(self,p,t,dgs):
#------------------------------------------------------------------------------

        f=self.f[p]

#------------------------------------------------------------------------------
#   delta: surface fluxes
#------------------------------------------------------------------------------

        np.matmul(dgs.P1,self.ax,out=self.a[0])
        np.matmul(dgs.P2,self.ax,out=self.a[1])

        np.matmul(dgs.P1,self.bx,out=self.b[0])
        np.matmul(dgs.P2,self.bx,out=self.b[1])

        np.matmul(dgs.P1,f,out=self.H[0])
        np.matmul(dgs.P2,f,out=self.H[1])

        self.F[0]=self.H[0].copy()
        self.F[1]=self.H[1].copy()

        np.multiply(self.a[0],self.F[0],out=self.F[0])
        np.multiply(self.a[1],self.F[1],out=self.F[1])

        self.F[0][self.a[0]>0.]=0.
        self.F[1][self.a[1]<0.]=0.

        self.F[1][0,:-1]=self.F[1][0,:-1]+self.F[0][0,1:]
        self.F[0][0,1:]=self.F[1][0,:-1]

        self.F[0][0,0]=self.F[0][0,0]+self.F[1][0,-1]
        self.F[1][0,-1]=self.F[0][0,0]

        np.matmul(dgs.Q2,f,out=self.G[0])
        np.matmul(dgs.Q1,f,out=self.G[1])

        np.multiply(self.b[0],self.G[0],out=self.G[0])
        np.multiply(self.b[1],self.G[1],out=self.G[1])

        self.G[1][0,:-1]=self.G[1][0,:-1]+self.G[0][0,1:]
        self.G[0][0,1:]=self.G[1][0,:-1]

        self.G[0][0,0]=self.G[0][0,0]+self.G[1][0,-1]
        self.G[1][0,-1]=self.G[0][0,0]

        np.subtract(self.F[0],self.G[0],out=self.F[0])
        np.subtract(self.F[1],self.G[1],out=self.F[1])

#------------------------------------------------------------------------------
#   delta: volume fluxes
#------------------------------------------------------------------------------

        np.multiply(self.ax,f,out=self.fx)
        np.multiply(self.bx,f,out=self.gx)

        np.matmul(dgs.C,self.fx,out=self.fx)

        np.matmul(dgs.U,self.gx,out=self.gx)
        np.matmul(dgs.C,self.gx,out=self.gx)

        np.subtract(self.fx,self.gx,out=self.df)

#------------------------------------------------------------------------------
#   delta: delta
#------------------------------------------------------------------------------

        np.matmul(dgs.F1,self.F[0],out=self.B[0])
        np.matmul(dgs.F2,self.F[1],out=self.B[1])

        np.subtract(self.B[0],self.B[1],out=self.B[0])

        np.add(self.B[0],self.df,out=self.df)

        np.multiply(t*self.g,self.df,out=self.df)    

#------------------------------------------------------------------------------
    def get_error(self,err,lgn):
#------------------------------------------------------------------------------

        L1=np.amax(abs(err))

        L2=np.matmul(lgn.C2,err)
        L2=np.sqrt(np.sum(L2*L2)*self.h)

        return L1,L2

#------------------------------------------------------------------------------
    def get_order(self,err):
#------------------------------------------------------------------------------

        k=np.sort(err.keys())

        odr=np.zeros(len(k)-1)

        for i in range(len(k)-1):
            odr[i]=np.log2(err[k[i]]/err[k[i+1]])

        return odr

#------------------------------------------------------------------------------
    def save_err(self,err1,err2,n,name):
#------------------------------------------------------------------------------

        pickle.dump(err1,open(name+"_n="+str(n)+"_L1.bin","wb"))
        pickle.dump(err2,open(name+"_n="+str(n)+"_L2.bin","wb"))

###############################################################################
