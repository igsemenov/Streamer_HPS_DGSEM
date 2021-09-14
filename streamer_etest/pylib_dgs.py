###############################################################################
" pylib_dgs_lgn "
###############################################################################
import numpy as np
###############################################################################
class dgs_edg:
#------------------------------------------------------------------------------
    def __init__(self,n,m):
#------------------------------------------------------------------------------

        self.m=m

        self.f=np.zeros((n*n,m))

        self.a=np.zeros((n,m))
        self.b=np.zeros((n,m))

        self.F=np.zeros((n,m))
        self.G=np.zeros((n,m))

    def refine(self,edg,dgs,p1,p2):

        np.matmul(dgs.T[p1],self.f,out=edg.f[:,0::2])
        np.matmul(dgs.T[p2],self.f,out=edg.f[:,1::2])

        np.matmul(dgs.T1,self.a,out=edg.a[:,0::2])
        np.matmul(dgs.T2,self.a,out=edg.a[:,1::2])

        np.matmul(dgs.T1,self.b,out=edg.b[:,0::2])
        np.matmul(dgs.T2,self.b,out=edg.b[:,1::2])

        np.multiply(2.,edg.b,out=edg.b)

    def split(self,edg1,edg2):

        m=self.f.shape[1]/2

        edg1.f[:,:]=self.f[:,0:m]
        edg2.f[:,:]=self.f[:,m::]

        edg1.a[:,:]=self.a[:,0:m]
        edg2.a[:,:]=self.a[:,m::]

        edg1.b[:,:]=self.b[:,0:m]
        edg2.b[:,:]=self.b[:,m::]

    def merge(self,edg1,edg2):

        m=self.F.shape[1]/2

        self.F[:,0:m]=edg1.F[:,:]
        self.F[:,m::]=edg2.F[:,:]

    def flux(self,edg,dgs):

        np.matmul(dgs.K1,edg.F[:,0::2],out=edg.F[:,0::2])
        np.matmul(dgs.K2,edg.F[:,1::2],out=edg.F[:,1::2])

        np.add(edg.F[:,0::2],edg.F[:,1::2],out=self.F)

###############################################################################
class dgs_rib:
#------------------------------------------------------------------------------
    def __init__(self,box,ind):
#------------------------------------------------------------------------------

        n=box[0].n

        self.i=ind

        self.e=None

        self.r1=None
        self.r2=None

        self.r1_=None
        self.r2_=None

        k=[[2,0],[1,3],[0,2],[2,0],[3,1],[1,3]]

        if ind>1:

            m=box[2].m[0]

            k1=k[ind][0]
            k2=k[ind][1]

            self.r1_=dgs_edg(n,m/2)
            self.r2_=dgs_edg(n,m/2)

            self.r1=dgs_edg(n,m)
            self.r2=dgs_edg(n,m)

            self.e=[box[0].e[k1],box[1].e[k1],box[2].e[k2]]

        else:

            k1=k[ind][0]
            k2=k[ind][1]

            self.e=[box[0].e[k1],box[1].e[k2]]

#------------------------------------------------------------------------------
    def flux(self,dgs):
#------------------------------------------------------------------------------

        if self.i==0:
            self.flxh(self.e[0],self.e[1],dgs)

        if self.i==1:
            self.flxv(self.e[0],self.e[1],dgs)

        if self.i==2:

            self.e[2].split(self.r1_,self.r2_)

            self.r1_.refine(self.r1,dgs,1,3)
            self.r2_.refine(self.r2,dgs,1,3)

            self.flxh(self.r1,self.e[0],dgs)
            self.flxh(self.r2,self.e[1],dgs)

            self.r1_.flux(self.r1,dgs)
            self.r2_.flux(self.r2,dgs)

            self.e[2].merge(self.r1_,self.r2_)

        if self.i==3:

            self.e[2].split(self.r1_,self.r2_)

            self.r1_.refine(self.r1,dgs,0,2)
            self.r2_.refine(self.r2,dgs,0,2)

            self.flxh(self.e[0],self.r1,dgs)
            self.flxh(self.e[1],self.r2,dgs)

            self.r1_.flux(self.r1,dgs)
            self.r2_.flux(self.r2,dgs)

            self.e[2].merge(self.r1_,self.r2_)

        if self.i==4:

            self.e[2].split(self.r1_,self.r2_)

            self.r1_.refine(self.r1,dgs,2,3)
            self.r2_.refine(self.r2,dgs,2,3)

            self.flxv(self.r1,self.e[0],dgs)
            self.flxv(self.r2,self.e[1],dgs)

            self.r1_.flux(self.r1,dgs)
            self.r2_.flux(self.r2,dgs)

            self.e[2].merge(self.r1_,self.r2_)

        if self.i==5:

            self.e[2].split(self.r1_,self.r2_)

            self.r1_.refine(self.r1,dgs,0,1)
            self.r2_.refine(self.r2,dgs,0,1)

            self.flxv(self.e[0],self.r1,dgs)
            self.flxv(self.e[1],self.r2,dgs)

            self.r1_.flux(self.r1,dgs)
            self.r2_.flux(self.r2,dgs)

            self.e[2].merge(self.r1_,self.r2_)

#------------------------------------------------------------------------------
    def flxv(self,e1,e2,dgs):
#------------------------------------------------------------------------------

        e1.a=0.5*(e1.a+e2.a)
        e2.a[:]=e1.a[:]

        e1.b=0.5*(e1.b+e2.b)
        e2.b[:]=e1.b[:]

        np.matmul(dgs.P2x,e1.f,out=e1.F)
        np.matmul(dgs.P1x,e2.f,out=e2.F)

        np.multiply(e1.a,e1.F,out=e1.F)
        np.multiply(e2.a,e2.F,out=e2.F)

        e1.F[e1.a<0.]=0.0
        e2.F[e2.a>0.]=0.0

        np.matmul(dgs.Q1x,e1.f,out=e1.G)
        np.matmul(dgs.Q2x,e2.f,out=e2.G)

        np.multiply(e1.b,e1.G,out=e1.G)
        np.multiply(e2.b,e2.G,out=e2.G)

        np.subtract(e1.F,e1.G,out=e1.F)
        np.subtract(e2.F,e2.G,out=e2.F)

        np.add(e1.F,e2.F,out=e1.F)

        e2.F[:]=e1.F[:]

#------------------------------------------------------------------------------
    def flxh(self,e1,e2,dgs):
#------------------------------------------------------------------------------

        e1.a=0.5*(e1.a+e2.a)
        e2.a[:]=e1.a[:]

        e1.b=0.5*(e1.b+e2.b)
        e2.b[:]=e1.b[:]

        np.matmul(dgs.P2y,e1.f,out=e1.F)
        np.matmul(dgs.P1y,e2.f,out=e2.F)

        np.multiply(e1.a,e1.F,out=e1.F)
        np.multiply(e2.a,e2.F,out=e2.F)

        e1.F[e1.a<0.]=0.
        e2.F[e2.a>0.]=0.

        np.matmul(dgs.Q1y,e1.f,out=e1.G)
        np.matmul(dgs.Q2y,e2.f,out=e2.G)

        np.multiply(e1.b,e1.G,out=e1.G)
        np.multiply(e2.b,e2.G,out=e2.G)

        np.subtract(e1.F,e1.G,out=e1.F)
        np.subtract(e2.F,e2.G,out=e2.F)

        np.add(e1.F,e2.F,out=e1.F)

        e2.F[:]=e1.F[:]

###############################################################################
class dgs_box2d:
#------------------------------------------------------------------------------
    def __init__(self,nod):
#------------------------------------------------------------------------------

        n=nod.n
        m=nod.m

        k=m[0]*m[1]

        self.n=nod.n
        self.m=nod.m
        self.h=nod.h

#------------------------------------------------------------------------------

        ij=np.arange(k,dtype=int).reshape(m[0],m[1])

        self.ie=[ij[0,:].flatten(),ij[:,-1].flatten(),
                 ij[-1,:].flatten(),ij[:,0].flatten()]

#------------------------------------------------------------------------------

        self.g=(1./nod.h[0],1./nod.h[1])

        self.f=[np.zeros_like(nod.x),np.zeros_like(nod.x),
                np.zeros_like(nod.x),np.zeros_like(nod.x),
                np.zeros_like(nod.x),np.zeros_like(nod.x),
                np.zeros_like(nod.x),np.zeros_like(nod.x)]

        self.df=[np.zeros_like(nod.x),np.zeros_like(nod.x)]

        self.ax=np.zeros_like(nod.x)
        self.ay=np.zeros_like(nod.y)

        self.bx=np.zeros_like(nod.x)
        self.by=np.zeros_like(nod.y)

        self.fx=np.zeros_like(nod.x)
        self.fy=np.zeros_like(nod.y)

        self.gx=np.zeros_like(nod.x)
        self.gy=np.zeros_like(nod.y)

        self.jx=np.zeros_like(nod.x)
        self.sx=np.zeros_like(nod.x)

        self.jx=1./nod.x

        self.a=[np.zeros((n,k)) for i in range(4)]
        self.b=[np.zeros((n,k)) for i in range(4)]

        self.F=[np.zeros((n,k)) for i in range(4)]
        self.G=[np.zeros((n,k)) for i in range(4)]

        self.B=[np.zeros_like(nod.x) for i in range(4)]

        self.e=[dgs_edg(n,m[1]),dgs_edg(n,m[0]),
                dgs_edg(n,m[1]),dgs_edg(n,m[0])]

        self.bx.fill(0.02*self.g[0])
        self.by.fill(0.02*self.g[1])

#------------------------------------------------------------------------------
    def  edge(self,p,dgs):
#------------------------------------------------------------------------------

        f=self.f[p]

        self.e[0].f[:,:]=f[:,self.ie[0]]
        self.e[2].f[:,:]=f[:,self.ie[2]]

        self.e[1].f[:,:]=f[:,self.ie[1]]
        self.e[3].f[:,:]=f[:,self.ie[3]]

        np.matmul(dgs.P2x,self.ax[:,self.ie[1]],out=self.e[1].a)
        np.matmul(dgs.P1x,self.ax[:,self.ie[3]],out=self.e[3].a)

        np.matmul(dgs.P2x,self.bx[:,self.ie[1]],out=self.e[1].b)
        np.matmul(dgs.P1x,self.bx[:,self.ie[3]],out=self.e[3].b)

        np.matmul(dgs.P1y,self.ay[:,self.ie[0]],out=self.e[0].a)
        np.matmul(dgs.P2y,self.ay[:,self.ie[2]],out=self.e[2].a)

        np.matmul(dgs.P1y,self.by[:,self.ie[0]],out=self.e[0].b)
        np.matmul(dgs.P2y,self.by[:,self.ie[2]],out=self.e[2].b)

        np.matmul(dgs.P1y,self.e[0].f,out=self.e[0].F)
        np.matmul(dgs.P2y,self.e[2].f,out=self.e[2].F)

        np.matmul(dgs.P2x,self.e[1].f,out=self.e[1].F)
        np.matmul(dgs.P1x,self.e[3].f,out=self.e[3].F)

        np.multiply(self.e[0].a,self.e[0].F,out=self.e[0].F)
        np.multiply(self.e[2].a,self.e[2].F,out=self.e[2].F)
        np.multiply(self.e[1].a,self.e[1].F,out=self.e[1].F)
        np.multiply(self.e[3].a,self.e[3].F,out=self.e[3].F)

        self.e[3].F[:,:]=0.0
        self.e[1].F[:,:]=0.0

        self.e[0].G[:,:]=0.0
        self.e[1].G[:,:]=0.0
        self.e[2].G[:,:]=0.0
        self.e[3].G[:,:]=0.0

#------------------------------------------------------------------------------
    def delta(self,p,dgs):
#------------------------------------------------------------------------------

        f=self.f[p]

        ix=dgs.ix
        iy=dgs.iy
        ie=dgs.ie

#------------------------------------------------------------------------------
#   delta: surface flux-x
#------------------------------------------------------------------------------

        np.matmul(dgs.P2x,self.ax,out=self.a[1])
        np.matmul(dgs.P1x,self.ax,out=self.a[3])

        np.matmul(dgs.P2x,self.bx,out=self.b[1])
        np.matmul(dgs.P1x,self.bx,out=self.b[3])

        self.a[1][:,ix[0]]=0.5*(self.a[1][:,ix[0]]+self.a[3][:,ix[1]])
        self.a[3][:,ix[1]]=self.a[1][:,ix[0]]

        self.b[1][:,ix[0]]=0.5*(self.b[1][:,ix[0]]+self.b[3][:,ix[1]])
        self.b[3][:,ix[1]]=self.b[1][:,ix[0]]

        np.matmul(dgs.P2x,f,out=self.F[1])
        np.matmul(dgs.P1x,f,out=self.F[3])

        np.multiply(self.a[1],self.F[1],out=self.F[1])
        np.multiply(self.a[3],self.F[3],out=self.F[3])

        self.F[1][self.a[1]<0.]=0.0
        self.F[3][self.a[3]>0.]=0.0

        self.F[1][:,ix[0]]=self.F[1][:,ix[0]]+self.F[3][:,ix[1]]
        self.F[3][:,ix[1]]=self.F[1][:,ix[0]]

        np.matmul(dgs.Q1x,f,out=self.G[1])
        np.matmul(dgs.Q2x,f,out=self.G[3])

        np.multiply(self.b[1],self.G[1],out=self.G[1])
        np.multiply(self.b[3],self.G[3],out=self.G[3])

        self.G[1][:,ix[0]]=self.G[1][:,ix[0]]+self.G[3][:,ix[1]]
        self.G[3][:,ix[1]]=self.G[1][:,ix[0]]

        np.subtract(self.F[1],self.G[1],out=self.F[1])
        np.subtract(self.F[3],self.G[3],out=self.F[3])

        self.F[3][:,ie[3]]=self.e[3].F[:,:]
        self.F[1][:,ie[1]]=self.e[1].F[:,:]

#------------------------------------------------------------------------------
#   delta: surface flux-y
#------------------------------------------------------------------------------

        np.matmul(dgs.P1y,self.ay,out=self.a[0])
        np.matmul(dgs.P2y,self.ay,out=self.a[2])

        np.matmul(dgs.P1y,self.by,out=self.b[0])
        np.matmul(dgs.P2y,self.by,out=self.b[2])

        self.a[2][:,iy[0]]=0.5*(self.a[2][:,iy[0]]+self.a[0][:,iy[1]])
        self.a[0][:,iy[1]]=self.a[2][:,iy[0]]

        self.b[2][:,iy[0]]=0.5*(self.b[2][:,iy[0]]+self.b[0][:,iy[1]])
        self.b[0][:,iy[1]]=self.b[2][:,iy[0]]

        np.matmul(dgs.P1y,f,out=self.F[0])
        np.matmul(dgs.P2y,f,out=self.F[2])

        np.multiply(self.a[0],self.F[0],out=self.F[0])
        np.multiply(self.a[2],self.F[2],out=self.F[2])

        self.F[0][self.a[0]>0.]=0.0
        self.F[2][self.a[2]<0.]=0.0

        self.F[2][:,iy[0]]=self.F[2][:,iy[0]]+self.F[0][:,iy[1]]
        self.F[0][:,iy[1]]=self.F[2][:,iy[0]]

        np.matmul(dgs.Q2y,f,out=self.G[0])
        np.matmul(dgs.Q1y,f,out=self.G[2])

        np.multiply(self.b[0],self.G[0],out=self.G[0])
        np.multiply(self.b[2],self.G[2],out=self.G[2])

        self.G[2][:,iy[0]]=self.G[2][:,iy[0]]+self.G[0][:,iy[1]]
        self.G[0][:,iy[1]]=self.G[2][:,iy[0]]

        np.subtract(self.F[0],self.G[0],out=self.F[0])
        np.subtract(self.F[2],self.G[2],out=self.F[2])

        self.F[2][:,ie[2]]=self.e[2].F[:,:]
        self.F[0][:,ie[0]]=self.e[0].F[:,:]

#------------------------------------------------------------------------------
#   delta: volume flux
#------------------------------------------------------------------------------

        np.multiply(self.ax,f,out=self.fx)
        np.multiply(self.ay,f,out=self.fy)

        np.matmul(dgs.Ux,f,out=self.gx)
        np.matmul(dgs.Uy,f,out=self.gy)

        np.multiply(self.bx,self.gx,out=self.gx)
        np.multiply(self.by,self.gy,out=self.gy)

        np.subtract(self.fx,self.gx,out=self.fx)
        np.subtract(self.fy,self.gy,out=self.fy)

        self.sx=self.jx*self.fx

        np.matmul(dgs.Cx,self.fx,out=self.fx)
        np.matmul(dgs.Cy,self.fy,out=self.fy)

#------------------------------------------------------------------------------
#   delta: delta
#------------------------------------------------------------------------------

        np.matmul(dgs.F1y,self.F[0],out=self.B[0])
        np.matmul(dgs.F2y,self.F[2],out=self.B[2])

        np.matmul(dgs.F1x,self.F[3],out=self.B[3])
        np.matmul(dgs.F2x,self.F[1],out=self.B[1])

        np.subtract(self.B[1],self.B[3],out=self.B[1])
        np.subtract(self.B[2],self.B[0],out=self.B[2])

        np.subtract(self.B[1],self.fx,out=self.fx)
        np.subtract(self.B[2],self.fy,out=self.fy)

        np.multiply(self.g[0],self.fx,out=self.fx)
        np.multiply(self.g[1],self.fy,out=self.fy)

        np.add(self.fx,self.fy,out=self.df[0])

        self.df[0]=-self.df[0]-self.sx

#        self.df[0]=-self.jx*self.df[0]

#------------------------------------------------------------------------------
    def split(self,p,dgs):
#------------------------------------------------------------------------------

        n=self.n
        m=self.m[0]

        iq=dgs.iq
        iv=dgs.iv

        q=np.zeros((n*n,4*m*m))

        q[:,iq[0]]=np.matmul(dgs.T[0],self.f[p])
        q[:,iq[1]]=np.matmul(dgs.T[1],self.f[p])
        q[:,iq[2]]=np.matmul(dgs.T[2],self.f[p])
        q[:,iq[3]]=np.matmul(dgs.T[3],self.f[p])

        f=[None,None,None,None]

        f[0]=q[:,iv[0]]
        f[1]=q[:,iv[1]]
        f[2]=q[:,iv[2]]
        f[3]=q[:,iv[3]]

        return f

#------------------------------------------------------------------------------
    def merge(self,p,f,dgs):
#------------------------------------------------------------------------------

        n=self.n
        m=self.m[0]

        iq=dgs.iq
        iv=dgs.iv

        q=np.zeros((n*n,4*m*m))

        q[:,iv[0]]=f[0][:,:]
        q[:,iv[1]]=f[1][:,:]
        q[:,iv[2]]=f[2][:,:]
        q[:,iv[3]]=f[3][:,:]

        self.f[p][:,:]=0.0

        self.f[p]=self.f[p]+np.matmul(dgs.S[0],q[:,iq[0]])
        self.f[p]=self.f[p]+np.matmul(dgs.S[1],q[:,iq[1]])
        self.f[p]=self.f[p]+np.matmul(dgs.S[2],q[:,iq[2]])
        self.f[p]=self.f[p]+np.matmul(dgs.S[3],q[:,iq[3]])

#------------------------------------------------------------------------------
    def getsol(self,p):
#------------------------------------------------------------------------------

        n=self.n

        mx=self.m[1]
        my=self.m[0]

        f=self.f[p]

        f=f.T.reshape((mx*my,n,n)).reshape((my,mx,n,n))

        f=np.concatenate(np.concatenate(f,axis=1),axis=1)

        return f

#------------------------------------------------------------------------------
    def getlev(self,p,l,lgn):
#------------------------------------------------------------------------------

        n=self.n

        mx=self.m[1]
        my=self.m[0]

        f=self.f[p]

        f=f.T.reshape((mx*my,n,n)).reshape((my,mx,n,n))

        f=np.matmul(lgn.T[l],f)
        f=np.matmul(f,lgn.T[l].T)

        f=np.concatenate(np.concatenate(f,axis=1),axis=1)

        return f

#------------------------------------------------------------------------------
    def geterr(self,err,lgn):
#------------------------------------------------------------------------------

        L1=np.amax(abs(err))

        L2=np.matmul(lgn.C2,err)

        L2=np.sum(np.sum(L2*L2,axis=0))*(self.h[0]*self.h[1])

        L2=np.sqrt(L2)

        return L1,L2

#------------------------------------------------------------------------------
    def geteps(self,lgn):
#------------------------------------------------------------------------------

        E=self.f[6]

        eps=(119.44+(4.3666e7)/(E*E*E))*np.exp(-273./E)

        eps=0.25*np.sum(lgn.I*eps,axis=0)*self.h[0]

        eps=np.amax(eps)

        return eps

###############################################################################
