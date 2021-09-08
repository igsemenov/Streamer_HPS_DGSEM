###############################################################################
" pylib_fvs "
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
###############################################################################
def limiter(a,b):

    m=np.zeros_like(a)

#    m=0.5*(np.sign(a)+np.sign(b))*np.minimum(abs(a),abs(b))

#    m1=0.5*(np.sign(a)+np.sign(b))*np.minimum(abs(2*a),abs(b))
#    m2=0.5*(np.sign(a)+np.sign(b))*np.minimum(abs(a),abs(2*b))
#    m=0.5*(np.sign(a)+np.sign(b))*np.maximum(abs(m1),abs(m2))

    m=np.zeros_like(a)

    aa=a*a
    ab=a*b

    c1=.25*ab
    c2=2.5*ab

    i1=(aa<=c1)
    i2=(aa>=c1)&(aa<c2)
    i3=(aa>=c2)&(ab>0.)

    m[i1]=a[i1]
    m[i2]=(a[i2]+a[i2]+b[i2])/6.
    m[i3]=b[i3]

    m=2.*m

    return m

###############################################################################
class fvs_edg:
#------------------------------------------------------------------------------
    def __init__(self,k):
#------------------------------------------------------------------------------

        self.k=k

        self.f=np.zeros((k,2))
        self.a=np.zeros((k,2))
        self.b=np.zeros((k,2))

        self.s=np.zeros(k)
        self.d=np.zeros(k)
        self.q=np.zeros(k)

        self.F=np.zeros(k)

#------------------------------------------------------------------------------
    def split(self,edg1,edg2):

        k=self.k/2

        edg1.f[:,:]=self.f[0:k,:]
        edg2.f[:,:]=self.f[k::,:]

        edg1.a[:,:]=self.a[0:k,:]
        edg2.a[:,:]=self.a[k::,:]

        edg1.b[:,:]=self.b[0:k,:]
        edg2.b[:,:]=self.b[k::,:]

#------------------------------------------------------------------------------
    def refine(self,edg,T):

        edg.f[1:-2:2,:]=0.75*self.f[0:-1,:]+0.25*self.f[1::,:]
        edg.f[2:-1:2,:]=0.25*self.f[0:-1,:]+0.75*self.f[1::,:]

        edg.f[0,:]=1.25*self.f[0,:]-0.25*self.f[1,:]
        edg.f[-1,:]=1.25*self.f[-1,:]-0.25*self.f[-2,:]

        np.matmul(edg.f,T,out=edg.f)

        edg.a[1:-2:2,:]=0.75*self.a[0:-1,:]+0.25*self.a[1::,:]
        edg.a[2:-1:2,:]=0.25*self.a[0:-1,:]+0.75*self.a[1::,:]

        edg.a[0,:]=1.25*self.a[0,:]-0.25*self.a[1,:]
        edg.a[-1,:]=1.25*self.a[-1,:]-0.25*self.a[-2,:]

        np.matmul(edg.a,T,out=edg.a)

        edg.b[1:-2:2,:]=0.75*self.b[0:-1,:]+0.25*self.b[1::,:]
        edg.b[2:-1:2,:]=0.25*self.b[0:-1,:]+0.75*self.b[1::,:]

        edg.b[0,:]=1.25*self.b[0,:]-0.25*self.b[1,:]
        edg.b[-1,:]=1.25*self.b[-1,:]-0.25*self.b[-2,:]

        np.matmul(edg.b,T,out=edg.b)

        edg.b=2.*edg.b

###############################################################################
class fvs_rib:
#------------------------------------------------------------------------------
    def __init__(self,box,ind):
#------------------------------------------------------------------------------

        k=box[0].m*box[0].s

        self.i=ind

        self.e=None

        self.a=np.zeros(k)
        self.b=np.zeros(k)
        self.s=np.zeros(k)

        self.r1=None
        self.r2=None

        self.r1_=None
        self.r2_=None

        self.T=[np.array([[1.25,-0.25],[0.75,0.25]]).T,
                np.array([[0.75,0.25],[-0.25,1.25]]).T]

        k_=[[2,0],[1,3],[0,2],[2,0],[3,1],[1,3]]

        if ind>1:

            k1=k_[ind][0]
            k2=k_[ind][1]

            self.e=[box[0].e[k1],box[1].e[k1],box[2].e[k2]]

            self.r1_=fvs_edg(k/2)
            self.r2_=fvs_edg(k/2)

            self.r1=fvs_edg(k)
            self.r2=fvs_edg(k)

        else:

            k1=k_[ind][0]
            k2=k_[ind][1]

            self.e=[box[0].e[k1],box[1].e[k2]]

#------------------------------------------------------------------------------
    def flux(self):
        
        if self.i==0:
            self.flux_(self.e[0],self.e[1])
        if self.i==1:
            self.flux_(self.e[0],self.e[1])

        if (self.i==2)|(self.i==4):

            self.e[2].split(self.r1_,self.r2_)

            self.r1_.refine(self.r1,self.T[1])
            self.r2_.refine(self.r2,self.T[1])

            self.flux_(self.r1,self.e[0])
            self.flux_(self.r2,self.e[1])

            self.r1_.F=0.5*self.r1.F[0::2]+0.5*self.r1.F[1::2]
            self.r2_.F=0.5*self.r2.F[0::2]+0.5*self.r2.F[1::2]

            self.e[2].F=np.hstack([self.r1_.F,self.r2_.F])

            self.r1_.q=(self.e[0].f[0::2,0]+self.e[0].f[1::2,0]+\
                        self.e[0].f[0::2,1]+self.e[0].f[1::2,1])*0.25

            self.r2_.q=(self.e[1].f[0::2,0]+self.e[1].f[1::2,0]+\
                        self.e[1].f[0::2,1]+self.e[1].f[1::2,1])*0.25

            self.e[2].q=np.hstack([self.r1_.q,self.r2_.q])

        if (self.i==3)|(self.i==5):

            self.e[2].split(self.r1_,self.r2_)

            self.r1_.refine(self.r1,self.T[0])
            self.r2_.refine(self.r2,self.T[0])

            self.flux_(self.e[0],self.r1)
            self.flux_(self.e[1],self.r2)

            self.r1_.F=0.5*self.r1.F[0::2]+0.5*self.r1.F[1::2]
            self.r2_.F=0.5*self.r2.F[0::2]+0.5*self.r2.F[1::2]

            self.e[2].F=np.hstack([self.r1_.F,self.r2_.F])

            self.r1_.q=(self.e[0].f[0::2,0]+self.e[0].f[1::2,0]+\
                        self.e[0].f[0::2,1]+self.e[0].f[1::2,1])*0.25

            self.r2_.q=(self.e[1].f[0::2,0]+self.e[1].f[1::2,0]+\
                        self.e[1].f[0::2,1]+self.e[1].f[1::2,1])*0.25

            self.e[2].q=np.hstack([self.r1_.q,self.r2_.q])

#------------------------------------------------------------------------------
    def flux_(self,e1,e2):

        self.a=0.5*(e1.a[:,1]+e2.a[:,0])
        self.b=0.5*(e1.b[:,1]+e2.b[:,0])

        e1.s=e1.f[:,1]-e1.f[:,0]
        e2.s=e2.f[:,1]-e2.f[:,0]

        self.s=e2.f[:,0]-e1.f[:,1]

        e1.d=limiter(e1.s,self.s)
        e2.d=limiter(self.s,e2.s)

        e1.q=e1.f[:,1]+0.5*e1.d
        e2.q=e2.f[:,0]-0.5*e2.d

        e1.F=e1.q*np.maximum(self.a,0.)+e2.q*np.minimum(self.a,0.)

        e1.F=e1.F-self.b*self.s

        e2.F[:]=e1.F[:]

        e1.q[:]=e2.f[:,0]
        e2.q[:]=e1.f[:,1]

###############################################################################
class fvs_box2d:
#------------------------------------------------------------------------------
    def __init__(self,c,h,fvs):
#------------------------------------------------------------------------------

        m=fvs.m

        self.m=fvs.m
        self.n=fvs.n
        self.s=fvs.s

        hx=h[0]/m
        hy=h[0]/m

        cx=np.arange(m)*(2.*hx)+(c[0]-h[0]+hx)
        cy=np.arange(m)*(2.*hy)+(c[1]-h[1]+hy)

        x=fvs.x[...,None]*hx+cx[None,...]
        y=fvs.x[...,None]*hy+cy[None,...]

        x=np.hstack(x.T)
        y=np.hstack(y.T)

        [x,y]=np.meshgrid(x,y)

        self.x=x
        self.y=y

        self.hx=hx*fvs.h
        self.hy=hy*fvs.h

        self.gx=1./self.hx
        self.gy=1./self.hy

#------------------------------------------------------------------------------

        k=x.shape[0]

        self.f=[np.zeros_like(x),np.zeros_like(x),
                np.zeros_like(x),np.zeros_like(x),
                np.zeros_like(x),np.zeros_like(x),
                np.zeros_like(x),np.zeros_like(x)]

        self.df=[np.zeros_like(x),np.zeros_like(x)]

        self.ax=np.zeros_like(x)
        self.ay=np.zeros_like(x)

        self.bx=np.zeros_like(x)
        self.by=np.zeros_like(x)

        self.ax_=0.5*self.ax[:,0:-1]+0.5*self.ax[:,1::]
        self.bx_=0.5*self.bx[:,0:-1]+0.5*self.bx[:,1::]

        self.ay_=0.5*self.ay[0:-1,:]+0.5*self.ay[1::,:]
        self.by_=0.5*self.by[0:-1,:]+0.5*self.by[1::,:]

        self.dx=np.zeros_like(x)
        self.dy=np.zeros_like(x)

        self.F=[np.zeros_like(x),np.zeros_like(x),
                np.zeros_like(x),np.zeros_like(x)]

        self.sx=np.zeros((k,k+1))
        self.sy=np.zeros((k+1,k))

        self.Fx=np.zeros((k,k+1))
        self.Fy=np.zeros((k+1,k))

#------------------------------------------------------------------------------

        self.e=[fvs_edg(k),fvs_edg(k),fvs_edg(k),fvs_edg(k)]

#------------------------------------------------------------------------------

#        self.ax.fill(0.)
#        self.ay.fill(0.)

        self.bx.fill(0.0219*self.gx)
        self.by.fill(0.0180*self.gy)

#------------------------------------------------------------------------------
    def  edge(self,p):
#------------------------------------------------------------------------------

        f=self.f[p]

        self.e[0].f[:,:]=f[0:2,:].T
        self.e[2].f[:,:]=f[-2:,:].T
        self.e[1].f[:,:]=f[:,-2:]
        self.e[3].f[:,:]=f[:,0:2]

        self.e[0].a[:,:]=self.ay[0:2,:].T
        self.e[2].a[:,:]=self.ay[-2:,:].T
        self.e[1].a[:,:]=self.ax[:,-2:]
        self.e[3].a[:,:]=self.ax[:,0:2]

        self.e[0].b[:,:]=self.by[0:2,:].T
        self.e[2].b[:,:]=self.by[-2:,:].T
        self.e[1].b[:,:]=self.bx[:,-2:]
        self.e[3].b[:,:]=self.bx[:,0:2]

        self.e[0].F=f[0,:]*self.ay[0,:]
        self.e[2].F=f[-1,:]*self.ay[-1,:]
        self.e[1].F=f[:,-1]*self.ax[:,-1]
        self.e[3].F=f[:,0]*self.ax[:,0]

        self.e[3].F[:]=0.0
        self.e[1].F[:]=0.0

#------------------------------------------------------------------------------
    def delta(self,p):

        f=self.f[p]

        np.add(self.ax[:,0:-1],self.ax[:,1::],out=self.ax_)
        np.add(self.bx[:,0:-1],self.bx[:,1::],out=self.bx_)

        np.multiply(self.ax_,0.5,out=self.ax_)
        np.multiply(self.bx_,0.5,out=self.bx_)

        np.add(self.ay[0:-1,:],self.ay[1::,:],out=self.ay_)
        np.add(self.by[0:-1,:],self.by[1::,:],out=self.by_)

        np.multiply(self.ay_,0.5,out=self.ay_)
        np.multiply(self.by_,0.5,out=self.by_)

        self.sx[:,1:-1]=f[:,1::]-f[:,:-1]
        self.sy[1:-1,:]=f[1::,:]-f[:-1,:]

        self.sx[:,0]=f[:,0]-self.e[3].q[:]
        self.sx[:,-1]=self.e[1].q[:]-f[:,-1]

        self.sy[0,:]=f[0,:]-self.e[0].q[:]
        self.sy[-1,:]=self.e[2].q[:]-f[-1,:]

        self.dx[:,:]=0.
        self.dy[:,:]=0.

        self.dx[:,:]=limiter(self.sx[:,0:-1],self.sx[:,1::])
        self.dy[:,:]=limiter(self.sy[0:-1,:],self.sy[1::,:])

        self.F[1]=f+0.5*self.dx
        self.F[3]=f-0.5*self.dx

        self.F[0]=f-0.5*self.dy
        self.F[2]=f+0.5*self.dy

        self.Fx[:,1:-1]=self.F[1][:,:-1]*np.maximum(self.ax_,0.)+\
                        self.F[3][:,1::]*np.minimum(self.ax_,0.)

        self.Fy[1:-1,:]=self.F[2][:-1,:]*np.maximum(self.ay_,0.)+\
                        self.F[0][1::,:]*np.minimum(self.ay_,0.)

        self.Fx[:,1:-1]=self.Fx[:,1:-1]-self.sx[:,1:-1]*self.bx_
        self.Fy[1:-1,:]=self.Fy[1:-1,:]-self.sy[1:-1,:]*self.by_

        self.Fx[:,0]=self.e[3].F[:]
        self.Fx[:,-1]=self.e[1].F[:]

        self.Fy[0,:]=self.e[0].F[:]
        self.Fy[-1,:]=self.e[2].F[:]

        self.df[0]=-(self.Fx[:,1::]-self.Fx[:,:-1])*self.gx-\
                    (self.Fy[1::,:]-self.Fy[:-1,:])*self.gy

        self.dx=0.5*(self.Fx[:,:-1]+self.Fx[:,1::])/self.x
        self.df[0]=self.df[0]-self.dx

#        self.df[0]=self.df[0]/self.x

#------------------------------------------------------------------------------
    def split(self,p,fvs):
#------------------------------------------------------------------------------

        f=self.f[p]

        m=self.m
        s=self.s

        g=np.zeros((2*m,2*m,s,s))

        f=np.asarray(np.split(f,m,axis=1))
        f=np.asarray(np.split(f,m,axis=1))

        g[0::2,0::2,:,:]=np.matmul(np.matmul(fvs.T1,f),fvs.T1.T)
        g[1::2,0::2,:,:]=np.matmul(np.matmul(fvs.T2,f),fvs.T1.T)
        g[0::2,1::2,:,:]=np.matmul(np.matmul(fvs.T1,f),fvs.T2.T)
        g[1::2,1::2,:,:]=np.matmul(np.matmul(fvs.T2,f),fvs.T2.T)

        q=[None,None,None,None]

        q[0]=np.concatenate(np.concatenate(g[0:m,0:m,:,:],axis=1),axis=1)
        q[1]=np.concatenate(np.concatenate(g[m::,0:m,:,:],axis=1),axis=1)
        q[2]=np.concatenate(np.concatenate(g[0:m,m::,:,:],axis=1),axis=1)
        q[3]=np.concatenate(np.concatenate(g[m::,m::,:,:],axis=1),axis=1)

        return q

#------------------------------------------------------------------------------
    def merge(self,p,f,fvs):
#------------------------------------------------------------------------------

        m=self.m

        f=np.hstack([np.vstack(f[0:2]),np.vstack(f[2:4])])

        f=np.asarray(np.split(f,2*m,axis=1))
        f=np.asarray(np.split(f,2*m,axis=1))

        g=np.matmul(fvs.S1,f[0::2,:])+np.matmul(fvs.S2,f[1::2,:])

        g=np.matmul(g[:,0::2],fvs.S1.T)+np.matmul(g[:,1::2],fvs.S2.T)

        g=np.concatenate(np.concatenate(g,axis=1),axis=1)

        self.f[p][:,:]=g[:,:]

#        return g

#------------------------------------------------------------------------------
    def getlev(self,f,l,lgn):
#------------------------------------------------------------------------------

        n=self.n

        mx=self.m
        my=self.m

        f=f.T.reshape((mx*my,n,n)).reshape((my,mx,n,n))

        f=np.matmul(lgn.T[l],f)
        f=np.matmul(f,lgn.T[l].T)

        f=np.concatenate(np.concatenate(f,axis=1),axis=1)

        return f

#------------------------------------------------------------------------------
    def geteps(self):

        E=self.f[6]

        Ki=(119.44+(4.3666e7)/(E*E*E))*np.exp(-273./E)

        eps1=Ki*self.hx

        eps2=E*E

        eps1=np.amax(eps1)
        eps2=np.amax(eps2)

        return eps1,eps2

###############################################################################

