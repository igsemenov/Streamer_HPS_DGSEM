###############################################################################
" pylib_lgn "
###############################################################################
import numpy as np
from scipy.linalg import block_diag as bdiag
###############################################################################
def getpol(x,n):

    P=np.zeros((len(x),n))

    P[:,0]=1.
    P[:,1]=x

    for k in range(2,P.shape[1]):
        P[:,k]=(2.*k-1.)*x*P[:,k-1]-(k-1.)*P[:,k-2]
        P[:,k]=P[:,k]/(k-0.)

    return P

###############################################################################
def permat(n,m):

    P=np.zeros((n*m,n*m))

    for j in range(m):
        for i in range(n):
            P[i*m+j,i+j*n]=1

    return P

###############################################################################
class lgn1d:
#------------------------------------------------------------------------------
    def __init__(self,n):
#------------------------------------------------------------------------------

        x=np.cos(np.pi*(4.*np.arange(n,0,-1)-1.)/(4.*n+2))

        P=np.zeros((n,n+2))
        G=np.zeros((n,n+2))
    
        for i in range(7):
    
            P[:,0]=1.
            P[:,1]=x
    
            G[:,0]=0.
            G[:,1]=1.
    
            for k in range(2,P.shape[1]):
     
                P[:,k]=(2.*k-1.)*x*P[:,k-1]-(k-1.)*P[:,k-2]
                G[:,k]=(2.*k-1.)*x*G[:,k-1]-(k-0.)*G[:,k-2]
    
                P[:,k]=P[:,k]/(k-0.)
                G[:,k]=G[:,k]/(k-1.)
    
            err=P[:,n]/G[:,n]
    
            x=x-err
    
        I=(2.)/((1.-x*x)*G[:,n]*G[:,n])
    
        k=np.arange(n)+0.5
    
        C=k[...,None]*(I*P[:,0:n].transpose())

#------------------------------------------------------------------------------

        U=G[:,0:n].dot(C)

        L=np.zeros((n,n))

        L[:,0]=0.0
        L[:,1]=0.0
        L[:,2]=3.0
    
        for k in range(3,L.shape[1]):
            L[:,k]=(2.*k-1.)*x*L[:,k-1]-(k+1.)*L[:,k-2]
            L[:,k]=L[:,k]/(k-2.)
    
        L=L[:,0:n].dot(C)
    
        Pb=np.ones(n+2);
        Pa=np.ones(n+2); Pa[1::2]=-1;

#------------------------------------------------------------------------------

        Iax=np.zeros((n,n+1))
        Ixb=np.zeros((n,n+1))
    
        Max=np.zeros((n,n))
        Mxb=np.zeros((n,n))
    
        Iax[:,0]=1.+x
        Ixb[:,0]=1.-x
    
        for k in range(1,n+1):
            Iax[:,k]=((P[:,k+1]-Pa[k+1])-(P[:,k-1]-Pa[k-1]))/(2.*k+1.)
            Ixb[:,k]=((Pb[k+1]-P[:,k+1])-(Pb[k-1]-P[:,k-1]))/(2.*k+1.)

        Max[:,0]=0.5*(x*x-1.)
        Mxb[:,0]=0.5*(1.-x*x)
    
        for k in range(1,n):
            Max[:,k]=((k+1.)*Iax[:,k+1]+k*Iax[:,k-1])/(2.*k+1.)
            Mxb[:,k]=((k+1.)*Ixb[:,k+1]+k*Ixb[:,k-1])/(2.*k+1.)

        Iax=Iax[:,0:-1].dot(C)
        Ixb=Ixb[:,0:-1].dot(C)

        Max=Max.dot(C)
        Mxb=Mxb.dot(C)

        P1=0.5*(x[...,None]-1.)
        P2=0.5*(x[...,None]+1.)

        G1=(Max+Iax)
        G2=(Mxb-Ixb)

        G=P1*G1+P2*G2
        H=0.5*(G1+G2)

        Pa=Pa[0:n].dot(C)
        Pb=Pb[0:n].dot(C)

#------------------------------------------------------------------------------

        x1=-0.5+0.5*x
        x2=+0.5+0.5*x
    
        T1=np.matmul(getpol(x1,n),C)
        T2=np.matmul(getpol(x2,n),C)

        y1=(x[0:n/2]+0.5)*2.
        y2=(x[n/2::]-0.5)*2.

        S1=np.matmul(getpol(y1,n),C)
        S2=np.matmul(getpol(y2,n),C)

#------------------------------------------------------------------------------

        self.n=n

        self.x=x

        self.C=C; self.I=I;
        self.U=U; self.L=L;
        self.G=G; self.H=H;

        self.Pa=Pa; self.Pb=Pb;

        self.T1=T1; self.T2=T2;
        self.S1=S1; self.S2=S2;

        self.C2=np.sqrt(2./(2.*np.arange(n)[...,None]+1.))*C

        if (n==4)|(n==6):
            self.lbt=lbt1d(self)

###############################################################################
class lbt1d:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------

        n=lgn.n

        if n==4:

            x=np.array([-1.,-1./np.sqrt(5.),1./np.sqrt(5.),+1])

            I=np.array([1./6.,5./6.,5./6.,1./6.])

        if n==6:

            x=np.array([-1.,-np.sqrt((7.+2.*np.sqrt(7.))/21.),
                             -np.sqrt((7.-2.*np.sqrt(7.))/21.),
                              np.sqrt((7.-2.*np.sqrt(7.))/21.),
                              np.sqrt((7.+2.*np.sqrt(7.))/21.),1.])
            
            I=np.array([1./15,(14.-np.sqrt(7.))/30.,
                        (14.+np.sqrt(7.))/30.,
                        (14.+np.sqrt(7.))/30.,
                        (14.-np.sqrt(7.))/30.,1./15])

#------------------------------------------------------------------------------
#       self.T_lg: gauss (lgn) -> lobatto (lbt) transform
#       self.T_gl: lobatto (lbt) -> gauss (lgn) transform
#------------------------------------------------------------------------------

        T_lg=np.dot(getpol(x,n),lgn.C)
        T_gl=np.linalg.solve(T_lg,np.eye(n))

        U=np.matmul(T_lg,np.matmul(lgn.U,T_gl))

#------------------------------------------------------------------------------

        T1=np.matmul(T_lg,np.matmul(lgn.T1,T_gl))
        T2=np.matmul(T_lg,np.matmul(lgn.T2,T_gl))

#------------------------------------------------------------------------------

        self.n=n; self.x=x;
        self.I=I; self.U=U;

        self.T_lg=T_lg
        self.T_gl=T_gl

        self.T1=T1
        self.T2=T2

###############################################################################
class lgn2d:
#------------------------------------------------------------------------------
    def __init__(self,n,m):
#------------------------------------------------------------------------------

        self.n=n
        self.m=m

        self.lgn=lgn1d(n)

        P=permat(n,n)

        T1=bdiag(*[self.lgn.T1 for i in range(n)])
        T2=bdiag(*[self.lgn.T2 for i in range(n)])

        self.T=[None,None,None,None]

        self.T[0]=np.matmul(P,np.matmul(T1,np.matmul(P,T1)))
        self.T[1]=np.matmul(P,np.matmul(T2,np.matmul(P,T1)))
        self.T[2]=np.matmul(P,np.matmul(T1,np.matmul(P,T2)))
        self.T[3]=np.matmul(P,np.matmul(T2,np.matmul(P,T2)))

        T=np.vstack([self.lgn.T1,self.lgn.T2])
        S=bdiag(self.lgn.S1,self.lgn.S2)

        self.T_=[np.eye(m*n),bdiag(*[T for i in range(m)]),
                             bdiag(*[S for i in range(m)])]

#        C2=bdiag(*[lgn.C2 for i in range(n)])
#        self.C2=np.matmul(P,np.matmul(C2,np.matmul(P,C2)))

#------------------------------------------------------------------------------

        self.lbt=lbt2d(self.lgn)

###############################################################################
class lbt2d:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------

        n=lgn.n

        P=permat(n,n)

        T=bdiag(*[lgn.lbt.T_lg for i in range(n)])
        S=bdiag(*[lgn.lbt.T_gl for i in range(n)])

        self.T_lg=np.matmul(P,np.matmul(T,np.matmul(P,T)))
        self.T_gl=np.matmul(P,np.matmul(S,np.matmul(P,S)))

        T1=bdiag(*[lgn.lbt.T1 for i in range(n)])
        T2=bdiag(*[lgn.lbt.T2 for i in range(n)])

        self.T=[None,None,None,None]

        self.T[0]=np.matmul(P,np.matmul(T1,np.matmul(P,T1)))
        self.T[1]=np.matmul(P,np.matmul(T2,np.matmul(P,T1)))
        self.T[2]=np.matmul(P,np.matmul(T1,np.matmul(P,T2)))
        self.T[3]=np.matmul(P,np.matmul(T2,np.matmul(P,T2)))

###############################################################################
class ope2d:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------
        
        n=lgn.n

        P=permat(n,n)
        T=P.T

        D=np.linalg.solve(lgn.G,np.eye(n))

        G=bdiag(*[lgn.G for i in range(n)])
        H=bdiag(*[lgn.H for i in range(n)])

        D=bdiag(*[D for i in range(n)])

        Q=np.matmul(P,np.matmul(G,np.matmul(P,G)))
        D=np.matmul(D,np.matmul(T,np.matmul(D,T)))

        self.Q=Q

        self.Ux=np.matmul(P,np.matmul(G,np.matmul(P,H)))
        self.Uy=np.matmul(P,np.matmul(H,np.matmul(P,G)))

        self.Ly=G
        self.Lx=np.matmul(P,np.matmul(G,P))

###############################################################################
class opb2d:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------

        n=lgn.n
    
        P=permat(n,n)
    
        p1=0.5*(1.-lgn.x)
        p2=0.5*(1.+lgn.x)
    
        P1x=bdiag(*[p1[...,None] for i in range(n)])
        P2x=bdiag(*[p2[...,None] for i in range(n)])
    
        Ix=bdiag(*[np.ones((n,1)) for i in range(n)])
    
        P1y=P.dot(P1x)
        P2y=P.dot(P2x)
    
        Iy=P.dot(Ix)
    
        G=lgn.G
        H=lgn.H

        S=np.linalg.solve(G,np.eye(n))
        N=H.dot(S)
    
        n1=N.dot(p1)
        n2=N.dot(p2)
    
        d1=S.dot(p1)
        d2=S.dot(p2)
    
#------------------------------------------------------------------------------
    
        Q=[0,0];
    
        Ux=[0,0]; Uy=[0,0];
        Lx=[0,0]; Ly=[0,0];
    
#------------------------------------------------------------------------------
#   Corner nodes:
#------------------------------------------------------------------------------
    
        Q[0]=[-P1y.dot(p1),-P1y.dot(p2),-P2y.dot(p2),-P2y.dot(p1)]
    
        Ux[0]=[-P1y.dot(n1),-P1y.dot(n2),-P2y.dot(n2),-P2y.dot(n1)]
        Uy[0]=[-P1x.dot(n1),-P2x.dot(n1),-P2x.dot(n2),-P1x.dot(n2)]
    
        Lx[0]=[-P1y.dot(d1),-P1y.dot(d2),-P2y.dot(d2),-P2y.dot(d1)]
        Ly[0]=[-P1x.dot(d1),-P2x.dot(d1),-P2x.dot(d2),-P1x.dot(d2)]
    
#------------------------------------------------------------------------------
#   Edge nodes:
#------------------------------------------------------------------------------
    
        Q[1]=[P1y,P2x,P2y,P1x]
    
        Ux[1]=[P1y.dot(N),0.5*Ix,P2y.dot(N),-0.5*Ix]
        Uy[1]=[-0.5*Iy,P2x.dot(N),0.5*Iy,P1x.dot(N)]
    
        Lx[1]=[P1y.dot(S),Ix*0.0,P2y.dot(S),Ix*0.0]
        Ly[1]=[Iy*0.0,P2x.dot(S),Iy*0.0,P1x.dot(S)]
    
#------------------------------------------------------------------------------
#   Exclusion of the corner nodes:
#------------------------------------------------------------------------------
    
        P=[lgn.Pa,lgn.Pb]; a=0; b=1;
    
        c=[[(0,a),(1,b)],[(1,a),(2,b)],[(3,a),(2,b)],[(0,a),(3,b)]]
    
        for k in range(4):
            for (i,j) in c[k]:
    
                Q[1][k]=Q[1][k]+0.5*np.outer(Q[0][i],P[j])
    
                Ux[1][k]=Ux[1][k]+0.5*np.outer(Ux[0][i],P[j])
                Uy[1][k]=Uy[1][k]+0.5*np.outer(Uy[0][i],P[j])
    
                Lx[1][k]=Lx[1][k]+0.5*np.outer(Lx[0][i],P[j])
                Ly[1][k]=Ly[1][k]+0.5*np.outer(Ly[0][i],P[j])
    
        self.Q=np.hstack(Q[1])
    
        self.Ux=np.hstack(Ux[1])
        self.Uy=np.hstack(Uy[1])
    
        self.Lx=np.hstack(Lx[1])
        self.Ly=np.hstack(Ly[1])

###############################################################################
class elm2d:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------

        n=lgn.n

        P=permat(n,n)

        self.n=n

        self.ope=ope2d(lgn)
        self.opb=opb2d(lgn)

        P1x=bdiag(*[lgn.Pa for i in range(n)])
        P2x=bdiag(*[lgn.Pb for i in range(n)])

        P1y=P1x.dot(P)
        P2y=P2x.dot(P)

        Ne1x=P1x.dot(self.ope.Ux)
        Ne2x=P2x.dot(self.ope.Ux)
        Nb1x=P1x.dot(self.opb.Ux)
        Nb2x=P2x.dot(self.opb.Ux)

        Ne1y=P1y.dot(self.ope.Uy)
        Ne2y=P2y.dot(self.ope.Uy)
        Nb1y=P1y.dot(self.opb.Uy)
        Nb2y=P2y.dot(self.opb.Uy)

        self.Ne=np.vstack([Ne1y,Ne2x,Ne2y,Ne1x])
        self.Nb=np.vstack([Nb1y,Nb2x,Nb2y,Nb1x])

###############################################################################
class dgs1d_lgn:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------

        M=lgn.I
        U=lgn.U

        Pa=lgn.Pa
        Pb=lgn.Pb

        Ua=np.dot(lgn.Pa,U)
        Ub=np.dot(lgn.Pb,U)

        Q=np.vstack([np.hstack([Pb,-Pa]),
                     np.hstack([Ub,-Ua])])

        u,s,v=np.linalg.svd(Q)

        Q=(v[2::,:].T).dot(v[2::,:])

        Q=np.hsplit(np.vsplit(Q,2)[0],2)

        J=1.0/M

        self.U=U
        self.C=J[...,None]*(U.T*M)

        self.P1=Pa
        self.P2=Pb

        self.F1=J*Pa
        self.F2=J*Pb

        self.Q1=Ub.dot(Q[0])
        self.Q2=Ub.dot(Q[1])

        self.P1=self.P1[None,...]
        self.P2=self.P2[None,...]

        self.F1=self.F1[...,None]
        self.F2=self.F2[...,None]

        self.Q1=self.Q1[None,...]
        self.Q2=self.Q2[None,...]

###############################################################################
class dgs1d_lbt:
#------------------------------------------------------------------------------
    def __init__(self,lbt):
#------------------------------------------------------------------------------

        n=lbt.n

        I=lbt.I 
        U=lbt.U

        J=1./I

        C=J[...,None]*(U.T*I)

#------------------------------------------------------------------------------

        Pa=np.eye(n)[:,0]
        Pb=np.eye(n)[:,-1]

        Ua=U[0,:]
        Ub=U[-1,:]

        Q=np.vstack([np.hstack([Pb,-Pa]),
                     np.hstack([Ub,-Ua])])

        u,s,v=np.linalg.svd(Q)

        Q=(v[2::,:].T).dot(v[2::,:])

        self.Q=Q

        Q=np.hsplit(np.vsplit(Q,2)[0],2)

#------------------------------------------------------------------------------

        self.n=n

        self.I=I
        self.J=J

        self.U=U 
        self.C=C

        self.P1=Pa
        self.P2=Pb

        self.F1=Pa*J[0]
        self.F2=Pb*J[-1]

        self.P1=self.P1[None,...]
        self.P2=self.P2[None,...]

        self.F1=self.F1[...,None]
        self.F2=self.F2[...,None]

        self.Q1=Ub.dot(Q[0])
        self.Q2=Ub.dot(Q[1])

        self.Q1=self.Q1[None,...]
        self.Q2=self.Q2[None,...]

###############################################################################
class dgs2d_lbt:
#------------------------------------------------------------------------------
    def __init__(self,lgn):
#------------------------------------------------------------------------------

        n=lgn.n

        P=permat(n,n)

        dgs=dgs1d_lbt(lgn)

        self.J=dgs.J

        self.Cx=bdiag(*[dgs.C for i in range(n)])
        self.Ux=bdiag(*[dgs.U for i in range(n)])

        self.Cy=np.matmul(P,np.matmul(self.Cx,P))
        self.Uy=np.matmul(P,np.matmul(self.Ux,P))

        self.Q1x=bdiag(*[dgs.Q1 for i in range(n)])
        self.Q2x=bdiag(*[dgs.Q2 for i in range(n)])

        self.Q1y=np.matmul(self.Q1x,P)
        self.Q2y=np.matmul(self.Q2x,P)

        self.P1x=bdiag(*[dgs.P1 for i in range(n)])
        self.P2x=bdiag(*[dgs.P2 for i in range(n)])

        self.F1x=bdiag(*[dgs.F1 for i in range(n)])
        self.F2x=bdiag(*[dgs.F2 for i in range(n)])

        self.P1y=np.matmul(self.P1x,P)
        self.P2y=np.matmul(self.P2x,P)

        self.F1y=np.matmul(P,self.F1x)
        self.F2y=np.matmul(P,self.F2x)

        self.T1=lgn.lbt.T1
        self.T2=lgn.lbt.T2

        self.K1=np.matmul(self.T1,np.eye(n)).T
        self.K2=np.matmul(self.T2,np.eye(n)).T

        self.K1=(0.5/lgn.lbt.I[...,None])*(self.K1*lgn.lbt.I)
        self.K2=(0.5/lgn.lbt.I[...,None])*(self.K2*lgn.lbt.I)

        T1=bdiag(*[self.T1 for i in range(n)])
        T2=bdiag(*[self.T2 for i in range(n)])

        self.T={}

        self.T[0]=np.matmul(P,np.matmul(T1,np.matmul(P,T1)))
        self.T[1]=np.matmul(P,np.matmul(T2,np.matmul(P,T1)))
        self.T[2]=np.matmul(P,np.matmul(T1,np.matmul(P,T2)))
        self.T[3]=np.matmul(P,np.matmul(T2,np.matmul(P,T2)))

###############################################################################
class nod1d:
#------------------------------------------------------------------------------
    def __init__(self,c,h,m,elm):
#------------------------------------------------------------------------------

        n=elm.n

        h_=h/m

        c=np.arange(m)*(2.*h_)+(c-h+h_)

        self.n=n
        self.m=m
        self.h=h_

        self.x=elm.x[...,None]*h_+c[None,...]

###############################################################################
class nod2d:
#------------------------------------------------------------------------------
    def __init__(self,c,h,m,elm):
#------------------------------------------------------------------------------

        n=elm.n

        mx=m
        my=m

        hx=h[0]/mx
        hy=h[1]/my

        cx_=np.arange(mx)*(2.*hx)+(c[0]-h[0]+hx)
        cy_=np.arange(my)*(2.*hy)+(c[1]-h[1]+hy)

        [cx,cy]=np.meshgrid(cx_,cy_)

        cx=cx.flatten()
        cy=cy.flatten()

        x,y=np.meshgrid(elm.x,elm.x)

        x=x.flatten()
        y=y.flatten()

        self.n=n

        self.m=(my,mx)
        self.h=(hx,hy)

        self.x=x[...,None]*hx+cx[None,...]
        self.y=y[...,None]*hy+cy[None,...]

#------------------------------------------------------------------------------

    def getmsh(self):

        n=self.n

        mx=self.m[1]
        my=self.m[0]

        x=self.x.T.reshape((mx*my,n,n)).reshape((my,mx,n,n))
        y=self.y.T.reshape((mx*my,n,n)).reshape((my,mx,n,n))

        x=np.concatenate(np.concatenate(x,axis=1),axis=1)
        y=np.concatenate(np.concatenate(y,axis=1),axis=1)

        return x,y

    def getedg(self):

        x=[q.T.flatten() for q in self.x_]
        y=[q.T.flatten() for q in self.y_]

        return x,y

###############################################################################
