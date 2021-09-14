###############################################################################
" pylib_psn"
###############################################################################
import numpy as np
###############################################################################
class psn_dtn:
#------------------------------------------------------------------------------
    def __init__(self):
#------------------------------------------------------------------------------

        self.i=1

        self.n=None

        self.A=None
        self.S=None
        self.Q=None

        self.Tab=None
        self.Tba=None

#------------------------------------------------------------------------------
    def mrgh(self,dtn1,dtn2):

        A=dtn1.A
        B=dtn2.A

        n=dtn1.n
        m=dtn2.n

        self.n=np.array([n[0]+m[0],m[1],n[2]+m[2],n[3]])

        self.S=np.linalg.solve(A[1,1]-B[3,3],np.eye(n[1]))

        Q=[-A[1,0],B[3,0],B[3,1],-A[1,2],B[3,2],-A[1,3]]

        Q=[np.matmul(self.S,q) for q in Q]

        C=[0,0,0,0,0,0]

        C[0]=[np.matmul(A[0,1],q) for q in Q]
        C[1]=[np.matmul(B[0,3],q) for q in Q]
        C[2]=[np.matmul(B[1,3],q) for q in Q]
        C[3]=[np.matmul(A[2,1],q) for q in Q]
        C[4]=[np.matmul(B[2,3],q) for q in Q]
        C[5]=[np.matmul(A[3,1],q) for q in Q]

        for [i,j] in [[0,0],[3,2],[5,3]]:

            C[i][0]=C[i][0]+A[j,0]
            C[i][3]=C[i][3]+A[j,2]
            C[i][5]=C[i][5]+A[j,3]

        for [i,j] in [[1,0],[2,1],[4,2]]:

            C[i][1]=C[i][1]+B[j,0]
            C[i][2]=C[i][2]+B[j,1]
            C[i][4]=C[i][4]+B[j,2]

        k=np.add.accumulate(self.n[0:-1])

        C=np.vstack([np.hstack(p) for p in C])
        C=[np.hsplit(p,k) for p in np.vsplit(C,k)]

        Q=np.hstack(Q)
        Q=np.hsplit(Q,k)

        self.A={(i,j): C[i][j] for i in range(4) for j in range(4)}
        self.Q={i: Q[i] for i in range(4)}

#------------------------------------------------------------------------------
    def mrgv(self,dtn1,dtn2):

        A=dtn1.A
        B=dtn2.A

        n=dtn1.n
        m=dtn2.n

        self.n=np.array([n[0],n[1]+m[1],m[2],n[3]+m[3]])

        self.S=np.linalg.solve(A[2,2]-B[0,0],np.eye(n[2]))

        Q=[-A[2,0],-A[2,1],B[0,1],B[0,2],-A[2,3],B[0,3]]

        Q=[np.matmul(self.S,q) for q in Q]

        C=[0,0,0,0,0,0]

        C[0]=[np.matmul(A[0,2],q) for q in Q]
        C[1]=[np.matmul(A[1,2],q) for q in Q]
        C[2]=[np.matmul(B[1,0],q) for q in Q]
        C[3]=[np.matmul(B[2,0],q) for q in Q]
        C[4]=[np.matmul(A[3,2],q) for q in Q]
        C[5]=[np.matmul(B[3,0],q) for q in Q]

        for [i,j] in [[0,0],[1,1],[4,3]]:

            C[i][0]=C[i][0]+A[j,0]
            C[i][1]=C[i][1]+A[j,1]
            C[i][4]=C[i][4]+A[j,3]

        for [i,j] in [[2,1],[3,2],[5,3]]:

            C[i][2]=C[i][2]+B[j,1]
            C[i][3]=C[i][3]+B[j,2]
            C[i][5]=C[i][5]+B[j,3]

        k=np.add.accumulate(self.n[0:-1])

        C=np.vstack([np.hstack(p) for p in C])
        C=[np.hsplit(p,k) for p in np.vsplit(C,k)]

        Q=np.hstack(Q)
        Q=np.hsplit(Q,k)

        self.A={(i,j): C[i][j] for i in range(4) for j in range(4)}
        self.Q={i: Q[i] for i in range(4)}

#------------------------------------------------------------------------------
    def mrgv_(self,dtn1,dtn2):

        A=dtn1.A
        B=dtn2.A

        n=dtn1.n
        m=dtn2.n

        self.n=np.array([n[0],n[1]+m[1],m[2],n[3]+m[3]])

        B[0,0]=np.matmul(B[0,0],self.Tba)
        B[0,0]=np.matmul(self.Tab,B[0,0])

        B[0,1]=np.matmul(self.Tab,B[0,1])
        B[0,2]=np.matmul(self.Tab,B[0,2])
        B[0,3]=np.matmul(self.Tab,B[0,3])

        B[1,0]=np.matmul(B[1,0],self.Tba)
        B[2,0]=np.matmul(B[2,0],self.Tba)
        B[3,0]=np.matmul(B[3,0],self.Tba)

        self.S=np.linalg.solve(A[2,2]-B[0,0],np.eye(n[2]))

        Q=[-A[2,0],-A[2,1],B[0,1],B[0,2],-A[2,3],B[0,3]]

        Q=[np.matmul(self.S,q) for q in Q]

        C=[0,0,0,0,0,0]

        C[0]=[np.matmul(A[0,2],q) for q in Q]
        C[1]=[np.matmul(A[1,2],q) for q in Q]
        C[2]=[np.matmul(B[1,0],q) for q in Q]
        C[3]=[np.matmul(B[2,0],q) for q in Q]
        C[4]=[np.matmul(A[3,2],q) for q in Q]
        C[5]=[np.matmul(B[3,0],q) for q in Q]

        for [i,j] in [[0,0],[1,1],[4,3]]:

            C[i][0]=C[i][0]+A[j,0]
            C[i][1]=C[i][1]+A[j,1]
            C[i][4]=C[i][4]+A[j,3]

        for [i,j] in [[2,1],[3,2],[5,3]]:

            C[i][2]=C[i][2]+B[j,1]
            C[i][3]=C[i][3]+B[j,2]
            C[i][5]=C[i][5]+B[j,3]

        k=np.add.accumulate(self.n[0:-1])

        C=np.vstack([np.hstack(p) for p in C])
        C=[np.hsplit(p,k) for p in np.vsplit(C,k)]

        Q=np.hstack(Q)
        Q=np.hsplit(Q,k)

        self.A={(i,j): C[i][j] for i in range(4) for j in range(4)}
        self.Q={i: Q[i] for i in range(4)}

#------------------------------------------------------------------------------
    def mrgh_(self,dtn1,dtn2):

        A=dtn1.A
        B=dtn2.A

        n=dtn1.n
        m=dtn2.n

        self.n=np.array([n[0]+m[0],m[1],n[2]+m[2],n[3]])

        B[3,3]=np.matmul(B[3,3],self.Tba)
        B[3,3]=np.matmul(self.Tab,B[3,3])

        B[3,0]=np.matmul(self.Tab,B[3,0])
        B[3,1]=np.matmul(self.Tab,B[3,1])
        B[3,2]=np.matmul(self.Tab,B[3,2])

        B[0,3]=np.matmul(B[0,3],self.Tba)
        B[1,3]=np.matmul(B[1,3],self.Tba)
        B[2,3]=np.matmul(B[2,3],self.Tba)

        self.S=np.linalg.solve(A[1,1]-B[3,3],np.eye(n[1]))

        Q=[-A[1,0],B[3,0],B[3,1],-A[1,2],B[3,2],-A[1,3]]

        Q=[np.matmul(self.S,q) for q in Q]

        C=[0,0,0,0,0,0]

        C[0]=[np.matmul(A[0,1],q) for q in Q]
        C[1]=[np.matmul(B[0,3],q) for q in Q]
        C[2]=[np.matmul(B[1,3],q) for q in Q]
        C[3]=[np.matmul(A[2,1],q) for q in Q]
        C[4]=[np.matmul(B[2,3],q) for q in Q]
        C[5]=[np.matmul(A[3,1],q) for q in Q]

        for [i,j] in [[0,0],[3,2],[5,3]]:

            C[i][0]=C[i][0]+A[j,0]
            C[i][3]=C[i][3]+A[j,2]
            C[i][5]=C[i][5]+A[j,3]

        for [i,j] in [[1,0],[2,1],[4,2]]:

            C[i][1]=C[i][1]+B[j,0]
            C[i][2]=C[i][2]+B[j,1]
            C[i][4]=C[i][4]+B[j,2]

        k=np.add.accumulate(self.n[0:-1])

        C=np.vstack([np.hstack(p) for p in C])
        C=[np.hsplit(p,k) for p in np.vsplit(C,k)]

        Q=np.hstack(Q)
        Q=np.hsplit(Q,k)

        self.A={(i,j): C[i][j] for i in range(4) for j in range(4)}
        self.Q={i: Q[i] for i in range(4)}

###############################################################################
class psn_elm:
#------------------------------------------------------------------------------
    def __init__(self,x,h,elm):
#------------------------------------------------------------------------------

        n=elm.n

        g=h[0]/x[...,None]

        Ae=elm.ope.Ly+elm.ope.Lx+g*elm.ope.Ux
        Ab=elm.opb.Ly+elm.opb.Lx+g*elm.opb.Ux

        Ae=np.linalg.solve(Ae,np.eye(Ae.shape[0]))

        G=np.matmul(Ae,-Ab)

        A=np.matmul(elm.Ne,G)+elm.Nb

        self.i=0

        self.n=np.array([n,n,n,n])

        self.Q=np.vstack([np.matmul(elm.ope.Q,G)+elm.opb.Q,
                          np.matmul(elm.ope.Ux,G)+elm.opb.Ux,
                          np.matmul(elm.ope.Uy,G)+elm.opb.Uy])

        self.R=np.vstack([np.matmul(elm.ope.Q,Ae),
                          np.matmul(elm.ope.Ux,Ae),
                          np.matmul(elm.ope.Uy,Ae)])

        self.A={(i,j): A[n*i:n*(i+1),n*j:n*(j+1)]/h[0] 
                for i in range(4) for j in range(4)}

        self.b=np.vsplit(np.matmul(elm.Ne,Ae)/h[0],4)

###############################################################################
class psn_rhs:
#------------------------------------------------------------------------------
    def __init__(self,dtn):
#------------------------------------------------------------------------------

        n=dtn.n

        self.n=n

        self.b=[np.zeros(n[0]),np.zeros(n[1]),
                np.zeros(n[2]),np.zeros(n[3])]

        self.f=[np.zeros(n[0]),np.zeros(n[1]),
                np.zeros(n[2]),np.zeros(n[3])]

        self.r=None
        self.q=None

        if dtn.i==1:

            m=dtn.S.shape[0]

            self.r=np.zeros(m)

            self.q=[np.zeros(m),np.zeros(m),np.zeros(m),np.zeros(m)]

#------------------------------------------------------------------------------
    def mrgh(self,rhs1,rhs2,dtn1,dtn2,dtn):

        k0=rhs1.n[0]
        k2=rhs1.n[2]

        a=rhs1.b
        b=rhs2.b

        A=dtn1.A
        B=dtn2.A

        np.dot(dtn.Tab,b[3],out=self.r)

        np.subtract(self.r,a[1],out=self.r)

        np.matmul(dtn.S,self.r,out=self.r)

        np.matmul(B[1,3],self.r,out=self.b[1])
        np.matmul(A[3,1],self.r,out=self.b[3])

        np.add(self.b[1],b[1],out=self.b[1])
        np.add(self.b[3],a[3],out=self.b[3])

        np.matmul(A[0,1],self.r,out=self.b[0][0:k0])
        np.matmul(B[0,3],self.r,out=self.b[0][k0::])

        np.add(self.b[0][0:k0],a[0],out=self.b[0][0:k0])
        np.add(self.b[0][k0::],b[0],out=self.b[0][k0::])

        np.matmul(A[2,1],self.r,out=self.b[2][0:k2])
        np.matmul(B[2,3],self.r,out=self.b[2][k2::])

        np.add(self.b[2][0:k2],a[2],out=self.b[2][0:k2])
        np.add(self.b[2][k2::],b[2],out=self.b[2][k2::])

#------------------------------------------------------------------------------
    def mrgv(self,rhs1,rhs2,dtn1,dtn2,dtn):

        k1=rhs1.n[1]
        k3=rhs1.n[3]

        a=rhs1.b
        b=rhs2.b

        A=dtn1.A
        B=dtn2.A

        np.dot(dtn.Tab,b[0],out=self.r)

        np.subtract(self.r,a[2],out=self.r)

        np.matmul(dtn.S,self.r,out=self.r)

        np.matmul(A[0,2],self.r,out=self.b[0])
        np.matmul(B[2,0],self.r,out=self.b[2])

        np.add(self.b[0],a[0],out=self.b[0])
        np.add(self.b[2],b[2],out=self.b[2])

        np.matmul(A[1,2],self.r,out=self.b[1][0:k1])
        np.matmul(B[1,0],self.r,out=self.b[1][k1::])

        np.add(self.b[1][0:k1],a[1],out=self.b[1][0:k1])
        np.add(self.b[1][k1::],b[1],out=self.b[1][k1::])

        np.matmul(A[3,2],self.r,out=self.b[3][0:k3])
        np.matmul(B[3,0],self.r,out=self.b[3][k3::])

        np.add(self.b[3][0:k3],a[3],out=self.b[3][0:k3])
        np.add(self.b[3][k3::],b[3],out=self.b[3][k3::])

#------------------------------------------------------------------------------
    def solv(self,dtn,rhs1,rhs2):

        k1=rhs1.n[1]
        k3=rhs1.n[3]

        np.matmul(dtn.Q[0],self.f[0],out=self.q[0])
        np.matmul(dtn.Q[1],self.f[1],out=self.q[1])
        np.matmul(dtn.Q[2],self.f[2],out=self.q[2])
        np.matmul(dtn.Q[3],self.f[3],out=self.q[3])

        np.add(self.q[0],self.q[1],out=self.q[1])
        np.add(self.q[2],self.q[3],out=self.q[3])
        np.add(self.q[1],self.q[3],out=self.q[1])

        np.add(self.q[1],self.r,out=self.r)

        rhs1.f[0][:]=self.f[0][:]
        rhs2.f[2][:]=self.f[2][:]

        rhs1.f[2][:]=self.r[:]

        np.dot(dtn.Tba,self.r,out=rhs2.f[0])

        rhs1.f[1][:]=self.f[1][0:k1]
        rhs2.f[1][:]=self.f[1][k1::]

        rhs1.f[3][:]=self.f[3][0:k3]
        rhs2.f[3][:]=self.f[3][k3::]

#------------------------------------------------------------------------------
    def solh(self,dtn,rhs1,rhs2):

        k0=rhs1.n[0]
        k2=rhs1.n[2]

        np.matmul(dtn.Q[0],self.f[0],out=self.q[0])
        np.matmul(dtn.Q[1],self.f[1],out=self.q[1])
        np.matmul(dtn.Q[2],self.f[2],out=self.q[2])
        np.matmul(dtn.Q[3],self.f[3],out=self.q[3])

        np.add(self.q[0],self.q[1],out=self.q[1])
        np.add(self.q[2],self.q[3],out=self.q[3])
        np.add(self.q[1],self.q[3],out=self.q[1])

        np.add(self.q[1],self.r,out=self.r)

        rhs1.f[1][:]=self.r[:]

        np.dot(dtn.Tba,self.r,out=rhs2.f[3])

        rhs1.f[3][:]=self.f[3][:]
        rhs2.f[1][:]=self.f[1][:]

        rhs1.f[0][:]=self.f[0][0:k0]
        rhs2.f[0][:]=self.f[0][k0::]

        rhs1.f[2][:]=self.f[2][0:k2]
        rhs2.f[2][:]=self.f[2][k2::]

###############################################################################
class psn_rhs_:
#------------------------------------------------------------------------------
    def __init__(self,dtn,dim):
#------------------------------------------------------------------------------

        self.n=dtn.n

        self.r=None
        self.q=None

        self.b=[None,None,None,None]
        self.f=[None,None,None,None]

        for i in range(4):

            self.b[i]=np.zeros((dtn.n[i],dim))

            self.f[i]=np.zeros((dtn.n[i],dim))

        if dtn.i==1:

            n=dtn.S.shape[0]

            self.r=np.zeros((n,dim))

            self.q=[np.zeros((n,dim)),np.zeros((n,dim)),
                    np.zeros((n,dim)),np.zeros((n,dim))]

#------------------------------------------------------------------------------
    def mrgh(self,rhs1,rhs2,dtn1,dtn2,dtn):

        k=rhs1.n[0]

        a=rhs1.b
        b=rhs2.b

        A=dtn1.A
        B=dtn2.A

        np.subtract(b[3],a[1],out=self.r)

        np.matmul(dtn.S,self.r,out=self.r)

        np.matmul(B[1,3],self.r,out=self.b[1])
        np.matmul(A[3,1],self.r,out=self.b[3])

        np.add(self.b[1],b[1],out=self.b[1])
        np.add(self.b[3],a[3],out=self.b[3])

        np.matmul(A[0,1],self.r,out=self.b[0][0:k,:])
        np.matmul(B[0,3],self.r,out=self.b[0][k::,:])

        np.add(self.b[0][0:k,:],a[0],out=self.b[0][0:k,:])
        np.add(self.b[0][k::,:],b[0],out=self.b[0][k::,:])

        np.matmul(A[2,1],self.r,out=self.b[2][0:k,:])
        np.matmul(B[2,3],self.r,out=self.b[2][k::,:])

        np.add(self.b[2][0:k,:],a[2],out=self.b[2][0:k,:])
        np.add(self.b[2][k::,:],b[2],out=self.b[2][k::,:])

#------------------------------------------------------------------------------
    def mrgv(self,rhs1,rhs2,dtn1,dtn2,dtn):

        k=rhs1.n[1]

        a=rhs1.b
        b=rhs2.b

        A=dtn1.A
        B=dtn2.A

        np.subtract(b[0][:,1::2],a[2][:,0::2],out=self.r)

        np.matmul(dtn.S,self.r,out=self.r)

        np.matmul(A[0,2],self.r,out=self.b[0])
        np.matmul(B[2,0],self.r,out=self.b[2])

        np.add(self.b[0],a[0][:,0::2],out=self.b[0])
        np.add(self.b[2],b[2][:,1::2],out=self.b[2])

        np.matmul(A[1,2],self.r,out=self.b[1][0:k,:])
        np.matmul(B[1,0],self.r,out=self.b[1][k::,:])

        np.add(self.b[1][0:k,:],a[1][:,0::2],out=self.b[1][0:k,:])
        np.add(self.b[1][k::,:],b[1][:,1::2],out=self.b[1][k::,:])

        np.matmul(A[3,2],self.r,out=self.b[3][0:k,:])
        np.matmul(B[3,0],self.r,out=self.b[3][k::,:])

        np.add(self.b[3][0:k,:],a[3][:,0::2],out=self.b[3][0:k,:])
        np.add(self.b[3][k::,:],b[3][:,1::2],out=self.b[3][k::,:])

#------------------------------------------------------------------------------
    def solv(self,dtn,rhs1,rhs2):

        k=rhs1.n[1]

        np.matmul(dtn.Q[0],self.f[0],out=self.q[0])
        np.matmul(dtn.Q[1],self.f[1],out=self.q[1])
        np.matmul(dtn.Q[2],self.f[2],out=self.q[2])
        np.matmul(dtn.Q[3],self.f[3],out=self.q[3])

        np.add(self.q[0],self.q[1],out=self.q[1])
        np.add(self.q[2],self.q[3],out=self.q[3])
        np.add(self.q[1],self.q[3],out=self.q[1])

        np.add(self.q[1],self.r,out=self.r)

        rhs1.f[0][:,0::2]=self.f[0][:,:]
        rhs2.f[2][:,1::2]=self.f[2][:,:]

        rhs1.f[2][:,0::2]=self.r[:,:]
        rhs2.f[0][:,1::2]=self.r[:,:]

        rhs1.f[1][:,0::2]=self.f[1][0:k,:]
        rhs2.f[1][:,1::2]=self.f[1][k::,:]

        rhs1.f[3][:,0::2]=self.f[3][0:k,:]
        rhs2.f[3][:,1::2]=self.f[3][k::,:]

#------------------------------------------------------------------------------
    def solh(self,dtn,rhs1,rhs2):

        k=rhs1.n[0]

        np.matmul(dtn.Q[0],self.f[0],out=self.q[0])
        np.matmul(dtn.Q[1],self.f[1],out=self.q[1])
        np.matmul(dtn.Q[2],self.f[2],out=self.q[2])
        np.matmul(dtn.Q[3],self.f[3],out=self.q[3])

        np.add(self.q[0],self.q[1],out=self.q[1])
        np.add(self.q[2],self.q[3],out=self.q[3])
        np.add(self.q[1],self.q[3],out=self.q[1])

        np.add(self.q[1],self.r,out=self.r)

        rhs1.f[1][:,:]=self.r[:,:]        
        rhs2.f[3][:,:]=self.r[:,:]        

        rhs1.f[3][:,:]=self.f[3][:,:]
        rhs2.f[1][:,:]=self.f[1][:,:]

        rhs1.f[0][:,:]=self.f[0][0:k,:]
        rhs2.f[0][:,:]=self.f[0][k::,:]

        rhs1.f[2][:,:]=self.f[2][0:k,:]
        rhs2.f[2][:,:]=self.f[2][k::,:]

###############################################################################
class psn_box:
#------------------------------------------------------------------------------
    def __init__(self,nod,elm):
#------------------------------------------------------------------------------

        n=nod.n

        mx=nod.m[1]
        my=nod.m[0]

        sx=int(np.log2(mx))+1
        sy=int(np.log2(my))

        self.n=nod.n
        self.m=nod.m

        self.k=np.arange(4)*(n*n)

        self.h=(nod.h[0]*nod.h[0],1./nod.h[0],1./nod.h[1])

        i_=np.split(np.arange(2**sx-1),
                    np.add.accumulate(2**np.arange(sx-1,0,-1)))

        j_=np.arange(sy)+i_[-1][0]+1

        self.e=np.arange(mx)

        self.i=list()
        self.j=list()

        for k in range(1,len(i_)):
            for i in range(len(i_[k])):
                self.i.append((i_[k][i],i_[k-1][2*i],i_[k-1][2*i+1]))

        self.j.append((j_[0],i_[-1][0],i_[-1][0]))

        for k in range(1,len(j_)):
            self.j.append((j_[k],j_[k-1],j_[k-1]))

        self.i_=self.i[::-1]
        self.j_=self.j[::-1]

        self.dtn=[None for i in range(2**sx-1+sy)]

        for i in self.e:
            self.dtn[i]=psn_elm(nod.x[:,i],nod.h,elm)

        for i in self.i:
            self.dtn[i[0]]=psn_dtn()
            self.dtn[i[0]].mrgh(self.dtn[i[1]],self.dtn[i[2]])

        for j in self.j:
            self.dtn[j[0]]=psn_dtn()
            self.dtn[j[0]].mrgv(self.dtn[j[1]],self.dtn[j[2]])

        self.rhs=[None for i in range(2**sx-1+sy)]

        for i in self.e:
            self.rhs[i]=psn_rhs_(self.dtn[i],my)

        for i in self.i:
            self.rhs[i[0]]=psn_rhs_(self.dtn[i[0]],my)

        ky=my
        for j in self.j:
            ky=ky/2
            self.rhs[j[0]]=psn_rhs_(self.dtn[j[0]],ky)

        self.q=np.zeros((1*n*n,my))
        self.g=np.zeros((3*n*n,my))

        self.r=np.zeros_like(nod.x)

        self.f=[np.zeros_like(nod.x),
                np.zeros_like(nod.x),
                np.zeros_like(nod.x)]

#------------------------------------------------------------------------------
    def up(self):
#------------------------------------------------------------------------------

        m=self.m[1]

        np.multiply(self.h[0],self.r,out=self.r)

        for i in self.e:

            self.q[:]=self.r[:,i::m]

            np.matmul(self.dtn[i].b[0],self.q,out=self.rhs[i].b[0])
            np.matmul(self.dtn[i].b[1],self.q,out=self.rhs[i].b[1])
            np.matmul(self.dtn[i].b[2],self.q,out=self.rhs[i].b[2])
            np.matmul(self.dtn[i].b[3],self.q,out=self.rhs[i].b[3])

        for i in self.i:

            self.rhs[i[0]].mrgh(self.rhs[i[1]],self.rhs[i[2]],
                                self.dtn[i[1]],self.dtn[i[2]],
                                self.dtn[i[0]])

        for i in self.j:

            self.rhs[i[0]].mrgv(self.rhs[i[1]],self.rhs[i[2]],
                                self.dtn[i[1]],self.dtn[i[2]],
                                self.dtn[i[0]])

#------------------------------------------------------------------------------
    def down(self):
#------------------------------------------------------------------------------

        m=self.m[1]
        k=self.k

        for i in self.j_:

            self.rhs[i[0]].solv(self.dtn[i[0]],
                                self.rhs[i[1]],self.rhs[i[2]])

        for i in self.i_:

            self.rhs[i[0]].solh(self.dtn[i[0]],
                                self.rhs[i[1]],self.rhs[i[2]])

        for i in self.e:

            dtn=self.dtn[i]
            rhs=self.rhs[i]

            np.matmul(dtn.Q,
                      np.vstack([rhs.f[0],rhs.f[1],rhs.f[2],rhs.f[3]]),
                      out=self.g)

            self.q[:]=self.r[:,i::m]

            self.g=self.g+np.matmul(self.dtn[i].R,self.q)

            self.f[0][:,i::m]=self.g[k[0]:k[1],:]
            self.f[1][:,i::m]=self.g[k[1]:k[2],:]
            self.f[2][:,i::m]=self.g[k[2]:k[3],:]

        np.multiply(self.h[1],self.f[1],out=self.f[1])
        np.multiply(self.h[2],self.f[2],out=self.f[2])

#------------------------------------------------------------------------------
    def geterr(self,err,lgn):
#------------------------------------------------------------------------------

        L1=np.amax(abs(err))

        L2=np.matmul(lgn.C2,err)

        L2=np.sum(np.sum(L2*L2,axis=0))*self.h[0]

        L2=np.sqrt(L2)

        return L1,L2

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
    def geteps(self,lgn):

        E=1809.0*self.f[2]+15.

        eps=np.matmul(lgn.Uy,E)*self.h[1]

        eps=np.amax(abs(eps))

        return eps

###############################################################################