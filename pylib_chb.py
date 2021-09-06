###############################################################################
" pylib_chb "
###############################################################################
import numpy as np
from scipy.linalg import block_diag as bdiag
###############################################################################
def permat(n,m):

    P=np.zeros((n*m,n*m))

    for j in range(m):
        for i in range(n):
            P[i*m+j,i+j*n]=1

    return P

###############################################################################
class chb1d:

    def __init__(self,n):

        t=(np.arange(n)/(n-1.)+1.)*np.pi

        x=np.cos(t)

        C=np.cos(np.outer(np.arange(n),t)).T

        C=np.linalg.solve(C,np.eye(n))

        P=np.zeros((n,n))
        U=np.zeros((n,n))

        P[:,0]=1.
        P[:,1]=2.*x

        for i in range(1,n-1):
            P[:,i+1]=2.*x*P[:,i]-P[:,i-1]

        U[:,1::]=P[:,0:-1]

        U=(U*np.arange(n)).dot(C)

        self.x=x

        self.C=C
        self.U=U

        self.L=np.matmul(U,U)

###############################################################################
class chb2d:

    def __init__(self,n):

        chb=chb1d(n)

        x=chb.x

        U=chb.U
        L=chb.L

        P=permat(n,n)

        self.x=x

        Ux=bdiag(*[U for i in range(n)])
        Uy=np.matmul(P,np.matmul(Ux,P))

        Lx=bdiag(*[L for i in range(n)])
        Ly=np.matmul(P,np.matmul(Lx,P))

        self.Ux=Ux
        self.Uy=Uy

        self.Lx=Lx
        self.Ly=Ly

###############################################################################
