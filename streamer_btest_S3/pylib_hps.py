###############################################################################
" pylib_hps "
###############################################################################
import pylib_lgn as lg
import pylib_psn as ps
import pylib_fvs as fv
import pylib_dgs as dg
#------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pickle
###############################################################################
def getKi(E):

    Ki=(1.64616*E)*np.exp(-197.6/E)

    return Ki

def mask_dgs_fvs(elm):

    flag=True

    return flag

###############################################################################
class sol:
    def __init__(self):

        self.x=None
        self.y=None

        self.Ne=None
        self.Ni=None

        self.Ex=None
        self.Ey=None

        self.E=None

        self.y_fld_axs=None
        self.y_lgn_axs=None

        self.Ne_axs=None
        self.Ni_axs=None

        self.E_axs=None
        self.U_axs=None
        self.h_axs=None

    def rho(self):
        plt.plot(self.y_axs,self.Ne_axs,"r-")
        plt.plot(self.y_axs,self.Ni_axs,"b-")

    def rho_log(self):
        plt.semilogy(self.y_axs,self.Ne_axs,"r-")
        plt.semilogy(self.y_axs,self.Ni_axs,"b-")

    def field(self):
        plt.plot(self.y_axs,self.E_axs,"-b")

    def plot2d(self,dat):

        fmin=np.amin(dat[~np.isnan(dat)])
        fmax=np.amax(dat[~np.isnan(dat)])

        lev=np.linspace(fmin,fmax,51)

        plt.contourf(+self.x,self.y,dat,levels=lev,cmap="jet")
        plt.contourf(-self.x,self.y,dat,levels=lev,cmap="jet")

        plt.axis("equal")

###############################################################################
#   class unit
##############################################################################
class unit:
#------------------------------------------------------------------------------
    def __init__(self,c,h):

        self.t=0.

        self.dim=0

        self.top=knot(c,h,[0,0,0])

        self.elm=[]

        self.rib_fvs=[]
        self.rib_dgs=[]

        self.brd=None
        self.edg=None
        self.ebd=None

        self.flx=np.array([0.,0.,0.,0.])

        self.sol=sol()

#------------------------------------------------------------------------------
    def add(self,mask):
         self.top.add(self.top,mask)

#------------------------------------------------------------------------------
    def fix(self):

        flag=False

        self.getelm()

        edg=[np.array([],dtype=int),np.array([],dtype=int)]

        self.top.getedg(self.top,edg,self.dim)

        i=np.where(abs(edg[0][1::2]-edg[1][1::2])>1)[0]

        if len(i):

            p=edg[0][2*i+1]-edg[1][2*i+1]

            j1=np.where(p<0)[0]
            j2=np.where(p>0)[0]

            p=np.append(edg[0][2*i[j1]],edg[1][2*i[j2]])

            for p_ in np.unique(p):
                self.elm[p_].split()

            flag=True

        if flag:
            self.getelm()

        return flag

#------------------------------------------------------------------------------
    def fix_amr(self,fvs,dgs,grn,lgn):

        flag=False

        self.getelm()

        edg=[np.array([],dtype=int),np.array([],dtype=int)]

        self.top.getedg(self.top,edg,self.dim)

        i=np.where(abs(edg[0][1::2]-edg[1][1::2])>1)[0]

        if len(i):

            p=edg[0][2*i+1]-edg[1][2*i+1]

            j1=np.where(p<0)[0]
            j2=np.where(p>0)[0]

            p=np.append(edg[0][2*i[j1]],edg[1][2*i[j2]])

            for p_ in np.unique(p):
                self.elm[p_].split_amr(fvs,dgs,grn,lgn)

            flag=True

        if flag:
            self.getelm()

        return flag

#------------------------------------------------------------------------------

    def setnod_lgn(self,lgn):
        self.top.setnod_lgn(self.top,lgn)

    def setpsn(self,grn):
        self.top.setpsn(self.top,grn)

    def setfld_fvs(self,fvs):
        self.top.setfld_fvs(self.top,fvs)

    def setfld_dgs(self,dgs):
        self.top.setfld_dgs(self.top,dgs)

#------------------------------------------------------------------------------

    def psndtn(self,grn):

        self.top.psndtn(self.top,grn)

        A1=np.hstack([self.top.dtn.A[1,1],self.top.dtn.A[1,3]])
        A3=np.hstack([self.top.dtn.A[3,1],self.top.dtn.A[3,3]])

        A=np.vstack([A1,A3])

        A=np.linalg.solve(A,np.eye(A.shape[0]))

        self.top.dtn.B=[None,None,None,None]

        n=self.top.dtn.n[1]

        self.top.dtn.B[0]=A[0:n,0:n]
        self.top.dtn.B[1]=A[0:n,n::]
        self.top.dtn.B[2]=A[n::,0:n]
        self.top.dtn.B[3]=A[n::,n::]

    def psnrhs(self):

        self.top.psnrhs(self.top)

        self.top.rhs.f[0][:]=0.
        self.top.rhs.f[2][:]=0.

        f1=np.matmul(self.top.dtn.B[0],self.top.rhs.b[1])
        f2=np.matmul(self.top.dtn.B[1],self.top.rhs.b[3])

        self.top.rhs.f[1][:]=-f1-f2

        f1=np.matmul(self.top.dtn.B[2],self.top.rhs.b[1])
        f2=np.matmul(self.top.dtn.B[3],self.top.rhs.b[3])

        self.top.rhs.f[3][:]=-f1-f2

    def psnsol(self):
        self.top.psnsol(self.top)

#------------------------------------------------------------------------------
    def getelm(self):

        self.popelm()

        self.top.getelm(self.top,self.elm)

        self.dim=max([elm_.p[2] for elm_ in self.elm])

    def popelm(self):
        for i in range(len(self.elm)):
            self.elm.pop()

#------------------------------------------------------------------------------
    def getrib(self):

        self.getelm()

        rib=[[],[],[],[],[],[]]

        rib[0]=np.array([[],[]],dtype=int)
        rib[1]=np.array([[],[]],dtype=int)

        rib[2]=np.array([[],[],[]],dtype=int)
        rib[3]=np.array([[],[],[]],dtype=int)
        rib[4]=np.array([[],[],[]],dtype=int)
        rib[5]=np.array([[],[],[]],dtype=int)

        self.top.getrib(self.top,rib,self.dim)

        m=len(self.elm)

        ind=[np.ones(m,dtype=int),np.ones(m,dtype=int),
             np.ones(m,dtype=int),np.ones(m,dtype=int)]

        ind[2][rib[0][0,:]]=0
        ind[0][rib[0][1,:]]=0

        ind[1][rib[1][0,:]]=0
        ind[3][rib[1][1,:]]=0

        ind[0][rib[2][0,:]]=0
        ind[0][rib[2][1,:]]=0
        ind[2][rib[2][2,:]]=0

        ind[2][rib[3][0,:]]=0
        ind[2][rib[3][1,:]]=0
        ind[0][rib[3][2,:]]=0

        ind[3][rib[4][0,:]]=0
        ind[3][rib[4][1,:]]=0
        ind[1][rib[4][2,:]]=0

        ind[1][rib[5][0,:]]=0
        ind[1][rib[5][1,:]]=0
        ind[3][rib[5][2,:]]=0

        brd=[np.where(ind[0]==1)[0],np.where(ind[1]==1)[0],
             np.where(ind[2]==1)[0],np.where(ind[3]==1)[0]]

        self.brd=[[],[],[],[]]

        for i in [0,1,2,3]:
            for e in brd[i]:
                self.brd[i].append(self.elm[e].j)

        return rib

#------------------------------------------------------------------------------
    def fvsrib(self):

        rib=self.getrib()

        self.poprib_fvs()

        for r in range(rib[0].shape[1]):

            i=rib[0][0,r]
            j=rib[0][1,r]

            box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs]

            self.rib_fvs.append(fv.fvs_rib(box,0))

        for r in range(rib[1].shape[1]):

            i=rib[1][0,r]
            j=rib[1][1,r]

            box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs]

            self.rib_fvs.append(fv.fvs_rib(box,1))

        for r in range(rib[2].shape[1]):

            i=rib[2][0,r]
            j=rib[2][1,r]
            k=rib[2][2,r]

            box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,self.elm[k].fld_fvs]

            self.rib_fvs.append(fv.fvs_rib(box,2))

        for r in range(rib[3].shape[1]):

            i=rib[3][0,r]
            j=rib[3][1,r]
            k=rib[3][2,r]

            box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,self.elm[k].fld_fvs]

            self.rib_fvs.append(fv.fvs_rib(box,3))

        for r in range(rib[4].shape[1]):

            i=rib[4][0,r]
            j=rib[4][1,r]
            k=rib[4][2,r]

            box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,self.elm[k].fld_fvs]

            self.rib_fvs.append(fv.fvs_rib(box,4))

        for r in range(rib[5].shape[1]):

            i=rib[5][0,r]
            j=rib[5][1,r]
            k=rib[5][2,r]

            box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,self.elm[k].fld_fvs]

            self.rib_fvs.append(fv.fvs_rib(box,5))

#------------------------------------------------------------------------------
    def dgsrib(self):

        rib=self.getrib()

        self.poprib_dgs()

        for r in range(rib[0].shape[1]):

            i=rib[0][0,r]
            j=rib[0][1,r]

            box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs]

            self.rib_dgs.append(dg.dgs_rib(box,0))

        for r in range(rib[1].shape[1]):

            i=rib[1][0,r]
            j=rib[1][1,r]

            box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs]

            self.rib_dgs.append(dg.dgs_rib(box,1))

        for r in range(rib[2].shape[1]):

            i=rib[2][0,r]
            j=rib[2][1,r]
            k=rib[2][2,r]

            box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,self.elm[k].fld_dgs]

            self.rib_dgs.append(dg.dgs_rib(box,2))

        for r in range(rib[3].shape[1]):

            i=rib[3][0,r]
            j=rib[3][1,r]
            k=rib[3][2,r]

            box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,self.elm[k].fld_dgs]

            self.rib_dgs.append(dg.dgs_rib(box,3))

        for r in range(rib[4].shape[1]):

            i=rib[4][0,r]
            j=rib[4][1,r]
            k=rib[4][2,r]

            box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,self.elm[k].fld_dgs]

            self.rib_dgs.append(dg.dgs_rib(box,4))

        for r in range(rib[5].shape[1]):

            i=rib[5][0,r]
            j=rib[5][1,r]
            k=rib[5][2,r]

            box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,self.elm[k].fld_dgs]

            self.rib_dgs.append(dg.dgs_rib(box,5))

#------------------------------------------------------------------------------
    def rib_dgs_fvs(self):

        rib=self.getrib()

        self.poprib_fvs()
        self.poprib_dgs()

        for r in range(rib[0].shape[1]):

            i=rib[0][0,r]
            j=rib[0][1,r]

            mask1=mask_dgs_fvs(self.elm[i])
            mask2=mask_dgs_fvs(self.elm[j])

            mask=[mask1,mask2]

            if mask==[True,True]:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,0))
            elif mask==[False,False]:
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,0))
            else:                
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,0))
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,0))

        for r in range(rib[1].shape[1]):

            i=rib[1][0,r]
            j=rib[1][1,r]

            mask1=mask_dgs_fvs(self.elm[i])
            mask2=mask_dgs_fvs(self.elm[j])

            mask=[mask1,mask2]

            if mask==[True,True]:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,1))
            elif mask==[False,False]:
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,1))
            else:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,1))        
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,1))

        for r in range(rib[2].shape[1]):

            i=rib[2][0,r]
            j=rib[2][1,r]
            k=rib[2][2,r]

            mask1=mask_dgs_fvs(self.elm[i])
            mask2=mask_dgs_fvs(self.elm[j])
            mask3=mask_dgs_fvs(self.elm[k])

            mask=[mask1,mask2,mask3]

            if mask==[True,True,True]:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,2))
            elif mask==[False,False,False]:                
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,2))
            else:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,2))
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,2))         

        for r in range(rib[3].shape[1]):

            i=rib[3][0,r]
            j=rib[3][1,r]
            k=rib[3][2,r]

            mask1=mask_dgs_fvs(self.elm[i])
            mask2=mask_dgs_fvs(self.elm[j])
            mask3=mask_dgs_fvs(self.elm[k])

            mask=[mask1,mask2,mask3]

            if mask==[True,True,True]:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,3))
            elif mask==[False,False,False]:                
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,3))
            else:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,3))
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,3))

        for r in range(rib[4].shape[1]):

            i=rib[4][0,r]
            j=rib[4][1,r]
            k=rib[4][2,r]

            mask1=mask_dgs_fvs(self.elm[i])
            mask2=mask_dgs_fvs(self.elm[j])
            mask3=mask_dgs_fvs(self.elm[k])

            mask=[mask1,mask2,mask3]

            if mask==[True,True,True]:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,4))
            elif mask==[False,False,False]:                
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,4))
            else:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,4))
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,4))

        for r in range(rib[5].shape[1]):

            i=rib[5][0,r]
            j=rib[5][1,r]
            k=rib[5][2,r]

            mask1=mask_dgs_fvs(self.elm[i])
            mask2=mask_dgs_fvs(self.elm[j])
            mask3=mask_dgs_fvs(self.elm[k])

            mask=[mask1,mask2,mask3]

            if mask==[True,True,True]:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,5))
            elif mask==[False,False,False]:                
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,5))
            else:
                box=[self.elm[i].fld_fvs,self.elm[j].fld_fvs,
                     self.elm[k].fld_fvs]
                self.rib_fvs.append(fv.fvs_rib(box,5))
                box=[self.elm[i].fld_dgs,self.elm[j].fld_dgs,
                     self.elm[k].fld_dgs]
                self.rib_dgs.append(dg.dgs_rib(box,5))

#------------------------------------------------------------------------------

    def poprib_fvs(self):
        for i in range(len(self.rib_fvs)):
            self.rib_fvs.pop()

    def poprib_dgs(self):
        for i in range(len(self.rib_dgs)):
            self.rib_dgs.pop()

#------------------------------------------------------------------------------
    def setedg(self):

        self.edg=[(),(),(),()]
        
        for i in [0,1,2,3]:

            x_=[]
            y_=[]

            for j in self.brd[i]:

                x,y=self.elm[j].nod_lgn.getedg()

                x_.append(x[i])
                y_.append(y[i])

            self.edg[i]=(np.hstack(x_),np.hstack(y_))

#------------------------------------------------------------------------------
    def getsol(self,l,lgn,fvs):

        sol=self.sol

        n=fvs.s
        m=fvs.m

        dim=self.dim

        hx=self.top.h[0]/(2**dim)
        hy=self.top.h[1]/(2**dim)

        imax=[]

        for elm in self.elm:
            if elm.p[2]>=l:

                [i,j,p]=elm.p

                k=2**(dim-p)

                imax.append(k*(i+1))

        imax=max(imax)
        jmax=2**dim

        ky=jmax
        kx=imax

        cx=self.top.c[0]-self.top.h[0]
        cy=self.top.c[1]-self.top.h[1]

        hx_=hx/m
        hy_=hy/m

        cx=np.arange(kx*m)*(2.*hx_)+cx+hx_
        cy=np.arange(ky*m)*(2.*hx_)+cy+hy_

        x=fvs.x[...,None]*hx_+cx[None,...]
        y=fvs.x[...,None]*hy_+cy[None,...]

        x=x.T.flatten()
        y=y.T.flatten()

        [x,y]=np.meshgrid(x,y)

        sol.x=x
        sol.y=y

        nm=n*m

        k=(ky*nm,kx*nm)

        sol.Ne=np.zeros(k)*np.nan
        sol.Ni=np.zeros(k)*np.nan

        sol.Ex=np.zeros(k)*np.nan
        sol.Ey=np.zeros(k)*np.nan

        for elm in self.elm:
            if elm.p[2]>=l:

                [i,j,p]=elm.p

                k=2**(dim-p)

                i=k*i
                j=k*j

                i=np.arange(nm*i,nm*(i+k))
                j=np.arange(nm*j,nm*(j+k))

                i,j=np.meshgrid(i,j)

                Ne=elm.fld_fvs.f[0]
                Ni=elm.fld_fvs.f[2]

                if p==dim:
                    sol.Ne[j,i]=Ne.copy()
                    sol.Ni[j,i]=Ni.copy()
                else:
                    sol.Ne[j,i]=fvs.refine(Ne,dim-p)
                    sol.Ni[j,i]=fvs.refine(Ni,dim-p)

                if elm.psn:                    

                    Ex=-1809.*elm.psn.f[1]
                    Ey=-1809.*elm.psn.f[2]-52.

                    Ex=fvs.g2f(Ex)
                    Ey=fvs.g2f(Ey)
                    
                    if p==dim:
                        sol.Ex[j,i]=Ex.copy()
                        sol.Ey[j,i]=Ey.copy()
                    else:
                        sol.Ex[j,i]=fvs.refine(Ex,dim-p)
                        sol.Ey[j,i]=fvs.refine(Ey,dim-p)

        sol.E=np.sqrt(sol.Ex*sol.Ex+sol.Ey*sol.Ey)

        sol.y_axs=[]

        sol.Ne_axs=[]
        sol.Ni_axs=[]

        sol.E_axs=[]

        for j in self.brd[3]:

            elm=self.elm[j]            

            sol.y_axs.append(elm.fld_fvs.y[:,0])

            sol.Ne_axs.append(elm.fld_fvs.f[0][:,0])
            sol.Ni_axs.append(elm.fld_fvs.f[2][:,0])

            E=fvs.g2f(elm.psn.f[2])

            E=-1809.*E[:,0]-52.

            sol.E_axs.append(E)

        sol.y_axs=np.hstack(sol.y_axs)

        sol.Ne_axs=np.hstack(sol.Ne_axs)
        sol.Ni_axs=np.hstack(sol.Ni_axs)

        sol.E_axs=np.hstack(sol.E_axs)

        return self.sol

#------------------------------------------------------------------------------
    def plot(self):

        self.getelm()

        for elm_ in self.elm:
            elm_.plot()
        plt.axis("equal")

#------------------------------------------------------------------------------
    def save(self,path):

        self.getelm()

        Ne=[]
        Ni=[]
        for elm in self.elm:
            Ne.append(elm.fld_fvs.f[0])
            Ni.append(elm.fld_fvs.f[2])

        pickle.dump([Ne,Ni],open(path+"/fld_fvs.bin","wb"))

        Ne=[]
        Ni=[]
        for elm in self.elm:
            Ne.append(elm.fld_dgs.f[0])
            Ni.append(elm.fld_dgs.f[2])

        pickle.dump([Ne,Ni],open(path+"/fld_dgs.bin","wb"))

        self.top.free(self.top)

        self.popelm()

        self.poprib_fvs()
        self.poprib_dgs()

        self.sol=None

        pickle.dump(self,open(path+"/hps.bin","wb"))

    def load(self,path,dgs,fvs,lgn):

        self.setnod_lgn(lgn)

        self.setfld_fvs(fvs)
        self.setfld_dgs(dgs)

        self.getelm()

        [Ne,Ni]=pickle.load(open(path+"/fld_fvs.bin","rb"))

        Ne=Ne[::-1]
        Ni=Ni[::-1]
        for elm in self.elm:
            elm.fld_fvs.f[0][:,:]=Ne.pop()
            elm.fld_fvs.f[2][:,:]=Ni.pop()

        [Ne,Ni]=pickle.load(open(path+"/fld_dgs.bin","rb"))

        Ne=Ne[::-1]
        Ni=Ni[::-1]
        for elm in self.elm:
            elm.fld_dgs.f[0][:,:]=Ne.pop()
            elm.fld_dgs.f[2][:,:]=Ni.pop()

        self.sol=sol()

#------------------------------------------------------------------------------
    def step_psn(self,lgn):

        for elm in self.elm:
            elm.psn.r=elm.fld_dgs.f[0]-elm.fld_dgs.f[2]

        self.psnrhs()
        self.psnsol()

#------------------------------------------------------------------------------
    def step_dgs_fvs(self,dt,nt,dgs,fvs,lgn):

        for j in range(nt):

            for elm in self.elm:
                elm.psn.r=elm.fld_dgs.f[0]-elm.fld_dgs.f[2]

            self.psnrhs()
            self.psnsol()

            for elm in self.elm:

                psn=elm.psn

                fld=elm.fld_dgs

                fld.f[4][:,:]=-1809.0*psn.f[1][:,:]
                fld.f[5][:,:]=-1809.0*psn.f[2][:,:]-52.

                fld.f[6]=np.sqrt(fld.f[4]*fld.f[4]+fld.f[5]*fld.f[5])

                fld.ax[:,:]=-0.038*fld.f[4]
                fld.ay[:,:]=-0.038*fld.f[5]

                fld.f[7]=getKi(fld.f[6])

                fld.edge(0,dgs)

                fld=elm.fld_fvs

                fld.f[4][:,:]=-1809.*fvs.g2f(psn.f[1])
                fld.f[5][:,:]=-1809.*fvs.g2f(psn.f[2])-52.

                fld.f[6]=np.sqrt(fld.f[4]*fld.f[4]+fld.f[5]*fld.f[5])

                fld.ax[:,:]=-0.038*fld.f[4]
                fld.ay[:,:]=-0.038*fld.f[5]

                fld.f[7]=getKi(fld.f[6])

                fld.edge(0)

#            self.flx=self.flx+0.5*dt*self.flux(lgn)

            for rib in self.rib_dgs:
                rib.flux(dgs)
            for rib in self.rib_fvs:
                rib.flux()

            for elm in self.elm:

                mask=mask_dgs_fvs(elm)

                if not mask:

                    fld=elm.fld_dgs

                    fld.delta(0,dgs)

                    fld.df[1]=fld.f[7]*fld.f[0] 

                    fld.f[1]=fld.f[0]+dt*fld.df[1]+dt*fld.df[0]
                    fld.f[3]=fld.f[2]+dt*fld.df[1]

                    elm.fld_fvs.f[1][:,:]=fvs.g2f(elm.fld_dgs.f[1])
                    elm.fld_fvs.f[3][:,:]=fvs.g2f(elm.fld_dgs.f[3])

                else:

                    fld=elm.fld_fvs

                    fld.delta(0)

                    fld.df[1]=fld.f[7]*fld.f[0] 

                    fld.f[1]=fld.f[0]+dt*fld.df[1]+dt*fld.df[0]
                    fld.f[3]=fld.f[2]+dt*fld.df[1]

                    elm.fld_dgs.f[1]=fvs.f2g(elm.fld_fvs.f[1])
                    elm.fld_dgs.f[3]=fvs.f2g(elm.fld_fvs.f[3])

            for elm in self.elm:
                elm.psn.r=elm.fld_dgs.f[1]-elm.fld_dgs.f[3]

            self.psnrhs()
            self.psnsol()

            for elm in self.elm:

                psn=elm.psn

                fld=elm.fld_dgs

                fld.f[4][:,:]=-1809.0*psn.f[1][:,:]
                fld.f[5][:,:]=-1809.0*psn.f[2][:,:]-52.

                fld.f[6]=np.sqrt(fld.f[4]*fld.f[4]+fld.f[5]*fld.f[5])

                fld.ax[:,:]=-0.038*fld.f[4]
                fld.ay[:,:]=-0.038*fld.f[5]

                fld.f[7]=getKi(fld.f[6])

                fld.edge(1,dgs)

                fld=elm.fld_fvs

                fld.f[4][:,:]=-1809.*fvs.g2f(psn.f[1])
                fld.f[5][:,:]=-1809.*fvs.g2f(psn.f[2])-52.

                fld.f[6]=np.sqrt(fld.f[4]*fld.f[4]+fld.f[5]*fld.f[5])

                fld.ax[:,:]=-0.038*fld.f[4]
                fld.ay[:,:]=-0.038*fld.f[5]

                fld.f[7]=getKi(fld.f[6])

                fld.edge(1)

#            for elm in self.elm:
#                elm.fld_dgs.edge(1,dgs)
#                elm.fld_fvs.edge(1)

#            self.flx=self.flx+0.5*dt*self.flux(lgn)

            for rib in self.rib_dgs:
                rib.flux(dgs)
            for rib in self.rib_fvs:
                rib.flux()

            for elm in self.elm:

                mask=mask_dgs_fvs(elm)

                if not mask:

                    fld=elm.fld_dgs

                    fld.delta(1,dgs)

                    fld.df[1]=fld.f[7]*fld.f[1]

                    fld.f[0]=0.5*(fld.f[0]+fld.f[1]+dt*fld.df[1]+dt*fld.df[0])
                    fld.f[2]=0.5*(fld.f[2]+fld.f[3]+dt*fld.df[1])

                    elm.fld_fvs.f[0][:,:]=fvs.g2f(elm.fld_dgs.f[0])
                    elm.fld_fvs.f[2][:,:]=fvs.g2f(elm.fld_dgs.f[2])

                else:

                    fld=elm.fld_fvs
    
                    fld.delta(1)

                    fld.df[1]=fld.f[7]*fld.f[1]

                    fld.f[0]=0.5*(fld.f[0]+fld.f[1]+dt*fld.df[1]+dt*fld.df[0])
                    fld.f[2]=0.5*(fld.f[2]+fld.f[3]+dt*fld.df[1])

                    elm.fld_dgs.f[0]=fvs.f2g(elm.fld_fvs.f[0])
                    elm.fld_dgs.f[2]=fvs.f2g(elm.fld_fvs.f[2])

        self.t=self.t+nt*dt

    def flux(self,lgn):

        n=lgn.n
        
        f0=0.
        f2=0.

        for i in self.brd[0]:
            fld=self.elm[i].fld_dgs

            x=self.elm[i].nod_lgn.x[0:n,0:1]
            f0=f0+np.sum(lgn.lgn.I[...,None]*x*fld.e[0].F)*fld.h[0]

#            f0=f0+np.sum(lgn.lgn.I[...,None]*fld.e[0].F)*fld.h[0]

        for i in self.brd[2]:
            fld=self.elm[i].fld_dgs

            x=self.elm[i].nod_lgn.x[0:n,0:1]
            f2=f2+np.sum(lgn.lgn.I[...,None]*x*fld.e[2].F)*fld.h[0]

#            f2=f2+np.sum(lgn.lgn.I[...,None]*fld.e[2].F)*fld.h[0]

        return np.array([f0,0.,f2,0.])

    def charge(self,lgn):

        rho_=0.0
        for elm in self.elm:

            fld=elm.fld_dgs
            nod=elm.nod_lgn

            rho=(fld.f[2]-fld.f[0])*nod.x

            rho=np.sum(lgn.I*rho)*(fld.h[0]*fld.h[0])

            rho_=rho_+rho

        return rho_

###############################################################################
#   class knot
###############################################################################
class knot:
#------------------------------------------------------------------------------
    def __init__(self,c,h,p):
#------------------------------------------------------------------------------

        self.flg=False
        self.key=False

        self.lft=None
        self.rgt=None

        self.c=c
        self.h=h
        self.p=p

        self.j=0

        self.t=None

        self.e=[None,None,None,None]

        self.nod_lgn=None

        self.psn=None

        self.fld_fvs=None
        self.fld_dgs=None

        self.dtn=None
        self.rhs=None

        self.eps1=0.
        self.eps2=0.

#------------------------------------------------------------------------------
    def split(self):

        self.flg=1

        [cx,cy]=self.c
        [hx,hy]=self.h
        [i,j,p]=self.p

        hx=0.5*hx

        self.lft=knot([cx-hx,cy],[hx,hy],[None,None,p])
        self.rgt=knot([cx+hx,cy],[hx,hy],[None,None,p])

        self.lft.flg=2
        self.rgt.flg=2

        hy=0.5*hy

        self.lft.lft=knot([cx-hx,cy-hy],[hx,hy],[2*i+0,2*j+0,p+1])
        self.lft.rgt=knot([cx-hx,cy+hy],[hx,hy],[2*i+0,2*j+1,p+1])

        self.rgt.lft=knot([cx+hx,cy-hy],[hx,hy],[2*i+1,2*j+0,p+1])
        self.rgt.rgt=knot([cx+hx,cy+hy],[hx,hy],[2*i+1,2*j+1,p+1])

#------------------------------------------------------------------------------
    def add(self,root,mask):
        if root:
            self.add(root.lft,mask)
            self.add(root.rgt,mask)
            if not root.flg:
                if mask(root):
                    root.split()

    def sub(self,root):
        if root:
            self.sub(root.lft)
            self.sub(root.rgt)
            if root.key:

                del root.lft
                del root.rgt

                root.lft=None
                root.rgt=None

                root.flg=False

    def keyup(self,root):
        if root:
            self.keyup(root.lft)
            self.keyup(root.rgt)

            root.key=False

            if root.flg==1:

                flag=[root.lft.lft.flg,root.lft.rgt.flg,
                      root.rgt.lft.flg,root.rgt.rgt.flg]

                if flag==[False,False,False,False]:
                    root.key=True

#------------------------------------------------------------------------------
    def geteps(self,root):
        if root:
            self.geteps(root.lft)
            self.geteps(root.rgt)
            if not root.flg:

                eps1,eps2=root.fld_fvs.geteps()
                
                root.eps1=eps1
                root.eps2=eps2

#------------------------------------------------------------------------------
    def free(self,root):
        if root:
            self.free(root.lft)
            self.free(root.rgt)

            root.t=None

            root.e=[None,None,None,None]

            root.nod_lgn=None

            root.fld_fvs=None;
            root.fld_dgs=None;

            root.psn=None;
            root.dtn=None; root.rhs=None;

#------------------------------------------------------------------------------
    def getelm(self,root,elm):
        if root:
            self.getelm(root.lft,elm)
            self.getelm(root.rgt,elm)
            if not root.flg:
                elm.append(root)      
                root.j=len(elm)-1

#------------------------------------------------------------------------------
    def getedg(self,root,edg,dim):
        if root:
            self.getedg(root.lft,edg,dim)
            self.getedg(root.rgt,edg,dim)

            if not root.flg:

                odr=2**(dim-root.p[2])

                edg_=np.array([root.j,root.p[2]]*odr,dtype=int)

                root.e[0]=edg_
                root.e[1]=edg_
                root.e[2]=edg_
                root.e[3]=edg_

            else:

                if root.flg==2:

                    root.e[0]=root.lft.e[0]
                    root.e[2]=root.rgt.e[2]

                    root.e[1]=np.append(root.lft.e[1],root.rgt.e[1])
                    root.e[3]=np.append(root.lft.e[3],root.rgt.e[3])

                    edg[0]=np.append(edg[0],root.lft.e[2])
                    edg[1]=np.append(edg[1],root.rgt.e[0])

                else:

                    root.e[0]=np.append(root.lft.e[0],root.rgt.e[0])
                    root.e[2]=np.append(root.lft.e[2],root.rgt.e[2])

                    root.e[1]=root.rgt.e[1]
                    root.e[3]=root.lft.e[3]

                    edg[0]=np.append(edg[0],root.lft.e[1])
                    edg[1]=np.append(edg[1],root.rgt.e[3])

#------------------------------------------------------------------------------
    def getrib(self,root,rib,dim):
        if root:
            self.getrib(root.lft,rib,dim)
            self.getrib(root.rgt,rib,dim)
            if not root.flg:

                odr=2**(dim-root.p[2])

                rib_=np.array([root.j,odr],dtype=int)

                root.e[0]=rib_
                root.e[1]=rib_
                root.e[2]=rib_
                root.e[3]=rib_

            else:

                if root.flg==2:

                    root.e[0]=root.lft.e[0]
                    root.e[2]=root.rgt.e[2]

                    root.e[1]=np.append(root.lft.e[1],root.rgt.e[1])
                    root.e[3]=np.append(root.lft.e[3],root.rgt.e[3])

                    p,t=root.getpar(root.lft.e[2],root.rgt.e[0])

                    rib[0]=np.append(rib[0],p[2],axis=1)
                    rib[2]=np.append(rib[2],p[0],axis=1)
                    rib[3]=np.append(rib[3],p[1],axis=1)

                    root.t=t

                else:

                    root.e[0]=np.append(root.lft.e[0],root.rgt.e[0])
                    root.e[2]=np.append(root.lft.e[2],root.rgt.e[2])

                    root.e[1]=root.rgt.e[1]
                    root.e[3]=root.lft.e[3]

                    p,t=root.getpar(root.lft.e[1],root.rgt.e[3])

                    rib[1]=np.append(rib[1],p[2],axis=1)
                    rib[4]=np.append(rib[4],p[0],axis=1)
                    rib[5]=np.append(rib[5],p[1],axis=1)

                    root.t=t

#------------------------------------------------------------------------------
    def getpar(self,a,b):

        p0=a[1::2]
        p1=b[1::2]

        m0=len(p0)
        m1=len(p1)

        n0=a[0::2].copy()
        n1=b[0::2].copy()

        m_=min(np.amin(p0),np.amin(p1))

        p0=p0/m_
        p1=p1/m_

        m_=np.sum(p0)

        T=np.zeros((m_,m1),dtype=int)
        S=np.zeros((m0,m_),dtype=int)

        i=np.arange(m_); j=np.repeat(np.arange(m1),p1); T[i,j]=1;
        i=np.repeat(np.arange(m0),p0); j=np.arange(m_); S[i,j]=1;

        Q=S.dot(T)

        ij=np.argwhere(Q)

        i=np.unique(ij[:,0],return_index=True,return_counts=True)
        j=np.unique(ij[:,1],return_index=True,return_counts=True)

        q1=np.array([[],[],[]],dtype=int)
        q2=np.array([[],[],[]],dtype=int)

        t1=np.ones(ij.shape[0],dtype=int)

        if np.sum(i[2])-m0:

            k=i[1][i[2]==2]

            t1[k]=3
            t1[k+1]=3

            k1=ij[k,0]
            k2=ij[k,1]
            k3=ij[k+1,1]

            q1=np.vstack([n1[k2],n1[k3],n0[k1]])

            n0[k1]=-1
            n1[k2]=-1
            n1[k3]=-1

        if np.sum(j[2])-m1:

            k=j[1][j[2]==2]

            t1[k]=2
            t1[k+1]=2

            k1=ij[k,1]
            k2=ij[k,0]
            k3=ij[k+1,0]

            q2=np.vstack([n0[k2],n0[k3],n1[k1]])

            n1[k1]=-1
            n0[k2]=-1
            n0[k3]=-1

        q3=np.vstack([n0[n0>-1],n1[n1>-1]])

        t1[np.where((t1==3)|(t1==2))[0][0::2]]=0
        t1=t1[t1!=0]

        t2=t1.copy()

        t2[t1==3]=2
        t2[t1==2]=3

        t=np.array([],dtype=int)

        if max(abs(t1-t2)):
            t=np.vstack([t1-1,t2-1])

        return (q1,q2,q3),t

#------------------------------------------------------------------------------

    def setnod_lgn(self,root,lgn):
        if root:
            self.setnod_lgn(root.lft,lgn)
            self.setnod_lgn(root.rgt,lgn)
            if not root.flg:
                root.nod_lgn=lg.nod2d(root.c,root.h,lgn)

    def setpsn(self,root,grn):
        if root:
            self.setpsn(root.lft,grn)
            self.setpsn(root.rgt,grn)
            if not root.flg:
                root.psn=ps.psn_box(root.nod_lgn,grn)

    def setfld_fvs(self,root,fvs):
        if root:
            self.setfld_fvs(root.lft,fvs)
            self.setfld_fvs(root.rgt,fvs)
            if not root.flg:
                root.fld_fvs=fv.fvs_box2d(root.c,root.h,fvs)

    def setfld_dgs(self,root,dgs):
        if root:
            self.setfld_dgs(root.lft,dgs)
            self.setfld_dgs(root.rgt,dgs)
            if not root.flg:
                root.fld_dgs=dg.dgs_box2d(root.nod_lgn)

#------------------------------------------------------------------------------
    def psndtn(self,root,grn):
        if root:
            self.psndtn(root.lft,grn)
            self.psndtn(root.rgt,grn)
            if not root.flg:

                root.dtn=ps.psn_dtn()

                root.dtn.i=0

                root.dtn.n=root.psn.dtn[-1].n

                root.dtn.A={}
                for (i,j) in root.psn.dtn[-1].A.keys():
                    root.dtn.A[i,j]=root.psn.dtn[-1].A[i,j].copy()
 
                root.rhs=ps.psn_rhs(root.dtn)

            else:
               
                root.dtn=ps.psn_dtn()

                if root.flg==2:
                    if root.t.shape[0]:

                        Tab=[grn.T[q] for q in root.t[0,:]]
                        Tba=[grn.T[q] for q in root.t[1,:]]

                        root.dtn.Tab=lg.bdiag(*Tab)
                        root.dtn.Tba=lg.bdiag(*Tba)

                        root.dtn.mrgv_(root.lft.dtn,root.rgt.dtn)

                        root.rhs=ps.psn_rhs(root.dtn)

                    else:

                        root.dtn.Tab=1.
                        root.dtn.Tba=1.

                        root.dtn.mrgv(root.lft.dtn,root.rgt.dtn)

                        root.rhs=ps.psn_rhs(root.dtn)

                else:
                    if root.t.shape[0]:

                        Tab=[grn.T[q] for q in root.t[0,:]]
                        Tba=[grn.T[q] for q in root.t[1,:]]

                        root.dtn.Tab=lg.bdiag(*Tab)
                        root.dtn.Tba=lg.bdiag(*Tba)

                        root.dtn.mrgh_(root.lft.dtn,root.rgt.dtn)

                        root.rhs=ps.psn_rhs(root.dtn)

                    else:

                        root.dtn.Tab=1.
                        root.dtn.Tba=1.

                        root.dtn.mrgh(root.lft.dtn,root.rgt.dtn)

                        root.rhs=ps.psn_rhs(root.dtn)

#------------------------------------------------------------------------------
    def psnrhs(self,root):
        if root:
            self.psnrhs(root.lft)
            self.psnrhs(root.rgt)
            if not root.flg:

                root.psn.up()

                root.rhs.b[0]=root.psn.rhs[-1].b[0][:,0]
                root.rhs.b[1]=root.psn.rhs[-1].b[1][:,0]
                root.rhs.b[2]=root.psn.rhs[-1].b[2][:,0]
                root.rhs.b[3]=root.psn.rhs[-1].b[3][:,0]

            else:

                if root.flg==2:
                   root.rhs.mrgv(root.lft.rhs,root.rgt.rhs,
                                 root.lft.dtn,root.rgt.dtn,root.dtn)
                else:
                    root.rhs.mrgh(root.lft.rhs,root.rgt.rhs,
                                  root.lft.dtn,root.rgt.dtn,root.dtn)

#------------------------------------------------------------------------------
    def psnsol(self,root):
        if root:
            if root.flg:
                if root.flg==2:
                    root.rhs.solv(root.dtn,root.lft.rhs,root.rgt.rhs)
                else:
                    root.rhs.solh(root.dtn,root.lft.rhs,root.rgt.rhs)
            else:

                root.psn.rhs[-1].f[0][:,0]=root.rhs.f[0]
                root.psn.rhs[-1].f[1][:,0]=root.rhs.f[1]
                root.psn.rhs[-1].f[2][:,0]=root.rhs.f[2]
                root.psn.rhs[-1].f[3][:,0]=root.rhs.f[3]

                root.psn.down()

            self.psnsol(root.lft)
            self.psnsol(root.rgt)

#------------------------------------------------------------------------------
    def split_amr(self,fvs,dgs,grn,lgn):

        self.split()

        elm=[self.lft.lft,self.lft.rgt,
             self.rgt.lft,self.rgt.rgt]

        for i in [0,1,2,3]:

            elm[i].nod_lgn=lg.nod2d(elm[i].c,elm[i].h,lgn)

            elm[i].fld_dgs=dg.dgs_box2d(elm[i].nod_lgn)

            elm[i].fld_fvs=fv.fvs_box2d(elm[i].c,elm[i].h,fvs)

            elm[i].psn=ps.psn_box(elm[i].nod_lgn,grn)

#        f=self.fld_fvs.split(0,fvs)
#        elm[0].fld_fvs.f[0][:,:]=f[0][:,:]
#        elm[1].fld_fvs.f[0][:,:]=f[1][:,:]
#        elm[2].fld_fvs.f[0][:,:]=f[2][:,:]
#        elm[3].fld_fvs.f[0][:,:]=f[3][:,:]
#        f=self.fld_fvs.split(2,fvs)
#        elm[0].fld_fvs.f[2][:,:]=f[0][:,:]
#        elm[1].fld_fvs.f[2][:,:]=f[1][:,:]
#        elm[2].fld_fvs.f[2][:,:]=f[2][:,:]
#        elm[3].fld_fvs.f[2][:,:]=f[3][:,:]

        f=self.fld_dgs.split(0,dgs)

        elm[0].fld_dgs.f[0][:,:]=f[0][:,:]
        elm[1].fld_dgs.f[0][:,:]=f[1][:,:]
        elm[2].fld_dgs.f[0][:,:]=f[2][:,:]
        elm[3].fld_dgs.f[0][:,:]=f[3][:,:]

        f=self.fld_dgs.split(2,dgs)

        elm[0].fld_dgs.f[2][:,:]=f[0][:,:]
        elm[1].fld_dgs.f[2][:,:]=f[1][:,:]
        elm[2].fld_dgs.f[2][:,:]=f[2][:,:]
        elm[3].fld_dgs.f[2][:,:]=f[3][:,:]

        elm[0].fld_fvs.f[0]=fvs.g2f(elm[0].fld_dgs.f[0])
        elm[1].fld_fvs.f[0]=fvs.g2f(elm[1].fld_dgs.f[0])
        elm[2].fld_fvs.f[0]=fvs.g2f(elm[2].fld_dgs.f[0])
        elm[3].fld_fvs.f[0]=fvs.g2f(elm[3].fld_dgs.f[0])

        elm[0].fld_fvs.f[2]=fvs.g2f(elm[0].fld_dgs.f[2])
        elm[1].fld_fvs.f[2]=fvs.g2f(elm[1].fld_dgs.f[2])
        elm[2].fld_fvs.f[2]=fvs.g2f(elm[2].fld_dgs.f[2])
        elm[3].fld_fvs.f[2]=fvs.g2f(elm[3].fld_dgs.f[2])

        elm[0].eps1=self.eps1*0.5
        elm[1].eps1=self.eps1*0.5
        elm[2].eps1=self.eps1*0.5
        elm[3].eps1=self.eps1*0.5

        elm[0].eps2=self.eps2*0.5
        elm[1].eps2=self.eps2*0.5
        elm[2].eps2=self.eps2*0.5
        elm[3].eps2=self.eps2*0.5

        self.nod_lgn=None

        self.fld_fvs=None
        self.fld_dgs=None

        self.psn=None

    def merge_amr(self,fvs,dgs,grn,lgn):

        self.nod_lgn=lg.nod2d(self.c,self.h,lgn)

        self.fld_fvs=fv.fvs_box2d(self.c,self.h,fvs)

        self.fld_dgs=dg.dgs_box2d(self.nod_lgn)

        self.psn=ps.psn_box(self.nod_lgn,grn)

#        f=[self.lft.lft.fld_fvs.f[0],self.lft.rgt.fld_fvs.f[0],
#           self.rgt.lft.fld_fvs.f[0],self.rgt.rgt.fld_fvs.f[0]]
#        self.fld_fvs.merge(0,f,fvs)
#        f=[self.lft.lft.fld_fvs.f[2],self.lft.rgt.fld_fvs.f[2],
#           self.rgt.lft.fld_fvs.f[2],self.rgt.rgt.fld_fvs.f[2]]
#        self.fld_fvs.merge(2,f,fvs)

        f=[self.lft.lft.fld_dgs.f[0],self.lft.rgt.fld_dgs.f[0],
           self.rgt.lft.fld_dgs.f[0],self.rgt.rgt.fld_dgs.f[0]]

        self.fld_dgs.merge(0,f,dgs)

        f=[self.lft.lft.fld_dgs.f[2],self.lft.rgt.fld_dgs.f[2],
           self.rgt.lft.fld_dgs.f[2],self.rgt.rgt.fld_dgs.f[2]]

        self.fld_dgs.merge(2,f,dgs)

        self.fld_fvs.f[0]=fvs.g2f(self.fld_dgs.f[0])
        self.fld_fvs.f[2]=fvs.g2f(self.fld_dgs.f[2])

        self.eps1=np.amax([self.lft.lft.eps1,self.lft.rgt.eps1,
                           self.rgt.lft.eps1,self.rgt.rgt.eps1])

        self.eps2=np.amax([self.lft.lft.eps2,self.lft.rgt.eps2,
                           self.rgt.lft.eps2,self.rgt.rgt.eps2])

        del self.lft
        del self.rgt

        self.lft=None
        self.rgt=None

        self.flg=False

#------------------------------------------------------------------------------
    def add_amr(self,root,mask,fvs,dgs,grn,lgn):
        if root:
            self.add_amr(root.lft,mask,fvs,dgs,grn,lgn)
            self.add_amr(root.rgt,mask,fvs,dgs,grn,lgn)
            if not root.flg:
                if mask(root):
                    root.split_amr(fvs,dgs,grn,lgn)

    def sub_amr(self,root,mask,fvs,dgs,grn,lgn):
        if root:
            self.sub_amr(root.lft,mask,fvs,dgs,grn,lgn)
            self.sub_amr(root.rgt,mask,fvs,dgs,grn,lgn)
            if root.key:

                root.eps1=max([root.lft.lft.eps1,root.lft.rgt.eps1,
                               root.rgt.lft.eps1,root.rgt.rgt.eps1])

                root.eps2=max([root.lft.lft.eps2,root.lft.rgt.eps2,
                               root.rgt.lft.eps2,root.rgt.rgt.eps2])

                if mask(root):
                    root.merge_amr(fvs,dgs,grn,lgn)

#------------------------------------------------------------------------------
    def plot(self):

        cx=self.c[0]
        cy=self.c[1]

        hx=self.h[0]
        hy=self.h[1]

        x=np.array([cx-hx,cx+hx,cx+hx,cx-hx,cx-hx])
        y=np.array([cy-hy,cy-hy,cy+hy,cy+hy,cy-hy])

        plt.plot(x,y,"-b")

        plt.text(cx,cy,str(self.j))

###############################################################################
