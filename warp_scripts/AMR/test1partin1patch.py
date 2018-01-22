from warp import *
from warp.field_solvers.AMR import *
top.nhist=100000

l_mr = 0   # switch to turn mesh refinement ON

ntransit=0 # number of transition (ghost) cells for reduction of spurious self-forces with MR

if not l_mr:
  bfact=2
else:
  bfact=1

EnableAll()

w3d.solvergeom=w3d.XZgeom
top.pboundxy = reflect

w3d.nz = w3d.nx = 64*bfact
w3d.zmmin=w3d.xmmin=-1.
w3d.zmmax=w3d.xmmax= 1.
w3d.dx=(w3d.xmmax-w3d.xmmin)/w3d.nx

# --- add electron
elec = Species(type=Electron,weight=1.e8)
elec.addpart(-6*w3d.dx*bfact,0,0.,0.,0.,0.,lallindomain=1)

top.dt = bfact*100*w3d.dx/clight

MRroot=MRBlock2D()
registersolver(MRroot)
if l_mr:
   MRroot.addchild(mins=[-0.25,0.,-0.25],maxs=[0.25,0.,0.25],nguard=[ntransit,0,ntransit])

winon()
package('w3d');generate()
# --- set shortcuts
pg = top.pgroup

loadrho()
fieldsol()

hx = AppendableArray()
def addhx():
    hx.append(getx()[0])
    if top.it%50==0:
        fma()
        MRroot.pfzx(contours=0,cellarray=1);MRroot.drawboxzx(withguards=1)
        ppzx(color=red,msize=5,titles=0)
        refresh()
  
installafterstep(addhx)

addhx()
step(2000)
winon(1);pla(hx[:]/w3d.dx+32);ptitles('X vs time','Time step','')

