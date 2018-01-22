from warp import *
from warp.lattice import *
from warp.field_solvers.AMR import *
from warp.envelope.env_match import *

case = 'lowres'
#case = 'highres'
#case = 'AMR'

if case=='lowres':
    l_amr   = 0
    nx = ny = 64
if case=='highres':
    l_amr   = 0
    nx = ny = 64*8
if case=='AMR': # uses 3 levels of refinement
    l_amr   = 1
    nx = ny = 64*2

l_test           = true
l_match          = true
l_plots          = 1
w3d.l4symtry     = 1
top.npmax        = 50000
freqplot         = 50#0

frz.l_get_fields_on_grid=0

winon()

w3d.solvergeom=w3d.XYgeom


# AMR parameters
if l_amr:
    w3d.AMRlevels             = 3
    w3d.AMRmaxlevel_density   = 1
    w3d.AMRmaxlevel_gradient  = w3d.AMRlevels
    w3d.AMRthreshold_gradient = 0.2
    w3d.AMRgenerate_periodicity = 2
    w3d.AMRtransit              = 2
    w3d.gchange()
    w3d.AMRuse_inactive_regions=True
    w3d.AMRcoalescing=0.6

# --- Set four-character run id, comment lines, user's name.
top.runid    = "test"
top.pline2   = "Timer with moments and plots"
top.pline1   = "Slice beam. 128x128"
top.runmaker = "David P. Grote"
# --- Invoke setup routine 
numticks = 4
if not l_test:
    setup()
else:
    window(0)

# --- Set input parameters describing the beam, 72 to 17.
# --- Parameters calculated with envelope code ignoring focusing effect
# --- of dipoles.
top.a0       = 15.358933450767e-3
top.b0       =  8.6379155933081e-3
top.ap0      = 0.e0
top.bp0      = 0.e0
top.ibeam    = 2.e-03
top.emit     = 51.700897052724e-6*0.01
top.vbeam    = 0.e0
top.ekin     = 80.e3
top.aion     = 39.1e0
top.zion     = 1.e0
top.lrelativ = true
derivqty()
top.vthz     = .5e0*top.vbeam*top.emit/sqrt(top.a0*top.b0) # Vthz ~ Vthperp
# +++ Set up arrays describing lattice.
# --- Set temp variables.
hlp = 36.0e-2      # half lattice period length
piperad = 3.445e-2 # pipe radius
quadlen = 11.e-2   # quadrupole length
gaplen  = 3.e-2
rodlen  = quadlen + gaplen
rodrad  = 8./7.*piperad
platewid = 1.e-2
dbdx = .949/quadlen
# --- Set general lattice variables.
top.tunelen   = 2.e0*hlp
env.zl        = 0.0
env.zu        = 2.0*hlp
env.dzenv     = top.tunelen/100.e0

# --- Set up quadrupoles
# --- This uses the MAD-like lattice input mechanism
dhalf = Drft(l=(hlp-quadlen)/2.,ap=piperad)
f1 = Quad(l=quadlen,de=-dbdx*top.vbeam,ap=piperad,rr=rodrad,rl=rodlen,
        gl=gaplen,gp=1,pw=platewid)
d1 = Quad(l=quadlen,de=+dbdx*top.vbeam,ap=piperad,rr=rodrad,rl=rodlen,
        gl=gaplen,gp=1,pw=platewid)

lattice = dhalf + d1 + 2*dhalf + f1 + dhalf
madtowarp(lattice)

top.zlatperi = 2.*hlp
top.zlatstrt = hlp/2.

# +++ Set input parameters describing the 3d simulation.
w3d.nx = nx;w3d.ny = ny
top.ibpush = 2           # accurate B advance
steps_p_perd = 25
top.dt = (top.tunelen/steps_p_perd)/top.vbeam
top.prwall = piperad
# --- Set to finite beam.
w3d.xmmin = -piperad*15./7. #*(1. + 1./(w3d.nx/2-1))
w3d.xmmax =  piperad*15./7. #*(1. + 1./(w3d.nx/2-1))
w3d.ymmin = -piperad*15./7. #*(1. + 1./(w3d.ny/2-1))
w3d.ymmax =  piperad*15./7. #*(1. + 1./(w3d.ny/2-1))
# --- Load Semi-Gaussian cigar beam.
w3d.distrbtn = "semigaus"
w3d.xrandom = "digitrev"
w3d.vtrandom = "digitrev"
w3d.vzrandom = "digitrev"
w3d.ldprfile = "polar"
# --- Set up some windows.
top.rwindows[:,1] = [0.e0,.005e0]
top.rwindows[:,2] = [0.e0,.01e0]
top.rwindows[:,3] = [0.e0,.02e0]
# --- Select plot intervals, etc.
top.nhist = 1
top.itmomnts[0:4]=[0,1000000,abs(top.nhist),0]
top.itplseldom[0:4]=[0,100000,steps_p_perd,0]
top.itplalways[0:4]=[0,100000,1,0]

timeplots = 0.
def myplots(l_force=0):
  global timeplots
  if top.it%freqplot!=0 and not l_force:return
  tstart = wtime()
  fma()
  pfxy(contours=10,filled=1)
  ppxy(chopped=0.05)
  window(1);fma();ppxy(chopped=0.5)
  if l_amr:
    AMRtree.draw(allmesh=1)
  window(0)
  pyg_pending()
  pyg_idler()
  timeplots += wtime()-tstart

# --- Run the envelope solver to provide data used to initialize particles.
package("env"); generate(); step()
if l_match:match1(100)
penv();fma()

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.).
package("wxy"); generate()

top.plottime+=timeplots
generatefstime = top.fstime
generateplottime = top.plottime

if l_plots:installafterstep(myplots)

t0=wtime()
step(250)

print('50 steps in ',wtime()-t0,'s')

# --- Kludge to remove generate field-solve time from the total field-solve time.
top.fstime = top.fstime - generatefstime
top.plottime = top.plottime - generateplottime

printtimers()

print()
print("Time for 50 steps) = ",top.steptime)
print("Plot time per step = ",top.plottime/top.it)
print("Field solve time per step = ",top.fstime/top.it)
print()
print("Number of particles = ",top.pgroup.nps)
print("Mesh dimensions = ",w3d.nx,"x",w3d.ny)

n = 64

window(0);fma();ppxy(color='density',nx=n,ny=n,contours=20,filled=1,view=3)
palette('heat.gp')

rho = zeros([n+1,n+1],'d')
setgrid2d(getn(),getx(),gety(),n,n,rho,min(getx()),max(getx()),min(gety()),max(gety()))

ppxvx(view=4)

hpepsx(color=red,width=3,plsysval=10);limits(0.,top.time,0.,7.e-6)

window(1);myplots(True)
window(1);draw_mesh(w3d.nx,w3d.ny,w3d.xmmin,w3d.ymmin,w3d.dx,w3d.dy,color='red');ppxy()
if l_amr:AMRtree.draw(allmesh=1)
limits(-0.02,0.02,-0.012,0.012)
try:
  from Opyndx import *
  DXMountainPlot(rho,colorbar=0,xmin=w3d.xmmin,ymin=w3d.ymmin,dx=w3d.dx,dy=w3d.dy,l_interactive=0)
except:
  pass