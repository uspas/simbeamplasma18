from warp import *
from warp.field_solvers.em3dsolver import *
from warp.run_modes.boosted_frame import *
#import PRpickle as PR
#import PWpickle as PW
import os
home=os.getenv('HOME')
# --- flags turning off unnecessary diagnostics (ignore for now)
top.ifzmmnt = 0
top.itmomnts = 0
top.itplps = 0
top.itplfreq = 0
top.zzmomnts = 0
top.zzplps = 0
top.zzplfreq = 0
top.nhist = top.nt
top.iflabwn = 0
w3d.lrhodia3d = false
w3d.lgetese3d = false
w3d.lgtlchg3d = false

#-------------------------------------------------------------------------------
# main parameters
#-------------------------------------------------------------------------------
#dim = "3d"                 # 3D calculation
dim = "2d"                 # 2D calculation 
#dim = "1d"                 # 1D calculation 
dpi=100                     # graphics resolution
l_test             = 0      # Will open output window on screen 
                            # and stop before entering main loop.
l_gist             = 1      # Turns gist plotting on/off
l_restart          = false  # To restart simulation from an old run (works?)
restart_dump       = ""     # dump file to restart from (works?)
l_moving_window    = 1      # on/off (Galilean) moving window
l_plasma           = 1      # on/off plasma
l_ions             = 1      # on/off plasma ions
l_beam             = 1      # on/off electron beam
l_injectplane      = 1      # on/off beam injection through plane
l_usesavedist      = 0      # if on, uses dump of beam particles distribution
savedist           = home+'/runs/warp/lhc/quasistatic/unitdist4sym300000' # beam initial distribution file
svstride           = 100    # loads only every svstride particles from dump
l_smooth           = 1      # on/off smoothing of current density
l_laser            = 1      # on/off laser
l_pdump            = 0      # on/off regular dump of beam data
stencil            = 1      # 0 = Yee; 1 = Yee-enlarged (Karkkainen) on EF,B; 2 = Yee-enlarged (Karkkainen) on E,F 
                            # use 0 or 1; 2 does not verify Gauss Law
if dim=="1d":stencil=0
dtcoef             = 1.     # coefficient to multiply default time step that is set at the EM solver CFL
top.depos_order    = 3      # particles deposition order (1=linear, 2=quadratic, 3=cubic)
top.efetch         = 4      # field gather type (1=from nodes "momentum conserving"; 4=from Yee mesh "energy conserving")
nzstations         = 50     # number of beam diag z-stations
nzfieldstations    = 20     # number of field diag z-stations
l_pselect          = 0      # on/off selection of particles (i.e. remove halo) for diagnostics
top.runid          = "lpa_basic"                         # run name
top.pline1         = "basic lpa"                         # comment line on plots
top.runmaker       = "J.-L. Vay,"                        # run makers
top.lrelativ       = true                                # on/off relativity (for particles push)
top.pgroup.lebcancel_pusher=true                         # flag for particle pusher (0=Boris pusher; 1=Vay PoP 08 pusher)
l_verbose          = 0                                   # verbosity level (0=off; 1=on)
l_correct_num_Cherenkov = True                           # flag for correction to numerical Cherenkov instability (NCI)

#-------------------------------------------------------------------------------
# diagnostics parameters + a few other settings
#-------------------------------------------------------------------------------
hist_freq          = 50    # frequency (in time steps) of beam history data saving
live_plot_freq     = 20   # frequency (in time steps) of live plots (off is l_test is off)

#-------------------------------------------------------------------------------
# boosted frame
#-------------------------------------------------------------------------------
gammafrm           = 6.
betafrm            = sqrt(1.-1./gammafrm**2)
if gammafrm>1.: # turns ON plasma ions if in boosted frame
  l_ions = 1
  l_moving_window = 1

#-------------------------------------------------------------------------------
# some units for convenience (& clarity)
#-------------------------------------------------------------------------------
microns            = 1.e-6
femtoseconds       = 1.e-15

#-------------------------------------------------------------------------------
# plasma density and max nb particles per cell 
#-------------------------------------------------------------------------------
# --- in lab frame
dfact             = 100.                                     # coefficient factor for plasma density (for scaled simulations)
dens0lab          = dfact*1.e23                            # plasma density
wplab             = sqrt(dens0lab*echarge**2/(eps0*emass)) # plasma frequency
kplab             = wplab/clight                           # plasma wavenumber
lambda_plasma_lab = 2.*pi/kplab                            # plasma wavelength
length_pramp_lab  = 1.e-10*0.05/(dfact**1.5)                      # plasma entrance ramp length
length_pramp_exit_lab = 0.05/(dfact**1.5)                  # plasma exit ramp length
Lplasma_lab       = 1./dfact**1.5                          # plasma length
# --- in boosted frame
dens0             = dens0lab*gammafrm                      # plasma density
length_pramp      = length_pramp_lab/gammafrm              # plasma ramp length
length_pramp_exit = length_pramp_exit_lab/gammafrm         # plasma ramp length
Lplasma           = Lplasma_lab/gammafrm                   # plasma length

#-------------------------------------------------------------------------------
# laser parameters
#-------------------------------------------------------------------------------
# --- in lab frame
KP_L               = 2.0                     # normalized length - standard gaussian form P=P0*exp(-2xi^2/L^2)  
KP_SIGMA           = 5.3                     # normalized transverse spot size - standard gaussian form I=I0*exp(-2r^2/SIGMA^2)  
laser_radius       = KP_SIGMA / kplab         
laser_waist        = laser_radius*2.354820   # radius -> FWHM
laser_length_lab   = KP_L/kplab              # laser length 
laser_duration_lab = laser_length_lab/clight # laser duration (FWHM)
lambda_laser_lab   = 0.8*microns             # wavelength 
laser_polangle     = pi/2                    # polarization (0=aligned with x; pi/2=aligned with y)
a0                 = 1.                      # normalized potential vector (amplitude)
k0lab              = 2.*pi/lambda_laser_lab
w0lab              = k0lab*clight
ZR                 = 0.5*k0lab*(laser_waist**2)   # Rayleigh length
# --- in boosted frame
lambda_laser       = lambda_laser_lab*gammafrm*(1.+betafrm)   # wavelength 
laser_duration     = laser_duration_lab*gammafrm*(1.+betafrm)
laser_length       = laser_duration/clight
k0                 = 2.*pi/lambda_laser
w0                 = k0*clight
Eamp               = a0*w0*emass*clight/echarge
Bamp               = Eamp/clight 
if l_laser==0:Eamp=Bamp=0.

#-------------------------------------------------------------------------------
# plasma cont'd
#-------------------------------------------------------------------------------
K                   = k0lab/kplab
Ld=LINEAR_DEPHASING = 0.5*lambda_plasma_lab**3/lambda_laser_lab**2             # linear dephasing length
BETAG_LINEAR_LAB    = sqrt(1-(1./K)*(1./K))                                    # linear beta of wake in lab
GAMMAG_LINEAR_LAB   = 1./sqrt(1-BETAG_LINEAR_LAB*BETAG_LINEAR_LAB)             # linear gamma of wake in lab
BETAG_LINEAR        = (BETAG_LINEAR_LAB-betafrm)/(1.-BETAG_LINEAR_LAB*betafrm) # linear beta of wake in simulation frame
GAMMAG_LINEAR       = 1./sqrt(1-BETAG_LINEAR*BETAG_LINEAR)                     # linear gamma of wake in simulation frame
kp                  = kplab*(gammafrm*(1.-BETAG_LINEAR_LAB*betafrm))
lambda_plasma       = 2.*pi/kp
densc               = emass*eps0*w0lab**2/echarge**2                           # critical density

#-------------------------------------------------------------------------------
# print some plasma parameters to the screen
#-------------------------------------------------------------------------------
print("the laser group velocity is: ")
print((BETAG_LINEAR*clight))
print("the laser spot size is: ")
print(laser_waist)
print("the Rayleigh length is: ")
print(ZR)
print("the laser wavelength is: ")
print(lambda_laser)
print("the plasma wavelength is: ")
print(lambda_plasma)
print("the plasma length is: ")
print(Lplasma)

#-------------------------------------------------------------------------------
# e-beam
#-------------------------------------------------------------------------------
# --- in lab frame
E_BEAM_GAMMA      = GAMMAG_LINEAR_LAB*1.5
E_BEAM_ENERGY_MEV = 0.511*(E_BEAM_GAMMA-1.)
E_BEAM_BETA       = sqrt(1.- 1./(E_BEAM_GAMMA*E_BEAM_GAMMA))
E_BEAM_U          = E_BEAM_GAMMA * E_BEAM_BETA
E_BEAM_RADIUS     = 0.825e-6/sqrt(dfact)
E_BEAM_LENGTH     = 0.85e-6/sqrt(dfact)

# --- transverse spread (RMS Gaussian)
GAMMAVXSIGMA = 0.
GAMMAVYSIGMA = 0.
# --- longitudinal spread (RMS Gaussian)
GAMMAVZSIGMA = 0.

E_BEAM_DENSITY_PEAK = 1.0e10
E_BEAM_PHASE        = 5.*pi/4.
E_BEAM_DISTANCE_BEHIND_LASER = E_BEAM_PHASE/kplab

#-------------------------------------------------------------------------------
# number of grid cells
#-------------------------------------------------------------------------------
# --- transverse
nx = 50
# --- longitudinal
nzplambda = 32

#-------------------------------------------------------------------------------
# number of plasma macro-particles/cell
#-------------------------------------------------------------------------------
nppcellx = 1
nppcelly = 1
nppcellz = 1

if dim=="2d":
  nppcelly = 1
if dim=="1d":
  nppcellx = nppcelly = 1

#-------------------------------------------------------------------------------
# grid dimensions, nb cells and BC
#-------------------------------------------------------------------------------
w3d.xmmax = 3.*laser_radius
w3d.xmmin = -w3d.xmmax
w3d.ymmax = w3d.xmmax
w3d.ymmin = -w3d.ymmax
w3d.nx = nx
w3d.ny = w3d.nx    
w3d.dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
if dim in ["1d"]:
    w3d.nx = 2
    w3d.xmmin = -float(w3d.nx)/2
    w3d.xmmax = float(w3d.nx)/2
if dim in ["1d","2d"]:
    w3d.ny = 2
    w3d.ymmin = -float(w3d.ny)/2
    w3d.ymmax = float(w3d.ny)/2
w3d.zmmin = -2.*lambda_plasma
w3d.zmmax = 2.*lambda_laser#/nzplambda
w3d.nz = int((w3d.zmmax-w3d.zmmin)*nzplambda/lambda_laser) 
w3d.dx = (w3d.xmmax-w3d.xmmin)/w3d.nx
w3d.dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
w3d.dz = (w3d.zmmax-w3d.zmmin)/w3d.nz

# --- enforces dx = dz (= dy) in boosted frame
if 0:#gammafrm>1:
  if w3d.dz<w3d.dx:
    # x (and y) dimensions are rescaled if needed
    w3d.dx = w3d.dz
    w3d.nx = 2*nint(0.5*(w3d.xmmax-w3d.xmmin)/w3d.dx)
    w3d.xmmax = w3d.nx*w3d.dx/2
    w3d.xmmin = -w3d.xmmax
    if dim=='3d':
      w3d.dy = w3d.dz
      w3d.ny = 2*nint(0.5*(w3d.ymmax-w3d.ymmin)/w3d.dy)
      w3d.ymmax = w3d.ny*w3d.dy/2
      w3d.ymmin = -w3d.ymmax
  elif w3d.dx<w3d.dz:
    # z dimensions are rescaled if needed
    w3d.dz = w3d.dx
    w3d.nz = nint((w3d.zmmax-w3d.zmmin)/w3d.dz)
    w3d.zmmax = w3d.zmmin+w3d.nz*w3d.dz

if gammafrm>10.:
  if dim=="1d":
    dt0 = w3d.dz/clight
  if dim=="2d":
    dt0 = 1./(clight*sqrt(1./w3d.dx**2+1./w3d.dz**2))
  if dim=="3d":
    dt0 = 1./(clight*sqrt(1./w3d.dx**2+1./w3d.dy**2+1./w3d.dz**2))
  mydt = w3d.dz/(sqrt(2.)*clight)
  dtcoef = min(1.,mydt/dt0)

# --- sets field boundary conditions
w3d.bound0  = w3d.boundnz = openbc
#w3d.bound0  = w3d.boundnz = -1 # reflective longitudinal BC
w3d.boundxy = -1                # reflective transverse BC

# --- sets particles boundary conditions
top.pboundxy = reflect#absorb
top.pbound0  = absorb
top.pboundnz = absorb

if dim=="1d":
  w3d.boundxy = periodic
  top.pboundxy = periodic

#-------------------------------------------------------------------------------
# Plasma channel
#-------------------------------------------------------------------------------
#  Matched spot size in channel (a0)
#  Note - to decrease channel dN at SIGMA, make matchspot > SIGMA
#  then dN/dN_matched = (SIGMA^4/matchspot^4), e.g. 5.7% increase in matchspot -> 20% detuning of dN 
#  tor: matchspot/sigma = (dN_matched/dN)^0.25
CHANNELFAC = 0.6 # factor by which density rise is less than matched density (to compensate self guide)
matchspot  = laser_radius*(1/CHANNELFAC)**0.25 # channel
max_parab_radius   = 0.8*w3d.xmmax
max_radius        = 0.96*w3d.xmmax
diff_density      = gammafrm*1.13e14/(matchspot**4)  #formula for the channel density without the radius squared term
if dim=="1d":diff_density=0.
norm_diff_density  = diff_density/dens0        #the normalized channel density without the radius squared term
max_diff_density  = norm_diff_density*(max_parab_radius**2) #the maximum normalized channel density

#-------------------------------------------------------------------------------
# set max time
#-------------------------------------------------------------------------------
tmaxlab = (Lplasma_lab)/clight
tmax = tmaxlab / gammafrm

#-------------------------------------------------------------------------------
# dump/plots intervals (in pico-seconds)
#-------------------------------------------------------------------------------
dump_intervals = tmax/20
beamdump_intervals = tmax/10
plot_intervals = tmax/10

#-------------------------------------------------------------------------------
# set graphics
#-------------------------------------------------------------------------------
if l_gist:
 if l_test:
  winon(0,dpi=dpi)
 else:
  setup()
  winon(0,dpi=dpi)
else:
 setup()
 
#-------------------------------------------------------------------------------
# set particles 
#-------------------------------------------------------------------------------
weight     = dens0*w3d.dx*w3d.dy*w3d.dz/(nppcellx*nppcelly*nppcellz) # weight of plasma macro-particles
weightbeam = 0. # needs to be fixed

# --- create e- beam species
if l_beam:
  beam = Species(type=Electron,weight=weightbeam)
# --- create plasma electron species
elec = Species(type=Electron,weight=weight)
# --- create plasma ion species
if l_ions: 
  ions = Species(type=Proton,charge_state=1,weight=weight)
#  ions.sm*=1.e20  # to give virtually infinite mass to ions
top.wpid = nextpid() # creates data space for variable weights (needed for plasma ramps)
top.depos_order[...] = top.depos_order[0,0] # sets deposition order of all species = those of species 0
top.efetch[...] = top.efetch[0] # same for field gathering
if dim in ["1d","2d"]:
  top.depos_order[1,:]=1
if dim=="1d":
  top.depos_order[0,:]=1

#-------------------------------------------------------------------------------
# set smoothing of current density
#-------------------------------------------------------------------------------
if l_smooth:
  # --- 1 time nilinear (0.25,0.5,0.25) + 1 time relocalization (-1, 3/2,-1.)
  npass_smooth = [[ 0 , 0 ],[ 0 , 0 ],[ 1 , 1 ]]
  alpha_smooth = [[ 0.5, 3.],[ 0.5, 3.],[0.5, 3./2]]
  stride_smooth = [[ 1 , 1 ],[ 1 , 1 ],[ 1 , 1 ]]
  if dim=='1d':
    for i in range(len(npass_smooth[0])):
      npass_smooth[0][i]=0
  if dim in ['1d','2d']:
    for i in range(len(npass_smooth[0])):
      npass_smooth[1][i]=0
else:
  npass_smooth = [[ 0 ],[ 0 ],[ 0 ]]
  alpha_smooth = [[ 1.],[ 1.],[ 1.]]
  stride_smooth = [[ 1 ],[ 1 ],[ 1 ]]

#-------------------------------------------------------------------------------
# initializes WARP
#-------------------------------------------------------------------------------
top.fstype = -1 # sets field solver to None (desactivates electrostatic solver)
package('w3d');generate()
#-------------------------------------------------------------------------------
# set a few shortcuts
#-------------------------------------------------------------------------------
pg = top.pgroup

#-------------------------------------------------------------------------------
# set input plasma density arrays
#-------------------------------------------------------------------------------
# --- set intial positions and weights
zpos = zstart0 = 0.
nppcell = nppcellx*nppcelly*nppcellz
dx = w3d.dx/nppcellx
dy = w3d.dy/nppcelly
dz = w3d.dz/nppcellz
nx = nppcellx*w3d.nx
ny = nppcelly*w3d.ny
nz = nppcellz

if dim in ["1d","2d"]:
  if dim=="1d":
    nx=1
    xp0,zp0 = getmesh2d(0.,dx,0,
                        dz/2,dz,nz-1)
    yp0 = xp0*0.
  else:
    xp0,zp0 = getmesh2d(w3d.xmmin+dx/2,dx,nx-1,
                        -w3d.dz+dz/2,dz,nz-1)
    yp0 = xp0*0.
else:
  xp0,yp0,zp0 = getmesh3d(w3d.xmmin+dx/2,dx,nx-1,
                       w3d.ymmin+dy/2,dy,ny-1,
                       -w3d.dz+dz/2,dz,nz-1)

zp0-=minnd(zp0) # ensures that zp0 starts at 0

# --- transform to 1D arrays
xp0=xp0.flatten()
yp0=yp0.flatten()
zp0=zp0.flatten()

# --- select particles within computational box of local processor
ii=compress((xp0>=w3d.xmminlocal) & (xp0<w3d.xmmaxlocal) & \
            (yp0>=w3d.ymminlocal) & (yp0<w3d.ymmaxlocal),arange(len(xp0)))
xp0=take(xp0,ii)
yp0=take(yp0,ii)
zp0=take(zp0,ii)

# --- select particles within max_radius
rp0=sqrt(xp0**2+yp0**2)
ii=compress(rp0<=max_radius,arange(len(rp0)))
xp0=take(xp0,ii)
yp0=take(yp0,ii)
zp0=take(zp0,ii)
rp0=sqrt(xp0**2+yp0**2)

# --- set the transverse profile
def plasma_trans_profile(r):
    wp = ones(shape(r)[0])
    slope = -(1+max_diff_density)/(max_radius-max_parab_radius)
    intercept = -max_radius*slope
    wp = where(r<max_parab_radius,1+norm_diff_density*(r**2),slope*r+intercept)
    wp = where(wp<0.,0.,wp)
    return wp

wp0=plasma_trans_profile(rp0)

def plasma_long_profile(wp0,z):
  wp = where(z<length_pramp,sin(0.5*pi*z/length_pramp)*wp0,wp0) # plasma density rises as half sin
  wp = where(z>Lplasma-length_pramp_exit,sin(0.5*pi*(Lplasma-z)/length_pramp_exit)*wp,wp) # plasma density falls as half sin
  return wp

def pldens():
  # --- plots longitudinal and transverse plasma profile
  nz = 1000
  za = arange(-0.,Lplasma,(Lplasma)/nz)
  wp = ones(shape(za)[0])
  for i,z in enumerate(za.tolist()):
    wp[i] = plasma_long_profile(1.,z)
  plsys(9)
  pla(wp,za*gammafrm*1000,width=3,color=red)
  limits(0,za[-1]*gammafrm*1000,0,1.1)
  ptitles('longitudinal density profile','z (mm)','')
  plsys(10)
  r = arange(0,w3d.xmmax,w3d.xmmax/nz)
  wp = plasma_trans_profile(r)
  pla(wp,r*1.e6,width=3,color=red)
  limits(0,max_radius*1.e6,0,1.1*max(wp))
  ptitles('radial density profile','r (microns)','')

# --- defines subroutine injecting plasma
def loadelec():
 global zstart0,zpos,xp0,yp0,zp0,wp0,l_moving_window
 while(zpos>=zstart0 and (zpos+betafrm*clight*top.time<Lplasma)): 
    z0 = zstart0+zp0
    # --- sets ramp by adjusting weight
    zi = z0+betafrm*clight*top.time
    wp = plasma_long_profile(wp0,zi)
    # --- sets velocity
    vx = 0.#001*clight*ranf(dz)
    vy = 0.#001*clight*ranf(dz)
    vz = -betafrm*clight#+0.001*clight*ranf(dz)
    # --- sets positions
    x = xp0.copy()
    y = yp0.copy()
    z = z0.copy()
    # --- inject electrons
    elec.addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=wp,lallindomain=false)#,lusespaceabove=true)
    # --- inject ions at same locations
    if l_ions:ions.addpart(x=x,y=y,z=z,vx=vx,vy=vy,vz=vz,w=wp,lallindomain=false)#,lusespaceabove=true)
    ladd=1
    zstart0+=w3d.dz
 if l_moving_window:
    zpos+=(top.vbeamfrm+em.vzgrid)*top.dt
#    zpos+=clight*top.dt # was not tested for a long time; is this correct?
 else:
    zpos+=clight*top.dt # was not tested for a long time; is this correct?
 zstart0-=betafrm*clight*top.dt

if l_plasma:installuserinjection(loadelec)

laser_total_duration=4.*laser_duration
laser_total_duration_lab=4.*laser_duration_lab
laser_total_length_lab=laser_total_duration_lab*clight
#-------------------------------------------------------------------------------
# set laser pulse shape
#-------------------------------------------------------------------------------
def laser_amplitude(time):
 global laser_total_length,Eamp,laser_duration
 l = laser_duration/2.354820
 return Eamp*exp(-0.5*((time-0.5*laser_total_duration)/l)**2)

def laser_profile(x,y):
  global laser_waist
  r2 = x**2 + y**2
  rw = laser_radius  
  return exp(-0.5*r2/rw**2)

#-------------------------------------------------------------------------------
# set laser amplitude by combining the pulse shape, laser profile, and laser phase
#-------------------------------------------------------------------------------
def laser_func(x,y,t):
  global laser_amplitude,laser_phase,laser_profile
  return laser_amplitude(t)*cos(k0*clight*t)*laser_profile(x,y)
    
def laser_func_circ(x,y,t):
#  --- test circular polarization
  global laser_amplitude,laser_phase,laser_profile
  E = laser_amplitude(t)*laser_profile(x,y)
  phase = k0*clight*t
  return [E*cos(phase),E*sin(phase)]

if l_laser:
  laser_source_z=0.
else:
  laser_func=laser_source_z=Eamp=None
  
#-------------------------------------------------------------------------------
# initializes main field solver block
#-------------------------------------------------------------------------------
em = EM3D(       laser_func=laser_func,
                 laser_source_z=laser_source_z,
                 laser_source_v=-betafrm*clight,
                 laser_polangle=laser_polangle,
                 laser_emax=Eamp,
                 stencil=stencil,
                 npass_smooth=npass_smooth,
                 alpha_smooth=alpha_smooth,
                 stride_smooth=stride_smooth,
                 l_2dxz=dim=="2d",
                 l_1dz=dim=="1d",
                 dtcoef=dtcoef,
                 l_correct_num_Cherenkov=l_correct_num_Cherenkov,
                 l_verbose=l_verbose)

#-------------------------------------------------------------------------------
# restarts from dump file
#-------------------------------------------------------------------------------
if l_restart:
  restore(dump_file)

# --- load diagnostics
exec(compile(open('lwfa_diags.py').read(), 'lwfa_diags.py', 'exec'))
if nzfieldstations>0:installafterstep(updateebstations)

#-------------------------------------------------------------------------------
# intializes e- beam
#-------------------------------------------------------------------------------
if l_beam:
  # --- add beam particles
  np_beam = 4000
  top.vbeam = E_BEAM_BETA*clight
  if me==0: # --- do only if processor 0
   if np_beam==1:
    # --- loads single test electron
    beam.addpart(x=array([0.]),
                 y=array([0.]),
                 z=array([0.]),
                 vx=array([0.]),
                 vy=array([0.]),
                 vz=array([E_BEAM_GAMMA*E_BEAM_BETA*clight]),
                 gi=array([1./E_BEAM_GAMMA]),
                         lmomentum=True,
                         lallindomain=True)
   else:
    # --- loads e- beam
     if l_usesavedist:
       # --- loads distribution from file
       try:
         ff=PR.PR(savedist+'.pdb')
       except:
         ff=PR.PR(savedist+'.pyp')
       ux = ff.xp[::svstride]*GAMMAVXSIGMA
       uy = ff.yp[::svstride]*GAMMAVYSIGMA
       uz = ff.dp[::svstride]*GAMMAVZSIGMA+E_BEAM_GAMMA*E_BEAM_BETA*clight
       gi = 1./sqrt(1.+(ux**2+uy**2+uz**2)/clight**2)
       beam.addpart(ff.x[::svstride]*E_BEAM_RADIUS*2,
                    ff.y[::svstride]*E_BEAM_RADIUS*2,
                    ff.z[::svstride]*E_BEAM_LENGTH,
                    ux,uy,uz,gi=gi,
                    lmomentum=True,
                    lallindomain=True)
     else:
       # --- loads gaussian electron beam
       beam.add_gaussian_dist(np=np_beam,
                              deltax=E_BEAM_RADIUS*2*1,
                              deltay=E_BEAM_RADIUS*2*1,
                              deltaz=E_BEAM_LENGTH,
                              vthx=GAMMAVXSIGMA*1,
                              vthy=GAMMAVYSIGMA*1,
                              vthz=GAMMAVZSIGMA,
                              zmean=0.,
                              vzmean=E_BEAM_GAMMA*E_BEAM_BETA*clight,
                              lmomentum=True,
                              lallindomain=True,
                              zdist='regular')
  np_beam = beam.getn()
  # --- sets e- beam macro-particles weight
  if dim=="1d":
    beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)*w3d.dx*w3d.dy*E_BEAM_LENGTH*(2.*pi)))/np_beam
  if dim=="2d":
    beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)*w3d.dy*E_BEAM_LENGTH*(2.*pi)))/np_beam
  if dim=="3d":
    beam.sw = (E_BEAM_DENSITY_PEAK*((E_BEAM_RADIUS*2)**2*E_BEAM_LENGTH*(2.*pi)*1.5))/np_beam
  # --- install e- beam diagnostic routine
  if sum(pg.nps)>0:
    # --- set beam position in lab frame
    pg.zp += zstart0-E_BEAM_DISTANCE_BEHIND_LASER-0.5*laser_total_length_lab
  # --- transform particle positions and velocity to boosted frame
  bf=Boosted_Frame(gammafrm,l_setselfb=0)
  bf.boost(beam,l_inject_plane=l_injectplane,lallindomain=1,l_rmzmean=0,zinject=-5.*w3d.dz)
  particleboundaries3d(top.pgroup,-1,False)

#-------------------------------------------------------------------------------
# register solver
#-------------------------------------------------------------------------------
print('register solver')
registersolver(em)
print('done')

#-------------------------------------------------------------------------------
# sets moving window velocity
#-------------------------------------------------------------------------------
if l_moving_window:
  top.vbeamfrm=0.
  em.vzgrid=BETAG_LINEAR*clight
  
def liveplots():
  if top.it%live_plot_freq==0:
    fma();


    if dim=="1d":
      pzxex(view=3,msize=2,titles=0,gridscale=1.e-12)
      ptitles('Ex [V/m]','z [um]')
      density=elec.get_density()
      plsys(4)
      pla(density)
      ptitles('n/n0','z [um]','X [um]')
      pzxez(view=5,msize=2,titles=0,gridscale=1.e-9)
      ptitles('Ez [V/m]','z [um]')
#      elec.ppzuz(view=6)
      plke(view=6)
      refresh()
    else:    
      pzxex(view=3,msize=2,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-12)
      ptitles('Ex [TV/m]','z [um]','X [um]')
#      pzxey(view=4,msize=2)
      density=elec.get_density()
      ppg(transpose(density),view=4,titles=0, \
          xmin=w3d.zmmin+top.zgrid,xmax=w3d.zmmax+top.zgrid, \
          ymin=w3d.xmmin,ymax=w3d.xmmax,
          xscale=1e6,yscale=1.e6,gridscale=1./dens0)

      ptitles('n/n0','z [um]','X [um]')
      pzxez(view=5,msize=2,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-9)
      ptitles('Ez [GV/m]','z [um]','X [um]')
#      elec.ppzuz(view=6)
      plke(view=6)
      refresh()

installafterstep(liveplots)

print('\nInitialization complete\n')

# if this is a test, then stop, else execute main loop
if l_test:
  print('<<< To execute n steps, type "step(n)" at the prompt >>>')
#  raise('')
else:
  # --- run until last beam particle has crossed the last beam z-station diagnostic
  doit = True
  while(doit):
    em.step(10)
    if getn()==0:
      doit = True
    else:
      doit = min(getz())<max(beamzstations/gammafrm-top.time*betafrm*clight)

savebeamstations()
if nzfieldstations>0:saveebstations()
