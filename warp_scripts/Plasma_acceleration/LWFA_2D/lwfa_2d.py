"""
This is a Warp input script for laser-wakefield acceleration.

Usage
-----
- Type "python -i lwfa_2d.py" in a terminal
"""
from warp.init_tools import *  # Import warp-specific packages

# --------------------------------------------------------
# Parameters (Modify the values below to suit your needs)
# --------------------------------------------------------

# General parameters
# ------------------
dim = "2d"  # Dimension of simulation ("3d" or "2d")
N_steps = 6000  # Total number of timesteps in the simulation
diag_period = 200  # Period of diagnostics (in number of timesteps)
cgmplotfreq = 200

# Simulation box
# --------------
# Number of grid cells in each direction
Nx = 40*2 ; Ny = 40 ; Nz = 110*2
# Dimension of the box in longitudinal direction (meters)
zmin = -11.e-6 ; zmax = 0.e-6
# Dimension of the box in transverse direction (meters)
xmax = ymax = 10e-6
# Field boundary conditions (longitudinal and transverse respectively)
f_boundz = f_boundxy = openbc
# Particles boundary conditions (longitudinal and transverse respectively)
p_boundz = p_boundxy = absorb

# Numerical parameters
# --------------------
stencil = 1  # Field solver (0:Yee, 1:Karkkainen on EF,B)
depos_order = 1  # Particle shape (1:linear, 2:quadratic, 3:cubic)

# Gaussian laser parameters
# -------------------------
laser_a0 = 2.5  # Laser amplitude at focus
laser_w0 = 4.e-6  # Waist at focus (meters)
laser_ctau = 1.5e-6  # Length of the pulse (meters)
laser_z0 = zmax - 2*laser_ctau  # Initial position of the centroid (meters)

# Plasma content and profile
# --------------------------
n_plasma = 2.5e25  # Plasma density (in number of particles per m^3)
# Number of macroparticles per cell in each direction
plasma_nx = 1 ; plasma_ny = 1 ; plasma_nz = 1
# Positions between which the plasma is initialized
plasma_zmin = 0.e-6 ; plasma_zmax = 1500.e-6
# Define your own profile and profile parameters below
ramp_start = 0.e-6 ; ramp_length = 20.e-6
def plasma_dens_func( x, y, z ):
    """User-defined function: density profile at position x, y, z"""
    # Allocate relative density
    n = ones_like(z)
    # Make linear ramp
    n = where( z<ramp_start+ramp_length, (z-ramp_start)/ramp_length, n )
    return(n)

# Relativistic beam
# -----------------
beam_tot_charge = 10.e-12  # Total charge (C)
beam_gamma = 1000.  # Lorentz factor
beam_Np = 1000  # Number of macroparticles
beam_zmean = -8.5e-6 # Mean position of the beam in z (meters)
# RMS width of the beam in each direction
beam_zrms = 0.1e-6 ; beam_xrms = 0.1e-6 ; beam_yrms = 0.1e-6

# -----------------------------------------------------------------------------
# Initialization of the simulation (Normal users should not modify this part.)
# -----------------------------------------------------------------------------

# Set some general options for warp
# ---------------------------------
set_diagnostics( interactive=False )
set_boundary_conditions( f_boundz, f_boundxy, p_boundz, p_boundxy )
set_simulation_box( Nz, Nx, Ny, zmin, zmax, xmax, ymax, dim )
set_moving_window( True, v_moving_window=clight )

# Creation of the species
# -----------------------
# Create the plasma species
elec_weight = prepare_weights( n_plasma, plasma_nx, plasma_ny, plasma_nz, dim, circ_m=0 )
elec = Species(type=Electron, weight=elec_weight, name='electrons')
# Create the beam
beam_weight = beam_tot_charge / beam_Np / echarge
beam = Species(type=Electron, weight=beam_weight, name='beam')
# Set the numerical parameters only now: they affect the newly created species
top.nhist = 5
set_numerics( depos_order, efetch=1, particle_pusher=1, dim=dim)

# Setup the field solver object
# -----------------------------
em = EM3D( stencil=stencil, 
           npass_smooth=[[1],[0],[1]],
           alpha_smooth=[[0.5],[0.],[0.5]],
           l_2dxz=(dim == "2d"), l_getrho=True )
registersolver(em)

# Load the laser and the electron beam (with initial space-charge fields)
# -----------------------------------------------------------------------
add_laser( em, dim, laser_a0, laser_w0,
           laser_ctau, laser_z0, theta_pol=pi/2, source_z=zmax-1.e-6 )
# Load the beam
beam.add_gaussian_dist( beam_Np, beam_xrms, beam_yrms, beam_zrms,
            zmean=beam_zmean, vzmean=clight*sqrt(1-1./beam_gamma**2) )
initialize_beam_fields( em, dim, beam, w3d, top )

# Introduce the plasma
# --------------------
# Create an object to store the information about plasma injection
plasma_injector = PlasmaInjector( elec, None, w3d, top, dim,
        plasma_nx, plasma_ny, plasma_nz, plasma_zmin,
        plasma_zmax, xmax, ymax, plasma_dens_func )
# Continuously inject the plasma, as the moving window progresses
installuserinjection( plasma_injector.continuous_injection )
        
# Setup the diagnostics
# ---------------------
remove_existing_directory( ['diags'] )
diag1 = FieldDiagnostic( period=diag_period, top=top, w3d=w3d, em=em,
    comm_world=comm_world, write_dir='diags' )
diag1.write()
installafterstep( diag1.write )
diag2 = ParticleDiagnostic( period=diag_period, top=top, w3d=w3d,
    species={ species.name : species for species in listofallspecies }, 
    comm_world=comm_world, write_dir='diags' )
diag2.write()
installafterstep( diag2.write )

# --- accumulation of history of beam kinetic energy and position
zh = AppendableArray()
ekh = AppendableArray()
@callfromafterstep
def accuhist():
    if top.it%10==0:
        zh.append(ave(beam.getz()))
        ekh.append(ave(beam.getke()))

# Setup up cgm plots
# ------------------
if cgmplotfreq:
    setup()
    winon()
    @callfromafterstep
    def cgmplots():
        if top.it%cgmplotfreq==0:
            fma()
            if dim=='3d':
                em.pfez(direction=1,view=9,l_transpose=1,contours=100,filled=1,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-9)
                beam.ppzx(view=9,color=red,titles=0,xscale=1e6,yscale=1.e6)
                ptitles('E_z^^ [GV/m]','Z [um]','X [um]')
                em.pfez(direction=0,view=10,l_transpose=1,contours=100,filled=1,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-9)
                beam.ppzy(view=10,color=red,titles=0,xscale=1e6,yscale=1.e6)
                ptitles('E_z^^ [GV/m]','Z [um]','Y [um]')
            if dim=='2d':
                d=elec.get_density()/n_plasma 
                ppgeneric(transpose(d),xmin=(w3d.zmmin+top.zgrid)*1.e6,
                            xmax=(w3d.zmmax+top.zgrid)*1.e6,
                            ymin=w3d.xmmin*1.e6,
                            ymax=w3d.xmmax*1.e6,
                            contours=100,filled=1,view=14)
                ptitles('n/n0','Z [um]','X [um]')
                em.pfez(l_transpose=1,contours=100,filled=1,titles=0,xscale=1e6,yscale=1.e6,gridscale=1.e-9,view=15)
                beam.ppzx(color=red,titles=0,xscale=1e6,yscale=1.e6,view=15)
                ptitles('E_z^^ [GV/m]','Z [um]','X [um]')
                ez=em.gatherez()
                plsys(16)
                pla(ez[int(w3d.nx/2),:]*1.e-9,1.e6*(w3d.zmesh+top.zgrid))
                ptitles('','Z [um]','E_z^^ on axis [GV/m]')
            refresh()

# -----------------------------------------------------------------------------
# Accumulation of CPU time and saving in output file
# -----------------------------------------------------------------------------
tottime = AppendableArray()
def savetottime():
    global tottime
    f=PW.PW('tottime.pdb')
    f.time=tottime[:]
    f.close()
@callfromafterstep
def accuttime():
    global tottime
    tottime.append(top.steptime)
    if me==0 and top.it%200==0:
        savetottime()

# -----------------------------------------------------------------------------
# Simulation loop (Normal users should not modify this part either.)
# -----------------------------------------------------------------------------
step(N_steps)
savetottime()

# --- plots e- beam kinetic energy history
window(3)
pla(ekh[...]*1.e-6,zh[...]*1.e6,color=red,width=3)
ptitles('','Z [um]','KE [MeV]')
pdf('KE_history_2d')
