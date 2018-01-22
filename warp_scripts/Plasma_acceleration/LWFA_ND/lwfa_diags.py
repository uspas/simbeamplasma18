from warp import *
# --- define shortcut
ppg=ppgeneric

def labdata(z,t,ux,uy,uz,gi,uzfrm):
  if me==0 and l_verbose:print('enter labdata')
  np = shape(z)[0]
  gammafrm = sqrt(1.+uzfrm**2/clight**2)
  zpr = gammafrm*z-uzfrm*t
  tpr = gammafrm*t-uzfrm*z/clight**2
  setu_in_uzboosted_frame3d(np,ux,uy,uz,gi,uzfrm,gammafrm)
  if me==0 and l_verbose:print('exit labdata')
  return zpr,tpr,ux,uy,uz

def pzx(msize=1,color=red,titles=1,xscale=1.,yscale=1.,view=1):
  """
  plots ZX projection of e- beam
  """
  if not l_beam:return
  ppzx(msize=msize,color=color,titles=titles,xscale=xscale,yscale=yscale,view=view)
  try:
    ppzx(pgroup=bf.pgroup,color=blue,msize=msize,titles=titles,xscale=xscale,yscale=yscale,view=view)
  except:
    pass
  
def pzy(msize=1,color=red,titles=1,xscale=1.,yscale=1.,view=1):
  """
  plots ZY projection of e- beam
  """
  if not l_beam:return
  ppzy(msize=msize,color=color,titles=titles,xscale=xscale,yscale=yscale,view=view)
  try:
    ppzy(pgroup=bf.pgroup,color=blue,msize=msize,titles=titles,xscale=xscale,yscale=yscale,view=view)
  except:
    pass
  
def pxy(msize=1,color=red,titles=1,xscale=1.,yscale=1.,view=1):
  """
  plots XY projection of e- beam
  """
  if not l_beam:return
  ppxy(msize=msize,color=color,titles=titles,xscale=xscale,yscale=yscale,view=view)
  try:
    ppxy(pgroup=bf.pgroup,color=blue,msize=msize,titles=titles,xscale=xscale,yscale=yscale,view=view)
  except:
    pass
  
def pzxex(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  if dim=='1d':
    em.pfex(l_transpose=dim!='1d',direction=1,view=view,titles=titles)
    z=getz()
    if me==0:
      ppg(z*0,z*xscale,color=red,msize=msize,view=view,titles=titles)
    try:
      z=getz(pgroup=bf.pgroup)
      if me==0:
        ppg(z*0,z*xscale,color=blue,msize=msize,view=view,titles=titles)
    except:
      pass
  else:  
    em.pfex(l_transpose=dim!='1d',direction=1,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
    pzx(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)

def pzxey(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  if dim=='1d':
    em.pfey(l_transpose=dim!='1d',direction=1,view=view,titles=titles)
    z=getz()
    if me==0:
      ppg(z*0,z*xscale,color=red,msize=msize,view=view,titles=titles,xscale=xscale,yscale=yscale)
    try:
      z=getz(pgroup=bf.pgroup)
      if me==0:
        ppg(z*0,z*xscale,color=blue,msize=msize,view=view,titles=titles,xscale=xscale,yscale=yscale)
    except:
      pass
  else:  
    em.pfey(l_transpose=dim!='1d',direction=1,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
    pzx(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)

def pzxez(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  if dim=='1d':
    em.pfez(l_transpose=dim!='1d',direction=1,view=view,titles=titles)
    z=getz()
    if me==0:
      ppg(z*0,z*xscale,color=red,msize=msize,view=view,titles=titles,xscale=xscale,yscale=yscale)
    try:
      z=getz(pgroup=bf.pgroup)
      if me==0:
        ppg(z*0,z*xscale,color=blue,msize=msize,view=view,titles=titles,xscale=xscale,yscale=yscale)
    except:
      pass
  else:  
    em.pfez(l_transpose=dim!='1d',direction=1,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
    pzx(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)
  
def pzyex(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  em.pfex(l_transpose=1,direction=0,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
  pzy(msize=msize,view=view,titles=0,xscale=xscale,yscale=yscale)

def pzyey(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  em.pfey(l_transpose=1,direction=0,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
  pzy(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)

def pzyez(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  em.pfez(l_transpose=1,direction=0,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
  pzy(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)
  
def pxyex(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  em.pfex(l_transpose=1,direction=2,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
  pxy(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)

def pxyey(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  em.pfey(l_transpose=1,direction=2,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
  pxy(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)

def pxyez(msize=1,titles=1,xscale=1.,yscale=1.,view=1,gridscale=1.):
  em.pfez(l_transpose=1,direction=2,view=view,titles=titles,xscale=xscale,yscale=yscale,gridscale=gridscale)
  pxy(msize=msize,titles=0,view=view,xscale=xscale,yscale=yscale)
  
zstart0lab=0.
dzstations=Lplasma_lab/nzstations
beamzstations = zstart0lab+arange(0.,Lplasma_lab,dzstations) # list of diag stations z-locations in lab frame

ekstations = zeros(shape(beamzstations),'d')
ppzstations = zeros(shape(beamzstations),'d')
xbarstations = zeros(shape(beamzstations),'d')
xpbarstations = zeros(shape(beamzstations),'d')
xsqstations = zeros(shape(beamzstations),'d')
xpsqstations = zeros(shape(beamzstations),'d')
xxpstations = zeros(shape(beamzstations),'d')
if dim == "3d":
  ybarstations = zeros(shape(beamzstations),'d')
  ypbarstations = zeros(shape(beamzstations),'d')
  ysqstations = zeros(shape(beamzstations),'d')
  ypsqstations = zeros(shape(beamzstations),'d')
  yypstations = zeros(shape(beamzstations),'d')
tbarstations = zeros(shape(beamzstations),'d')
tsqstations = zeros(shape(beamzstations),'d')
ekstationstime = zeros(shape(beamzstations),'d')
ekstationscnt = zeros(shape(beamzstations),'d')
ekstationscnt2 = zeros(shape(beamzstations),'d')
npz = 1001
pzbeamstations = (200.e6/dfact*arange(npz))/(npz-1)
pzstations = zeros((shape(beamzstations)[0],npz),'d')

top.zoldpid=nextpid()
def updatebeamstations():
  global timestart
#  if top.it%10<>0:return
  if me==0 and l_verbose:print('enter updatebeamstations')
  # --- compute beta*gamma*c
  uzfrm=-betafrm*gammafrm*clight
  # --- get nb particles on each CPU
  np = getn(gather=0,bcast=0) 
  if np>0:
    # --- get z on each CPU
    z=getz(gather=0,bcast=0).copy()
    zold=getpid(id=top.zoldpid-1,gather=0,bcast=0)
    zoldlab = gammafrm*zold-uzfrm*(top.time-top.dt)
    # --- get z, time and velocities in lab frame
    zlab,tlab,uxlab,uylab,uzlab = labdata(z,
                                          top.time,
                                          getux(gather=0,bcast=0).copy(),
                                          getuy(gather=0,bcast=0).copy(),
                                          getuz(gather=0,bcast=0).copy(),
                                          getgaminv(gather=0,bcast=0).copy(),
                                          uzfrm=uzfrm)
    w = abs(zlab-zoldlab)/dzstations
    # --- get x,y on each CPU
    x = getx(gather=0,bcast=0).copy()
    y = gety(gather=0,bcast=0).copy()
    # --- compute gamma in lab frame
    myglab = sqrt(1.+(uxlab**2+uylab**2+uzlab**2)/clight**2)
    # --- compute kinetic energy in lab frame
    mykelab = beam.sm*(myglab-1.)*clight**2/echarge
    # --- defines cutoffs if particle selection is ON
    if l_pselect:
      # --- set threshold on radius
      XYcutoff = E_BEAM_RADIUS*5.
      # --- set threshold on longitudinal velocity
      UZcutoff = 0.95*E_BEAM_GAMMA*E_BEAM_BETA*clight
      # --- set threshold on energy
      KEcutoff = 1.e6 # eV
#      KEcutoff = None
#      XYcutoff = None
#      UZcutoff = None
    else:
      XYcutoff = None
      UZcutoff = None
      KEcutoff  = None
    if XYcutoff is not None:
      # --- select particle based on radius
      if dim=="3d":
        r2 = x*x+y*y
        XYcutoff2 = XYcutoff**2
        ii = compress(r2<XYcutoff2,arange(np))
      else:
        ii = compress(abs(x)<XYcutoff,arange(np))
      # --- get # of selected particles
      np = len(ii)
      # --- get weight, position, time, velocity and energy of selected particles
      w = take(w,ii)
      x = take(x,ii)
      y = take(y,ii)
      zlab = take(zlab,ii)
      tlab = take(tlab,ii)
      uxlab = take(uxlab,ii)
      uylab = take(uylab,ii)
      uzlab = take(uzlab,ii)
      mykelab = take(mykelab,ii)
    if UZcutoff is not None:
      # --- select particle based on longitudinal velocity
      ii = compress(uzlab>UZcutoff,arange(np))
      # --- get # of selected particles
      np = len(ii)
      # --- get weight, position, time, velocity and energy of selected particles
      w = take(w,ii)
      x = take(x,ii)
      y = take(y,ii)
      zlab = take(zlab,ii)
      tlab = take(tlab,ii)
      uxlab = take(uxlab,ii)
      uylab = take(uylab,ii)
      uzlab = take(uzlab,ii)
      mykelab = take(mykelab,ii)
    if KEcutoff is not None:
      # --- select particle based on energy
      ii = compress(mykelab>KEcutoff,arange(np))
      # --- get # of selected particles
      np = len(ii)
      # --- get weight, position, time, velocity and energy of selected particles
      w = take(w,ii)
      x = take(x,ii)
      y = take(y,ii)
      zlab = take(zlab,ii)
      tlab = take(tlab,ii)
      uxlab = take(uxlab,ii)
      uylab = take(uylab,ii)
      uzlab = take(uzlab,ii)
      mykelab = take(mykelab,ii)
    if np>0:
      xplab = uxlab/clight # normalized (gamma*beta*xp)
      yplab = uylab/clight # normalized (gamma*beta*yp) 
      nz = shape(ekstations)[0]
      deposgrid1dw(1,np,zlab,mykelab,w,nz-1,ekstations,ekstationscnt,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,beam.sm*uzlab*clight/echarge,w,nz-1,ppzstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,x,w,nz-1,xbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,x**2,w,nz-1,xsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,xplab,w,nz-1,xpbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,xplab**2,w,nz-1,xpsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,x*xplab,w,nz-1,xxpstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      if dim == "3d":
        deposgrid1dw(1,np,zlab,y,w,nz-1,ybarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,y**2,w,nz-1,ysqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,yplab,w,nz-1,ypbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,yplab**2,w,nz-1,ypsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
        deposgrid1dw(1,np,zlab,y*yplab,w,nz-1,yypstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,tlab,w,nz-1,tbarstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      deposgrid1dw(1,np,zlab,tlab**2,w,nz-1,tsqstations,ekstationscnt2,beamzstations[0],beamzstations[-1])
      setgrid2dw(np,zlab,uzlab*beam.sm*clight/echarge,w,nz-1,npz-1,pzstations,
                  beamzstations[0],beamzstations[-1],pzbeamstations[0],pzbeamstations[-1])
  if top.it%hist_freq==0:
    savebeamstations()
  if me==0 and l_verbose:print('exit updatebeamstations')

if l_beam:installafterstep(updatebeamstations)

def savebeamstations():
 if me==0 and l_verbose:print('enter savebeamstations')
 pnums = parallelsum(ekstationscnt)
 pnum = where(pnums==0.,1.,pnums)
 ekst = parallelsum(ekstations)/pnum
 ppzst = parallelsum(ppzstations)/pnum
 xbar = parallelsum(xbarstations)/pnum
 xsq = parallelsum(xsqstations)/pnum
 xpbar = parallelsum(xpbarstations)/pnum
 xpsq = parallelsum(xpsqstations)/pnum
 xxp = parallelsum(xxpstations)/pnum
 if dim == "3d":
   ybar = parallelsum(ybarstations)/pnum
   ysq = parallelsum(ysqstations)/pnum
   ypbar = parallelsum(ypbarstations)/pnum
   ypsq = parallelsum(ypsqstations)/pnum
   yyp = parallelsum(yypstations)/pnum
 tbar = parallelsum(tbarstations)/pnum
 tsq = parallelsum(tsqstations)/pnum
 wti  = parallelsum(ekstationstime)/pnum
 pzst = parallelsum(pzstations)#*beam.sw
 if me==0:
  os.system('mv -f ebeamstations.pdb ebeamstationsold.pdb')
  f = PW.PW('ebeamstations.pdb')
  f.ekstations = ekst
  f.ppzstations = ppzst
  f.xbarstations = xbar
  f.xrmsstations = sqrt(xsq)
  f.xpbarstations = xpbar
  f.xprmsstations = sqrt(xpsq)
  f.xemitnstations = sqrt((xsq-xbar*xbar)*(xpsq-xpbar*xpbar)-(xxp-xbar*xpbar)**2)
  if dim == "3d":
    f.ybarstations = ybar
    f.yrmsstations = sqrt(ysq)
    f.ypbarstations = ypbar
    f.yprmsstations = sqrt(ypsq)
    f.yemitnstations = sqrt((ysq-ybar*ybar)*(ypsq-ypbar*ypbar)-(yyp-ybar*ypbar)**2)
  f.tbarstations = tbar
  f.trmsstations = sqrt(tsq-tbar*tbar)
  f.ekstationstime = wti
  f.pzstations = pzst
  f.beamzstations = beamzstations-zstart0lab
  f.pzbeamstations = pzbeamstations
  f.pnumstations = pnums
  f.nx = w3d.nx
  f.ny = w3d.ny
  f.nz = w3d.nz
  f.time = top.time
  f.dt = top.dt
  f.it = top.it
  f.stencil=stencil
  f.dim=dim
  f.close()
  os.system('rm -f ebeamstationsold.pdb')
 if me==0 and l_verbose:print('exit savebeamstations')
  
def plke(view=1):
 global kelab,pxlab,pylab,pzlab,zhlab
 if me==0 and l_verbose:print('enter plke')
 ekcnt = parallelsum(ekstationscnt)
 ekcnt = where(ekcnt==0.,1.,ekcnt)
 ekst = parallelsum(ekstations)/ekcnt
 if me==0:
    plsys(view)
    pla(ekst*1.e-6,beamzstations*1.e3,color=red)
    ptitles('Energy (MeV)','Z (mm)','')
 if me==0 and l_verbose:print('exit plke')

if nzfieldstations>0:
 zstations = arange(0.,Lplasma_lab,Lplasma_lab/nzfieldstations) # list of field diag stations z-locations in lab frame
 exstations = []
 eystations = []
 ezstations = []
 bxstations = []
 bystations = []
 bzstations = []
 tstations = []
 for i in range(shape(zstations)[0]):
    exstations.append(AppendableArray(typecode='d'))
    eystations.append(AppendableArray(typecode='d'))
    ezstations.append(AppendableArray(typecode='d'))
    bxstations.append(AppendableArray(typecode='d'))
    bystations.append(AppendableArray(typecode='d'))
    bzstations.append(AppendableArray(typecode='d'))
    tstations.append(AppendableArray(typecode='d'))

 def updateebstations():
  #  --- routine for accumulating EM field value at z locations (array zstations)
  global em,hist_freq
  if me==0 and l_verbose:print('enter updateebstations')
  # --- compute z in calculation (boosted) frame
  zcalc = zstations/gammafrm-betafrm*clight*top.time

  # --- select z that are within grid range
  n = shape(zcalc)[0]
  ilist = compress((zcalc>=(w3d.zmmin+top.zgrid)) & (zcalc<(w3d.zmmax+top.zgrid)),arange(n))
  zcalc = take(zcalc,ilist)
  n = shape(zcalc)[0]
  if n==0:return

  # --- gather EM fields in calculation frame
  x=zcalc*0.
  y=zcalc*0.
  ex,ey,ez,bx,by,bz = em.getfieldsfrompositions(x,y,zcalc)

  if me==0:
    uxf = 0.
    uyf = 0.
    uzf = -gammafrm*betafrm*clight
    gammaf = gammafrm
    # --- convert EM fields to lab frame
    seteb_in_boosted_frame(n,ex,ey,ez,bx,by,bz,uxf,uyf,uzf,gammaf)
    # --- compute time in lab frame
    t = gammafrm*(top.time+betafrm*zcalc/clight)
    # --- store field and time values in appendable arrays
    for j,i in enumerate(ilist):
        tstations[i].append(t[j])
        exstations[i].append(ex[j])
        eystations[i].append(ey[j])
        ezstations[i].append(ez[j])
        bxstations[i].append(bx[j])
        bystations[i].append(by[j])
        bzstations[i].append(bz[j])

    # --- save data every hist_freq time steps
    if top.it%hist_freq==0:
      saveebstations()
  if me==0 and l_verbose:print('exit updateebstations')

 def  saveebstations():
  if me>0:return
  if me==0 and l_verbose:print('enter saveebstations')
  fname='ebhist.pdb'
  os.system('mv -f ebhist.pdb ebhistold.pdb')
  f=PW.PW(fname)
  ex = []
  ey = []
  ez = []
  bx = []
  by = []
  bz = []
  t = []
  for i in range(len(exstations)):
    ex.append(exstations[i][:])
    ey.append(eystations[i][:])
    ez.append(ezstations[i][:])
    bx.append(bxstations[i][:])
    by.append(bystations[i][:])
    bz.append(bzstations[i][:])
    t.append(tstations[i][:])
  f.ex=ex
  f.ey=ey
  f.ez=ez
  f.bx=bx
  f.by=by
  f.bz=bz
  f.t=t
  f.z=zstations

  f.close()
  os.system('rm -f ebhistold.pdb')
  if me==0 and l_verbose:print('exit saveebstations')


tottime = AppendableArray()
def accuttime():
  global tottime
  tottime.append(time.clock())
  if me==0 and top.it%200==0:
    f=PW.PW('tottime.pdb')
    f.time=tottime[:]
    f.close()

installafterstep(accuttime)


