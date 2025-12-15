 # input

debug = False
analyze = False
nstep_analysis = 10

ne = 8
tdoc_frozen = 5

tbasis = {'Ar': 'ccpvdz' }
tecp = {  }
tgeom = "Ar 0 0 0.0; "
tcharge = 0
tspin = 0

pbasis = {'H': 'sto-3g'}
pecp = {  }
elp = "H"
xp = 0
yp = 0
zp = -1000
pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
pcharge = +1
pspin = 0

i_init = 1
dtime = 0.05

zmax = 200.0
ngrid = 60
gridtype = 'exp'  #lin or exp
vproj = 1.0
bmin =  0.5
bmax =  8.5
nbb = 1
xmlfile = 'csfs.xml'


