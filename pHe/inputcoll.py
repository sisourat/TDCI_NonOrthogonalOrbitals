 # input

debug = True
analyze = False
nstep_analysis = 10

ne = 2
tdoc_frozen = 0

tbasis = {'He': 'aug-ccpvdz'}
tgeom = "He 0 0 0.0 ; "
tcharge = 0
tspin = 0

pbasis = {'H': 'aug-ccpvdz'}
elp = "H"
xp = 0
yp = 0
zp = -1000
pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
pcharge = +1
pspin = 0

#i_init = 1
dtime = 0.05

zmax = 100.0
ngrid = 100
gridtype = 'exp'  #lin or exp
vproj = 2.0
bmin =  0.5
bmax =  6.5
nbb = 18
xmlfile = 'phe.xml'


