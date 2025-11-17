 # input

debug = True
ne = 2

tbasis = {'He': 'sto-3g'}
tgeom = "He 0 0 0;"
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

i_init = 1

zmax = -100.0
ngrid = 100
vproj = 0.4
bmin =  0.5
bmax =  8.5
nbb = 16
xmlfile = 'phe.xml'

