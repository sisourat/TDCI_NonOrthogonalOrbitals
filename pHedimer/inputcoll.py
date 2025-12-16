 # input

debug = True
analyze = True
nstep_analysis = 10

ne = 4

tbasis = {'He': 'ccpvdz', 'ghost': [[0, [38.36, 0.023809], [5.77, 0.154891], [1.24, 0.469987]], [0, [0.2976, 1.0]], [1, [1.275, 1.0]]]
}
tecp = {  }
tgeom = "He 0 0 -12.0 ; He  0 0.0 12.0; "
tcharge = 0
tspin = 0

pbasis = {'H': 'aug-ccpvdz'}
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

zmax = 60.0
ngrid = 150
gridtype = 'lin'  #lin or exp
vproj = 0.6
bmin =  3.5
bmax =  6.5
nbb = 1
xmlfile = 'phedimer.xml'


