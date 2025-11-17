 # input

debug = True
ne = 4

tbasis = {'He': 'ccpvdz', 'ghost': [[0, [38.36, 0.023809], [5.77, 0.154891], [1.24, 0.469987]], [0, [0.2976, 1.0]], [1, [1.275, 1.0]]]
}
tgeom = "He 0 0 -4; He 0 0 4; "
tcharge = 0
tspin = 0

pbasis = {'H': 'augccpvdz'}
elp = "H"
xp = 0
yp = 0
zp = -1000
pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
pcharge = +1
pspin = 0

i_init = 1

zmax = 100.0
ngrid = 300
gridtype = 'lin'  #lin or exp
vproj = 0.4
bmin =  0.5
bmax =  6.5
nbb = 24
xmlfile = 'phe.xml'

