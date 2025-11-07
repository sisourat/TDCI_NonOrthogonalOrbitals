from pyscf import gto
import numpy as np
import sys
import os
import pathlib

from libcollision import *
#from inputcoll import *
from numpy import linalg as LA
from generate_csfs import *
from cimat import *

if __name__ == "__main__":

  pdir = pathlib.Path().resolve()
  sys.path.append(pdir)
  from inputcoll import *
  sgeom = tgeom + pgeom
  sbasis = tbasis | pbasis
  scharge = tcharge + pcharge
  sspin = tspin + pspin

  tmol = gto.M(atom=tgeom,basis=tbasis,charge=tcharge,spin=tspin)
  tmol.build( unit = 'Bohr')
  pmol = gto.M(atom=pgeom,basis=pbasis,charge=pcharge,spin=pspin)
  pmol.build( unit = 'Bohr')
  mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
  mol.build( unit = 'Bohr')

# Computes target and projectile HF orbitals and ""model"" orbitals
  print("TARGET")
  tmo, tmo_e = system(tmol,debug)
  ntmo = len(tmo)

  print()

  print("PROJECTILE")
  pmo, pmo_e = system(pmol,debug)
  npmo = len(pmo)

  csfs = process_xml_csf(xmlfile)
  ncsfs = len(csfs)

  if(sys.argv[1]=='0'):
      sys.exit()

# Collision dynamics
  rlist = np.linspace(0.5,zmax,ngrid)

  nmo = ntmo + npmo

  h1emo = np.zeros((nmo,nmo))*1.0j
  ovmo = np.zeros((nmo,nmo))*1.0j

  smo = np.zeros((nmo,nmo))
  smo[0:ntmo,0:ntmo] = tmo
  smo[ntmo:nmo,ntmo:nmo] = pmo

  for r in rlist:
     phase = np.ones(ncsfs)
     xp = 0
     yp = 0
     zp = r
     pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
     sgeom = tgeom + pgeom
     mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
     mol.build( unit = 'Bohr')
     nucrep = mol.energy_nuc()

     ovl, kin, pot = hcore(mol,smo)
     h1e = kin + pot
     eri = twoeints(mol,smo)

     hmat, smat = cimat(ne,nmo,ovl,h1e,eri,csfs,phase)
     # Print S^-1M in a file and call Prop
     inv = np.linalg.inv(smat)
     mat = np.matmul(inv, hmat)

     eigenvalues, eigenvectors = LA.eig(mat)
     sorted_indices = np.argsort(eigenvalues)
     e = eigenvalues[sorted_indices].real+nucrep
     print(r,' '.join(map(str,e)))
     #print(hmat)

     #print(r,end=" ")
     #for i in range(len(ovl)-1):
     #    print(' '.join(map(str,np.abs(ovl[i])**2)),end=" ")
     #print(' '.join(map(str,np.abs(ovl[-1])**2)))

     #print(r,end=" ")
     #for i in range(len(smat)-1):
     #    print(' '.join(map(str,np.abs(smat[i])**2)),end=" ")
     #print(' '.join(map(str,np.abs(smat[-1])**2)))






