from pyscf import gto
import numpy as np
import sys
import os
import pathlib

from libcollision import *
from libdyn import *
#from inputcoll import *
from generate_csfs import *
from cimat import *
from scipy.interpolate import interp1d

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

  if(sys.argv[1]=='0'):
      sys.exit()

# Collision dynamics
  blist = np.linspace(bmin,bmax,nbb)
  if(gridtype=='exp'):
    zlist =  np.sort(np.concatenate((-np.logspace(-18, 0, base=2, num=ngrid)*zmax,np.logspace(-18, 0, base=2, num=ngrid)*zmax)))
  elif(gridtype=='lin'):
    zlist =  np.linspace(-zmax,zmax,ngrid)
  else:
    raise NotImplementedError('Exp or Lin Grid only')

  nmo = ntmo + npmo

  h1emo = np.zeros((nmo,nmo))*1.0j
  ovmo = np.zeros((nmo,nmo))*1.0j

  smo = np.zeros((nmo,nmo))
  smo[0:ntmo,0:ntmo] = tmo
  smo[ntmo:nmo,ntmo:nmo] = pmo

  # Checks the asymptotic energies
  xp = 0
  yp = 0
  zp = -10000.0
  pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
  sgeom = tgeom + pgeom
  mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
  mol.build( unit = 'Bohr')

  csfs = process_xml_csf(xmlfile)
  ncsfs = len(csfs)

  # Asymptotic Energies
  phase = np.ones(ncsfs)
  ovl, kin, pot = hcore(mol,smo)
  h1e = kin + pot
  eri = twoeints(mol,smo)

  hmat, smat = cimat(ne,nmo,ovl,h1e,eri,csfs,phase)

  esta = []
  net = []
  nep = []
  print()
  print("Asymptotic energies")
  print()
  for i, csf in enumerate(csfs):
    _, alpe, betae =  csf.terms[0]
    nte = int(np.count_nonzero(np.array(alpe)<ntmo) + np.count_nonzero(np.array(betae)<ntmo))
    npe = len(alpe)+len(betae)-nte
    net.append(nte)
    nep.append(npe)
    esta.append(hmat[i,i].real)
    print(i, esta[i], "          ", alpe, betae, nte, npe)
    #print(i, esta[i], "          ", csf)

  for b in blist:
   tmat = []
   ovl_mo = []
   for zproj in zlist:
     time = zproj/vproj
     phase = []
     for i, csf in enumerate(csfs):
       if(nep[i]==0):
         phase.append(1.0)
       elif(nep[i]==1):
         phase.append(np.exp(-vproj*zproj*1.0j)*np.exp(+0.5*vproj**2*time*1.0j))
       elif(nep[i]==2):
         phase.append(np.exp(-vproj*zproj*1.0j)**2*np.exp(+vproj**2*time*1.0j))

     xp = b
     yp = 0
     zp = zproj
     pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
     sgeom = tgeom + pgeom
     mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
     mol.build( unit = 'Bohr')

     ovl, kin, pot = hcore(mol,smo)
     ovl_mo.append(ovl)
     h1e = kin + pot
     eri = twoeints(mol,smo)

     hmat, smat = cimat(ne,nmo,ovl,h1e,eri,csfs,phase)
     inv = np.linalg.inv(smat)
     mat = np.matmul(inv, hmat)
     tmat.append(mat)

   # running dynamics for bproj = b
   tlist = zlist/vproj
   hmat = interp1d(tlist, tmat , axis=0)
   psi0 = np.zeros(ncsfs)
   psi0[i_init-1] = 1.0
   ntime = int(2.0*zmax/(vproj*dtime))
   t_grid = np.linspace(zlist[0]/vproj,zlist[-1]/vproj,ntime)
   wf_t = solve_tdse_numba(hmat, psi0, t_grid, hbar=1.0, method='auto', normalize=True)
   prob = np.abs(wf_t[-1])**2
   print(b,' '.join(map(str,prob)),np.sum(prob))

   #analyzing the dynamics for bproj = b
   if analyze == True:
    fout = open('prob_time_'+str(vproj)+'_'+str(b)+'.txt','w')
    for it, time in enumerate(t_grid):
       zproj = vproj * time
       if(it%10==0):
        amp2 = np.abs(wf_t[it])**2
        print(zproj,' '.join(map(str,amp2)),np.sum(amp2),file=fout)
    fout.close()
    # can do 1rdm analysis here

