from pyscf import gto
import numpy as np
import sys
import os
import pathlib
import matplotlib.pyplot as plt

from libcollision import *
from libdyn import *
from libanalysis import *
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
  zp = -zmax
  pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
  sgeom = tgeom + pgeom
  mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
  mol.build( unit = 'Bohr')

  csfs = process_xml_csf(xmlfile)
  ncsfs = len(csfs)

  for i in range(ncsfs):
      print(i,csfs[i])

  if(sys.argv[1]=='0'):
      sys.exit()

  # Asymptotic Energies
  phase = np.ones(ncsfs)
  ovl, kin, pot = hcore(mol,smo)
  eri = twoeints(mol,smo)
  eecore = 2.0*np.trace(eri[0:tdoc_frozen,0:tdoc_frozen,:,:])
  h1e = kin + pot + eecore

  hmat, smat = cimat(ne,nmo,ovl,h1e,eri,csfs,phase)
  eig, eigv = np.linalg.eig(hmat)
  idx = eig.argsort()[::-1]
  eig = eig[idx]
  eigv = eigv[:,idx]
  nsta = len(eig)


  esta = []
  net = []
  nep = []
  print()
  print("Asymptotic energies")
  print()

  # Get the absolute values of the eigenvectors
  abs_eigv = np.abs(eigv)
  # Find the index of the largest component for each eigenvector (column)
  largest_component_indices = np.argmax(abs_eigv, axis=0)
  # Find the values of the largest components
  largest_component_values = abs_eigv[largest_component_indices, range(eigv.shape[1])]


  for i in range(nsta):
    _, alpe, betae =  csfs[largest_component_indices[i]].terms[0]
    nte = int(np.count_nonzero(np.array(alpe)<ntmo) + np.count_nonzero(np.array(betae)<ntmo))
    npe = len(alpe)+len(betae)-nte
    net.append(nte)
    nep.append(npe)
    #esta.append(hmat[i,i].real)
    esta.append(eig[i])
    print(i, esta[i].real)#, "          ", alpe, betae, nte, npe, largest_component_indices[i],largest_component_values[i])
    print(*[ (j, eigv[j,i].real) for j in range(ncsfs)])
    print()
    #print(i, esta[i], "          ", csf)

  i_init = int(input('Enter the initial state : '))


  for b in blist:

   tmat = []
   for zproj in zlist:

     time = zproj/vproj
     phase = []
     for i, csf in enumerate(csfs):
       if(nep[i]==0):
         phase.append(1.0)
       elif(nep[i]==1):
         #phase.append(1.0)
         phase.append(np.exp(-vproj*zproj*1.0j)*np.exp(+0.5*vproj**2*time*1.0j))
       elif(nep[i]==2):
         #phase.append(1.0)
         phase.append(np.exp(-vproj*zproj*1.0j)**2*np.exp(+vproj**2*time*1.0j))
       else:
         raise NotImplementedError("Wrong number of projectile electrons")

     xp = b
     yp = 0
     zp = zproj
     pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
     sgeom = tgeom + pgeom
     mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
     mol.build( unit = 'Bohr')

     ovl, kin, pot = hcore(mol,smo)
     eri = twoeints(mol,smo)
     eecore = 2.0*np.trace(eri[0:tdoc_frozen,0:tdoc_frozen,:,:])
     h1e = kin + pot + eecore

     matH, matS = cimat(ne,nmo,ovl,h1e,eri,csfs,phase)

     hmat = np.linalg.inv(eigv) @ matH @ eigv
     smat = np.linalg.inv(eigv) @ matS @ eigv
     inv = np.linalg.inv(smat)
     mat = np.matmul(inv, hmat)
     tmat.append(mat)

   # running dynamics for bproj = b
   tlist = zlist/vproj
   hmat_interp = interp1d(tlist, tmat , axis=0)
   psi0 = np.zeros(ncsfs, dtype=complex)
   psi0[i_init] = 1.0
   ntime = int(2.0*zmax/(vproj*dtime))
   t_grid = np.linspace(zlist[0]/vproj,zlist[-1]/vproj,ntime)
   wf_t = solve_tdse(hmat_interp, psi0, t_grid)
   prob = np.abs(wf_t[-1,:])**2
   formatted_string = "  ".join([f"{num:.6f}" for num in prob])
   print(b, formatted_string,' ', f"{np.sum(prob):.6f}")
   #print(b,*prob,np.sum(prob))

   #analyzing the dynamics for bproj = b
   if analyze == True:
    fout = open('prob_time_'+str(vproj)+'_'+str(b)+'.txt','w')
    for it, time in enumerate(t_grid):
       zproj = vproj * time
       if(it%nstep_analysis==0):
        amp2 = np.abs(wf_t[it])**2
        print(zproj,' '.join(map(str,amp2)),np.sum(amp2),file=fout)
    fout.close()
    # can do 1rdm analysis here

    ztime = []
    stime = []
    for it, time in enumerate(t_grid):
       zproj = vproj * time
       if(it%nstep_analysis==0):

        ztime.append(zproj)
        xp = b
        yp = 0
        zp = zproj
        pgeom = elp + " " + str(xp) + " " + str(yp) + " " + str(zp)
        sgeom = tgeom + pgeom
        mol = gto.M(atom=sgeom,basis=sbasis,charge=scharge,spin=sspin)
        mol.build( unit = 'Bohr')

        ovl, kin, pot = hcore(mol,smo)
        cicoeffs = wf_t[it]

        #rdm1, rdm2 =  one_rdm_nonorth(csfs, cicoeffs, ovl)
        rdm1 =  one_rdm_nonorth(csfs, cicoeffs, ovl)
        rrdm1 = np.real(rdm1)
        #rrdm2 = np.real(rdm2)

        eigenvalues, eigenvectors = np.linalg.eigh(rdm1)
        # Compute von Neumann entropy
        entropy = -np.sum(eigenvalues * np.log(eigenvalues + 1e-12))  # Add small value to avoid log(0)

        print(zproj,entropy)
        #stime.append(entropy)
        #print(zproj,rrdm1[0,0],rrdm1[1,1],rrdm1[3,3],rrdm1[4,4],rrdm1[2,2],rrdm1[5,5])
        #print(np.real(rdm1))
        #print(np.trace(np.real(rdm1)))
        #print()

        #print("Diagonal elements of the 2-RDM (pair populations):")
        #for p in range(rdm2.shape[0]):
        # for q in range(rdm2.shape[1]):
        #  print(f"Γ({p},{q},{p},{q}) = {rdm2[p, q, p, q]:.3f}")
        # Sum over p for each q
        #sum_over_p = np.zeros(2*nmo)  # Initialize for each q
        #for q in range(2*nmo):
        #    for p in range(0,2):
        #        sum_over_p[q] += rdm2[p, q, p, q]

        #print(zproj,sum_over_p[2],sum_over_p[5])
        #print(zproj,rrdm2[3,2,3,2],rrdm2[4,2,4,2],rrdm1[3,3],rrdm1[4,4],rrdm1[2,2])
        #print(zproj,' '.join(map(str,np.abs(rdm2[0,2,4,:]))))


    #plt.plot(ztime, stime)
    #plt.show()

