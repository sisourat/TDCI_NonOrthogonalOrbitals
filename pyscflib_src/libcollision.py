import sys
import numpy as np

from pyscf import gto, scf, ao2mo, lo

def hcore(mol,mo):
   """ Computes the one-electron terms for given mol and mo """
   T = mol.intor('int1e_ovlp')
   ovl = mo.T @ T @ mo
   T = mol.intor('int1e_kin')
   kin = mo.T @ T @ mo
   T = mol.intor('int1e_nuc')
   pot = mo.T @ T @ mo
   return ovl, kin, pot

def twoeints(mol,mo):
    """ Compute the two electron integrals in MO basis """

# saves the two-electron integrals in the file ftmp.name
    ao2mo.kernel(mol, mo, erifile = 'hf.h5', dataname = 'test')
# load 2e integrals by filename and dataname
    with ao2mo.load('hf.h5', 'test') as eri:
      erimo = ao2mo.restore(1, np.asarray(eri), mo.shape[1])

    return erimo

def system(mol,debug=False):
  '''
    Computes target or projectile HF orbitals
  '''
  conv, e, mo_e, mo, mo_occ = scf.hf.kernel(scf.hf.SCF(mol), dm0=np.eye(mol.nao_nr()))
  nmo = len(mo_e)

 # mf = scf.RHF(mol).run()
 # lomo = lo.orth_ao(mf, 'nao')

  print()
  print(e)
  print()
  for i in range(len(mo_e)):
     print(i,mo_e[i])
  print()

  return mo, mo_e
