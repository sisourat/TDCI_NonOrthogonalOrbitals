from generate_csfs import *
import numpy as np
import sys

from lowdin_nb import *
from numba import jit, njit

# generate the CSFs

def cimat(ne,nmo,ovmo,h1emo,r12mo,csfs,phase):

 ncsfs = len(csfs)
 hmat = np.zeros((ncsfs,ncsfs), dtype=complex)
 smat = np.zeros((ncsfs,ncsfs), dtype=complex)

 imat1 = 0
 for i1, csf1 in enumerate(csfs):
    nterm1 = csf1.nterms
    imat2 = 0
    #print(csf1)
    for i2, csf2 in enumerate(csfs):
        nterm2 = csf2.nterms
        #print(csf2)
        for j1 in range(nterm1):
            c1, alp, bet = csf1.terms[j1]
            c1 = c1 * phase[i1]
            det1 = Sdeterminant(len(alp), len(bet), np.array(alp), np.array(bet))
            #det1 = Sdeterminant(len(alp), len(bet), alp, bet)
            for j2 in range(nterm2):
                 c2, alp, bet = csf2.terms[j2]
                 c2 = c2 * phase[i2]
                 det2 = Sdeterminant(len(alp), len(bet), np.array(alp), np.array(bet))
                 #det2 = Sdeterminant(len(alp), len(bet), alp, bet)
                 ov, h1e, r12 = lowdin(ne, nmo, ovmo, h1emo, r12mo, det1, det2)
                 hmat[imat1,imat2] += c1*np.conj(c2)*(h1e + r12)
                 smat[imat1,imat2] += c1*np.conj(c2)*ov
        imat2+=1
    imat1+=1

 #sys.exit()
 return hmat, smat
