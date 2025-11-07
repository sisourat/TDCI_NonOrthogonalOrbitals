import sys
import numpy as np


dat = np.loadtxt(sys.argv[1])
bimp = dat[:,0]
nbimp = len(bimp)
nsta = len(dat[0,:])-2
bprob = np.zeros((nbimp,nsta))
for ib in range(nbimp):
    for i in range(nsta):
        bprob[ib,i] = bimp[ib]*(dat[ib,1+i]/dat[ib,-1])

sig = []
for i in range(nsta):
 sig.append(np.trapz(bprob[:,i],bimp)*2.0*np.pi*0.28)
 print(i,sig[i])

print()
print(np.sum(sig[1:5]))

