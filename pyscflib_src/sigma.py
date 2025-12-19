import sys
import numpy as np


dat = np.loadtxt(sys.argv[1])
bimp = dat[:,0]
nbimp = len(bimp)
nsta = len(dat[0,:])-2
bprob = np.zeros((nbimp,nsta))
bpbiem = np.zeros((nbimp,nsta))
for ib in range(nbimp):
    for i in range(nsta):
        bprob[ib,i] = bimp[ib]*(dat[ib,1+i]/dat[ib,-1])
        p = dat[ib,1+i]/dat[ib,-1]
        piem = 2.0*p*(1.0-p)
        bpbiem[ib,i] = bimp[ib]*piem

sig = []
for i in range(nsta):
 sig.append(np.trapz(bprob[:,i],bimp)*2.0*np.pi*0.28)
 print(i,sig[i])

print()
print(sig[13],np.sum(sig[10:13]),np.sum(sig[5:10]))
#print(sig[5],np.sum(sig[6:9]),np.sum(sig[9:14]))
#print(sig[9]+sig[14])
#print(np.sum(sig[1:5]))
#print(np.sum(sig[9:18]))
sys.exit()

print()
print("IEM Model")
print()
sig = []
for i in range(nsta):
 sig.append(np.trapz(bpbiem[:,i],bimp)*2.0*np.pi*0.28)
 print(i,sig[i])

print()
print(sig[9]+sig[14])
#print(np.sum(sig[1:5]))
#print(np.sum(sig[9:18]))

