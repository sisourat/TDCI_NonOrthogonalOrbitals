import numpy as np
from dataclasses import dataclass
import sys

from numba import jit

@dataclass
class Sdeterminant:
    nalpha: int
    nbeta: int
    alpha: list
    beta: list

@jit
def compute_det(mat):
    return np.linalg.det(mat)

def lowdin(ne, nmo, ovmo, h1emo, r12mo, det1: Sdeterminant, det2: Sdeterminant):
    ovmat = np.zeros((ne, ne), dtype=complex)
    ovstore = np.zeros_like(ovmat)

    # Overlap matrix between determinants
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            ovmat[j, i] = ovmo[ja, ia]
            #print('s alpha',j,i,ia,ja,ovmo[ja, ia])

    for i in range(det1.nbeta):
        ib = det1.beta[i]
        for j in range(det2.nbeta):
            jb = det2.beta[j]
            #print('s beta',j + det2.nalpha, i + det1.nalpha,ovmo[jb, ib])
            ovmat[j + det2.nalpha, i + det1.nalpha] = ovmo[jb, ib]

    ovstore[:, :] = ovmat[:, :]
    ov = compute_det(ovmat)

    # One-electron term
    h1e = 0 + 0j
    comat = np.zeros((ne - 1, ne - 1), dtype=complex)

    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            ovmat[:, :] = ovstore[:, :]
            comat[:, :] = np.delete(np.delete(ovmat, i, axis=1), j, axis=0)
            h1e += (-1) ** (i + j) * h1emo[ja, ia] * compute_det(comat)
            #print("alpha-alpha 1e",i,j,h1emo[ja, ia] , compute_det(comat), (-1) ** (i + j))

    for i in range(det1.nbeta):
        ib = det1.beta[i]
        for j in range(det2.nbeta):
            jb = det2.beta[j]
            ovmat[:, :] = ovstore[:, :]
            comat[:, :] = np.delete(np.delete(ovmat, i+det1.nalpha, axis=1), j+det2.nalpha, axis=0)
            h1e += (-1) ** (det1.nalpha + i + det2.nalpha + j) * h1emo[jb, ib] * compute_det(comat)
            #print("beta-beta 1e",i,j,h1emo[jb, ib] , compute_det(comat), (-1) ** (det1.nalpha + i + det2.nalpha + j))

    # Two-electron term
    r12 = 0 + 0j
    if(ne<2):
        return ov, h1e, 0.0
    comat2 = np.zeros((ne - 2, ne - 2), dtype=complex)

    # Alpha-alpha terms
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            for k in range(i+1,det1.nalpha):
                ka = det1.alpha[k]
                for l in range(j+1,det2.nalpha):
                    la = det2.alpha[l]
                    ovmat[:, :] = ovstore[:, :]
                    comat2[:, :] = np.delete(np.delete(ovmat, [i, k], axis=1),
                                              [j, l], axis=0)
                    r12 += (r12mo[la, ka, ja, ia]-r12mo[la, ja, ka, ia]) * (-1) ** (i + j + k + l) * compute_det(comat2)
                    #print("alpha-alpha 2e",r12mo[la, ka, ja, ia],r12mo[la, ja, ka, ia],compute_det(comat2),(ia , ja , ka , la))

    # Alpha-beta terms
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            for k in range(det1.nbeta):
                kb = det1.beta[k]
                for l in range(det2.nbeta):
                    lb = det2.beta[l]
                    ovmat[:, :] = ovstore[:, :]
                    #1comat2[:, :] = np.delete(np.delete(ovmat, [i, k+det1.nalpha], axis=0),
                    #1                          [j, l+det2.nalpha], axis=1)
                    comat2[:, :] = np.delete(np.delete(ovmat, [i, k+det1.nalpha], axis=1),
                                              [j, l+det2.nalpha], axis=0)
                    r12 += r12mo[lb, kb, ja, ia] * (-1) ** (i + j + k + l + det1.nalpha + det2.nalpha) * compute_det(comat2)

                    #print("alpha-beta 2e",r12mo[lb, kb, ja, ia], compute_det(comat2), (-1) ** (i + j + k + l + det1.nalpha + det2.nalpha), lb, kb, ja, ia)

    # Beta-beta terms
    for i in range(det1.nbeta):
        ib = det1.beta[i]
        for j in range(det2.nbeta):
            jb = det2.beta[j]
            for k in range(i+1,det1.nbeta):
                kb = det1.beta[k]
                for l in range(j+1,det2.nbeta):
                    lb = det2.beta[l]
                    ovmat[:, :] = ovstore[:, :]
                    comat2[:, :] = np.delete(np.delete(ovmat, [k+det1.nalpha, i+det1.nalpha], axis=1),
                                              [l+det2.nalpha,j+det2.nalpha], axis=0)
                    r12 += (r12mo[lb, kb, jb, ib] - r12mo[lb, jb, kb, ib]) * (-1) ** (i + j + k + l + 2 * det1.nalpha + 2 * det2.nalpha) * compute_det(comat2)
                    #print("beta-beta 2e",r12mo[lb, kb, jb, ib], r12mo[lb, jb, kb, ib], compute_det(comat2), (-1) ** (i + j + k + l + 2 * det1.nalpha + 2 * det2.nalpha))

    #print(h1e,r12)
    #print()
    return ov, h1e, r12

if __name__ == "__main__":
    # Example test case
    ne = 4
    nmo = 5
    nalpha = 2
    nbeta = 2

    print("\n===== Test 1: Random Non-Orthogonal Orbitals =====")
    ovmo = np.random.rand(nmo, nmo) + 1j * np.random.rand(nmo, nmo)
    h1emo = np.random.rand(nmo, nmo) + 1j * np.random.rand(nmo, nmo)
    r12mo = np.random.rand(nmo, nmo, nmo, nmo) + 1j * np.random.rand(nmo, nmo, nmo, nmo)

    det1 = Sdeterminant(nalpha, nbeta, [0, 4], [2, 3])
    det2 = Sdeterminant(nalpha, nbeta, [0, 3], [2, 3])

    ov, h1e, r12 = lowdin(ne, nmo, ovmo, h1emo, r12mo, det1, det2)

    print("Overlap:", ov)
    print("One-electron term:", h1e)
    print("Two-electron term:", r12)

    print("\n===== Test 2: Orthogonal Orbitals =====")
    # Make orthogonal orbitals (identity matrix for simplicity)
    ovmo = np.eye(nmo, dtype=complex)
    h1emo = np.eye(nmo, dtype=complex)
    r12mo = np.zeros((nmo, nmo, nmo, nmo), dtype=complex)

    for i in range(nmo):
        r12mo[i, i, i, i] = 1.0 + 0j  # simple test value

    det1 = Sdeterminant(nalpha, nbeta, [0, 1], [2, 3])
    det2 = Sdeterminant(nalpha, nbeta, [0, 1], [2, 3])

    ov, h1e, r12 = lowdin(ne, nmo, ovmo, h1emo, r12mo, det1, det2)

    print("Overlap:", ov)
    print("One-electron term:", h1e)
    print("Two-electron term:", r12)
