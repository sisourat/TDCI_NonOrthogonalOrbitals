import numpy as np
from numba import njit, int64
from numba.experimental import jitclass

# -------------------- 1. Jitclass --------------------
spec = [
    ('nalpha', int64),
    ('nbeta', int64),
    ('alpha', int64[:]),
    ('beta', int64[:]),
]

@jitclass(spec)
class Sdeterminant:
    def __init__(self, nalpha, nbeta, alpha, beta):
        self.nalpha = nalpha
        self.nbeta = nbeta
        self.alpha = alpha
        self.beta = beta

# -------------------- 2. Determinant --------------------
@njit
def compute_det(mat):
    return np.linalg.det(mat)

# -------------------- 3. Direct copy function for deleting row/col --------------------
@njit
def copy_excluding(mat, out, row_indices, col_indices):
    """
    Copy mat into out, skipping rows and columns in row_indices and col_indices.
    row_indices and col_indices must be sorted ascending.
    """
    n_out, m_out = out.shape
    row_skip = set(row_indices)
    col_skip = set(col_indices)

    ri = 0
    for i in range(mat.shape[0]):
        if i in row_skip:
            continue
        ci = 0
        for j in range(mat.shape[1]):
            if j in col_skip:
                continue
            out[ri, ci] = mat[i, j]
            ci += 1
        ri += 1

# -------------------- 4. Lowdin function --------------------
@njit
def lowdin(ne, nmo, ovmo, h1emo, r12mo, det1, det2):
    ovmat = np.zeros((ne, ne), dtype=np.complex128)
    ovstore = np.zeros_like(ovmat)

    # Overlap matrix
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            ovmat[j, i] = ovmo[ja, ia]

    for i in range(det1.nbeta):
        ib = det1.beta[i]
        for j in range(det2.nbeta):
            jb = det2.beta[j]
            ovmat[j + det2.nalpha, i + det1.nalpha] = ovmo[jb, ib]

    ovstore[:, :] = ovmat[:, :]
    ov = compute_det(ovmat)

    # One-electron term
    h1e = 0.0 + 0.0j
    comat = np.zeros((ne - 1, ne - 1), dtype=np.complex128)

    # Alpha-alpha
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            ovmat[:, :] = ovstore[:, :]
            copy_excluding(ovmat, comat, np.array([j], dtype=np.int64), np.array([i], dtype=np.int64))
            h1e += ((-1) ** (i + j)) * h1emo[ja, ia] * compute_det(comat)

    # Beta-beta
    for i in range(det1.nbeta):
        ib = det1.beta[i]
        for j in range(det2.nbeta):
            jb = det2.beta[j]
            ovmat[:, :] = ovstore[:, :]
            copy_excluding(ovmat, comat,
                            np.array([j + det2.nalpha], dtype=np.int64),
                            np.array([i + det1.nalpha], dtype=np.int64))
            h1e += ((-1) ** (det1.nalpha + det2.nalpha + i + j)) * h1emo[jb, ib] * compute_det(comat)

    # Two-electron term
    r12 = 0.0 + 0.0j
    if ne < 2:
        return ov, h1e, r12

    comat2 = np.zeros((ne - 2, ne - 2), dtype=np.complex128)

    # Alpha-alpha 2e
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            for k in range(i+1, det1.nalpha):
                ka = det1.alpha[k]
                for l in range(j+1, det2.nalpha):
                    la = det2.alpha[l]
                    ovmat[:, :] = ovstore[:, :]
                    copy_excluding(ovmat, comat2,
                                    np.array([j, l], dtype=np.int64),
                                    np.array([i, k], dtype=np.int64))
                    r12 += (r12mo[la, ka, ja, ia] - r12mo[la, ja, ka, ia]) * \
                           ((-1) ** (i + j + k + l)) * compute_det(comat2)

    # Alpha-beta 2e
    for i in range(det1.nalpha):
        ia = det1.alpha[i]
        for j in range(det2.nalpha):
            ja = det2.alpha[j]
            for k in range(det1.nbeta):
                kb = det1.beta[k]
                for l in range(det2.nbeta):
                    lb = det2.beta[l]
                    ovmat[:, :] = ovstore[:, :]
                    copy_excluding(ovmat, comat2,
                                    np.array([j, l + det2.nalpha], dtype=np.int64),
                                    np.array([i, k + det1.nalpha], dtype=np.int64))
                    r12 += r12mo[lb, kb, ja, ia] * \
                           ((-1) ** (i + j + k + l + det1.nalpha + det2.nalpha)) * compute_det(comat2)

    # Beta-beta 2e
    for i in range(det1.nbeta):
        ib = det1.beta[i]
        for j in range(det2.nbeta):
            jb = det2.beta[j]
            for k in range(i+1, det1.nbeta):
                kb = det1.beta[k]
                for l in range(j+1, det2.nbeta):
                    lb = det2.beta[l]
                    ovmat[:, :] = ovstore[:, :]
                    copy_excluding(ovmat, comat2,
                                    np.array([j + det2.nalpha, l + det2.nalpha], dtype=np.int64),
                                    np.array([i + det1.nalpha, k + det1.nalpha], dtype=np.int64))
                    r12 += (r12mo[lb, kb, jb, ib] - r12mo[lb, jb, kb, ib]) * \
                           ((-1) ** (i + j + k + l + 2*det1.nalpha + 2*det2.nalpha)) * compute_det(comat2)

    return ov, h1e, r12

# -------------------- 5. Example usage --------------------
if __name__ == "__main__":
    ne = 4
    nmo = 5
    nalpha = 2
    nbeta = 2

    ovmo = np.random.rand(nmo, nmo) + 1j * np.random.rand(nmo, nmo)
    h1emo = np.random.rand(nmo, nmo) + 1j * np.random.rand(nmo, nmo)
    r12mo = np.random.rand(nmo, nmo, nmo, nmo) + 1j * np.random.rand(nmo, nmo, nmo, nmo)

    det1 = Sdeterminant(nalpha, nbeta,
                        np.array([0, 4], dtype=np.int64),
                        np.array([2, 3], dtype=np.int64))
    det2 = Sdeterminant(nalpha, nbeta,
                        np.array([0, 3], dtype=np.int64),
                        np.array([2, 3], dtype=np.int64))

    ov, h1e, r12 = lowdin(ne, nmo, ovmo, h1emo, r12mo, det1, det2)

    print("Overlap:", ov)
    print("One-electron term:", h1e)
    print("Two-electron term:", r12)

