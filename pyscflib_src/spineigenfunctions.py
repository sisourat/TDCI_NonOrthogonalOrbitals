class CSF:
    def __init__(self, nterms, terms):
        """
        nterms: int — number of terms in the CSF
        terms: list of tuples — each term is (coefficient: float, spin_a: list[int], spin_b: list[int])
        """
        self.nterms = nterms
        self.terms = terms

    def __repr__(self):
        return f"CSF(nterms={self.nterms}, terms={self.terms})"


def ne0spin1(paired, unpaired):
    if len(unpaired) != 0:
        raise ValueError("error in ne0spin1, incorrect unpaired length")

    terms = [
        (1.0, paired.copy(), paired.copy())
    ]
    return CSF(1, terms)


def ne1spin2(paired, unpaired):
    if len(unpaired) != 1:
        raise ValueError("error in ne1spin2, incorrect unpaired length")

    terms = [
        (1.0, paired.copy() + unpaired, paired.copy())
    ]
    return CSF(1, terms)


def ne2spin1(paired, unpaired):
    if len(unpaired) != 2:
        raise ValueError("error in ne2spin1, incorrect unpaired length")

    terms = [
        (0.70710678118654746, paired.copy() + [unpaired[0]], paired.copy() + [unpaired[1]]),
        (0.70710678118654746, paired.copy() + [unpaired[1]], paired.copy() + [unpaired[0]])
    ]
    return CSF(2, terms)


def ne2spin3(paired, unpaired):
    if len(unpaired) != 2:
        raise ValueError("error in ne2spin3, incorrect unpaired length")

    terms = [
        (1.0, paired.copy() + unpaired, paired.copy())
    ]
    return CSF(1, terms)


def ne3spin2(paired, unpaired):
    if len(unpaired) != 3:
        raise ValueError("error in ne3spin2, incorrect unpaired length")

    terms = [
        (0.70710678118654746, paired.copy() + [unpaired[0], unpaired[1]], paired.copy() + [unpaired[2]]),
        (-0.70710678118654746, paired.copy() + [unpaired[0], unpaired[2]], paired.copy() + [unpaired[1]]),
        (0.70710678118654746, paired.copy() + [unpaired[1], unpaired[2]], paired.copy() + [unpaired[0]])
    ]
    return CSF(3, terms)


def ne3spin4(paired, unpaired):
    if len(unpaired) != 3:
        raise ValueError("error in ne3spin4, incorrect unpaired length")

    terms = [
        (1.0, paired.copy() + unpaired, paired.copy())
    ]
    return CSF(1, terms)


def ne4spin1(paired, unpaired):
    if len(unpaired) != 4:
        raise ValueError("error in ne4spin1, incorrect unpaired length")

    terms = [
        (0.5, paired.copy() + [unpaired[0], unpaired[1]], paired.copy() + [unpaired[2], unpaired[3]]),
        (0.5, paired.copy() + [unpaired[0], unpaired[2]], paired.copy() + [unpaired[1], unpaired[3]]),
        (0.5, paired.copy() + [unpaired[0], unpaired[3]], paired.copy() + [unpaired[1], unpaired[2]]),
        (0.5, paired.copy() + [unpaired[1], unpaired[2]], paired.copy() + [unpaired[0], unpaired[3]]),
        (0.5, paired.copy() + [unpaired[1], unpaired[3]], paired.copy() + [unpaired[0], unpaired[2]]),
        (0.5, paired.copy() + [unpaired[2], unpaired[3]], paired.copy() + [unpaired[0], unpaired[1]])
    ]
    return CSF(6, terms)


def ne4spin3(paired, unpaired):
    if len(unpaired) != 4:
        raise ValueError("error in ne4spin3, incorrect unpaired length")

    terms = [
        (0.70710678118654746, paired.copy() + [unpaired[0], unpaired[1], unpaired[2]], paired.copy() + [unpaired[3]]),
        (-0.70710678118654746, paired.copy() + [unpaired[0], unpaired[1], unpaired[3]], paired.copy() + [unpaired[2]]),
        (0.70710678118654746, paired.copy() + [unpaired[0], unpaired[2], unpaired[3]], paired.copy() + [unpaired[1]]),
        (-0.70710678118654746, paired.copy() + [unpaired[1], unpaired[2], unpaired[3]], paired.copy() + [unpaired[0]])
    ]
    return CSF(4, terms)


def getCSF(spin, paired, unpaired):
    if spin == 1 and len(unpaired) == 0:
        return ne0spin1(paired, unpaired)
    elif spin == 2 and len(unpaired) == 1:
        return ne1spin2(paired, unpaired)
    elif spin == 1 and len(unpaired) == 2:
        return ne2spin1(paired, unpaired)
    elif spin == 3 and len(unpaired) == 2:
        return ne2spin3(paired, unpaired)
    elif spin == 2 and len(unpaired) == 3:
        return ne3spin2(paired, unpaired)
    elif spin == 4 and len(unpaired) == 3:
        return ne3spin4(paired, unpaired)
    elif spin == 1 and len(unpaired) == 4:
        return ne4spin1(paired, unpaired)
    elif spin == 3 and len(unpaired) == 4:
        return ne4spin3(paired, unpaired)
    else:
        return None

if __name__ == "__main__":
# Example test case
    csf = getCSF(1, [1, 2], [])
    print(csf)
    print(csf.nterms)
    print(csf.terms)

    csf2 = getCSF(1, [1, 2], [3, 4])
    print(csf2)
    print(csf.nterms)
    print(csf.terms)
