import numpy as np


def zShift(seq, pos):
    """Return components of Z curve shift.

    zCurve[0] = (A+G)-(C+T) # purine/pyrimidine
    zCurve[1] = (A+C)-(G+T) # amino/keto
    zCurve[2] = (A+T)-(G+C) # weak/strong
    """
    if seq[pos] == "A":
        return np.array([1, 1, 1])
    if seq[pos] == "G":
        return np.array([1, -1, -1])
    if seq[pos] == "C":
        return np.array([-1, 1, -1])
    if seq[pos] == "T":
        return np.array([-1, -1, 1])


def zCurve(seq):
    """Return 3-dimensional Z curve corresponding to sequence.

    zcurve[n] = zcurve[n-1] + zShift[n]
    """
    zcurve = np.zeros((len(seq), 3), dtype=int)
    zcurve[0] = zShift(seq, 0)
    for pos in range(1, len(seq)):
        zcurve[pos] = np.add(zcurve[pos - 1], zShift(seq, pos))
    return zcurve
