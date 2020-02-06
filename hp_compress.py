def hpCompress(seq):
    """Return hompolymer compressed string.
    """
    hpcompressed = seq[0]
    for pos in range(1,len(seq)):
        if seq[pos] != hpcompressed[-1]:
            hpcompressed += seq[pos]
    return hpcompressed
