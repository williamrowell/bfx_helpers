#!/usr/bin/env python3
"""Identify all RE sites in reference.

Find all RE sites for ENZYME within REFERENCE and output in bed format.
"""

import argparse
import sys
import os.path

from Bio import Restriction
from Bio import SeqIO

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('reference', help='reference fasta', type=str)
parser.add_argument('enzyme', help='restriction enzyme', type=str)
args = parser.parse_args()


def main(args):
    if getattr(Restriction, args.enzyme):
        enzyme = getattr(Restriction, args.enzyme)
    else:
        sys.exit("%s not found" % args.enzyme)
    offsets = (enzyme.charac[0] + 1, enzyme.size - enzyme.charac[0] - 1)
    with open(('.').join([os.path.basename(args.reference), args.enzyme, 'bed']), 'w') as f:
        for record in SeqIO.parse(args.reference, 'fasta'):
            for site in enzyme.search(record.seq):
                f.write('%s\t%s\t%s\t%s\n' % (record.id, site - offsets[0], site + offsets[1], args.enzyme))


if __name__ == '__main__':
    main(args)

