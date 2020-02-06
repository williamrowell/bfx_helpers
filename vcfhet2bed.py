#!/usr/bin/env python3

"""Convert 1KGP VCF to BED of het var intervals
"""

import sys
import re
import argparse
import csv
from collections import defaultdict


GT_DELIM = re.compile(r'/|\|')


def isHet(genotype):
    """Return True if genotype is heterozygous.
    
    >>> isHet('1|0')
    True
    >>> isHet('0/1')
    True
    >>> isHet('0|0')
    False
    >>> isHet('0')
    False
    """
    if re.search(GT_DELIM, genotype):
        haplotypes = re.split(GT_DELIM, genotype)
        return not haplotypes[0] == haplotypes[1]
    else:
        # genotypes without / or | are hemizygous
        return False


def parse_vcf(vcffile):
    """Return dict of heterozygous variant positions from vcf."""
    samples = []  # list of sample names
    hetdict = defaultdict(list)  # dictionary of het var positions
    for row in open(vcffile, 'rb'):
        fields = row.split()
        if not fields[0].startswith('#') and fields[2].startswith('rs'):
            # ignore metadata headers and CNVs
            genotypes = fields[9:]
            # add positions of heterozygous variants to lists per sample
            for i, g in enumerate(genotypes):
                if isHet(g):
                    hetdict[samples[i]].append(int(fields[1]))
        elif fields[0].startswith('#CHROM'):
            # sample names start at column 10 of header row
            samples = fields[9:]
            next
    return hetdict


def write_bed(bedname, hetlist, chrom):
    """Write bed of inter-het intervals.
    chromosome  start   end length
    """
    with open(bedname, 'wb') as bedfile:
        csvwriter = csv.writer(bedfile, delimiter='\t', lineterminator='\n')
        for x, y in zip(hetlist[0:-1], hetlist[1:]):
            csvwriter.writerow([chrom, x-1, y, y-x])


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help="1KGP vcf", type=str)
    parser.add_argument('chromosome', help="chromosome", type=str)
    args = parser.parse_args(arguments)

    hetdict = parse_vcf(args.vcf)
    for sample in hetdict.keys():
        if hetdict[sample]:
            write_bed(('.').join([args.chromosome, sample, 'bed']),
                                  hetdict[sample], args.chromosome)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

