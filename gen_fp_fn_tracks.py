#!/usr/bin/env python3
"""Split hap.py annotated vcf into variant types and false call types.
"""

import os
import sys
import argparse
import vcf


def filter_vcf(invcf):
    """Filter input vcf by type (SNP/INDEL) and false call type (FP/FN)."""
    vcf_reader = vcf.Reader(open(invcf, 'r'))
    vcf_writer_FP_INDEL = vcf.Writer(open('FP_INDEL.vcf', 'w'), vcf_reader)
    vcf_writer_FN_INDEL = vcf.Writer(open('FN_INDEL.vcf', 'w'), vcf_reader)
    vcf_writer_FP_SNP = vcf.Writer(open('FP_SNP.vcf', 'w'), vcf_reader)
    vcf_writer_FN_SNP = vcf.Writer(open('FN_SNP.vcf', 'w'), vcf_reader)
    for record in vcf_reader:
        if record.genotype('QUERY')['BVT'] == 'INDEL' and record.genotype('QUERY')['BD'] == 'FP':
            vcf_writer_FP_INDEL.write_record(record)
        elif record.genotype('TRUTH')['BVT'] == 'INDEL' and record.genotype('TRUTH')['BD'] == 'FN':
            vcf_writer_FN_INDEL.write_record(record)
        elif record.genotype('QUERY')['BVT'] == 'SNP' and record.genotype('QUERY')['BD'] == 'FP':
            vcf_writer_FP_SNP.write_record(record)
        elif record.genotype('TRUTH')['BVT'] == 'SNP' and record.genotype('TRUTH')['BD'] == 'FN':
            vcf_writer_FN_SNP.write_record(record)
    vcf_writer_FN_SNP.close()
    vcf_writer_FP_SNP.close()
    vcf_writer_FN_INDEL.close()
    vcf_writer_FP_INDEL.close()


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__)
    parser.add_argument('invcf', help="hap.py vcf with annotated FP and FN", type=str)
    args = parser.parse_args(arguments)
    filter_vcf(args.invcf)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

