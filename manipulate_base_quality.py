#!/usr/bin/env python

"""Manipulate fastq quality records.

For each record in a fastq, shuffle the quality scores, effectively
randomizing per-base quality while holding mean read quality constant.

Alternately, set all base qualities to a fixed value.
"""

import argparse
import random
import sys

from Bio import SeqIO


def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("fastq", help="file to process", type=argparse.FileType("r"))
    parser.add_argument("-s", "--seed", help="random seed (23)", type=int, default=23)
    parser.add_argument(
        "-f",
        "--fixed_score",
        help="Instead of randomizing, used fixed score (phred33)",
        type=int,
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file",
        default=sys.stdout,
        type=argparse.FileType("w"),
    )

    args = parser.parse_args(arguments)

    if args.fixed_score and args.fixed_score not in range(1, 99):
        sys.exit("Fixed score must be an integer in [1,99].")

    random.seed(args.seed)

    for record in SeqIO.parse(args.fastq, "fastq"):
        if args.fixed_score:
            record.letter_annotations["phred_quality"] = [args.fixed_score] * len(
                record.letter_annotations["phred_quality"]
            )
        else:
            random.shuffle(record.letter_annotations["phred_quality"])
        SeqIO.write(record, args.outfile, "fastq")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
