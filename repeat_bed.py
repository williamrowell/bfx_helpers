#!/usr/bin/env python3
"""Identify low-complexity regions in a fasta file.

Outputs a bed file with annotated homopolymers and dinucleotide streches.
"""
import argparse
import csv
import re
from multiprocessing import Pool, cpu_count

from Bio import SeqIO
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("reference", help="path for input fasta", type=str)
parser.add_argument(
    "--minhrun", help="minimum homopolymer run length (0=disabled)", type=int, default=0
)
parser.add_argument(
    "--minrepeats",
    help="minimum dinucleotide repeats (0=disabled)",
    type=int,
    default=0,
)
parser.add_argument(
    "--threads",
    help="number of threads (1)",
    type=int,
    default=1,
    choices=range(1, cpu_count()),
)
args = parser.parse_args()

if not (args.minhrun or args.minrepeats):
    parser.error("No action requested, add --minhrun or --minrepeats")


# generate regex
BASES = IUPAC.IUPACUnambiguousDNA.letters
DINUC_REPEATS = [x + y for x in BASES for y in BASES if x != y]
PATTERN = ""
if args.minrepeats:
    PATTERN = ("|").join(["(%s){%s,}" % (x, args.minrepeats) for x in DINUC_REPEATS])
if args.minhrun:
    if PATTERN:
        PATTERN += "|"
    PATTERN = ("|").join(["(%s){%s,}" % (x, args.minhrun) for x in BASES])
print(PATTERN)


def find_repeats(seq_record, pattern=PATTERN):
    """Find and record all occurrences of PATTERN in seq_record.
    Case-insensitive.

    Outputs list of 4-tuples.
        Example:
        [('chr1', 11019, 11026, 'Gx7'),
        ('chr1', 36793, 36803, 'TGx5'),...]
    """
    repeat_list = list()
    for repeat in re.finditer(pattern, str(seq_record.seq), flags=re.I):
        unit = [x for x in repeat.groups() if x][0].upper()  # 'AC'
        size = len(repeat.group()) / len(unit)  # len('ACACAC')/len('AC')
        repeat_list.append(
            (seq_record.id, repeat.start(), repeat.end(), unit + "x" + str(size))
        )
    return repeat_list


p = Pool(args.threads)
repeat_list = p.map(find_repeats, SeqIO.parse(args.reference, "fasta"))
repeat_list = [item for items in repeat_list for item in items]  # unpack

outfile = (
    args.reference
    + ".dinuc_repeats_ge"
    + str(args.minrepeats)
    + ".hrun_ge"
    + str(args.minhrun)
    + ".bed"
)

with open(outfile, "wb") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerows(repeat_list)
