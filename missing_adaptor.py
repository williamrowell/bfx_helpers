#!/usr/local/env python3
"""Calculate alignment of sequences flanking adaptor and plot distribution.

CCS reads from SMRTbell templates missing one adaptor should result in
palindromic reads.  This script attempts to identify these reads by
aligning bases that flank the adaptor.  The outputs are a histogram of
normalized alignment score (so you can pick a threshold that fits your
data) and a table with alignment scores and other useful per-read stats,
so that you can easily create whitelists.

My early attempts at aligning the full CCS read to its reverse complement
1) resulted in distributions that were difficult to explain and 2) took
too long for large datasets.  Comparing the first and last N bases in the
CCS read solved both of those problems.  The tradeoff is that this method
might be less robust to biological repeats.  One can always change the
number of flanking bases used if desired.

The constants at the top are empirically chosen for bimodal separation
after testing on a couple datasets. YMMV.
"""

import os.path
import sys
import argparse
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import pairwise2, SeqIO


matplotlib.use("Agg")

FLANK_LEN = 100  # number of adaptor-flanking bases to align
MATCH = 5  # match score
MISMATCH = -4  # mismatch penalty
GAP_OPEN = -10  # gap open penalty
GAP_EXTEND = -0.5  # gap extend penalty


def calc_identity(record):
    """Calculate the percent identity of adaptor flanking sequences.

    refer to http://biopython.org/DIST/docs/api/Bio.pairwise2-module.html
    """
    leftflank = record.seq[0:FLANK_LEN]
    rightflank = record.seq[-FLANK_LEN:].reverse_complement()
    alignment = pairwise2.align.globalms(
        leftflank, rightflank, MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND, score_only=True
    )
    return alignment, (float(alignment) / (MATCH * FLANK_LEN)) * 100


def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="input fasta", type=str)
    parser.add_argument("--nohist", help="suppress histogram", action="store_false")
    args = parser.parse_args(arguments)

    outcsv = os.path.splitext(args.fasta)[0] + "_self_alignment.csv"

    data = defaultdict(list)
    for record in SeqIO.parse(args.fasta, "fasta"):
        data["name"].append(record.id)
        data["zmw"].append(record.id.split("/")[1])
        data["length"].append(len(record))
        score, percent = calc_identity(record)
        data["score"].append(score)
        data["percent"].append(percent)
    df = pd.DataFrame(data=data)
    df = df[["name", "zmw", "length", "score", "percent"]]

    fig, axarr = plt.subplots(2)
    fig.set_size_inches(8, 6)
    axarr[0].set_title(args.fasta)
    bins = np.linspace(df.percent.min(), df.percent.max(), 50)

    if args.nohist:
        # plot alignment distribution
        df.percent.hist(ax=axarr[0], bins=bins, log=True)
        axarr[0].set_ylabel("log count")

        # plot cumulative distribution
        df.percent.hist(ax=axarr[1], bins=bins, cumulative=True, density=True)
        axarr[1].set_ylabel("cum. freq.")
        axarr[-1].set_xlabel("percent alignment")

        plt.tight_layout()
        outfig = os.path.splitext(args.fasta)[0] + "_self_alignment_distribution.png"
        plt.savefig(outfig)

    # save table
    df.to_csv(outcsv, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
