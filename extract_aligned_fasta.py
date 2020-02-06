#!/usr/bin/env python

"""Extract subsequences mapped to contig:start-stop.
Given an aligned bam and coordinates, return the reference-oriented
subsequences of alignments that completely overlap [start, stop].
"""

import os.path
import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import pysam


direction = {True: 'rc', False: 'f'}


def get_seq_record(record, start, stop, description):
    """Return a SeqRecord for query between start and stop.
    Given a sam record, find query sequences that cover the [start, stop]
    interval completely and create a SeqRecord object with sequence in the
    reference orientation.
    """
    # get the query positions of the bases mapped to start, stop
    # TODO: handle cases where no base is mapped (None)
    positions = record.get_aligned_pairs(matches_only=True)
    first_position = [item[0] for item in positions if item[1]==start]
    last_position = [item[0] for item in positions if item[1]==stop]
    # fetch and reorient sequence
    if first_position and last_position:
        name = record.query_name
        if not record.is_reverse:
            seq = Seq(record.get_forward_sequence()[first_position[0]:last_position[0]],
                      IUPAC.IUPACUnambiguousDNA())
            direction='f'
        else:
            length = record.query_length
            seq = Seq(record.get_forward_sequence()[length - last_position[0]:length - first_position[0]],
                      IUPAC.IUPACUnambiguousDNA()).reverse_complement()
            direction = 'rc'
        return SeqRecord(seq, id=name, description=('|').join([description, direction[record.is_reverse]]))


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inbam', help="input bam", type=argparse.FileType('r'))
    parser.add_argument('outfasta', help="output fasta", type=argparse.FileType('w'))
    parser.add_argument('contig', help="reference contig", type=str)
    parser.add_argument('start', help="start position on reference contig", type=int)
    parser.add_argument('stop', help="stop position on reference contig", type=int)

    args = parser.parse_args(arguments)

    description = args.contig + ':' + str(args.start) + '-' + str(args.stop)
    my_records = []
    with pysam.AlignmentFile(args.inbam, 'rb') as samfile:
        for record in samfile.fetch(contig=args.contig, start=args.start, stop=args.stop):
            my_records.append(get_seq_record(record, args.start, args.stop, description))
    my_records = [record for record in my_records if record is not None]
    if my_records:
        SeqIO.write(sorted(my_records, key=len), args.outfasta, 'fasta')


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

