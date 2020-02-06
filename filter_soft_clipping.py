#!/usr/bin/env python3
"""Mark duplicate reads in input BAM and write to output BAM
example usage:
python markdup_stream.py --inbam consensusalignments.bam \
                         --outbam markdup.consensusalignments.bam \
                         --stats dedup.stats.txt \
                         --dup_list dupreads.txt
"""
import os
import sys
import argparse
import pysam
from collections import deque


def _parse_args(cmdl):
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--inbam", "-i", help="input BAM (/dev/stdin)",
                        type=argparse.FileType("r"), default=sys.stdin)
    parser.add_argument("--outbam", "-o", help="output BAM (/dev/null)",
                        type=argparse.FileType("w"), default=os.devnull)
    parser.add_argument("--stats", "-s", help="output stats (/dev/stderr)",
                        type=argparse.FileType("w"), default=sys.stderr)
    parser.add_argument("--rmdup", "-r", help="remove duplicates", action='store_true')
    parser.add_argument("--aln_wiggle", "-a", help="bp alignment wiggle at ends (2)",
                        type=int, default=2)
    parser.add_argument("--len_wiggle", "-l", help="fraction read length wiggle (0.01)",
                        type=int, default=0.01)
    parser.add_argument("--dup_list", "-d", help="output duplicate read names",
                        type=argparse.FileType("w"))
    return parser.parse_args(cmdl)


def ref_match(curr, prev):
    """curr is aligned to same reference id as prev"""
    return curr.reference_id == prev.reference_id


def start_match(curr, prev, aln_wiggle):
    """curr starts within aln_wiggle bp of prev"""
    return abs(curr.reference_start - prev.reference_start) <= aln_wiggle


def end_match(curr, prev, aln_wiggle):
    """curr ends within aln_wiggle bp of prev"""
    return abs(curr.reference_end - prev.reference_end) <= aln_wiggle


def length_match(curr, prev, len_wiggle):
    """curr infer_query_length within len_wiggle bp of prev"""
    return abs(curr.infer_query_length() - prev.infer_query_length()) <= \
        len_wiggle * curr.infer_query_length()


def is_duplicate(curr, buffer, aln_wiggle, len_wiggle):
    """Return True if current alignment is a duplicate of alignments in buffer.
    Reference id must be identical.
    Start and end coordinates must be similar.
    Length must be similar.
    """
    for prev in buffer:
        if all([ref_match(curr, prev),
                start_match(curr, prev, aln_wiggle),
                end_match(curr, prev, aln_wiggle),
                length_match(curr, prev, len_wiggle)]):
            return True
    return False


def flush_buffer(curr, buffer, aln_wiggle):
    """Remove reads that are out of range of curr.
    buffer is modified
    removed reads are returned
    """
    # check oldest item in buffer until it's within aln_wiggle
    to_write = list()
    while buffer and not (ref_match(curr, buffer[0]) and start_match(curr, buffer[0], aln_wiggle)):
        to_write.append(buffer.popleft())
    return to_write


def write_alignments(alignments, rmdup, outfile):
    """Write alignments to outfile.
    If rmdup flag is set, don't write duplicate reads.
    """
    for alignment in alignments:
        if not (rmdup and alignment.is_duplicate):
            outfile.write(alignment)


def write_stats(duplicate_set, total_set, args):
    """Write stats to args.stats
    """
    values = [('input', args.inbam.name),
              ('output', args.outbam.name),
              ('alignment wiggle', '{} bp'.format(args.aln_wiggle)),
              ('length wiggle', '{:.02}'.format(args.len_wiggle)),
              ('total reads', len(total_set)),
              ('PCR duplicate reads', len(duplicate_set)),
              ('PCR duplicate frequency', '{:.04f}%'.format(100 * len(duplicate_set)/float(len(total_set))))]
    if args.stats:
        args.stats.writelines(f'{name}	{value}\n' for name, value in values)


def main(arguments):
    """Flag duplicates in input BAM and write to output BAM.
    """
    args = _parse_args(arguments)
    with pysam.AlignmentFile(args.inbam, "rb") as infile, \
            pysam.AlignmentFile(args.outbam, "wb", template=infile) as outfile:
        buffer = deque()  # buffer of previous reads within aln_wiggle
        duplicate_set = set()  # set of query_name for all duplicate reads
        total_set = set()  # set of query_name for all reads

        for curr in infile:
            # flush output buffer of any alignments out of range of curr and write to disk
            write_alignments(flush_buffer(curr, buffer, args.aln_wiggle), args.rmdup, outfile)

            # clear curr is_duplicate flag
            curr.is_duplicate = False

            # mark the current alignment as duplicate if it is:
            #   1. from a known duplicat read, or
            #   2. a duplicate of an alignment in the buffer
            if curr.query_name in duplicate_set or \
                    is_duplicate(curr, buffer, args.aln_wiggle, args.len_wiggle):
                curr.is_duplicate = True
                duplicate_set.add(curr.query_name)

            # add curr to buffer and set of all query names
            buffer.append(curr)
            total_set.add(curr.query_name)
        # write alignments remaining in output buffer
        write_alignments(buffer, args.rmdup, outfile)

        # print the duplicate rate to stdout
        write_stats(duplicate_set, total_set, args)

        # write list of duplicate reads
        if args.dup_list:
            args.dup_list.writelines("%s\n" % read_name for read_name in list(duplicate_set))


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

