#!/usr/bin/env python3
"""
Add nhomalt to VCF.
"""

__author__ = "William Rowell"
__version__ = "0.1.0"
__license__ = "MIT"

import vcfpy


def nhomalt(record):
    """"""
    nhomalt = 0
    for call in record:
        if call.gt_alleles == [1, 1]:
            nhomalt += 1
    return nhomalt


def main():
    # read from stdin
    reader = vcfpy.Reader.from_path("/dev/stdin")
    # add 'nhomalt' to header
    reader.header.add_info_line(
        vcfpy.OrderedDict(
            [
                ("ID", "nhomalt"),
                ("Number", "A"),
                ("Type", "Integer"),
                (
                    "Description",
                    "The number of individuals that are called homozygous for the alternate allele.",
                ),
            ]
        )
    )
    # write to stdout with modified header
    with vcfpy.Writer.from_path("/dev/stdout", reader.header) as writer:
        for record in reader:
            record.INFO["nhomalt"] = [nhomalt(record)]
            writer.write_record(record)


if __name__ == "__main__":
    """This is executed when run from the command line"""
    main()
