#!/usr/bin/env python

"""Add binned mq tags to BAM.
"""

import pysam


with pysam.AlignmentFile("-", "r") as infile, \
     pysam.AlignmentFile("-", "w", template=infile) as outfile:
    for s in infile:
        s.set_tag('mq', s.mapping_quality / 10, value_type='i')
        outfile.write(s)

