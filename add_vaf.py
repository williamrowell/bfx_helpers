#!/usr/bin/env python3
import vcfpy

# read from input file
reader = vcfpy.Reader.from_path("/dev/stdin")
# add 'VAF' to header
reader.header.add_format_line(
    vcfpy.OrderedDict(
        [
            ("ID", "VAF"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description", "Depth of variant alleles as a fraction of total depth."),
        ]
    )
)
# write to stdout using modified header
writer = vcfpy.Writer.from_path("/dev/stdout", reader.header)
# TODO: Handle multisample VCFs
for record in reader:
    vaf = sum(record.calls[0].data["AD"][1:]) / record.calls[0].data["DP"]
    record.add_format("VAF", round(vaf, 2))
    writer.write_record(record)
