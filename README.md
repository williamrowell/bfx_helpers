# bfx_helpers

A junk drawer of one-off scripts that don't deserve their own repo yet.

## repeat_bed.py

Identify low-complexity sequence in a fasta file and output a bed file with annotated homopolymers and dinucleotide repeats.

```text
usage: repeat_bed.py [-h] [--minhrun MINHRUN] [--minrepeats MINREPEATS]
                     [--threads {1,2,3,4,5,6,7}]
                     reference

Identify low-complexity regions in a fasta file. Outputs a bed file with
annotated homopolymers and dinucleotide streches.

positional arguments:
  reference             path for input fasta

optional arguments:
  -h, --help            show this help message and exit
  --minhrun MINHRUN     minimum homopolymer run length (0=disabled)
  --minrepeats MINREPEATS
                        minimum dinucleotide repeats (0=disabled)
  --threads {1,2,3,4,5,6,7}
                        number of threads (1)
```

## add_mq_tag.py

Add tag `mq` to reads in BAM, with a binned representation of the MAPQ (MAPQ/10).

```text
usage: samtools view my.bam | python3 add_mq_tag.py > my.mq.bam
```

## add_vaf.py

Add FORMAT field for variant allele fraction to pbsv VCF.

```text
usage: cat pbsv.vcf | python3 add_vaf.py > pbsv.vaf.vcf
```

## depth_mean_stddev.py

For a mosdepth bed file provided to stdin, determine the mean and stddev as
if the data were normally distributed, and report these values.

```text
usage:  tabix ${BAM%.*}.median.regions.bed.gz <list of autosomes> | python depth_mean_stddev.py
```

## extract_aligned_fasta.py

```text
usage: extract_aligned_fasta.py [-h] inbam outfasta contig start stop

Extract subsequences mapped to contig:start-stop.
Given an aligned bam and coordinates, return the reference-oriented
subsequences of alignments that completely overlap [start, stop].

positional arguments:
  inbam       input bam
  outfasta    output fasta
  contig      reference contig
  start       start position on reference contig
  stop        stop position on reference contig

optional arguments:
  -h, --help  show this help message and exit
```

## extract_aligned_fastq.py

```text
usage: extract_aligned_fastq.py [-h] inbam outfastq contig start stop

Extract subsequences mapped to contig:start-stop.
Given an aligned bam and coordinates, return the subsequences of alignments
that completely overlap [start, stop].

positional arguments:
  inbam       input bam
  outfastq    output fastq
  contig      reference contig
  start       start position on reference contig
  stop        stop position on reference contig

optional arguments:
  -h, --help  show this help message and exit
```

## gen_fp_fn_tracks.py

```text
usage: gen_fp_fn_tracks.py [-h] invcf

Split hap.py annotated vcf into variant types and false call types.

positional arguments:
  invcf       hap.py vcf with annotated FP and FN

optional arguments:
  -h, --help  show this help message and exit
```

## markdup_stream.py

Mark duplicate reads in input BAM and write to output BAM. Reference based.

```text
usage:
python markdup_stream.py --inbam consensusalignments.bam \
                         --outbam markdup.consensusalignments.bam \
                         --stats dedup.stats.txt \
                         --dup_list dupreads.txt
```

## REbed.py

Mark duplicate reads in input BAM and write to output BAM

```text
usage:
python markdup_stream.py --inbam consensusalignments.bam \
                         --outbam markdup.consensusalignments.bam \
                         --stats dedup.stats.txt \
                         --dup_list dupreads.txt
```

## missing_adaptors.py

Calculate alignment of sequences flanking adaptor and plot distribution. CCS
reads from SMRTbell templates missing one adaptor should result in palindromic
reads. This script attempts to identify these reads by aligning bases that
flank the adaptor. The outputs are a histogram of normalized alignment score
(so you can pick a threshold that fits your data) and a table with alignment
scores and other useful per-read stats, so that you can easily create
whitelists. My early attempts at aligning the full CCS read to its reverse
complement 1) resulted in distributions that were difficult to explain and 2)
took too long for large datasets. Comparing the first and last N bases in the
CCS read solved both of those problems. The tradeoff is that this method might
be less robust to biological repeats. One can always change the number of
flanking bases used if desired. The constants at the top are empirically
chosen for bimodal separation after testing on a couple datasets. YMMV.

```text
usage: missing_adaptor.py [-h] [--nohist] fasta
```

## vcf2hetbed.py

Given a vcf, create a bed of intervals between consecutive heterozygous variants.

```text
usage: vcfhet2bed.py [-h] vcf chromosome

Convert 1KGP VCF to BED of het var intervals

positional arguments:
  vcf         1KGP vcf
  chromosome  chromosome

optional arguments:
  -h, --help  show this help message and exit
```

## hp_compress.py & zCurve.py

Helper functions for compressing homopolymers and calculating a DNA zCurve.
