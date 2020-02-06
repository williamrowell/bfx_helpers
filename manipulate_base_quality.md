Manipulate the per-base quality in fastq files to test the effect of bq
on a downstream process like variant calling.

```bash
$ rand_fq_qual.py --help
usage: rand_fq_qual.py [-h] [-s SEED] [-f FIXED_SCORE] [-o OUTFILE] fastq

Manipulate fastq quality records.

For each record in a fastq, shuffle the quality scores, effectively
randomizing per-base quality while holding mean read quality constant.

Alternately, set all bases to a fixed value.

positional arguments:
  fastq                 file to process

optional arguments:
  -h, --help            show this help message and exit
  -s SEED, --seed SEED  random seed (23)
  -f FIXED_SCORE, --fixed_score FIXED_SCORE
                        Instead of randomizing, used fixed score (phred33)
  -o OUTFILE, --outfile OUTFILE
                        Output file
```

Shuffle per-base quality score.
```bash
rand_fq_qual.py -o shuffle.reads.fastq reads.fastq
```

Shuffle per-base quality score and align reads.
```bash
rand_fq_qual.py --seed 23 reads.fastq | \
  minimap2 -t 24 -ax map-pb \
    -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k \
    -R "@RG\tSM:samplename" \
    reference.mmi - | \
  samtools sort -@ 24 - > shuffle.seed23.aln.bam
```

Change all per-base quality scores to fixed value and align reads.
```bash
rand_fq_qual.py --fixed 20 reads.fastq | \
  minimap2 -t 24 -ax map-pb \
    -O 5,56 -E 4,1 -B 5 -z 400,50 -r 2k \
    -R "@RG\tSM:samplename" \
    reference.mmi - | \
  samtools sort -@ 24 - > fixed.q20.aln.bam
```