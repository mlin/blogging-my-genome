# Freebayes Variant Caller

## What does this app do?

This app calls variants (SNPs, indels, and other events) using FreeBayes v0.9.9.

## What are typical use cases for this app?

This app is used when you have mappings and need to identify variants (SNPs and indels).
Variation calling is performed with the Freebayes package.

## What data are required for this app to run?

This app requires coordinate-sorted mappings in BAM format (`*.bam`), and the associated reference
genome sequence in gzipped fasta format (`*.fasta.gz`, `*.fa.gz`).

Ideally, the mappings should contain associated read sample information (in technical language, each
read should have an RG tag pointing to a read group with an SM tag). Read groups must be unique across
samples. When multiple samples are present in the input, Freebayes will perform population variation
calling, and will distinguish among samples in the input file(s) based on the SM tag.

## What does this app output?

This app outputs a gzipped VCF file (`*.vcf.gz`) with the called variants. If multiple samples are
present in the input, the output will be a "multi-sample" VCF file.

## How does this app work?

This app runs 'freebayes'. For more information, consult the Freebayes manual at:

https://github.com/ekg/freebayes/blob/master/README
