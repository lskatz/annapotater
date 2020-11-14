# Annapotater

A quick and dirty SNP annotator.
Only annotates SYNONYMOUS, NONSYNONYMOUS, or STOP.

Does not work on multisnp sites or indels.

# Installation

    cd annapotater
    perl Makefile.PL
    make
    make test

# Usage

    perl scripts/annapotater.pl t/data/MN908947.3.gbk t/out.vcf.gz.masked.vcf.gz

Give it a genbank file with CDS features and sequence,
and the next parameter is a gzipped VCF file.

# Output

Prints vcf to stdout.

Adds a notation onto the INFO field with the key `EFF`.
E.g.,
`EFF=NONSYNONYMOUS`

# Requirements

* perl and perl modules:
  * BioPerl
  * Array::IntSpan

# Etymology

Tomato, tomato  
Potato, potato  
Anna, anna  

