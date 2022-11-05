# Annapotater

A quick and dirty SNP annotator.
Only annotates SYNONYMOUS, NONSYNONYMOUS, STOP, STOPLOST, or NA.
NA stands for "not assessed."

Does not work on multisnp sites or indels.

# Installation

Step 1: download the latest release and decompress it somewhere,
e.g., $HOME/bin/.../
and then `cd` into the annapotater folder.

    cd annapotater
    perl Makefile.PL
    cpanm --installdeps .
    make
    make test
    # Make install is optional, if you want to put the script in $HOME/bin
    make install

# Usage

    annapotater.pl t/data/MN908947.3.gbk t/data/SRR12894848.vcf.gz

Give it a genbank file with CDS features and sequence,
and the next parameter is a gzipped VCF file.
Use `--prot` to include protein translations.
Use `--all` to output all sites, useful as a passthrough.

# Output

Prints vcf to stdout.
Only prints sites that are annotated.

Adds a notation onto the INFO field with the key `EFF`.
E.g.,
`EFF=NONSYNONYMOUS`

If you give `--prot`, will also have a field `AA=`.

# Requirements

* perl and perl modules:
  * BioPerl
  * Array::IntSpan

# Etymology

Sounds like "annotator" but:  
Tomato, tomato  
Potato, potato  
Anna, anna  

And also, in some places, 'tater' is slang for 'potato.'

