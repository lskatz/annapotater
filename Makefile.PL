use strict;
use warnings;
 
use 5.010;
 
use ExtUtils::MakeMaker;

WriteMakefile(
  NAME         => 'Annopotater',
  VERSION_FROM => 'scripts/annopotater.pl',
  AUTHOR => q{Lee S. Katz <gzu2@cdc.gov>},
  EXE_FILES    => [
    'scripts/annopotater.pl',
  ],
  PREREQ_PM    => {
    # Core modules
    'Getopt::Long' => 0,
    'File::Basename'=>0,
    'Data::Dumper' =>0,

    'Array::IntSpan' => 0,
    'Bio::SeqIO' => 0,

  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
  MIN_PERL_VERSION => '5.16.1',
  LICENSE          => "gpl_3",
  dist=>{
    COMPRESS => "gzip -9f",
    SUFFIX   => "gz",
    PREOP    => $preop,
  },
  META_MERGE => {
    'meta-spec' => { version => 2 },
    resources => {
      repository => {
        type => 'git',
        url => 'git://github.com/lskatz/annopotater.git',
        web => 'https://github.com/lskatz/annopotater',
      },
    },
  },
);

