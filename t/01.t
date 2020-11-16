#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>1;
use File::Basename qw/basename/;
use FindBin qw/$RealBin/;

$ENV{PATH}="$RealBin/../scripts:$ENV{PATH}";

subtest 'annotate' => sub{
  my $vcf = "$RealBin/data/SRR12894848.vcf.gz";
  my $gbk = "$RealBin/data/MN908947.3.gbk";
  my @out = `annapotater.pl $gbk $vcf`;
  if($?){
    BAIL_OUT("ERROR: could not run annapotater.pl: $!");
  }
  chomp(@out);

  my %expectedEff = (
    3037  => "SYNONYMOUS",
    14408 => "SYNONYMOUS",
    16647 => "NONSYNONYMOUS",
    23401 => "SYNONYMOUS",
    23403 => "NONSYNONYMOUS",
    29849 => "STOPLOST",
    29850 => "STOP",
  );

  plan tests => scalar(keys(%expectedEff)) + 1;

  my $vcfFormatLine = shift(@out);
  is($vcfFormatLine, '##fileformat=VCFv4.2', "VCF format line, first line");

  for my $line(@out){
    next if($line =~ /^#/);

    my @F = split(/\t/, $line);
    my $pos = $F[1];
    my @keyValue = split(/;/, $F[9]);
    next if(!$expectedEff{$pos});
    for (@keyValue){
      my($key, $value) = split(/=/);
      if($key eq 'EFF'){
        is($value, $expectedEff{$pos}, "Expected effect for $pos ($expectedEff{$pos})");
      }
    }
  }

};
