#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

use Array::IntSpan;
use Bio::SeqIO;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help)) or die $!;
  usage() if($$settings{help} || !@ARGV);

  my($gbk, @vcf) = @ARGV;

  my ($sequence, $cds) = readGbk($gbk, $settings);

  for my $vcf(@vcf){
    my $snps = readVcf($vcf, $settings);

    # Decide which SNPs need to get annotated better
    for my $chrom(sort keys(%$snps)){
      for my $pos(sort{$a<=>$b} keys(%{$$snps{$chrom}})){
        my $annotatedSnp = annotateSnp($$snps{$chrom}{$pos}, $sequence, $cds, $settings);
      }
    }

  }

  return 0;
}

sub readGbk{
  my($gbk, $settings)=@_;

  my(%seq, %ranges);
  my $numCds = 0;

  my $in=Bio::SeqIO->new(-file=>$gbk, -format=>"genbank");
  while(my $seq = $in->next_seq){
    # Cache these two variables for some speed and readability
    my $seqid    = $seq->id;
    my $sequence = $seq->seq;

    $seq{$seqid} = $sequence;
    $ranges{$seqid} = Array::IntSpan->new();

    for my $feat($seq->get_SeqFeatures){
      next if($feat->primary_tag ne 'CDS');

      my $start = $feat->location->start;
      my $end   = $feat->location->end;
      my $name  = "hypothetical";
      if($feat->has_tag('gene')){
        $name = (sort $feat->get_tag_values('gene'))[0];
      }

      $name .= "~~~$start~~~$end";
      $ranges{$seqid}->set_range($start, $end, $name);
      logmsg join("\t", $name, $start, $end) if($$settings{debug});
    }
  }

  return(\%seq, \%ranges);
}

sub readVcf{
  my($vcf, $settings) = @_;

  my %snp;

  my @header = qw(chrom pos id ref alt qual filter info format formatValues);

  open(my $fh, "zcat $vcf |") or die "ERROR: could not read $vcf: $!";
  while(<$fh>){
    next if(/^#/);

    chomp;

    # Get a hash of values with header=>value
    my %s = ();
    my @F = split(/\t/, $_);
    @s{@header} = @F;

    # Bioperl strips the version and so we will too
    $s{chrom} =~ s/\.\d+$//;

    if(scalar(@F) > scalar(@header)){
      logmsg "WARNING: there might be multiple samples in this vcf $vcf";
    }

    # Don't care if there is no alt
    if($s{alt} eq '.'){
      next;
    }
    # Don't care if the references _is_ the alt
    if($s{ref} eq $s{alt}){
      next;
    }

    $snp{$s{chrom}}{$s{pos}} = \%s;
  }
  close $fh;

  return \%snp;
}

# args:
#   $snp: the hash of snps with keys chrom pos id ref alt...
#   $sequence: hash of seqid=>sequence
#   $cds: hash of Array::IntSpan objects with CDS names
sub annotateSnp{
  my($snp, $sequence, $cds, $settings) = @_;

  # Get some basic info on the SNP
  my $chrom = $$snp{chrom};
  my $pos   = $$snp{pos};
  my $cdsInfo = $$cds{$chrom}->lookup($pos);
  my ($cdsName, $start, $end) = split(/~~~/, $cdsInfo);
  my $refSequence = $$sequence{$chrom};

  #$pos = 2; $start=2; $$snp{alt}="C"; $$snp{ref}=substr($refSequence,$pos-1,1); logmsg "DEBUG $$snp{ref}$pos$$snp{alt}";

  # The sequence of the CDS
  #my $refSequence = substr($$sequence{$chrom}, $start, $end - $start +1);

  my $altSequence = $refSequence;
  substr($altSequence, $pos - 1, 1) = $$snp{alt};

  # position within the CDS
  my $posWithinCds = $pos - $start + 1; # one-based int
  # 1, 2, or 3
  my $codonPos = $posWithinCds % 3 + 1;

  # Where the codon starts on the reference sequence
  my $codonStart = $posWithinCds - $codonPos + 1;
  my $refCodonSequence = substr($refSequence, $codonStart - 1 - $codonPos, 3);
  my $altCodonSequence = substr($altSequence, $codonStart - 1 - $codonPos, 3);

  my $refAA = translate($refCodonSequence);
  my $altAA = translate($altCodonSequence);
  #logmsg "substr(refSequence, $codonStart - 1, 3)";
  #logmsg "refCodonSequence $refCodonSequence $refAA";
  #logmsg "altCodonSequence $altCodonSequence $altAA";
  #die;

  logmsg substr($refSequence, 0, 60);
  logmsg "posWithinCds $posWithinCds = $pos - $start + 1";
  logmsg "  $$snp{ref}$pos$$snp{alt}";
  logmsg "  codonPos     $codonPos = $posWithinCds % 3 + 1";
  logmsg "  codonStart   $codonStart = $posWithinCds - $codonPos + 1";
  logmsg "    refCodonSequence $refCodonSequence $refAA";
  logmsg "    altCodonSequence $altCodonSequence $altAA";

  #die;
}

# Translate a DNA string
sub translate{
  my($sequence) = @_;
  my $seq = Bio::Seq->new(-id=>"tmp", -seq=>$sequence);
  my $prot= $seq->translate();
  return $prot->seq;
}

sub usage{
  print "$0: annotates a VCF to stdout
  Usage: $0 [options] annotations.gbk 1.vcf.gz [2.vcf.gz...]
  --help   This useful help menu
";
  exit 0;
}
