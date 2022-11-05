#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;

use Array::IntSpan;
use Bio::SeqIO;

our $VERSION = 0.5;

my @vcfHeader = qw(chrom pos id ref alt qual filter info format formatValues);

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(prot|protein all version help)) or die $!;
  if($$settings{version}){
    print "$0 v$VERSION\n";
    exit 0;
  }
  usage() if($$settings{help} || !@ARGV);

  my($gbk, @vcf) = @ARGV;

  my ($seqs, $cds) = readGbk($gbk, $settings);

  for my $vcf(@vcf){
    my ($vcfHeader, $snps) = readVcf($vcf, $settings);
    print $vcfHeader;

    # Decide which SNPs need to get annotated better
    for my $chrom(sort keys(%$snps)){
      for my $pos(sort{$a<=>$b} keys(%{$$snps{$chrom}})){
        print "$chrom\t$pos";
        my $thisSnp = $$snps{$chrom}{$pos};
        my $effect = snpEffect($thisSnp, $seqs, $cds, $settings);

        my $effString = "EFF=$$effect{effect}";
        if($$settings{prot}){
          $effString.= ";AA=$$effect{altAA}";
        }
        my $info = $$thisSnp{info};
        if($info eq '.'){
          $info  = $effString;
        }
        else {
          $info .= ";$effString";
        }

        $$thisSnp{info}=$info;

        for my $h(@vcfHeader){
          print "\t$$thisSnp{$h}";
        }
        print "\n";

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

    $seq{$seqid} = $seq;
    $ranges{$seqid} = Array::IntSpan->new();

    # Start the counter at -1 so that we can increment at
    # the beginning of the loop before the CDS filter.
    # Therefore any usage of $feature_count will let us use
    # zero-based integers.
    my $feature_count = -1;
    for my $feat($seq->get_SeqFeatures){
      $feature_count++;
      next if($feat->primary_tag ne 'CDS');
      #delete($$feat{_gsf_seq}); delete($$feat{_gsf_tag_hash}{translation}); die Dumper $feat;

      my $name  = "hypothetical";
      if($feat->has_tag('gene')){
        $name = (sort $feat->get_tag_values('gene'))[0];
      }

      # Look at each location for the CDS.
      # A CDS can be split and so it's best to loop through
      # each location.
      for my $location($feat->location->each_Location){
        my $start = $location->start;
        my $end   = $location->end;
        my $strand= $location->strand;
        my $loadedName = join("~~~",$name,$feature_count,$start,$end,$strand);
        # Mark the location of this CDS
        $ranges{$seqid}->set_range($start, $end, $loadedName);
      }

    }
  }

  return(\%seq, \%ranges);
}

sub readVcf{
  my($vcf, $settings) = @_;

  my ($headerStr, %snp);

  my @header = qw(chrom pos id ref alt qual filter info format formatValues);

  open(my $fh, "zcat $vcf |") or die "ERROR: could not read $vcf: $!";
  while(<$fh>){
    if(/^#/){
      $headerStr.=$_;
      next;
    }
    next if(/^\s*$/);

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
      if(!$$settings{all}){
        next;
      }
    }
    # Don't care if the references _is_ the alt
    if($s{ref} eq $s{alt}){
      if(!$$settings{all}){
        next;
      }
    }

    $snp{$s{chrom}}{$s{pos}} = \%s;
  }
  close $fh;

  return ($headerStr, \%snp);
}

# args:
#   $snp: the hash of snps with keys chrom pos id ref alt...
#   $sequence: hash of seqid=>sequence
#   $cds: hash of Array::IntSpan objects with CDS names
sub snpEffect{
  my($snp, $seqs, $cds, $settings) = @_;

  my %effect = (
    effect  => "NA", 
    altAA   => "",
  );

  # Get some basic info on the SNP
  my $chrom = $$snp{chrom};
  my $pos   = $$snp{pos};
  my $cdsInfo = $$cds{$chrom}->lookup($pos);
  if(!$cdsInfo){
    return \%effect;
  }

  # If this isn't really a variant and the alt==ref, then no effect
  if($$snp{alt} eq '.' || $$snp{alt} eq $$snp{ref}){
    return \%effect;
  }

  my ($cdsName, $feature_idx, $start, $end, $strand) = split(/~~~/, $cdsInfo);
  my $refSeq = $$seqs{$chrom};
  my $altSeq = $refSeq->clone;

  # mutate the altSeq with a substr() trick.
  # With substr(), positions are in base-0 and not in
  # base-1 like in bioperl.
  my $altSequence = $altSeq->seq;
  substr($altSequence, $pos - 1, 1) = $$snp{alt};
  $altSeq->seq($altSequence);

  # Reference protein sequence
  my @refFeat = $refSeq->get_SeqFeatures;
  my $refFeature = $refFeat[$feature_idx];
  my $refAA = $refFeature->seq->translate->seq;

  # Alt protein sequence
  my @altFeat = $altSeq->get_SeqFeatures;
  my $altFeature = $altFeat[$feature_idx];
  my $altAA = $altFeature->seq->translate->seq;
  $effect{altAA} = $altAA;

  my $snpName = "$$snp{ref}$pos$$snp{alt}";

  #logmsg "DEBUG"; $$settings{debug}=1;
  if($$settings{debug}){
    $refAA = substr($refAA, 600, 60);
    $altAA = substr($altAA, 600, 60);
    logmsg "$snpName";
    logmsg "  $refAA";
    logmsg "  $altAA";
  }

  # If nothing changed, then by definition, the mutation
  # was synonymous
  if($refAA eq $altAA){
    $effect{effect} = "SYNONYMOUS";
    return \%effect;
  }
  
  # If there are more stops than the ref, then there was
  # a stop codon
  my $numRefStops = () = $refAA =~ /\*/g;
  my $numAltStops = () = $altAA =~ /\*/g;
  if($numRefStops < $numAltStops){
    $effect{effect} = "STOP";
    return \%effect;
  }
  if($numRefStops > $numAltStops){
    $effect{effect} = "STOPLOST";
    return \%effect;
  }

  $effect{effect} = "NONSYNONYMOUS";
  return \%effect;
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
  --prot   Add the protein sequence translation for any variant features
  --all    Output all VCF lines regardless of there being an effect. Warning: slow
  --help   This useful help menu
";
  exit 0;
}
