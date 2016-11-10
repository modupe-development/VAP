#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR VariantdetectionAnalysis.pl

=pod

=head1 NAME

$0 -- Comprehensive pipeline : Variant detection using PICARD, GATK and produces and output folder.

=head1 SYNOPSIS

VariantdetectionAnalysis.pl -a configfile [--help] [--manual]

=head1 DESCRIPTION

Accepts all folders from frnakenstein output.
 
=head1 OPTIONS

=over 3

=item B<-a, -c, --config>=FILE

Configuration file (a template can be found @ .. ).  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   Getopt::Long
   Pod::Usage
use strict;
use File::Basename;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use DateTime;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage
=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR 
use strict;
use File::Basename;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use DateTime;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

 
#ARGUMENTS
my($help,$manual,$config,%CONFIGURE, $bamfile, $flag);
GetOptions (
                                "i|c|a|config=s"    =>      \$config,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument(s) are not found.\n", -exitval => 2, -verbose => 1)  if (! $config);
#file path for configuration file
open(CONFIG, "<", $config) or die "Configuration File \"$config\" can't be found\nTERMINATED!\n";
while (<CONFIG>){
  chomp;
  if ($_ =~ /\=/){
    my @ALL = split"\="; $ALL[1] =~ s/\[.+\]//g; #removing comments 
    foreach my $i (0..$#ALL) {$ALL[$i] =~ s/^\s+|\s+$//g;} #removing whitespace
    $CONFIGURE{$ALL[0]} = $ALL[1];
  }
}close CONFIG;
#INDEPENDENT PROGRAMS TO RUN
my $TOPHAT=$CONFIGURE{"TOPHAT"};
my $BOWTIE2=$CONFIGURE{"BOWTIE2"};
my $PICARDDIR=$CONFIGURE{"PICARDDIR"}."/picard.jar";
my $GATKDIR=$CONFIGURE{"GATKDIR"}."/GenomeAnalysisTK.jar";
my $VEP=$CONFIGURE{"VEP"} . "/scripts/variant_effect_predictor/variant_effect_predictor.pl";
my $SNPdat=$CONFIGURE{"SNPDAT"};
my $SNPEFF=$CONFIGURE{"SNPEFF"};
my $REF = $CONFIGURE{"REFERENCE"};
my $ANN = $CONFIGURE{"ANN"};
my $outputfolder = $CONFIGURE{"FOLDER"};
my $specie = $CONFIGURE{"ORGANISM"};
my $variantsname = "$outputfolder/all_variants.vcf";
my ($email, $notify); #notification options
#CLEAN UP OPTIONS
unless ($CONFIGURE{"TYPEOFDATA"} =~ /^(DNA|RNA)$/i) {
  die "TYPEOFDATA can only be DNA/RNA\t Sorry!\n";
}
unless ($CONFIGURE{"TYPEOFINPUT"} =~ /^(BAM|FASTQ)$/i) {
  die "TYPEOFINPUT can only be = BAM/FASTQ\t Sorry!\n";
}

#CREATE OUTPUT DIRECTORY
`mkdir $outputfolder`;

#EMAILS
my $subject = 'VAP-'.`date +%m-%d-%y_%T`; chomp $subject; my $notification = $outputfolder."/".$subject.'.log';
if (length($CONFIGURE{"SUBJECT"}) >= 1) { $subject = $CONFIGURE{"SUBJECT"}; }
if (length ($CONFIGURE{"EMAIL"}) >= 1) {
  $email = $CONFIGURE{"EMAIL"};
  $notify = "yes";
my $welcome = <<"ENDWELCOME";
Welcome to the VAP. 
You've subscribed for email update notifications and will arrive momentarily.
Currently you are running the VAP with the following options:
    1. Type of Data : $CONFIGURE{"TYPEOFDATA"}
    2. Type of Input : $CONFIGURE{"TYPEOFINPUT"}
    3. Reference file : $CONFIGURE{"REFERENCE"}
    4. Annotation file : $CONFIGURE{"ANN"}
    5. Output folder : $CONFIGURE{"FOLDER"}
    6. Filtering : $CONFIGURE{"FILTER"}
ENDWELCOME
  if ($CONFIGURE{"FILTER"} =~ /YES/i){
    $welcome .= "    7. Filtering options\n\t\tDP : $CONFIGURE{'DP'}\n";
  }
  my $create = ".welcome";
  open (WELCOME, ">$create");
  print WELCOME "$welcome\n";
  close WELCOME;
  system "mail -s \"VAP - $subject\" $email < $create";
  system "rm -rf $create"; 
}

#PROCESSING
#FILEPARSE NAME 
my $outgatk = fileparse($REF, qr/\.[^.]*(\.gz)?$/);
SORTGATK();
if ($CONFIGURE{"TYPEOFINPUT"} =~ /^BAM$/i) {
  if ($notify) { NOTIFICATION("Processing Variants from Bam File");  }
  $bamfile = $CONFIGURE{"FILENAME"};
  VARIANTS();
}
elsif ($CONFIGURE{"TYPEOFINPUT"} =~ /^FASTQ$/i) {
  if ($notify) { NOTIFICATION("Genome Assembly using TOPHAT");  }
  $bamfile = "$outputfolder/accepted_hits.bam"; 
  ASSEMBLY();
}
#NOT IN EFFECT
#ANNOTATION TOOL OPTIONS
#will be removed
my ($vepformat);
if ($CONFIGURE{"TOOL"} =~ /vep/i) {
  if ($CONFIGURE{"VEPFORMAT"}) {
    $vepformat= "--".lc($CONFIGURE{"VEPFORMAT"})."--everything on --terms ensembl";
  }
  else {
    $vepformat = "--database --terms ensembl"; 
  }
}
elsif ($CONFIGURE{"TOOL"} =~ /snpeff/i) {
  #working progress
}


##SUBROUTINES
sub ASSEMBLY {
  #TOPHAT assembly
  if ($notify) { NOTIFICATION("Assembly of the genome");  }
  #building index
  #`$BOWTIE2-build $REF $outputfolder/$outgatk`;
  if ($ANN){
    `$TOPHAT --library-type fr-unstranded --no-coverage-search -G $ANN -p 24 -o $outputfolder $outputfolder/$outgatk $CONFIGURE{"FILENAME"}`;
  }else {
    `$TOPHAT --library-type fr-unstranded --no-coverage-search -p 24 -o $outputfolder $outputfolder/$outgatk $CONFIGURE{"FILENAME"}`;
  }
  VARIANTS();
}

sub VARIANTS{
  my $DICT = $outputfolder."/".fileparse($REF, qr/(\..*)?$/).".dict";

  #CREATE DICTIONARY for gatk
  my $java = "java -jar $PICARDDIR CreateSequenceDictionary R=$REF O=$DICT"; `$java`;
  if ($notify) { NOTIFICATION("Sequence Dictionary complete");  }
  
  #QUALITY SCORE DISTRIBUTION
  $java = "java -jar $PICARDDIR QualityScoreDistribution INPUT=$bamfile OUTPUT=$outputfolder/qualityscores.txt CHART=$outputfolder/qualityscores.chart"; print $java;
  #CHECK QUALITY SCORE DISTRIBUTION
  open(CHECK,"<$outputfolder/qualityscores.txt");
  while (<CHECK>) { if (((split("\t",$_, 2))[0]) > 59){ $flag = "--fix_misencoded_quality_scores"; } }
  close CHECK;
  
  #SORT BAM
  if ($specie =~ /human/i || $specie =~ /homo/i){ #human samples can be reordered.
    `java -jar $PICARDDIR ReorderSam INPUT=$bamfile OUTPUT=$outputfolder/aln_sorted.bam REFERENCE=$REF`;
  } else {
    `java -jar $PICARDDIR SortSam INPUT=$bamfile OUTPUT=$outputfolder/aln_sorted.bam SO=coordinate`;
  }
  if ($notify) { NOTIFICATION("Sort Bam complete");  }      
  #ADDREADGROUPS
  my $addreadgroup = "java -jar $PICARDDIR AddOrReplaceReadGroups INPUT=$outputfolder/aln_sorted.bam OUTPUT=$outputfolder/aln_sorted_add.bam SO=coordinate RGID=Label RGLB=Label RGPL=illumina RGPU=Label RGSM=Label";
  `$addreadgroup`;
  if ($notify) { NOTIFICATION("Add read groups complete");  }
  
  #MARKDUPLICATES
  my $markduplicates = "java -jar $PICARDDIR MarkDuplicates INPUT=$outputfolder/aln_sorted_add.bam OUTPUT=$outputfolder/aln_sorted_mdup.bam M=$outputfolder/aln_sorted_mdup.metrics CREATE_INDEX=true";
  `$markduplicates`;
  if ($notify) { NOTIFICATION("Mark duplicates complete");  }
  
  if ($CONFIGURE{"TYPEOFDATA"} =~ /^RNA$/i) {
    #SPLIT&TRIM
    my $splittrim = "java -jar $GATKDIR -T SplitNCigarReads $flag -R $REF -I $outputfolder/aln_sorted_mdup.bam -o $outputfolder/aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar";
    `$splittrim`;
    if ($notify) { NOTIFICATION("SplitNCigars complete");  }
    
    #GATK
    my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_sorted_split.bam -o $variantsname";
    `$gatk`;
    if ($notify) { NOTIFICATION("Haplotype caller complete");  }
  }
  else {
    my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_sorted_mdup.bam -o $variantsname";
    `$gatk`;
    if ($notify) { NOTIFICATION("Haplotype caller complete");  }
  }
  #perl to select DP > 5 & get header information
  if ($CONFIGURE{"FILTER"} =~ /YES/){
    FILTERING($outputfolder, "$outputfolder/$variantsname");
    if ($notify) { NOTIFICATION("Fitering complete");  }
    $variantsname = "$outputfolder/all_filtered.vcf";
  }
  system "rm -rf $notification";
#  ANNOTATION($CONFIGURE{"TOOL"});
  if ($notify) { NOTIFICATION("Job Completed!!");  }
}

sub FILTERING {
  my $input = $_[1];
  my $wkdir = $_[0];
  unless(open(FILE,$input)){
    print "File \'$input\' doesn't exist\n";
    exit;
  }
  my $out = fileparse($input, qr/(\.vcf)?$/);
  my $output = "all_filtered.vcf";
  open(OUT,">$wkdir/$output");

  my @file = <FILE>; chomp @file; close (FILE);
  foreach my $chr (@file){
    unless ($chr =~ /^\#/){
      my @chrdetails = split('\t', $chr);
      my $chrIwant = $chrdetails[7];
      my @morechrsplit = split(';', $chrIwant);
      foreach my $Imptchr (@morechrsplit){
        if ($Imptchr =~ m/^DP/) {
          my @addchrsplit = split('=', $Imptchr);
          if ($addchrsplit[1] > $CONFIGURE{"DP"}){print OUT "$chr\n";}
        }
      }
    }
    else {
      print OUT "$chr\n";
    }
  }
  close (OUT);
}

#sub ANNOTATION {
#  if ($_[0] =~ /vep/i){
#    #ANNOTATIONS : running VEP
#    if (length $specie > 1 ){
#      my $veptxt = "perl $VEP -i $variantsname --fork 24 --species $specie $vepformat -o $outputfolder/all_VEP.txt";
#      `$veptxt`;
#      if ($notify) { NOTIFICATION("Vep \"TXT\" version complete");  }
#      my $vepvcf = "perl $VEP -i $variantsname --fork 24 --species $specie  --dir /home/modupeore17/.vep/ --cache --vcf --merged --everything on --terms ensembl --output_file $outputfolder/all_VEP.vcf";
#      `$vepvcf`;
#      if ($notify) { NOTIFICATION("Vep \"VCF\" version complete");  }
#    }
#  }
#  elsif ($_[0] =~ /snpeff/i){
#    #nothing yet
#  }
#}

sub NOTIFICATION {
  open (NOTE, ">>$notification");
  print NOTE "\nStatus is : $_[0]";
  close NOTE;
  system "mail -s \"$subject\" $email < $notification";
}

sub SORTGATK {
  $/ = "\>";
  my $zipgatk;
  if ($REF =~ /\.gz$/){$zipgatk=1;}
  # FILE HANDLES
  my ($GDATA,$GOUT1,$GOUT2,%GSEQ, %GSEQnum, %GSEQheader);

  # OPEN INPUT FILE(s)(in1 &/ in2)
  if($zipgatk) { $GDATA = IO::Uncompress::Gunzip->new( $REF ) or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; }
  else { open ($GDATA,$REF) or die $!;}

  #OPEN OUTPUT FILE(s)
  open($GOUT1, "> $outputfolder/$outgatk.fa") or die $!;
  open($GOUT2, "> $outputfolder/$outgatk.fa.fai") or die $!;
  my @gatkref = <$GDATA>;
  shift(@gatkref);
  foreach my $entry (@gatkref){
    my @pieces = split(/\n/, $entry);
    my $seq = '';
    foreach my $num (1.. $#pieces-1){
      $seq .= $pieces[$num];
    }
    my @allnuc = split('',$pieces[$#pieces]);
    if ($allnuc[$#allnuc] =~ /\>/){ $seq .= substr($pieces[$#pieces],0,-1); }
    else { $seq .= $pieces[$#pieces]; }
    $GSEQ{$pieces[0]} = $seq;
    $GSEQnum{$pieces[0]} = length($seq);
    $GSEQheader{$pieces[0]} = length($pieces[0]);
  }
  my ($check, $start, $newstart, $last); 
  unless ($specie =~ /human/i || $specie =~ /homo/i){
    foreach my $header (sort keys %GSEQ){
      if (length($header) >= 1) {
        print $GOUT1 ">$header\n$GSEQ{$header}\n";
        unless ($check){
          $start = $GSEQheader{$header}+2;
          $last = $GSEQnum{$header}+1;
          print $GOUT2 "$header\t$GSEQnum{$header}\t$start\t$GSEQnum{$header}\t$last\n";
          $check = "yes";
        }
        else {
          $newstart = $GSEQheader{$header}+2+$last+$start;
          $start = $newstart;
          $last = $GSEQnum{$header}+1;
          print $GOUT2 "$header\t$GSEQnum{$header}\t$start\t$GSEQnum{$header}\t$last\n";
          $check = "yes";
        }
      }
    }
  }
  else {
    my %TEMPLATE = ''; my $prefix; my $header; my @chr_others = ("X","Y","M");
    foreach my $pullheader (keys %GSEQ){
      if (length($pullheader) >= 1) {
        if ($pullheader =~ /^([A-Za-z]*)(\d+)$/) {
          $prefix = $1;
          $TEMPLATE{$2} = $pullheader;
        }
      }
    }
    foreach my $testheader (sort { $a <=> $b } keys %TEMPLATE){
      if (length($testheader) >= 1) {
        $header = $TEMPLATE{$testheader};
        print $GOUT1 ">$header\n$GSEQ{$header}\n";
        unless ($check){
          $start = $GSEQheader{$header}+2;
          $last = $GSEQnum{$header}+1;
          print $GOUT2 "$header\t$GSEQnum{$header}\t$start\t$GSEQnum{$header}\t$last\n";
          $check = "yes";
        }
        else {
          $newstart = $GSEQheader{$header}+2+$last+$start;
          $start = $newstart;
          $last = $GSEQnum{$header}+1;
          print $GOUT2 "$header\t$GSEQnum{$header}\t$start\t$GSEQnum{$header}\t$last\n";
          $check = "yes";
        }
      }
      delete $GSEQ{$header};
    }
    foreach (@chr_others){
      my $header = $prefix.$_;
      print $GOUT1 ">$header\n$GSEQ{$header}\n";
      unless ($check){
        $start = $GSEQheader{$header}+2;
        $last = $GSEQnum{$header}+1;
        print $GOUT2 "$header\t$GSEQnum{$header}\t$start\t$GSEQnum{$header}\t$last\n";
        $check = "yes";
      }
      else {
        $newstart = $GSEQheader{$header}+2+$last+$start;
        $start = $newstart;
        $last = $GSEQnum{$header}+1;
        print $GOUT2 "$header\t$GSEQnum{$header}\t$start\t$GSEQnum{$header}\t$last\n";
        $check = "yes";
      }
      delete $GSEQ{$header};
    }
  }
  close $GDATA; close $GOUT1; close $GOUT2;
  $/ = "\n";
  $REF = "$outputfolder/$outgatk.fa";
}
