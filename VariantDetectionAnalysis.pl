#!/usr/bin/perl
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

our ($TOPHAT, $BOWTIE2, $PICARDDIR, $GATKDIR, $VEP, $SNPdat, $SNPEFF, $ANNOVAR, $REF, $ANN, $outputfolder, $specie, $variantsname);
my ($email, $notify); #notification options


my $date = `date +%m-%d-%y_%T`; chomp $date;
my $std_out = "VAP-$date.log";
my $std_err = "VAP-$date.err";

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
@ARGV == 0 || pod2usage("ERROR! Additional comments '@ARGV' not required\n");

configureinput(); #configure input options
#stdout and stderr
`mkdir -p $outputfolder`;
open(STDOUT, '>', "$outputfolder/$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$outputfolder/$std_err") or die "Error file doesn't exist";
 
#CREATE OUTPUT DIRECTORY
my @foldercontent = split("\n", `find $outputfolder`); #get details of the folder

#EMAILS
my $subject = "VAP-$date"; my $notification = $outputfolder."/".$subject.'.log';
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
    $welcome .= "    7. Filtering options\n";
    if ($CONFIGURE{"SEPARATEVARIANTS"} =~ /YES/i) {
      $welcome .= " \t#SNP\n\t\tDP : $CONFIGURE{'SNP-DP'}\n";
      $welcome .= "\t\tQD : $CONFIGURE{'SNP-QD'}\n\t\tMQ : $CONFIGURE{'SNP-MQ'}\n\t\tMQRankSum : $CONFIGURE{'SNP-MQRankSum'}\n\t\tFS : $CONFIGURE{'SNP-FS'}\n\t\tReadPosRankSum : $CONFIGURE{'SNP-ReadPosRankSum'}\n";
      $welcome .= " \t#INDEL\n\t\tDP : $CONFIGURE{'INDEL-DP'}\n";
      $welcome .= "\t\tQD : $CONFIGURE{'INDEL-QD'}\n\t\tFS : $CONFIGURE{'INDEL-FS'}\n\t\tReadPosRankSum : $CONFIGURE{'INDEL-ReadPosRankSum'}\n";
    } else {
      $welcome .= " \t\tDP : $CONFIGURE{'SNP-DP'}\n";
      $welcome .= "\t\tQD : $CONFIGURE{'SNP-QD'}\n\t\tMQ : $CONFIGURE{'SNP-MQ'}\n\t\tMQRankSum : $CONFIGURE{'SNP-MQRankSum'}\n\t\tFS : $CONFIGURE{'SNP-FS'}\n\t\tReadPosRankSum : $CONFIGURE{'SNP-ReadPosRankSum'}\n";
    }
  }
  my $create = ".welcome";
  open (WELCOME, ">$create");
  print WELCOME "$welcome\n";
  close WELCOME;
  system "mail -s \"VAP - $subject\" $email < $create";
  `mv $create $outputfolder/welcome-$date.log`;
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
sub configureinput {
  my ($prefix, $check);
  open(CONFIG, "<", $config) or die "Configuration File \"$config\" can't be found\nTERMINATED!\n";
  while (<CONFIG>){
    chomp;
    if (($_ =~ /\#(SNP)/) || ($_ =~ /\#(INDEL)/)){
      $prefix = $1;
      $check = 1;
    }  elsif ($_ =~ /^\#/) { $check = 2;}; #undo prefix addendum
    if ($_ =~ /\=/){
      my @ALL = split"\="; $ALL[1] =~ s/\[.+\]//g; #removing comments 
      foreach my $i (0..$#ALL) {$ALL[$i] =~ s/^\s+|\s+$//g;} #removing whitespace
      if ($check == 1) {
        my $newcheck = "$prefix-$ALL[0]";
        $CONFIGURE{$newcheck} = $ALL[1];
      } else {
        $CONFIGURE{$ALL[0]} = $ALL[1];
      }
    }
  } close CONFIG;
  #INDEPENDENT PROGRAMS TO RUN
  $TOPHAT=$CONFIGURE{"TOPHAT"};
  $BOWTIE2=$CONFIGURE{"BOWTIE2"};
  $PICARDDIR=$CONFIGURE{"PICARDDIR"}."/picard.jar";
  $GATKDIR=$CONFIGURE{"GATKDIR"}."/GenomeAnalysisTK.jar";
  $VEP=$CONFIGURE{"VEP"} . "/scripts/variant_effect_predictor/variant_effect_predictor.pl";
  $SNPdat=$CONFIGURE{"SNPDAT"};
  $SNPEFF=$CONFIGURE{"SNPEFF"};
  $ANNOVAR = $CONFIGURE{"ANNOVAR"} . "/table_annovar.pl";
  $REF = $CONFIGURE{"REFERENCE"};
  $ANN = $CONFIGURE{"ANN"};
  $outputfolder = $CONFIGURE{"FOLDER"};
  $specie = $CONFIGURE{"ORGANISM"};
  $variantsname = "$outputfolder/all_variants.vcf";
  #CLEAN UP OPTIONS
  unless ($CONFIGURE{"TYPEOFDATA"} =~ /^(DNA|RNA)$/i) {
    die "TYPEOFDATA can only be DNA/RNA\t Sorry!\n";
  }
  unless ($CONFIGURE{"TYPEOFINPUT"} =~ /^(BAM|FASTQ)$/i) {
    die "TYPEOFINPUT can only be = BAM/FASTQ\t Sorry!\n";
  }
}

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
  my $doesexist = (grep /\.dict/, @foldercontent)[0];
  unless ($doesexist) {
    `java -jar $PICARDDIR CreateSequenceDictionary R=$REF O=$DICT`;
    if ($notify) { NOTIFICATION("Sequence Dictionary complete");  }
  } else {
    if ($notify) { NOTIFICATION("Sequence Dictionary previously completed");  }
    print "\nVAPSTATUS:\tSequence Dictionary previously completed\n";
  }
  #QUALITY SCORE DISTRIBUTION
  $doesexist = (grep /qualityscores.txt/, @foldercontent)[0];
  unless ($doesexist) {
    `java -jar $PICARDDIR QualityScoreDistribution INPUT=$bamfile OUTPUT=$outputfolder/qualityscores.txt CHART=$outputfolder/qualityscores.chart`;
  }
  
  #CHECK QUALITY SCORE DISTRIBUTION
  open(CHECK,"<$outputfolder/qualityscores.txt");
  while (<CHECK>) { if (((split("\t",$_, 2))[0]) > 59){ $flag = "--fix_misencoded_quality_scores"; } }
  close CHECK;
  
  #SORT BAM
  $doesexist = (grep /aln_sorted.bam/, @foldercontent)[0];
  unless ($doesexist) {
   `java -jar $PICARDDIR SortSam INPUT=$bamfile OUTPUT=$outputfolder/aln_sorted.bam SO=coordinate`;
  #`java -jar $PICARDDIR ReorderSam INPUT=$bamfile OUTPUT=$outputfolder/aln_sorted.bam REFERENCE=$REF`;
    if ($notify) { NOTIFICATION("Sort Bam complete");  }
  } else {
    if ($notify) { NOTIFICATION("Sort Bam previously completed"); }
    print "\nVAPSTATUS:\tSort Bam previously completed\n"
  }
  
  #ADDREADGROUPS
  $doesexist = (grep /aln_sorted_add.bam/, @foldercontent)[0];
  unless ($doesexist) {
    my $addreadgroup = "java -jar $PICARDDIR AddOrReplaceReadGroups INPUT=$outputfolder/aln_sorted.bam OUTPUT=$outputfolder/aln_sorted_add.bam SO=coordinate RGID=Label RGLB=Label RGPL=illumina RGPU=Label RGSM=Label";
    `$addreadgroup`;
    if ($notify) { NOTIFICATION("Add read groups complete");  }
  } else {
    if ($notify) { NOTIFICATION("Add read groups previously completed");  }
    print "\nVAPSTATUS:\tAdd read groups previously completed\n";
  }
  
  #MARKDUPLICATES
  $doesexist = (grep /aln_sorted_mdup.bam/, @foldercontent)[0];
  unless ($doesexist) {
    my $markduplicates = "java -jar $PICARDDIR MarkDuplicates INPUT=$outputfolder/aln_sorted_add.bam OUTPUT=$outputfolder/aln_sorted_mdup.bam M=$outputfolder/aln_sorted_mdup.metrics CREATE_INDEX=true";
    `$markduplicates`;
    if ($notify) { NOTIFICATION("Mark duplicates complete");  }
  } else {
    if ($notify) { NOTIFICATION("Mark duplicates previously completed");  }
    print "\nVAPSTATUS:\tMark duplicates previously completed\n";
  }
  
  #REORDER SAM
  $doesexist = (grep /aln_resorted_mdup.bam/, @foldercontent)[0];
  unless ($doesexist) {
    `java -jar $PICARDDIR ReorderSam INPUT=$outputfolder/aln_sorted_mdup.bam OUTPUT=$outputfolder/aln_resorted_mdup.bam REFERENCE=$REF CREATE_INDEX=TRUE`;
    if ($notify) { NOTIFICATION("Resorted Mark duplicates complete");  }
  } else {
    if ($notify) { NOTIFICATION("Resorted Mark duplicates previously completed");  }
    print "\nVAPSTATUS:\tResorted Mark duplicates previously completed\n";
  }
  #specified into RNAseq or DNAseq
  if ($CONFIGURE{"TYPEOFDATA"} =~ /^RNA$/i) {
    #SPLIT&TRIM
    $doesexist = (grep /aln_sorted_split.bam/, @foldercontent)[0];
    unless ($doesexist) {
      my $splittrim = "java -jar $GATKDIR -T SplitNCigarReads $flag -R $REF -I $outputfolder/aln_sorted_mdup.bam -o $outputfolder/aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar";
      `$splittrim`;
      if ($notify) { NOTIFICATION("SplitNCigars complete");  }
    } else {
      if ($notify) { NOTIFICATION("SplitNCigars previously completed");  }
      print "\nVAPSTATUS:\tSplitNCigars previously completed\n";
    }
    
    #GATK
    $doesexist = (grep /all_variants.vcf/, @foldercontent)[0];
    unless ($doesexist) {
      my $gatk = "java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_sorted_split.bam -o $variantsname";
      `$gatk`;
      if ($notify) { NOTIFICATION("Haplotype caller complete");  }
    } else {
      if ($notify) { NOTIFICATION("Haplotype caller previously completed");  }
      print "\nVAPSTATUS:\tHaplotype caller previously completed\n";
    }
  }
  else {
    #working with DNA   
    #GATK
    $doesexist = (grep /all_variants.vcf/, @foldercontent)[0];
    unless ($doesexist) {
    `java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_resorted_mdup.bam -o $outputfolder/all_variants.vcf`;
    `java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_resorted_mdup.bam --emitRefConfidence GVCF -o $outputfolder/all_variants_emit.vcf`;
    `java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_resorted_mdup.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $outputfolder/all_variants_genotype.vcf`;
    `java -jar $GATKDIR -T HaplotypeCaller -R $REF -I $outputfolder/aln_resorted_mdup.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF -o $outputfolder/all_variants_genotypeemit.vcf`;
    if ($notify) { NOTIFICATION("Haplotype caller complete");  }
    } else {
      if ($notify) { NOTIFICATION("Haplotype caller previously completed");  }
      print "\nVAPSTATUS:\tHaplotype caller previously completed\n";
    }
  }
  
  #filter based on criteria
  if ($CONFIGURE{"FILTER"} =~ /YES/){
    my $filteryes = 0;
    my (@snpfilteryes, @indelfilteryes);
    my $snpfilterexpression = " --filterExpression \"";
    my $indelfilterexpression = " --filterExpression \"";
    #separate SNPS and INDELS
    if ($CONFIGURE{"SEPARATEVARIANTS"} =~ /YES/i) {
      `java -jar $GATKDIR -T SelectVariants -R $REF -V $outputfolder/all_variants.vcf -selectType SNP -o $outputfolder/raw_snps.vcf`;
      `java -jar $GATKDIR -T SelectVariants -R $REF -V $outputfolder/all_variants.vcf -selectType INDEL -o $outputfolder/raw_indels.vcf`;
      #filter based on defaults or specified criteria
      my $snpfilter = "java -jar $GATKDIR -T VariantFiltration -R $REF -V $outputfolder/raw_snps.vcf --filterName \"FAILED\" -o $outputfolder/filtered_snps.vcf";
      my $indelfilter = "java -jar $GATKDIR -T VariantFiltration -R $REF -V $outputfolder/raw_indels.vcf --filterName \"FAILED\" -o $outputfolder/filtered_indels.vcf";
      if ($CONFIGURE{"SNP-DP"} =~ /\d/){ $filteryes = 1;  push @snpfilteryes, "DP < $CONFIGURE{'SNP-DP'}"; }
      if ($CONFIGURE{"INDEL-DP"} =~ /\d/){ $filteryes = 1; push @indelfilteryes, "DP < $CONFIGURE{'INDEL-DP'}"; }
      if ($CONFIGURE{"SNP-QD"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "QD < $CONFIGURE{'SNP-QD'}"; }
      if ($CONFIGURE{"INDEL-QD"} =~ /\d/){ $filteryes = 1; push @indelfilteryes, "QD < $CONFIGURE{'INDEL-QD'}"; }
      if ($CONFIGURE{"SNP-MQ"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "MQ < $CONFIGURE{'SNP-MQ'}"; }
      if ($CONFIGURE{"SNP-MQRankSum"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "MQRankSum < $CONFIGURE{'SNP-MQRankSum'}"; }
      if ($CONFIGURE{"SNP-FS"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "FS < $CONFIGURE{'SNP-FS'}"; }
      if ($CONFIGURE{"INDEL-FS"} =~ /\d/){ $filteryes = 1; push @indelfilteryes, "FS < $CONFIGURE{'INDEL-FS'}"; }
      if ($CONFIGURE{"SNP-ReadPosRankSum"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "ReadPosRankSum < $CONFIGURE{'SNP-ReadPosRankSum'}"; }
      if ($CONFIGURE{"INDEL-ReadPosRankSum"} =~ /\d/){ $filteryes = 1; push @indelfilteryes, "ReadPosRankSum < $CONFIGURE{'INDEL-ReadPosRankSum'}"; }
      if ($filteryes == 1){ foreach ( @snpfilteryes ){ $snpfilterexpression .=  $_.' || ';} } $snpfilterexpression = substr($snpfilterexpression,0,-4); $snpfilterexpression .= "\""; $snpfilter .= $snpfilterexpression;
      if ($filteryes == 1){ foreach ( @indelfilteryes ){ $indelfilterexpression .=  $_.' || ';} } $indelfilterexpression = substr($indelfilterexpression,0,-4); $indelfilterexpression .= "\""; $indelfilter .= $indelfilterexpression;
      `$snpfilter`;
      `$indelfilter`;
      if ($notify) { NOTIFICATION("Filtering complete");  }
      FILTERING($outputfolder, "$outputfolder/filtered_snps.vcf"); #separating failed variants
      FILTERING($outputfolder, "$outputfolder/filtered_indels.vcf");
      if ($notify) { NOTIFICATION("Separating failed variants complete");  }
      #  `java -jar $GATKDIR -T VariantFiltration -R $REF -V $outputfolder/raw_snps.vcf \
      #    --filterExpression \"QD < $CONFIGURE{"QD"} || MQ < $CONFIGURE{"MQ"} || DP < $CONFIGURE{"DP"}  || MQRankSum < $CONFIGURE{"MQRankSum"} || FS > $CONFIGURE{"FS"} || ReadPosRankSum < $CONFIGURE{"ReadPosRankSum"}\" \
      #    --filterName "FAILED" \
      #    -o $outputfolder/filtered_snps.vcf`;
      #`java -jar $GATKDIR -T VariantFiltration -R $REF -V $outputfolder/raw_indels.vcf \
      #    --filterExpression \"QD < $CONFIGURE{"QD"} || DP < $CONFIGURE{"DP"} || FS >$CONFIGURE{"FS"} || ReadPosRankSum < $CONFIGURE{"ReadPosRankSum"}\" \
      #    --filterName "FAILED" \
      #    -o $outputfolder/filtered_indels.vcf`;
    } else {
      my $filter = "java -jar $GATKDIR -T VariantFiltration -R $REF -V $outputfolder/all_variants.vcf --filterName \"FAILED\" -o $outputfolder/all_filtered.vcf";
      if ($CONFIGURE{"SNP-DP"} =~ /\d/){ $filteryes = 1;  push @snpfilteryes, "DP < $CONFIGURE{'SNP-DP'}"; }
      if ($CONFIGURE{"SNP-QD"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "QD < $CONFIGURE{'SNP-QD'}"; }
      if ($CONFIGURE{"SNP-MQ"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "MQ < $CONFIGURE{'SNP-MQ'}"; }
      if ($CONFIGURE{"SNP-MQRankSum"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "MQRankSum < $CONFIGURE{'SNP-MQRankSum'}"; }
      if ($CONFIGURE{"SNP-FS"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "FS < $CONFIGURE{'SNP-FS'}"; }
      if ($CONFIGURE{"SNP-ReadPosRankSum"} =~ /\d/){ $filteryes = 1; push @snpfilteryes, "ReadPosRankSum < $CONFIGURE{'SNP-ReadPosRankSum'}"; }
      if ($filteryes == 1){ foreach ( @snpfilteryes ){ $snpfilterexpression .=  $_.' || ';} } $snpfilterexpression = substr($snpfilterexpression,0,-4); $snpfilterexpression .= "\""; $filter .= $snpfilterexpression;
      `$filter`;
      if ($notify) { NOTIFICATION("Filtering complete");  }
      FILTERING($outputfolder, "$outputfolder/all_filtered.vcf");
      if ($notify) { NOTIFICATION("Separating failed variants complete");  }
    }
    
    #FILTERING($outputfolder, "$outputfolder/$variantsname");
    #if ($notify) { NOTIFICATION("Separating failed variants complete");  }
    #$variantsname = "$outputfolder/all_filtered.vcf";
  }
  system "rm -rf $notification";
#  ANNOTATION($CONFIGURE{"TOOL"});
  if ($notify) { NOTIFICATION("Job Completed!!");  }
}

sub FILTERING {
  open(FILE,$_[1]) or die "File \'$_[1]\' doesn't exist\n";
  #my $out = fileparse($_[1], qr/(\.vcf)?$/);
  my $output = $_[0]."/".fileparse($_[1], qr/(\..*)?$/)."_passed.vcf";
  open(OUT,">$_[0]/$output");

  my @file = <FILE>; chomp @file; close (FILE);
  foreach my $chr (@file){
    unless ($chr =~ /^\#/){
      my @chrdetails = split('\t', $chr);
      my $chrIwant = $chrdetails[6];
      if ($chrIwant =~ /PASS/) {
        print OUT $chr, "\n";
      }
      #my @morechrsplit = split(';', $chrIwant);
      #foreach my $Imptchr (@morechrsplit){
      #  if ($Imptchr =~ m/^DP/) {
       #   my @addchrsplit = split('=', $Imptchr);
        #  if ($addchrsplit[1] > $CONFIGURE{"DP"}){print OUT "$chr\n";}
        #}
      #}
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
  `cat $notification >> $outputfolder/welcome-$date.log`;
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

