#!/usr/bin/perl
# CODE FOR 
use strict;
use File::Basename;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use threads;
use Thread::Queue;

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#tools
our ($TOPHAT, $BOWTIE, $PICARD, $GATK, $HISAT, $FASTQC, $STAR, $SAMTOOLS, $BCFTOOLS, $BWA);

#files
our ($value, %CONFIGURE, %FILE, @content, @arrayfile, @threads);

#inputs
our ($REF, $SAMREF, $GATKREF, $ANN, $ANNGTF, $outputfolder, $THREADS);
our ($sample, $reads, $doesexist, @eachread, $codes, $outted, $newout);
our ($email, $notify, $syntax); #notification options
our ($date, $pwd, $std_out, $std_err, $help, $manual, $config, $bamfile, $flag);

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - C O N F I G U R A T I O N - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

initialize_workflow(); 
configureinput(); #configure input options
INPUTFILES(); #sort inputfiles

#standard error and output file
`mkdir -p $outputfolder`;
chdir ("$outputfolder"); #working from the outputfolder;
open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";

#Notification configuration
my $subject = "VAP-$date"; my $notification = $outputfolder."/".$subject.'.log';
unless ((length($CONFIGURE{"SUBJECT"}) <= 1) || ($CONFIGURE{"SUBJECT"} =~ /false/)) { $subject = $CONFIGURE{"SUBJECT"}; }
unless ((length ($CONFIGURE{"EMAIL"}) <= 1) || ($CONFIGURE{"EMAIL"} =~ /false/)) {
  $email = $CONFIGURE{"EMAIL"};
  $notify = "yes";
my $welcome = <<"ENDWELCOME";
Welcome to the VAP. 
You've subscribed for email update notifications and will arrive momentarily.
  Good luck!
  VAP.
ENDWELCOME

  my $create = ".welcome";
  open (WELCOME, ">$create");
  print WELCOME "$welcome\n";
  close WELCOME;
  system "mail -s 'VAP - $subject' $email < $create";
  `mv $create $outputfolder/welcome-$date.log`;
  system "rm -rf $create"; 
} # end if email notification is selected
`mkdir -p "$outputfolder/tmp"`;
my $parsedemail = $email; $parsedemail =~ s/\@/\\\@/g;
#notification subroutine to all the different files
my $endnotify= <<"ENDNOTIFY";

sub NOTIFICATION {
  \$email = '$parsedemail';
  open (NOTE, ">>\$0.notice.log");
ENDNOTIFY
$endnotify .= <<'ENDNOTIFY';
  print NOTE "\nStatus is : $_[0]";
  close NOTE;
ENDNOTIFY
$endnotify .= <<"ENDNOTIFY";
  system "mail -s '$subject' \$email < \$0.notice.log";
}
ENDNOTIFY

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - W O R K F L O W - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#PROCESSING
if ( (exists $CONFIGURE{"FASTQ"}) && ($CONFIGURE{"FASTQ"} ne "false") ) {
  my $testindex = `ls $CONFIGURE{"FASTQ"} | head -n 1`; chomp $testindex; 
  if (-e $testindex) {
    foreach (keys %FILE){
      if ($CONFIGURE{"RUNFASTQC"} eq "true") { 
        $newout = "$outputfolder/tmp/$_-fastqc.txt";
        push @arrayfile, $newout; open (OUT, ">$newout");
        print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
        FASTQC($_, $FILE{$_}, $newout);
      } #fastqc
      if ($CONFIGURE{"SAMPLERNA"} eq "true") { #rna
        if ($CONFIGURE{"RUNTOPHAT"} eq "true"){
          $newout = "$outputfolder/tmp/$_-tophat.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          TOPHAT($_, $FILE{$_}, $newout);
        } #tophat
        if ($CONFIGURE{"RUNHISAT"} eq "true"){
          $newout = "$outputfolder/tmp/$_-hisat.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          HISAT($_, $FILE{$_}, $newout);
        } #hisat
        if ($CONFIGURE{"RUNSTAR"} eq "true"){
          $newout = "$outputfolder/tmp/$_-star.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          STAR($_, $FILE{$_}, $newout);
        } #star
      } elsif ($CONFIGURE{"SAMPLEDNA"} eq "true") { #dna
        if ($CONFIGURE{"RUNBOWTIE"} eq "true"){
          $newout = "$outputfolder/tmp/$_-bowtie.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          BOWTIE($_, $FILE{$_}, $newout);
        } #bowtie
        if ($CONFIGURE{"RUNBWA"} eq "true"){
          $newout = "$outputfolder/tmp/$_-bwa.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          BWA($_, $FILE{$_}, $newout);
        } #bwa
      } #end sequence type
    } #end file
  } #end if file found
} #end if fastq file
elsif ( (exists $CONFIGURE{"BAM"}) || (exists $CONFIGURE{"SAM"}) ) { #if bam or sam file
  my $nac = "RNA"; #default is RNA
  if ($CONFIGURE{"SAMPLEDNA"} eq "true") { $nac = "DNA"; } #if dna
  foreach (keys %FILE){ 
    if ($CONFIGURE{"RUNFASTQC"} eq "true"){
      $newout = "$outputfolder/tmp/$_-fastqc.txt";
      push @arrayfile, $newout;
      open (OUT, ">$newout"); print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
      FASTQC($_, $FILE{$_}, $newout);
    } #fastqc

    if ($CONFIGURE{"RUNVAP"} eq "true"){ #run variant analysis pipeline
      if ((exists $CONFIGURE{"SAM"}) && ($CONFIGURE{"SAM"} ne "false")) { #if sam file
        my $testindex = `ls $CONFIGURE{"SAM"} | head -n 1`; chomp $testindex;
        if (-e $testindex) {
          $newout = "$outputfolder/tmp/$_-sam.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          SAMBAM($_, $FILE{$_}, $newout, $nac);
        } 
      } elsif ((exists $CONFIGURE{"BAM"}) && ($CONFIGURE{"BAM"} ne "false")) { #if bam file  
       my $testindex = `ls $CONFIGURE{"BAM"} | head -n 1`; chomp $testindex;
        if (-e $testindex) {
          $newout = "$outputfolder/tmp/$_-bam.txt";
          push @arrayfile, $newout;
          open (OUT, ">$newout");
          print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
          VAP($_, $FILE{$_}, $newout,$nac, $_);
        } #end if file exist
      } #end if sam/bam 
    } #end if run vap
  } #end file
} #end if bam or sam file

my $queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..5; #execute 10 threads
$builder->join; #join threads
foreach (@threads){$_->join;}

#end of job
`echo "Job Completed" >>  $outputfolder/welcome-$date.log` if $notify;
system "mail -s \"VAP - $subject\" $email < $outputfolder/welcome-$date.log" if $notify;


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - -  S U B R O U T I N E S  - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


sub initialize_workflow {
  #set standard output and error files 
  $date = `date +%m-%d-%y_%T`; chomp $date;
  $std_out = "VAP-$date.log"; $std_err = "VAP-$date.err";
  $pwd = `pwd`; chomp $pwd;
  #ARGUMENTS
  GetOptions ( "i|c|a|config=s"=>\$config, "h|help"=>\$help, "man|manual"=>\$manual );

  # VALIDATE ARGS
  pod2usage( -verbose => 2 )  if ($manual);
  pod2usage( -verbose => 1 )  if ($help);
  pod2usage( -msg  => "ERROR!  Required argument '-c <config_file>' not found.\n", -exitval => 2, -verbose => 1)  if (! $config);
  @ARGV == 0 || pod2usage("ERROR! Additional comments '@ARGV' not required\n");
} #end of subroutine: initialize_workflow


sub main { #beginning thread
  foreach my $count (0..$#arrayfile) {
    while(1) {
      if ($queue->pending() <100) {
        $queue->enqueue($arrayfile[$count]);
        last;
      }
    }
  }
  foreach(1..5) { $queue-> enqueue(undef); }
} #end of subroutine: main


sub processor { #queueing the thread
  my $query;
  while ($query = $queue->dequeue()){
    parseinput($query);
  }
} #end of subroutine: processor


sub parseinput{ #working with each thread
  print "Submitted\t$_[0]\t",`date +%x-%X`;
  system("perl $_[0] 1>$_[0].log 2>$_[0].err");
  system("echo 'from $_[0].log' >> $outputfolder/$std_out");
  system("echo 'from $_[0].err' >> $outputfolder/$std_err");
  system("cat $_[0].log >> $outputfolder/$std_out");
  system("cat $_[0].err >> $outputfolder/$std_err");
  system("cat $_[0].notice.log >> $outputfolder/welcome-$date.log");
  print "Completed\t$_[0]\t",`date +%x-%X`;
} #end of subroutine: parseinput


sub configureinput {
  #read the config file
  open(CONFIG, "<", $config) or die "Configuration File \"$config\" can't be found\nTERMINATED!\n";
  while (<CONFIG>){
    chomp;
    unless ($_ =~ /^#/){
      if ($_ =~ /\=/){
        $_ =~ tr/"//d; $_ =~ tr/'//d;
        $_ =~ /(\S+)\s*=\s*(\S+)/;
        my ($tool, $info) = (uc($1), $2);
        $CONFIGURE{$tool} = $info;
      }
    }
  } close CONFIG;

  #INDEPENDENT PROGRAMS TO RUN
  $TOPHAT=$CONFIGURE{"TOPHAT"};
  $BOWTIE=$CONFIGURE{"BOWTIE"};
  $PICARD=$CONFIGURE{"PICARD"};
  $FASTQC=$CONFIGURE{"FASTQC"};
  $SAMTOOLS=$CONFIGURE{"SAMTOOLS"};
  $BCFTOOLS=$CONFIGURE{"BCFTOOLS"};
  $HISAT=$CONFIGURE{"HISAT"};
  $BWA=$CONFIGURE{"BWA"};
  $STAR=$CONFIGURE{"STAR"};
  $GATK=$CONFIGURE{"GATK"};
  $REF=$CONFIGURE{"GENOME"}; 
  $ANN=$CONFIGURE{"GFF"};
  $ANNGTF=$CONFIGURE{"GTF"};
  if ($CONFIGURE{"THREADS"} > 1 ){ $THREADS = $CONFIGURE{"THREADS"};} else {$THREADS = 1;}
  if ($CONFIGURE{"OUTPUTDIR"}) { $outputfolder=$CONFIGURE{"OUTPUTDIR"};} else {$outputfolder="$pwd/VAPOUT-$date";}
} #end of subroutine: configureinput


sub INPUTFILES { #sorting the input files
  #FASTQ files and get basename if exist
  if ( (exists $CONFIGURE{"FASTQ"}) && ($CONFIGURE{"FASTQ"} ne "false") ) {
    @content = sort {$a <=> $b || $a cmp $b} (split("\n", `ls $CONFIGURE{"FASTQ"}`)); #get details of the folder
    foreach (@content) {
      my $file = basename($_);
      if($file =~ /.+[_\.][Rr][12].+/) {
        $file =~ /(.+)_.*[_\.][Rr][12].+/; $value = $1;
      }
      elsif($file =~ /.+[_\.][12].+/) {
        $file =~ /(.+)_.*[_\.][12].+/; $value = $1;
      }
      elsif($file =~ /.+[_\.]pe[12]/) {
        $file =~ /(.+)_.*[_\.]pe[12].+/; $value = $1;
      }
      elsif($file =~ /.+[_\.]PE[12]/) {
        $file =~ /(.+)_.*[_\.]PE[12].+/; $value = $1;
      }
      else { $value = $1; }
      
      if (exists $FILE{$value}){ $FILE{$value} = "$FILE{$value} $_"; }
      else { $FILE{$value} = $_; }
    }
  } #end if fastq is specified
  #SAM/BAM files and get basename if exist
  elsif ( (exists $CONFIGURE{"BAM"}) && ($CONFIGURE{"BAM"} ne "false") ) {
    @content = sort {$a <=> $b || $a cmp $b} (split("\n", `ls $CONFIGURE{"BAM"}`)); #get details of the folder
    foreach (@content) {
      my $file = basename($_);
      if($file =~ /(.+)[\.]sam/) {
         $value = $1;
      } elsif($file =~ /(.+)[\.]bam/) {
        $value = $1;
      }
      $FILE{$value} = $_;
    }
  } # end if bam specified
  elsif ( (exists $CONFIGURE{"SAM"}) && ($CONFIGURE{"SAM"} eq "false") ) {
    @content = sort {$a <=> $b || $a cmp $b} (split("\n", `ls $CONFIGURE{"SAM"}`)); #get details of the folder
    foreach (@content) {
      my $file = basename($_);
      if($file =~ /(.+)[\.]sam/) {
        $value = $1;
      } elsif($file =~ /(.+)[\.]bam/) {
        $value = $1;
      }
      $FILE{$value} = $_;
    }
  } # end if sam specified
} #end of subroutine: INPUTFILES


sub FASTQC { #run fastqc
  ($sample, $reads, $outted) = @_;
  chdir("$outputfolder"); 
  `mkdir -p $sample/fastqc`;
  chdir ("$sample/fastqc"); 

#printing instructions in fastqc job file
$codes = <<"ENDCODES";
#FastQC $sample
chdir("$outputfolder");
`mkdir -p $sample/fastqc`;
chdir ("$sample/fastqc");
ENDCODES
  open (OUT, ">>$outted");
  print OUT "$codes\n";

  #parsing name of eventual fastqc file
  @eachread = split(/\s/, $reads);
  foreach my $single (@eachread) { 
    my $singlefolder = fileparse($single);
    if ($singlefolder =~ /fastq\.gz/) { $singlefolder =~ s/\.fastq\.gz/_fastqc/g; }
    elsif ($singlefolder =~ /fastq$/) { $singlefolder =~ s/\.fastq/_fastqc/g; }
    elsif ($singlefolder =~ /sam$/) { $singlefolder =~ s/\.sam/_fastqc/g; }
    elsif ($singlefolder =~ /bam$/) { $singlefolder =~ s/\.bam/_fastqc/g; }
    elsif ($singlefolder =~ /fq\.gz/) { $singlefolder =~ s/\.fq\.gz/_fastqc/g; }
    elsif ($singlefolder =~ /fq$/) { $singlefolder =~ s/\.fq/_fastqc/g; }
    else {die "Can't find fastq file type\n";}
    $doesexist = (grep /$singlefolder/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      $singlefolder = fileparse($single);

#printing instructions in fastqc job file
$codes = <<"ENDCODES";
`cp $reads ./`;
`$FASTQC $singlefolder`;
ENDCODES

      if ($singlefolder =~ /fastq\.gz/) { $singlefolder =~ s/\.fastq\.gz/_fastqc/g; }
      elsif ($singlefolder =~ /fastq$/) { $singlefolder =~ s/\.fastq/_fastqc/g; }
      elsif ($singlefolder =~ /sam$/) { $singlefolder =~ s/\.sam/_fastqc/g; }
      elsif ($singlefolder =~ /bam$/) { $singlefolder =~ s/\.bam/_fastqc/g; }
      elsif ($singlefolder =~ /fq\.gz/) { $singlefolder =~ s/\.fq\.gz/_fastqc/g; }
      elsif ($singlefolder =~ /fq$/) { $singlefolder =~ s/\.fq/_fastqc/g; }
      else {die "Can't find fastq file type\n";}

      print OUT "$codes\n";
$codes = <<"ENDCODES";
`unzip $singlefolder.zip`;
`cp $singlefolder/summary.txt $singlefolder.txt`;
`rm -rf $singlefolder`;
ENDCODES

      print OUT "$codes\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - FastQC complete\");\n"; }
    } # end unless ($doesexist)
  } # end foreach (@eachread)
  close OUT;
} #end of subroutine: FASTQC


sub TOPHAT { #runtophat
  ($sample, $reads, $outted) = @_;
  my $testindex = $CONFIGURE{'GENOMEDIR'}.'/*.rev.1.bt2'; $testindex = `ls $testindex`; chomp $testindex;
  unless (-e $CONFIGURE{'GENOMEDIR'}) { die "GENOMEDIR FOLDER specified doesn't exist\n"; }
  unless (-e $testindex) { die "Genome Indexes for TOPHAT2 dont exist in GENOMEDIR specified\n$testindex\n"; }
  my $newgenomeindex = `ls -1 $testindex`; $newgenomeindex = (split('.rev.1.bt2',$newgenomeindex))[0];
  unless (-e $TOPHAT) { print "ERROR: $TOPHAT path specified doesn't exist\n"; exit; }
  chdir("$outputfolder");
  `mkdir -p $sample/tophat`;
  chdir ("$sample/tophat"); 

$codes = <<"ENDCODES";
#TopHAT $sample
chdir("$outputfolder");
`mkdir -p $sample/tophat`;
chdir ("$sample/tophat");
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";

  $doesexist = (grep /tophat\.bam/, (split("\n", `find ./`)))[0]; #checking if file already exists
  unless ($doesexist) {  
$codes = <<"ENDCODES";
`$TOPHAT -p $THREADS --library-type fr-firststrand --no-coverage-search -G $ANN -o ./ $newgenomeindex $reads`;
`mv accepted_hits.bam $sample.tophat.bam`;
ENDCODES

    print OUT "$codes\n";

    if ($notify) { print OUT "NOTIFICATION(\"$sample - TOPHAT complete\");\n";  }
  }
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$outputfolder/$sample/tophat/$sample.tophat.bam", $outted, "RNA", "$sample/tophat");
  }
  
} #end of subroutine: TOPHAT


sub STAR { #runstar
  ($sample, $reads, $outted) = @_;
  unless (-e $REF) { die "GENOME specified doesn't exist\n"; }
  unless (-e $STAR) { die "STAR tool path is incorrect, make sure the correct path is specified\n"; }
  chdir("$outputfolder");
  `mkdir -p $sample/star`;
  chdir ("$sample/star");

$codes = <<"ENDCODES";
#STAR $sample
chdir("$outputfolder");
`mkdir -p $sample/star`;
chdir ("$sample/star");
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";

  $doesexist = (grep /star\.bam/, (split("\n", `find ./`)))[0];
  unless ($doesexist) {
    @eachread = split(/\s/, $reads);
    my $addGTF=""; #making sure the GTF file is specified
    if ($CONFIGURE{'GTF'}){$addGTF="--sjdbGTFfile $ANNGTF";} 
    my $readlength=`zcat $eachread[0] | head -n 2 | tail -n 1 | awk '{print length}'`; $readlength = $readlength - 1;
    my $stargen = "mkdir -p STARref; $STAR --runThreadN $THREADS --runMode genomeGenerate --genomeDir STARref --genomeFastaFiles $REF";
    my $starstat = "mkdir -p 1PASS; cd 1PASS; $STAR --runThreadN $THREADS --genomeDir $outputfolder/$sample/star/STARref --readFilesCommand zcat --sjdbOverhang $readlength --outFileNamePrefix $sample. $addGTF --readFilesIn ";
    my $star2ref = "mkdir -p STARref2; $STAR  --runThreadN $THREADS --runMode genomeGenerate --genomeDir STARref2 --genomeFastaFiles $REF --sjdbFileChrStartEnd $outputfolder/$sample/star/1PASS/SJ.out.tab";
    my $star2pass = "mkdir -p 2PASS; cd 2PASS; $STAR --runThreadN $THREADS --genomeDir $outputfolder/$sample/star/STARref2 --readFilesCommand zcat --sjdbOverhang $readlength --outFileNamePrefix $sample.2. $addGTF --readFilesIn ";
    foreach my $single (@eachread) {
      $star2pass .= "$single ";
      $starstat .= "$single ";
    }

$codes = <<"ENDCODES";
#`cp $reads ./`;
#`gunzip *gz`;
`$stargen`;
`$starstat; cd ..`;
`$star2ref`;
`$star2pass; cd ..`;
`cp 1PASS/$sample.Aligned.out.sam ./`;
`cp 2PASS/$sample.2.Aligned.out.sam ./`;
`$SAMTOOLS view -bS $sample.2.Aligned.out.sam -o $sample.2.star.bam`;
`$SAMTOOLS view -bS $sample.Aligned.out.sam -o $sample.star.bam`;
ENDCODES

    print OUT "$codes\n";

    if ($notify) { print OUT "NOTIFICATION(\"$sample - STAR complete\");\n";  }
  } # end unless(doesexist)
  close OUT;
  if ($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$sample.2.star.bam", $outted, "RNA","$sample/star2");
    VAP($_, "$sample.star.bam", $outted, "RNA","$sample/star");
  } # end if (parse to VAP)
} #end of subroutine: STAR


sub HISAT { #runhisat
  ($sample, $reads, $outted) = @_;
  my $testindex = $CONFIGURE{'GENOMEDIR'}.'/*.1.ht2'; $testindex = `ls $testindex`; chomp $testindex;
  unless (-e $CONFIGURE{'GENOMEDIR'}) { die "GENOMEDIR FOLDER specified doesn't exist\n"; }
  unless (-e $testindex) { die "Genome Indexes for HiSAT2 dont exist in GENOMEDIR specified\n"; }
  my $newgenomeindex = `ls -1 $testindex`; $newgenomeindex = (split('.1.ht2',$newgenomeindex))[0];  
  unless (-e $HISAT) { die "HISAT tool path is incorrect, make sure the correct path is specified\n"; }

  chdir("$outputfolder");
  `mkdir -p $sample/hisat`;
  chdir ("$sample/hisat");

$codes = <<"ENDCODES";
#HiSAT $sample
chdir("$outputfolder");
`mkdir -p $sample/hisat`;
chdir ("$sample/hisat");
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";
 
  
  $doesexist = (grep /hisat\.bam/, (split("\n", `find ./`)))[0];
  my $hisatstat;
  unless ($doesexist) {
    @eachread = split(/\s/, $reads);
    if ($#eachread >= 1){ 
      $hisatstat = "$HISAT -p $THREADS -x $newgenomeindex -S $sample.hisat.sam -1 $eachread[0] -2 $eachread[1] 2> $sample"."_align.txt";
    } else {
      $hisatstat = "$HISAT -p $THREADS -x $newgenomeindex -S $sample.hisat.sam -U $eachread[0] 2> $sample"."_align.txt";
    }

$codes = <<"ENDCODES";
`$hisatstat`;
`$SAMTOOLS view -bS $sample.hisat.sam -o $sample.hisat.bam`;
ENDCODES

    print OUT "$codes\n";
    if ($notify) { print OUT "NOTIFICATION(\"$sample - HiSAT complete\");\n";  }
  } # end unless(doesexist)
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$outputfolder/$sample/hisat/$sample.hisat.bam", $outted, "RNA", "$sample/hisat");
  } # end if (parse to VAP)
} #end of subroutine: HISAT


sub BWA { #run bwa mem
  ($sample, $reads, $outted) = @_;
  unless (-e $REF) { die "GENOME specified doesn't exist\n"; }
  my $testindex = $REF.'.pac';
  unless (-e $testindex) { die "Genome Indexes for BWA dont exist in GENOME folder or with GENOME prefix\n"; }
  unless (-e $BWA) { die "BWA tool path is incorrect, make sure the correct path is specified\n"; }
  chdir("$outputfolder");
  `mkdir -p $sample/bwa`;
  chdir ("$sample/bwa");

$codes = <<"ENDCODES";
#BWA $sample
chdir("$outputfolder");
`mkdir -p $sample/bwa`;
chdir ("$sample/bwa");
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";

  $doesexist = (grep /bwa\.bam/, (split("\n", `find ./`)))[0];
  unless ($doesexist) {

$codes = <<"ENDCODES";    
`$BWA mem -t $THREADS $REF $reads > $sample.bwa.sam`;
`$SAMTOOLS view -bS $sample.bwa.sam -o $sample.bwa.bam`;
ENDCODES

    print OUT "$codes\n";
    if ($notify) { print OUT "NOTIFICATION(\"$sample -BWA complete\");\n";  }
  } # end unless(doesexist)
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$sample.bwa.bam", $outted, "DNA", "$sample/bwa");
  } #end if (parse to VAP)
} #end of subroutine: BWA


sub BOWTIE { #run bowtie
  ($sample, $reads, $outted) = @_;
  my $testindex = $CONFIGURE{'GENOMEDIR'}.'/*.rev.1.bt2'; $testindex = `ls $testindex`; chomp $testindex;
  unless (-e $CONFIGURE{'GENOMEDIR'}) { die "GENOMEDIR FOLDER specified doesn't exist\n"; }
  unless (-e $testindex) { die "Genome Indexes for BOWTIE2 dont exist in GENOMEDIR specified\n"; }
  my $newgenomeindex = `ls -1 $testindex`; $newgenomeindex = (split('.rev.1.bt2',$newgenomeindex))[0];
  unless (-e $BOWTIE) { die "BOWTIE2 tool path is incorrect, make sure the correct path is specified\n"; }

  chdir("$outputfolder");
  `mkdir -p $sample/bowtie`;
  chdir ("$sample/bowtie");

$codes = <<"ENDCODES";
#BOWTIE $sample
chdir("$outputfolder");
`mkdir -p $sample/bowtie`;
chdir ("$sample/bowtie");
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";

  $doesexist = (grep /bowtie\.bam/, (split("\n", `find ./`)))[0];
  unless ($doesexist) {
    @eachread = split(/\s/, $reads);
    my $bowtiestat = "$BOWTIE -p $THREADS -x $newgenomeindex -S $sample.bowtie.sam -1 $eachread[0] -2 $eachread[1]";

$codes = <<"ENDCODES";
`$bowtiestat`;
`$SAMTOOLS view -bS $sample.bowtie.sam -o $sample.bowtie.bam`;
ENDCODES

    print OUT "$codes\n";
    if ($notify) { print OUT "NOTIFICATION(\"$sample -BOWTIE complete\");\n";  }
  } # end unless(doesexist)
  close OUT;

  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$sample.bowtie.bam", $outted, "DNA", "$sample/bowtie");
  } # end if (parse to VAP)
} #end of subroutine: BOWTIE



sub SAMBAM { #run samtobam
  unless (-e $SAMTOOLS) { die "SAMTOOLS tool path is incorrect, make sure the correct path is specified\n"; }
  ($sample, $reads, $outted) = @_[0..2];
  chdir("$outputfolder");
  `mkdir -p $sample`;
  chdir ("$sample"); 

$codes = <<"ENDCODES";
#SAMtoBAM $sample
chdir("$outputfolder");
`mkdir -p $sample`;
chdir ("$sample");
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";
  $doesexist = (grep /$sample\.bam/, (split("\n", `find ./`)))[0];
  unless ($doesexist) {
    print OUT "`$SAMTOOLS view -bS $reads -o $sample.bam`;\n";
    if ($notify) { print OUT "NOTIFICATION(\"$sample - SAMtoBAM complete\");\n";  }
  } #end unless(doesexist)
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($sample, "$sample.bam", $outted, $_[3], $sample);
  } #end if (parse to VAP)
} #end of subroutine: SAMBAM


sub VAP { #run VAP
  ($sample, $reads, $outted) = @_[0..2];
  chdir("$outputfolder");
  `mkdir -p $_[4]`;
  chdir ("$_[4]");

$codes = <<"ENDCODES";
#VAP $sample
chdir("$outputfolder");
`mkdir -p $_[4]`;
chdir ("$_[4]");
my \$locale=`pwd`; chomp \$locale;
ENDCODES

  open (OUT, ">>$outted");
  print OUT "$codes\n";
 
  if($CONFIGURE{"RUNGATK"} eq "true") {
    unless (-e $PICARD) { die "PICARD tool path is incorrect, make sure the correct path is specified\n"; }
    unless (-e $GATK) { die "GATK tool path is incorrect, make sure the correct path is specified\n"; }
    my $dict = $outputfolder.'/*.dict'; $dict = `ls $dict`; chomp $dict;
    $GATKREF = basename($REF); #( split('/',$REF) )[-1];
    my $testindex = $CONFIGURE{'GENOMEDIR'}.'/*.dict'; $testindex = `ls $testindex`; chomp $testindex;
    unless (-e $CONFIGURE{'GENOMEDIR'}) { die "GENOMEDIR FOLDER specified doesn't exist\n"; }
    unless (-e $testindex) {
      unless (-e $dict) {
        print "Genome Indexes for GATK don't exist in GENOMEDIR specified\nCreating temporary index\n";
        chdir ("$outputfolder");
        $dict = (( split("\.fa",$GATKREF) )[0]).".dict";
        `cp $REF ./; java -jar $PICARD CreateSequenceDictionary R=$GATKREF O=$dict;`;
      }
    } else {
      $GATKREF = $REF;
    }

    print OUT '`mkdir -p $locale/variants/GATK`;',"\n",'chdir("$locale/variants/GATK");',"\n";
		
    #QUALITY SCORE DISTRIBUTION
    $doesexist = (grep /GATK\/qualityscores.txt/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`java -jar $PICARD QualityScoreDistribution INPUT=$reads OUTPUT=qualityscores.txt CHART=qualityscores.chart`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Quality Score Distribution complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Quality Score Distribution previously completed\");\n";  }
    } # end unless+else (doesexist)
		
    #SORT BAM
    $doesexist = (grep /GATK\/aln_sorted.bam/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`java -jar $PICARD SortSam INPUT=$_[1] OUTPUT=aln_sorted.bam SO=coordinate`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Sort Bam complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Sort Bam previously completed\");\n";  }
    } # end unless+else (doesexist)
		
    #ADDREADGROUPS
    $doesexist = (grep /GATK\/aln_sorted_add.bam/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`java -jar $PICARD AddOrReplaceReadGroups INPUT=aln_sorted.bam OUTPUT=aln_sorted_add.bam SO=coordinate RGID=Label RGLB=Label RGPL=illumina RGPU=Label RGSM=Label`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Add read groups complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Add read groups previously completed\");\n";  }
    }
		
    #MARKDUPLICATES
    $doesexist = (grep /GATK\/aln_sorted_mdup.bam/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`java -jar $PICARD MarkDuplicates INPUT=aln_sorted_add.bam OUTPUT=aln_sorted_mdup.bam M=aln_sorted_mdup.metrics CREATE_INDEX=true`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Mark duplicates complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Mark duplicates previously completed\");\n";  }
    }
		
    #REORDER SAM
    $doesexist = (grep /GATK\/aln_resorted_mdup.bam/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`java -jar $PICARD ReorderSam INPUT=aln_sorted_mdup.bam OUTPUT=aln_resorted_mdup.bam REFERENCE=$GATKREF CREATE_INDEX=TRUE`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Resorted Mark duplicates complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Resorted Mark duplicates previously completed\");\n";  }
    }
		
    #specified into RNAseq or DNAseq
    if ($_[3] eq "RNA") {		
      #SPLIT&TRIM
      $doesexist = (grep /GATK\/aln_sorted_split.bam/, (split("\n", `find ./`)))[0];
      unless ($doesexist) {
			
$codes = <<'ENDCODES';
my $file = `tail -n2 qualityscores.txt | head -n 1 | awk -F" " '{print \$1}'`;
if ($file >= 59) {
ENDCODES
$codes .= <<"ENDCODES";
`java -jar $GATK -T SplitNCigarReads --fix_misencoded_quality_scores -R $GATKREF -I aln_resorted_mdup.bam -o aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`;
} else {
`java -jar $GATK -T SplitNCigarReads -R $GATKREF -I aln_resorted_mdup.bam -o aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`;
}
ENDCODES

        print OUT "$codes\n";
        if ($notify) { print OUT "NOTIFICATION(\"$sample - SplitNCigars complete\");\n";  }
      } else {
        if ($notify) { print OUT "NOTIFICATION(\"$sample - SplitNCigars previously completed\");\n";  }
      }
			
      #GATK
      $doesexist = (grep /GATK\/$_[0]_all.vcf/, (split("\n", `find ./`)))[0];
      unless ($doesexist) {
        print OUT "`java -jar $GATK -T HaplotypeCaller -R $GATKREF -I aln_sorted_split.bam -o $_[0]_all.vcf`;\n";
        if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller complete\");\n";  }
      } else {
        if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller previously completed\");\n";  }
      }
    } #end if RNA
    else {  #working with DNA   
      #GATK
      $doesexist = (grep /GATK\/$_[0]_all.vcf/, (split("\n", `find ./`)))[0];
      unless ($doesexist) {
        print OUT "`java -jar $GATK -T HaplotypeCaller -R $GATKREF -I aln_resorted_mdup.bam -o $_[0]_all.vcf`;\n";
        print OUT "`java -jar $GATK -T HaplotypeCaller -R $GATKREF -I aln_resorted_mdup.bam --emitRefConfidence GVCF -o $_[0]_all_emit.vcf`;\n";
        if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller complete\");\n";  }
      } else {
        if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller previously completed\");\n";  }
      }
    } #end else DNA
  } #end unless GATK is false
  if($CONFIGURE{"RUNSAMTOOLS"} eq "true") {
    unless (-e $SAMTOOLS) { die "SAMTOOLS tool path is incorrect, make sure the correct path is specified\n"; }
    unless (-e $BCFTOOLS) { die "BCFTOOLS tool path is incorrect, make sure the correct path is specified\n"; }
    my $dict = $outputfolder.'/*.fai'; $dict = `ls $dict`; chomp $dict;
    $SAMREF = basename($REF); #( split('/',$REF) )[-1];
    my $testindex = $CONFIGURE{'GENOMEDIR'}.'/*.fai'; $testindex = `ls $testindex`; chomp $testindex;
    unless (-e $CONFIGURE{'GENOMEDIR'}) { die "GENOMEDIR FOLDER specified doesn't exist\n"; }
    unless (-e $testindex) {
      unless (-e $dict) {
        print "Genome Indexes for GATK don't exist in GENOMEDIR specified\nCreating temporary index\n";
        chdir ("$outputfolder");
        `cp $REF ./; $SAMTOOLS faidx $SAMREF;`;
      }
    } else {
      $SAMREF = $REF;
    }

    print OUT '`mkdir -p $locale/variants/samtools`;',"\n",'chdir ("$locale/variants/samtools");',"\n";
		
    #SORT BAM
    $doesexist = (grep /samtools\/aln_sorted.bam/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`$SAMTOOLS sort $_[1] -o aln_sorted.bam`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Sort Bam complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Sort Bam previously completed\");\n";  }
    }
		
    #variant call to BCF output
    $doesexist = (grep /samtools\/$_[0]_all.bcf/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`$SAMTOOLS mpileup -f $SAMREF -g aln_sorted.bam > $_[0]_all.bcf`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Variant call to BCF complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - Variant call to BCF previously completed\");\n";  }
    }
		
    #convert to VCF
    $doesexist = (grep /samtools\/$_[0]_all.vcf/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      print OUT "`$BCFTOOLS view -vcg $_[0]_all.bcf > $_[0]_all.vcf`;\n";
      if ($notify) { print OUT "NOTIFICATION(\"$sample - convert BCF to VCF complete\");\n";  }
    } else {
      if ($notify) { print OUT "NOTIFICATION(\"$sample - convert BCF to VCF  previously completed\");\n";  }
    }
  } #end unless SAMTOOLS is false
} #end of subroutine: VAP


## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - V A P - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

=head1 NAME

 -- Variant Analysis Pipeline : For Reference mapping and Variant detection. --

=head1 SYNOPSIS

 VariantdetectionAnalysis.pl -a <configfile> [--help] [--manual]

=head1 DESCRIPTION

 Comprehensive pipeline for variant analysis using a wide suite of bioinformatics tools 
  from reference assembly using HiSAT, TopHAT, STAR, BWA, BOWTIE 
  and variant detection using PICARD, GATK and/or SAMTOOLS using default settings.

=head1 OPTIONS

=over 8

=item B<-a, -c, --config>=FILE

Configuration file (a template can be found below ).  (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 NITTY GRITTY

=over 8

=back 

=head2 Configuration File Variables

 Here is a list of all configuration file parameters currently accepted, 
tools utilized in the pipeline should be specified, else specify false.

 ######============== TOOL LOCATIONS ==============######
	TOPHAT = /path/to/tophat2
	BOWTIE = /path/to/bowtie2
	STAR = /path/to/star
	HISAT = /path/to/hisat2
	BWA = /path/to/bwa/
	PICARD = /path/to/picard.jar
	GATK = /path/to/gatk.jar
 	SAMTOOLS = /path/to/samtools
	BCFTOOLS = /path/to/bcftools
	FASTQC = /path/to/fastqc
 ######============ INPUT FILE OPTIONS ============######
	GFF = /path/to/gff_file
	GTF = /path/to/gtf_file
	GENOME = /path/to/genome_fasta
	GENOMEDIR = /path/to/genome_indexes_directory
	FASTQ = /path/to/fastq_files/*gz [or false]
	SAM = /path/to/sam_files/*sam [or false]
	BAM = /path/to/bam_files/*bam [or false]
 ######================== OPTIONS ==================###### 
	THREADS = 16    #number of threads for alignment process
 ######=============== OUTPUT FOLDER ===============######
 	OUTPUTDIR = /path/to/output_folder
 ######================== WORKFLOW =================######
	sampleRNA = true [or false] #is sample RNA
	sampleDNA = false [or true] #is sample DNA
	#Alignment options
	runTopHat = true [or false] #must specify sampleRNA=true
	runHISAT = true [or false] #must specify sampleRNA=true
	runSTAR = true [or false] #must specify sampleRNA=true
	runBWA = false [or true] #must specify sampleDNA=true
	runBOWTIE = false [or true] #must specify sampleDNA=true
	#Housekeeping + VariantAnalysis 
	runFastqc = true [or false]
	runVAP = false [or false]
	#Variant Caller options
	runGATK = true [or false]
	runSAMTOOLS = true [or false]
 ######=============== NOTIFICATION ===============######
	EMAIL= youremail@address.com
	SUBJECT = Title_of_Email #space sensitive

=over 8

=back

=head2 LOG & ERROR files
 Detail reports for the log and error files for each step of the workflow are stored
  in the <OUTPUTDIR>/tmp folder.

=over 8

=back

=head1 DEPENDENCIES

 Requires the following Perl libraries (all standard in most Perl installs).
  Getopt::Long
  Pod::Usage
  File::Basename
  Time::localtime
  Time::Piece
  File::stat
  threads
  Thread::Queue

=head1 AUTHOR

 Written by Modupeore O. Adetunji, 
 Animal and Food Sciences Department, Schmidt Lab., University of Delaware.
 Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

 Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

 Copyright 2019 MOA.  
 License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
 This is free software: you are free to change and redistribute it.  
 There is NO WARRANTY, to the extent permitted by law.  

 Please acknowledge author and affiliation in published work arising from this script's usage
=cut

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - T H E - E N D - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

