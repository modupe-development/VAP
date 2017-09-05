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

#tools
our ($TOPHAT, $BOWTIE, $PICARD, $GATK, $HISAT, $FASTQC, $STAR, $SAMTOOLS, $BCFTOOLS, $BWA);

#files
our ($value, %CONFIGURE, %FILE, @content, @arrayfile, @threads);
#inputs
our ($REF, $ANN, $outputfolder, $THREADS);
our ($sample, $reads, $doesexist, @eachread, $codes, $outted, $newout);
my ($email, $notify, $syntax); #notification options


my $date = `date +%m-%d-%y_%T`; chomp $date;
my $std_out = "VAP-$date.log";
my $std_err = "VAP-$date.err";

#ARGUMENTS
my($help,$manual,$config,$bamfile, $flag);
GetOptions (
                                "i|c|a|config=s"    =>      \$config,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# VALIDATE ARGS
pod2usage( -verbose => 2 )  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument '-c <config_file>' not found.\n", -exitval => 2, -verbose => 1)  if (! $config);
@ARGV == 0 || pod2usage("ERROR! Additional comments '@ARGV' not required\n");

configureinput(); #configure input options
INPUTFILES();
#stdout and stderr
`mkdir -p $outputfolder`;
chdir ("$outputfolder"); #working from the outputfolder;
open(STDOUT, '>', "$outputfolder/$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$outputfolder/$std_err") or die "Error file doesn't exist";
 

#EMAILS
my $subject = "VAP-$date"; my $notification = $outputfolder."/".$subject.'.log';
if (length($CONFIGURE{"SUBJECT"}) >= 1) { $subject = $CONFIGURE{"SUBJECT"}; }
if (length ($CONFIGURE{"EMAIL"}) >= 1) {
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
}
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

#PROCESSING
if ($CONFIGURE{"FASTQ"} ne "false") {
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
} #end if fastq file
elsif ( ($CONFIGURE{"BAM"} ne "false") || ($CONFIGURE{"SAM"} ne "false") ) { #if bam or sam file
  my $nac = "RNA"; #default is RNA
  if ($CONFIGURE{"SAMPLEDNA"} eq "true") { $nac = "DNA"; } #if dna
  elsif ($CONFIGURE{"SAMPLERNA"} eq "true") { $nac = "RNA"; } #if rna
  foreach (keys %FILE){
    if ($CONFIGURE{"RUNFASTQC"} eq "true"){
			$newout = "$outputfolder/tmp/$_-fastqc.txt";
			push @arrayfile, $newout;
			open (OUT, ">$newout"); print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
			FASTQC($_, $FILE{$_}, $newout);
		} #fastqc
    if ($CONFIGURE{"RUNVAP"} eq "true"){ #run variant analysis pipeline
      if ($CONFIGURE{"SAM"} ne "false") { #if sam file
				$newout = "$outputfolder/tmp/$_-sam.txt";
				push @arrayfile, $newout;
				open (OUT, ">$newout");
				print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
				SAMBAM($_, $FILE{$_}, $newout, $nac);
			} else {
				$newout = "$outputfolder/tmp/$_-bam.txt";
				push @arrayfile, $newout;
				open (OUT, ">$newout");
				print OUT "#!/usr/bin/perl\n$endnotify\n"; close (OUT);
				VAP($_, $FILE{$_}, $newout,$nac, $_);
			} #end if sam or bam file
    } #end if run vap
  } # end file
} #end if bam or sam file

my $queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..5; #execute 10 threads
$builder->join; #join threads
foreach (@threads){$_->join;}

#end of job
`echo "Job Completed" >>  $outputfolder/welcome-$date.log` if $notify;
system "mail -s \"VAP - $subject\" $email < $outputfolder/welcome-$date.log" if $notify;


##SUBROUTINES
sub main { #main thread routine
  foreach my $count (0..$#arrayfile) {
    while(1) {
      if ($queue->pending() <100) {
        $queue->enqueue($arrayfile[$count]);
        last;
      }
    }
	}
	foreach(1..5) { $queue-> enqueue(undef); }
}

sub processor { #queueing the thread
	my $query;
	while ($query = $queue->dequeue()){
		parseinput($query);
	}
}

sub parseinput{ #what is done in each thread
  print "Submitted\t$_[0]\t",`date +%x-%X`;
  system("perl $_[0] 1>$_[0].log 2>$_[0].err");
	system("echo 'from $_[0].log' >> $outputfolder/$std_out");
  system("echo 'from $_[0].err' >> $outputfolder/$std_err");
  system("cat $_[0].log >> $outputfolder/$std_out");
  system("cat $_[0].err >> $outputfolder/$std_err");
	system("cat $_[0].notice.log >> $outputfolder/welcome-$date.log");
	print "Completed\t$_[0]\t",`date +%x-%X`;
}

sub configureinput {
  #read the config file
  open(CONFIG, "<", $config) or die "Configuration File \"$config\" can't be found\nTERMINATED!\n";
  while (<CONFIG>){
    chomp;
    unless ($_ =~ /^#/){
      if ($_ =~ /\=/){
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
  $REF = $CONFIGURE{"GENOME"};
  $ANN = $CONFIGURE{"GFF"};
  if ($CONFIGURE{"THREADS"} > 1 ){ $THREADS = $CONFIGURE{"THREADS"};} else {$THREADS = 1;}
  $outputfolder = $CONFIGURE{"OUTPUTDIR"};
}

sub INPUTFILES { #sorting the input files
  #FASTQ files
  unless ($CONFIGURE{"FASTQ"} eq "false") {
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
  }
 if ( ($CONFIGURE{"BAM"} ne "false") || ($CONFIGURE{"SAM"} ne "false") ) {
    if ($CONFIGURE{"BAM"} eq "false") {
      @content = sort {$a <=> $b || $a cmp $b} (split("\n", `ls $CONFIGURE{"SAM"}`)); #get details of the folder
    } else {
      @content = sort {$a <=> $b || $a cmp $b} (split("\n", `ls $CONFIGURE{"BAM"}`)); #get details of the folder
    }
    foreach (@content) {
      my $file = basename($_);
      if($file =~ /(.+)[\.]sam/) {
        $value = $1;
      } elsif($file =~ /(.+)[\.]bam/) {
        $value = $1;
      }
      $FILE{$value} = $_;
    }
  }
}

sub FASTQC { #run fastqc
  ($sample, $reads, $outted) = @_;
	chdir("$outputfolder");
	`mkdir -p $sample/fastqc`;
	chdir ("$sample/fastqc"); 
$codes = <<"ENDCODES";
#FastQC $sample
chdir("$outputfolder");
`mkdir -p $sample/fastqc`;
chdir ("$sample/fastqc");
ENDCODES

open (OUT, ">>$outted");
print OUT "$codes\n";

  @eachread = split(/\s/, $reads);
  foreach my $single (@eachread) {
    my $singlefolder = fileparse($single);
    if ($singlefolder =~ /fastq\.gz/) { $singlefolder =~ s/\.fastq\.gz/_fastqc/g; }
		elsif ($singlefolder =~ /fastq$/) { $singlefolder =~ s/\.fastq/_fastqc/g; }
		elsif ($singlefolder =~ /sam$/) { $singlefolder =~ s/\.sam/_fastqc/g; }
		elsif ($singlefolder =~ /bam$/) { $singlefolder =~ s/\.bam/_fastqc/g; }
		else {die "Can't find fastq file type\n";}
    $doesexist = (grep /$singlefolder/, (split("\n", `find ./`)))[0];
    unless ($doesexist) {
      $singlefolder = fileparse($single);
      
$codes = <<"ENDCODES";
`cp $reads ./`;
`$FASTQC $singlefolder`;
ENDCODES
if ($singlefolder =~ /fastq\.gz/) { $singlefolder =~ s/\.fastq\.gz/_fastqc/g; }
elsif ($singlefolder =~ /fastq$/) { $singlefolder =~ s/\.fastq/_fastqc/g; }
elsif ($singlefolder =~ /sam$/) { $singlefolder =~ s/\.sam/_fastqc/g; }
elsif ($singlefolder =~ /bam$/) { $singlefolder =~ s/\.bam/_fastqc/g; }
else {die "Can't find fastq file type\n";}
print OUT "$codes\n";
$codes = <<"ENDCODES";
`unzip $singlefolder.zip`;
`cp $singlefolder/summary.txt $singlefolder.txt`;
`rm -rf $singlefolder`;
ENDCODES

print OUT "$codes\n";

      if ($notify) { print OUT "NOTIFICATION(\"$sample - FastQC complete\");\n"; }
    }
  }
  close OUT;
}

sub TOPHAT { #runtophat
  ($sample, $reads, $outted) = @_;
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

  $doesexist = (grep /tophat\.bam/, (split("\n", `find ./`)))[0];
  unless ($doesexist) {
    
$codes = <<"ENDCODES";
`$TOPHAT -p $THREADS --library-type fr-firststrand --no-coverage-search -G $ANN -o ./ $CONFIGURE{"GENOMEINDEX"} $reads`;
`mv accepted_hits.bam $sample.tophat.bam`;
ENDCODES

print OUT "$codes\n";

    if ($notify) { print OUT "NOTIFICATION(\"$sample - TOPHAT complete\");\n";  }
  }
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$outputfolder/$sample/tophat/$sample.tophat.bam", $outted, "RNA", "$sample/tophat");
  }
  
}

sub STAR { #runstar
  ($sample, $reads, $outted) = @_;
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
    my $starstat = "$STAR --runThreadN $THREADS --genomeDir $CONFIGURE{'GENOMEDIR'} --outFileNamePrefix $sample. --readFilesIn ";
    foreach my $single (@eachread) {
      my $singlefolder = fileparse($single);
      $singlefolder =~ s/q\.gz/q/g;
      $starstat .= "$singlefolder ";
    }

$codes = <<"ENDCODES";
`cp $reads ./`;
`gunzip *gz`;
`$starstat`;
`$SAMTOOLS view -bS $sample.Aligned.out.sam -o $sample.star.bam`;
ENDCODES

print OUT "$codes\n";

    if ($notify) { print OUT "NOTIFICATION(\"$sample - STAR complete\");\n";  }
  }
  close OUT;
  if ($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$outputfolder/$sample/star/$sample.star.bam", $outted, "RNA","$sample/star");
  }
  
}

sub HISAT { #runhisat
  ($sample, $reads, $outted) = @_;
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
			$hisatstat = "$HISAT -p $THREADS -x $CONFIGURE{'GENOMEINDEX'} -S $sample.hisat.sam -1 $eachread[0] -2 $eachread[1] 2> $sample"."_align.txt";
		} else {
			$hisatstat = "$HISAT -p $THREADS -x $CONFIGURE{'GENOMEINDEX'} -S $sample.hisat.sam -U $eachread[0] 2> $sample"."_align.txt";
		}
$codes = <<"ENDCODES";
`$hisatstat`;
`$SAMTOOLS view -bS $sample.hisat.sam -o $sample.hisat.bam`;
ENDCODES

print OUT "$codes\n";

    if ($notify) { print OUT "NOTIFICATION(\"$sample - HiSAT complete\");\n";  }
  }
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$outputfolder/$sample/hisat/$sample.hisat.bam", $outted, "RNA", "$sample/hisat");
  }
}

sub BWA { #run bwa
  ($sample, $reads, $outted) = @_;
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
  }
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$sample.bwa.bam", $outted, "DNA", "$sample/bwa");
  }
}       

sub BOWTIE { #run bowtie
  ($sample, $reads, $outted) = @_;
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
    my $bowtiestat = "$BOWTIE -p $THREADS -x $CONFIGURE{'GENOMEINDEX'} -S $sample.bowtie.sam -1 $eachread[0] -2 $eachread[1]";
$codes = <<"ENDCODES";
`$bowtiestat`;
`$SAMTOOLS view -bS $sample.bowtie.sam -o $sample.bowtie.bam`;
ENDCODES

print OUT "$codes\n";
    if ($notify) { print OUT "NOTIFICATION(\"$sample -BOWTIE complete\");\n";  }
  }
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($_, "$sample.bowtie.bam", $outted, "DNA", "$sample/bowtie");
  }
}       


sub SAMBAM { #run samtobam
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
  }
  close OUT;
  if($CONFIGURE{"RUNVAP"} eq "true") {
    VAP($sample, "$sample.bam", $outted, $_[3], $sample);
  }
}

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
		print OUT '`mkdir -p $locale/variants/GATK`;',"\n",'chdir("$locale/variants/GATK");',"\n";
		
		#QUALITY SCORE DISTRIBUTION
		$doesexist = (grep /GATK\/qualityscores.txt/, (split("\n", `find ./`)))[0];
		unless ($doesexist) {
			print OUT "`java -jar $PICARD QualityScoreDistribution INPUT=$reads OUTPUT=qualityscores.txt CHART=qualityscores.chart`;\n";
			if ($notify) { print OUT "NOTIFICATION(\"$sample - Quality Score Distribution complete\");\n";  }
		} else {
			if ($notify) { print OUT "NOTIFICATION(\"$sample - Quality Score Distribution previously completed\");\n";  }
		}
		
		#SORT BAM
		$doesexist = (grep /GATK\/aln_sorted.bam/, (split("\n", `find ./`)))[0];
		unless ($doesexist) {
		print OUT "`java -jar $PICARD SortSam INPUT=$_[1] OUTPUT=aln_sorted.bam SO=coordinate`;\n";
			if ($notify) { print OUT "NOTIFICATION(\"$sample - Sort Bam complete\");\n";  }
		} else {
			if ($notify) { print OUT "NOTIFICATION(\"$sample - Sort Bam previously completed\");\n";  }
		}
		
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
			print OUT "`java -jar $PICARD ReorderSam INPUT=aln_sorted_mdup.bam OUTPUT=aln_resorted_mdup.bam REFERENCE=$REF CREATE_INDEX=TRUE`;\n";
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
`java -jar $GATK -T SplitNCigarReads --fix_misencoded_quality_scores -R $REF -I aln_resorted_mdup.bam -o aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`;
} else {
`java -jar $GATK -T SplitNCigarReads -R $REF -I aln_resorted_mdup.bam -o aln_sorted_split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`;
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
				print OUT "`java -jar $GATK -T HaplotypeCaller -R $REF -I aln_sorted_split.bam -o $_[0]_all.vcf`;\n";
				if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller complete\");\n";  }
			} else {
				if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller previously completed\");\n";  }
			}
		}
		else {  #working with DNA   
			#GATK
			$doesexist = (grep /GATK\/$_[0]_all.vcf/, (split("\n", `find ./`)))[0];
			unless ($doesexist) {
				print OUT "`java -jar $GATK -T HaplotypeCaller -R $REF -I aln_resorted_mdup.bam -o $_[0]_all.vcf`;\n";
				print OUT "`java -jar $GATK -T HaplotypeCaller -R $REF -I aln_resorted_mdup.bam --emitRefConfidence GVCF -o $_[0]_all_emit.vcf`;\n";
				if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller complete\");\n";  }
			} else {
				if ($notify) { print OUT "NOTIFICATION(\"$sample - Haplotype caller previously completed\");\n";  }
			}
		}
	} #end unless GATK is false
	if($CONFIGURE{"RUNSAMTOOLS"} eq "true") {
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
		print OUT "`$SAMTOOLS mpileup -f $REF -g aln_sorted.bam > $_[0]_all.bcf`;\n";
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
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
## MANUAL FOR VariantdetectionAnalysis.pl
#
#=pod
#
#=head1 NAME
#
#$0 -- Comprehensive pipeline : Variant detection using PICARD, GATK and produces and output folder.
#
#=head1 SYNOPSIS
#
#VariantdetectionAnalysis.pl -a configfile [--help] [--manual]
#
#=head1 DESCRIPTION
#
#Accepts all folders from frnakenstein output.
# 
#=head1 OPTIONS
#
#=over 3
#
#=item B<-a, -c, --config>=FILE
#
#Configuration file (a template can be found @ .. ).  (Required)
#
#=item B<-h, --help>
#
#Displays the usage message.  (Optional) 
#
#=item B<-man, --manual>
#
#Displays full manual.  (Optional) 
#
#=back
#
#=head1 DEPENDENCIES
#
#Requires the following Perl libraries (all standard in most Perl installs).
#   Getopt::Long
#   Pod::Usage
#use strict;
#use File::Basename;
#use Getopt::Long;
#use Time::localtime;
#use Pod::Usage;
#use Time::Piece;
#use File::stat;
#use DateTime;
#use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
#
#=head1 AUTHOR
#
#Written by Modupe Adetunji, 
#Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.
#
#=head1 REPORTING BUGS
#
#Report bugs to amodupe@udel.edu
#
#=head1 COPYRIGHT
#
#Copyright 2015 MOA.  
#License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
#This is free software: you are free to change and redistribute it.  
#There is NO WARRANTY, to the extent permitted by law.  
#
#Please acknowledge author and affiliation in published work arising from this script's usage
#=cut
#
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
