#!/usr/bin/nextflow
//TOOLS
//RNA
tophat = params.tophat
star = params.star
hisat = params.hisat

//DNA
bwa = params.bwa
bowtie = params.bowtie

//VARIANTS
picard = params.picard
gatk = params.gatk

//OTHERS
samtools = params.samtools
fastqc = params.fastqc

//INPUT FILES
GFF = params.gff
genome = params.genome
genomeDir = params.genomeDir
genomeIndex = params.genomeIndex
params.fastq = "/home/folder/*gz"
params.sam = "/home/folder/sam.sam"
params.bam = "/home/folder/bam.bam"
params.threads = 1

//PARENT OUTPUT FOLDER
params.finalDir = "/home/final"

//PARAMETERS TO RUN
params.runFastqc = true
params.runTopHat = false
params.runHISAT = false
params.runSTAR = false
params.runBWA = false
params.runBOWTIE = false
params.sampleRNA = false
params.sampleDNA = false
params.runVAP = false

//PROCESS INPUT FILE
if (!(params.fastq || params.sam || params.bam)){
	error 'First No input reads specified. Please specify parameters "fastq" or "sam" or "bam"'
}
else if ( params.fastq ) {
  demuxed = Channel.fromPath(params.fastq)
  println "Executing pipeline for \"${params.fastq}\""
	groupedDemux = demuxed
		.map {it -> [sample(it), it ]}
		.groupTuple (sort:true)
		.map{ sample, reads -> tuple (sample, reads[0], reads[1]) }

	groupedDemux.tap{Sample1}
	            .tap{Sample2}
							.tap{Sample3}
							.tap{Sample4}
						
	process fastqc {
	//Objective : run FastQC
	  tag { sample }
	  publishDir "${params.finalDir}/${sample}/fastqc", mode: 'link'
	
	  input:
	    set val(sample), file(read1), file(read2) from Sample1

	  when:
	    params.runFastqc

	  output:
	    set sample, file ('*.html'), file ('*.zip'),
				file ( "${read1.name.replaceFirst(/\.fastq.gz/, '.txt')}"),
				file ( "${read2.name.replaceFirst(/\.fastq.gz/, '.txt')}") into fastqcOut

	  script:
	    """
	    $fastqc ${read1} ${read2}
			unzip ${read1.name.replaceFirst(/\.fastq\.gz/, '_fastqc.zip')}
			unzip ${read2.name.replaceFirst(/\.fastq\.gz/, '_fastqc.zip')}
			cp ${read1.name.replaceFirst(/\.fastq\.gz/, '_fastqc')}/summary.txt ${read1.name.replaceFirst(/\.fastq.gz/, '.txt')}
			cp ${read2.name.replaceFirst(/\.fastq\.gz/, '_fastqc')}/summary.txt ${read2.name.replaceFirst(/\.fastq.gz/, '.txt')}
	    """
	}
	if ( params.sampleRNA ) {
		process alignTopHAT {
			tag { sample }
			publishDir "${params.finalDir}/${sample}/tophat", mode: 'link'
			
			input:
				set val(sample), file(read1), file(read2) from Sample2
	
			output:
				file '*' into alignmentOut1
				set val(sample), file ('*tophat.bam'), val('tophat') into BamOut1
				set val(sample), val("${sample}-tophat"), val('tophat') into BamStat1
			
			when:
				params.runTopHat
			
			shell:
				"""
				$tophat \
					-p ${threads} \
					--library-type fr-firststrand \
					--no-coverage-search \
					-G $GFF \
					-o ./ \
					$genomeIndex \
					${read1} ${read2}
				mv accepted_hits.bam ${sample}.tophat.bam
				
				"""
		}
	
		process alignSTAR {
			tag { sample }
			publishDir "${params.finalDir}/${sample}/star", mode: 'link'
			
			input:
				set val(sample), file(read1), file(read2) from Sample3
	
			output:
				file '*' into alignmentOut2
				set val(sample), file ('*Aligned.out.bam'), val('star') into BamOut2
				set val(sample), val("${sample}-star"), val('star') into BamStat2
				
			
			when:
				params.runSTAR
			
			shell:
				"""
				cp ${read1} ${read1.name.replaceFirst(/\.gz/, '2.gz')}
				cp ${read2} ${read2.name.replaceFirst(/\.gz/, '2.gz')}
				gunzip ${read1.name.replaceFirst(/\.gz/, '2.gz')}
				gunzip ${read2.name.replaceFirst(/\.gz/, '2.gz')}
                                mv ${read1.name.replaceFirst(/\.gz/, '2')} ${read1.name.replaceFirst(/\.gz/, '')}
                                mv ${read2.name.replaceFirst(/\.gz/, '2')} ${read2.name.replaceFirst(/\.gz/, '')}
				$star \
					--runThreadN 16 \
					--genomeDir $genomeDir \
					--outFileNamePrefix ${sample}. \
					--readFilesIn ${read1.name.replaceFirst(/\.gz/, '')} \
					${read2.name.replaceFirst(/\.gz/,'')}
				$samtools view -bS ${sample}.Aligned.out.sam -o ${sample}.Aligned.out.bam
				"""
		}
	
		process alignHiSAT {
			tag { sample }
			publishDir "${params.finalDir}/${sample}/hisat", mode: 'link'
			
			input:
				set val(sample), file(read1), file(read2) from Sample4
	
			output:
				file '*' into alignmentOut3
				set val(sample), file ('*hisat.bam'), val('hisat') into BamOut3
				set val(sample), val("${sample}-hisat"), val('hisat') into BamStat3
			
			when:
				params.runHISAT
		
			shell:
				"""
				$hisat \
					-p ${threads} \
					-x $genomeIndex \
					-1 ${read1} -2 ${read2} \
					-S ${sample}.hisat.sam \
					&> ${sample}_align.txt
					$samtools view -bS ${sample}.hisat.sam -o ${sample}.hisat.bam
				"""
		}
		BamOut = BamOut1.mix( BamOut2 )
						.mix( BamOut3 )
		BamStat = BamStat1.mix( BamStat2 )
						  .mix( BamStat3 )
	}
	else if ( params.sampleDNA ) {
		process alignBWA {
			tag { sample }
			publishDir "${params.finalDir}/${sample}/BWA", mode: 'link'
			
			input:
				set val(sample), file(read1), file(read2) from Sample2
	
			output:
				file '*' into alignmentOut1
				set val(sample),file ('*bwa.bam'),val('BWA') into BamOut1
				set val(sample), val("${sample}-bwa"), val('BWA') into BamStat1
			
			when:
				params.runBWA
			
			shell:
				"""
				$bwa \
					mem \
					-t ${threads} \
					$genome \
					${read1} ${read2} \
					> ${sample}_bwa.sam
					$samtools view -bS ${sample}_bwa.sam -o ${sample}_bwa.bam
				"""
		}
		
		process alignBOWTIE {
			tag { sample }
			publishDir "${params.finalDir}/${sample}/BOWTIE", mode: 'link'
			
			input:
				set val(sample), file(read1), file(read2) from Sample3
	
			output:
				file '*' into alignmentOut2
				set val(sample), file ('*bowtie.bam'), val('BOWTIE') into BamOut2
				set val(sample), val("${sample}-bowtie"), val('BOWTIE') into BamStat2
			
			when:
				params.runBOWTIE
			
			shell:
				"""
				$bowtie \
					-p ${threads} \
					-x $genomeIndex \
					-1 ${read1} -2 ${read2} \
					-S ${sample}_bowtie.sam
					$samtools view -bS ${sample}_bowtie.sam -o ${sample}_bowtie.bam
				"""
		}
		BamOut = BamOut1.mix( BamOut2 )
		BamStat = BamStat1.mix( BamStat2 )
	}
} 
else if ((params.sam) || (params.bam)) { //Dealing with a bam or sam file
	if (params.sam) {
		demuxed = Channel.fromPath(params.sam)
		println "Executing pipeline for \"${params.sam}\""
		SamMuxed = demuxed
			.map {it -> [samheader(it), it ]}
			.groupTuple (sort:true)
			.map{ sample, reads -> tuple (sample, reads[0]) }
		
		process convertBAM {
		//Objective : convert sam to bam file
			tag { sample }
			publishDir "${params.finalDir}/${sample}/", mode: 'link'
		
			input:
				set val(sample), file(reads) from SamMuxed
	
			output:
				set val(sample), file ('*bam'), val('.'), into ConvertBam
				set val(sample), val(sample), val('.') into BamStat				
	
			script:
				"""
				$samtools view -bS ${reads} -o ${sample}.bam
				"""
		}
	} else if (params.bam) {
		demuxed = Channel.fromPath(params.bam)
		demuxed.tap{demuxed1}
			   .tap{demuxed2}
		println "Executing pipeline for \"${params.bam}\""
		ConvertBam = demuxed1
			.map {it -> [samheader(it), it ]}
			.groupTuple (sort:true)
			.map{ sample, reads -> tuple (sample, reads[0]) }
			.spread(['.'])
		BamStat = demuxed2
			.map {it -> [samheader(it), it ]}
			.groupTuple (sort:true)
			.map{ sample, reads -> tuple (sample, sample) }
			.spread(['.'])
	}
	
	ConvertBam.tap{Convert}
			  .tap{BamOut}
							 
	process fastqc {
	//Objective : run FastQC on Bam file 
	  tag { sample }
	  publishDir "${params.finalDir}/${sample}/${folder}/fastqc", mode: 'link'
	
	  input:
	    set val(sample), file(reads), val(folder) from Convert

	  when:
	    params.runFastqc

	  output:
	    set sample, file ('*.html'), file ('*.zip') into fastqcOut

	  script:
	    """
	    $fastqc ${reads}
            unzip ${reads.name.replaceFirst(/\.bam/, '_fastqc.zip')}
	    cp ${reads.name.replaceFirst(/\.bam/, '_fastqc')}/summary.txt ${reads.name.replaceFirst(/\.fastq.gz/, '.txt')}
	    """
	}
} else {
	error 'Please specify parameters "fastq" or "bam" or "sam"'
}
								
if ( params.runVAP ) {
	BamOut.tap{ NewBam1 }
			.tap{ NewBam2 }
	BamStat.tap{ NewBamStat1 }
			.tap{ NewBamStat2 }
	
	process QualityScore {
	//Objective : create sequence dictionary
		tag {stat}
		publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
	
	  input:
	    set val(sample), file(bam), val(folder) from NewBam1
		set val(sample), val(stat), val(folder) from NewBamStat1

	  output:
	    set sample, file ('*') into QualityScoreOut
		set sample, folder, file('qualityscores.txt') into CheckQuality

	  script:
	    """
	    java -jar $picard \
				QualityScoreDistribution \
				INPUT=${bam} \
				OUTPUT=qualityscores.txt \
				CHART=qualityscores.chart
	    """
	}
	
	process SortBam {
	//Objective : sort bam file in co-ordinate
		tag {stat}
		publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
	
	  input:
	    set val(sample), file(bam), val(folder) from NewBam2
		set val(sample), val(stat), val(folder) from NewBamStat2

	  output:
	    set sample, file('*') into SortBamAll
		set sample, stat, folder, file ('aln_sorted.bam') into SortBamOut

	  script:
	    """
	    java -jar $picard \
				SortSam \
				INPUT=${bam} \
				OUTPUT=aln_sorted.bam \
				SO=coordinate
	    """
	}
  
	process AddReadGroup {
	//Objective : add default read groups
		tag {stat}
		publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
	
	  input:
	    set val(sample), val(stat), val(folder), file(bam) from SortBamOut

	  output:
	    set sample, file ('*') into AddReadAll
			set sample, stat, folder, file ('aln_sorted_add.bam') into AddReadOut

	  script:
	    """
	    java -jar $picard \
				AddOrReplaceReadGroups \
				INPUT=${bam} \
				OUTPUT=aln_sorted_add.bam \
				SO=coordinate \
				RGID=Label RGLB=Label RGPL=illumina RGPU=Label RGSM=Label
	    """
	}
	
	process MarkDuplicates {
	//Objective : mark duplicates
		tag {stat}
		publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
	
	  input:
	    set val(sample), val(stat), val(folder), file(bam) from AddReadOut

	  output:
	    set sample, file ('*') into MarkDupAll
			set sample, stat, folder, file ('aln_sorted_mdup.bam'), file ('aln_sorted_mdup.bai') into MarkDupOut

	  script:
	    """
	    java -jar $picard \
				MarkDuplicates \
				INPUT=${bam} \
				OUTPUT=aln_sorted_mdup.bam \
				M=aln_sorted_mdup.metrics \
				CREATE_INDEX=true
	    """
	}
	
	process ReOrderSam {
	//Objective : reordering bam file to match with reference file
		tag {stat}
		publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
	
	  input:
	    set val(sample), val(stat), val(folder), file(bam), file(bai) from MarkDupOut

	  output:
	    set sample, file ('*') into ReOrderAll
			set sample, stat, folder, file ('aln_resorted_mdup.bam'), file ('aln_resorted_mdup.bai') into ReOrderOut

	  script:
	    """
	    java -jar $picard \
			ReorderSam \
			INPUT=${bam} \
			OUTPUT=aln_resorted_mdup.bam \
			REFERENCE=$genome \
			CREATE_INDEX=true
	    """
	}
	
	ReOrderOut.tap{ReOrderOut1}
						.tap{ReOrderOut2}
						
	if (params.sampleRNA) {
		process SplitTrim {
		//Objective : split and trim reads
			tag {stat}
			publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
		
			input:
				set sample, folder, file(quality) from CheckQuality
				set val(sample), val(stat), val(folder), file(bam), file(bai) from ReOrderOut1
	
			output:
				set sample, file ('*') into SplitTrimAll
				set sample, stat, folder, file ('aln_sorted_split.bam'), file ('aln_sorted_split.bai') into SplitTrimOut
	
			shell:
				'''
				for i in `tail -n2 !{quality} | head -n 1 | awk -F" " '{print $1}'`; 
					do
						if [ $i -ge 59 ]; then
							`java -jar !{gatk} \
								-T SplitNCigarReads \
								--fix_misencoded_quality_scores \
								-R !{genome} -I !{bam} \
								-o aln_sorted_split.bam \
								-rf ReassignOneMappingQuality \
								-RMQF 255 -RMQT 60 --filter_reads_with_N_cigar`
						else 
							`java -jar !{gatk} \
							-T SplitNCigarReads \
							-R !{genome} -I !{bam} \
							-o aln_sorted_split.bam \
							-rf ReassignOneMappingQuality \
							-RMQF 255 -RMQT 60 --filter_reads_with_N_cigar` 
						fi
					done
				'''
		}
		
		process HaplotypeCallRNA {
		//Objective : call variants
			tag {stat}
			publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
		
			input:
				set val(sample), val(stat), val(folder), file(bam), file(bai) from SplitTrimOut
	
			output:
				set sample, file ('*') into VariantsAll
				
			script:
				"""
				java -jar $gatk \
					-T HaplotypeCaller \
					-R $genome \
					-I ${bam} \
					-o ${sample}_all.vcf
				"""
		}
	}
	else if (params.sampleDNA) {
		process HaplotypeCallDNA {
		//Objective : call variants
			tag {stat}
			publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
		
			input:
				set val(sample), val(stat), val(folder), file(bam), file(bai) from ReOrderOut1

			output:
				set sample, file ('*') into VariantsAll
				
			script:
				"""
				java -jar $gatk \
					-T HaplotypeCaller \
					-R $genome \
					-I ${bam} \
					-o ${sample}_all.vcf
				"""
		}
		process HC_DNAEmit {
		//Objective : call variants; emit reference confidence
			tag {sample}
			publishDir "${params.finalDir}/${sample}/${folder}", mode: 'link'
		
			input:
				set val(sample), val(stat), val(folder), file(bam), file(bai) from ReOrderOut2

			output:
				set sample, file ('*') into VariantsAllEmit
				
			script:
				"""
				java -jar $gatk \
					-T HaplotypeCaller \
					-R $genome \
					-I ${bam} \
					--emitRefConfidence GVCF \
					-o ${sample}_all_emit.vcf
				"""
		}
	} else {
		error 'Specify if your sample is DNA or RNA in job file'
	}
}

/*=================================================================================================================

//============================================ WORKFLOW ROUTINES ==================================================*/
////DEFINITIONS
def sample(fileproc) {
  def procfile = fileproc.getFileName().toString()
  if(procfile =~ /.+[_\.][Rr][12].+/) {
    value = procfile =~ /(.+)_.*[_\.][Rr][12].+/
  }
  else if(procfile =~ /.+[_\.][12].+/) {
    value = procfile =~ /(.+)_.*[_\.][12].+/
  }
  else if(procfile =~ /.+[_\.]pe[12]/) {
    value = procfile =~ /(.+)_.*[_\.]pe[12].+/
  }
  else if(procfile =~ /.+[_\.]PE[12]/) {
    value = procfile =~ /(.+)_.*[_\.]PE[12].+/
  }
  return value[0][1]
}
def samheader(fileproc) { //get the prefix before output
  def procfile = fileproc.getFileName().toString()
  if(procfile =~ /.+[\.]sam/) {
    value = procfile =~ /(.+)[\.]sam/
  }
	else if(procfile =~ /.+[\.]bam/) {
    value = procfile =~ /(.+)[\.]bam/
  }
  return value[0][1]
}

////ON COMPLETION
workflow.onComplete {
    println ""
    println "WORKFLOW SUMMARY"
    println "Pipeline started at: $workflow.start"
    println "Pipeline completed at: $workflow.complete"
    println "Pipeline duration: $workflow.duration"
    println "Results Directory: ${params.finalDir}"
    println "Error message: $workflow.errorMessage"
    println "Temporary Directory : $workflow.workDir"
    println "Execution status: ${ workflow.success ? 'SUCCESS' : 'failed' }"
}

/*=================================================================================================================*/
