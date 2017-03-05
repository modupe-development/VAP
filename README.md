#VAP
Thanks for your interest in using the Variant Analysis Pipeline.

**N.B.** : parameters of all tools are set to default.
	Contact maintainer to make custom changes to the different tools if needed.

##Indexes of assembly and variant tools
Before running the pipeline. Create indexes for the different assemblers specified
**SYNTAXS:**
- FASTA INDEX :
```
	samtools faidx <reference.fa>
```
- GATK : 
```
	java -jar <picard directory>/picard.jar CreateSequenceDictionary R=<reference.fa> O=<dictionary.dict>
```
- HISAT :
```
	hisat2-build <reference.fa> <path to reference.fa>/<index_name>
```
- BOWTIE/TOPHAT :
```	
	bowtie2-build <reference.fa> <path to reference.fa>/<index_name>
```
- BWA :
```
	bwa index <reference.fa> 
```
- STAR :
```
	STAR --runMode genomeGenerate --genomeDir <path to reference.fa> --genomeFastaFiles <reference.fa> [--sjdbGTFfile <annotation [gff|gtf] file>] [--sjdbOverhang <reads average length if other than 100>]
```
_**N.B.** All ```<index_name>``` should be the same and stored in the ```<reference.fa>``` directory_

##Job File
change **_custom_job.config_** file with settings as required.
If parameters are not needed, they must be changed to *false*
Needed workflows (prefix: run) must be change to *true* (case-sensitive)
e.g : SAM = false (this means there is no sam file) or input the file directory
e.g : runTopHAT = true (this means the pipeline should run TopHAT

##To run workflow
perl VariantAnalysisPipeline.pl -c custom_job.config

