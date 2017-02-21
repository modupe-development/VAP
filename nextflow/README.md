#Nextflow version of VAP
Thanks for your continued interest in using the Variant Analysis Pipeline.

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
change **_custom_job.txt_** file with settings as required.
If parameters are not needed, they must be changed to *false*
e.g : params.sam = false (this means there is no sam file) or input the file directory


##To run workflow
bash nextVAP.sh custom_job.txt


