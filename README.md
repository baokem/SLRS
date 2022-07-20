# SLRS
## Introduction
SLRS is a scaffolding method based on synthetic long reads. It utilizes co-barcoding information to assemble contigs into longer sequences.
## User Guide
### 1.Before installing and running
```
    Please install BWA from https://github.com/lh3/bwa.
    Please install Samtools from https://sourceforge.net/projects/samtools/files/samtools/.
    Please build and install Bamtools from https://github.com/pezmaster31/bamtools.
```
### 2.Installing
```
    SLRS should run on Linux operating sysetm with gcc. We test SLHSD using gcc9.4.0 on Ubuntu.
        git clone https://github.com/baokem/SLRS.git
	cd SLRS
	export BAMTOOLS_HOME_INCLUDE=~/bamtools/build/usr/local/include/bamtools
	export BAMTOOLS_HOME_LIB=~/bamtools/build/usr/local/lib
	make all
```
### 3.Running
```
    step 1: bwa index contigs.fasta
    step 2: bwa mem -a contigs.fasta R1_reads.fasta > align1.sam
    step 3: samtools view -Sb align1.sam > align1.bam
    step 4: bwa mem -a contigs.fasta R2_reads.fasta > align2.sam
    step 5: samtools view -Sb align2.sam > align2.bam
    step 6: ./SLRS -c contigs.fasta -l align1.bam -r align2.bam -o resultOutPutDirectory [other options] 
    
    -c contig-file: Input file with fasta format
    -l R1-reads-bamfile: R1 reads aligned to contigs
    -r R2-reads-bamfile: R2 reads aligned to contigs
    -a alignmentResultFile: Collected alignment information files
    -o result-file: Result file
    -t contigLengthThreshold: Contigs with length less than this parameter are ignored (default 2000)\n
    -e endlength: Contig head/tail length for collecting aligned reads (default 30000)\n
    -m readsnum:The coefficient of variation is calculated when the number of reads exceeds this parameter (default 7)\n
    -n readsnumperbarcode: The number of reads each barcode has (default 7)\n
    -f endreadsnum: The number of reads aligned to contig end (default 3)\n
    -s sortbarcode: Sort barcodes\n
    
```
### 4.Output
```
    The output file "scaffold_set.fa" is the scaffolding result. 
```
