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
	cd SLHSD
	export BAMTOOLS_HOME_INCLUDE=/path_bamtools_include_api_shared/
	export BAMTOOLS_HOME_LIB=/path_bamtools_lib_libbamtools.a/
	make all
```
