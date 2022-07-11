#!/bin/zsh

### SRA Toolkit can be used to automatically download sequencing files from NCBI SRA webpage
### The files can also be downloaded manually

# download one set of fastq files using SRA toolkit
fasterq-dump --split-files SRR15852399