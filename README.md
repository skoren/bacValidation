# bacValidation
BAC-based validation of assemblies

## Requirements
- minimap2 >= 2.17
- R >= 3.5
- samtools >= 1.9
- GCC >= 4.8

These must all be available in your path.

## License

All code developed by me is released into the public domain. Other code included in the repo may have different licenses, see the header present at the top of each file for details.

## Installation
Clone the repo and compile the sam parsing code by running make. It should make a single binary named samToErrorRate. If that exists, the installation succeeded.


## BAC libraries
This pipeline assumes you have a file named bacs.fasta in the folder where you are running from your genome to validate. There are some BAC libraries in NCBI for several human genomes, for example <a href="">CHM13</a>, NA12878, and HG0733.

## Usage

There is a single shell script, getStats.sh which takes one argument, a fasta file for your assembly. The script will align the bacs to your assembly, convert the sam to text, and report some statistics.

## Optional usage

If you want to validate on a subset of BACs, e.g. ones coming from the unique regions of the genome or you have higher confidence in, create a file named goodBacs listing the IDs from the fasta file (up to the first space). The same getStats.sh script above will report stats only on this subset of BACs.
