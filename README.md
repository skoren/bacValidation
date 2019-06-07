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
This pipeline assumes you have a file named bacs.fasta in the folder where you are running from your genome to validate. There are some BAC libraries in NCBI for several human genomes, for example <a href="https://www.ncbi.nlm.nih.gov/nuccore/?term=VMRC59+and+complete">CHM13</a>, <a href="https://www.ncbi.nlm.nih.gov/nuccore/?term=VMRC53+and+complete">NA12878</a>, and <a href="https://www.ncbi.nlm.nih.gov/nuccore/?term=VMRC62+and+complete">HG0733</a>.

## Usage

There is a single shell script, getStats.sh which takes one argument, a fasta file for your assembly. The script will align the bacs to your assembly, convert the sam to text, and report some statistics. Remove the sam and txt files to force the pipeline to re-map the bacs.

## Advanced usage

If you want to validate on a subset of BACs, e.g. ones coming from the unique regions of the genome or you have higher confidence in, create a file named goodBacs listing the IDs from the fasta file (up to the first space). The same getStats.sh script above will report stats only on this subset of BACs. You can do this after running the pipeline, the mappings won't be re-generated.

You can also adjust the default break length when mapping (which defaults to 2kb) to any other value by passing it as the second parameter (e.g. sh getStats.sh asm.fasta 5000 will use a break length of 5kb). This will map more BACs in single pieces but will decrease the identity since it will tolerage more noise in the alignment.

## Example output for <a href="https://github.com/nanopore-wgs-consortium/CHM13">T2T v0.6</a> assembly

```
******************* BAC SUMMARY ******************
 TOTAL    : 341
 BP       : 51532183
************** Statistics for: chm13.draft_v0.6.fasta ****************
BACs closed: 280 (82.1114) BASES 42501309 (82.4753)
Median:		 99.98025 
MedianQV:	 37.04433 
Mean:		 99.80281 
MeanQV:		 27.05109 
***** STATS IGNORING INDELS ********************
Median:		 100 
MedianQV:	 Inf 
Mean:		 99.96021 
MeanQV:		 34.00199 
```

