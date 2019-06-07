#!/bin/bash

ASM=$1
BREAK_LEN=2000
if [ $# -ge 2 ]; then
   BREAK_LEN=$2
fi

PREFIX=`echo $1 |sed s/.fasta//g`
PREFIX="${PREFIX}_${BREAK_LEN}"
SCRIPT_PATH=$BASH_SOURCE
SCRIPT_PATH=`dirname $SCRIPT_PATH`

SAM=$PREFIX.sam
cores=`grep -c ^processor /proc/cpuinfo`

if [ x$ASM == "x" ]; then
   echo "Error: please provide an input assembly to validate"
   exit
fi

if [ ! -e bacs.fasta ]; then
  echo "Error: expect a file named bacs.fasta in current folder!"
  exit
fi

if [ ! -e $ASM ]; then
   echo "Error: cannot read assembly $ASM"
   exit
fi

if [ ! -e $PREFIX.txt ]; then
   echo -ne "Starting mapping for asm $ASM using a break len of $BREAK_LEN...."
   minimap2 --secondary=no -t $cores -ax asm20 -r $BREAK_LEN $ASM bacs.fasta > $SAM
   # do this bam back to sam to add header in case asm is >4gb
   samtools view -b -T $ASM -o $PREFIX.tmp $SAM
   samtools view -h $PREFIX.tmp > $SAM
   rm -f $PREFIX.tmp
   # now get the txt
   $SCRIPT_PATH/samToErrorRate $SAM $ASM > $PREFIX.txt
   echo "Done"
fi

samtools faidx bacs.fasta
total=`grep -c ">" bacs.fasta`
bp=`cat bacs.fasta.fai |awk '{SUM+=$2; } END {print SUM}'`

echo "******************* BAC SUMMARY ******************"
echo " TOTAL    : $total"
echo " BP       : $bp"

echo "************** Statistics for: $ASM ****************"
cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) print $1}'|sort |uniq |awk '{print $1}' > $PREFIX.resolved 
DONE=`cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) print $1}'|sort |uniq |awk '{print $1}' |wc -l`
DONEBP=`cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) print $1" "$8}'|sort |uniq |awk '{SUM+=$NF; } END { print SUM}'`
echo "$DONE $total $DONEBP $bp" |awk '{print "BACs closed: "$1" ("$1/$2*100") BASES "$3" ("$3/$4*100")"}'
   
echo "" > Rstats
echo "dat = read.table('$PREFIX.R')" >> Rstats
echo "cat('Median:\t\t',  median(dat[,1]), '\n')" >> Rstats
echo "cat('MedianQV:\t', -10*log10(1-(median(dat[,1])/100)), '\n')" >> Rstats
echo "cat('Mean:\t\t', mean(dat[,1]), '\n')" >> Rstats
echo "cat('MeanQV:\t\t', -10*log10(1-(mean(dat[,1])/100)), '\n')" >> Rstats

cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) { print $0}}' |awk '{print $4}' > $PREFIX.R
Rscript Rstats
echo "***** STATS IGNORING INDELS ********************"
cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) { print $0}}' |awk '{print $(NF-1)}' > $PREFIX.R
Rscript Rstats

# now unique bacs
if [ -e goodBacs ]; then
   echo "****** STATS FOR goodBacs ONLY ***********" 
   cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) { print $0}}' |grep -f goodBacs |awk '{print $4}' > $PREFIX.R
   Rscript Rstats
fi
echo "**********************************************"
