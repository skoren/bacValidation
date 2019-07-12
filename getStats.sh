#!/bin/bash

# sanity check inputs/prereq
if [ $# -lt 1 ]; then
   echo "Error: I need an assembly as the first parameter"
   exit
fi

haveR=`which R 2>/dev/null |wc -l`
if [ $haveR -eq 0 ]; then
   echo "Error: no R found in path, please add R to your path before running"
   exit
fi

haveBed=`which bedtools 2>/dev/null |wc -l`
if [ $haveBed -eq 0 ]; then
   echo "Error: no bedtools found in path, please add bedtools to your path before running"
   exit
fi

haveSam=`which samtools 2>/dev/null |wc -l`
if [ $haveSam -eq 0 ]; then
   echo "Error: no samtools found in path, please add samtools to your path before running"
   exit
fi

haveMinimap=`which minimap2 2>/dev/null |wc -l`
if [ $haveMinimap -eq 0 ]; then
   echo "Error: no minimap2 found in path, please add minimap2 to your path before running"
   exit
fi

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
# figure out attempted, we have a bit more tolerance than the 995 for resolved here, just 95 to account for indels in the attempted reconstruction
# first we convert the paf to a BED with tig name and BAC length included
# the groupby will give a list of mapping positions on each contig, like this:
# AC275637.1 tig00000007     0,151967,152735,57006   49467,153339,173266,139925      173266
# then we just go through and sum up the intervals and divide
cat $PREFIX.txt |awk '{print $1"\t"$6"\t"$7"\t"$7-$6"\t"$2"\t"$8}'|sort -gk1 |sort -k1,1 -k 5,5 |bedtools groupby -g 1,5 -c 2,3,6 -o collapse,collapse,max -i - |awk -F "\t" '{s=split($3, S, ","); e=split($4, E, ","); if (e != s) { print "Error: non-matching intervals"; } else { SUM=0; for (i = 1; i <= s; i++) { SUM+=E[i]-S[i]; } if (SUM/$NF > 0.95) print $1}}'|sort |uniq > $PREFIX.attempted

echo "******************* BAC SUMMARY ******************"
echo " TOTAL    : $total"
echo " BP       : $bp"

echo "************** Statistics for: $ASM ****************"
cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) print $1}'|sort |uniq |awk '{print $1}' > $PREFIX.resolved 
ATTEMPTED=`cat $PREFIX.attempted |wc -l`
FAILED=`cat $PREFIX.attempted |grep -v -f $PREFIX.resolved |wc -l`
DONE=`cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) print $1}'|sort |uniq |awk '{print $1}' |wc -l`
DONEBP=`cat $PREFIX.txt |awk '{if (($7-$6)/$8 > 0.995) print $1" "$8}'|sort |uniq |awk '{SUM+=$NF; } END { print SUM}'`
echo "$DONE $total $DONEBP $bp $FAILED" |awk '{print "BACs closed: "$1" ("$1/$2*100"); BACs attempted: "$1+$NF" %good = "($1/($1+$NF))*100"; BASES "$3" ("$3/$4*100")"}'
   
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
