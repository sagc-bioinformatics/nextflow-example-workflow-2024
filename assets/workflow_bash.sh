#!/usr/bin/bash

# stop script when error encountered
set -e  


# ---------------------------------
# Requires:
# - minimap2
# - samtools
# - fastp 
# - multiqc

# ---------------------------------



# define variables
samples="A B C"
InDir="data/samples"
OutDir="output/aligned"
genome="data/genome.fa"
QCDir="output/QC"
LogDir="output/Logs"

threads=2

mkdir -p $OutDir
mkdir -p $QCDir
mkdir -p $LogDir

# iterate through samples
for sample in $samples
do
  echo "Processing sample ${sample}"

  # define input and output
  reads=${InDir}/${sample}.fastq
  samfile=${OutDir}/${sample}.sam
  bamfile=${OutDir}/${sample}.bam

  # run pre-mapping QC
  echo "Processing pre-mapping QC"
  fastp -i ${reads} \
    --json ${QCDir}/${sample}.fastp.json \
    --html ${QCDir}/${sample}.fastp.html
  
  # map reads
  echo "mapping reads"
  echo "minimap2 -t ${threads} -a -x sr ${genome} ${reads} > ${samfile} 2> ${LogDir}/${sample}.minimap2.log"
  minimap2 -t ${threads} -a -x sr ${genome} ${reads} > ${samfile} 2> ${LogDir}/${sample}.minimap2.log

  # convert SAM to BAM
  echo "convert SAM to BAM"
  samtools view -bh -o ${bamfile} ${samfile}

  # delete temporary data
  echo "delete temp data"
  rm ${samfile}
  
  # post-mapping QC
  samtools flagstat ${bamfile} > ${bamfile}.flagstat.txt
done


# combine QC data into a single multiQC report
multiqc -c multiqc.config -n ${QCDir}/multiqc_report ${QCDir} ${OutDir}



