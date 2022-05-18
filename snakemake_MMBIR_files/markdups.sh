#module purge
#module load stack/2019.1
#module load openjdk/1.8.0_202-b08_gcc-9.1.0




INPUT_BAM=$1
DEDUPED_BAM_OUT=${INPUT_BAM%.bam}.deduped.bam
METRICS=${INPUT_BAM%.bam}.dedup_metrics.txt




#java -jar ~/bin/picard.jar MarkDuplicates \
#    REMOVE_DUPLICATES=true \
#    I=$INPUT_BAM \
#    O=$DEDUPED_BAM_OUT \
#    M=$METRICS






module load jdk/8u121

~/workspace/dna/tools/gatk-4.2.5.0/gatk MarkDuplicates \
  --INPUT $INPUT_BAM \
  --OUTPUT $DEDUPED_BAM_OUT \
  --METRICS_FILE $METRICS \
  --CREATE_INDEX TRUE \
  --VALIDATION_STRINGENCY LENIENT \
  --REMOVE_DUPLICATES TRUE



#Count Reads

echo $DEDUPED_BAM_OUT >> /Users/$USER/exp_shared/trial2/deduped_counts.txt
samtools view -c -@ 56 $DEDUPED_BAM_OUT >> /Users/$USER/exp_shared/trial2/deduped_counts.txt
