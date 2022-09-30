
ref="/Users/$USER/exp_shared/GRCh38.p14/GRCh38.p14.fa"
reads=$1
OUT=${reads%.UM_fq}

bwa mem -t 56 $ref $1 | samtools sort -@56 -m 2G -O bam -T temp.sort -o $2


