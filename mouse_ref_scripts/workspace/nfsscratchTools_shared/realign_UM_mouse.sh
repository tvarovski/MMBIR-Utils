
ref="/Users/$USER/MMBSearch/exp_environment/Mus_musculus.GRCm39.dna.primary_assembly/Mus_musculus.GRCm39.dna.primary_assembly.fa"
reads=$1
OUT=${reads%.UM_fq}

bwa mem -t 56 $ref $1 | samtools sort -@56 -m 2G -O bam -T temp.sort -o $2


