if [ $# != "3" ] || [ $1 = "-h" ] || [ $1 = "--help" ]; then
        echo "$0 input.bam output.bam output.fq"
        echo "Obtains all reads with I,D, or S in the alignment CIGAR string"
        echo ""
        exit 1
fi

samtools view -h $1 \
| awk '$1 ~ "^@" || $6 ~ "I|D|S"' \
| samtools view -b - > $2

samtools bam2fq $2 > $3
