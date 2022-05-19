if [ $# != "3" ] || [ $1 = "-h" ] || [ $1 = "--help" ]; then
        echo "$0 full.bam unmapped.bam unmapped.fq"
        echo "full.bam = bam file obtained from TCGA (big file!)"
        echo "unmapped.bam = a file name for the bam slice containing unmapped reads only"
	echo "unmapped.fq = a file name for the fq file generated from all unmapped reads"
        exit 1
fi

samtools view -f 4 -b $1 > $2
samtools bam2fq $2 > $3
