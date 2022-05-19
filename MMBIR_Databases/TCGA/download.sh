#$ -q UI,TELOMER
#$ -pe smp 16
#$ -N downloadTCGA-slice-0
#$ -cwd


#make sure these are right (also change job name)

MANIFEST="manifest_slices/TCGA-BRCA-WXS-BAM-manifest-slice-0.txt"
DATA_DUMP_DIR="/nfsscratch/twarowski/TCGA-BRCA/slice-0"

USER_TOKEN="gdc-user-token.2022-04-20T15_28_32.414Z.txt"

./gdc-client download -m $MANIFEST -t $USER_TOKEN -d $DATA_DUMP_DIR
