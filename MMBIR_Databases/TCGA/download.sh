#$ -q TELOMER
#$ -pe smp 8
#$ -N download-TCGA-OV-slice-0
#$ -cwd


#make sure these are right (also change job name)

MANIFEST="manifest_slices/TCGA-OV-WXS-BAM-manifest-slice-0.txt"
DATA_DUMP_DIR="/nfsscratch/$USER/TCGA-OV/slice-0"

USER_TOKEN="gdc-user-token.txt"

./gdc-client download -m $MANIFEST -t $USER_TOKEN -d $DATA_DUMP_DIR
