#$ -q all.q
# -q TELOMER,all.q,UI,BIOLOGY
#$ -pe smp 4
#$ -t 1-25:1
# -t 1
#$ -r y
#$ -ckpt user
#$ -j y
# -o /dev/null

bash
module load openmpi
module load gcc

PATH=$PATH:$HOME/bin:$PATH:$HOME/exp_shared
export PATH

#Enter sample info for TUMOR sample analysis

current_sample='trial2\/BASE'
minCluster='minConsolidate=5'
maxCluster='maxConsolidate=200'
minInsLen='minBirLength=10'
TemMatch='tempToBirPercentage=0.8'
Reference='hg38_chr\/hg38_chr.fa'

########################################################################################################
#Common commands and variables are below. DO NOT CHANGE.

common_dir='trial2\/PROJECT_DIR'
commonCluster='minConsolidate=3'
commonMxCluster='maxConsolidate=500'
commonInsLen='minBirLength=10'
commonMatch='tempToBirPercentage=0.8'
commonRef='hg38_chr\/hg38_chr.fa'

cd $HOME/exp_shared/config_Tchr/
cp config1404-chr${SGE_TASK_ID}.txt configTemp-chr${SGE_TASK_ID}.txt

sed -i -e "s/$common_dir/$current_sample/g" configTemp-chr${SGE_TASK_ID}.txt
sed -i -e "s/$commonCluster/$minCluster/g" configTemp-chr${SGE_TASK_ID}.txt
sed -i -e "s/$commonMxCluster/$maxCluster/g" configTemp-chr${SGE_TASK_ID}.txt
sed -i -e "s/$commonMatch/$TemMatch/g" configTemp-chr${SGE_TASK_ID}.txt
sed -i -e "s/$commonRef/$Reference/g" configTemp-chr${SGE_TASK_ID}.txt
sed -i -e "s/$commonInsLen/$minInsLen/g" configTemp-chr${SGE_TASK_ID}.txt

cd $HOME/exp_shared

./mmbirfinder0918 config_Tchr/configTemp-chr${SGE_TASK_ID}.txt > Log2/runTMerged${SGE_TASK_ID}.txt


cd $HOME/exp_shared/config_Tchr/
rm configTemp-chr${SGE_TASK_ID}.txt

