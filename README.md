# MMBIR-Utils
For running the MMBIR pipeline using TCGA data and snakemake


Install Mambaforge
```bash
mkdir -p ~/workspace/sm/tools

curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh 

bash Mambaforge-Linux-x86_64.sh

```
log out and log back in. Check installation with:

```bash
which conda
```
should show something like this:
```bash
~/mambaforge/bin/conda
```

Install Snakemake
```bash
mkdir -p ~/workspace/sm/snakemake-tutorial/run
cd ~/workspace/sm/snakemake-tutorial/run

wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz

tar --wildcards -xf v5.24.1.tar.gz --strip 1 "*/data" "*/environment.yaml"
```

Activate Environment
```bash
conda activate base

mamba env create --name snakemake-tutorial --file environment.yaml

conda activate snakemake-tutorial
```

Check if Snakemake works

```bash
snakemake --help
```
