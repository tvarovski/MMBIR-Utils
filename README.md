# MMBIR-Utils
For running the MMBIR pipeline using TCGA data and snakemake

![Pipeline Diagram](https://github.com/tvarovski/MMBIR-Utils/blob/main/MMBIR-pipeline.png)

## Environment Set Up
To run MMBIR pipeline efficiently we utilized Snakemake to create our pipeline. To set it up, follow the instructions below.

### Install Mambaforge

To get mambaforge run following commands:
```bash
mkdir -p ~/workspace/sm/tools

curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh 

bash Mambaforge-Linux-x86_64.sh

```
Next, log out and log back in. Check installation with:

```bash
which conda
```
should show something like this:
```bash
~/mambaforge/bin/conda
```

### Install Snakemake
After getting Mambaforge, it is time to get Snakemake. Follow the commands below:

```bash
mkdir -p ~/workspace/sm/snakemake-tutorial/run
cd ~/workspace/sm/snakemake-tutorial/run

wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz

tar --wildcards -xf v5.24.1.tar.gz --strip 1 "*/data" "*/environment.yaml"
```

Activate Environment and Check if Snakemake works:
```bash
conda activate base

mamba env create --name snakemake-tutorial --file environment.yaml

conda activate snakemake-tutorial
```

```bash
snakemake --help
```


### Install GATK

```bash
mkdir ~/workspace/dna/tools/
cd ~/workspace/dna/tools
wget https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip
unzip gatk-4.2.5.0.zip
```

### Install BioAid

```bash
pip install bio-aid
```
### Get GDC-client

section in progress

