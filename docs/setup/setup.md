# Pipeline Setup

This section provides detailed instructions for setting up the cosigt pipeline on your system.

## Prerequisites

Cosigt is built on version 7 of the popular [Snakemake](https://snakemake.github.io) workflow management system. Detailed instructions for installing Snakemake are available in the [official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Below is a minimal installation script using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

::: warning
We have tested the cosigt pipeline with Snakemake version 7.32.4. At the time of writing, Snakemake has migrated to versions 8/9, introducing changes that we have not yet evaluated - therefore, versions >=8 are not currently supported.
:::

```bash
# Install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# Create environment with Snakemake, Apptainer, and other dependencies
micromamba create \
    -n smk7324app132 \
    bioconda::snakemake=7.32.4 \
    conda-forge::apptainer=1.3.2 \
    conda-forge::cookiecutter=2.6.0 \
    conda-forge::gdown #only used to download the test dataset
# Activate the environment
micromamba activate smk7324app132
# Confirm successful installation
snakemake --version
# 7.32.4
singularity --version
# apptainer version 1.3.2
```

::: tip
The code above should automatically install `PyYAML` in the environment. If this package is missing, you can install it with `pip install PyYAML`
:::

Cosigt combines numerous tools across different workflow branches, which can present setup challenges for some users. To simplify deployment and enhance reproducibility, we provide [Docker](https://www.docker.com/) containers with pre-compiled binaries for all required software. These containers are automatically managed by the pipeline when a working [Apptainer](https://apptainer.org/) (formerly Singularity) installation is available. Alternatively, the pipeline can be run using dedicated [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environments, which are provided for each rule in the workflow. For users who prefer not to use Singularity or Conda, all [required tools](#tools) must be manually installed and made available in the system `$PATH`.

## Configuration

Once the [prerequisites](#prerequisites) are in place, clone the cosigt pipeline repository:

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
```

### Job Submission on HPC

For users working with high-performance computing (HPC) clusters, we recommend creating a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to manage job submission. For example, a cookiecutter profile for the SLURM workload manager is available [here](https://github.com/Snakemake-Profiles/slurm) and can be configured as follows:

```bash
template="gh:Snakemake-Profiles/slurm"
cookiecutter \
    --output-dir config \
    $template
```

A sample configuration is shown below. Note that you may need to adjust these settings based on your specific cluster configuration:

```txt
[1/17] profile_name (slurm):
[2/17] Select use_singularity
    1 - False
    2 - True
    Choose from [1/2] (1): 2
[3/17] Select use_conda
    1 - False
    2 - True
    Choose from [1/2] (1): 1
[4/17] jobs (500):
[5/17] restart_times (0): 3
[6/17] max_status_checks_per_second (10):
[7/17] max_jobs_per_second (10):
[8/17] latency_wait (5): 30
[9/17] Select print_shell_commands
    1 - False
    2 - True
    Choose from [1/2] (1): 2
[10/17] sbatch_defaults (): partition=cpuq #other options like account, qos, etc. can be added here
[11/17] cluster_sidecar_help (Use cluster sidecar. NB! Requires snakemake >= 7.0! Enter to continue...):
[12/17] Select cluster_sidecar
    1 - yes
    2 - no
    Choose from [1/2] (1): 1
[13/17] cluster_name ():
[14/17] cluster_jobname (%r_%w):
[15/17] cluster_logpath (logs/slurm/%r/%j):
[16/17] cluster_config_help (The use of cluster-config is discouraged. Rather, set snakemake CLI options in the profile configuration file (see snakemake documentation on best practices). Enter to continue...):
[17/17] cluster_config ():
```

This process will generate the following directory structure:

```txt
config/slurm
├── CookieCutter.py
├── config.yaml
├── settings.json
├── slurm-jobscript.sh
├── slurm-sidecar.py
├── slurm-status.py
├── slurm-submit.py
└── slurm_utils.py
```

::: tip
The configuration above enables running the cosigt pipeline using Docker containers through Singularity. If you prefer a Conda-based solution, change your answers in steps 2/17 and 3/17 accordingly, or modify the resulting `config/slurm/config.yaml` file with `use-singularity: "False"` and `use-conda: "True"`.
:::

Profiles for other cluster management systems can be found online (such as [this one](https://github.com/Christian-Heyer/snakemake-lsf) for LSF), though we have not tested them with this workflow.

### Organization of Pipeline Input

We provide a [Python script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py) to automate the organization of folders and files used by the pipeline (see the [→ Usage](/usage/usage.html) section for a detailed description and [→ Use Cases](/usecases/usecases.html) for examples). It is **strongly recommended** to use this script during setup, as our workflow depends on the specific file structure it generates. Deviations from this structure may cause the pipeline to malfunction.

### Running Cosigt

The [Python script](#organization-of-pipeline-input) generates a ready-to-use Bash script (`cosigt_smk.sh`) to run the cosigt pipeline through Singularity (or Conda). To run cosigt without using Singularity or Conda, ensure all necessary [tools](#tools) are available in your `$PATH` and execute `snakemake cosigt -j <cores>`.

## Tools

Below is a list of tools and their versions used across all branches of the pipeline (in alphabetical order). These versions correspond to the latest release of the pipeline. For guidance on manual installation of most tools, refer to the Dockerfiles in [this repository](https://github.com/davidebolo1993/cosigt_containers/tree/main).


| Tool            | Version/Commit          |
| --------------- | ---------------- |
| [bedtools](https://github.com/arq5x/bedtools2)  | v2.31.0 |
| [cosigt](https://github.com/davidebolo1993/cosigt)  | v0.1.4 |
| [gafpack](https://github.com/pangenome/gafpack)  |  v0.1.2 |
| [gfainject](https://github.com/AndreaGuarracino/gfainject)  |  v0.2.0 |
| [impg](https://github.com/pangenome/impg)  |  v0.2.3 |
| [odgi](https://github.com/pangenome/odgi)  | v0.9.2 |
| [pggb](https://github.com/pangenome/pggb)  | v0.7.4 |
| [samtools](https://github.com/samtools/samtools)  | v1.21 |
| [wally](https://github.com/tobiasrausch/wally)  | v0.7.1 |
| [wfmash](https://github.com/waveygang/wfmash)  | v0.14.0 |

The reads-to-assemblies alignment step uses branch-specific tools:

| Tool            | Version/Commit         | Branch          |
| --------------- | ---------------- |---------------- |
| [bwa](https://github.com/lh3/bwa)  | v0.7.18 | short-reads, ancient genomes |
| [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)  | v2.2.1 | short-reads, modern genomes  |

Various calculations and visualizations are performed in R using multiple libraries. A complete list of required R packages can be found [here](https://github.com/davidebolo1993/cosigt_containers/blob/877bd1ecac86de8c1f079d08884c434c823786c7/renv/4.3.3/Dockerfile#L21-L50).

## Installation Test

We provide a small test dataset to verify your installation, available for download [here](https://drive.google.com/file/d/1RClqLk7pObNXmTNns3cyHtCwC_i3tnLP/view?usp=sharing). This dataset includes:

1. A small reference sequence (`genome/genome.fa`): a 2 Mbp fragment extracted from chromosome 22 of the human GRCh38 reference genome
2. Five assemblies (`assemblies/assemblies.fa`): each containing a single deletion relative to the reference, ranging from 10 to 50 Kbp
3. Synthetic reads (`sample/sample.bam`): 50% derived from the reference sequence and 50% from the assembly with a 30 Kbp deletion
4. A BED file (`region/region.bed`) defining the region of interest

To run the test:

```bash
# Working directory: cosigt/cosigt_smk
micromamba activate smk7324app132
# Download test data using gdown
gdown "https://drive.google.com/file/d/1RClqLk7pObNXmTNns3cyHtCwC_i3tnLP/view?usp=sharing" -O ../cosigt.test.pack.tar.gz
tar -xvf ../cosigt.test.pack.tar.gz -C ../cosigt.test.pack
# Create minimal input mapping
asm=$(find ../cosigt.test.pack -name "assemblies.fa" -exec readlink -f {} \;)
echo -e "chr22\t$asm" > asm_map.tsv
# Run organize.py to prepare directory structure
python workflow/scripts/organize.py \
    -a asm_map.tsv \
    -g ../cosigt.test.pack/genome/genome.fa \
    -r ../cosigt.test.pack/sample \
    -b ../cosigt.test.pack/region/region.bed \
    -o cosigt_test \
    --threads 5
# Run cosigt pipeline
sh cosigt_smk.sh
# Container pulling happens only once and won't be repeated in future runs
```

After the initial pull of Singularity containers (a one-time operation), the entire test should take approximately 1 minute to complete. The genotype results will be available in `cosigt_test/cosigt/sample/chr22/chr22_1000000_1060000/cosigt_genotype.tsv`, showing the genotype of `sample` in region `chr22:1000000-1060000`. The expected genotype reported by cosigt should be:

```txt
haplotype.1                haplotype.2
del#3#30k:1000001-1029999  ref#chr22:1000001-1060000
```