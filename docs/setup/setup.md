# Pipeline setup

In this section we provide detailed instructions on how to setup cosigt.

## Prerequisites

Cosigt comes with a pipeline built on version 7(.32.4) of the popular [snakemake](https://snakemake.github.io/#) workflow language. One can find detailed instructions on how to install snakemake [in the official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). A minimal installation snippet based on [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) is shown below.

::: warning
At the time of writing, snakemake migrated to version 8 and introduced some changes we have not explored yet - therefore version 8 is not supported at the moment.
:::

```bash
# install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# create snakemake + apptainer environment
micromamba create \
    -n smk7324app132 \
    bioconda::snakemake=7.32.4 \
    conda-forge::apptainer=1.3.2 \
    conda-forge::cookiecutter=2.6.0 \
    conda-forge::gdown #only used to download the test set
# activate
micromamba activate smk7324app132
# confirm snakemake has been installed
snakemake --version
# 7.32.4
singularity --version
# apptainer version 1.3.2
```

::: tip
The code above should also install `PyYAML` in the environment. If this is not the case, this can be installed for instance with `pip` (`pip install PyYAML`)
:::

Cosigt pipeline combines a large number of tools across the different [branches](/docs/introduction/introduction.md), which may cause troubles at the setup stage for some users. Therefore, we find it more convenient and reproducible to provide [docker](https://www.docker.com/) containers with precompiled binaries exposing all the required softwares. These containers are internally managed by the pipeline, given a working [Apptainer](https://apptainer.org/) - formerly Singularity - installation available. For the same reason, the pipeline can be also run using dedicated [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environments, which are available for each rule of the workflow. For users not using singularity nor conda, [tools](#tools) need to be installed manually and made available in `$PATH`. 

## Configuration

Once [prerequisites](#prerequisites) are met, cosigt pipeline can be cloned. 

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
```

### Job submission on HPC

For users working with high-performance computing architectures we suggest building a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to manage job submission to the cluster. For instance, a cookiecutter profile for the SLURM Workload Manager is provided [here](https://github.com/Snakemake-Profiles/slurm) and can be configured as follows.

```bash
template="gh:Snakemake-Profiles/slurm"
cookiecutter \
    --output-dir config \
    $template
```

A possible configuration is shown below. Note that other SLURM users may need to adjust the following based on their cluster specification.

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
[10/17] sbatch_defaults (): partition=cpuq #other stuff goes here, like account, qos, etc.
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

This process will result in the following folder structure for `config/slurm`:

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
The above configuration strategy allows to run cosigt pipeline using the dedicated docker containers through singularity. Users that prefer a solution based on conda need to modify answers in 2/17 and 3/17 accordingly or modify the resulting `config/slurm/config.yaml` file with `use-singularity: "False"`and `use-conda: "True"`.
:::

Profiles for other clusters can be found online (for instance [here](https://github.com/Christian-Heyer/snakemake-lsf) for the LSF workload management solution) - though they were not tested with this workflow. 

### Organization of pipeline input

We provide a [python script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py) to automatize the organization of folders and files used by the pipeline (see the [usage](/docs/usage/usage.md) section for a detailed description and [use cases](/docs/usecases/usecases.md#prepare-cosigt-workflow) for an example). It is **strongly recommended** to use that script in the setup stage, as cosigt workflow strictly depends on the resulting organization of files and differences with respect to this structure may cause the pipeline not to work properly.

### Run cosigt

The [python script](#organization-of-pipeline-input)  builds a ready-to-use `bash` (`cosigt_smk.sh`) script to run cosigt pipeline through singularity (or conda). 
Running cosigt without relying on singularity/conda requires all the necessary [tools](#tools) to be available in `$PATH` and it's done with `snakemake cosigt --cores <cores>`.

## Tools

A list of tools/versions used across all the branches of the pipeline is provided below (alphabetical order). Versions here refer to the latest version of the pipeline. Dockerfiles in [this repo](https://github.com/davidebolo1993/cosigt_containers/tree/main) we mantain have snippets that can guide users through manual installation of the majority of the tools.


| Tool            | Version/Commit          |
| --------------- | ---------------- |
| [bedtools](https://github.com/arq5x/bedtools2)  | v2.31.0 |
| [cosigt](https://github.com/davidebolo1993/cosigt)  | v0.1.4 |
| [gafpack](https://github.com/pangenome/gafpack)  |  v0.1.2 |
| [gfainject](https://github.com/AndreaGuarracino/gfainject)  |  v0.2.0 |
| [impg](https://github.com/pangenome/impg)  |  v0.2.3 |
| [wally](https://github.com/tobiasrausch/wally)  | v0.7.1 |
| [odgi](https://github.com/pangenome/odgi)  | v0.9.2 |
| [pggb](https://github.com/pangenome/pggb)  | v0.7.4 |
| [samtools](https://github.com/samtools/samtools)  | v1.21 |
| [wfmash](https://github.com/waveygang/wfmash)  | v0.14.0 |

The reads-to-assemblies alignment step is currently branch-specific.

| Tool            | Version/Commit         | Branch          |
| --------------- | ---------------- |---------------- |
| [bwa](https://github.com/lh3/bwa)  | v0.7.18 | short-reads, ancient genomes |
| [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)  | v2.2.1 | short-reads, modern genomes  |

Some calculations and plotting are done in R and use different R libraries. Lines [here](https://github.com/davidebolo1993/cosigt_containers/blob/877bd1ecac86de8c1f079d08884c434c823786c7/renv/4.3.3/Dockerfile#L21-L50) list packages we require.

## Installation test

A small installation test is available to be downloaded from [here](https://drive.google.com/file/d/1RClqLk7pObNXmTNns3cyHtCwC_i3tnLP/view?usp=sharing). This includes:
1. a small reference sequence (`genome/genome.fa`; a 2 Mbp-long piece extracted from chr22 of then human GRCh38 reference genome)
2. 5 assemblies (`assemblies/assemblies.fa`), each exhibiting a single deletion with respect to the reference, ranging from 10 to 50 Kbp
3. synthetic reads (`sample/sample.bam`), with 50% of them coming from the reference sequence and 50% coming from the assembly with a 30 Kbp deletion 
4. a bed file (`region/region.bed`) with the region of interest

```bash
#working directory: cosigt/cosigt_smk
micromamba activate smk7324app132
#easy download with gdown
gdown "https://drive.google.com/file/d/1RClqLk7pObNXmTNns3cyHtCwC_i3tnLP/view?usp=sharing" -O ../cosigt.test.pack.tar.gz
tar -xvf ../cosigt.test.pack.tar.gz -C ../cosigt.test.pack
#minimal input organization
asm=$(find ../cosigt.test.pack -name "assemblies.fa" -exec readlink -f {} \;)
echo -e "chr22\t$asm" > asm_map.tsv
#run organize.py
python workflow/scripts/organize.py \
    -a asm_map.tsv \
    -g ../cosigt.test.pack/genome/genome.fa \
    -r ../cosigt.test.pack/sample \
    -b ../cosigt.test.pack/region/region.bed \
    -o cosigt_test \
    --threads 5
#run cosigt pipeline
sh cosigt_smk.sh
#it will take some time to pull containers but it's done once and won't be repeated for future cosigt iterations
```
After having the singularity containers pulled (this is done once), the enire run should take $\sim$ 1 minute. `cosigt_test/cosigt/sample/chr22/chr22_1000000_1060000/cosigt_genotype.tsv` shows the genotype of the sample `sample` in region `chr22:1000000-1060000` in `chr22`. The genotype reported by cosigt should be:
```txt
haplotype.1                haplotype.2
del#3#30k:1000001-1029999  ref#chr22:1000001-1060000
```

Next section [usage](/docs/usage/usage.md) contains an in-depth description of the entire workflow. [Use cases](/docs/usecases/usecases.md) has some practice on how to run cosigt pipeline. 