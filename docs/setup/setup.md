# Setup

## Cosigt pipeline

### Prerequisites

Cosigt comes with a pipeline built on version 7 of the popular [snakemake](https://snakemake.github.io/#) workflow language. One can find detailed instructions on how to install snakemake [in the official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). A minimal installation snippet based on [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) is shown below.

::: warning
At the time of writing, snakemake migrated to version 8 and introduced some changes we have not explored yet - therefore version 8 is not supported at the moment.
:::

```bash
#install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
#create snakemake + apptainer environment
micromamba create \
    -n smk7324app132 \
    bioconda::snakemake=7.32.4 \
    conda-forge::apptainer=1.3.2 \
    conda-forge::cookiecutter=2.6.0
#activate
micromamba activate smk7324app132
#confirm snakemake has been installed
snakemake --version
#7.32.4
singularity --version
#apptainer version 1.3.2
```

::: tip
The code above should also install `PyYAML`and `pandas` in the environment. If this is not the case, those can be installed for instance with `pip` (`pip install pandas==2.2.2 PyYAML==6.0.2`)
:::

Cosigt pipeline combines a large number of tools across the different [branches](/docs/introduction/introduction.md), which may cause troubles at the setup stage for some users. Therefore, we find it more convenient and reproducible to provide [docker](https://www.docker.com/) containers with precompiled binaries exposing all the required softwares. These containers are internally managed by the pipeline, given a working [Apptainer](https://apptainer.org/) - formerly Singularity - installation available. For the same reason, the pipeline can be also run using dedicated [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environments, which are available for each rule of the workflow. For users not using singularity nor conda, [tools](#tools) need to be installed manually and made available in `$PATH`. 

### Pipeline setup

Once [prerequisites](#prerequisites) are met, cosigt pipeline can be cloned. 

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
```

For those working with high-performance computing architectures we suggest building a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to manage job submission to the cluster. For instance, a cookiecutter profile for the SLURM Workload Manager is provided [here](https://github.com/Snakemake-Profiles/slurm) and can be automatically configured as follows:

```bash
template="gh:Snakemake-Profiles/slurm"
cookiecutter \
    --output-dir config \
    $template
```

A possible configuration is shown below. Note that users may need to adjust the following based on their SLURM cluster specification. 

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
[10/17] sbatch_defaults (): mem=1000 partition=cpuq time=1
[11/17] cluster_sidecar_help (Use cluster sidecar. NB! Requires snakemake >= 7.0! Enter to continue...):
[12/17] Select cluster_sidecar
    1 - yes
    2 - no
    Choose from [1/2] (1): 2
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
The above configuration strategy will run the pipeline using the dedicated docker containers through singularity. Users that prefer a solution based on conda need to modify answers in 2/17 and 3/17 accordingly or modify the resulting `config/slurm/config.yaml` file with `use-singularity: "False"`and `use-conda: "True"`.
:::

We further provide a [python script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py) to automatize the organization of folders and files used by the pipeline (see the [usage](/docs/usage/usage.md) section for a detailed description and the [use cases](/docs/usecases/usecases.md#prepare-cosigt-workflow) for an example).

```bash
python workflow/scripts/organize.py \
    --assemblies <assemblies.fasta> \
    -r <reference.fasta> \
    -a <alignment.folder> \
    --roi <regions.bed> \
    --output <output.folder>
```
The code above should result in the following folder structure in `config`

```txt
config/
|-- config.yaml
|-- samples.tsv
`-- slurm
```

and the following folder structure in `resources`

```txt
resources/
|-- alignment
|   |-- sample.{bam,cram} -> </path/to/sample.{bam,cram}>
|   |-- sample.{bam,cram}.{bai,csi,crai}  -> </path/to/sample.{bam,cram}.{bai,csi,crai}>
|   `-- ...
|-- extra
|   `-- blacklist.txt
|-- assemblies
|   `-- assemblies.fasta -> </path/to/assemblies.fasta>
|-- reference
|   `-- reference.fasta -> </path/to/reference.fasta>
`-- regions
    |-- region.bed -> </path/to/region.bed>
    `-- ...
```

This script also builds ready-to-use commands to run cosigt pipeline through singularity (`snakemake.singularity*run.sh`) or conda (`snakemake.conda*run.sh`). Running cosigt without relying on singularity/conda requires all the necessary [tools](#tools) to be available in `$PATH` and it's done with `snakemake cosigt -j <threads>`.


### Tools

A list of tools/versions used across all the branches of the pipeline is provided below (alphabetical order).

| Tool            | Version/Commit          |
| --------------- | ---------------- |
| [bedtools](https://github.com/arq5x/bedtools2)  | v2.31.0 |
| [cosigt](https://github.com/davidebolo1993/cosigt)  | v0.1.2 |
| [gafpack](https://github.com/pangenome/gafpack)  |  v0.1.0 |
| [gfainject](https://github.com/AndreaGuarracino/gfainject)  |  v0.1.0 |
| [impg](https://github.com/pangenome/impg)  |  v0.2.1 |
| [megadepth](https://github.com/ChristopherWilks/megadepth)  | v1.2.0 |
| [odgi](https://github.com/pangenome/odgi)  | v0.9.0 |
| [pggb](https://github.com/pangenome/pggb)  | v0.7.1 |
| [samtools](https://github.com/samtools/samtools)  | v1.19.2 |
| [wfmash](https://github.com/waveygang/wfmash)  | v0.14.0 |

The re-alignment step is currently branch-specific.

| Tool            | Version/Commit         | Branch          |
| --------------- | ---------------- |---------------- |
| [bwa](https://github.com/lh3/bwa)  | v0.7.18 | ancient |
| [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)  | v2.2.1 | master | 
| [minimap2](https://github.com/lh3/minimap2)  | v2.28 | long |

Structural clustering of the haplotypes and its visualisation require the following R libraries. 

| Tool            | Version/Commit          |
| --------------- | ---------------- |
| [data.table](https://cran.r-project.org/web/packages/data.table/index.html)  | v1.15.4 |
| [dendextend](https://cran.r-project.org/web/packages/dendextend/index.html)| v1.18.1 |
| [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)| v3.5.1 |
| [NbClust](https://cran.r-project.org/web/packages/NbClust/index.html)| v3.0.1 |
| [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)| v1.4.4 | 
| [rjson](https://cran.r-project.org/web/packages/rjson/index.html) | v0.2.23 |

Additional R libraries are used in companion scripts.

| Tool            | Version/Commit          |
| --------------- | ---------------- |
| [gggenes](https://cran.r-project.org/web/packages/gggenes/index.html) | v0.5.1 |
| [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html) | v1.62.0 |
| [scales](https://cran.r-project.org/web/packages/scales/index.html) | v1.3.0 |


### Installation test

A small installation test is available in `cosigt/test/cosigt.test.pack.tar.gz`. This includes:
1. a small reference sequence (`genome/genome.fa`; a 2 Mbp-long piece extracted from chr22 of then human GRCh38 reference genome)
2. 5 assemblies (`assemblies/assemblies.fa`), each exhibiting a single deletion with respect to the reference, ranging from 10 to 50 Kbp
3. synthetic reads (`sample/sample.bam`), with 50% of them coming from the reference sequence and 50% coming from the haplotype with a 30 Kbp deletion 
4. a bed file (`region/region.bed`) with the region of interest

```bash
#wd: cosigt/cosigt_smk
#conda activate smk7324app132
tar -xvf ../test/cosigt.test.pack.tar.gz -C ../test
#organize input
python workflow/scripts/organize.py \
    --assemblies ../test/cosigt.test.pack/assemblies/assemblies.fa \
    -r ../test/cosigt.test.pack/genome/genome.fa \
    -a ../test/cosigt.test.pack/sample/ \
    --roi ../test/cosigt.test.pack/region/region.bed \
    --output cosigt.test \
    --profile None \
    --threads 8
#run with singularity
sh snakemake.singularity.run.sh
#it will take some time to pull containers, but it's done once
```

`cosigt.test/cosigt/sample/ref_chr22_1000000_1060000/cosigt_genotype.tsv` shows results from the pipeline. As expected, the genotype predicted by cosigt is:
```txt
haplotype.1                haplotype.2
del#3#30k:1000001-1029999  ref#chr22:1000001-1060000
```