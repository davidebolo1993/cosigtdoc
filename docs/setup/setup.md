# Setup


## Cosigt pipeline


### Prerequisites

Cosigt comes with a pipeline built on the popular [Snakemake](https://snakemake.github.io/#) workflow language. One can find detailed installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and a minimal installation snippet based on [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) below.

```bash
#install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
#create snakemake environment. Add a couple of additional packages useful for preprocessing
micromamba create \
    -c conda-forge \
    -c bioconda \
    -c anaconda \
    -n snakemake \
    snakemake cookiecutter pyyaml pandas
#activate
micromamba activate snakemake
#confirm snakemake has been installed
snakemake --version
#8.6.0
```

Cosigt pipeline combines a large number of tools across the different [branches](/introduction/introduction.md). We provide a [docker](https://www.docker.com/) container with precompiled binaries for all the required softwares. This [container](https://hub.docker.com/r/davidebolo1993/graph_genotyper) is internally managed by our pipeline, given a working [Apptainer](https://apptainer.org/) - formerly Singularity - installation available. Alternatively, [tools](#tools) need to be installed manually (or through conda-like solutions where possible) and made available in `$PATH`. 

### Pipeline setup

Once [prerequisites](#prerequisites) are met, cosigt pipeline can be cloned. 

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
```

For those working with high-performance computing architectures we suggest building a [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to manage jobs. For instance, a cookiecutter profile for the SLURM Workload Manager is provided [here](https://github.com/Snakemake-Profiles/slurm) and can be automatically configured as follows:

```bash
mkdir -p config
template="gh:Snakemake-Profiles/slurm"
cookiecutter \
    --output-dir config \
    $template
```

A possible configuration is shown below:

```txt
[1/17] profile_name (slurm):
[2/17] Select use_singularity
    1 - False
    2 - True
    Choose from [1/2] (1): 2
[3/17] Select use_conda
    1 - False
    2 - True
    Choose from [1/2] (1): 2
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
```

which should result in the following folder structure:

```txt
tree config/slurm

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
We further provide a [python script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py) to automatize the organization of folders and files used by the pipeline (see the [usage](/usage/usage.md) section for a detailed description).

```bash
python workflow/scripts/organize.py \
    --paf <pairwise.alignment.paf> \
    --fasta <assemblies.fasta> \
    -r <reference.fasta> \
    -a <alignment.folder> \
    --roi <regions.bed> \
    --output <output.folder>
```
which should result in the following folder structure:

```txt
tree config

config/
|-- config.yaml
|-- samples.tsv
`-- slurm
    |-- config.yaml
    |-- CookieCutter.py
    |-- __pycache__
    |   |-- CookieCutter.cpython-310.pyc
    |   `-- slurm_utils.cpython-310.pyc
    |-- settings.json
    |-- slurm-jobscript.sh
    |-- slurm-sidecar.py
    |-- slurm-status.py
    |-- slurm-submit.py
    `-- slurm_utils.py

tree resources

resources/
|-- alignment
|   |-- sample1.{bam,cram} -> </path/to/sample1.{bam,cram}>
|   |-- sample1.{bam,cram}.{bai,csi,crai}  -> </path/to/sample1.{bam,cram}.{bai,csi,crai}>
|   |-- sample2.{bam,cram} -> </path/to/sample2.{bam,cram}>
|   `-- sample2.{bam,cram}.{bai,csi,crai}  -> </path/to/sample2.{bam,cram}.{bai,csi,crai}>
|-- extra
|   `-- blacklist.txt
|-- paf
|   `-- assemblies.paf -> </path/to/assemblies.paf>
|-- assemblies
|   `-- assemblies.fasta -> </path/to/assemblies.fasta>
|-- reference
|    -- reference.fasta -> </path/to/reference.fasta>
`-- regions
    |-- region1.bed -> </path/to/region1.bed>
    `-- region2.bed -> </path/to/region2.bed>
```

The previous command also allows users to specify resources used by some of the different tools - those that usually requires most of the resources -, define additional flags for some of the tools, identify the path in the graph that should be used as reference (this assumes haplotypes in the variation graph to follow the [PanSN-spec](https://github.com/pangenome/PanSN-spec) and potentially bind additional folders into the different singularity containers used by the pipeline. Running the above python script will also result in a ready-to-use command to run cosigt pipeline through singularity (`snakemake.run.sh`).

### Tools

