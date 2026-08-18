# Pipeline Setup

This section provides detailed instructions for setting up the cosigt pipeline on your system.

## Prerequisites

Cosigt is a [Snakemake](https://snakemake.github.io) workflow. Detailed instructions for installing Snakemake are available in the [official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). Below is a minimal installation using [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html).

::: warning
Cosigt now requires **Snakemake 8 or newer** and has been tested with 9.23.1. Earlier releases of the pipeline targeted Snakemake 7; that version is no longer supported, because the workflow relies on the executor-plugin architecture introduced in Snakemake 8.
:::

```bash
# Install micromamba
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# Create the environment
micromamba create \
    -n cosigt \
    -c conda-forge -c bioconda \
    snakemake=9.23.1 \
    apptainer
# Activate it
micromamba activate cosigt
# Confirm
snakemake --version
# 9.23.1
apptainer --version
# apptainer version 1.3.x
```

Cosigt combines numerous tools, which can present setup challenges. To simplify deployment and enhance reproducibility we provide [Docker](https://www.docker.com/) containers with pre-compiled binaries for all required software. These are managed automatically when a working [Apptainer](https://apptainer.org/) (formerly Singularity) installation is available, and this is the default. Alternatively the pipeline can run with [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environments, provided per rule. If you prefer neither, all [required tools](#tools) must be installed and on your `$PATH`.

### Cluster execution

Running on a cluster needs the matching Snakemake executor plugin. Install only the one your site uses:

```bash
pip install snakemake-executor-plugin-slurm
# or
pip install snakemake-executor-plugin-lsf
# or
pip install snakemake-executor-plugin-cluster-generic
```

`make check` verifies that the plugin required by your chosen profile is present and prints the exact install command if it is missing, so you do not have to work this out in advance.

::: tip
Profiles for `local`, `slurm`, `lsf` and `cluster-generic` ship with the pipeline in `cosigt_smk/profiles/`. The cookiecutter-based profile generation described in earlier versions of this documentation is no longer needed.
:::

## Getting started

Clone the repository and run the three commands that drive the whole pipeline:

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt

make init     # create config/ from the examples
make check    # validate the environment and your configuration
make run      # run the pipeline
```

### make init

Copies the example files into `cosigt_smk/config/` and writes a `.cosigt.mk` holding your default settings. Existing files are never overwritten, so it is safe to re-run.

```txt
cosigt_smk/config/
├── config.yaml         # the main configuration
├── samples.tsv         # sample -> alignment
├── regions.bed         # regions to genotype
├── assemblies.tsv      # chromosome -> assembly FASTA
├── alleles.tsv         # region -> allele FASTA   (allele_source: custom)
└── truth_graphs.tsv    # region -> truth graph    (benchmark_mode: leave_all_out)
```

Edit these before running `make check`. See the [→ Configuration](/configuration/configuration.html) section for every available key.

### make check

`check` validates everything needed for the settings you give it and fails with a specific message rather than part-way through a run:

- Snakemake is on `$PATH`.
- The executor plugin required by `PROFILE` is installed; for `local` it reports the cores it will use.
- The runtime required by `SOFTWARE` is available — `apptainer` (accepting either `apptainer` or `singularity`), `conda` (including the 24.7.1 minimum Snakemake needs), or, for `none`, that every tool the current configuration will invoke is on `$PATH`. That list is derived from your config, so it accounts for `read_mode`, `vcf`, `gtf` and the visualisation switches.
- The config file, sample table, region BED, assembly or allele tables, and every referenced file and index.
- Finally it composes the Apptainer flags this run needs.

```txt
cosigt check  (profile=local software=apptainer target=cosigt cores=32)

environment:
  ok    snakemake (9.23.1)
  ok    local execution, 32 of 32 detected cores
  ok    apptainer (apptainer version 1.3.6)

configuration and inputs:
  ok    config, sample table, regions, indexes and input files
  ok    region metadata and flagger blacklist written

apptainer flags (used automatically by 'make run'):
  -B /group/soranzo,/project/ham,/tmp -e

All checks passed. Run the pipeline with: make run
```

Those flags are the bind mounts covering every configured input and output location, collapsed to the shortest set of parent directories, plus `-e` (`--cleanenv`), which pggb requires. `make run` picks them up automatically, so bind paths no longer have to be worked out by hand.

### make run

Runs the selected target. It re-reads the flags written by `check`, so run `check` first, and again after changing any input path.

### Settings

All three commands take the same variables, either on the command line or persisted in `.cosigt.mk` by `make init`:

| Variable   | Default        | Meaning |
| ---------- | -------------- | ------- |
| `PROFILE`  | `local`        | `local`, `slurm`, `lsf`, `cluster-generic`, a path, or `none` |
| `SOFTWARE` | `apptainer`    | `apptainer`, `conda`, or `none` |
| `TARGET`   | `cosigt`       | `cosigt`, `refine`, or `benchmark` |
| `CORES`    | all detected   | passed to `--cores` |
| `SMK_ARGS` | empty          | extra Snakemake arguments, e.g. `-n` for a dry run |

```bash
make run PROFILE=slurm                     # cluster run
make run SMK_ARGS=-n                       # dry run
make run TARGET=benchmark                  # benchmarking instead of genotyping
make run SOFTWARE=conda                    # conda environments instead of containers
make run PROFILE=cluster-generic SMK_ARGS='--cluster-generic-submit-cmd "qsub"'
```

::: tip
Site-specific Apptainer flags belong in `apptainer_extra` in `config.yaml`, not in `SMK_ARGS`, which Snakemake itself interprets.
:::

## Tools

Below are the tools used across the workflow and the versions pinned by the conda environments. For manual installation, refer to the Dockerfiles in [this repository](https://github.com/davidebolo1993/cosigt_containers/tree/main).

| Tool | Version |
| ---- | ------- |
| [bcftools](https://github.com/samtools/bcftools) | 1.23.1 |
| [bedtools](https://github.com/arq5x/bedtools2) | 2.31.1 |
| [bwa](https://github.com/lh3/bwa) | 0.7.18 |
| [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) | 2.2.1 |
| [cosigt](https://github.com/davidebolo1993/cosigt) | 0.2 |
| [gafpack](https://github.com/pangenome/gafpack) | 0.1.3 |
| [gfainject](https://github.com/AndreaGuarracino/gfainject) | 0.2.0 |
| [impg](https://github.com/pangenome/impg) | 0.3.3 |
| [kfilt](https://github.com/davidebolo1993/kfilt) | 0.1.1 |
| [meryl](https://github.com/marbl/meryl) | 1.4.1 |
| [minimap2](https://github.com/lh3/minimap2) | 2.31 |
| [miniprot](https://github.com/lh3/miniprot) | 0.18 |
| [odgi](https://github.com/pangenome/odgi) | 0.9.3 |
| [pangene](https://github.com/lh3/pangene) | r231 (v1.1) |
| [panplexity](https://github.com/AndreaGuarracino/panplexity) | 0.1.1 |
| [pggb](https://github.com/pangenome/pggb) | 0.7.4 |
| [samtools](https://github.com/samtools/samtools) | 1.23.1 |
| [svim-asm](https://github.com/eldariont/svim-asm) | 1.0.3 |
| [wally](https://github.com/tobiasrausch/wally) | 0.7.1 |

Which of these are actually needed depends on your configuration: `bwa` only for `read_mode: ancient`, `kfilt`/`meryl` only for `read_mode: short`, `bcftools` only with `vcf: true`, `pangene`/`miniprot` only with `gtf`/`proteins`, and so on. `make check SOFTWARE=none` reports exactly the set your configuration requires.

Calculations and plots are performed in R. The package list is in [the renv Dockerfile](https://github.com/davidebolo1993/cosigt_containers/blob/main/renv/4.3.3/Dockerfile).

::: warning
`r-dbscan` is pinned to **1.2.2** deliberately. Releases from 1.2.3 onward crash with a segmentation fault inside `dbcv()`, which the clustering step calls, and because a segfault terminates the R process it cannot be caught. Do not relax this pin.
:::

::: tip
Two conda pins differ slightly from the containers, because the versions the containers build are not on bioconda: `gfainject` (0.2.0 vs 0.2.1) and `bedtools` (2.31.1 vs 2.31.0, since 2.31.0 cannot coexist with samtools 1.23.1). The `benchmark` target is container-only, as its `compute_qv` helper has no conda package.
:::
