# Configuration

Everything the pipeline needs is declared in `cosigt_smk/config/config.yaml` and the small tables it points to. `make init` creates them from the shipped examples; `make check` validates them.

Earlier versions of cosigt kept the three genotyping strategies on separate git branches and generated the run layout with a helper script. Both are gone: there is one branch, and the choices that used to be branches are now configuration keys.

| Old branch | Now |
| ---------- | --- |
| `master` | `read_mode: short` |
| `ancient_dna` | `read_mode: ancient` |
| `custom_alleles` | `allele_source: custom` |

## Input tables

### samples.tsv

The samples to genotype. Alignments must be BAM or CRAM against the same reference given in `reference`, and must be indexed.

```txt
sample     alignment
HG00438    /data/reads/HG00438.cram
HG01952    /data/reads/HG01952.cram
```

### regions.bed

The regions to genotype. Columns 1–3 are required; the rest are optional.

```txt
chr1  103998686  104408600  AMY1A
chr6  31972057   32055418   C4A       .  2
chr4  144109304  144140761  GYPA      chr4:144917483-144946100
```

| Column | Meaning |
| ------ | ------- |
| 1–3 | chromosome, start, end (0-based, half-open) |
| 4 | region label, reported as `gene_name` in benchmark output |
| 5 | comma-separated alternative intervals (`chrom:start-end`) whose reads are also collected |
| 6 | ploidy, default `2` |

Ploidy is per region, so a haploid locus is simply `1`. Values above 2 are genotyped normally, but the SVbyEye, svim-asm and benchmarking steps are diploid-only and are skipped for those regions with a message.

### assemblies.tsv

Used when `allele_source: assemblies`. Maps each chromosome appearing in `regions.bed` to a FASTA of assembled contigs for that chromosome. Contigs must follow the [PanSN](https://github.com/pangenome/PanSN-spec) specification, `sample#haplotype#contig`, and the FASTA must be indexed.

```txt
chromosome  fasta
chr1        /data/assemblies/chr1.fa.gz
chr6        /data/assemblies/chr6.fa.gz
```

### alleles.tsv

Used when `allele_source: custom`. Supplies the alleles for each region directly, skipping assembly alignment and liftover.

```txt
region              fasta
chr1_103998686_104408600  /data/alleles/amy.fa.gz
```

### truth_graphs.tsv

Only for `benchmark_mode: leave_all_out`. Maps each region to a pre-built odgi graph (`.og`) that contains the genotyped samples' own haplotypes. See [→ Benchmarking](/benchmarking/benchmarking.html).

## Core keys

```yaml
output: results
reference: /path/to/reference.fa.gz
samples: config/samples.tsv
regions: config/regions.bed

allele_source: assemblies      # or custom
assemblies: config/assemblies.tsv
alleles: config/alleles.tsv

read_mode: short               # short | ancient | long:<preset>
pansn_prefix: grch38#1#
tmpdir: /tmp
```

`read_mode` selects the realignment strategy:

| Value | Aligner | Notes |
| ----- | ------- | ----- |
| `short` | bwa-mem2 `mem` | modern short reads |
| `ancient` | bwa `aln`/`samse` | ancient DNA parameters |
| `long:<preset>` | minimap2 | preset is one of `map-pb`, `map-hifi`, `map-ont`, `map-iclr`, `lr:hq` |

::: warning
Unmapped-read rescue runs in `short` mode only. There, reads that failed to map to the linear reference are recovered by matching region-specific 31-mers and realigned alongside the mapped ones. It is deliberately skipped for `long:<preset>` and `ancient`, where per-read error rates make exact 31-mer matching unreliable. In those modes the `meryl` and `kfilt` resource blocks are unused.
:::

## Optional outputs

```yaml
gtf:                 # gene annotation, GTF
proteins:            # protein-coding translations, FASTA
flagger_blacklist:   # regions to exclude, BED with PanSN contig names

odgi_viz: true       # node-coverage plot per region
pangene_viz: true    # gene-order plot (needs gtf + proteins)
wally_viz: false     # read pileups against each allele
svbyeye_viz: false   # all-vs-all plot of predicted haplotypes
sv_calling: false    # svim-asm calls from predicted haplotypes
vcf: false           # merged multi-sample VCF of the genotypes
```

`gtf` and `proteins` must be given together or not at all. The VCF is produced for every `read_mode`.

## Containers

```yaml
apptainer_cleanenv: true
apptainer_extra:
```

`make check` derives the bind mounts from the paths you configured and writes them to `cosigt_smk/.cosigt/apptainer.args`, which `make run` then uses. `-e` (`--cleanenv`) is included because pggb fails without it; set `apptainer_cleanenv: false` to drop it. Anything in `apptainer_extra` is appended verbatim — this is where site-specific flags belong.

## Resources

Every rule group has a `threads` / `mem_mb` / `runtime` block. `runtime` is in minutes. Memory and runtime are multiplied by the attempt number on retry, so a job that runs out of either gets more on the next try.

```yaml
pggb:
  threads: 24
  mem_mb: 20000
  runtime: 40
  tmpdir: /tmp
  params: "-c 2 -k 101"
```

Two blocks deserve attention.

**`bwa-mem2.max_occ`** is bwa-mem2's `-c`, the seed-occurrence cap.

```yaml
bwa-mem2:
  max_occ: 500
  min_max_occ: 50
  retries: 3
```

bwa-mem2 aborts with `assert failed for seqPair size` when too many seeds share a position. An allele panel provokes this by construction, since every locus repeats once per haplotype, so seed counts scale with panel size. Lowering `-c` is the documented workaround, so each retry halves it down to `min_max_occ` — 500, 250, 125, 62 — rather than requiring a manual edit.

::: warning
A run that succeeds on a later attempt used a more restrictive `-c` than the first, and therefore discarded more high-occurrence seeds. The rule logs the value it used. Set `retries: 0` if you would rather it fail loudly and tune `max_occ` yourself.
:::

**`cluster`** exposes the tuning options of the clustering step. All values shown are the defaults.

```yaml
cluster:
  similarity_threshold: automatic   # or a fixed number
  eps_selection: quality            # or legacy
  eps_min: 0
  eps_max: 0.30
  eps_step: 0.01
  score_use_dbcv: true
  giant_cluster_fraction: 0.85
  ...
```

## Job grouping on clusters

The per-sample, per-region part of the workflow is a chain of about seven rules. Submitted individually that is `7 × samples × regions` cluster jobs, and the total is roughly:

```txt
jobs ≈ 7·S·R + 15·R + S
```

so 500 samples over 100 regions is around 350,000 submissions. Those rules are therefore assigned to a Snakemake group called `genotype`. Snakemake submits each connected component of a group as a single cluster job and runs its rules in order, cutting those submissions by 7× with no change to the DAG, to rule granularity, or to resumability. Each rule still uses its own container or conda environment inside the group job.

`--group-components genotype=N` packs N independent (sample, region) chains into one submission, giving a dial between full parallelism and a few large jobs:

```bash
make run PROFILE=slurm SMK_ARGS='--group-components genotype=8'
```

Group resources are derived automatically. Rules that run in series combine by taking the maximum, except `runtime` which is summed; rules that could run in parallel have their resources summed and `runtime` maximised. Raising `--group-components` multiplies the packed chains, so increase walltime accordingly if jobs start hitting limits.

## Reruns

The workflow keeps reusable graph, allele FASTA, aligner index, region metadata, sample unmapped-read FASTA, and final genotype outputs. Large per-sample, per-region transport files — mapped FASTA slices, filtered FASTA, realigned CRAMs, GAF, gafpack coverage — are temporary. Adding samples therefore reuses existing region graphs and indexes, and adding regions reuses chromosome-level assembly mappings, reference k-mer data, and sample-level unmapped reads.
