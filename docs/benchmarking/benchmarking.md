# Benchmarking

Cosigt ships a benchmarking target that measures how well the pipeline recovers haplotypes it can be scored against. It is not part of a normal genotyping run:

```bash
make run TARGET=benchmark
```

The principle is simple. Take the haplotypes cosigt predicted, take the sample's true haplotypes, align them with [edlib](https://github.com/Martinsos/edlib), and report the result as a phred-scaled quality value. What differs between the two modes is where the truth comes from, and therefore what a good score means.

::: warning
Benchmarking is diploid-only. Regions whose ploidy is not 2 are emitted with `NA` values rather than dropped. The target is also container-only: the `compute_qv` helper has no conda package, so run it with `SOFTWARE=apptainer`.
:::

## The two modes

### leave_zero_out

```yaml
benchmark_mode: leave_zero_out
```

The genotyped samples' own haplotypes are part of the pangenome. An exact hit is therefore possible, and QV measures whether cosigt finds it. This is the natural sanity check and an **upper bound** on performance: if cosigt cannot recover a haplotype that is sitting in the graph, something is wrong upstream.

### leave_all_out

```yaml
benchmark_mode: leave_all_out
truth_graphs: config/truth_graphs.tsv
```

The graph cosigt genotypes against contains **no haplotype from any genotyped sample**, so it must reconstruct each sample from other people's haplotypes. This is the realistic setting — it is what happens for any sample not already in your panel.

You build that graph by pointing `assemblies` at a panel that excludes the genotyped samples. The truth then has to come from somewhere else, which is what `truth_graphs` supplies: a TSV mapping each region to a pre-built odgi graph that *does* contain the samples' own haplotypes.

```txt
region                  graph
chr1_103998686_104408600  /data/truth/amy.truth.og
chr6_31972057_32055418    /data/truth/c4.truth.og
```

::: tip
Both predicted and true sequences are read from that one graph, so they are guaranteed to share an orientation frame. Graphs must be odgi `.og`; convert a GFA with `odgi build -g in.gfa -o out.og`. `make check` verifies the table is present, covers every region, and that each graph exists and has the right extension.
:::

## QVfrac: separating two kinds of failure

Under `leave_all_out` a perfect answer is impossible by construction — the exact haplotype simply is not available. A low QV could therefore mean either of two very different things:

1. cosigt chose badly from what it had, or
2. the panel contained nothing close, and cosigt did as well as anyone could.

Raw QV cannot distinguish them, so each true haplotype is additionally compared against **every** haplotype the genotyping graph did offer, and the best match for each is taken:

```txt
QV_sum_best = best(hap_1_true vs panel) + best(hap_2_true vs panel)
QVfrac      = QV_sum_pred / QV_sum_best
```

`QVfrac` near 1 means cosigt picked about the best pair available to it, so a low absolute QV reflects an incomplete panel rather than a genotyping error. `QVfrac` well below 1 means there was a better answer in the graph that cosigt did not take — that one is on the genotyper.

::: warning
This scan is the expensive part of `leave_all_out`: one edlib alignment per (true haplotype, panel haplotype) pair, so cost scales with panel size × samples. Raise `benchmark.runtime` for large panels, and `benchmark.threads` to parallelise it.
:::

## Output

Per-region tables are collected into `benchmark/benchmark.qv.tsv`:

| Column | Meaning |
| ------ | ------- |
| `sample`, `region`, `gene_name` | identifiers; `gene_name` is column 4 of the regions BED |
| `hap_1_pred`, `hap_2_pred` | predicted haplotypes |
| `cluster_1_pred`, `cluster_2_pred` | their clusters, as assigned during [structural clustering](/workflow/workflow.html#structural-clustering) |
| `hap_1_true`, `hap_2_true` | the sample's own haplotypes, in the matched order |
| `QV_1_pred`, `QV_2_pred`, `QV_sum_pred` | phred-scaled quality per haplotype, and their sum |
| `error_rate_1_pred`, `error_rate_2_pred` | edit distance over alignment length |
| `hap_1_best`, `hap_2_best` | best-matching haplotype available in the graph (`leave_all_out`) |
| `QV_1_best`, `QV_2_best`, `QV_sum_best` | the QV those would have achieved |
| `QVfrac` | `QV_sum_pred / QV_sum_best` |

Both assignments of the two predictions to the two truths are scored and the better one reported, so haplotype order does not matter. Rows are emitted with `NA` rather than dropped when a sample's own haplotypes are not in the truth graph exactly twice, or when a region's ploidy is not 2.

::: tip
`compute_qv` floors the edit distance at 0.5, so an exact match yields a length-dependent ceiling rather than infinity — roughly QV 36 for a 2 kb haplotype and 53 for 100 kb. Treat the top of the scale as "indistinguishable", not as an absolute.
:::

## Plot

`benchmark/benchmark.qv.png` summarises the result per gene as a stacked bar, genes ordered by the share in the best band, with the number of samples printed above each bar. What is banded depends on the mode:

- **`leave_zero_out`** bands `QV_pred` into four quality bands — `<= 17`, `17–23`, `23–33`, `> 33` — ordered by the share reaching `> 33`.
- **`leave_all_out`** bins `QVfrac` into quintiles, ordered by the share in `Q5: 0.8–1.0`.

`benchmark.max_bars_per_row` (default 30) splits the plot into rows when there are many regions.

## Interpreting the two together

Running both modes on the same cohort is more informative than either alone:

| leave_zero_out | leave_all_out | Reading |
| -------------- | ------------- | ------- |
| high QV | high QVfrac | working as intended |
| high QV | low QVfrac | the genotyper is missing better candidates present in the panel |
| low QV | — | something upstream is wrong: liftover, graph construction, or read recovery |

A high `QV_pred` with a low `QVfrac` is the most actionable combination, since it points at the deconvolution step rather than at panel composition.
