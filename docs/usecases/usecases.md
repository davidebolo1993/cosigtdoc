# Use Cases

This section provides a comprehensive end-to-end example of how to run the [cosigt](https://github.com/davidebolo1993/cosigt) pipeline. The rationale behind each main step is explained in greater detail in the [usage section](/docs/usage/usage.md).

## Data Acquisition

We can create a minimal test dataset to demonstrate cosigt functionality using publicly available data. Let's start by creating a test directory:

```bash
mkdir cosigt_test
cd cosigt_test
```

### Required Tools

Several tools are necessary to complete the following sections. The most straightforward approach is to install them using [conda](https://docs.conda.io/projects/conda/en/latest/index.html) or similar package managers, as shown below (similar to the approach described in [setup documentation](/docs/setup/setup.md)):

```bash
cd cosigt_test
micromamba create \
    -p use_cases_env \
    -c bioconda \
    -c conda-forge \
    samtools \
    odgi \
    agc \
    mash
micromamba activate $PWD/use_cases_env
```

### Sequencing Reads

The cosigt pipeline requires sequencing reads previously aligned to a reference genome in standard [bam](https://samtools.github.io/hts-specs/SAMv1.pdf) or [cram](https://samtools.github.io/hts-specs/CRAMv3.pdf) format. We can obtain suitable alignment files from the [1000 Genomes Project](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38):

```bash
mkdir alignments
cd alignments
# Download the complete sample index
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index
# Select two specific samples
grep -E "HG00438|HG01952" 1000G_698_related_high_coverage.sequence.index | cut -f 1 > 1000G.selected.txt
# Add their index files
sed 's/$/.crai/' 1000G.selected.txt >> 1000G.selected.txt
# Download the selected files
wget -i 1000G.selected.txt
cd ..
```

### Reference Genome

A reference genome is required for two purposes: extracting region-specific reads from the original alignments (particularly for .cram format files) and serving as one of the assemblies against which samples are compared. In this example, reads are aligned to [GRCh38](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/):

```bash
mkdir reference
cd reference
# Download the reference genome
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
# Download its index
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
cd ..
```

### Genome Assemblies

The cosigt pipeline begins with individual contigs to construct pangenome graphs of regions targeted for genotyping. Our workflow expects contigs to be grouped by their corresponding reference chromosome. Below are two alternative strategies for generating such chromosome-specific assemblies, which users can choose based on their available data.

#### Option 1: Starting from a Chromosome-Specific Graph

One approach is to extract contig lists for specific chromosomes from existing chromosome-specific pangenome graphs. For example, [Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies) of the [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/) provides [such data](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/pggb/chroms/). We'll limit our analysis to `chr6` and `chr22` since these contain our [regions of interest](#regions-of-interest):

```bash
mkdir assemblies_a
cd assemblies_a
# Download chromosome-specific pangenome graphs
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr6.hprc-v1.0-pggb.gfa.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr22.hprc-v1.0-pggb.gfa.gz
# Download AGC (Assembly Graph Container)
wget https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1 -O HPRC-yr1.agc 
# Process each chromosome to subset HG00438, HG01952, plus 8 additional samples
for chr in chr6 chr22; do
    # Extract all contig IDs
    zgrep "^P" "$chr.hprc-v1.0-pggb.gfa.gz" | cut -f2 | sort -u > "${chr}.all_ids.txt"
    # Subset to our samples of interest
    grep -E "^HG00438#|^HG01952#" "${chr}.all_ids.txt" > "${chr}.subset.txt"
    # Add 8 additional samples for diversity
    grep -v -E "^HG00438#|^HG01952#" "${chr}.all_ids.txt" \
        | cut -d'#' -f1 \
        | sort -u \
        | head -n 8 \
        | while read id; do grep "^$id#" "${chr}.all_ids.txt"; done >> "${chr}.subset.txt"
    rm "${chr}.all_ids.txt"
    # Extract FASTA sequences
    cat "${chr}.subset.txt" | while read f; do
        agc getctg HPRC-yr1.agc $f | bgzip >> "${chr}.subset.fasta.gz"
    done
    rm "${chr}.subset.txt"
    # Index the FASTA file
    samtools faidx "${chr}.subset.fasta.gz"
done
cd ..
```

::: tip
This represents just one possible strategy for subsetting assemblies from chromosome-specific pangenome graphs. Alternatively, you can subset the initial .gfa to paths of interest using `odgi paths -i GFA -K PATHS.TXT -o SUBSET.OG` and then extract the assemblies as FASTA from the subgraph with `odgi paths -i SUBSET.OG -f > SUBSET.FASTA`. 
:::

#### Option 2: Starting from Individual Contigs

An alternative approach is to determine which reference chromosome each sample-specific contig most likely belongs to using [Mash distances](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x). This method doesn't require chromosome-specific graphs and can start directly from individual assemblies:

```bash
mkdir assemblies_b
cd assemblies_b
# Reference to AGC from previous step
agc="../assemblies_a/HPRC-yr1.agc"
# Extract individual IDs
ids=$(cut -d "#" -f 1,2 ../assemblies_a/*fai | sed 's/#/./g' | sort | uniq)
# Download corresponding FASTA files
for id in $ids; do
    agc getset $agc $id | bgzip -@ 8 > $id.fasta.gz
    samtools faidx $id.fasta.gz
done
# Create a sketch of the reference
mash sketch -s 100000 -p 8 -i ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
for id in $ids; do
    # Calculate distances (each contig vs reference chromosomes)
    mash dist ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.msh "$id.fasta.gz" -i -s 100000 -p 8 > "$id.tsv"
    # Identify contigs with chr6/chr22 as best matches
    cut -f 2 "$id.tsv" | sort -u | while read f; do
        best=$(grep -w "$f" "$id.tsv" | sort -k3n | head -1)
        chr=$(echo "$best" | awk '{print $1}')
        hap=$(echo "$best" | awk '{print $2}')
        if [[ "$chr" == "chr6" || "$chr" == "chr22" ]]; then
            samtools faidx "$id.fasta.gz" "$hap" | bgzip -@ 8 >> "${chr}.subset.fasta.gz"
        fi
    done
    rm "$id.tsv"
done
# Index the resulting chromosome-specific files
for chr in chr6 chr22; do
    samtools faidx "${chr}.subset.fasta.gz"
done
```

### Regions of Interest

The final piece of information required for the pipeline is a list of target regions located on the chromosomes for which we have assemblies (chromosomes 6 and 22 in this example). Notable examples include the [C4](https://en.wikipedia.org/wiki/Complement_component_4) and [CYP2D6](https://en.wikipedia.org/wiki/CYP2D6) loci:

```bash
mkdir regions
cd regions
echo -e "chr6\t31972057\t32055418\tC4" > roi.bed
echo -e "chr22\t42077656\t42253758\tCYP2D6" >> roi.bed
cd ..
```

## Pipeline Execution

### Clone cosigt Repository

First, clone the cosigt repository from [GitHub](https://github.com/davidebolo1993/cosigt). Then follow the instructions in the [setup section](/docs/setup/setup.md) to prepare your working environment:

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
micromamba activate smk7324app132
# Or: micromamba activate /path/to/smk7324app132
```

### Configure the Workflow

With all necessary data prepared, we can now configure the pipeline using the [dedicated setup script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py):

```bash
# Map each chromosome to its assembly file (required)
for c in $(ls ../../assemblies_a/*.fasta.gz); do fasta=$(basename $c) && id=$(echo $fasta | cut -d "." -f 1) && echo -e "$id\t$c" >> asm_map.tsv; done
# Create a mapping of alignment filenames to sample IDs (optional, but produces cleaner sample names)
for s in $(ls ../../alignments/*.cram); do cram=$(basename $s) && id=$(echo $cram | cut -d "." -f 1) && echo -e "$s\t$id" >> sample_map.tsv ; done
# Download gene annotations for the regions of interest (optional)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
# Organize input for cosigt
python workflow/scripts/organize.py \
    -a asm_map.tsv \
    -g ../../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -r ../../alignments \
    -b ../../regions/roi.bed \
    -o cosigt_test \
    --map sample_map.tsv \
    --gtf gencode.v47.annotation.gtf.gz \
    --tmp /tmp \
    --threads 8
```

### Run the Workflow

The `cosigt_smk.sh` file generated by the previous command contains a ready-to-use command that will execute cosigt using the Docker containers provided with the pipeline:

```bash
sh cosigt_smk.sh
```

## Exploring Results

### Genotyping Results

Cosigt identifies the combination of haplotypes in the pangenome that best describes the structure of the sequenced sample (from short-reads) in the region of interest. This process is called `haplotype deconvolution` (or `genotyping`), as detailed in the [usage section](/docs/usage/usage.md).

After running the cosigt pipeline as described above, results are stored in `cosigt_test/cosigt/{sample}/{chromosome}/{region}`. The directory structure should resemble the following:

```txt
tree cosigt_test/cosigt
cosigt_test/cosigt
|-- HG00438
|   |-- chr22
|   |   `-- chr22_42077656_42253758
|   |       |-- cosigt_genotype.tsv
|   |       `-- sorted_combos.tsv
|   `-- chr6
|       `-- chr6_31972057_32055418
|           |-- cosigt_genotype.tsv
|           `-- sorted_combos.tsv
`-- HG01952
    |-- chr22
    |   `-- chr22_42077656_42253758
    |       |-- cosigt_genotype.tsv
    |       `-- sorted_combos.tsv
    `-- chr6
        `-- chr6_31972057_32055418
            |-- cosigt_genotype.tsv
            `-- sorted_combos.tsv
```

Each subdirectory contains a `cosigt_genotype.tsv` file that records the haplotype pair (`haplotype.1`, `haplotype.2`) assigned by cosigt to the sample (`sample.id`). The selected haplotype combination maximizes the [cosine similarity](https://en.wikipedia.org/wiki/Cosine_similarity) between the sample and the pangenome. Detailed information about how this metric is calculated is available in the [usage section](/docs/usage/usage.md). Cosine similarity scores for all alternative combinations are stored in the `sorted_combos.tsv` file.

### Structural Clustering

XXX Add details on structural clustering here XXX