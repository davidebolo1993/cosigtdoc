# Use cases

In this section, we provide an end-to-end example on how to run [cosigt](https://github.com/davidebolo1993/cosigt). The rationale behind the individual rules is better detailed in the [usage section](/docs/usage/usage.md).

## Get the data

We construct a minimal input dataset to demonstrate cosigt usage starting from publicly available data. We start by creating the test data folder:

```bash
mkdir cosigt_test
cd cosigt_test
```

### Tools

A bunch of tools are needed to complete the sections below. The easiest way is likely to install them through [conda](https://docs.conda.io/projects/conda/en/latest/index.html)-like solutions, similarly to what is presented [here](/docs/setup/setup.md#prerequisites):

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

### Reads

Our pipeline requires sequencing reads previously aligned to a reference genome in the usual [bam](https://samtools.github.io/hts-specs/SAMv1.pdf) or [cram](https://samtools.github.io/hts-specs/CRAMv3.pdf) format. We can get a couple of these alignment files from [1000G](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).

```bash
mkdir alignments
cd alignments
#get the list of all the samples
wget https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/1000G_2504_high_coverage/additional_698_related/1000G_698_related_high_coverage.sequence.index
#select a couple
grep -E "HG00438|HG01952" 1000G_698_related_high_coverage.sequence.index | cut -f 1 > 1000G.selected.txt
#also add their indexes
sed 's/$/.crai/' 1000G.selected.txt >> 1000G.selected.txt
#download the files
wget -i 1000G.selected.txt
#all done
cd ..
```

### Reference

A reference genome is needed to fetch region-specific reads from the original alignments (if these are in .cram format). The same reference is also used as one of the assemblies samples are compared to.
In this example, reads are aligned to [grch38](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/).

```bash
mkdir reference
cd reference
#get the actual reference
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
#and its index
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai
#all done
cd ..
```

### Assemblies

In cosigt pipeline, we start from individual contigs to construct pangenome graphs of the regions we want to genotype. Our workflow expects the contigs to be grouped by the reference chromosome they belong to. A couple of tipical strategies one can adopt to generate such chromosome-specific assemblies are detailed below and users can decide which one to choose based on the available data.

#### Starting from a chromosome-specific graph

One can extract a list of contigs for a certain chromosome starting from chromosome-specific pangenome graphs. [Year 1](https://github.com/human-pangenomics/HPP_Year1_Assemblies) of [HPRC](https://humanpangenome.org/) provides [such data](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/pggb/chroms/), for instance. We will limit our analysis to `chr6` and `chr22` in this example - since these are the chromosomes our [regions of interest](/docs/usecases/usecases.md#region-of-interest) lie on.

```bash
mkdir assemblies_a
cd assemblies_a
#get the chromosome-specific pangenome graphs
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr6.hprc-v1.0-pggb.gfa.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/pggb/chroms/chr22.hprc-v1.0-pggb.gfa.gz
#get AGC
wget https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1 -O HPRC-yr1.agc 
#subset both to HG00438,HG01952 plus additional 8 samples
for chr in chr6 chr22; do
    #all the contigs
    zgrep "^P" "$chr.hprc-v1.0-pggb.gfa.gz" | cut -f2 | sort -u > "${chr}.all_ids.txt"
    #subset to HG00438,HG01952
    grep -E "^HG00438#|^HG01952#" "${chr}.all_ids.txt" > "${chr}.subset.txt"
    #add additional 8
    grep -v -E "^HG00438#|^HG01952#" "${chr}.all_ids.txt" \
        | cut -d'#' -f1 \
        | sort -u \
        | head -n 8 \
        | while read id; do grep "^$id#" "${chr}.all_ids.txt"; done >> "${chr}.subset.txt"
    rm "${chr}.all_ids.txt"
    #extract fasta
    cat "${chr}.subset.txt" | while read f; do
        agc getctg HPRC-yr1.agc $f | bgzip >> "${chr}.subset.fasta.gz"
    done
    rm "${chr}.subset.txt"
    #index
    samtools faidx "${chr}.subset.fasta.gz"
done
cd ..
```

::: tip
This is just one of the possible strategies to subset assemblies from a chromosome-specific pangenome graph. One can also subset the initial .gfa to paths of interest with `odgi paths -i GFA -K PATHS.TXT -o SUBSET.OG` and then get the assemblies .fasta out of the subgraph with `odgi paths -i SUBSET.OG -f > SUBSET.FASTA`. 
:::

#### Starting from individual contigs

One can define which is the reference chromosome each of the sample-specific contig most-likely belong to using [mash distances](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x). With this approach, we don't need chromosome-specific graphs but one can start from individual assemblies.

```bash
mkdir assemblies_b
cd assemblies_b
#AGC
agc="../assemblies_a/HPRC-yr1.agc"
#individual ids
ids=$(cut -d "#" -f 1,2 ../assemblies_a/*fai | sed 's/#/./g' | sort | uniq)
#download correspoding fasta
for id in $ids; do
    agc getset $agc $id | bgzip -@ 8 > $id.fasta.gz
    samtools faidx  $id.fasta.gz
done
#sketch the reference
mash sketch -s 100000 -p 8 -i ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
for id in $ids; do
    # distance (each contig vs reference chromosomes)
    mash dist ../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.msh "$id.fasta.gz" -i -s 100000 -p 8 > "$id.tsv"
    # find those having chr6/chr22 as best matches
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
#index
for chr in chr6 chr22; do
    samtools faidx "${chr}.subset.fasta.gz"
done
```

### Region of interest

One last piece of information we need to provide to the pipeline is a list of regions of interest laying on the chromosomes we have assemblies for (human chromosomes 6 and 22, in this example). Examples of such regions are the [C4](https://en.wikipedia.org/wiki/Complement_component_4) and [CYP2D6](https://en.wikipedia.org/wiki/CYP2D6) loci.

```bash
mkdir regions
cd regions
echo -e "chr6\t31972057\t32055418\tC4" > roi.bed
echo -e "chr22\t42077656\t42253758\tCYP2D6" >> roi.bed
cd ..
```

## Run the pipeline

### Clone cosigt

First step is to clone cosigt from [this GitHub page](https://github.com/davidebolo1993/cosigt). One can then follow istructions available in the [setup section](/docs/setup/setup.md) to prepare a working environment.

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
micromamba activate smk7324app132
#or micromamba activate /path/to/smk7324app132
```

### Prepare cosigt workflow

Having all the data at hand, we can now setup the pipeline using the [dedicated script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py).

```bash
#map each chromosome to its assembly - required
for c in $(ls ../../assemblies_a/*.fasta.gz); do fasta=$(basename $c) && id=$(echo $fasta | cut -d "." -f 1) && echo -e "$id\t$c" >> asm_map.tsv; done
#create a map of alignment name to id - optional, but leads to nicer sample names
for s in $(ls ../../alignments/*.cram); do cram=$(basename $s) && id=$(echo $cram | cut -d "." -f 1) && echo -e "$s\t$id" >> sample_map.tsv ; done
#also add some annotations for the region of interest - optional
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz
#organize input for cosigt
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

### Run cosigt workflow

The `cosigt_smk.sh` file generated with the above command contains a ready-to-use command that will run cosigt using the docker containers provided within the pipelines for executables

```bash
sh cosigt_smk.sh
```

## Output exploration

### Genotyping results

Cosigt identifies the combination of haplotypes in the pangenome describing the structure of the sequenced sample (*i.e.*, the short-reads) in the region of interest. We refer to this as `haplotype deconvolution` (or `genotyping`), which is further detailed in [usage section](/docs/usage/usage.md)
Having run cosigt pipeline as described above, results are in `cosigt_test/cosigt/{sample}/{chromosome}/{region}`. Structure should look like the following for each cosigt run:

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

In each sub-directory there is a `cosigt_genotype.tsv` file, with the haplotype pair (`haplotype.1`, `haplotype.2`) assigned by cosigt to the sample (`sample.id`). The chosen haplotype combination is the one that maximise the [cosine similarity](https://en.wikipedia.org/wiki/Cosine_similarity) between the sample and the pangenome. Details on how this metric is computed is available in the [usage section](/docs/usage/usage.md). Cosine similarity scores for all the other combinations are stored in `sorted_combos.tsv`

### Structural clustering


