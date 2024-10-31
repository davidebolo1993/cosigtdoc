# Use cases

In this section, we provide an end-to-end example on how to run [cosigt](https://github.com/davidebolo1993/cosigt). The rationale behind single rules is better detailed in the [usage](/docs/usage/usage.md) section. 

## Get the data

Following, we construct a minimal input dataset to demonstrate cosigt usage starting from publicly available data. We start by creating the test data folder:

```bash
mkdir cosigt_test
cd cosigt_test
```

### Tools

A couple of tools is needed to complete this workflow.

```bash
mkdir tools
cd tools
#samtools
wget https://github.com/samtools/samtools/releases/download/1.21/samtools-1.21.tar.bz2
tar -jxvf samtools-1.21.tar.bz2
rm samtools-1.21.tar.bz2
cd samtools-1.21
./configure
make
cd ..
#agc binary for Linux
wget https://github.com/refresh-bio/agc/releases/download/v1.1/agc-1.1_x64-linux.tar.gz 
tar -xvzf agc-1.1_x64-linux.tar.gz
rm agc-1.1_x64-linux.tar.gz
```

### Aligned reads

Our pipeline requires sequencing reads previously aligned to a reference genome in the usual [bam](https://samtools.github.io/hts-specs/SAMv1.pdf) or [cram](https://samtools.github.io/hts-specs/CRAMv3.pdf) format. We can get a couple of these from [1000G](https://www.internationalgenome.org/data-portal/data-collection/30x-grch38).

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

### Reference genome

A reference genome is needed to fetch reads from the original alignments if in cram format. In this example, reads are aligned to [grch38](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/).

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

::: warning
Currently, cosigt pipeline does not distinguish between cram and bam input, so the reference the original reads are aligned to must be provided regradless of the input file format. This behaviour may change in the future.
:::


### Assemblies

In cosigt pipeline, we use individual assemblies to construct local pangenome graphs of the regions we want to genotype. In this example, we retrieve assemblies spanning a certain region of interest (see [below](/docs/usecases/usecases.md#region-of-interest) - including those for the chosen individuals (see [above](/docs/usecases/usecases.md#aligned-reads)) - from [year 1 data](https://github.com/human-pangenomics/HPP_Year1_Assemblies) of the [Human Pangenome Reference Consortium](https://humanpangenome.org/).

```bash
mkdir assemblies
cd assemblies
#download year 1 assemblies in agc format
wget https://zenodo.org/record/5826274/files/HPRC-yr1.agc?download=1 -O HPRC-yr1.agc 
#get ids of HPRCv1 assemblies covering chromosome 6
wget https://raw.githubusercontent.com/davidebolo1993/cosigt/refs/heads/master/docs/chr6.hprcy1.frompggb.small.txt
#extract corresponding fasta sequence
while read -r line; do
    ../tools/agc-1.1_x64-linux/agc getctg HPRC-yr1.agc $line >> chr6.y1.fa
done < chr6.hprcy1.frompggb.small.txt
#index
../tools/samtools-1.21/samtools faidx chr6.y1.fa
#all done
cd ..
```

::: tip
`chr6.hprcy1.frompggb.small.txt` contains only 88 haplotypes spanning the region of interest. One can also use the full set of 1408 contigs spanning chromosome 6 in `chr6.hprcy1.frompggb.txt`, which are available using `wget https://raw.githubusercontent.com/davidebolo1993/cosigt/refs/heads/master/docs/chr6.hprcy1.frompggb.txt`. The remaining code above should be adjusted accordingly.
:::

### Region of interest

One last piece of information we need to provide to our pipeline is one (or more) region of interest laying on the chromosome we have assemblies for (human chromosome 6, in this example). One such region is the [C4 locus](https://en.wikipedia.org/wiki/Complement_component_4).

```bash
mkdir regions
cd regions
echo -e "grch38#chr6\t31972057\t32055418" > roi.bed
#all done
cd ..
```

::: warning
While is possible in principle to run cosigt pipeline w/ assemblies and loci belonging to many different chromosomes all at once, this has not been properly tested at the moment and is discouraged in favour of a chromosome-specific approach. This behaviour may change for a more  
:::

## Run the pipeline

### Clone cosigt

First step is to clone cosigt from [this GitHub page](https://github.com/davidebolo1993/cosigt). One can then follow istructions available in the [setup section](/docs/setup/setup.md#pipeline-setup) to prepare a working environment. Here we use the environment built as described [here](/docs/setup/setup.md#prerequisites).

```bash
git clone https://github.com/davidebolo1993/cosigt
cd cosigt/cosigt_smk
conda activate smk7324app132
```

### Prepare cosigt workflow

As better detailed in the [usage](/docs/usage/usage.md) section, we provide a [dedicated script](https://github.com/davidebolo1993/cosigt/blob/master/cosigt_smk/workflow/scripts/organize.py) to automatically structure our cosigt workflow following [Snakemake's best practices for distribution and reproducibility](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility).

```bash
#create a map of alignment name to id - optional
for s in $(ls ../../alignments/*.cram); do cram=$(basename $s) && id=$(echo $cram | cut -d "." -f 1) && echo -e "$cram\t$id" >> sample.map.tsv ; done
#also add some annotations for the region of interest - optional
wget https://raw.githubusercontent.com/davidebolo1993/cosigt/refs/heads/master/docs/chr6.c4.annotation.bed
#organize input for cosigt
python workflow/scripts/organize.py \
    -a ../../alignments \
    -r ../../reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
    --assemblies ../../assemblies/chr6.y1.fa \
    --roi ../../regions/roi.bed \
    --output cosigtC4 \
    --wfmash_tmpdir /localscratch \
    --pggb_tmpdir /localscratch \
    --samplemap sample.map.tsv \
    --annotation chr6.c4.annotation.bed
```

### Run cosigt workflow

The `snakemake.singularity.profile.run.sh` file generated with the above command contains a ready-to-use command that will run cosigt on a slurm cluster, using the docker containers provided within the pipelines for executables.

```bash
sh snakemake.singularity.profile.run.sh
```

### Initial exploration: cosigt output

Cosigt identifies the combination of haplotypes in the pangenome describing the structure of the sequenced sample in the region of interest. It also performs clustering of the haplotypes so that structurally-similar haplotypes are likely to fall in the same cluster. This is better illustrated in the [usage](/docs/usage/usage.md) section.
Having run cosigt pipeline as described above, results are in `cosigtC4/cosigt/{sample}/{region}/cosigt_genotype.tsv`.  For instance, looking at `cosigtC4/cosigt/HG00438/grch38_chr6_31972057_32055418/cosigt_genotype.tsv`, the 2 haplotypes that best describe sample `HG00438` are, as expected, `HG00438#1#JAHBCB010000040.1` and `HG00438#2#JAHBCA010000042.1`, both belonging to the same structural cluster. The cosine similarity is also reported in the last column.

### Additional exploration: pipeline output

In addition to ... continues from here ...

