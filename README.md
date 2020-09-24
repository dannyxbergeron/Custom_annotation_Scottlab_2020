# Creation of a custom annotation
__Author__ : Danny Bergeron

__Email__ :  _<danny.bergeron@usherbrooke.ca>_

## Description
The goal of this snakemake project was to create a new custom annotation base on the newest Ensembl release `Ensembl_V101` and to include all missing genes from `RefSeq`, `snoDB`, `gtRNAdb` and our `blockbuster` data.

## Processing
#### tRNA
Since there is no tRNA in the base Ensembl annotation, gtf entries (gene, transcript, exon) for all tRNA from gtRNAdb were generated (619 tRNA). The bed12 file was taken from: `http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/Hsapi38-gene-list.html` by going to the `Download tRNAscan-SE Results`.

#### rRNA
There is already a lot of 5S and 5.8S in the base Ensembl annotation, but no 18S, 28S or 45S. A bedtools intersect (pybedtools) is used to get the rRNA present in Refseq but not in Ensembl. It turns out that the only missing rRNA genes were 3X(18S and 28S). By looking at the Refseq rRNA genes, I noticed that there was also 45S rRNA (that were not caught by bedtools intersect since they were overlapping with 5.8S genes). Since the 45S is ovelapping with the 5.8S, 18S and 28S, I decided to add the 3X45S as gene and transcript and add the 5.8S, 18S and 28S as exon (only the exon_id represent the refseq id of those genes). I took all the information from refseq for those rRNA genes and I removed the overlapping 5.8S genes from the base Ensembl annotation.

#### blockbuster_genes
For the blockbuster genes, there was 308 genes present in our last custom annotation that were not updated for the new cutoff that was used in the final article. For this new annotation I kept only the 111 genes presented in the article `"Reducing the structure bias of RNA-Seq reveals a large
number of non-annotated non-coding RNA" (PMID: 31980822)`. I did a bedtools intersect to be sure that they had the same id, and since they did, I added the entries from our previous custom annotation (with their weird names) for those 111 genes. If the final biotype in the table was not "Unknown", I replaced the gene_biotype with the final biotype (the box type was not kept for the snoRNA).

#### snoRNA
I first checked if there was overlapping snoRNA/scaRNA in the Ensembl and Refseq annotation. It turns out that there were none for Refseq and 19 for Ensembl. To decide which one to remove from the base annotation, I kept preferentially the ones that the name was starting with either 'SNOR' or 'U' (U13 for example). If they had both the same name or were both 'SNOR' or 'U' or both not, I removed the shorter one. I did a bedtools intersect to get all the snoRNA in Refseq that were not in ensembl (113 snoRNA were missing). Then I pooled Ensembl, refseq and blockbuster snoRNA/scaRNA and did a bedtools intersect with snoDB to get the missing snoRNA. I manually removed 'TERC' and the snoRNA that were in scafold chromosomes (since we don't use them) from the resulting missing snoDB snoRNA. Since there were no gene_id for those genes, I used the concatenation of 'snoDB' and the 'id' column to create the gene_id. I also removed duplicates for the remaining snoRNA (with the same method mentionned previoulsy) and create the gtf entries (gene, transcript and exon) for the Refseq and snoDB missing snoRNA (381 snoRNA from snoDB, a lot of U3).


## Installation
For Linux users :
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Answer `yes` to `Do you wish the installer to initialize Miniconda3?`


To create the Snakemake environment used to launch Snakemake, run the following. The `conda create` command can appear to be stuck on `Solving environment`. While we are actually arguably [never going to solve the environment](https://www.ipcc.ch/sr15/chapter/spm/), the command is probably not stuck. Just be patient.

```bash
exec bash
conda config --set auto_activate_base False
conda create --name smake -c bioconda -c conda-forge snakemake=5.4.5
```

Before running Snakemake, you have to initialize the environment
```bash
conda activate smake
```

## Run
To run the workflow locally simply run the following command in the Snakemake conda environment, where `$CORES` is the number of available cores.
```bash
snakemake --use-conda --cores='$CORES'
```
