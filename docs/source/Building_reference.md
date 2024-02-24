# Building Reference Genomes

First, we need to get the reference genomes. Here I'm just using human and mouse as examples. For other species, it should be similar. The thing is that you should talk to people in your community and choose a build that most people are using in your field.

## Download Reference

For the reference genomes, we are going to use the standard `mm10` build for mouse and the [analysis set sequence files](https://genome.ucsc.edu/FAQ/FAQdownloads.html#downloadAnalysis) of `hg38` for human. In addition, some benchmarking experiments use a mixture of the human and mouse cells together. Therefore, it is also useful for us to create a species mixing genome.

We download them from the [UCSC Genome Browser](https://hgdownload.soe.ucsc.edu/downloads.html). Create a directory where you want to hold the data and reference. You can call this directory `scg_prep_test`. Everything we do will be inside this directory. Do the following:

```bash
# download and extract human fasta to the hg38 directory
mkdir hg38
wget -P hg38/ -c 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz'
gunzip hg38/hg38.analysisSet.fa.gz

# download and extract mouse fasta to the mm10 directory
mkdir mm10
wget -P mm10/ -c 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz'
gunzip mm10/mm10.fa.gz

# create the species mixing genome fasta
mkdir mix_hg38_mm10
cat <(sed 's/chr/hg38_chr/' hg38/hg38.analysisSet.fa) \
    <(sed 's/chr/mm10_chr/' mm10/mm10.fa) \
    > mix_hg38_mm10/genome.fa
```

Now, we also need the annotation files. Here, we use the annotation from [GENCODE](https://www.gencodegenes.org).

```bash
# download human from https://www.gencodegenes.org/human/
# We need the Comprehensive Gene Annotation GTF file
wget -P hg38/ -c 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz'
gunzip hg38/gencode.v41.annotation.gtf.gz

# donwload mouse from https://www.gencodegenes.org/mouse/
# We need the Comprehensive Gene Annotation GTF file
wget -P mm10/ -c 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz'
gunzip mm10/gencode.vM25.annotation.gtf.gz

# create the species mixing annotation file
cat <(sed 's/^chr/hg38_chr/' hg38/gencode.v41.annotation.gtf) \
    <(sed 's/^chr/mm10_chr/' mm10/gencode.vM25.annotation.gtf) \
    > mix_hg38_mm10/annotation.gtf
```

## Making Indices

For scRNA-seq, we will use `starsolo`, so we need to build the index using `star` using 4 cores (change `--runThreadN` accordingly if using more or less cores):

```bash
# build human
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir hg38/star_index \
     --genomeFastaFiles hg38/hg38.analysisSet.fa \
     --sjdbGTFfile hg38/gencode.v41.annotation.gtf

# build mouse
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir mm10/star_index \
     --genomeFastaFiles mm10/mm10.fa \
     --sjdbGTFfile mm10/gencode.vM25.annotation.gtf

# build species mixed genome
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir mix_hg38_mm10/star_index \
     --genomeFastaFiles mix_hg38_mm10/genome.fa \
     --sjdbGTFfile mix_hg38_mm10/annotation.gtf
```

For scATAC-seq, we will use `chromap`, so we index the genome again using `chromap` with 4 cores (change `-t` accordingly if using more or less cores):

```bash
# build human
mkdir -p hg38/chromap_index
chromap -i -t 4 -r hg38/hg38.analysisSet.fa -o hg38/chromap_index/genome.index

# build mouse
mkdir -p mm10/chromap_index
chromap -i -t 4 -r mm10/mm10.fa -o mm10/chromap_index/genome.index

# build species mixed genome
mkdir -p mix_hg38_mm10/chromap_index
chromap -i -t 4 -r mix_hg38_mm10/genome.fa -o mix_hg38_mm10/chromap_index/genome.index
```

Once all the above steps are finished without errors, the `scg_prep_test` directory should look like this:

```console
scg_prep_test/
├── hg38
│   ├── chromap_index
│   │   └── genome.index
│   ├── gencode.v41.annotation.gtf
│   ├── hg38.analysisSet.fa
│   └── star_index
│       ├── chrLength.txt
│       ├── chrNameLength.txt
│       ├── chrName.txt
│       ├── chrStart.txt
│       ├── exonGeTrInfo.tab
│       ├── exonInfo.tab
│       ├── geneInfo.tab
│       ├── Genome
│       ├── genomeParameters.txt
│       ├── Log.out
│       ├── SA
│       ├── SAindex
│       ├── sjdbInfo.txt
│       ├── sjdbList.fromGTF.out.tab
│       ├── sjdbList.out.tab
│       └── transcriptInfo.tab
├── mix_hg38_mm10
│   ├── annotation.gtf
│   ├── chromap_index
│   │   └── genome.index
│   ├── genome.fa
│   └── star_index
│       ├── chrLength.txt
│       ├── chrNameLength.txt
│       ├── chrName.txt
│       ├── chrStart.txt
│       ├── exonGeTrInfo.tab
│       ├── exonInfo.tab
│       ├── geneInfo.tab
│       ├── Genome
│       ├── genomeParameters.txt
│       ├── Log.out
│       ├── SA
│       ├── SAindex
│       ├── sjdbInfo.txt
│       ├── sjdbList.fromGTF.out.tab
│       ├── sjdbList.out.tab
│       └── transcriptInfo.tab
└── mm10
    ├── chromap_index
    │   └── genome.index
    ├── gencode.vM25.annotation.gtf
    ├── mm10.fa
    └── star_index
        ├── chrLength.txt
        ├── chrNameLength.txt
        ├── chrName.txt
        ├── chrStart.txt
        ├── exonGeTrInfo.tab
        ├── exonInfo.tab
        ├── geneInfo.tab
        ├── Genome
        ├── genomeParameters.txt
        ├── Log.out
        ├── SA
        ├── SAindex
        ├── sjdbInfo.txt
        ├── sjdbList.fromGTF.out.tab
        ├── sjdbList.out.tab
        └── transcriptInfo.tab

9 directories, 57 files
```

Now, we are ready to go.