# ISSAAC-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html) to see how __ISSAAC-seq__ libraries are generated experimentally. In this method, open chromatin DNA and RNA/DNA hybrid after reverse transcription were tagged by the transposase Tn5. Then single nuclei are captured either through FACS or by a droplet system with Nextera sequence on the bead. In this documentation, we are only demonstrating the droplet-based workflow with 10x Genomics Single Cell ATAC system. Soon, we will update a new version of the FACS-based workflow, which will unify both FACS- and droplet-based reagents.

## For Your Own Experiments

If you follow the protocol, you will see that __ISSAAC-seq__ leverage the __10x Genomics Single Cell ATAC__ kit for the joint detection of both gene expression (RNA) and chromatin accessibility (ATAC) from the same cell. Therefore, both the RNA and the ATAC libraries have the same configuration, and they are the same as the __10x Genomics Single Cell ATAC__ library configuration.

If you sequence your data via your core facility or a company, you will need to provide the sample and modality index sequence, which is basically the Illumina Nextera `N7xx` primer, to them and ask they to sequence the library as if they are __10x Genomics Single Cell ATAC__ libraries. They will know what to do and demultiplex for you.

If you sequence by yourself, you need to run `bcl2fastq` by yourself with a `SampleSheet.csv` in a very specific way. Here is an example of `SampleSheet.csv` of a NextSeq run with two different samples using the Illumina Nextera Indexing primers of `N701`, `N702`, `N703` and `N704`:

```text
[Header],,,,,,,,,,,
IEMFileVersion,5,,,,,,,,,,
Date,17/12/2019,,,,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,,,
Application,NextSeq FASTQ Only,,,,,,,,,,
Instrument Type,NextSeq/MiniSeq,,,,,,,,,,
Assay,AmpliSeq Library PLUS for Illumina,,,,,,,,,,
Index Adapters,AmpliSeq CD Indexes (384),,,,,,,,,,
Chemistry,Amplicon,,,,,,,,,,
,,,,,,,,,,,
[Reads],,,,,,,,,,,
151,,,,,,,,,,,
151,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample1_ATAC,,,,,,N701,TAAGGCGA,,,,
Sample1_RNA,,,,,,N702,CGTACTAG,,,,
Sample2_ATAC,,,,,,N703,AGGCAGAA,,,,
Sample2_RNA,,,,,,N704,TCCTGAGC,,,,
```

You can see the `i7` sequence in the `N7xx` primer serves as the index to discriminate both the sample (__Sample1__ vs __Sample2__) and the modality (__ATAC__ vs __RNA__).

```{eval-rst}
.. important::
  
  Make sure you understand how sequencing is done for **ISSAAC-seq** by checking `this GitHub page <https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html>`_. For convenience, we always sequence the RNA library and the ATAC library in the same lane together. Very often, we also mix them with other **10x Genomics Single Cell ATAC** libraries. For both ISSAAC-ATAC and ISSAAC-RNA, there are a total of four reads in this order:

  =====  ================  =======  ===============================================================================
  Order  Read              Cycle    Description
  =====  ================  =======  ===============================================================================
    1    Read 1            > 50     This normally yields ``R1_001.fastq.gz``, RNA cDNA reads or ATAC genomic insert
    2    Index 1 (**i7**)  8        This normally yields ``I1_001.fastq.gz``, Sample and modality index
    3    Index 2 (**i5**)  16       This normally yields ``I2_001.fastq.gz``, Cell barcodes
    4    Read 2            > 50     This normally yields ``R2_001.fastq.gz``, RNA UMI or ATAC genomic insert
  =====  ================  =======  ===============================================================================
```

Now let's look at the order of the sequencing read configuration above, as you can see, the first (`R1`), the 3rd (`I2`) and the 4th (`R2`) reads are all important for us. Therefore, we would like to get all of them for each sample based on sample and modality index, that is, the 2nd read (`I1`). To do this, you should run `bcl2fastq` in the following way:

```console
bcl2fastq --use-bases-mask=Y151,I8,Y16,Y151 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is `--use-bases-mask=Y151,I8,Y16,Y151`. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y151` at the first position indicates "use the cycle as a real read", so you will get 151-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `I8` at the second position indicates "use the cycle as an index read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read.
3. `Y16` at the third position indicates "use the cycle as a real read", so you will get 16-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read, though it is the 3rd read overall.
4. `Y151` at the fourth position indicates "use the cycle as a real read", so you will get 151-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

Therefore, you will get four fastq file per sample per modality. Using the examples above, these are the files you should get:

```bash
# files for Sample1_ATAC

Sample1_ATAC_S1_I1_001.fastq.gz # 8 bp: sample and modality index
Sample1_ATAC_S1_R1_001.fastq.gz # 151 bp: genomic insert
Sample1_ATAC_S1_R2_001.fastq.gz # 16 bp: cell barcodes
Sample1_ATAC_S1_R3_001.fastq.gz # 151 bp: genomic insert 

# files for Sample1_RNA

Sample1_RNA_S2_I1_001.fastq.gz # 8 bp: sample and modality index
Sample1_RNA_S2_R1_001.fastq.gz # 151 bp: cDNA reads
Sample1_RNA_S2_R2_001.fastq.gz # 16 bp: cell barcodes
Sample1_RNA_S2_R3_001.fastq.gz # 151 bp: The first 10 bp are UMI, the rest (poly-T) are ignored

# files for Sample2_ATAC

Sample2_ATAC_S3_I1_001.fastq.gz # 8 bp: sample and modality index
Sample2_ATAC_S3_R1_001.fastq.gz # 151 bp: genomic insert
Sample2_ATAC_S3_R2_001.fastq.gz # 16 bp: cell barcodes
Sample2_ATAC_S3_R3_001.fastq.gz # 151 bp: genomic insert

# files for Sample2_RNA

Sample2_RNA_S4_I1_001.fastq.gz # 8 bp: sample and modality index
Sample2_RNA_S4_R1_001.fastq.gz # 151 bp: cDNA reads
Sample2_RNA_S4_R2_001.fastq.gz # 16 bp: cell barcodes
Sample2_RNA_S4_R3_001.fastq.gz # 151 bp: The first 10 bp are UMI, the rest (poly-T) are ignored
```

We can safely ignore the `I1` files, but the naming here is really different from our normal usage. The `R1` files are good. The `R2` files here actually mean `I2` in our normal usage. The `R3` files here actually mean `R2` in our normal usage. Anyway, __DO NOT get confused__. You are ready to go from here.

## Public Data

The data is from the following paper from our group:

```{eval-rst}
.. note::
  Xu W, Yang W, Zhang Y, Chen Y, Zhang Q, Wang X, Song K, Jin W, Chen X (2022) **ISSAAC-seq enables sensitive and flexible multimodal profiling of chromatin accessibility and gene expression in single cells.** *bioRxiv* 2022.01.16.476488. https://doi.org/10.1101/2022.01.16.476488
```

where we developed the __ISSAAC-seq__ method for the first time. The data is deposited to ArrayExpress under the accession code [E-MTAB-11264](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11264/). There are quite a few samples there, but we will just use the __mCortex_Droplet_rep1__ sample.

```console
mkdir -p ISSAAC-seq/data
wget -P ISSAAC-seq/data -c \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847057/mCortex_rep1_Droplet_ATAC_S1_L001_I2_001.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847057/mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847057/mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847058/mCortex_rep1_Droplet_RNA_S1_L001_I2_001.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847058/mCortex_rep1_Droplet_RNA_S1_L001_R1_001.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR984/ERR9847058/mCortex_rep1_Droplet_RNA_S1_L001_R2_001.fastq.gz
```

These are the files after download:

```console
scg_prep_test/ISSAAC-seq/data/
├── mCortex_rep1_Droplet_ATAC_S1_L001_I2_001.fastq.gz
├── mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz
├── mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz
├── mCortex_rep1_Droplet_RNA_S1_L001_I2_001.fastq.gz
├── mCortex_rep1_Droplet_RNA_S1_L001_R1_001.fastq.gz
└── mCortex_rep1_Droplet_RNA_S1_L001_R2_001.fastq.gz

0 directories, 6 files
```

As you can see, we did not upload the `I1` file, because it is not important. In addition, instead of using the `R1`, `R2` and `R3` names from the `bcl2fastq`, we renamed the files during the submission. The current naming is consistent with our normal usage: `R1`, `I2` and `R2`, where `I2` is the cell barcodes.

For the RNA data, we need an extra step. To use `starsolo`, we need to stitch the cell barcode (`I2`) and the UMI (the first 10bp of `R2`) into a single `fastq`:

```bash
paste <(zcat ISSAAC-seq/data/mCortex_rep1_Droplet_RNA_S1_L001_I2_001.fastq.gz) \
      <(zcat ISSAAC-seq/data/mCortex_rep1_Droplet_RNA_S1_L001_R2_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 substr($2,1,10)} }' | \
      gzip > ISSAAC-seq/data/mCortex_rep1_Droplet_RNA_S1_L001_CB_UMI.fastq.gz
```

That yields a new fastq file `mCortex_rep1_Droplet_RNA_S1_L001_CB_UMI.fastq.gz` with 26bp in length. The first 16 bp are cell barcodes and the last 10 bp are UMI. Now you ready to go.

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the __10x Chromium Single Cell ATAC__ kit, because that is the kit used in the mCortex data set. Both the RNA library and the ATAC library are using this whitelist. If you have `cellranger-atac` in your computer, you will find a file called `737K-cratac-v1.txt.gz` in the `lib/python/atac/barcodes` directory. If you don't have `cellranger-atac`, I have prepared the file for you:

```console
# download the whitelist

wget -P ISSAAC-seq/data https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-cratac-v1.txt.gz
gunzip ISSAAC-seq/data/737K-cratac-v1.txt.gz

# reverse complement the whitelist

cat ISSAAC-seq/data/737K-cratac-v1.txt | \
    rev | tr 'ACGT' 'TGCA' > \
    ISSAAC-seq/data/737K-cratac-v1_rc.txt
```

### Explain Whitelist

You may wonder what is the reverse complementary step about. The cell barcodes in the `737K-cratac-v1.txt` are the sequences on the gel beads. These 16 bp cell barcodes are in the `i5` index location, that is, between Illumina P5 and the Nextera Read 1 sequence. It means they will be sequenced as `Index 2` (`I2`). How `i5` or `Index 2` is sequenced depends on the machine. Previously, __MiSeq__, __HiSeq 2000__, __HiSeq 2500__, __MiniSeq (Rapid)__ and __NovaSeq 6000 (v1.0)__ use the bottom strand as the template, so the index reads will be the same as the barcodes in the `737K-cratac-v1.txt`. However, more recent machines and chemistries, like __iSeq 100__, __MiniSeq (Standard)__, __NextSeq__, __HiSeq X__, __HiSeq 3000__, __HiSeq 4000__ and __NovaSeq 600 (v1.5)__, use the top strand as the template, so the index reads will be reverse complementary to the barcodes in the `737K-cratac-v1.txt`. Therefore, we need to create a reverse complementary file as the whitelist for some data.

```{eval-rst}
.. tip::
  
  Whenever you are dealing with reads that come from ``index 2``, you should:

  1. Check the sequencing machine used to generate the data
  2. Make sure you are familiar with different sequencing modes from different Illumina machines by looking at `this page <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_.
  3. Extract some sequences from ``index 2``, compare them to the whitelist or the reverse complementary to the whitelist. 
```

## From FastQ To Count Matrices

Now we are ready to map the reads to the genome using `starsolo` for the RNA library and `chromap` for the ATAC library:

```console
mkdir -p ISSAAC-seq/star_outs
mkdir -p ISSAAC-seq/chromap_outs

# process the RNA library using starsolo

STAR --runThreadN 4 \
     --genomeDir mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix ISSAAC-seq/star_outs/ \
     --readFilesIn ISSAAC-seq/data/mCortex_rep1_Droplet_RNA_S1_L001_R1_001.fastq.gz ISSAAC-seq/data/mCortex_rep1_Droplet_RNA_S1_L001_CB_UMI.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
     --soloCBwhitelist ISSAAC-seq/data/737K-cratac-v1_rc.txt \
     --clip3pNbases 116 \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate

# process the ATAC library using chromap

## map and generate the fragment file

chromap -t 4 --preset atac \
        -x mm10/chromap_index/genome.index \
        -r mm10/mm10.fa \
        -1 ISSAAC-seq/data/mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz \
        -2 ISSAAC-seq/data/mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz \
        -b ISSAAC-seq/data/mCortex_rep1_Droplet_ATAC_S1_L001_I2_001.fastq.gz \
        --barcode-whitelist ISSAAC-seq/data/737K-cratac-v1_rc.txt \
        -o ISSAAC-seq/chromap_outs/fragments.tsv

## compress and index the fragment file

bgzip ISSAAC-seq/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed ISSAAC-seq/chromap_outs/fragments.tsv.gz
```

After this stage, we are done with the RNA library. The count matrix and other useful information can be found in the `star_outs` directory. For the ATAC library, two new files `fragments.tsv.gz` and `fragments.tsv.gz.tbi` are generated. They will be useful and sometimes required for other programs to perform downstream analysis. There are still some extra work.

### Explain star and chromap

If you understand the __ISSAAC-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html), the commands above should be straightforward to understand.

#### Explain star

`--runThreadN 4`
  
> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mm10/star_index`

> Pointing to the directory of the star index. The public data we are analysing is from the cerebral cortex of an adult mouse.

`--readFilesCommand zcat`

> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix ISSAAC-seq/star_outs/`

> We want to keep everything organised. This directs all output files inside the `ISSAAC-seq/star_outs` directory.

`--readFilesIn`

> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. In __ISSAAC-seq__, cDNA reads come from Read 1, and the cell barcode and UMI come from `CB_UMI.fastq.gz` file we just prepared before . Check [the ISSAAC-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html) if you are not sure. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--soloType CB_UMI_Simple`

> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10`

> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist ISSAAC-seq/data/737K-cratac-v1_rc.txt`

> The plain text file containing all possible valid cell barcodes, one per line. For this data set, the __ISSAAC-seq__ droplet workflow used the __10x Chromium Single Cell ATAC__ kit as the single cell capture platform. This is a commercial platform. The whitelist is taken from their commercial software `cellranger-atac`. This experiment was sequenced on a __NovaSeq (v1.5)__, so we should use the reverse complementary version of the original whitelist. In other cases, you might want to use the original one.

`--clip3pNbases 116`

> We sequenced this library together with libraries from other people that requires 151 bp pair end sequencing. Therefore, we have 151 bp in the cDNA reads (`R1_001.fastq.gz`) which is unnecessarily long. The 3' of the read may contain adaptor sequences, so we just used the first 35 bp of Read 1 for the mapping, which is sufficient. This option remove 116 bp from the 3' end. Change this parameter or drop this option accordingly for your own data.

`--soloCellFilter EmptyDrops_CR`

> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Forward`

> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __ISSAAC-seq__, the cDNA reads are from the Read 1 file. During the experiment, the transposed cDNA molecules in the RNA/DNA hybrid are captured by barcoded primer containing Illumina Nextera Read 1 sequence. The reads come from the coding strand. Therefore, use `Forward` for __ISSAAC-seq__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the ISSAAC-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/ISSAAC-seq.html) if you are not sure.

`--outSAMattributes CB UB`

> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

> We want sorted `BAM` for easy handling by other programs.

#### Explain chromap

`-t 4`

> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`-x mm10/chromap_index/genome.index`

> The `chromap` index file. The public data We are analysing is from the cerebral cortex of an adult mouse.

`-r mm10/mm10.fa`

> Reference genome sequence in `fasta` format. This is basically the file which you used to create the `chromap` index file.

`-1`, `-2` and `-b`

> They are Read 1 (genomic), Read 2 (genomic) and cell barcode read, respectively. For ATAC-seq, the sequencing is usually done in pair-end mode. Therefore, you normally have two genomic reads for each genomic fragment: Read 1 and Read 2. For the reason described previously, `R1` is the genomic Read 1 and should be passed to `-1`; `R2` is the genomic Read 2 and should be passed to `-2`; `I2` is the cell barcode read and should be passed to `-b`. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--barcode-whitelist ISSAAC-seq/data/737K-cratac-v1_rc.txt`

> The plain text file containing all possible valid cell barcodes, one per line. For this data set, the __ISSAAC-seq__ droplet workflow used the __10x Chromium Single Cell ATAC__ kit as the single cell capture platform. This is a commercial platform. The whitelist is taken from their commercial software `cellranger-atac`. This experiment was sequenced on a __NovaSeq (v1.5)__, so we should use the reverse complementary version of the original whitelist. In other cases, you might want to use the original one.

`-o ISSAAC-seq/chromap_outs/fragments.tsv`

> Direct the mapped fragments to a file. The format is described in the [10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments).

### From ATAC Fragments To Reads

The fragment file is the following format:

| Column number | Meaning                               |
|---------------|---------------------------------------|
| 1             | fragment chromosome                   |
| 2             | fragment start                        |
| 3             | fragment end                          |
| 4             | cell barcode                          |
| 5             | Number of read pairs of this fragment |

It is very useful, but we often need the peak-by-cell matrix for the downstream analysis. Therefore, we need to perform a peak calling process to identify open chromatin regions. We need to convert the fragment into reads. For each fragment, we will have two reads, a forward read and a reverse read. The read length is not important, but we could generate a 50-bp read pair for each fragment.

First, we need to create a genome size file, which is a tab delimited file with only two columns. The first column is the chromosome name and the second column is the length of the chromosome in bp:

```console
# we also sort the output by chromosome name
# which will be useful later

faSize -detailed mm10/mm10.fa | \
    sort -k1,1 > mm10/mm10.chrom.sizes
```

This is the first 5 lines of `mm10/mm10.chrom.sizes`:

```
chr1	195471971
chr10	130694993
chr11	122082543
chr12	120129022
chr13	120421639
```

Now let's generate the reads from fragments:

```bash
# we use bedClip to remove reads outside the chromosome boundary
# we also remove reads mapped to the mitochondrial genome (chrM)

zcat ISSAAC-seq/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin mm10/mm10.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > ISSAAC-seq/chromap_outs/reads.bed.gz
```

Note we also sort the output reads by `sort -k1,1 -k2,2n`. In this way, the order of chromosomes in the `reads.bed.gz` is the same as that in `mm10.chrom.sizes`, which makes downstream processes easier. The output `reads.bed.gz` are the reads in `bed` format, with the 4th column holding the cell barcodes.

### Peak Calling By MACS2

Now we can use the newly generated read file for the peak calling using `MACS2`:

```console
macs2 callpeak -t ISSAAC-seq/chromap_outs/reads.bed.gz \
               -g mm -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir ISSAAC-seq/chromap_outs \
               -n aggregate
```

#### Explain MACS2

The reasons of choosing those specific parameters are a bit more complicated. I have dedicated a post for this a while ago. Please have a look at [__this post__](https://dbrg77.github.io/posts/2020-12-09-atac-seq-peak-calling-with-macs2/) if you are still confused. The following output files are particularly useful:

| File                       | Description                                                               |
|----------------------------|---------------------------------------------------------------------------|
| aggregate_peaks.narrowPeak | Open chromatin peak locations in the narrowPeak format                    |
| aggregate_peaks.xls        | More information about peaks                                              |
| aggregate_treat_pileup.bdg | Signal tracks. Can be used to generate the bigWig file for visualisation |

### Getting The Peak-By-Cell Count Matrix

Now that we have the peak and reads files, we can compute the number of reads in each peak for each cell. Then we could get the peak-by-cell count matrix. There are different ways of doing this. The following is the method I use.

#### Find Reads In Peaks Per Cell

First, we use the `aggregate_peaks.narrowPeak` file. We only need the first 4 columns (chromosome, start, end, peak ID). You can also remove the peaks that overlap [the black list regions](https://www.nature.com/articles/s41598-019-45839-z). The black list is not available for every species and every build, so I'm not doing it here. We also need to sort the peak to make sure the order of the chromosomes in the peak file is the same as that in the `mm10.chrom.sizes` and `reads.bed.gz` files. Then we could find the overlap by `bedtools`. We need to do this in a specific way to get the number of reads in each peak from each cell:

```bash
# format and sort peaks

cut -f 1-4 ISSAAC-seq/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > ISSAAC-seq/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a ISSAAC-seq/chromap_outs/aggregate_peaks_sorted.bed \
    -b ISSAAC-seq/chromap_outs/reads.bed.gz \
    -wo -sorted -g mm10/mm10.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > ISSAAC-seq/chromap_outs/peak_read_ov.tsv.gz
```

##### Explain Finding Reads In Peaks Per Cell

We start with the command before the first pipe, that is, the intersection part. If you read the manual of the `bedtools intersect`, it should be straightforward to understand. The `-wo` option will output the records in both `-a` file and `-b` file. Since the `reads.bed.gz` file has the cell barcode information at the 4th column, we would get an output with both peak and cell information for the overlap. The `-sorted -g mm10/mm10.chrom.sizes` options make the program use very little memory. Here is an example (top 5 lines) of the output of this part:

```console
chr1	3012629	3012836	aggregate_peak_54	chr1	3012465	3012665	GATGCATTGAATCGAA	.	+	36
chr1	3012629	3012836	aggregate_peak_54	chr1	3012465	3012665	GTTCTGGCTTGACTCA	.	+	36
chr1	3012629	3012836	aggregate_peak_54	chr1	3012540	3012740	CCTCCCTTGTTCGTTT	.	+	111
chr1	3012629	3012836	aggregate_peak_54	chr1	3012562	3012762	CAGGAGCCTCCTAGGC	.	-	133
chr1	3012629	3012836	aggregate_peak_54	chr1	3012577	3012777	AGGTAGCGACGTACAT	.	+	148
```

We see that the 8th column holds the cell barcode and we want to group them using `bedtools groupby`. Therefore, we need to sort by this column, that is the `sort -k8,8`. When we group by the 8th column, we are interested in how many times each peak appear per group, so we could gather the information of the peak ID (4th column). That is the `-g 8 -c 4 -o freqdesc`. The `-o freqdesc` option returns a `value:frequency` pair in descending order. Here are some records from `peak_read_ov.tsv.gz`:

```console
AAACAACGAAAACTGA	aggregate_peak_109919:2,aggregate_peak_200603:2
AAACAACGAAAAGCTA	aggregate_peak_57301:2
AAACAACGAAAAGGTT	aggregate_peak_41947:1
```

In a way, that is sort of a count matrix in an awkward format. For example:

- The first line means that in cell `AAACAACGAAAACTGA`, the peak `aggregate_peak_109919` has 2 counts and the peak `aggregate_peak_200603` has 2 counts. All the rest peaks not mentioned here have 0 counts in this cell.
- The second line means that in cell `AAACAACGAAAAGCTA`, the peak `aggregate_peak_57301` has 2 counts. All the rest peaks not mentioned here have 0 counts in this cell.

#### Output The Peak-By-Cell Matrix

At this stage, we pretty much have all the things needed. Those two files `aggregate_peaks_sorted.bed` and `peak_read_ov.tsv.gz` contain all information for a peak-by-cell count matrix. We just need a final touch to make the output in a standard format: a [market exchange format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). Since most downstream software takes the input from the __10x Genomics Single Cell ATAC__ results, we are going to generate the MEX and the associated files similar to the output from 10x Genomics.

Here, I'm using a python script for this purpose. You don't have to do this. Choose whatever works for you. The point here is to just generate similar files as the __peak-barcode matrix__ described from [the 10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

First, let's make a directory to hold the output files and generate the `peaks.bed` and `barcodes.tsv` files, which are easy to do:

```bash
# create dirctory
mkdir -p ISSAAC-seq/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 ISSAAC-seq/chromap_outs/aggregate_peaks_sorted.bed > \
    ISSAAC-seq/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat ISSAAC-seq/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    ISSAAC-seq/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
```

The slightly more difficult file to generate is `matrix.mtx`. This is the python script `generate_csc_mtx.py` for this purpose:

```python
# import helper packages
# most entries in the count matrix is 0, so we are going to use a sparse matrix
# since we need to keep updating the sparse matrix, we use lil_matrix from scipy
import sys
import gzip
from scipy.io import mmwrite
from scipy.sparse import lil_matrix

# the unique peak ID is a good renference
# generate a dictionary with peak_id : index_in_the_file
# sys.argv[1] is the 4-column bed file aggregate_peaks_sorted.bed
peak_idx = {}
with open(sys.argv[1]) as fh:
    for i, line in enumerate(fh):
        _, _, _, peak_name = line.strip().split('\t')
        peak_idx[peak_name] = i

# determine and create the dimension of the output matrix
# that is, to calculate the number of peaks and the number of barcodes
# sys.argv[2] is barcodes.tsv
num_peaks = len(peak_idx.keys())
num_cells = len(open(sys.argv[2]).readlines())
mtx = lil_matrix((num_peaks, num_cells), dtype = int)

# read the information from peak_read_ov.tsv.gz
# update the counts into the mtx object
# sys.argv[3] is peak_read_ov.tsv.gz
with gzip.open(sys.argv[3], 'rt') as fh:
    for i, line in enumerate(fh):
        col_idx = i # each column is a cell barcode
        count_info = line.strip().split('\t')[1]
        items = count_info.split(',')
        for pn_count in items:
            pn, count = pn_count.split(':')
            row_idx = peak_idx[pn] # each row is a peak
            mtx[row_idx, col_idx] = int(count)

# get a CSC sparse matrix, which is the same as the 10x Genomics matrix.mtx
mtx = mtx.tocsc()

# sys.argv[4] is the path to the output directory
mmwrite(sys.argv[4] + '/matrix.mtx', mtx, field='integer')
```

Run that script in the terminal:

```bash
python generate_csc_mtx.py \
    ISSAAC-seq/chromap_outs/aggregate_peaks_sorted.bed \
    ISSAAC-seq/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    ISSAAC-seq/chromap_outs/peak_read_ov.tsv.gz \
    ISSAAC-seq/chromap_outs/raw_peak_bc_matrix
```

After that, you should have the `matrix.mtx` in the `ISSAAC-seq/chromap_outs/raw_peak_bc_matrix` directory.

#### Cell Calling (Filter Cell Barcodes)

Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, we could use `starolo`. It is used for scRNA-seq in general, but it does have a cell calling function that takes a directory containing raw mtx and associated files, and return the filtered ones. Since `starsolo` looks for the following three files in the input directory: `matrix.mtx`, `features.tsv` and `barcodes.tsv`. Those are the output from the 10x Genomics scRNA-seq workflow. In this case, we can use `peaks.bed` as our `features.tsv`:

```console
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed ISSAAC-seq/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     ISSAAC-seq/chromap_outs/raw_peak_bc_matrix \
     ISSAAC-seq/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv ISSAAC-seq/chromap_outs/filtered_peak_bc_matrix/peaks.bed
```

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/ISSAAC-seq/
├── chromap_outs
│   ├── aggregate_control_lambda.bdg
│   ├── aggregate_peaks.narrowPeak
│   ├── aggregate_peaks_sorted.bed
│   ├── aggregate_peaks.xls
│   ├── aggregate_summits.bed
│   ├── aggregate_treat_pileup.bdg
│   ├── filtered_peak_bc_matrix
│   │   ├── barcodes.tsv
│   │   ├── features.tsv
│   │   ├── matrix.mtx
│   │   └── peaks.bed -> features.tsv
│   ├── fragments.tsv.gz
│   ├── fragments.tsv.gz.tbi
│   ├── peak_read_ov.tsv.gz
│   ├── raw_peak_bc_matrix
│   │   ├── barcodes.tsv
│   │   ├── features.tsv -> peaks.bed
│   │   ├── matrix.mtx
│   │   └── peaks.bed
│   └── reads.bed.gz
├── data
│   ├── 737K-cratac-v1_rc.txt
│   ├── 737K-cratac-v1.txt
│   ├── mCortex_rep1_Droplet_ATAC_S1_L001_I2_001.fastq.gz
│   ├── mCortex_rep1_Droplet_ATAC_S1_L001_R1_001.fastq.gz
│   ├── mCortex_rep1_Droplet_ATAC_S1_L001_R2_001.fastq.gz
│   ├── mCortex_rep1_Droplet_RNA_S1_L001_CB_UMI.fastq.gz
│   ├── mCortex_rep1_Droplet_RNA_S1_L001_I2_001.fastq.gz
│   ├── mCortex_rep1_Droplet_RNA_S1_L001_R1_001.fastq.gz
│   └── mCortex_rep1_Droplet_RNA_S1_L001_R2_001.fastq.gz
└── star_outs
    ├── Aligned.sortedByCoord.out.bam
    ├── Log.final.out
    ├── Log.out
    ├── Log.progress.out
    ├── SJ.out.tab
    └── Solo.out
        ├── Barcodes.stats
        └── Gene
            ├── Features.stats
            ├── filtered
            │   ├── barcodes.tsv
            │   ├── features.tsv
            │   └── matrix.mtx
            ├── raw
            │   ├── barcodes.tsv
            │   ├── features.tsv
            │   └── matrix.mtx
            ├── Summary.csv
            └── UMIperCellSorted.txt

9 directories, 42 files
```