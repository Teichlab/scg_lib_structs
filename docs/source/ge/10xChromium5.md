# 10x Genomics Single Cell 5' V2

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html) to see how __10x Genomics Single Cell 5' V2__ libraries are generated experimentally. This is a droplet-based method, where cells are captured inside droplets. At the same time, gel beads with barcoded Templates Switch Oligo (TSO) primer containing UMIs are also captured inside the droplet. Reverse transcription happens inside the droplet, and the mRNA molecules are captured at the 5' end by the TSO. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead.

In the 5' strategy, the sequencing configuration of the `V2` chemistry is the same as the `V1` chemistry. In the 5' kit, both gene expression library and V(D)J library can be generated. We focused on the gene expression preprocessing in this documentation.

## For Your Own Experiments

You sequencing read configuration is like this:

| Order | Read             | Cycle           | Description                                                                     |
|-------|------------------|-----------------|---------------------------------------------------------------------------------|
| 1     | Read 1           | At least 26     | `R1_001.fastq.gz`, 16 bp cell barcodes + 10 bp UMI + TTTCTTATATGGG + 5' of cDNA |
| 2     | Index 1 (__i7__) | 8 or 10         | `I1_001.fastq.gz`, Sample index                                                 |
| 3     | Index 2 (__i5__) | 8 or 10 or None | `I2_001.fastq.gz`, Sample index (if any)                                        |
| 4     | Read 2           | >50             | `R2_001.fastq.gz`, cDNA reads                                                   |

Most people just do 26 cycles for Read 1, but sometimes it is good to sequence longer because you will get the transcription start sites information.

If you sequence your data via your core facility or a company, you will need to provide the sample index sequence, which is the primer (__PN-1000213__) taken from the commercial kit from 10x Genomics, to them and they will demultiplex for you. You will get two `fastq` files per sample. Read 1 contains the cell barcodes and UMI and Read 2 contains the reads from cDNA.

If you sequence by yourself, you need to run `bcl2fastq` by yourself with a `SampleSheet.csv`. Here is an example of `SampleSheet.csv` of a NextSeq run with two different samples using the indexing primers from the A1 and B1 wells, respectively:

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
26,,,,,,,,,,,
98,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample01,,,,,,SI-GA-A1_1,GGTTTACT,,,,
Sample01,,,,,,SI-GA-A1_2,CTAAACGG,,,,
Sample01,,,,,,SI-GA-A1_3,TCGGCGTC,,,,
Sample01,,,,,,SI-GA-A1_4,AACCGTAA,,,,
Sample02,,,,,,SI-GA-B1_1,GTAATCTT,,,,
Sample02,,,,,,SI-GA-B1_2,TCCGGAAG,,,,
Sample02,,,,,,SI-GA-B1_3,AGTTCGGC,,,,
Sample02,,,,,,SI-GA-B1_4,CAGCATCA,,,,
```

You can see each sample actually has four different index sequences. This is because each well from the plate __PN-1000213__ actually contain four different indices for base balancing. You can also use dual index (__PN-1000215__), and you should add that to the `SampleSheet.csv` if you use that. After the step is done, for each sample, you will have `R1_001.fastq.gz` and `R2_001.fastq.gz`. You are good to go from here.

## Public Data

For the purpose of demonstration, we will use the __10x Genomics Single Cell 5' V2 Gene Expression__ data from the following paper:

```{eval-rst}
.. note::
  Masuda K, Kornberg A, Miller J, Lin S, Suek N, Botella T, Secener KA, Bacarella AM, Cheng L, Ingham M, Rosario V, Al-Mazrou AM, Lee-Kong SA, Kiran RP, Stoeckius M, Smibert P, Portillo AD, Oberstein PE, Sims PA, Yan KS, Han A(2020) **Multiplexed single-cell analysis reveals prognostic and nonprognostic T cell types in human colorectal cancer.** *JCI Insight* 7:e154646. https://doi.org/10.1172/jci.insight.154646
```

where the authors investigated quite a lot of T cells from colorectal cancer. We are going to use the data from the __Aug13_sample1__ sample. You can download the `fastq` file from [this ArrayExpress page](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-9455/samples/?query=10x+chromium+5%27).

```console
mkdir -p masuda2022/10x5p
wget -P masuda2022/10x5p -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR466/006/ERR4667456/ERR4667456_1.fastq.gz
wget -P masuda2022/10x5p -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR466/006/ERR4667456/ERR4667456_2.fastq.gz
```

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the `V2` chemistry. If you have `cellranger` in your computer, you will find a file called `737K-august-2016.txt` in the `lib/python/cellranger/barcodes/` directory. If you don't have `cellranger`, I have prepared the file for you:

```console
# download the whitelist 
wget -P masuda2022/10x5p https://teichlab.github.io/scg_lib_structs/data/737K-august-2016.txt.gz
gunzip masuda2022/10x5p/737K-august-2016.txt.gz
```

## From FastQ To Count Matrix

Now we could start the preprocessing by simply doing:

```console
STAR --runThreadN 40 \
     --genomeDir hg38/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix masuda2022/star_outs/ \
     --readFilesIn masuda2022/10x5p/ERR4667456_2.fastq.gz masuda2022/10x5p/ERR4667456_1.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
     --soloBarcodeReadLength 0 \
     --soloCBwhitelist masuda2022/10x5p/737K-august-2016.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Reverse \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __10x Genomics Single Cell 5' V2__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir hg38/star_index`

>>> Pointing to the directory of the star index. The data is from human samples.

`--readFilesCommand zcat`

>>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix masuda2022/star_outs/`

>>> We want to keep everything organised. This directs all output files inside the `masuda2022/star_outs` directory.

`--readFilesIn masuda2022/10x5p/ERR4667456_2.fastq.gz masuda2022/10x5p/ERR4667456_1.fastq.gz`

>>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. In __10x Genomics Single Cell 5' V2__, cDNA reads come from Read 2, and the cell barcode and UMI come from Read 1. Check [the 10x Genomics Single Cell 5' V2 GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html) if you are not sure.

`--soloType CB_UMI_Simple`

>>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10`

>>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloBarcodeReadLength 0`

>>> The length of the cell barcode + UMI (Read 1) is 16 + 10 = 26 bp. Therefore, `star` by default makes sure that reads from the Read 1 file is 26-bp long. However, the data we are analysing have 28 bp in length. The last two bp of Read 1 is all `TT`. This option turns off the length check and make sure the program runs without throwing an error.

`--soloCBwhitelist masuda2022/10x5p/737K-august-2016.txt`

>>> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell 5' V2__ is a commercial platform. The whitelist is taken from their commercial software `cellranger`.

`--soloCellFilter EmptyDrops_CR`

>>> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Reverse`

>>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __10x Genomics Single Cell 5' V2__, the cDNA reads are from the Read 2 file. During the experiment, the mRNA molecules are captured at the 5' end by the TSO with an Illumina Read 1 sequence. Therefore, Read 1 consists of cell barcodes and UMI comes from the coding strand. Read 2 comes from the first strand, complementary to the coding strand. Therefore, use `Reverse` for __10x Genomics Single Cell 5' V2__ data. Check [the 10x Genomics Single Cell 5' V2 GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html) if you are not sure.

`--outSAMattributes CB UB`

>>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/masuda2022/
├── 10x5p
│   ├── 737K-august-2016.txt
│   ├── ERR4667456_1.fastq.gz
│   └── ERR4667456_2.fastq.gz
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

6 directories, 18 files
```