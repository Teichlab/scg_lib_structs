# 10x Genomics Single Cell 3' V3

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) to see how __10x Genomics Single Cell 3' V3__ libraries are generated experimentally. This is a droplet-based method, where cells are captured inside droplets. At the same time, gel beads with barcoded oligo-dT primer containing UMIs are also captured inside the droplet. Reverse transcription happens inside the droplet. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead.

The `V3` chemistry is a significant improvement over the previous `V2` version of the kit in terms of sensitivity: many more detected genes per cell. The library structure is extremely similar to the `V2` version.

## For Your Own Experiments

Your sequencing read configuration is like this:

| Order | Read             | Cycle           | Description                                           |
|-------|------------------|-----------------|-------------------------------------------------------|
| 1     | Read 1           | 28              | `R1_001.fastq.gz`, 16 bp cell barcodes + 12 bp UMI    |
| 2     | Index 1 (__i7__) | 8 or 10         | `I1_001.fastq.gz`, Sample index                       |
| 3     | Index 2 (__i5__) | 8 or 10 or None | `I2_001.fastq.gz`, Sample index (if using dual index) |
| 4     | Read 2           | >50             | `R2_001.fastq.gz`, cDNA reads                         |

If you sequence your data via your core facility or a company, you will need to provide the sample index sequence, which is the primer (__PN-1000213/PN-2000240__) taken from the commercial kit from 10x Genomics, to them and they will demultiplex for you. You will get two `fastq` files per sample. Read 1 contains the cell barcodes and UMI and Read 2 contains the reads from cDNA.

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
28,,,,,,,,,,,
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

You can see each sample actually has four different index sequences. This is because each well from the plate __PN-1000213/PN-2000240__ actually contain four different indices for base balancing. You can also use dual index (__PN-1000215/PN-3000431__), and you should add that to the `SampleSheet.csv` if you use that. After the step is done, for each sample, you will have `R1_001.fastq.gz` and `R2_001.fastq.gz`. You are good to go from here.

## Public Data

For the purpose of demonstration, we will use the __10x Genomics Single Cell 3' V3__ data from the following paper:

```{eval-rst}
.. note::
  Mereu E, Lafzi A, Moutinho C, Ziegenhain C, McCarthy DJ, Álvarez-Varela A, Batlle E, Sagar, Grün D, Lau JK, Boutet SC, Sanada C, Ooi A, Jones RC, Kaihara K, Brampton C, Talaga Y, Sasagawa Y, Tanaka K, Hayashi T, Braeuning C, Fischer C, Sauer S, Trefzer T, Conrad C, Adiconis X, Nguyen LT, Regev A, Levin JZ, Parekh S, Janjic A, Wange LE, Bagnoli JW, Enard W, Gut M, Sandberg R, Nikaido I, Gut I, Stegle O, Heyn H (2020) **Benchmarking single-cell RNA-sequencing protocols for cell atlas projects.** *Nat Biotechnol* 38:747–755. https://doi.org/10.1038/s41587-020-0469-4
```

where the authors benchmarked quite a few different scRNA-seq methods using a standardised sample: a mixture of different human, mouse and dog cells. We are going to use the data from the __10x Genomics Single Cell 3' V3__ method. You can download the `fastq` file from [this ENA page](https://www.ebi.ac.uk/ena/browser/view/PRJNA593571?show=reads). There are two runs, but I'm just downloading the first run for the demonstration.

```console
mkdir -p mereu2020/10xV3
wget -P mereu2020/10xV3 -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10587809/SRR10587809_1.fastq.gz
wget -P mereu2020/10xV3 -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR105/009/SRR10587809/SRR10587809_2.fastq.gz
```

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the `V3` chemistry. If you have `cellranger` in your computer, you will find a file called `3M-february-2018.txt.gz` in the `lib/python/cellranger/barcodes/` directory. If you don't have `cellranger`, I have prepared the file for you:

```console
# download the whitelist
wget -P mereu2020/10xV3 wget https://teichlab.github.io/scg_lib_structs/data/3M-february-2018.txt.gz
gunzip mereu2020/10xV3/3M-february-2018.txt.gz
```

## From FastQ To Count Matrix

Now we could start the preprocessing by simply doing:

```console
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix mereu2020/star_outs/ \
     --readFilesIn mereu2020/10xV3/SRR10587809_2.fastq.gz mereu2020/10xV3/SRR10587809_1.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
     --soloCBwhitelist mereu2020/10xV3/3M-february-2018.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __10x Genomics Single Cell 3' V3__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mix_hg38_mm10/star_index`

>>> Pointing to the directory of the star index. The public data from the above paper was produced using the HCA reference sample, which consists of human PBMCs (60%), and HEK293T (6%), mouse colon (30%), NIH3T3 (3%) and dog MDCK cells (1%). Therefore, we need to use the species mixing reference genome. We also need to add the dog genome, but the dog cells only take 1% of all cells, so I did not bother in this documentation.

`--readFilesCommand zcat`

>>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix mereu2020/star_outs/`

>>> We want to keep everything organised. This directs all output files inside the `mereu2020/star_outs` directory.

`--readFilesIn mereu2020/10xV3/SRR10587809_2.fastq.gz mereu2020/10xV3/SRR10587809_1.fastq.gz`

>>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. In __10x Genomics Single Cell 3' V3__, cDNA reads come from Read 2, and the cell barcode and UMI come from Read 1. Check [the 10x Genomics Single Cell 3' V3 GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) if you are not sure.

`--soloType CB_UMI_Simple`

>>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12`

>>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist mereu2020/10xV3/3M-february-2018.txt`

>>> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell 3' V3__ is a commercial platform. The whitelist is taken from their commercial software `cellranger`.

`--soloCellFilter EmptyDrops_CR`

>>> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Forward`

>>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __10x Genomics Single Cell 3' V3__, the cDNA reads are from the Read 2 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing UMI and the Illumina Read 1 sequence. Therefore, Read 1 consists of cell barcodes and UMI comes from the first strand, complementary to the coding strand. Read 2 comes from the coding strand. Therefore, use `Forward` for __10x Genomics Single Cell 3' V3__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the 10x Genomics Single Cell 3' V3 GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) if you are not sure.

`--outSAMattributes CB UB`

>>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/mereu2020/
├── 10xV3
│   ├── 3M-february-2018.txt
│   ├── SRR10587809_1.fastq.gz
│   └── SRR10587809_2.fastq.gz
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