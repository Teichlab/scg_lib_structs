# 10x Genomics Single Cell 3' V2

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) to see how __10x Genomics Single Cell 3' V2__ libraries are generated experimentally. This is a droplet-based method, where cells are captured inside droplets. At the same time, gel beads with barcoded oligo-dT primer containing UMIs are also captured inside the droplet. Reverse transcription happens inside the droplet. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead.

The `V2` chemistry is a significant improvement over the original `V1` version of the kit. The library structure is also quite different from the `V1` version. However, with the introduction of the `V3` chemistry, the `V2` kit will gradually become obsolete.

## For Your Own Experiments

Your sequencing read configuration is like this:

| Order | Read             | Cycle           | Description                                                       |
|-------|------------------|-----------------|-------------------------------------------------------------------|
| 1     | Read 1           | 26              | This yields `R1_001.fastq.gz`, 16 bp cell barcodes + 10 bp UMI    |
| 2     | Index 1 (__i7__) | 8 or 10         | This yields `I1_001.fastq.gz`, Sample index                       |
| 3     | Index 2 (__i5__) | 8 or 10 or None | This yields `I2_001.fastq.gz`, Sample index (if using dual index) |
| 4     | Read 2           | >50             | This yields `R2_001.fastq.gz`, cDNA reads                         |

If you sequence your data via your core facility or a company, you will need to provide the sample index sequence, which is the primer (__PN-120262/PN-220103__) taken from the commercial kit from 10x Genomics, to them and they will demultiplex for you. You will get two `fastq` files per sample. Read 1 contains the cell barcodes and UMI and Read 2 contains the reads from cDNA.

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

You can see each sample actually has four different index sequences. This is because each well from the plate __PN-120262/PN-220103__ actually contain four different indices for base balancing. Simply run `bcl2fastq` like this:

```console
bcl2fastq --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

After this, you will have `R1_001.fastq.gz` and `R2_001.fastq.gz` for each sample:

```bash
# sample01
Sample01_S1_R1_001.fastq.gz # 26 bp: cell barcode + UMI
Sample01_S1_R2_001.fastq.gz # cDNA reads

# sample02
Sample02_S2_R1_001.fastq.gz # 26 bp: cell barcode + UMI
Sample02_S2_R2_001.fastq.gz # cDNA reads
```

You are ready to go from here.

## Public Data

For the purpose of demonstration, we will use the __10x Genomics Single Cell 3' V2__ data from the following paper:

```{eval-rst}
.. note::
  Setty M, Kiseliovas V, Levine J, Gayoso A, Mazutis L, Pe'er D (2019) **Characterization of cell fate probabilities in single-cell data with Palantir.** *Nat Biotechnol* 37:451-460. https://doi.org/10.1038/s41587-019-0068-4

```

where the authors developed a computational method called `Palantir` to perform trajectory analysis on scRNA-seq data. They used the method on human bone marrow scRNA-seq to study haematopoietic differentiation. The library prepration method is __10x Genomics Single Cell 3' V2__. There are quite a few samples in this study, and you can find the raw `FASTQ` files via the accession code [PRJEB37166](https://www.ebi.ac.uk/ena/browser/view/PRJEB37166) from **ENA**. The full metadata can be obtained from the [Human Cell Atlas data portal](https://explore.data.humancellatlas.org/projects/091cf39b-01bc-42e5-9437-f419a66c8a45/project-metadata). Note that the `FASTQ` files are also available from the Human Cell Atlas website, but I found it is easier to download from the **ENA** webpage. Here, for the demonstration, we will just use the `HS_BM_P1_cells_1` sample from the donor `HS_BM_P1`. We could download them as follows:

```console
mkdir -p setty2019/data
wget -P setty2019/data -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR736/ERR7363162/Run4_SI-GA-H11_R1.fastq.gz
wget -P setty2019/data -c ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR736/ERR7363162/Run4_SI-GA-H11_R2.fastq.gz
```

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the `V2` chemistry. If you have `cellranger` in your computer, you will find a file called `737K-august-2016.txt` in the `lib/python/cellranger/barcodes/` directory. If you don't have `cellranger`, I have prepared the file for you:

```console
# download the whitelist 
wget -P setty2019/data https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz
gunzip setty2019/data/737K-august-2016.txt.gz
```

## From FastQ To Count Matrix

Now we could start the preprocessing by simply doing:

```console
STAR --runThreadN 4 \
     --genomeDir hg38/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix setty2019/star_outs/ \
     --readFilesIn setty2019/data/Run4_SI-GA-H11_R2.fastq.gz setty2019/data/Run4_SI-GA-H11_R1.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 \
     --soloBarcodeReadLength 0 \
     --soloCBwhitelist setty2019/data/737K-august-2016.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __10x Genomics Single Cell 3' V2__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir hg38/star_index`

> Pointing to the directory of the star index. The public data from the above paper was produced using CD34+ cells from bone marrow sorted by FACS from human donors. Therefore, we are using the human reference.

`--readFilesCommand zcat`

> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix setty2019/star_outs/`

> We want to keep everything organised. This directs all output files inside the `setty2019/star_outs/` directory.

`--readFilesIn setty2019/data/Run4_SI-GA-H11_R2.fastq.gz setty2019/data/Run4_SI-GA-H11_R1.fastq.gz`

> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. In __10x Genomics Single Cell 3' V2__, cDNA reads come from Read 2, and the cell barcode and UMI come from Read 1. Check [the 10x Genomics Single Cell 3' V2 GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) if you are not sure.

`--soloType CB_UMI_Simple`

> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10`

> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloBarcodeReadLength 0`

> Normally, when we specify the positions and lengths of cell barcodes and UMIs using `--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10`, the program checks if the read length matches the input. In this case, we have 16 bp cell barcodes and 10 bp UMIs. Therefore, the program expects the reads in the files containing the cell barcodes and UMIs (in this case: `R1.fastq.gz`) are 16 + 10 = 26 bp in length. However, in the study we are using, it seems the authors used 50 bp PE mode. It is possible the library was sequenced together with other libraries that require 50 bp PE. In this case, the program will throw an error and stop because 50 != 26. To prevent this from happening, we need to specify `--soloBarcodeReadLength 0`, which turns off the length check.

`--soloCBwhitelist setty2019/data/737K-august-2016.txt`

> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell 3' V2__ is a commercial platform. The whitelist is taken from their commercial software `cellranger`.

`--soloCellFilter EmptyDrops_CR`

> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Forward`

> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __10x Genomics Single Cell 3' V2__, the cDNA reads are from the Read 2 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing UMI and the Illumina Read 1 sequence. Therefore, Read 1 consists of cell barcodes and UMI comes from the first strand, complementary to the coding strand. Read 2 comes from the coding strand. Therefore, use `Forward` for __10x Genomics Single Cell 3' V2__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the 10x Genomics Single Cell 3' V2 GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html) if you are not sure.

`--outSAMattributes CB UB`

> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/setty2019/
├── data
│   ├── 737K-august-2016.txt
│   ├── Run4_SI-GA-H11_R1.fastq.gz
│   └── Run4_SI-GA-H11_R2.fastq.gz
├── filereport_read_run_PRJEB37166_tsv.txt
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