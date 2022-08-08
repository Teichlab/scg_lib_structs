# mcSCRB-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/SCRB-seq.html) to see how __SCRB-seq/mcSCRB-seq__ libraries are generated experimentally. Those are plate-based methods, where single cells are sorted into 96- or 384-well plates in a one-cell-per-well manner. Each well contains a barcoded RT primers to label mRNA from each well. After reverse transcription, cDNA from all wells are pooled together, and one single library is prepared per plate. A plate barcode is added at the library amplification stage using the Illumina Nextera Index primers. The cells can be identified using the combination of the well barcode and the plate barcode.

## For Your Own Experiments

You sequencing read configuration is like this:

| Order | Read             | Cycle     | Description                                          |
|-------|------------------|-----------|------------------------------------------------------|
| 1     | Read 1           | 16        | `R1_001.fastq.gz`, 6 bp cell barcodes + 10 bp UMI    |
| 2     | Index 1 (__i7__) | 8         | `I1_001.fastq.gz`, Plate barcode                     |
| 3     | Index 2 (__i5__) | 8 or None | `I2_001.fastq.gz`, Plate index (if using dual index) |
| 4     | Read 2           | >50       | `R2_001.fastq.gz`, cDNA reads                        |

If you sequence your data via your core facility or a company, you will need to provide the plate index sequence, which is the `N7xx` primer from the Illumina Nextera Index king, to them and they will demultiplex for you. You will get one two `fastq` files per plate. The cell barcode is basically the first 6 bp of Read 1 (`R1_001.fastq.gz`).

If you sequence by yourself, you need to run `bcl2fastq` by yourself. You have two choices:

__Option 1:__ Prepare a `SampleSheet.csv` with index information for each plate, and get the `fastq` for each plate by running `bcl2fastq`. Here is an example of `SampleSheet.csv` of a NextSeq run with five different plates:

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
16,,,,,,,,,,,
75,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Plate01,,,,,,N701,TAAGGCGA,,,,
Plate02,,,,,,N702,CGTACTAG,,,,
Plate03,,,,,,N703,AGGCAGAA,,,,
Plate04,,,,,,N704,TCCTGAGC,,,,
Plate05,,,,,,N705,GGACTCCT,,,,
```

Then, for each plate, you will get tw `fastq` files:

```
R1_001.fastq.gz
R2_001.fastq.gz
```

You can process each plate independently. See the later section on how to prepare a cell barcode `whitelist`.

__Option 2:__ I actually prefer this way: run `bcl2fastq` without a `SampleSheet.csv` like this:

```console
bcl2fastq --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

In this case, you will get three `fastq` files for each experiment:

```
Undetermined_S0_I1_001.fastq.gz
Undetermined_S0_R1_001.fastq.gz
Undetermined_S0_R2_001.fastq.gz
```

Your plate barcodes are in the `Undetermined_S0_I1_001.fastq.gz` and the well barcodes are the first 6 bp of `Undetermined_S0_R1_001.fastq.gz`. The cell barcode is basically the combination of the plate barcode and the well barcode. To run `starsolo`, we need to prepare cell barcode and UMI read files by stitching `I1` and `R1`, that is, put the sequence from `I1` in front of `R1`:

```bash
paste <(zcat Undetermined_S0_I1_001.fastq.gz) <(zcat Undetermined_S0_R1_001.fastq.gz) | \
    awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
    gzip > Undetermined_S0_CB_UMI.fastq.gz
```

The resulting `Undetermined_S0_CB_UMI.fastq.gz` file contains reads with 24 bp in length. The first 14 bp are the cell barcodes (plate barcode + well barcode), and the rest 10 bp are UMI. See the later section on how to prepare a cell barcode `whitelist`.

## Public Data

For the purpose of demonstration, we will use the __mcSCRB-seq__ data from the following paper:

```{eval-rst}
.. note::
  Mereu E, Lafzi A, Moutinho C, Ziegenhain C, McCarthy DJ, Álvarez-Varela A, Batlle E, Sagar, Grün D, Lau JK, Boutet SC, Sanada C, Ooi A, Jones RC, Kaihara K, Brampton C, Talaga Y, Sasagawa Y, Tanaka K, Hayashi T, Braeuning C, Fischer C, Sauer S, Trefzer T, Conrad C, Adiconis X, Nguyen LT, Regev A, Levin JZ, Parekh S, Janjic A, Wange LE, Bagnoli JW, Enard W, Gut M, Sandberg R, Nikaido I, Gut I, Stegle O, Heyn H (2020) **Benchmarking single-cell RNA-sequencing protocols for cell atlas projects.** *Nat Biotechnol* 38:747–755. https://doi.org/10.1038/s41587-020-0469-4
```

where the authors benchmarked quite a few different scRNA-seq methods using a standardised sample: a mixture of different human, mouse and dog cells. We are going to use the data from the __mcSCRB-seq__ method. You can download the `fastq` file from [this ENA page](https://www.ebi.ac.uk/ena/browser/view/PRJNA551755?show=reads). There are 8 plates in total, but I'm just downloading one plate for the demonstration.

```console
mkdir -p mereu2020/mcscrb-seq
wget -P mereu2020/mcscrb-seq -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR962/007/SRR9621767/SRR9621767_1.fastq.gz
wget -P mereu2020/mcscrb-seq -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR962/007/SRR9621767/SRR9621767_2.fastq.gz
```

## Prepare Whitelist

The data we get from the above paper is already demultiplexed based on plate, which is the same as the output from __Option 1__ desribed above. In this case, the cell barcode is essentially the barcode in the RT primer in the well. The barcoded oligo-dT can be download from [the protocols.io page of mcSCRB-seq](https://www.protocols.io/view/mcscrb-seq-protocol-5qpvo4d7g4o1/v2).

```console
# Download the oligo-dT file from protocols.io
# If this does not work, go to the protocols.io page and download manually
wget -O mereu2020/mcscrb-seq/mcSCRB-seq_oligodT.txt https://s3.amazonaws.com/pr-journal/uhyjf6e.txt

# the 3rd column is the well barcode
tail -n +2 mereu2020/mcscrb-seq/mcSCRB-seq_oligodT.txt | \
    cut -f 3 > mereu2020/mcscrb-seq/whitelist.txt
```

You should get 384 lines (well barcodes) in `mereu2020/mcscrb-seq/whitelist.txt`, and here I'm just showing the first five lines:

```
AAAACT
AAAATC
AAAGTT
AAATAC
AAATTG
```

If you perform mcSCRT-seq on your own and generated the `fastq` files with __Option 2__ described above, the `whitelist` should be a combination of the plate barcode and the well barcode. Just using the five plates (`N701`, `N702`, `N703`, `N704` and `N705`) described above as an example, you should put the each of the `N7xx` barcode in front of the 384 well barcode to get a total of 5 * 384 = 1920 cell barcodes. You can do this in the following way:

```bash
for pb in TAAGGCGA CGTACTAG AGGCAGAA TCCTGAGC GGACTCCT;
    do tail -n +2 mereu2020/mcscrb-seq/mcSCRB-seq_oligodT.txt | \
        awk -v PLATE_BARCODE=${pb} '{print PLATE_BARCODE $3}';
    done > mereu2020/mcscrb-seq/whitelist2.txt
```

You should get 1920 lines (well barcodes) in `mereu2020/mcscrb-seq/whitelist2.txt`, and here I'm jsut showing the first five lines:

```
TAAGGCGAAAAACT
TAAGGCGAAAAATC
TAAGGCGAAAAGTT
TAAGGCGAAAATAC
TAAGGCGAAAATTG
```

## From FastQ To Count Matrix

Now we could start the preprocessing by simply doing:

```console
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix mereu2020/star_outs/ \
     --readFilesIn mereu2020/mcscrb-seq/SRR9621767_2.fastq.gz mereu2020/mcscrb-seq/SRR9621767_1.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10 \
     --soloCBwhitelist mereu2020/mcscrb-seq/whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

The above scenario is from __Option 1__. If you are doing the experiments by yourself and using __Option 2__, you should do:

```console
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix /path/to/star_outs/ \
     --readFilesIn /path/to/Undetermined_S0_R2_001.fastq.gz /path/to/Undetermined_S0_CB_UMI.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10 \
     --soloCBwhitelist mereu2020/mcscrb-seq/whitelist2.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __SCRB-seq/mcSCRB-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/SCRB-seq.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mix_hg38_mm10/star_index`

>>> Pointing to the directory of the star index. The public data from the above paper was produced using the HCA reference sample, which consists of human PBMCs (60%), and HEK293T (6%), mouse colon (30%), NIH3T3 (3%) and dog MDCK cells (1%). Therefore, we need to use the species mixing reference genome. We also need to add the dog genome, but the dog cells only take 1% of all cells, so I did not bother in this documentation.

`--readFilesCommand zcat`

>>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix mereu2020/star_outs/`

>>> We want to keep everything organised. This directs all output files inside the `mereu2020/star_outs` directory.

`--readFilesIn mereu2020/mcscrb-seq/SRR9621767_2.fastq.gz mereu2020/mcscrb-seq/SRR9621767_1.fastq.gz` or `--readFilesIn /path/to/Undetermined_S0_R2_001.fastq.gz /path/to/Undetermined_S0_CB_UMI.fastq.gz`

>>> If you check the manual, we sould put two files here. The first file is the reads that come from cDNA, and the second file should contain cell barcode and UMI. In __SCRB-seq/mcSCRB-seq__, cDNA reads come from Read 2, and the cell barcode and UMI come from either Read 1 or the `CB_UMI.fastq.gz` file we created as described above. Check [the mcSCRB-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/SCRB-seq.html) if you are not sure.

`--soloType CB_UMI_Simple`

>>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter. This is what we did for __Option 2__ as described above.

`--soloCBstart 1 --soloCBlen 6 --soloUMIstart 7 --soloUMIlen 10` or `--soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10`

>>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist mereu2020/mcscrb-seq/whitelist.txt`

>>> The plain text file containing all possible valid cell barcodes, one per line. The content of this file varies depending on the strategy (__Option 1__ vs __Option 2__) and the primer used during the experiment. We just demonstrated the standard whitelist in the paper.

`--soloStrand Forward`

>>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __SCRB-seq__ and __mcSCRB-seq__, the cDNA reads are from the Read 2 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing the Illumina Read 1 sequence. Therefore, Read 1 comes from the first strand, complementary to the coding strand. Read 2 comes from the coding strand. Therefore, use `Forward` for __SCRB-seq/mcSCRB-seq__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear.

`--outSAMattributes CB UB`

>>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/mereu2020/
├── mcscrb-seq
│   ├── mcSCRB-seq_oligodT.txt
│   ├── SRR9621767_1.fastq.gz
│   ├── SRR9621767_2.fastq.gz
│   ├── whitelist2.txt
│   └── whitelist.txt
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

6 directories, 20 files
```