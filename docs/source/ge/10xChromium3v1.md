# 10x Genomics Single Cell 3' V1

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html) to see how __10x Genomics Single Cell 3' V1__ libraries are generated experimentally. This is a droplet-based method, where cells are captured inside droplets. At the same time, gel beads with barcoded oligo-dT primer containing UMIs are also captured inside the droplet. Reverse transcription happens inside the droplet. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead.

## For Your Own Experiments

The `V1` chemistry is already obsolete, but I'm still providing the preprocessing pipeline for the sake of keeping a record. Although it is highly unlikely that you will do this on your own in future, but just in case, this is the configuration:

| Order | Read             | Cycle | Description                                           |
|-------|------------------|-------|-------------------------------------------------------|
| 1     | Read 1           | >50   | This normally yields `R1_001.fastq.gz`, cDNA reads    |
| 2     | Index 1 (__i7__) | 14    | This normally yields `I1_001.fastq.gz`, Cell barcodes |
| 3     | Index 2 (__i5__) | 8     | This normally yields `I2_001.fastq.gz`, Sample index  |
| 4     | Read 2           | 10    | This normally yields `R2_001.fastq.gz`, UMI           |

Look at the order of the sequencing, as you can see, the first (`R1`), the 2nd (`I1`) and the 4th (`R2`) reads are all important for us. Therefore, you would like to get all of them for each sample based on sample index, that is, the 3rd read (`I2`). You could prepare a `SampleSheet.csv` with the sample index information. Here is an example of `SampleSheet.csv` of a NextSeq run with a sample using standard `i5` indexing primers:

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
75,,,,,,,,,,,
10,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample01,,,,,,,,SI-GA-A1_1,AGGCTGGT,,
Sample01,,,,,,,,SI-GA-A1_2,CACAACTA,,
Sample01,,,,,,,,SI-GA-A1_3,GTTGGTCC,,
Sample01,,,,,,,,SI-GA-A1_4,TTGTAAGA,,
```

You can see each sample actually has four different index sequences. This is because each well from the index plate actually contains four different indices for base balancing. To get the reads you need, you should run `bcl2fastq` in the following way:

```console
bcl2fastq --use-bases-mask=Y75,Y14,I8,Y10 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is `--use-bases-mask=Y75,Y14,I8,Y10`. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y75` at the first position indicates "use the cycle as a real read", so you will get 75-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `Y14` at the second position indicates "use the cycle as a real read", so you will get 14-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read.
3. `I8` at the third position indicates "use the cycle as an index read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read, though it is the 3rd read overall.
4. `Y10` at the fourth position indicates "use the cycle as a real read", so you will get 10-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

Therefore, you will get four fastq file per sample. Using the examples above, these are the files you should get:

```bash
Sample01_S1_I1_001.fastq.gz # 8 bp: sample index
Sample01_S1_R1_001.fastq.gz # 75 bp: cDNA reads
Sample01_S1_R2_001.fastq.gz # 14 bp: cell barcodes
Sample01_S1_R3_001.fastq.gz # 10 bp: UMI
```

We can safely ignore the `I1` files, but the naming here is really different from our normal usage. The `R1` files are good. The `R2` files here actually mean `I1` in our normal usage. The `R3` files here actually mean `R2` in our normal usage. Anyway, __DO NOT get confused__.

To run `starsolo`, we need to get the cell barcodes and the UMI into the same fastq file. This can be simply achieved by stitching `R2` and `R3` together:

```bash
paste <(zcat Sample01_S1_R2_001.fastq.gz) \
      <(zcat Sample01_S1_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > Sample01_S1_CB_UMI.fastq.gz
```

After that, you are ready to go.

## Public Data

The data is from the following paper:

```{eval-rst}
.. note::
  Zheng GXY, Terry JM, Belgrader P, Ryvkin P, Bent ZW, Wilson R, Ziraldo SB, Wheeler TD, McDermott GP, Zhu J, Gregory MT, Shuga J, Montesclaros L, Underwood JG, Masquelier DA, Nishimura SY, Schnall-Levin M, Wyatt PW, Hindson CM, Bharadwaj R, Wong A, Ness KD, Beppu LW, Deeg HJ, McFarland C, Loeb KR, Valente WJ, Ericson NG, Stevens EA, Radich JP, Mikkelsen TS, Hindson BJ, Bielas JH (2017) **Massively parallel digital transcriptional profiling of single cells.** *Nat Commun* 8:14049. https://doi.org/10.1038/ncomms14049
```

where the technology was officially described for the first time. Data can be access from the [10x website](https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5Bfacets%5D%5B0%5D=chemistryVersionAndThroughput&configure%5Bfacets%5D%5B1%5D=pipeline.version&configure%5BhitsPerPage%5D=500). Here, we are using the human HEK293T + mouse NIH3T3 mixture data, which contains roughly 1000 cells:

```console
mkdir -p zheng2017/data
wget -P zheng2017/data -c https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t_3t3/293t_3t3_fastqs.tar
tar xf zheng2017/data/293t_3t3_fastqs.tar -C zheng2017/data/
```
After the extraction, you should see the following files:

```console
scg_prep_test/zheng2017/data/
├── 293t_3t3_fastqs.tar
└── fastqs
    ├── flowcell1
    │   ├── read-I1_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
    │   ├── read-I1_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
    │   ├── read-I1_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
    │   ├── read-I1_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
    │   ├── read-I1_si-CACAACTA_lane-001-chunk-001.fastq.gz
    │   ├── read-I1_si-CACAACTA_lane-002-chunk-000.fastq.gz
    │   ├── read-I1_si-CACAACTA_lane-003-chunk-003.fastq.gz
    │   ├── read-I1_si-CACAACTA_lane-004-chunk-002.fastq.gz
    │   ├── read-I1_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
    │   ├── read-I1_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
    │   ├── read-I1_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
    │   ├── read-I1_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
    │   ├── read-I1_si-TCATCAAG_lane-001-chunk-001.fastq.gz
    │   ├── read-I1_si-TCATCAAG_lane-002-chunk-000.fastq.gz
    │   ├── read-I1_si-TCATCAAG_lane-003-chunk-003.fastq.gz
    │   ├── read-I1_si-TCATCAAG_lane-004-chunk-002.fastq.gz
    │   ├── read-I2_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
    │   ├── read-I2_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
    │   ├── read-I2_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
    │   ├── read-I2_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
    │   ├── read-I2_si-CACAACTA_lane-001-chunk-001.fastq.gz
    │   ├── read-I2_si-CACAACTA_lane-002-chunk-000.fastq.gz
    │   ├── read-I2_si-CACAACTA_lane-003-chunk-003.fastq.gz
    │   ├── read-I2_si-CACAACTA_lane-004-chunk-002.fastq.gz
    │   ├── read-I2_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
    │   ├── read-I2_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
    │   ├── read-I2_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
    │   ├── read-I2_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
    │   ├── read-I2_si-TCATCAAG_lane-001-chunk-001.fastq.gz
    │   ├── read-I2_si-TCATCAAG_lane-002-chunk-000.fastq.gz
    │   ├── read-I2_si-TCATCAAG_lane-003-chunk-003.fastq.gz
    │   ├── read-I2_si-TCATCAAG_lane-004-chunk-002.fastq.gz
    │   ├── read-RA_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
    │   ├── read-RA_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
    │   ├── read-RA_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
    │   ├── read-RA_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
    │   ├── read-RA_si-CACAACTA_lane-001-chunk-001.fastq.gz
    │   ├── read-RA_si-CACAACTA_lane-002-chunk-000.fastq.gz
    │   ├── read-RA_si-CACAACTA_lane-003-chunk-003.fastq.gz
    │   ├── read-RA_si-CACAACTA_lane-004-chunk-002.fastq.gz
    │   ├── read-RA_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
    │   ├── read-RA_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
    │   ├── read-RA_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
    │   ├── read-RA_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
    │   ├── read-RA_si-TCATCAAG_lane-001-chunk-001.fastq.gz
    │   ├── read-RA_si-TCATCAAG_lane-002-chunk-000.fastq.gz
    │   ├── read-RA_si-TCATCAAG_lane-003-chunk-003.fastq.gz
    │   └── read-RA_si-TCATCAAG_lane-004-chunk-002.fastq.gz
    └── flowcell2
        ├── read-I1_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
        ├── read-I1_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
        ├── read-I1_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
        ├── read-I1_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
        ├── read-I1_si-CACAACTA_lane-001-chunk-001.fastq.gz
        ├── read-I1_si-CACAACTA_lane-002-chunk-000.fastq.gz
        ├── read-I1_si-CACAACTA_lane-003-chunk-003.fastq.gz
        ├── read-I1_si-CACAACTA_lane-004-chunk-002.fastq.gz
        ├── read-I1_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
        ├── read-I1_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
        ├── read-I1_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
        ├── read-I1_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
        ├── read-I1_si-TCATCAAG_lane-001-chunk-001.fastq.gz
        ├── read-I1_si-TCATCAAG_lane-002-chunk-000.fastq.gz
        ├── read-I1_si-TCATCAAG_lane-003-chunk-003.fastq.gz
        ├── read-I1_si-TCATCAAG_lane-004-chunk-002.fastq.gz
        ├── read-I2_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
        ├── read-I2_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
        ├── read-I2_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
        ├── read-I2_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
        ├── read-I2_si-CACAACTA_lane-001-chunk-001.fastq.gz
        ├── read-I2_si-CACAACTA_lane-002-chunk-000.fastq.gz
        ├── read-I2_si-CACAACTA_lane-003-chunk-003.fastq.gz
        ├── read-I2_si-CACAACTA_lane-004-chunk-002.fastq.gz
        ├── read-I2_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
        ├── read-I2_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
        ├── read-I2_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
        ├── read-I2_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
        ├── read-I2_si-TCATCAAG_lane-001-chunk-001.fastq.gz
        ├── read-I2_si-TCATCAAG_lane-002-chunk-000.fastq.gz
        ├── read-I2_si-TCATCAAG_lane-003-chunk-003.fastq.gz
        ├── read-I2_si-TCATCAAG_lane-004-chunk-002.fastq.gz
        ├── read-RA_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
        ├── read-RA_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
        ├── read-RA_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
        ├── read-RA_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
        ├── read-RA_si-CACAACTA_lane-001-chunk-001.fastq.gz
        ├── read-RA_si-CACAACTA_lane-002-chunk-000.fastq.gz
        ├── read-RA_si-CACAACTA_lane-003-chunk-003.fastq.gz
        ├── read-RA_si-CACAACTA_lane-004-chunk-002.fastq.gz
        ├── read-RA_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
        ├── read-RA_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
        ├── read-RA_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
        ├── read-RA_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
        ├── read-RA_si-TCATCAAG_lane-001-chunk-001.fastq.gz
        ├── read-RA_si-TCATCAAG_lane-002-chunk-000.fastq.gz
        ├── read-RA_si-TCATCAAG_lane-003-chunk-003.fastq.gz
        └── read-RA_si-TCATCAAG_lane-004-chunk-002.fastq.gz

3 directories, 97 files
```

In reality, it is better to run `bcl2fastq` with the `--create-fastq-for-index-reads` flag without a `SampleSheet.csv`. You should get four fastq files per experiment:

```console
Undetermined_S0_I1_001.fastq.gz    # cell barcodes (14 bp)
Undetermined_S0_I2_001.fastq.gz    # sample index (8 bp)
Undetermined_S0_R1_001.fastq.gz    # cDNA reads (98 bp)
Undetermined_S0_R2_001.fastq.gz    # UMI (10 bp)
```

However, the files from the 10x website are __NOT__ like that because they demultiplexed the sample based on `I2`. They used different sample indices even though there is only one sample. The sample was also split into different flow cells and lanes. That is why there are so many files, but essentially, they are all from the same sample.

We can safely ignore all the `I2` files, and just look at the `I1` (cell barcodes) and `RA` (cDNA + UMI) files. If you look at the content of any `RA` file, you will realise that they are interleaved `fastq` files, containing cDNA and UMI reads next to each other. For example, these are the first 16 lines (4 records) of `flowcell1/read-RA_si-AGGCTGGT_lane-001-chunk-001.fastq.gz`:

```console
@NB500915:156:HYKFKBGXX:1:11101:14387:1086 1:N:0:0
TTCCTGGCCGCCAGAAGATCCACATCTCAAAGAAGTGGGGCTTCACCAAGTTCAATGCTGATGAATTTGAAGACATGGTGGCTGAAAAGCGGCTCATC
+
/AAAAEEE/EEAEEAEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEAEEEEEEEEEEEEAEE/EEEEEEEEAEA/EEEEEEEEA/EA
@NB500915:156:HYKFKBGXX:1:11101:14387:1086 4:N:0:0
GCACGNGNTN
+
A//AA#A#E#
@NB500915:156:HYKFKBGXX:1:11101:25884:1109 1:N:0:0
GACCTTTTGGCATGGCCCAGACTGGGGTGCCCTTTGGGGAAGTAAGCATGGTCCGGGACTGGTTGGGCATTGTGGGGCGTGTGCTGACCCATACCCAA
+
AAA/AEEEEEEAEEEAE/EE/EEEEEEEEEAEEEEEE/EEEEEEEEEEEEA</AEE</AEE<A<EEEEEEEAEAA/E/////AE/E/E6E/E/<////
@NB500915:156:HYKFKBGXX:1:11101:25884:1109 4:N:0:0
GTAGTTTTGG
+
A////AEEEE
```

## Reformat FastQ Files

To use `starsolo`, we need to prepare `fastq` files into a file containing cDNA reads and a file with cell barcode + UMI. To get the cDNA reads, we need every other read from the `RA` file:

```bash
mkdir -p zheng2017/data/combined_fastqs

# for cDNA reads
# get the lines whose line number (NR) mod 8 is between 1 and 4
zcat zheng2017/data/fastqs/flowcell1/read-RA_si-*.gz \
     zheng2017/data/fastqs/flowcell2/read-RA_si-*.gz | \
     awk 'NR%8>=1&&NR%8<=4' | \
     gzip > zheng2017/data/combined_fastqs/cDNA_reads.fastq.gz
```

Then, we should append the UMI from the `RA` file to the cell barcode `I1`. This can be achieved using a one liner, but if you are not comfortable, you can split it into different steps for readability.

```bash
# for UMI reads
# get the lines whose line number (NR) mod 8 is 5, 6, 7 or 0
paste <(zcat zheng2017/data/fastqs/flowcell1/read-I1_si-*.gz \
             zheng2017/data/fastqs/flowcell2/read-I1_si-*.gz) \
      <(zcat zheng2017/data/fastqs/flowcell1/read-RA_si-*.gz \
             zheng2017/data/fastqs/flowcell2/read-RA_si-*.gz | \
             awk 'NR%8==5||NR%8==6||NR%8==7||NR%8==0') | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > zheng2017/data/combined_fastqs/CB_UMI_reads.fastq.gz
```

The files `cDNA_reads.fastq.gz` and `CB_UMI_reads.fastq.gz` are just what we need.

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the `V1` chemistry. If you have `cellranger` in your computer, you will find a file called `737K-april-2014_rc.txt` in the `lib/python/cellranger/barcodes/` directory. If you don't have `cellranger`, I have prepared the file for you:

```console
# download the whitelist 
wget -P zheng2017/data/ https://teichlab.github.io/scg_lib_structs/data/737K-april-2014_rc.txt.gz
gunzip zheng2017/data/737K-april-2014_rc.txt.gz
```

## From FastQ To Count Matrix

Now we could start the preprocessing by simply doing:

```console
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix zheng2017/star_outs/ \
     --readFilesIn zheng2017/data/combined_fastqs/cDNA_reads.fastq.gz zheng2017/data/combined_fastqs/CB_UMI_reads.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10 \
     --soloCBwhitelist zheng2017/data/737K-april-2014_rc.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __10x Genomics Single Cell 3' V1__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mix_hg38_mm10/star_index`

>> Pointing to the directory of the star index. The public data from the 10x website is human HEK293T + mouse NIH3T3 cell mixtures. Therefore, we need to use the species mixing reference genome.

`--readFilesCommand zcat`

>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix zheng2017/star_outs/`

>> We want to keep everything organised. This directs all output files inside the `zheng2017/star_outs` directory.

`--readFilesIn zheng2017/data/combined_fastqs/cDNA_reads.fastq.gz zheng2017/data/combined_fastqs/CB_UMI_reads.fastq.gz`

>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. We have gone through all the trouble to generate those files using the procedures described above.

`--soloType CB_UMI_Simple`

>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter. That is why went through those procedures to reformat the `fastq` files.

`--soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10`

>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist zheng2017/data/737K-april-2014_rc.txt`

>> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell 3' V1__ is a commercial platform. The whitelist is taken from their commercial software `cellranger`.

`--soloCellFilter EmptyDrops_CR`

>> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter. 

`--soloStrand Forward`

>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __10x Genomics Single Cell 3' V1__, the cDNA reads are from the Read 1 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing the Illumina Read 2 sequence. Therefore, Read 2 comes from the first strand, complementary to the coding strand. Read 1 comes from the coding strand. Therefore, use `Forward` for __10x Genomics Single Cell 3' V1__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear.

`--outSAMattributes CB UB`

>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/zheng2017
├── data
│   ├── 293t_3t3_fastqs.tar
│   ├── 737K-april-2014_rc.txt
│   ├── combined_fastqs
│   │   ├── CB_UMI_reads.fastq.gz
│   │   └── cDNA_reads.fastq.gz
│   └── fastqs
│       ├── flowcell1
│       │   ├── read-I1_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
│       │   ├── read-I1_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
│       │   ├── read-I1_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
│       │   ├── read-I1_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
│       │   ├── read-I1_si-CACAACTA_lane-001-chunk-001.fastq.gz
│       │   ├── read-I1_si-CACAACTA_lane-002-chunk-000.fastq.gz
│       │   ├── read-I1_si-CACAACTA_lane-003-chunk-003.fastq.gz
│       │   ├── read-I1_si-CACAACTA_lane-004-chunk-002.fastq.gz
│       │   ├── read-I1_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
│       │   ├── read-I1_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
│       │   ├── read-I1_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
│       │   ├── read-I1_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
│       │   ├── read-I1_si-TCATCAAG_lane-001-chunk-001.fastq.gz
│       │   ├── read-I1_si-TCATCAAG_lane-002-chunk-000.fastq.gz
│       │   ├── read-I1_si-TCATCAAG_lane-003-chunk-003.fastq.gz
│       │   ├── read-I1_si-TCATCAAG_lane-004-chunk-002.fastq.gz
│       │   ├── read-I2_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
│       │   ├── read-I2_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
│       │   ├── read-I2_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
│       │   ├── read-I2_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
│       │   ├── read-I2_si-CACAACTA_lane-001-chunk-001.fastq.gz
│       │   ├── read-I2_si-CACAACTA_lane-002-chunk-000.fastq.gz
│       │   ├── read-I2_si-CACAACTA_lane-003-chunk-003.fastq.gz
│       │   ├── read-I2_si-CACAACTA_lane-004-chunk-002.fastq.gz
│       │   ├── read-I2_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
│       │   ├── read-I2_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
│       │   ├── read-I2_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
│       │   ├── read-I2_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
│       │   ├── read-I2_si-TCATCAAG_lane-001-chunk-001.fastq.gz
│       │   ├── read-I2_si-TCATCAAG_lane-002-chunk-000.fastq.gz
│       │   ├── read-I2_si-TCATCAAG_lane-003-chunk-003.fastq.gz
│       │   ├── read-I2_si-TCATCAAG_lane-004-chunk-002.fastq.gz
│       │   ├── read-RA_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
│       │   ├── read-RA_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
│       │   ├── read-RA_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
│       │   ├── read-RA_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
│       │   ├── read-RA_si-CACAACTA_lane-001-chunk-001.fastq.gz
│       │   ├── read-RA_si-CACAACTA_lane-002-chunk-000.fastq.gz
│       │   ├── read-RA_si-CACAACTA_lane-003-chunk-003.fastq.gz
│       │   ├── read-RA_si-CACAACTA_lane-004-chunk-002.fastq.gz
│       │   ├── read-RA_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
│       │   ├── read-RA_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
│       │   ├── read-RA_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
│       │   ├── read-RA_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
│       │   ├── read-RA_si-TCATCAAG_lane-001-chunk-001.fastq.gz
│       │   ├── read-RA_si-TCATCAAG_lane-002-chunk-000.fastq.gz
│       │   ├── read-RA_si-TCATCAAG_lane-003-chunk-003.fastq.gz
│       │   └── read-RA_si-TCATCAAG_lane-004-chunk-002.fastq.gz
│       └── flowcell2
│           ├── read-I1_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
│           ├── read-I1_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
│           ├── read-I1_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
│           ├── read-I1_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
│           ├── read-I1_si-CACAACTA_lane-001-chunk-001.fastq.gz
│           ├── read-I1_si-CACAACTA_lane-002-chunk-000.fastq.gz
│           ├── read-I1_si-CACAACTA_lane-003-chunk-003.fastq.gz
│           ├── read-I1_si-CACAACTA_lane-004-chunk-002.fastq.gz
│           ├── read-I1_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
│           ├── read-I1_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
│           ├── read-I1_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
│           ├── read-I1_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
│           ├── read-I1_si-TCATCAAG_lane-001-chunk-001.fastq.gz
│           ├── read-I1_si-TCATCAAG_lane-002-chunk-000.fastq.gz
│           ├── read-I1_si-TCATCAAG_lane-003-chunk-003.fastq.gz
│           ├── read-I1_si-TCATCAAG_lane-004-chunk-002.fastq.gz
│           ├── read-I2_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
│           ├── read-I2_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
│           ├── read-I2_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
│           ├── read-I2_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
│           ├── read-I2_si-CACAACTA_lane-001-chunk-001.fastq.gz
│           ├── read-I2_si-CACAACTA_lane-002-chunk-000.fastq.gz
│           ├── read-I2_si-CACAACTA_lane-003-chunk-003.fastq.gz
│           ├── read-I2_si-CACAACTA_lane-004-chunk-002.fastq.gz
│           ├── read-I2_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
│           ├── read-I2_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
│           ├── read-I2_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
│           ├── read-I2_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
│           ├── read-I2_si-TCATCAAG_lane-001-chunk-001.fastq.gz
│           ├── read-I2_si-TCATCAAG_lane-002-chunk-000.fastq.gz
│           ├── read-I2_si-TCATCAAG_lane-003-chunk-003.fastq.gz
│           ├── read-I2_si-TCATCAAG_lane-004-chunk-002.fastq.gz
│           ├── read-RA_si-AGGCTGGT_lane-001-chunk-001.fastq.gz
│           ├── read-RA_si-AGGCTGGT_lane-002-chunk-000.fastq.gz
│           ├── read-RA_si-AGGCTGGT_lane-003-chunk-003.fastq.gz
│           ├── read-RA_si-AGGCTGGT_lane-004-chunk-002.fastq.gz
│           ├── read-RA_si-CACAACTA_lane-001-chunk-001.fastq.gz
│           ├── read-RA_si-CACAACTA_lane-002-chunk-000.fastq.gz
│           ├── read-RA_si-CACAACTA_lane-003-chunk-003.fastq.gz
│           ├── read-RA_si-CACAACTA_lane-004-chunk-002.fastq.gz
│           ├── read-RA_si-GTTGGTCC_lane-001-chunk-001.fastq.gz
│           ├── read-RA_si-GTTGGTCC_lane-002-chunk-000.fastq.gz
│           ├── read-RA_si-GTTGGTCC_lane-003-chunk-003.fastq.gz
│           ├── read-RA_si-GTTGGTCC_lane-004-chunk-002.fastq.gz
│           ├── read-RA_si-TCATCAAG_lane-001-chunk-001.fastq.gz
│           ├── read-RA_si-TCATCAAG_lane-002-chunk-000.fastq.gz
│           ├── read-RA_si-TCATCAAG_lane-003-chunk-003.fastq.gz
│           └── read-RA_si-TCATCAAG_lane-004-chunk-002.fastq.gz
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

10 directories, 115 files
```