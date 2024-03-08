# inDrop

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html) to see how __inDrop__ libraries are generated experimentally. This is a droplet-based method, where cells are captured inside droplets. At the same time, gel beads with barcoded oligo-dT primer containing UMIs are also captured inside the droplet. Reverse transcription happens inside the droplet. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead. The barcodes on the gel beads are generated using a split-pool based combinatorial indexing strategy.

```{eval-rst}
.. important::
  
  1. Be aware that there are different versions of **inDrop**. There is also a commercial version from `1CellBio <https://1cell-bio.com>`_. They have different adaptor sequences and final library structures. Therefore, they require different parameters for the preprocessing steps. Check the `inDrop GitHub Page <https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html>`_ and `the inDrop GitHub repository <https://github.com/indrops/indrops>`_ for some details.
  2. The real situation is more complicated. If you do the experiments by yourself, you probably know the sequences. If you use a commercial inDrop platform or analyse published data, make sure you check the kit manual or the authors of the publication about the details of the sequences used in the protocol.
  3. **Some helpful threads:** `STARsolo issue #605 <https://github.com/alexdobin/STAR/issues/605>`_, `STARsolo issue #785 <https://github.com/alexdobin/STAR/issues/785>`_ and `indrops issue #32 <https://github.com/indrops/indrops/issues/32>`_
  4. If you find errors, which is highly likely, please do let me know via email ``chenx9@sustech.edu.cn`` or `raised an issue <https://github.com/Teichlab/scg_lib_structs/issues>`_. Thank you in advance.

```

## For Your Own Experiments

The read configuration is a bit more complicated and depends on the version.

### The V1 Configuration

| Order | Read             | Cycle    | Description                                                           |
|-------|------------------|----------|-----------------------------------------------------------------------|
| 1     | Read 1           | >=47     | This yields `R1_001.fastq.gz`, Cell barcodes, linker sequence and UMI |
| 2     | Index 1 (__i7__) | 6 or 8   | This yields `I1_001.fastq.gz`, sample index                           |
| 3     | Index 2 (__i5__) | Optional | This yields `I2_001.fastq.gz`, not really used but can be present     |
| 4     | Read 2           | >=35     | This yields `R2_001.fastq.gz`, cDNA reads                             |

The content of __Read 1__ is like this:

| Length | Sequence (5' -> 3')                                                                         |
|--------|---------------------------------------------------------------------------------------------|
| >46    | 8 - 11 bp + __Barcode1__ GAGTGATTGCTTGTGACGCCTT + 8 bp __Barcode2__ + 6 bp __UMI__ + poly-T |

### The V2 Configuration

| Order | Read             | Cycle    | Description                                                           |
|-------|------------------|----------|-----------------------------------------------------------------------|
| 1     | Read 1           | >=35     | This yields `R1_001.fastq.gz`, cDNA reads                             |
| 2     | Index 1 (__i7__) | 6 or 8   | This yields `I1_001.fastq.gz`, sample index                           |
| 3     | Index 2 (__i5__) | Optional | This yields `I2_001.fastq.gz`, not really used but can be present     |
| 4     | Read 2           | >=47     | This yields `R2_001.fastq.gz`, Cell barcodes, linker sequence and UMI |

The content of __Read 2__ is like this:

| Length | Sequence (5' -> 3')                                                                         |
|--------|---------------------------------------------------------------------------------------------|
| >46    | 8 - 11 bp + __Barcode1__ GAGTGATTGCTTGTGACGCCTT + 8 bp __Barcode2__ + 6 bp __UMI__ + poly-T |

### The V3 Configuration

This configuration is more complicated and the naming of the output files does not really follow our normal convention. __DO NOT__ get confused.

| Order | Read             | Cycle  | Description                                                          |
|-------|------------------|--------|----------------------------------------------------------------------|
| 1     | Read 1           | >50    | This normally yields `R1_001.fastq.gz`, cDNA reads                   |
| 2     | Index 1          | 8      | This normally yields `I1_001.fastq.gz`, __Barcode1__ of the gel bead |
| 3     | Index 2 (__i7__) | 6 or 8 | This normally yields `I2_001.fastq.gz`, sample index                 |
| 4     | Read 2           | 14     | This normally yields `R2_001.fastq.gz`, __Barcode2 + UMI__           |

The content of __Read 2__ is like this:

| Length | Sequence (5' -> 3')              |
|--------|----------------------------------|
| 14     | 8 bp __Barcode2__ + 6 bp __UMI__ |

You can think of the combination of 8 - 11 bp __Barcode1__ and 8 bp __Barcode2__ as the cell barcodes. If you use this method, you have to sequence the library on your own, you need to get the `fastq` files by running `bcl2fastq` by yourself. In this case it is better to write a `SampleSheet.csv` with `i7` indices for each sample. Note that the `i7` index is the 2nd read (Index 1) in __V1__ and __V2__, but the 3rd read (Index 2) in __V3__, so you need to adjust that accordingly. Here are examples of `SampleSheet.csv` of NextSeq runs with two samples using some standard index:

### The V1 & V2 SampleSheet

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
50,,,,,,,,,,,
50,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample01,,,,,,BC1,CAGATC,,,,
Sample02,,,,,,BC2,ACTTGA,,,,
```

### The V3 SampleSheet

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
50,,,,,,,,,,,
14,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample01,,,,,,BC1,CAGATC,,,,
Sample02,,,,,,BC2,ACTTGA,,,,
```

You need to run `bcl2fastq` differently based on the configuration like this:

```bash
# for V1 and V2 configuration

bcl2fastq --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4

# for V3 configuration

bcl2fastq --use-bases-mask=Y50,Y8,I6,Y14 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is `--use-bases-mask=Y50,Y8,I6,Y14` in the __V3__ configuration. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y50` at the first position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `Y8` at the second position indicates "use the cycle as a real read", so you will get 8-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read.
3. `I6` at the third position indicates "use the cycle as an index read", so you will get 6-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read, though it is the 3rd read overall.
4. `Y14` at the fourth position indicates "use the cycle as a real read", so you will get 14-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

After that, you will get two files per sample in the __V1__ and __V2__ configurations, and four files per sample in the __V3__ configuration:

```bash
# V1 configuration
Sample01_S1_R1_001.fastq.gz # 50 bp: barcodes, UMI and linkers, poly-T
Sample01_S1_R2_001.fastq.gz # 50 bp: cDNA reads
Sample02_S2_R1_001.fastq.gz # 50 bp: barcodes, UMI and linkers, poly-T
Sample02_S2_R2_001.fastq.gz # 50 bp: cDNA reads

# V2 configuration
Sample01_S1_R1_001.fastq.gz # 50 bp: cDNA reads
Sample01_S1_R2_001.fastq.gz # 50 bp: barcodes, UMI and linkers, poly-T
Sample02_S2_R1_001.fastq.gz # 50 bp: cDNA reads
Sample02_S2_R2_001.fastq.gz # 50 bp: barcodes, UMI and linkers, poly-T

# V3 configuration
Sample01_S1_I1_001.fastq.gz # 6 bp: sample index, can be ignored
Sample01_S1_R1_001.fastq.gz # 50 bp: cDNA reads
Sample01_S1_R2_001.fastq.gz # 8 bp: Barcode1
Sample01_S1_R3_001.fastq.gz # 14 bp: Barcode2 + UMI
Sample02_S2_I1_001.fastq.gz # 6 bp: sample index, can be ignored
Sample02_S2_R1_001.fastq.gz # 50 bp: cDNA reads
Sample02_S2_R2_001.fastq.gz # 8 bp: Barcode1
Sample01_S2_R3_001.fastq.gz # 14 bp: Barcode2 + UMI
```

For the __V1__ and __V2__ configurations, you are ready to go. However, since the cell barcodes and UMI are distributed in different reads, we need to collect them into one `fastq` file in order to use `starsolo`. This can be done by simple stitching the reads:

```bash
# Sample01
paste <(zcat Sample01_S1_R2_001.fastq.gz) \
      <(zcat Sample01_S1_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > Sample01_S1_CB_UMI.fastq.gz

# Sample02
paste <(zcat Sample02_S2_R2_001.fastq.gz) \
      <(zcat Sample02_S2_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > Sample02_S2_CB_UMI.fastq.gz
```

After that, you are ready to go.

## Public Data

For the purpose of demonstration, we will use the __inDrop__ data from the following paper:

```{eval-rst}
.. note::

  Mereu E, Lafzi A, Moutinho C, Ziegenhain C, McCarthy DJ, Álvarez-Varela A, Batlle E, Sagar, Grün D, Lau JK, Boutet SC, Sanada C, Ooi A, Jones RC, Kaihara K, Brampton C, Talaga Y, Sasagawa Y, Tanaka K, Hayashi T, Braeuning C, Fischer C, Sauer S, Trefzer T, Conrad C, Adiconis X, Nguyen LT, Regev A, Levin JZ, Parekh S, Janjic A, Wange LE, Bagnoli JW, Enard W, Gut M, Sandberg R, Nikaido I, Gut I, Stegle O, Heyn H (2020) **Benchmarking single-cell RNA-sequencing protocols for cell atlas projects.** *Nat Biotechnol* 38:747–755. https://doi.org/10.1038/s41587-020-0469-4

```

where the authors benchmarked quite a few different scRNA-seq methods using a standardised sample: a mixture of different human, mouse and dog cells. We are going to use the data from the __inDrop__ method, which uses the commercial __1CellBio__ platform. You can download the `fastq` file from [this ENA page](https://www.ebi.ac.uk/ena/browser/view/PRJNA551759?show=reads):

```console
# get fastq files
mkdir -p mereu2020/indrop
wget -P mereu2020/indrop -c \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR962/004/SRR9621794/SRR9621794_1.fastq.gz \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR962/004/SRR9621794/SRR9621794_2.fastq.gz
```

## Prepare Whitelist

The full oligo sequences can be found in the [Supplementary Table S2](https://teichlab.github.io/scg_lib_structs/data/inDrop/41596_2017_BFnprot2016154_MOESM456_ESM.xlsx) and [Supplementary Table S3](https://teichlab.github.io/scg_lib_structs/data/inDrop/41596_2017_BFnprot2016154_MOESM457_ESM.xlsx) from the [__inDrop Nature Protocols paper__](https://www.nature.com/articles/nprot.2016.154). As you can see, there are a total of 384 different __Barcode1__ with 8 - 11 bp in length, and 384 different __Barcode2__. The oligos are added to the gel beads by primer extension during the split-pool procedures. The cell barcodes are basically the combination of __Barcode1__ and __Barcode2__. There will be a total of __384 * 384 = 147456__ possible cell barcodes. I have organised the oligo information into two tables here (only showing 5 records):

__Barcode1__

| Name         | Sequence    | Reverse complement |
|--------------|-------------|:------------------:|
| W1-bc1.1-PE1 | AAACAAAC    |      GTTTGTTT      |
| W1-bc1.2-PE1 | AAACACGGT   |      ACCGTGTTT     |
| W1-bc1.3-PE1 | AAACACTATC  |      GATAGTGTTT    |
| W1-bc1.4-PE1 | AAACCGCCTCA |      TGAGGCGGTTT   |

__Barcode2__

| Name             | Sequence | Reverse complement |
|------------------|----------|:------------------:|
| BA19-N6-bc2.1-W1 | AAACAAAC |      GTTTGTTT      |
| BA19-N6-bc2.2-W1 | AAACACGG |      CCGTGTTT      |
| BA19-N6-bc2.3-W1 | AAACACTA |      TAGTGTTT      |
| BA19-N6-bc2.4-W1 | AAACCGCC |      GGCGGTTT      |

__You should notice that `Barcode1` has variable lengths, but the first 8 bp are exactly the same as `Barcode2`__. I have prepared the full tables in `csv` format for you to download:

[inDrop_Barcode1.csv](https://teichlab.github.io/scg_lib_structs/data/inDrop/inDrop_Barcode1.csv)  
[inDrop_Barcode2.csv](https://teichlab.github.io/scg_lib_structs/data/inDrop/inDrop_Barcode2.csv)

Let's download them to generate the whitelist:

```console
wget -P mereu2020/indrop \
    https://teichlab.github.io/scg_lib_structs/data/inDrop/inDrop_Barcode1.csv \
    https://teichlab.github.io/scg_lib_structs/data/inDrop/inDrop_Barcode2.csv
```

Now we need to generate the whitelist of those two sets of barcodes. Read very carefully of the [__inDrop GitHub page__](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html). Pay attention to the oligo orientation. The barcode sequences that we get from the [__inDrop Nature Protocols paper__](https://www.nature.com/articles/nprot.2016.154) are the sequences in the adaptors, which are used to generate the bead oligos. Therefore, the sequences on the bead oligos are reverse complement to the actual barcodes. Now, you can see that in the __V1__ and __V2__ configuration, __Barcode1__ and __Barcode2__ are in the same read and in the same direction of the bead oligo. Therefore, we should use the reverse complement of the barcode sequences for the whitelists. In the __V3__ configuration, __Barcode1__ is sequenced in the opposite direction of the bead oligo with only 8 cycles, so we need to use the first 8 bp of __Barcode1__ as they are. __Barcode2__ is sequenced in the same direction of the bead oligo, so we should take the reverse complement of the barcode sequence. In addition, since we stitch __Barcode1__, __Barcode2__ and __UMI__ together into the `CB_UMI.fastq.gz`, we should generate all possible combinations of the __Barcode1_8bp + Barcode2 rc__ as the whitelist. Here is how you could do this:

```bash
# for V1 and V2, prepare two plain lists from
# the reverse complement of Barcode1 and Barcode2
tail -n +2 mereu2020/indrop/inDrop_Barcode1.csv | 
    cut -f 3 -d, > mereu2020/indrop/V1and2_BC1.txt

tail -n +2 mereu2020/indrop/inDrop_Barcode2.csv | 
    cut -f 3 -d, > mereu2020/indrop/V1and2_BC2.txt

# for V3, we need to get all possible combination the first 8bp of the Barcode1
# to the reverse complement of Barcode2 sequence
for x in $(tail -n +2 mereu2020/indrop/inDrop_Barcode1.csv | cut -f 2 -d, | cut -c 1-8); do
    for y in $(tail -n +2 mereu2020/indrop/inDrop_Barcode2.csv | cut -f 3 -d,); do
        echo "${x}${y}"
        done
    done > mereu2020/indrop/V3_whitelist.txt
```

## From FastQ To Count Matrix

The public data was produced using the 1CellBio inDrop platform. Only have two fastq files are there. If you `grep` the sequence GAGTGATTGCTTGTGACGCCTT, it appears in __Read 2__ more often, so I assume it was based on __V2__. However, the perfect match only appears in ~24% of __Read 2__. I'm not entirely sure what's going on, maybe the commercial __1CellBio__ platform has different mechanism. Anyway, let's move on as if it is okay for the sake of demonstration.

```console
# map and generate the count matrix

STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix mereu2020/star_outs/ \
     --readFilesIn mereu2020/indrop/SRR9621794_1.fastq.gz mereu2020/indrop/SRR9621794_2.fastq.gz \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
     --soloAdapterMismatchesNmax 4 \
     --soloCBposition 0_0_2_-1 3_1_3_8 \
     --soloUMIposition 3_9_3_14 \
     --soloCBwhitelist mereu2020/indrop/V1and2_BC1.txt mereu2020/indrop/V1and2_BC2.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __inDrop__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mix_hg38_mm10/star_index`

>> Pointing to the directory of the star index. The public data from the above paper was produced using the HCA reference sample, which consists of human PBMCs (60%), and HEK293T (6%), mouse colon (30%), NIH3T3 (3%) and dog MDCK cells (1%). Therefore, we need to use the species mixing reference genome. We also need to add the dog genome, but the dog cells only take 1% of all cells, so I did not bother in this documentation.

`--readFilesCommand zcat`

>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix mereu2020/star_outs/`

>> We want to keep everything organised. This parameter directs all output files into the `mereu2020/star_outs/` directory.

`--readFilesIn`

>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second file should contain cell barcode and UMI. In __inDrop V2__, cDNA reads come from Read 1, and the cell barcode and UMI come from Read 2. Check [the inDrop GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html) if you are not sure.

`--soloType CB_UMI_Complex`

>> Since Read 2 not only has cell barcodes and UMI, the common linker sequences are also there. The cell barcodes are non-consecutive, separated by the linker sequences. In this case, we have to use the `CB_UMI_Complex` option. Of course, we could also extract them upfront into a new `fastq` file, but that's slow. It is better to use this option.

`--soloAdapterSequence GAGTGATTGCTTGTGACGCCTT`

>> The 8 - 11 bp variable length of __Barcode1__ at the beginning of __Read 2__ makes the situation complicated, because the absolute positions of __Barcode2__ and __UMI__ in each read will vary. However, by specifying an adaptor sequence, we could use this sequence as an anchor, and tell the program where cell barcodes and UMI are located relatively to the anchor. See below.

`--soloAdapterMismatchesNmax 3`

>> The number of mismatches are tolerated during the adapter finding. The adapter here is a bit long, so I want a bit relaxed matching, but you may want to try a few different options, like 1 (the default) or 2.

`--soloCBposition` and `--soloUMIposition`

>> These options specify the locations of cell barcode and UMI in the 2nd fastq files we passed to `--readFilesIn`. In this case, it is __Read 2__. Read the [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for more details. I have drawn a picture to help myself decide the exact parameters. There are some freedom here depending on what you are using as anchors. in __inDrop V1 & 2__, the __Barcode1__ has variable lengths, the absolute positions of __Barcode2__ and __UMI__ are variable. Therefore, using Read start as anchor will not work for them. We need to use the adaptor as the anchor, and specify the positions relative to the anchor. See the image:

![](https://teichlab.github.io/scg_lib_structs/data/inDrop/Star_CB_UMI_Complex_inDrop.jpg)

`--soloCBwhitelist`

>> Since the real cell barcodes consists of two non-consecutive parts: two sets of barcodes. The whitelist here is the combination of those two lists. We should provide them separately in the specified order and `star` will take care of the combinations.

`--soloCBmatchWLtype 1MM`

>> How stringent we want the cell barcode reads to match the whitelist. The default option (`1MM_Multi`) does not work here. We choose this one here for simplicity, but you might want to experimenting different parameters to see what the difference is.

`--soloCellFilter EmptyDrops_CR`

>> Experiments are never perfect. Even for barcodes that do not capture the molecules inside the cells, you may still get some reads due to various reasons, such as ambient RNA or DNA and leakage. In general, the number of reads from those cell barcodes should be much smaller, often orders of magnitude smaller, than those barcodes that come from real cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Forward`

>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __inDrop V2__, the cDNA reads are from the Read 1 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing UMI and Read 2 sequencing primer. Therefore, Read 2 consists of cell barcodes and UMI. They come from the first strand, complementary to the coding strand. Read 1 comes from the coding strand. Therefore, use `Forward` for __inDrop V2__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the inDrop GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html) if you are not sure.

`--outSAMattributes CB UB`

>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/mereu2020/
├── indrop
│   ├── inDrop_Barcode1.csv
│   ├── inDrop_Barcode2.csv
│   ├── SRR9621794_1.fastq.gz
│   ├── SRR9621794_2.fastq.gz
│   ├── V1and2_BC1.txt
│   ├── V1and2_BC2.txt
│   └── V3_whitelist.txt
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

6 directories, 22 files
```

## Some Extra Work

If you check the content of the `mereu2020/star_outs/Solo.out/Gene/Summary.csv` file, you will see that only ~24% of reads have valid barcodes. The commercial platform is difficult to crack. Let's use the sample data provided by the [indrops GitHub repository](https://github.com/indrops/indrops) to see if the above procedures we build are good or not.

```bash
# first, clone the indrops github
git clone https://github.com/indrops/indrops.git

# process V1 data, remember Read 2 is the cDNA in V1
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix indrops/v1_outs/ \
     --readFilesIn indrops/test/seq_runs/run_v1_single_file/R2.fastq.gz indrops/test/seq_runs/run_v1_single_file/R1.fastq.gz \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
     --soloAdapterMismatchesNmax 3 \
     --soloCBposition 0_0_2_-1 3_1_3_8 \
     --soloUMIposition 3_9_3_14 \
     --soloCBwhitelist mereu2020/indrop/V1and2_BC1.txt mereu2020/indrop/V1and2_BC2.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate

# process V2 data, remember Read 1 is the cDNA in V2
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix indrops/v2_outs/ \
     --readFilesIn indrops/test/seq_runs/run_v2_single_file/R1.fastq.gz indrops/test/seq_runs/run_v2_single_file/R2.fastq.gz \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GAGTGATTGCTTGTGACGCCTT \
     --soloAdapterMismatchesNmax 3 \
     --soloCBposition 0_0_2_-1 3_1_3_8 \
     --soloUMIposition 3_9_3_14 \
     --soloCBwhitelist mereu2020/indrop/V1and2_BC1.txt mereu2020/indrop/V1and2_BC2.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

If you check the `Summary.csv` files, more than 80% of reads contain valid barcodes. We know we are good.

In terms of __V3__, the file names are different from what we described at the beginning of the documentation. There are four types of fastq files in the repository directory `indrops/test/seq_runs/run_v3`. To put them into context in this documentation:

| Files in the indrops repo       | Files in this documentation | Description            |
|---------------------------------|-----------------------------|------------------------|
| Undetermined_S0_R1_001.fastq.gz | R1_001.fastq.gz             | cDNA reads             |
| Undetermined_S0_R2_001.fastq.gz | R2_001.fastq.gz             | __Barcode1__           |
| Undetermined_S0_R3_001.fastq.gz | I1_001.fastq.gz             | Sample index, ignored! |
| Undetermined_S0_R4_001.fastq.gz | R3_001.fastq.gz             | __Barcode2__ + __UMI__ |

```bash
# the files are split by lanes, so let's merge them
cat indrops/test/seq_runs/run_v3/*_R1_001.fastq.gz > indrops/test/seq_runs/run_v3/cDNA.fastq.gz

# stitch Barcode1 and Barcode2+UMI file
paste <(zcat indrops/test/seq_runs/run_v3/*_R2_001.fastq.gz) \
      <(zcat indrops/test/seq_runs/run_v3/*_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > indrops/test/seq_runs/run_v3/CB_UMI.fastq.gz

# do the mapping
# Note R4 contains some extra T after Barcode2+UMI
# so we need --soloBarcodeReadLength 0 to turn off length check
STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix indrops/v3_outs/ \
     --readFilesIn indrops/test/seq_runs/run_v3/cDNA.fastq.gz indrops/test/seq_runs/run_v3/CB_UMI.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 6 \
     --soloBarcodeReadLength 0 \
     --soloCBwhitelist mereu2020/indrop/V3_whitelist.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

If you check the `Summary.csv`, it seems >60% reads have valid barcodes. This might be okay, but I expect to have higher fraction.