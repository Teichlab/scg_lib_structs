# sci-RNA-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html) to see how __sci-RNA-seq__ libraries are generated experimentally. This is a split-pool based combinatorial indexing strategy, where fixed cells are used as the reaction chamber. mRNA molecules are marked by oligo-dT primer with distinct barcodes in 96 or 384 minibulk reactions in the plate format (the first plate). Then all cells are pooled and randomly distributed into a new 96- or 384-well plate (the second plate). Library preparation is performed using the Tn5-based Illumina Nextera strategy to add __i5__ and __i7__ indices. Single cells can be identified by the combination of the RT barcode and __i5 + i7__. In addition, another level of barcode can be added during the tagmentation by barcoded Tn5, but this documentation will just focus on two-level barcodes, without the Tn5 index.

## For Your Own Experiments

The read configuration is the same as a standard library:

| Order | Read             | Cycle   | Description                                                   |
|-------|------------------|---------|---------------------------------------------------------------|
| 1     | Read 1           | 18      | This yields `R1_001.fastq.gz`, 8 bp UMI + 10 bp RT barcode    |
| 2     | Index 1 (__i7__) | 8 or 10 | This yields `I1_001.fastq.gz`, well barcode for the 2nd plate |
| 3     | Index 2 (__i5__) | 8 or 10 | This yields `I2_001.fastq.gz`, well barcode for the 2nd plate |
| 4     | Read 2           | >50     | This yields `R2_001.fastq.gz`, cDNA reads                     |

You can think of the 10 bp __RT barcode__ as the well barcode for the 1st plate, and __i7 + i5__ are the well barcode for the 2nd plate. For a cell, it can go into a well in the 1st plate, then another well in the 2nd plate. Different cells have very low chance of going through the same combination of wells in the two plates. Therefore, if reads have the same combination of well barcodes (__RT barcode + i7 + i5__), we can safely think they are from the same cell.

If you sequence the library via your core facility or a company, you need to provide the `i5` and `i7` index sequence you used during the library PCR. Like mentioned previously, they are basically the well barcode for the 2nd plate. Then you will get two `fastq` files (`R1` and `R2`) per well. The total file number will depend on how many wells in the 2nd plate you are processing.

If you sequence the library on your own, you need to get the `fastq` files by running `bcl2fastq` by yourself. In this case you have two choices. One choice is to provide a `SampleSheet.csv` with `i7` and `i5` indices for each well in the 2nd plate. This will yield the `fastq` files similar to those from your core facility or the company. The other choice is to run `bcl2fastq` without a `SampleSheet.csv` but with the `--create-fastq-for-index-reads` flag, where all the reads from a run will be simply dumped into four simple `fastq` files. There are pro and cons for both choices. Let's look into more details.

### Run bcl2fastq With Strategy 1

I actually prefer this way. Here is an example of the `SampleSheet.csv` from a NextSeq run with a full 96-well plate (2nd plate) using some standard indices:

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
18,,,,,,,,,,,
52,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
A1,,,,,,N701,TAAGGCGA,S502,ATAGAGAG,,
A2,,,,,,N702,CGTACTAG,S502,ATAGAGAG,,
A3,,,,,,N703,AGGCAGAA,S502,ATAGAGAG,,
A4,,,,,,N704,TCCTGAGC,S502,ATAGAGAG,,
A5,,,,,,N705,GGACTCCT,S502,ATAGAGAG,,
A6,,,,,,N706,TAGGCATG,S502,ATAGAGAG,,
A7,,,,,,N707,CTCTCTAC,S502,ATAGAGAG,,
A8,,,,,,N710,CGAGGCTG,S502,ATAGAGAG,,
A9,,,,,,N711,AAGAGGCA,S502,ATAGAGAG,,
A10,,,,,,N712,GTAGAGGA,S502,ATAGAGAG,,
A11,,,,,,N714,GCTCATGA,S502,ATAGAGAG,,
A12,,,,,,N715,ATCTCAGG,S502,ATAGAGAG,,
B1,,,,,,N701,TAAGGCGA,S503,AGAGGATA,,
B2,,,,,,N702,CGTACTAG,S503,AGAGGATA,,
B3,,,,,,N703,AGGCAGAA,S503,AGAGGATA,,
B4,,,,,,N704,TCCTGAGC,S503,AGAGGATA,,
B5,,,,,,N705,GGACTCCT,S503,AGAGGATA,,
B6,,,,,,N706,TAGGCATG,S503,AGAGGATA,,
B7,,,,,,N707,CTCTCTAC,S503,AGAGGATA,,
B8,,,,,,N710,CGAGGCTG,S503,AGAGGATA,,
B9,,,,,,N711,AAGAGGCA,S503,AGAGGATA,,
B10,,,,,,N712,GTAGAGGA,S503,AGAGGATA,,
B11,,,,,,N714,GCTCATGA,S503,AGAGGATA,,
B12,,,,,,N715,ATCTCAGG,S503,AGAGGATA,,
C1,,,,,,N701,TAAGGCGA,S505,CTCCTTAC,,
C2,,,,,,N702,CGTACTAG,S505,CTCCTTAC,,
C3,,,,,,N703,AGGCAGAA,S505,CTCCTTAC,,
C4,,,,,,N704,TCCTGAGC,S505,CTCCTTAC,,
C5,,,,,,N705,GGACTCCT,S505,CTCCTTAC,,
C6,,,,,,N706,TAGGCATG,S505,CTCCTTAC,,
C7,,,,,,N707,CTCTCTAC,S505,CTCCTTAC,,
C8,,,,,,N710,CGAGGCTG,S505,CTCCTTAC,,
C9,,,,,,N711,AAGAGGCA,S505,CTCCTTAC,,
C10,,,,,,N712,GTAGAGGA,S505,CTCCTTAC,,
C11,,,,,,N714,GCTCATGA,S505,CTCCTTAC,,
C12,,,,,,N715,ATCTCAGG,S505,CTCCTTAC,,
D1,,,,,,N701,TAAGGCGA,S506,TATGCAGT,,
D2,,,,,,N702,CGTACTAG,S506,TATGCAGT,,
D3,,,,,,N703,AGGCAGAA,S506,TATGCAGT,,
D4,,,,,,N704,TCCTGAGC,S506,TATGCAGT,,
D5,,,,,,N705,GGACTCCT,S506,TATGCAGT,,
D6,,,,,,N706,TAGGCATG,S506,TATGCAGT,,
D7,,,,,,N707,CTCTCTAC,S506,TATGCAGT,,
D8,,,,,,N710,CGAGGCTG,S506,TATGCAGT,,
D9,,,,,,N711,AAGAGGCA,S506,TATGCAGT,,
D10,,,,,,N712,GTAGAGGA,S506,TATGCAGT,,
D11,,,,,,N714,GCTCATGA,S506,TATGCAGT,,
D12,,,,,,N715,ATCTCAGG,S506,TATGCAGT,,
E1,,,,,,N701,TAAGGCGA,S507,TACTCCTT,,
E2,,,,,,N702,CGTACTAG,S507,TACTCCTT,,
E3,,,,,,N703,AGGCAGAA,S507,TACTCCTT,,
E4,,,,,,N704,TCCTGAGC,S507,TACTCCTT,,
E5,,,,,,N705,GGACTCCT,S507,TACTCCTT,,
E6,,,,,,N706,TAGGCATG,S507,TACTCCTT,,
E7,,,,,,N707,CTCTCTAC,S507,TACTCCTT,,
E8,,,,,,N710,CGAGGCTG,S507,TACTCCTT,,
E9,,,,,,N711,AAGAGGCA,S507,TACTCCTT,,
E10,,,,,,N712,GTAGAGGA,S507,TACTCCTT,,
E11,,,,,,N714,GCTCATGA,S507,TACTCCTT,,
E12,,,,,,N715,ATCTCAGG,S507,TACTCCTT,,
F1,,,,,,N701,TAAGGCGA,S508,AGGCTTAG,,
F2,,,,,,N702,CGTACTAG,S508,AGGCTTAG,,
F3,,,,,,N703,AGGCAGAA,S508,AGGCTTAG,,
F4,,,,,,N704,TCCTGAGC,S508,AGGCTTAG,,
F5,,,,,,N705,GGACTCCT,S508,AGGCTTAG,,
F6,,,,,,N706,TAGGCATG,S508,AGGCTTAG,,
F7,,,,,,N707,CTCTCTAC,S508,AGGCTTAG,,
F8,,,,,,N710,CGAGGCTG,S508,AGGCTTAG,,
F9,,,,,,N711,AAGAGGCA,S508,AGGCTTAG,,
F10,,,,,,N712,GTAGAGGA,S508,AGGCTTAG,,
F11,,,,,,N714,GCTCATGA,S508,AGGCTTAG,,
F12,,,,,,N715,ATCTCAGG,S508,AGGCTTAG,,
G1,,,,,,N701,TAAGGCGA,S510,ATTAGACG,,
G2,,,,,,N702,CGTACTAG,S510,ATTAGACG,,
G3,,,,,,N703,AGGCAGAA,S510,ATTAGACG,,
G4,,,,,,N704,TCCTGAGC,S510,ATTAGACG,,
G5,,,,,,N705,GGACTCCT,S510,ATTAGACG,,
G6,,,,,,N706,TAGGCATG,S510,ATTAGACG,,
G7,,,,,,N707,CTCTCTAC,S510,ATTAGACG,,
G8,,,,,,N710,CGAGGCTG,S510,ATTAGACG,,
G9,,,,,,N711,AAGAGGCA,S510,ATTAGACG,,
G10,,,,,,N712,GTAGAGGA,S510,ATTAGACG,,
G11,,,,,,N714,GCTCATGA,S510,ATTAGACG,,
G12,,,,,,N715,ATCTCAGG,S510,ATTAGACG,,
H1,,,,,,N701,TAAGGCGA,S511,CGGAGAGA,,
H2,,,,,,N702,CGTACTAG,S511,CGGAGAGA,,
H3,,,,,,N703,AGGCAGAA,S511,CGGAGAGA,,
H4,,,,,,N704,TCCTGAGC,S511,CGGAGAGA,,
H5,,,,,,N705,GGACTCCT,S511,CGGAGAGA,,
H6,,,,,,N706,TAGGCATG,S511,CGGAGAGA,,
H7,,,,,,N707,CTCTCTAC,S511,CGGAGAGA,,
H8,,,,,,N710,CGAGGCTG,S511,CGGAGAGA,,
H9,,,,,,N711,AAGAGGCA,S511,CGGAGAGA,,
H10,,,,,,N712,GTAGAGGA,S511,CGGAGAGA,,
H11,,,,,,N714,GCTCATGA,S511,CGGAGAGA,,
H12,,,,,,N715,ATCTCAGG,S511,CGGAGAGA,,
```

Simply run `bcl2fastq` like this:

```console
bcl2fastq --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

After this, you will have `R1_001.fastq.gz` and `R2_001.fastq.gz` for each well:

```bash
A1_S1_R1_001.fastq.gz # 18 bp: 8 bp UMI + 10 bp RT barcode
A1_S1_R2_001.fastq.gz # 52 bp: cDNA
A2_S2_R1_001.fastq.gz # 18 bp: 8 bp UMI + 10 bp RT barcode
A2_S2_R2_001.fastq.gz # 52 bp: cDNA
...
...
...
H11_S95_R1_001.fastq.gz # 18 bp: 8 bp UMI + 10 bp RT barcode
H11_S95_R2_001.fastq.gz # 52 bp: cDNA
H12_S96_R1_001.fastq.gz # 18 bp: 8 bp UMI + 10 bp RT barcode
H12_S96_R2_001.fastq.gz # 52 bp: cDNA
```

That's it. You are ready to go from here using `starsolo`. You can and should treat each well as separate experiments, and single cell can be identified by the 10 bp RT barcode in `R1`. Each well needs to be processed independently as if they are from different experiments. For example, the RT barcode `GCGTTGGAGC` in the well __A1__ and the same RT barcode `GCGTTGGAGC` in the well __A2__ represent different cells. Therefore, we need to generate count matrix for each well separately, and combine them in the downstream analysis.

The advantage of this choice is that we actually divide each experiment into small chunks, and use the exact the same procedures for each chunk independently. In addition, the whitelist will simply be the 10 bp RT barcode for all the analysis. We will demonstrate this point in the later section using public data.

### Run bcl2fastq With Strategy 2

Alternatively, you can run `bcl2fastq` without a `SampleSheet.csv` with the `--create-fastq-for-index-reads` flag like this:

```console
bcl2fastq --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

In this case, you will have four `fastq` files per experiment:

```bash
Undetermined_S0_I1_001.fastq.gz  # 8 bp, i7
Undetermined_S0_I2_001.fastq.gz  # 8 bp, i5
Undetermined_S0_R1_001.fastq.gz  # 18 bp, 8 bp UMI + 10 bp RT barcode
Undetermined_S0_R2_001.fastq.gz  # 52 bp, cDNA
```

Sometimes, you want to process the whole experiments altogether. This choice is good for you if you want that. In this case, the cell barcodes will be the combination of the 10bp RT barcode, `i7`  and `i5`. Once we get those four `fastq` files, we need to put the cell barcode and UMI into the same `fastq` files so that we could use `starsolo`. To this end, we need to stitch `I1` (__i7__), `I2` (__i5__) and `R1` together like this:

```bash
paste <(zcat Undetermined_S0_R1_001.fastq.gz) \
      <(zcat Undetermined_S0_I1_001.fastq.gz) \
      <(zcat Undetermined_S0_I2_001.fastq.gz) | \
    awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2 $3} }' | \
    gzip > Undetermined_S0_CB_UMI.fastq.gz
```

After that, you are ready to throw `Undetermined_S0_R1_001.fastq.gz`, `Undetermined_S0_R2_001.fastq.gz` and `Undetermined_S0_CB_UMI.fastq.gz` into `starsolo`. See the later section.

## Public Data

For the purpose of demonstration, we will use the __sci-RNA-seq__ data from the following paper:

```{eval-rst}
.. note::

  Cao J, Packer JS, Ramani V, Cusanovich DA, Huynh C, Daza R, Qiu X, Lee C, Furlan SN, Steemers FJ, Adey A, Waterston RH, Trapnell C, Shendure J (2017) **Comprehensive single-cell transcriptional profiling of a multicellular organism.** *Science* 357:661–667. https://doi.org/10.1126/science.aam8940

```

where the authors developed __sci-RNA-seq__ for the first time and created a single cell atlas of an entire organism (C. elegans). The data is in GEO under the accession code [GSE98561](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98561). We are NOT going to use the C. elegans data. Instead, we will just use the __HEK293T, HeLa S3 and NIH/3T3 cell and nuclei mixture (96 x 96 sci-RNA-seq)__ sample data from the species mixing experiment.

As the sample name suggested, in this experiment the authors used 96 RT barcodes (1st plate) and 96 `i5 + i7` barcodes (2nd plate) for the experiment. You can go to [__SRX2784960__](https://www.ncbi.nlm.nih.gov/sra?term=SRX2784960) to see the data associated with the sample. As you can see, there are 96 run accessions, meaning that they already demultiplexed the experiment for us based on the well in the 2nd plate. Now, each `SRR` accession represent data from a well in the 2nd plate. We could get single cells by just take the 10 bp RT barcode in each of them. For the demonstration, I'm not going to do all 96 wells. Let's just use the first four wells:

```console
# get fastq files
mkdir -p sci-rna-seq/data
fastq-dump --split-files \
           --origfmt \
           --outdir sci-rna-seq/data \
           --defline-seq '@rd.$si:$sg:$sn' \
           SRR5509659 SRR5509660 SRR5509661 SRR5509662

# unfortunately, fastq-dump cannot generate gzipped file on the fly
# compress to save space
gzip sci-rna-seq/data/*.fastq
```

The reason of using the specific options above is a bit complicated. I wrote [__a post__](https://dbrg77.github.io/posts/2022-07-26-getting-index-reads-from-sra/) about this, and you can have a look to see the reason. Once the program finishes running, you will have two `fastq` per accession:

```
scg_prep_test/sci-rna-seq
└── data
    ├── SRR5509659_1.fastq.gz
    ├── SRR5509659_2.fastq.gz
    ├── SRR5509660_1.fastq.gz
    ├── SRR5509660_2.fastq.gz
    ├── SRR5509661_1.fastq.gz
    ├── SRR5509661_2.fastq.gz
    ├── SRR5509662_1.fastq.gz
    └── SRR5509662_2.fastq.gz

1 directory, 8 files
```

Now, let's have a look at the first 5 reads from `SRR5509659_1.fastq.gz`:

```text
@rd.1:TCGGATTCGG:1
GTTCTGGGATGGCGGATC
+1
AAAAAEEEEEEEEEEEEE
@rd.2:TCGGATTCGG:2
GTGTTATTACCAGCGCAG
+2
AAAAAEEEEEEEEEEEEE
@rd.3:TCGGATTCGG:3
TTTGAGTACCGCTGCTTC
+3
AAAAAEEEEEEEEEEEEE
@rd.4:TCGGATTCGG:4
GTTTGTTTCGGTCGTTAA
+4
6AAAAEEEEEEEEEEEAE
@rd.5:TCGGATTCGG:5
GTACACGCACCAGCGCAG
+5
AAAAAEEE/EEAEEEEEE
```

As you can see, those reads are 18 bp in length. The first 8 bp are UMI and the rest 10 bp are the RT barcode. If you look at the header, using `:` as the field separator, the second field is the well barcode, which should be `i7 + i5`. However, you don't have to do this. Anything that can index the individual well will work. In this case, you can see 10-bp sequence, which is basically the __i7__ barcode. This is because the PCR primer used by the authors have 96 different 10-bp __i7__ sequences to index each well in the second plate, so they don't really need the combination of __i5__ and __i7__. Just `i7` is enough. All reads in this particular well have the same __i7__ well barcode: `TCGGATTCGG`.

## Prepare Whitelist

To generate the whitelist, you need the 10-bp RT barcodes, the __i7__ and __i5__ indices. Generate a combination of them as the pool of all possible cell barcodes.

Unfortunately, in the [__sci-RNA-seq paper__](http://science.sciencemag.org/content/357/6352/661), I cannot seem to find the information of those oligos. However, in the [__sci-RNA-seq3 paper__](https://www.nature.com/articles/s41586-019-0969-x) which is an updated version of the original one, I can find 384 different 10-bp RT barcodes, 96 different 10-bp `i5` index and 96 different 10-bp `i7` index from the [Supplementary Table S11](https://teichlab.github.io/scg_lib_structs/data/41586_2019_969_MOESM3_ESM.xlsx) of the paper. The __sci-RNA-seq__ seem to use the same barcodes. We could collect the index sequences as tables as follows, and the names of the oligos are directly taken from the paper to be consistent (showing only 5 of the table to save space):

__RT Barcodes (10 bp)__

| Name             | Sequence   | Reverse complement |
|------------------|------------|:------------------:|
| sc_ligation_RT_1 | TCCTACCAGT |     ACTGGTAGGA     |
| sc_ligation_RT_2 | GCGTTGGAGC |     GCTCCAACGC     |
| sc_ligation_RT_3 | GATCTTACGC |     GCGTAAGATC     |
| sc_ligation_RT_4 | CTGATGGTCA |     TGACCATCAG     |
| sc_ligation_RT_5 | CCGAGAATCC |     GGATTCTCGG     |

__i7 Barcodes (10 bp)__

| Name | Sequence   | Reverse complement |
|------|------------|:------------------:|
| P7-1 | CCGAATCCGA |     TCGGATTCGG     |
| P7-2 | ATAAGCCGGA |     TCCGGCTTAT     |
| P7-3 | CCGGCGGCGA |     TCGCCGCCGG     |
| P7-4 | GGCTTGCCAA |     TTGGCAAGCC     |
| P7-5 | CCGCTAGCTG |     CAGCTAGCGG     |

__i5 Barcodes (10 bp)__

| Name | Sequence   | Reverse complement |
|------|------------|:------------------:|
| P5-1 | CTCCATCGAG |     CTCGATGGAG     |
| P5-2 | TTGGTAGTCG |     CGACTACCAA     |
| P5-3 | GGCCGTCAAC |     GTTGACGGCC     |
| P5-4 | CCTAGACGAG |     CTCGTCTAGG     |
| P5-5 | TCGTTAGAGC |     GCTCTAACGA     |

I have put those three tables into `csv` files and you can download them to have a look:

[sci-RNA-seq3_RT_bc.csv](https://teichlab.github.io/scg_lib_structs/data/sci-RNA-seq3_RT_bc.csv)  
[sci-RNA-seq3_p7.csv](https://teichlab.github.io/scg_lib_structs/data/sci-RNA-seq3_p7.csv)  
[sci-RNA-seq3_p5.csv](https://teichlab.github.io/scg_lib_structs/data/sci-RNA-seq3_p5.csv)  

Let's download them:

```console
wget -P sci-rna-seq/data \
    https://teichlab.github.io/scg_lib_structs/data/sci-RNA-seq3_RT_bc.csv \
    https://teichlab.github.io/scg_lib_structs/data/sci-RNA-seq3_p7.csv \
    https://teichlab.github.io/scg_lib_structs/data/sci-RNA-seq3_p5.csv
```

If you use the full capacity of those oligos, you could have a capacity of __384 * 96 * 96 = 3,538,944__ barcodes.

### Whitelist For Strategy 1

In this strategy, you are going to process the data for each well separately and independently, so you only need one whitelist for all wells. That is, the 10 bp RT barcodes. The RT barcode is in the same direction of the Illumina TruSeq Read 1 sequence, so we should take the sequences as they are:

```
tail -n +2 sci-rna-seq/data/sci-RNA-seq3_RT_bc.csv | \
    cut -f 2 -d, > sci-rna-seq/data/whitelist_strategy_1.txt
```

### Whitelist For Strategy 2

In this strategy, you are going to process the data for all wells in an experiment or multiple experiments. The cells will be identified by the combination of __RT barcode + i7 + i5__. The sequence of `i7` and `i5` depends on the primers you used. In this case for the public data, we only need the `i7`, because that is the index used to index each well. Therefore, we need to generate all combinations of __RT barcode + i7__ for this specific data set. Again, the RT barcode is in the same direction of the Illumina TruSeq Read 1 sequence, so we should take the sequences as they are. However, the `i7` index is always sequenced using the bottom strand as the template, so we need to take the reverse complement of the sequence. Check the [__sci-RNA-seq GitHub page__](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html) if you are still confused:

```bash
for x in $(tail -n +2 sci-rna-seq/data/sci-RNA-seq3_RT_bc.csv | cut -f 2 -d,); do
    for y in $(tail -n +2 sci-rna-seq/data/sci-RNA-seq3_p7.csv | cut -f 3 -d,); do
        echo "${x}${y}"
        done
    done > sci-rna-seq/data/whitelist_strategy_2.txt
```

## From FastQ To Count Matrix

### Strategy 1

It is relatively straightforward in this case. Using `SRR5509659` as an example:

```console
mkdir -p sci-rna-seq/star_outs/SRR5509659

# map and generate the count matrix

STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix sci-rna-seq/star_outs/SRR5509659/ \
     --readFilesIn sci-rna-seq/data/SRR5509659_2.fastq.gz sci-rna-seq/data/SRR5509659_1.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 9 --soloCBlen 10 --soloUMIstart 1 --soloUMIlen 8 \
     --soloCBwhitelist sci-rna-seq/data/whitelist_strategy_1.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

Once that is finished, you can do the exact the same thing with all the rest wells. In practice, you can do this via a loop or a pipeline. They can be run independently in parallel.

### Strategy 2

In this strategy, the cell barcodes are the combination of __RT barcode + i7__. We need to put the cell barcodes and UMI in one `fastq` files. The `i7` index is normally in the `I1_001.fastq.gz` file. However, this file is not available in [__SRA__](https://www.ncbi.nlm.nih.gov/sra). We could extract the `i7` sequence from the `fastq` header and generate fake quality strings to produce a fake `I1` file. Then we could stitch the `R1` and `I1` file to make a single `fastq` file with cell barcodes and UMI. Also, we need to merge all of them to mimic the output we directly get from `bcl2fastq`:

```bash
# merge R1 and R2
cat sci-rna-seq/data/*_1.fastq.gz > sci-rna-seq/data/Undetermined_S0_R1_001.fastq.gz
cat sci-rna-seq/data/*_2.fastq.gz > sci-rna-seq/data/Undetermined_S0_R2_001.fastq.gz

# generate fake I1 using the header of the fastq file
zcat sci-rna-seq/data/Undetermined_S0_R1_001.fastq.gz | \
    awk -F ':' '(NR%4==1){print $0 "\n" $2 "\n+\n" "IIIIIIIIII"}' | \
    gzip > sci-rna-seq/data/Undetermined_S0_I1_001.fastq.gz
```

Now, the three "Undetermined" `fastq` files are the starting point of strategy 2 as if you are sequencing on your own. Let's get the `CB_UMI` file as previously described:

```bash
paste <(zcat sci-rna-seq/data/Undetermined_S0_R1_001.fastq.gz) \
      <(zcat sci-rna-seq/data/Undetermined_S0_I1_001.fastq.gz) | \
    awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
    gzip > sci-rna-seq/data/Undetermined_S0_CB_UMI.fastq.gz
```

With that, we could start the mapping using similar settings to those used in strategy 1

```console
mkdir -p sci-rna-seq/star_outs/combined

STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix sci-rna-seq/star_outs/combined/ \
     --readFilesIn sci-rna-seq/data/Undetermined_S0_R2_001.fastq.gz sci-rna-seq/data/Undetermined_S0_CB_UMI.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 9 --soloCBlen 20 --soloUMIstart 1 --soloUMIlen 8 \
     --soloCBwhitelist sci-rna-seq/data/whitelist_strategy_2.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __sci-RNA-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mix_hg38_mm10/star_index`

>>> Pointing to the directory of the star index. The public data from the above paper was produced using the human and mouse mixture sample.

`--readFilesCommand zcat`

>>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix`

>>> We want to keep everything organised. This parameter directs all output files inside the specified directory.

`--readFilesIn`

>>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second file should contain cell barcode and UMI. In __sci-RNA-seq__, cDNA reads come from Read 2, and the cell barcode and UMI come from Read 1 or the `CB_UMI` file you just prepared. Check [the sci-RNA-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html) if you are not sure.

`--soloType CB_UMI_Simple`

>>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 9 --soloCBlen 10 or 20 --soloUMIstart 1 --soloUMIlen 8`

>>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist`

>>> The plain text file containing all possible valid cell barcodes, one per line. In the previous section, we prepared two versions of whitelist. The choice depends on the strategy you use.

`--soloCellFilter EmptyDrops_CR`

>>> Experiments are never perfect. Even for barcodes that do not capture the molecules inside the cells, you may still get some reads due to various reasons, such as ambient RNA or DNA and leakage. In general, the number of reads from those barcodes should be much smaller, often orders of magnitude smaller, than those barcodes from real cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

```{eval-rst}
.. important::

  You may want to take a closer look at the filtered results from **starsolo**. The count matrix in the ``Solo.out/Gene/filtered`` directory is the count matrix for "true" cells called by **starsolo**. In the **sci-RNA-seq** paper, it seems only 12 cells are sorted into each well in the 2nd plate. Therefore, we would expect to see ~12 cells in each accession in **Strategy 1**, and ~48 cells in total in **Strategy 2**. However, many more cells are called by **starsolo**. Maybe this ``EmptyDrops_CR`` parameter is **NOT** optimal for combinatorial indexing methods. Maybe we should try ``TopCells`` option. Open for discussion.
```

`--soloStrand Forward`

>>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __sci-RNA-seq__, the cDNA reads are from the Read 2 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing UMI and the Illumina Read 1 sequence. Therefore, Read 1 consists of RT barcodes and UMI. They come from the first strand, complementary to the coding strand. Read 2 comes from the coding strand. Therefore, use `Forward` for __sci-RNA-seq__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the sci-RNA-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html) if you are not sure.

`--outSAMattributes CB UB`

>>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/sci-rna-seq/
├── data
│   ├── sci-RNA-seq3_p5.csv
│   ├── sci-RNA-seq3_p7.csv
│   ├── sci-RNA-seq3_RT_bc.csv
│   ├── SRR5509659_1.fastq.gz
│   ├── SRR5509659_2.fastq.gz
│   ├── SRR5509660_1.fastq.gz
│   ├── SRR5509660_2.fastq.gz
│   ├── SRR5509661_1.fastq.gz
│   ├── SRR5509661_2.fastq.gz
│   ├── SRR5509662_1.fastq.gz
│   ├── SRR5509662_2.fastq.gz
│   ├── Undetermined_S0_CB_UMI.fastq.gz
│   ├── Undetermined_S0_I1_001.fastq.gz
│   ├── Undetermined_S0_R1_001.fastq.gz
│   ├── Undetermined_S0_R2_001.fastq.gz
│   ├── whitelist_strategy_1.txt
│   └── whitelist_strategy_2.txt
└── star_outs
    ├── combined
    │   ├── Aligned.sortedByCoord.out.bam
    │   ├── Log.final.out
    │   ├── Log.out
    │   ├── Log.progress.out
    │   ├── SJ.out.tab
    │   └── Solo.out
    │       ├── Barcodes.stats
    │       └── Gene
    │           ├── Features.stats
    │           ├── filtered
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── raw
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── Summary.csv
    │           └── UMIperCellSorted.txt
    ├── SRR5509659
    │   ├── Aligned.sortedByCoord.out.bam
    │   ├── Log.final.out
    │   ├── Log.out
    │   ├── Log.progress.out
    │   ├── SJ.out.tab
    │   └── Solo.out
    │       ├── Barcodes.stats
    │       └── Gene
    │           ├── Features.stats
    │           ├── filtered
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── raw
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── Summary.csv
    │           └── UMIperCellSorted.txt
    ├── SRR5509660
    │   ├── Aligned.sortedByCoord.out.bam
    │   ├── Log.final.out
    │   ├── Log.out
    │   ├── Log.progress.out
    │   ├── SJ.out.tab
    │   └── Solo.out
    │       ├── Barcodes.stats
    │       └── Gene
    │           ├── Features.stats
    │           ├── filtered
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── raw
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── Summary.csv
    │           └── UMIperCellSorted.txt
    ├── SRR5509661
    │   ├── Aligned.sortedByCoord.out.bam
    │   ├── Log.final.out
    │   ├── Log.out
    │   ├── Log.progress.out
    │   ├── SJ.out.tab
    │   └── Solo.out
    │       ├── Barcodes.stats
    │       └── Gene
    │           ├── Features.stats
    │           ├── filtered
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── raw
    │           │   ├── barcodes.tsv
    │           │   ├── features.tsv
    │           │   └── matrix.mtx
    │           ├── Summary.csv
    │           └── UMIperCellSorted.txt
    └── SRR5509662
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

27 directories, 92 files
```