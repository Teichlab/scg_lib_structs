# STRT-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/STRT-seq_family.html) to see how __STRT-seq__ libraries are generated experimentally. This is one of the few methods that sequence the 5' of the transcript. This is useful for both gene expression quantification and transcription start site (TSS) identification. The method is plate-based where single cells are sorted into each well in a one-cell-per-well manner.

```{eval-rst}
.. important::
  
  1. Be aware that there are three versions of **STRT-seq**, which are quite different in terms of experimental procedures and library structures. Make sure you check the `STRT-seq GitHub Page <https://teichlab.github.io/scg_lib_structs/methods_html/STRT-seq_family.html>`_ for some details.
  2. In this documentation, we will go through the procedures to process all versions just for the sake of record keeping.

```

## For Your Own Experiments

The read configuration varies greatly depending on the version.

### The Original Version

| Order | Read             | Cycle    | Description                                                       |
|-------|------------------|----------|-------------------------------------------------------------------|
| 1     | Read 1           | >=50     | This yields `R1_001.fastq.gz`, Cell barcodes + __GGG__ + cDNA     |
| 2     | Index 1 (__i7__) | 6 or 8   | This yields `I1_001.fastq.gz`, sample index                       |
| 3     | Index 2 (__i5__) | Optional | This yields `I2_001.fastq.gz`, not really used but can be present |
| 4     | Read 2           | >=50     | This yields `R2_001.fastq.gz`, cDNA                               |

__Read 1__ is the only required reads and the content is like this:

| Length | Sequence (5' -> 3')                        |
|--------|--------------------------------------------|
| >=50   |  6 bp __Cell barcodes__ + GGG + 5' of cDNA |

### The C1 Version

| Order | Read             | Cycle    | Description                                                                       |
|-------|------------------|----------|-----------------------------------------------------------------------------------|
| 1     | Read 1           | >=50     | This yields `R1_001.fastq.gz`, UMI + __GGG__ + cDNA                               |
| 2     | Index 1 (__i7__) | 8        | This yields `I1_001.fastq.gz`, Tn5 barcode which serves as the cell barcode index |
| 3     | Index 2 (__i5__) | Optional | This yields `I2_001.fastq.gz`, not really used but can be present                 |
| 4     | Read 2           | >=50     | This yields `R2_001.fastq.gz`, cDNA                                               |

The __Read 1__ and __Index 1__ are the only required reads, and the content of __Read 1__ is like this:

| Length  | Sequence (5' -> 3')             |
|---------|---------------------------------|
| >=50    | 5 bp __UMI__ + GGG + 5' of cDNA |

### The 2i Version

This configuration is more complicated and the naming of the output files does not really follow our normal convention. __DO NOT__ get confused.

| Order | Read             | Cycle    | Description                                              |
|-------|------------------|----------|----------------------------------------------------------|
| 1     | Read 1           | >=50     | This normally yields `R1_001.fastq.gz`, UMI + cDNA reads |
| 2     | Index 1          | 8        | This normally yields `I1_001.fastq.gz`, Subarray barcode |
| 3     | Index 2 (__i7__) | 5        | This normally yields `I2_001.fastq.gz`, Well barcode     |
| 4     | Read 2           | Optional | This normally yields `R2_001.fastq.gz`, cDNA reads       |

The content of __Read 1__ is like this:

| Length   | Sequence (5' -> 3')             |
|----------|---------------------------------|
| >=50     | 6 bp __UMI__ + GGG + 5' of cDNA |

In all cases, the pair-end sequencing mode can be used, but the original publications only used single-end reads. If you use this method, you have to sequence the library on your because custom sequencing primers are used, but that can be modified. You need to get the `fastq` files by running `bcl2fastq` by yourself. In the original version, it is better to write a `SampleSheet.csv` with `i7` indices for each sample. In the __C1__ and __2i__ versions, it is better just run `bcl2fastq` without a `SampleSheet.csv`. You will see the reason later. Here is an example of `SampleSheet.csv` of NextSeq runs with two samples using some standard index with the original version of __STRT-seq__:

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

You need to run `bcl2fastq` differently based on the versions like this:

```bash
# for the original version with the above SampleSheet.csv

bcl2fastq --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4

# for the C1 version without a SampleSheet.csv

bcl2fastq --use-bases-mask=Y50,I8,Y50 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4

# for the 2i version without a SampleSheet.csv

bcl2fastq --use-bases-mask=Y50,I8,I5,Y50 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is the `--use-bases-mask` flag in the __C1__ and __2i__ versions. Using the __2i__ version as an example, we have four reads in this case, and that parameter specifies how we treat each read in the stated order:

1. `Y50` at the first position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `I8` at the second position indicates "use the cycle as a real read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read, though it is the 2nd read overall.
3. `I5` at the third position indicates "use the cycle as an index read", so you will get 5-nt sequences, output as `I2_001.fastq.gz`, because this is the 2nd index read, though it is the 3rd read overall.
4. `Y50` at the fourth position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read, though it is the 4th read overall.

After that, you will get two files per sample for the original version, three files per run for the __C1__ version and four files per run for the __2i__ version:

```bash
# Original version
Sample01_S1_R1_001.fastq.gz # 50 bp: 6 bp cell barcodes + GGG + 5' cDNA
Sample01_S1_R2_001.fastq.gz # 50 bp: cDNA reads
Sample02_S2_R1_001.fastq.gz # 50 bp: 6 bp cell barcodes + GGG + 5' cDNA
Sample02_S2_R2_001.fastq.gz # 50 bp: cDNA reads

# C1 version
Undetermined_S0_R1_001.fastq.gz # 50 bp: 5 bp UMI + GGG + 5' cDNA
Undetermined_S0_I1_001.fastq.gz # 8 bp: cell barcodes
Undetermined_S0_R2_001.fastq.gz # 50 bp: cDNA reads

# 2i version
Undetermined_S0_R1_001.fastq.gz # 50 bp: 6bp UMI + GGG + 5' cDNA
Undetermined_S0_I1_001.fastq.gz # 8 bp: Subarray barcodes
Undetermined_S0_I2_001.fastq.gz # 5 bp: Well barcodes
Undetermined_S0_R2_001.fastq.gz # 50 bp: cDNA reads
```

There are no UMIs in the original version. For those types of data, we should demultiplex the `fastq` files based on the cell barcodes (the first 6 bp in __Read 1__), making one (for single-end) or two (for pair-end) files per cell. This can be achieved using any demultiplex programs, but we will use [cutadapt](https://cutadapt.readthedocs.io/en/stable/) as the demonstration later. For the __C1__ and __2i__ versions, the cell barcodes and UMI are distributed in different reads. We need to collect them into one `fastq` file in order to use `starsolo`. This can be done by simple stitching the reads:

```bash
# C1 version
paste <(zcat Undetermined_S0_I1_001.fastq.gz) \
      <(zcat Undetermined_S0_R1_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $2} else {print $1 $2} }' | \
      gzip > Undetermined_S0_CB_UMI.fastq.gz

# 2i version
paste <(zcat Undetermined_S0_I1_001.fastq.gz) \
      <(zcat Undetermined_S0_I2_001.fastq.gz) \
      <(zcat Undetermined_S0_R1_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $3} else {print $1 $2 $3} }' | \
      gzip > Undetermined_S0_CB_UMI_R1.fastq.gz
```

After that, you are ready to go.

## Public Data

For the purpose of demonstration, we are using the data from the following publications:

```{eval-rst}
.. note::
  
  **Original**

  Islam S, Kjällquist U, Moliner A, Zajac P, Fan J-B, Lönnerberg P, Linnarsson S (2011) **Characterization of the single-cell transcriptional landscape by highly multiplex RNA-seq.** *Genome Res* 21:1160–1167. https://doi.org/10.1101/gr.110882.110

  **C1**

  Islam S, Zeisel A, Joost S, Manno GL, Zajac P, Kasper M, Lönnerberg P, Linnarsson S (2014) **Quantitative single-cell RNA-seq with unique molecular identifiers.** *Nat Methods* 11:163–166. https://doi.org/10.1038/nmeth.2772

  **2i**

  Hochgerner H, Lönnerberg P, Hodge R, Mikes J, Heskol A, Hubschle H, Lin P, Picelli S, Manno GL, Ratz M, Dunne J, Husain S, Lein E, Srinivasan M, Zeisel A, Linnarsson S (2017) **STRT-seq-2i: dual-index 5ʹ single cell and nucleus RNA-seq on an addressable microwell array.** *Sci Rep*-uk 7:16327. https://doi.org/10.1038/s41598-017-16546-4

```

where the authors developed those methods for the first time.

### The Original Version

The raw data for the original version can be found from [__the PRJNA140307 ENA page__](https://www.ebi.ac.uk/ena/browser/view/PRJNA140307?show=reads). I have prepared the read information, and you can [__download here__](https://teichlab.github.io/scg_lib_structs/data/filereport_read_run_PRJNA140307.tsv). The authors already demultiplexed for us. They were using single-end sequencing mode, so there is one file per cell. To mimic what we get directly from the machine, we could merge all of them into one file.

```bash
# get individual fastq files and merge into one file
mkdir -p strt-seq/data
wget -P strt-seq/data https://teichlab.github.io/scg_lib_structs/data/filereport_read_run_PRJNA140307.tsv
wget -i <(cut -f 8 strt-seq/data/filereport_read_run_PRJNA140307.tsv | tail -n +2 | awk '{print "ftp://" $0}') \
     -O /dev/stdout >> trt-seq/data/STRT-seq.fastq.gz
```

Now we need to demultiplex the `fastq` file into individual files based on the first 6 bp. In this way, each cell has one file. Here, we use `cutadapt`. The cell barcode information can be found in this [__Supplementary Information__](https://teichlab.github.io/scg_lib_structs/data/STRT_GenomeRes_2011_SI.pdf) from the Genome Res. paper. We need the barcode in `fasta` format:

```
>bc01
TTTAGG
>bc02
ATTCCA
>bc03
GCTCAA
>bc04
CATCCC
>bc05
TTGGAC
. . .
```

I have already prepared the `fasta` file and you can [__download from here__](https://teichlab.github.io/scg_lib_structs/data/STRT_bc.fa), and pass the `fasta` to `cutadapt`:

```console
wget -P strt-seq/data https://teichlab.github.io/scg_lib_structs/data/STRT_bc.fa
cutadapt -j 4 -g ^file:strt-seq/data/STRT_bc.fa \
         --no-indels \
         -o "strt-seq/data/demul-{name}.fastq.gz" \
         strt-seq/data/STRT-seq.fastq.gz
```

It should finish without any problem, and we should have 97 more files under `strt-seq/data`. They are named as `demul-bc{01..96}.fastq.gz` and `demul-unknown.fastq.gz`. The size of the "unknown" file should be very small. We are ready to go from here for the original version.

### The C1 Version

The raw data for the C1 version can be found from [__the PRJNA203208 ENA page__](https://www.ebi.ac.uk/ena/browser/view/PRJNA203208?show=reads). I have prepared the read information as a TSV file including the barcode as the last column, and you can [__download here__](https://teichlab.github.io/scg_lib_structs/data/filereport_read_run_PRJNA203208.tsv). Again, the authors already demultiplexed for us. They were using single-end sequencing mode, so there is one file per cell.

```bash
mkdir -p strt-seq-c1/data
wget -P strt-seq-c1/data \
    https://teichlab.github.io/scg_lib_structs/data/filereport_read_run_PRJNA203208.tsv

# there are two types of libraries
# the one with the string "single" in the cell name is the regular one
# the one with the string "amplified" has extra 9 cycles of amplification
# we just use the regular one here
wget -P strt-seq-c1/data/ \
     -i <(tail -n +2 strt-seq-c1/data/filereport_read_run_PRJNA203208.tsv | grep '_single' | cut -f 8 | awk '{print "ftp://" $0}') 
```

Since the authors already demultiplexed the data to one file per cell, we need to add the cell barcode with fake quality scores in front of the reads and merge them into one `fastq` file. This mimics the `Undetermined_S0_CB_UMI_R1.fastq.gz` that we will get by ourselves. To this end, we do:

```bash
tail -n +2 strt-seq-c1/data/filereport_read_run_PRJNA203208.tsv | \
    grep '_single' | cut -f 4,10 | \
    while read -r line; do
        srr=$(echo "${line}" | cut -f 1)
        bc=$(echo "${line}" | cut -f 2)
        zcat strt-seq-c1/data/${srr}.fastq.gz | \
            awk -v BARCODE="${bc}" '{ if(NR%4==1||NR%4==3) {print $0} if(NR%4==2) {print BARCODE $0} if(NR%4==0) {print "IIIIIIII" $0} }' | \
            gzip >> strt-seq-c1/data/CB_UMI_R1.fastq.gz
    done
```

We are ready to go from here for the C1 version.

### The 2i Version

The raw data for the 2i version can be found from [__the PRJNA394919 ENA page__](https://www.ebi.ac.uk/ena/browser/view/PRJNA394919?show=reads). To preprocess the data, we need the oligo sequences. I have not got this information from the paper. I will update once I get them.

## Prepare Whitelist

### The Original Version

There are no UMIs in this version, and each cell has been demultiplexed into individual files. Therefore, we do not need a whitelist, but we do need to prepare a manifest, pointing the files to `starsolo`:

```bash
for i in $(ls strt-seq/data/demul-bc*.gz); do
    cell=$(echo ${i} | cut -f 3 -d '-')
    echo -e "${i}\t-\t${cell%.fastq.gz}"
done > islam2011_manifest.tsv
```

### The C1 Version

In this version, cDNA from individual cells are tagmented by barcoded Tn5 separately. The Tn5 barcode serves as the cell barcode. You can find the full sequence from the [__Supplementary Table 2__](https://teichlab.github.io/scg_lib_structs/data/41592_2014_BFnmeth2772_MOESM268_ESM.xlsx) from the Isalm2014 paper in Nature Methods. There are 96 different 8-bp Tn5 barcodes:

| Name      | Sequence | Reverse complement |
|-----------|----------|--------------------|
| C1-TN5-1  | CGTCTAAT | ATTAGACG           |
| C1-TN5-2  | AGACTCGT | ACGAGTCT           |
| C1-TN5-3  | GCACGTCA | TGACGTGC           |
| C1-TN5-4  | TCAACGAC | GTCGTTGA           |
| C1-TN5-5  | ATTTAGCG | CGCTAAAT           |
| C1-TN5-6  | ATACAGAC | GTCTGTAT           |
| C1-TN5-7  | TGCGTAGG | CCTACGCA           |
| C1-TN5-8  | TGGAGCTC | GAGCTCCA           |
| C1-TN5-9  | TGAATACC | GGTATTCA           |
| C1-TN5-10 | TCTCACAC | GTGTGAGA           |
| C1-TN5-11 | TACTGGTA | TACCAGTA           |
| C1-TN5-12 | ACGATAGG | CCTATCGT           |
| C1-TN5-13 | GATGTCGA | TCGACATC           |
| C1-TN5-14 | TTACGGGT | ACCCGTAA           |
| C1-TN5-15 | CACAGCAT | ATGCTGTG           |
| C1-TN5-16 | CTTTGACA | TGTCAAAG           |
| C1-TN5-17 | CCTTCAAG | CTTGAAGG           |
| C1-TN5-18 | GAGTCCTG | CAGGACTC           |
| C1-TN5-19 | CACACTGA | TCAGTGTG           |
| C1-TN5-20 | GTTACAGG | CCTGTAAC           |
| C1-TN5-21 | GGACCTTT | AAAGGTCC           |
| C1-TN5-22 | TTCCGTTC | GAACGGAA           |
| C1-TN5-23 | ACTGTTTG | CAAACAGT           |
| C1-TN5-24 | AAGTGGCT | AGCCACTT           |
| C1-TN5-25 | CTGTACAA | TTGTACAG           |
| C1-TN5-26 | CGCAAAGT | ACTTTGCG           |
| C1-TN5-27 | GTGCATGA | TCATGCAC           |
| C1-TN5-28 | GTCATTAG | CTAATGAC           |
| C1-TN5-29 | AGCTCCTT | AAGGAGCT           |
| C1-TN5-30 | TCACCCGA | TCGGGTGA           |
| C1-TN5-31 | GTTGCCAC | GTGGCAAC           |
| C1-TN5-32 | TGTACCAA | TTGGTACA           |
| C1-TN5-33 | AACGAGGT | ACCTCGTT           |
| C1-TN5-34 | AGCCACCA | TGGTGGCT           |
| C1-TN5-35 | GGTAATCA | TGATTACC           |
| C1-TN5-36 | CCAGTCCA | TGGACTGG           |
| C1-TN5-37 | ACCTCAGC | GCTGAGGT           |
| C1-TN5-38 | GGTGGACT | AGTCCACC           |
| C1-TN5-39 | GACAAACC | GGTTTGTC           |
| C1-TN5-40 | TAACTCCG | CGGAGTTA           |
| C1-TN5-41 | ACACCGTG | CACGGTGT           |
| C1-TN5-42 | GTAGAACG | CGTTCTAC           |
| C1-TN5-43 | GGATTGAC | GTCAATCC           |
| C1-TN5-44 | ACGTATCC | GGATACGT           |
| C1-TN5-45 | TTCGGAAA | TTTCCGAA           |
| C1-TN5-46 | AGTTGTGT | ACACAACT           |
| C1-TN5-47 | AAGCACAT | ATGTGCTT           |
| C1-TN5-48 | CTGTCATT | AATGACAG           |
| C1-TN5-49 | GTCCTATA | TATAGGAC           |
| C1-TN5-50 | CTACGCTG | CAGCGTAG           |
| C1-TN5-51 | GGGATTGT | ACAATCCC           |
| C1-TN5-52 | TGATGTAG | CTACATCA           |
| C1-TN5-53 | TTCGCTGT | ACAGCGAA           |
| C1-TN5-54 | GAAGACTT | AAGTCTTC           |
| C1-TN5-55 | TCTGGGCA | TGCCCAGA           |
| C1-TN5-56 | CAACTAGA | TCTAGTTG           |
| C1-TN5-57 | CCATGGGA | TCCCATGG           |
| C1-TN5-58 | ATGCGACG | CGTCGCAT           |
| C1-TN5-59 | GAGGGTAG | CTACCCTC           |
| C1-TN5-60 | CGGGTGAA | TTCACCCG           |
| C1-TN5-61 | GCCATCTT | AAGATGGC           |
| C1-TN5-62 | GCATAATC | GATTATGC           |
| C1-TN5-63 | TCTATGGT | ACCATAGA           |
| C1-TN5-64 | AGGACTTA | TAAGTCCT           |
| C1-TN5-65 | CGTGATTC | GAATCACG           |
| C1-TN5-66 | ACTAGCGA | TCGCTAGT           |
| C1-TN5-67 | GTAACTCC | GGAGTTAC           |
| C1-TN5-68 | CGGAAGTG | CACTTCCG           |
| C1-TN5-69 | CCGAGTAC | GTACTCGG           |
| C1-TN5-70 | GACGCAAT | ATTGCGTC           |
| C1-TN5-71 | ACCTGGAG | CTCCAGGT           |
| C1-TN5-72 | CATGGGTT | AACCCATG           |
| C1-TN5-73 | ATTCCTAG | CTAGGAAT           |
| C1-TN5-74 | AATCATGC | GCATGATT           |
| C1-TN5-75 | GCTTCCCT | AGGGAAGC           |
| C1-TN5-76 | AGGTAAAG | CTTTACCT           |
| C1-TN5-77 | CCACAACT | AGTTGTGG           |
| C1-TN5-78 | ACAGGCAT | ATGCCTGT           |
| C1-TN5-79 | TTTGTGTC | GACACAAA           |
| C1-TN5-80 | TGAGCATA | TATGCTCA           |
| C1-TN5-81 | TTAGACGC | GCGTCTAA           |
| C1-TN5-82 | CGCTTGCT | AGCAAGCG           |
| C1-TN5-83 | AGTCTGCC | GGCAGACT           |
| C1-TN5-84 | CATAGTCG | CGACTATG           |
| C1-TN5-85 | TCTTGCTG | CAGCAAGA           |
| C1-TN5-86 | GGGACAAC | GTTGTCCC           |
| C1-TN5-87 | ATATTCCC | GGGAATAT           |
| C1-TN5-88 | TGTTAAGC | GCTTAACA           |
| C1-TN5-89 | TACGCCTC | GAGGCGTA           |
| C1-TN5-90 | CACTTATC | GATAAGTG           |
| C1-TN5-91 | ACCGCTAA | TTAGCGGT           |
| C1-TN5-92 | TAAGGTCC | GGACCTTA           |
| C1-TN5-93 | GAAAGGTG | CACCTTTC           |
| C1-TN5-94 | ACGTTGTA | TACAACGT           |
| C1-TN5-95 | GCAGAGAA | TTCTCTGC           |
| C1-TN5-96 | GCATTTGG | CCAAATGC           |

I have prepared the full tables in `csv` format for you to download:

[STRT-seq_C1_bc.csv](https://teichlab.github.io/scg_lib_structs/data/STRT-seq_C1_bc.csv)  

If we check carefully about the oligo orientation in the [__STRT-seq C1 GitHub page__](https://teichlab.github.io/scg_lib_structs/methods_html/STRT-seq_family.html#STRT-seq-C1), we can see that the Tn5 barcodes are sequenced using the bottom strand as the template. Therefore, the barcode reads are actually reverse complement to the primer sequence. We should use the reverse complement as the whitelist:

```console
wget -P strt-seq-c1/data \
    https://teichlab.github.io/scg_lib_structs/data/STRT-seq_C1_bc.csv

tail -n +2 strt-seq-c1/data/STRT-seq_C1_bc.csv | \
    cut -f 3 -d, > strt-seq-c1/data/whitelist.txt
```

### The 2i Version

From the __Table 1__ of the Hochgerner2017 in Scientific Reports, there should be 32 different well barcodes (`DI-P1A-idx[1–32]-P1B`) and 96 different subarray barcodes (`STRT-Tn5-Idx[1–96]`). The cell barcodes are basically the combination of the subarray and well barcodes. Therefore, we should generate all combinations of the 96 subarray barcodes and 32 well barcodes for a total of __96 x 32 = 3072__ barcodes as whitelist. However, the sequences are not available from the paper. I will update once I get them.

## From FastQ To Count Matrix

Since we have already generated the manifest for the original version and the whitelist for the C1 version, it is now very easy to just run `starsolo`:

```console
# for the original version

STAR --runThreadN 4 \
     --genomeDir mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix strt-seq/star_outs/ \
     --readFilesManifest islam2011_manifest.tsv \
     --soloType SmartSeq \
     --clip5pNbases 3 \
     --soloUMIdedup Exact NoDedup \
     --soloStrand Forward \
     --outSAMtype BAM SortedByCoordinate

# for the C1 version

STAR --runThreadN 4 \
     --genomeDir mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix strt-seq-c1/star_outs/ \
     --readFilesIn strt-seq-c1/data/CB_UMI_R1.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 5 \
     --soloBarcodeMate 1 \
     --clip5pNbases 16 \
     --soloCBwhitelist strt-seq-c1/data/whitelist.txt \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __STRT-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/STRT-seq_family.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mm10/star_index`

>> Pointing to the directory of the star index. The public data from the above paper was from mouse embryonic stem cells (mESC).

`--readFilesCommand zcat`

>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix`

>> We want to keep everything organised. This parameter directs all output files into the `star_outs` directory within each method.

`--readFilesManifest` and `--readFilesIn`

>> For the original version, we need to provide the manifest here. For the C1 version, we provide the prepared read files containing cell barcodes, UMIs and the 5' of cDNA.

`--soloType`

>> The original version has no UMIs, and each cell has its own file, It is in the same situation of __SMART-seq__, so we put `SmartSeq` here. For the C1 version, we have prepared the files with cell barcodes and UMIs, so we use `CB_UMI_Simple` here.

`--soloCBstart 1 --soloCBlen 8 --soloUMIstart 9 --soloUMIlen 5`

>> This is for the C1 version. The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloBarcodeMate 1`

>> This is for the C1 version. This option is designed for the 5' sequencing methods, where one of the read contains not only cell barcodes + UMI, but useful cDNA as well. It tells the program that cell barcodes + UMI are in the first file in `--readFilesIn`. In this case, the public data is in single-end mode, so we only have one file.

`--clip5pNbases`

>> This option remove certain number of bases from the 5' of the read. In the original version, the cell barcodes are removed during the `cutadapt` demultiplexing step, but there are still a __GGG__ at the 5'. We need to ignore that. In the C1 version, the 5' of the read is 8 bp cell barcodes, 5 bp __UMI__ and __GGG__. Therefore, we need to remove __8 + 5 + 3 = 16 bp__.

`--soloUMIdedup Exact NoDedup`

The original version does not have UMI in the reads. `Exact` means perform the deduplication using the genomic coordinates, that is, fragments with the exact same starts and ends will be treated as duplicates. `NoDedup` means do not perform deduplication. In ChIP-seq, deduplication is standard. In non-UMI RNA-seq, it seems deduplication is not always enforced (I might be wrong). I'm not sure if this makes a huge difference. Putting both options here will generated two versions of count matrices, one with and one without deduplication.

`--soloCBwhitelist`

>> The plain text file containing all possible valid cell barcodes, one per line. We have prepared this file in the previous section. This is for the C1 version.

`--soloStrand Forward`

>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In all versions of __STRT-seq__, the cDNA reads from the Read 1 file are in the same direction of the mRNA, i.e. the coding strand. Therefore, use `Forward` for all __STRT-seq__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the STRT-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/STRT-seq_family.html) if you are not sure.

`--outSAMattributes CB UB`

>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
# The Original Version
scg_prep_test/strt-seq
├── data
│   ├── demul-bc01.fastq.gz
│   ├── demul-bc02.fastq.gz
│   ├── demul-bc03.fastq.gz
│   ├── demul-bc04.fastq.gz
│   ├── demul-bc05.fastq.gz
│   ├── demul-bc06.fastq.gz
│   ├── demul-bc07.fastq.gz
│   ├── demul-bc08.fastq.gz
│   ├── demul-bc09.fastq.gz
│   ├── demul-bc10.fastq.gz
│   ├── demul-bc11.fastq.gz
│   ├── demul-bc12.fastq.gz
│   ├── demul-bc13.fastq.gz
│   ├── demul-bc14.fastq.gz
│   ├── demul-bc15.fastq.gz
│   ├── demul-bc16.fastq.gz
│   ├── demul-bc17.fastq.gz
│   ├── demul-bc18.fastq.gz
│   ├── demul-bc19.fastq.gz
│   ├── demul-bc20.fastq.gz
│   ├── demul-bc21.fastq.gz
│   ├── demul-bc22.fastq.gz
│   ├── demul-bc23.fastq.gz
│   ├── demul-bc24.fastq.gz
│   ├── demul-bc25.fastq.gz
│   ├── demul-bc26.fastq.gz
│   ├── demul-bc27.fastq.gz
│   ├── demul-bc28.fastq.gz
│   ├── demul-bc29.fastq.gz
│   ├── demul-bc30.fastq.gz
│   ├── demul-bc31.fastq.gz
│   ├── demul-bc32.fastq.gz
│   ├── demul-bc33.fastq.gz
│   ├── demul-bc34.fastq.gz
│   ├── demul-bc35.fastq.gz
│   ├── demul-bc36.fastq.gz
│   ├── demul-bc37.fastq.gz
│   ├── demul-bc38.fastq.gz
│   ├── demul-bc39.fastq.gz
│   ├── demul-bc40.fastq.gz
│   ├── demul-bc41.fastq.gz
│   ├── demul-bc42.fastq.gz
│   ├── demul-bc43.fastq.gz
│   ├── demul-bc44.fastq.gz
│   ├── demul-bc45.fastq.gz
│   ├── demul-bc46.fastq.gz
│   ├── demul-bc47.fastq.gz
│   ├── demul-bc48.fastq.gz
│   ├── demul-bc49.fastq.gz
│   ├── demul-bc50.fastq.gz
│   ├── demul-bc51.fastq.gz
│   ├── demul-bc52.fastq.gz
│   ├── demul-bc53.fastq.gz
│   ├── demul-bc54.fastq.gz
│   ├── demul-bc55.fastq.gz
│   ├── demul-bc56.fastq.gz
│   ├── demul-bc57.fastq.gz
│   ├── demul-bc58.fastq.gz
│   ├── demul-bc59.fastq.gz
│   ├── demul-bc60.fastq.gz
│   ├── demul-bc61.fastq.gz
│   ├── demul-bc62.fastq.gz
│   ├── demul-bc63.fastq.gz
│   ├── demul-bc64.fastq.gz
│   ├── demul-bc65.fastq.gz
│   ├── demul-bc66.fastq.gz
│   ├── demul-bc67.fastq.gz
│   ├── demul-bc68.fastq.gz
│   ├── demul-bc69.fastq.gz
│   ├── demul-bc70.fastq.gz
│   ├── demul-bc71.fastq.gz
│   ├── demul-bc72.fastq.gz
│   ├── demul-bc73.fastq.gz
│   ├── demul-bc74.fastq.gz
│   ├── demul-bc75.fastq.gz
│   ├── demul-bc76.fastq.gz
│   ├── demul-bc77.fastq.gz
│   ├── demul-bc78.fastq.gz
│   ├── demul-bc79.fastq.gz
│   ├── demul-bc80.fastq.gz
│   ├── demul-bc81.fastq.gz
│   ├── demul-bc82.fastq.gz
│   ├── demul-bc83.fastq.gz
│   ├── demul-bc84.fastq.gz
│   ├── demul-bc85.fastq.gz
│   ├── demul-bc86.fastq.gz
│   ├── demul-bc87.fastq.gz
│   ├── demul-bc88.fastq.gz
│   ├── demul-bc89.fastq.gz
│   ├── demul-bc90.fastq.gz
│   ├── demul-bc91.fastq.gz
│   ├── demul-bc92.fastq.gz
│   ├── demul-bc93.fastq.gz
│   ├── demul-bc94.fastq.gz
│   ├── demul-bc95.fastq.gz
│   ├── demul-bc96.fastq.gz
│   ├── demul-unknown.fastq.gz
│   ├── STRT_bc.fa
│   └── STRT-seq.fastq.gz
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
            │   └── umiDedup-Exact.mtx
            ├── raw
            │   ├── barcodes.tsv
            │   ├── features.tsv
            │   ├── umiDedup-Exact.mtx
            │   └── umiDedup-NoDedup.mtx
            ├── Summary.csv
            └── UMIperCellSorted.txt

6 directories, 115 files

# The C1 Version
scg_prep_test/strt-seq-c1/
├── data
│   ├── CB_UMI_R1.fastq.gz
│   ├── filereport_read_run_PRJNA203208.tsv
│   ├── SRR1043197.fastq.gz
│   ├── SRR1043198.fastq.gz
│   ├── SRR1043199.fastq.gz
│   ├── SRR1043200.fastq.gz
│   ├── SRR1043201.fastq.gz
│   ├── SRR1043202.fastq.gz
│   ├── SRR1043203.fastq.gz
│   ├── SRR1043204.fastq.gz
│   ├── SRR1043205.fastq.gz
│   ├── SRR1043206.fastq.gz
│   ├── SRR1043207.fastq.gz
│   ├── SRR1043208.fastq.gz
│   ├── SRR1043209.fastq.gz
│   ├── SRR1043210.fastq.gz
│   ├── SRR1043211.fastq.gz
│   ├── SRR1043212.fastq.gz
│   ├── SRR1043213.fastq.gz
│   ├── SRR1043214.fastq.gz
│   ├── SRR1043215.fastq.gz
│   ├── SRR1043216.fastq.gz
│   ├── SRR1043217.fastq.gz
│   ├── SRR1043218.fastq.gz
│   ├── SRR1043219.fastq.gz
│   ├── SRR1043220.fastq.gz
│   ├── SRR1043221.fastq.gz
│   ├── SRR1043222.fastq.gz
│   ├── SRR1043223.fastq.gz
│   ├── SRR1043224.fastq.gz
│   ├── SRR1043225.fastq.gz
│   ├── SRR1043226.fastq.gz
│   ├── SRR1043227.fastq.gz
│   ├── SRR1043228.fastq.gz
│   ├── SRR1043229.fastq.gz
│   ├── SRR1043230.fastq.gz
│   ├── SRR1043231.fastq.gz
│   ├── SRR1043232.fastq.gz
│   ├── SRR1043233.fastq.gz
│   ├── SRR1043234.fastq.gz
│   ├── SRR1043235.fastq.gz
│   ├── SRR1043236.fastq.gz
│   ├── SRR1043237.fastq.gz
│   ├── SRR1043238.fastq.gz
│   ├── SRR1043239.fastq.gz
│   ├── SRR1043240.fastq.gz
│   ├── SRR1043241.fastq.gz
│   ├── SRR1043242.fastq.gz
│   ├── SRR1043243.fastq.gz
│   ├── SRR1043244.fastq.gz
│   ├── SRR1043245.fastq.gz
│   ├── SRR1043246.fastq.gz
│   ├── SRR1043247.fastq.gz
│   ├── SRR1043248.fastq.gz
│   ├── SRR1043249.fastq.gz
│   ├── SRR1043250.fastq.gz
│   ├── SRR1043251.fastq.gz
│   ├── SRR1043252.fastq.gz
│   ├── SRR1043253.fastq.gz
│   ├── SRR1043254.fastq.gz
│   ├── SRR1043255.fastq.gz
│   ├── SRR1043256.fastq.gz
│   ├── SRR1043257.fastq.gz
│   ├── SRR1043258.fastq.gz
│   ├── SRR1043259.fastq.gz
│   ├── SRR1043260.fastq.gz
│   ├── SRR1043261.fastq.gz
│   ├── SRR1043262.fastq.gz
│   ├── SRR1043263.fastq.gz
│   ├── SRR1043264.fastq.gz
│   ├── SRR1043265.fastq.gz
│   ├── SRR1043266.fastq.gz
│   ├── SRR1043267.fastq.gz
│   ├── SRR1043268.fastq.gz
│   ├── SRR1043269.fastq.gz
│   ├── SRR1043270.fastq.gz
│   ├── SRR1043271.fastq.gz
│   ├── SRR1043272.fastq.gz
│   ├── SRR1043273.fastq.gz
│   ├── SRR1043274.fastq.gz
│   ├── SRR1043275.fastq.gz
│   ├── SRR1043276.fastq.gz
│   ├── SRR1043277.fastq.gz
│   ├── SRR1043278.fastq.gz
│   ├── SRR1043279.fastq.gz
│   ├── SRR1043280.fastq.gz
│   ├── SRR1043281.fastq.gz
│   ├── SRR1043282.fastq.gz
│   ├── SRR1043283.fastq.gz
│   ├── SRR1043284.fastq.gz
│   ├── SRR1043285.fastq.gz
│   ├── SRR1043286.fastq.gz
│   ├── SRR1043287.fastq.gz
│   ├── SRR1043288.fastq.gz
│   ├── SRR1043289.fastq.gz
│   ├── SRR1043290.fastq.gz
│   ├── SRR1043291.fastq.gz
│   ├── SRR1043292.fastq.gz
│   ├── SRR1043293.fastq.gz
│   ├── SRR1043294.fastq.gz
│   ├── SRR1043295.fastq.gz
│   ├── SRR1043296.fastq.gz
│   ├── SRR1043297.fastq.gz
│   ├── SRR1043298.fastq.gz
│   ├── SRR1043299.fastq.gz
│   ├── SRR1043300.fastq.gz
│   ├── SRR1043301.fastq.gz
│   ├── SRR1043302.fastq.gz
│   ├── SRR1043303.fastq.gz
│   ├── SRR1043304.fastq.gz
│   ├── SRR1043305.fastq.gz
│   ├── SRR1043306.fastq.gz
│   ├── SRR1043307.fastq.gz
│   ├── SRR1043308.fastq.gz
│   ├── SRR1043309.fastq.gz
│   ├── SRR1043310.fastq.gz
│   ├── SRR1043311.fastq.gz
│   ├── SRR1043312.fastq.gz
│   ├── SRR1043313.fastq.gz
│   ├── SRR1043314.fastq.gz
│   ├── SRR1043315.fastq.gz
│   ├── SRR1043316.fastq.gz
│   ├── SRR1043317.fastq.gz
│   ├── SRR1043318.fastq.gz
│   ├── SRR1043319.fastq.gz
│   ├── SRR1043320.fastq.gz
│   ├── SRR1043321.fastq.gz
│   ├── SRR1043322.fastq.gz
│   ├── SRR1043323.fastq.gz
│   ├── SRR1043324.fastq.gz
│   ├── SRR1043325.fastq.gz
│   ├── SRR1043326.fastq.gz
│   ├── SRR1043327.fastq.gz
│   ├── SRR1043328.fastq.gz
│   ├── SRR1043329.fastq.gz
│   ├── SRR1043330.fastq.gz
│   ├── SRR1043331.fastq.gz
│   ├── SRR1043332.fastq.gz
│   ├── SRR1043333.fastq.gz
│   ├── SRR1043334.fastq.gz
│   ├── SRR1043335.fastq.gz
│   ├── SRR1043336.fastq.gz
│   ├── SRR1043337.fastq.gz
│   ├── SRR1043338.fastq.gz
│   ├── SRR1043339.fastq.gz
│   ├── SRR1043340.fastq.gz
│   ├── SRR1043341.fastq.gz
│   ├── SRR1043342.fastq.gz
│   ├── SRR1043343.fastq.gz
│   ├── SRR1043344.fastq.gz
│   ├── SRR1043345.fastq.gz
│   ├── SRR1043346.fastq.gz
│   ├── SRR1043347.fastq.gz
│   ├── SRR1043348.fastq.gz
│   ├── SRR1043349.fastq.gz
│   ├── SRR1043350.fastq.gz
│   ├── SRR1043351.fastq.gz
│   ├── SRR1043352.fastq.gz
│   ├── SRR1043353.fastq.gz
│   ├── SRR1043354.fastq.gz
│   ├── SRR1043355.fastq.gz
│   ├── SRR1043356.fastq.gz
│   ├── SRR1043357.fastq.gz
│   ├── SRR1043358.fastq.gz
│   ├── SRR1043359.fastq.gz
│   ├── SRR1043360.fastq.gz
│   ├── SRR1043361.fastq.gz
│   ├── SRR1043362.fastq.gz
│   ├── SRR1043363.fastq.gz
│   ├── SRR1043364.fastq.gz
│   ├── SRR1043365.fastq.gz
│   ├── SRR1043366.fastq.gz
│   ├── SRR1043367.fastq.gz
│   ├── SRR1043368.fastq.gz
│   ├── SRR1043369.fastq.gz
│   ├── SRR1043370.fastq.gz
│   ├── SRR1043371.fastq.gz
│   ├── SRR1043372.fastq.gz
│   ├── SRR1043373.fastq.gz
│   ├── SRR1043374.fastq.gz
│   ├── SRR1043375.fastq.gz
│   ├── SRR1043376.fastq.gz
│   ├── SRR1043377.fastq.gz
│   ├── SRR1043378.fastq.gz
│   ├── SRR1043379.fastq.gz
│   ├── SRR1043380.fastq.gz
│   ├── SRR1043381.fastq.gz
│   ├── SRR1043382.fastq.gz
│   ├── SRR1043383.fastq.gz
│   ├── SRR1043384.fastq.gz
│   ├── SRR1043385.fastq.gz
│   ├── SRR1043386.fastq.gz
│   ├── SRR1043387.fastq.gz
│   ├── SRR1043388.fastq.gz
│   ├── SRR1043389.fastq.gz
│   ├── SRR1043390.fastq.gz
│   ├── SRR1043391.fastq.gz
│   ├── SRR1043392.fastq.gz
│   ├── SRR1043393.fastq.gz
│   ├── SRR1043394.fastq.gz
│   ├── SRR1043395.fastq.gz
│   ├── SRR1043396.fastq.gz
│   ├── SRR1043397.fastq.gz
│   ├── SRR1043398.fastq.gz
│   ├── SRR1043399.fastq.gz
│   ├── SRR1043400.fastq.gz
│   ├── SRR1043401.fastq.gz
│   ├── SRR1043402.fastq.gz
│   ├── SRR1043403.fastq.gz
│   ├── SRR1043404.fastq.gz
│   ├── SRR1043405.fastq.gz
│   ├── SRR1043406.fastq.gz
│   ├── SRR1043407.fastq.gz
│   ├── SRR1043408.fastq.gz
│   ├── SRR1043409.fastq.gz
│   ├── SRR1043410.fastq.gz
│   ├── SRR1043411.fastq.gz
│   ├── SRR1043412.fastq.gz
│   ├── SRR1043413.fastq.gz
│   ├── SRR1043414.fastq.gz
│   ├── SRR1043415.fastq.gz
│   ├── SRR1043416.fastq.gz
│   ├── SRR1043417.fastq.gz
│   ├── SRR1043418.fastq.gz
│   ├── SRR1043419.fastq.gz
│   ├── SRR1043420.fastq.gz
│   ├── SRR1043421.fastq.gz
│   ├── SRR1043422.fastq.gz
│   ├── SRR1043423.fastq.gz
│   ├── SRR1043424.fastq.gz
│   ├── SRR1043425.fastq.gz
│   ├── SRR1043426.fastq.gz
│   ├── SRR1043427.fastq.gz
│   ├── SRR1043428.fastq.gz
│   ├── SRR1043429.fastq.gz
│   ├── SRR1043430.fastq.gz
│   ├── SRR1043431.fastq.gz
│   ├── SRR1043432.fastq.gz
│   ├── SRR1043433.fastq.gz
│   ├── SRR1043434.fastq.gz
│   ├── SRR1043435.fastq.gz
│   ├── SRR1043436.fastq.gz
│   ├── SRR1043437.fastq.gz
│   ├── SRR1043438.fastq.gz
│   ├── SRR1043439.fastq.gz
│   ├── SRR1043440.fastq.gz
│   ├── SRR1043441.fastq.gz
│   ├── SRR1043442.fastq.gz
│   ├── SRR1043443.fastq.gz
│   ├── SRR1043444.fastq.gz
│   ├── SRR1043445.fastq.gz
│   ├── SRR1043446.fastq.gz
│   ├── SRR1043447.fastq.gz
│   ├── SRR1043448.fastq.gz
│   ├── SRR1043449.fastq.gz
│   ├── SRR1043450.fastq.gz
│   ├── SRR1043451.fastq.gz
│   ├── SRR1043452.fastq.gz
│   ├── SRR1043453.fastq.gz
│   ├── SRR1043454.fastq.gz
│   ├── SRR1043455.fastq.gz
│   ├── SRR1043456.fastq.gz
│   ├── SRR1043457.fastq.gz
│   ├── SRR1043458.fastq.gz
│   ├── SRR1043459.fastq.gz
│   ├── SRR1043460.fastq.gz
│   ├── SRR1043461.fastq.gz
│   ├── SRR1043462.fastq.gz
│   ├── SRR1043463.fastq.gz
│   ├── SRR1043464.fastq.gz
│   ├── SRR1043465.fastq.gz
│   ├── SRR1043466.fastq.gz
│   ├── SRR1043467.fastq.gz
│   ├── SRR1043468.fastq.gz
│   ├── SRR1043469.fastq.gz
│   ├── SRR1043470.fastq.gz
│   ├── SRR1043471.fastq.gz
│   ├── SRR1043472.fastq.gz
│   ├── SRR1043473.fastq.gz
│   ├── SRR1043474.fastq.gz
│   ├── SRR1043475.fastq.gz
│   ├── SRR1043476.fastq.gz
│   ├── SRR1043477.fastq.gz
│   ├── SRR1043478.fastq.gz
│   ├── SRR1043479.fastq.gz
│   ├── SRR1043480.fastq.gz
│   ├── SRR1043481.fastq.gz
│   ├── SRR1043482.fastq.gz
│   ├── SRR1043483.fastq.gz
│   ├── SRR1043484.fastq.gz
│   ├── STRT-seq_C1_bc.csv
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

6 directories, 307 files
```