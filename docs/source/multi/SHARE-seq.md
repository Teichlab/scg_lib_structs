# SHARE-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/SHARE-seq.html) to see how __SHARE-seq__ libraries are generated experimentally. This is a split-pool based combinatorial indexing method, where open chromatin DNA are transposed by regular Tn5 and mRNA molecules are reverse transcribed by UMI-containing oligo-dT primers labelled by biotin. Then three rounds of ligation are performed to added three 8-bp barcodes. After that, ATAC and RNA molecules are separated by Streptavidin pull down, and the two libraries are prepared separately. Single cells can be identified by the combination of the three 8-bp barcodes.

## For Your Own Experiments

If you follow the protocol from the paper, you should have two libraries per sample: one for RNA and the other for ATAC. Normally, you prepare them independently and sequence them separately. If you use this assay, you have to run the sequencing for both RNA and ATAC by yourself using a custom sequencing recipe or ask your core facility to do this for you. See below the sequencing read configurations.

```{eval-rst}
.. important::
  
  Make sure you understand how sequencing is done for **SHARE-seq** by checking `this GitHub page <https://teichlab.github.io/scg_lib_structs/methods_html/SHARE-seq.html>`_. For both SHARE-ATAC and SHARE-RNA, there are a total of four reads in this order:

  =====  ================  =======  =======================================================================================
  Order  Read              Cycle    Description
  =====  ================  =======  =======================================================================================
    1    Read 1            > 30     This normally yields ``R1_001.fastq.gz``, RNA cDNA reads or ATAC genomic insert
    2    Index 1 (**i7**)  99       This normally yields ``I1_001.fastq.gz``, Cell barcodes
    3    Index 2 (**i5**)  8        This normally yields ``I2_001.fastq.gz``, Sample and modality index
    4    Read 2            > 30     This normally yields ``R2_001.fastq.gz``, RNA UMI (10bp) + polyT or ATAC genomic insert
  =====  ================  =======  =======================================================================================
```

The 2nd read (__i7__) has the following information, which will be used to identify single cells (__NOTE:__ `LB` means __Ligation Barcode__, see later section about whitelist, then you will understand):

| Length | Sequence (5' -> 3')                                                                                                 |
|--------|---------------------------------------------------------------------------------------------------------------------|
| 99 bp  | TCGGACGATCATGGG + 8 bp `LB` + CAAGTATGCAGCGCGCTCAAGCACGTGGAT + 8 bp `LB` AGTCGTACGCCGATGCGAAACATCGGCCAC + 8 bp `LB` |

After sequencing, you need to run `bcl2fastq` by yourself with a `SampleSheet.csv`. Here is an example of `SampleSheet.csv` of a NextSeq run with two different samples using the `Ad1.09` and `Ad1.17` primers from the [Supplementary Table 1](https://teichlab.github.io/scg_lib_structs/data/1-s2.0-S0092867420312538-mmc1.xlsx) from the [__SHARE-seq__ paper](https://www.sciencedirect.com/science/article/pii/S0092867420312538) as sample and modality index:

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
Sample01_ATAC,,,,,,,,Ad1.09,GATTTCCA,,
Sample01_RNA,,,,,,,,Ad1.17,TGCACGAA,,
```

If you look at the order of the sequencing read configuration above, as you can see, the first (`R1`), the 2nd (`I1`) and the 4th (`R2`) reads are all important for us. Therefore, we would like to get all of them for each sample based on sample index, that is, the 3rd read (`I2`). To do this, you should run `bcl2fastq` in the following way:

```console
bcl2fastq --use-bases-mask=Y50,Y99,I8,Y50 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is `--use-bases-mask=Y50,Y99,I8,Y50`. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y50` at the first position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `Y99` at the second position indicates "use the cycle as a real read", so you will get 99-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read.
3. `I8` at the third position indicates "use the cycle as an index read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read, though it is the 3rd read overall.
4. `Y50` at the fourth position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

Therefore, you will get four fastq files per sample per modality. Using the examples above, these are the files you should get:

```bash
# ATAC

Sample01_ATAC_S1_I1_001.fastq.gz # 8 bp: sample and modality index
Sample01_ATAC_S1_R1_001.fastq.gz # 50 bp: genomic insert
Sample01_ATAC_S1_R2_001.fastq.gz # 99 bp: barcodes and linker
Sample01_ATAC_S1_R3_001.fastq.gz # 50 bp: genomic insert

# RNA

Sample01_RNA_S2_I1_001.fastq.gz # 8 bp: sample and modality index
Sample01_RNA_S2_R1_001.fastq.gz # 50 bp: cDNA reads
Sample01_RNA_S2_R2_001.fastq.gz # 99 bp: barcodes and linker
Sample01_RNA_S2_R3_001.fastq.gz # 50 bp: 10 bp UMI + polyT
```

The naming here is really different from our normal usage. The `I1` files here actually mean `I2` in our normal usage, and we can safely ignore them. The `R1` files are the same as our normal usage. The `R2` files here actually mean `I1` in our normal usage. The `R3` files here actually means `R2` in our normal usage. Anyway, __DO NOT get confused__. You are ready to go from here.

## Public Data

We will use the data from the following paper:

```{eval-rst}
.. note::

  Ma S, Zhang B, LaFave LM, Earl AS, Chiang Z, Hu Y, Ding J, Brack A, Kartha VK, Tay T, Law T, Lareau C, Hsu Y-C, Regev A, Buenrostro JD (2020) **Chromatin Potential Identified by Shared Single-Cell Profiling of RNA and Chromatin.** *Cell* 183:1103-1116.e20. https://doi.org/10.1016/j.cell.2020.09.056

```

where the authors developed __SHARE-seq__ for the first time, and found chromatin accessibility changes precede gene expression during lineage commitment. The data is deposited to GEO under the accession code [GSE140203](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140203). However, SRA messed up the index read and the `fastq` header. Therefore, we are going to use the `fastq` files from the species mixing experiment provided by the author in the GoogleDrive. [Download from here (accessed 2022-Aug-08)](https://drive.google.com/drive/folders/19HdjJuWrpRJz8OeB6YNUMSojPyTMV7OP?usp=sharing):

```bash
mkdir -p share-seq/data

# then put the four fastq files into the directory:

scg_prep_test/share-seq/data/
├── sp.atac.R1.fastq.gz
├── sp.atac.R2.fastq.gz
├── sp.rna.R1.fastq.gz
└── sp.rna.R2.fastq.gz

0 directories, 4 files
```

Some explanation about the names of the files is needed before we proceed. The `R1` and `R2` here follow our normal naming convention described in the table at the beginning of the section, meaning they are biological reads. The index reads are not provided, but we could find the information in the `fastsq` header. For example, this is the top 20 lines (5 reads) from `sp.atac.R1.fastq.gz`:

```
@NB501804:939:HFCHMBGXG:1:11101:25826:1059_1:N:0:TCGGACGATCATGGGTTCAGCTCCAAGTATGCAGCGCGCTCAAGCACGTGGATGTAAGGTGAGTCGTACGCCGATGCGAAACATCGGCCANTNTCTTNA_NATTTCCA
GCTATGCCTCAGTTTTCTTCCCACATGAAA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE
@NB501804:939:HFCHMBGXG:1:11101:24470:1060_1:N:0:TCGGACGATCATGGGACAGTGGTCAAGTATGCAGCGCGCTCAAGCACGTGGATTCATTGAGAGTCGTACGCCGATGCGAAACATCGGCCANTNATCCNA_NACTCCTT
GGCTGGTGTCGTCTTCGGTGCGCGCCGGCG
+
AAAAAEEEEEEEEEEEEEEA<<A/EEEEAE
@NB501804:939:HFCHMBGXG:1:11101:5494:1060_1:N:0:TCGGACGATCATGGGGGTGAGTTCAAGTATGCAGCGCGCTCAAGCACGTGGATGCATGGCTAGTCGTACGCCGATGCGAAACATCGGCCACGNTGAGNT_NTTCATCA
CCCCTTTGGCGAGCTCGCGCGAGGACGTGC
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE
@NB501804:939:HFCHMBGXG:1:11101:26565:1060_1:N:0:TCGGACGATCATGGGTGGTTGTTCAAGTATGCAGCGCGCTCAAGCACGTGGATTCATTGAGAGTCGTACGCCGATGCGAACATCGGCCACTTNCATTNA_NTCATGTT
AATCATATCTCCAAACTTCCATTCAGTGCT
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE
@NB501804:939:HFCHMBGXG:1:11101:13273:1060_1:N:0:TCGGACGATCATGGGACAGTGGTCAAGTATGCAGCGCGCTCAAGCACGTGGATTCGTTAGCAGTCGTACGCCGATGCGAAACATCGGCCACTNTGCTNT_NTCATGTT
CTCAGACACACACCAGAAGAGGGCATTGGA
+
AAAAAEEEEEEEEEEEEEEEEEEEEEEEEE
```

As you can see from the header, using `:` as the delimiter, the last field is in the format of __99 bp i7__ + `_` + __8 bp i5__. We don't really need `i5`, so we could just get the __99 bp i7__ into a separate `fastq` file with fake quality strings. You can name the new file as `I1_001.fastq.gz` in the normal naming convention, or as `R2_001.fastq.gz` using the naming rules from `bcl2fastq`. To remove confusion, I'm just going to name it `CB`, which stands for __Cell Barcode__:

```bash
# output CB reads for ATAC
zcat share-seq/data/sp.atac.R1.fastq.gz | \
    awk -F ':' '(NR%4==1){print $0 "\n" substr($10,1,99) "\n+\n" "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' | \
    gzip > share-seq/data/sp.atac.CB.fastq.gz

# output CB reads for RNA
zcat share-seq/data/sp.rna.R1.fastq.gz | \
    awk -F ':' '(NR%4==1){print $0 "\n" substr($10,1,99) "\n+\n" "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' | \
    gzip > share-seq/data/sp.rna.CB.fastq.gz
```

Now, the files you have created kind of mimic the real output as if you are doing the experiment by yourself. To put those files in the context:

| Your own files after running __bcl2fastq__ | Species mixing data equivalent (__SHARE-seq__ paper) |
|--------------------------------------------|------------------------------------------------------|
| `Sample01_ATAC_S1_I1_001.fastq.gz`         | N/A, and not needed!                                 |
| `Sample01_ATAC_S1_R1_001.fastq.gz`         | __sp.atac.R1.fastq.gz__                              |
| `Sample01_ATAC_S1_R2_001.fastq.gz`         | __sp.atac.CB.fastq.gz__                              |
| `Sample01_ATAC_S1_R3_001.fastq.gz`         | __sp.atac.R2.fastq.gz__                              |
| `Sample01_RNA_S2_I1_001.fastq.gz`          | N/A, and not needed!                                 |
| `Sample01_RNA_S2_R1_001.fastq.gz`          | __sp.rna.R1.fastq.gz__                               |
| `Sample01_RNA_S2_R2_001.fastq.gz`          | __sp.rna.CB.fastq.gz__                               |
| `Sample01_RNA_S2_R3_001.fastq.gz`          | __sp.rna.R2.fastq.gz__                               |

I hope the above table gives you a rough idea about each file. __DO NOT__ get confused about the file names.

For the RNA library, the 10-bp UMI is located in the first 10 bp of `sp.rna.R1.fastq.gz`. Therefore, we need to put them together with the cell barcode in a new fastq file:

```bash
paste <(zcat share-seq/data/sp.rna.CB.fastq.gz) \
      <(zcat share-seq/data/sp.rna.R2.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print substr($1,16,8) substr($1,54,8) substr($1,92,8) substr($2,1,10)} }' | \
      gzip > share-seq/data/sp.rna.CB_UMI.fastq.gz
```

After that, we are ready to begin the preprocessing.

## Prepare Whitelist

There are three rounds of ligation. Each round will add 8-bp __Ligation Barcode__ to the molecules. There are 96 different __Ligation Barcodes__ in each round. The same set of 96 __Ligation Barcodes__ are used in each round. Single cells can be identified by the combination of themselves. Here is the information from the [Supplementary Table 1](https://teichlab.github.io/scg_lib_structs/data/1-s2.0-S0092867420312538-mmc1.xlsx) from the [__SHARE-seq__ paper](https://www.sciencedirect.com/science/article/pii/S0092867420312538):

| WellPosition | Name          | Sequence | Reverse complement |
|--------------|---------------|----------|:------------------:|
| A1           | Round1/2/3_01 | AACGTGAT |      ATCACGTT      |
| B1           | Round1/2/3_02 | AAACATCG |      CGATGTTT      |
| C1           | Round1/2/3_03 | ATGCCTAA |      TTAGGCAT      |
| D1           | Round1/2/3_04 | AGTGGTCA |      TGACCACT      |
| E1           | Round1/2/3_05 | ACCACTGT |      ACAGTGGT      |
| F1           | Round1/2/3_06 | ACATTGGC |      GCCAATGT      |
| G1           | Round1/2/3_07 | CAGATCTG |      CAGATCTG      |
| H1           | Round1/2/3_08 | CATCAAGT |      ACTTGATG      |
| A2           | Round1/2/3_09 | CGCTGATC |      GATCAGCG      |
| B2           | Round1/2/3_10 | ACAAGCTA |      TAGCTTGT      |
| C2           | Round1/2/3_11 | CTGTAGCC |      GGCTACAG      |
| D2           | Round1/2/3_12 | AGTACAAG |      CTTGTACT      |
| E2           | Round1/2/3_13 | AACAACCA |      TGGTTGTT      |
| F2           | Round1/2/3_14 | AACCGAGA |      TCTCGGTT      |
| G2           | Round1/2/3_15 | AACGCTTA |      TAAGCGTT      |
| H2           | Round1/2/3_16 | AAGACGGA |      TCCGTCTT      |
| A3           | Round1/2/3_17 | AAGGTACA |      TGTACCTT      |
| B3           | Round1/2/3_18 | ACACAGAA |      TTCTGTGT      |
| C3           | Round1/2/3_19 | ACAGCAGA |      TCTGCTGT      |
| D3           | Round1/2/3_20 | ACCTCCAA |      TTGGAGGT      |
| E3           | Round1/2/3_21 | ACGCTCGA |      TCGAGCGT      |
| F3           | Round1/2/3_22 | ACGTATCA |      TGATACGT      |
| G3           | Round1/2/3_23 | ACTATGCA |      TGCATAGT      |
| H3           | Round1/2/3_24 | AGAGTCAA |      TTGACTCT      |
| A4           | Round1/2/3_25 | AGATCGCA |      TGCGATCT      |
| B4           | Round1/2/3_26 | AGCAGGAA |      TTCCTGCT      |
| C4           | Round1/2/3_27 | AGTCACTA |      TAGTGACT      |
| D4           | Round1/2/3_28 | ATCCTGTA |      TACAGGAT      |
| E4           | Round1/2/3_29 | ATTGAGGA |      TCCTCAAT      |
| F4           | Round1/2/3_30 | CAACCACA |      TGTGGTTG      |
| G4           | Round1/2/3_31 | GACTAGTA |      TACTAGTC      |
| H4           | Round1/2/3_32 | CAATGGAA |      TTCCATTG      |
| A5           | Round1/2/3_33 | CACTTCGA |      TCGAAGTG      |
| B5           | Round1/2/3_34 | CAGCGTTA |      TAACGCTG      |
| C5           | Round1/2/3_35 | CATACCAA |      TTGGTATG      |
| D5           | Round1/2/3_36 | CCAGTTCA |      TGAACTGG      |
| E5           | Round1/2/3_37 | CCGAAGTA |      TACTTCGG      |
| F5           | Round1/2/3_38 | CCGTGAGA |      TCTCACGG      |
| G5           | Round1/2/3_39 | CCTCCTGA |      TCAGGAGG      |
| H5           | Round1/2/3_40 | CGAACTTA |      TAAGTTCG      |
| A6           | Round1/2/3_41 | CGACTGGA |      TCCAGTCG      |
| B6           | Round1/2/3_42 | CGCATACA |      TGTATGCG      |
| C6           | Round1/2/3_43 | CTCAATGA |      TCATTGAG      |
| D6           | Round1/2/3_44 | CTGAGCCA |      TGGCTCAG      |
| E6           | Round1/2/3_45 | CTGGCATA |      TATGCCAG      |
| F6           | Round1/2/3_46 | GAATCTGA |      TCAGATTC      |
| G6           | Round1/2/3_47 | CAAGACTA |      TAGTCTTG      |
| H6           | Round1/2/3_48 | GAGCTGAA |      TTCAGCTC      |
| A7           | Round1/2/3_49 | GATAGACA |      TGTCTATC      |
| B7           | Round1/2/3_50 | GCCACATA |      TATGTGGC      |
| C7           | Round1/2/3_51 | GCGAGTAA |      TTACTCGC      |
| D7           | Round1/2/3_52 | GCTAACGA |      TCGTTAGC      |
| E7           | Round1/2/3_53 | GCTCGGTA |      TACCGAGC      |
| F7           | Round1/2/3_54 | GGAGAACA |      TGTTCTCC      |
| G7           | Round1/2/3_55 | GGTGCGAA |      TTCGCACC      |
| H7           | Round1/2/3_56 | GTACGCAA |      TTGCGTAC      |
| A8           | Round1/2/3_57 | GTCGTAGA |      TCTACGAC      |
| B8           | Round1/2/3_58 | GTCTGTCA |      TGACAGAC      |
| C8           | Round1/2/3_59 | GTGTTCTA |      TAGAACAC      |
| D8           | Round1/2/3_60 | TAGGATGA |      TCATCCTA      |
| E8           | Round1/2/3_61 | TATCAGCA |      TGCTGATA      |
| F8           | Round1/2/3_62 | TCCGTCTA |      TAGACGGA      |
| G8           | Round1/2/3_63 | TCTTCACA |      TGTGAAGA      |
| H8           | Round1/2/3_64 | TGAAGAGA |      TCTCTTCA      |
| A9           | Round1/2/3_65 | TGGAACAA |      TTGTTCCA      |
| B9           | Round1/2/3_66 | TGGCTTCA |      TGAAGCCA      |
| C9           | Round1/2/3_67 | TGGTGGTA |      TACCACCA      |
| D9           | Round1/2/3_68 | TTCACGCA |      TGCGTGAA      |
| E9           | Round1/2/3_69 | AACTCACC |      GGTGAGTT      |
| F9           | Round1/2/3_70 | AAGAGATC |      GATCTCTT      |
| G9           | Round1/2/3_71 | AAGGACAC |      GTGTCCTT      |
| H9           | Round1/2/3_72 | AATCCGTC |      GACGGATT      |
| A10          | Round1/2/3_73 | AATGTTGC |      GCAACATT      |
| B10          | Round1/2/3_74 | ACACGACC |      GGTCGTGT      |
| C10          | Round1/2/3_75 | ACAGATTC |      GAATCTGT      |
| D10          | Round1/2/3_76 | AGATGTAC |      GTACATCT      |
| E10          | Round1/2/3_77 | AGCACCTC |      GAGGTGCT      |
| F10          | Round1/2/3_78 | AGCCATGC |      GCATGGCT      |
| G10          | Round1/2/3_79 | AGGCTAAC |      GTTAGCCT      |
| H10          | Round1/2/3_80 | ATAGCGAC |      GTCGCTAT      |
| A11          | Round1/2/3_81 | ATCATTCC |      GGAATGAT      |
| B11          | Round1/2/3_82 | ATTGGCTC |      GAGCCAAT      |
| C11          | Round1/2/3_83 | CAAGGAGC |      GCTCCTTG      |
| D11          | Round1/2/3_84 | CACCTTAC |      GTAAGGTG      |
| E11          | Round1/2/3_85 | CCATCCTC |      GAGGATGG      |
| F11          | Round1/2/3_86 | CCGACAAC |      GTTGTCGG      |
| G11          | Round1/2/3_87 | CCTAATCC |      GGATTAGG      |
| H11          | Round1/2/3_88 | CCTCTATC |      GATAGAGG      |
| A12          | Round1/2/3_89 | CGACACAC |      GTGTGTCG      |
| B12          | Round1/2/3_90 | CGGATTGC |      GCAATCCG      |
| C12          | Round1/2/3_91 | CTAAGGTC |      GACCTTAG      |
| D12          | Round1/2/3_92 | GAACAGGC |      GCCTGTTC      |
| E12          | Round1/2/3_93 | GACAGTGC |      GCACTGTC      |
| F12          | Round1/2/3_94 | GAGTTAGC |      GCTAACTC      |
| G12          | Round1/2/3_95 | GATGAATC |      GATTCATC      |
| H12          | Round1/2/3_96 | GCCAAGAC |      GTCTTGGC      |

Since during each ligation round, the same set of __Ligation Barcodes__ (96) are used. Therefore, the whitelist is basically the combination of those 96 barcodes themselves for three times: a total of __96 * 96 * 96 = 884736__ barcodes. Since the barcodes are sequenced as the `i7` index, which uses the bottom strand as the template, we should use the reverse complement to construct the whitelist. Again, if you are confused, check the [SHARE-seq GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/SHARE-seq.html). I have put the above table into a `csv` file so that you can download by [__click here__](https://teichlab.github.io/scg_lib_structs/data/share-seq_ligationBC.csv).

```bash
# download the ligation barcode file
wget -P share-seq/data https://teichlab.github.io/scg_lib_structs/data/share-seq_ligationBC.csv

# generate whitelist
for x in $(tail -n +2 share-seq/data/share-seq_ligationBC.csv | cut -f 4 -d,); do
    for y in $(tail -n +2 share-seq/data/share-seq_ligationBC.csv | cut -f 4 -d,); do
        for z in $(tail -n +2 share-seq/data/share-seq_ligationBC.csv | cut -f 4 -d,); do
            echo "${x}${y}${z}"
            done
        done
    done > share-seq/data/whitelist.txt
```

## From FastQ To Count Matrices

Let's map the reads to the genome using `starsolo` for the RNA library and `chromap` for the ATAC library:

```console
mkdir -p share-seq/star_outs
mkdir -p share-seq/chromap_outs

# process the RNA library using starsolo

STAR --runThreadN 4 \
     --genomeDir mix_hg38_mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix share-seq/star_outs/ \
     --readFilesIn share-seq/data/sp.rna.R1.fastq.gz share-seq/data/sp.rna.CB_UMI.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 24 --soloUMIstart 25 --soloUMIlen 10 \
     --soloCBwhitelist share-seq/data/whitelist.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate

# process the ATAC library using chromap

## map and generate the fragment file

chromap -t 4 --preset atac \
        -x mix_hg38_mm10/chromap_index/genome.index \
        -r mix_hg38_mm10/genome.fa \
        -1 share-seq/data/sp.atac.R1.fastq.gz \
        -2 share-seq/data/sp.atac.R2.fastq.gz \
        -b share-seq/data/sp.atac.CB.fastq.gz \
        --read-format bc:15:22,bc:53:60,bc:91:98 \
        --barcode-whitelist share-seq/data/whitelist.txt \
        -o share-seq/chromap_outs/fragments.tsv

## compress and index the fragment file

bgzip share-seq/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed share-seq/chromap_outs/fragments.tsv.gz
```

After this stage, we are done with the RNA library. The count matrix and other useful information can be found in the `star_outs` directory. For the ATAC library, two new files `fragments.tsv.gz` and `fragments.tsv.gz.tbi` are generated. They will be useful and sometimes required for other programs to perform downstream analysis. There are still some extra work.

### Explain star and chromap

If you understand the __SHARE-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/SHARE-seq.html), the commands above should be straightforward to understand.

#### Explain star

`--runThreadN 4`
  
>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mix_hg38_mm10/chromap_index/genome.index`

>> Pointing to the directory of the star index. The public data we are analysing is from human + mouse species mixing experiments.

`--readFilesCommand zcat`

>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix share-seq/star_outs/`

>> We want to keep everything organised. This directs all output files inside the `share-seq/star_outs` directory.

`--readFilesIn share-seq/data/sp.rna.R1.fastq.gz share-seq/data/sp.rna.CB_UMI.fastq.gz`

>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. In __SHARE-seq__, cDNA reads come from Read 1, and the cell barcode and UMI come from the file we just prepared previously. Check [the SHARE-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/SHARE-seq.html) if you are not sure. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--soloType CB_UMI_Simple`

>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 1 --soloCBlen 24 --soloUMIstart 25 --soloUMIlen 10`

>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist share-seq/data/whitelist.txt`

>> The plain text file containing all possible valid cell barcodes, one per line. __SHARE-seq__ uses the combination of the three 8-bp barcodes. Put the file we just prepared in the previous section.

`--soloCellFilter EmptyDrops_CR`

>> Experiments are never perfect. Even for barcodes that do not capture the molecules inside the cells, you may still get some reads due to various reasons, such as ambient RNA or DNA and leakage. In general, the number of reads from those cell barcodes should be much smaller, often orders of magnitude smaller, than those barcodes that come from real cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Forward`

>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __SHARE-seq__, the cDNA reads are from the Read 1 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing UMI and the Read 2 sequence. Therefore, Read 2 consists of cell barcodes and UMI come from the first strand, complementary to the coding strand. Read 1 comes from the coding strand. Therefore, use `Forward` for __SHARE-seq__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the SHARE-seq GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/SHARE-seq.html) if you are not sure.

`--outSAMattributes CB UB`

>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>> We want sorted `BAM` for easy handling by other programs.

#### Explain chromap

`-t 4`

>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`-x mix_hg38_mm10/chromap_index/genome.index`

>> The `chromap` index file. The public data we are analysing is from human mouse species mixing experiments.

`-r mix_hg38_mm10/genome.fa`

>> Reference genome sequence in `fasta` format. This is basically the file which you used to create the `chromap` index file.

`-1`, `-2` and `-b`

>> They are Read 1 (genomic), Read 2 (genomic) and cell barcode read, respectively. For ATAC-seq, the sequencing is usually done in pair-end mode. Therefore, you normally have two genomic reads for each genomic fragment: Read 1 and Read 2. For the reason described previously, `R1` is the genomic Read 1 and should be passed to `-1`; `R2` is the genomic Read 2 and should be passed to `-2`; `CB` is the cell barcode read we just prepared in the previous section and should be passed to `-b`. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--read-format bc:15:22,bc:53:60,bc:91:98`

>> Note that `sp.atac.CB.fastq.gz` contains the cell barcodes, which is the combination of three 8-bp barcodes separated by defined linker regions. Barcodes can be non-consecutive segments, which can be specified multiple times separated by `,`. The reads are 99-bp long, the first 15 bp are __TCGGACGATCATGGG__, position 16 - 23 are the first 8-bp barcode, position 24 - 53 are __CAAGTATGCAGCGCGCTCAAGCACGTGGAT__, position 54 - 61 are the second 8-bp barcode, position 62 - 91 are __AGTCGTACGCCGATGCGAAACATCGGCCAC__ and position 92 - 99 are the third 8-bp barcode. Therefore, we tell `chromap` to only use those positions (`15:22`, `53:60` and `91:98`) of the barcode file (`bc`) as the cell barcode. Be aware that the position is 0-based (the first base of the read is 0, __NOT__ 1). Check the `chromap` manual if you are not sure.

`--barcode-whitelist share-seq/data/whitelist.txt`

>> The plain text file containing all possible valid cell barcodes, one per line. This is the whitelist we just prepared in the previous section. It contains all possible combination of the 96 8-bp barcodes for three times. A total of 96 * 96 * 96 = 884736 barcodes are in this file.

`-o share-seq/chromap_outs/fragments.tsv`

>> Direct the mapped fragments to a file. The format is described in the [10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments).

```{eval-rst}
.. important::

  After this stage, we normally should look at the percentage of the reads mapped to each species. Make the decision if the cell is human or mouse or mix. After that, separate the human and mouse cells, and perform analysis on cells from the same species. However, for the demonstation, I'm doing all of cells together.

```

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

faSize -detailed mix_hg38_mm10/genome.fa | \
    sort -k1,1 > mix_hg38_mm10/genome.chrom.sizes
```

This is the first 5 lines of `mix_hg38_mm10/genome.chrom.sizes`:

```
hg38_chr1	248956422
hg38_chr10	133797422
hg38_chr11	135086622
hg38_chr11_KI270721v1_random	100316
hg38_chr12	133275309
```

Now let's generate the reads from fragments:

```bash
# we use bedClip to remove reads outside the chromosome boundary
# we also remove reads mapped to the mitochondrial genome (chrM)

zcat share-seq/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin mix_hg38_mm10/genome.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > share-seq/chromap_outs/reads.bed.gz
```

Note we also sort the output reads by `sort -k1,1 -k2,2n`. In this way, the order of chromosomes in the `reads.bed.gz` is the same as that in `genome.chrom.sizes`, which makes downstream processes easier. The output `reads.bed.gz` are the reads in `bed` format, with the 4th column holding the cell barcodes.

### Peak Calling By MACS2

Now we can use the newly generated read file for the peak calling using `MACS2`:

```console
macs2 callpeak -t share-seq/chromap_outs/reads.bed.gz \
               -g 4.57e9 -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir share-seq/chromap_outs \
               -n aggregate
```

#### Explain MACS2

The reasons of choosing those specific parameters are a bit more complicated. I have dedicated a post for this a while ago. Please have a look at [__this post__](https://dbrg77.github.io/posts/2020-12-09-atac-seq-peak-calling-with-macs2/) if you are still confused. Note the `-g`, which is the genome size parameter, is basically the sum of human and mouse. The following output files are particularly useful:

| File                       | Description                                                               |
|----------------------------|---------------------------------------------------------------------------|
| aggregate_peaks.narrowPeak | Open chromatin peak locations in the narrowPeak format                    |
| aggregate_peaks.xls        | More information about peaks                                              |
| aggregate_treat_pileup.bdg | Signal tracks. Can be used to generate the bigWig file for visualisation |

### Getting The Peak-By-Cell Count Matrix

Now that we have the peak and reads files, we can compute the number of reads in each peak for each cell. Then we could get the peak-by-cell count matrix. There are different ways of doing this. The following is the method I use.

#### Find Reads In Peaks Per Cell

First, we use the `aggregate_peaks.narrowPeak` file. We only need the first 4 columns (chromosome, start, end, peak ID). You can also remove the peaks that overlap [the black list regions](https://www.nature.com/articles/s41598-019-45839-z). The black list is not available for every species and every build, so I'm not doing it here. We also need to sort the peak to make sure the order of the chromosomes in the peak file is the same as that in the `genome.chrom.sizes` and `reads.bed.gz` files. Then we could find the overlap by `bedtools`. We need to do this in a specific way to get the number of reads in each peak from each cell:

```bash
# format and sort peaks

cut -f 1-4 share-seq/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > share-seq/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a share-seq/chromap_outs/aggregate_peaks_sorted.bed \
    -b share-seq/chromap_outs/reads.bed.gz \
    -wo -sorted -g mix_hg38_mm10/genome.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > share-seq/chromap_outs/peak_read_ov.tsv.gz
```

##### Explain Finding Reads In Peaks Per Cell

We start with the command before the first pipe, that is, the intersection part. If you read the manual of the `bedtools intersect`, it should be straightforward to understand. The `-wo` option will output the records in both `-a` file and `-b` file. Since the `reads.bed.gz` file has the cell barcode information at the 4th column, we would get an output with both peak and cell information for the overlap. The `-sorted -g mix_hg38_mm10/genome.chrom.sizes` options make the program use very little memory. Here is an example (top 5 lines) of the output of this part:

```console
hg38_chr1	181021	181929	aggregate_peak_1	hg38_chr1	180979	181029	TTCCATTGTTGGAGGTGATAGAGG	.	+	8
hg38_chr1	181021	181929	aggregate_peak_1	hg38_chr1	180982	181032	ACAGTGGTTAAGCGTTGTTAGCCT	.	-	11
hg38_chr1	181021	181929	aggregate_peak_1	hg38_chr1	180991	181041	TACTTCGGGGTCGTGTTCATTGAG	.	+	20
hg38_chr1	181021	181929	aggregate_peak_1	hg38_chr1	181074	181124	TCGAGCGTTTACTCGCTAGCTTGT	.	-	50
hg38_chr1	181021	181929	aggregate_peak_1	hg38_chr1	181079	181129	GACGGATTTTACTCGCTAGCTTGT	.	-	50
```

We see that the 8th column holds the cell barcode and we want to group them using `bedtools groupby`. Therefore, we need to sort by this column, that is the `sort -k8,8`. When we group by the 8th column, we are interested in how many times each peak appear per group, so we could gather the information of the peak ID (4th column). That is the `-g 8 -c 4 -o freqdesc`. The `-o freqdesc` option returns a `value:frequency` pair in descending order. Here are some records from `peak_read_ov.tsv.gz`:

```console
ACAGTGGTACAGTGGTACTTGATG	aggregate_peak_1909:2,aggregate_peak_29080:2,aggregate_peak_34079:2,aggregate_peak_145661:1,aggregate_peak_35738:1
ACAGTGGTACAGTGGTATCACGTT	aggregate_peak_159342:1
ACAGTGGTACAGTGGTCAGATCTG	aggregate_peak_105556:2,aggregate_peak_118790:2,aggregate_peak_126632:2,aggregate_peak_135617:2,aggregate_peak_140048:2,aggregate_peak_20789:2,aggregate_peak_3278:2,aggregate_peak_68156:2,aggregate_peak_106897:1,aggregate_peak_152567:1
```

In a way, that is sort of a count matrix in an awkward format. For example:

- The first line means that in cell `ACAGTGGTACAGTGGTACTTGATG`, the peak `aggregate_peak_1909` has 2 counts, the peak `aggregate_peak_29080` has 2 counts, the peak `aggregate_peak_34079` has 2 counts and the peak `aggregate_peak_145661` has 1 count. All the rest peaks not mentioned here have 0 counts in this cell.
- The second line means that in cell `ACAGTGGTACAGTGGTATCACGTT`, the peak `aggregate_peak_159342` has 1 count. All the rest peaks not mentioned here have 0 counts in this cell.

#### Output The Peak-By-Cell Matrix

At this stage, we pretty much have all the things needed. Those two files `aggregate_peaks_sorted.bed` and `peak_read_ov.tsv.gz` contain all information for a peak-by-cell count matrix. We just need a final touch to make the output in a standard format: a [market exchange format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). Since most downstream software takes the input from the __10x Genomics Single Cell ATAC__ results, we are going to generate the MEX and the associated files similar to the output from 10x Genomics.

Here, I'm using a python script for this purpose. You don't have to do this. Choose whatever works for you. The point here is to just generate similar files as the __peak-barcode matrix__ described from [the 10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

First, let's make a directory to hold the output files and generate the `peaks.bed` and `barcodes.tsv` files, which are easy to do:

```bash
# create dirctory
mkdir -p share-seq/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 share-seq/chromap_outs/aggregate_peaks_sorted.bed > \
    share-seq/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat share-seq/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    share-seq/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
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
    share-seq/chromap_outs/aggregate_peaks_sorted.bed \
    share-seq/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    share-seq/chromap_outs/peak_read_ov.tsv.gz \
    share-seq/chromap_outs/raw_peak_bc_matrix
```

After that, you should have the `matrix.mtx` in the `share-seq/chromap_outs/raw_peak_bc_matrix` directory.

#### Cell Calling (Filter Cell Barcodes)

Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, we could use `starolo`. It is used for scRNA-seq in general, but it does have a cell calling function that takes a directory containing raw mtx and associated files, and return the filtered ones. Since `starsolo` looks for the following three files in the input directory: `matrix.mtx`, `features.tsv` and `barcodes.tsv`. Those are the output from the 10x Genomics scRNA-seq workflow. In this case, we can use `peaks.bed` as our `features.tsv`:

```console
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed share-seq/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     share-seq/chromap_outs/raw_peak_bc_matrix \
     share-seq/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv share-seq/chromap_outs/filtered_peak_bc_matrix/peaks.bed
```

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/share-seq/
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
│   ├── share-seq_ligationBC.csv
│   ├── sp.atac.CB.fastq.gz
│   ├── sp.atac.R1.fastq.gz
│   ├── sp.atac.R2.fastq.gz
│   ├── sp.rna.CB.fastq.gz
│   ├── sp.rna.CB_UMI.fastq.gz
│   ├── sp.rna.R1.fastq.gz
│   ├── sp.rna.R2.fastq.gz
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

9 directories, 42 files
```