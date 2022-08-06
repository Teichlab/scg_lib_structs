# snATAC-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/snATAC-seq.html) to see how __snATAC-seq__ libraries are generated experimentally. The library preparation steps are the same as the __sci-ATAC-seq__ method. It is a split-pool based combinatorial indexing strategy, where nuclei are transposed in mini-bulk with indexed transposase Tn5 (__T5__ + __T7__). Then all nuclei are pooled and a fraction of the nuclei is randomly distributed into new wells. Library preparation is performed using a PCR primers with different indices (__i5__ + __i7__). Single cells can be identified by the combination of __T7 + i7 + T5 + i5__. The developers named the Tn5 barcode at the __i5__ side `T5`, and the __i7__ side `T7`. In this documentation, I will keep the name consistent with the developers.

```{eval-rst}
.. important::

   The library of the **snATAC-seq** is the same as **sci-ATAC-seq**, with minor sequence difference. It is kind of complicated, especially for beginners. First, make sure you are familiar with different sequencing modes from different Illumina machines by looking at `this page <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_. Then make sure you understand how the actual sequencing is done for **snATAC-seq** by checking `this GitHub page <https://teichlab.github.io/scg_lib_structs/methods_html/snATAC-seq.html>`_ and the associated publications.

```

## For Your Own Experiments

If you use this assay, you have to run the sequencing by yourself using a custom sequencing recipe or ask your core facility to do this for you. We don't have the proper setup to do this experiment here, so I haven't tried __snATAC-seq__ by myself. This is the educational guess of the sequencing configuration. Spike-in libraries are used in this assay, so that dark cycles are not needed when sequencing the index.

Using more recent machines and chemistries, like __iSeq 100__, __MiniSeq__, __NextSeq__, __HiSeq X__, __HiSeq 3000__, __HiSeq 4000__ and __NovaSeq 600 (v1.5)__, it should be (I call this __Configuration 1__):

| Order | Read             | Cycle | Description                                               |
|-------|------------------|-------|-----------------------------------------------------------|
| 1     | Read 1           | >50   | `R1_001.fastq.gz`, Genomic insert                         |
| 2     | Index 1 (__i7__) | 43    | `I1_001.fastq.gz`, 8-bp T7 + 27-bp linker + 8-bp i7 index |
| 3     | Index 2 (__i5__) | 37    | `I2_001.fastq.gz`, 8-bp T5 + 21-bp linker + 8-bp i5 index |
| 4     | Read 2           | >50   | `R2_001.fastq.gz`, Genomic insert                         |

Using older machines and chemistries like __MiSeq__, __HiSeq 2000__, __HiSeq 2500__, __NovaSeq 6000 (v1.0)__, it should be (I call this __Configuration 2__):

| Order | Read             | Cycle | Description                                               |
|-------|------------------|-------|-----------------------------------------------------------|
| 1     | Read 1           | >50   | `R1_001.fastq.gz`, Genomic insert                         |
| 2     | Index 1 (__i7__) | 43    | `I1_001.fastq.gz`, 8-bp T7 + 27-bp linker + 8-bp i7 index |
| 3     | Index 2 (__i5__) | 37    | `I2_001.fastq.gz`, 8-bp i5 + 21-bp linker + 8-bp T5 index |
| 4     | Read 2           | >50   | `R2_001.fastq.gz`, Genomic insert                         |

The two configurations are very similar. The only difference is how `I2` is sequenced and generated. Again, if you are not sure, check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/snATAC-seq.html) to see the reason.

After the sequencing is done, you could just run the `bcl2fastq` without a `SampleSheet.csv`, like this:

```console
bcl2fastq --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

In this case, you will have four `fastq` files per run:

```bash
Undetermined_S0_I1_001.fastq.gz
Undetermined_S0_I2_001.fastq.gz
Undetermined_S0_R1_001.fastq.gz # >50 bp, genomic insert
Undetermined_S0_R2_001.fastq.gz # >50 bp, genomic insert
```

The `R1` and `R2` files contain genomic inserts, which can be used directly. We don't need to do anything about them.

This is the content of `Undetermined_S0_I1_001.fastq.gz` regardless of the machine configuration:

| Length | Sequence (5' -> 3')                                 |
|--------|-----------------------------------------------------|
| 43 bp  | 8 bp `T7` + GGACAGGGACAGCCGAGCCCACGAGAC + 8 bp `i7` |

This is the content of `Undetermined_S0_I1_001.fastq.gz`:

| Configuration   | Length | Sequence (5' -> 3')                           |
|-----------------|--------|-----------------------------------------------|
| Configuration 1 | 37 bp  | 8 bp `T5` + GCGTGGAGACGCTGCCGACGA + 8 bp `i5` |
| Configuration 2 | 37 bp  | 8 bp `i5` + TCGTCGGCAGCGTCTCCACGC + 8 bp `T5` |

As you can see, your cell barcodes will be the first 8bp and the last 8 bp of `I1` + the first 8bp and the last 8 bp of `I2`. We could stitch them together to get the 32-bp cell barcodes to a single `fastq` file:

```bash
paste <(zcat Undetermined_S0_I1_001.fastq.gz) <(zcat Undetermined_S0_I2_001.fastq.gz) | \
    awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print substr($1,1,8) substr($1,36,8) substr($2,1,8) substr($2,30,8)} }' | \
    gzip > Undetermined_S0_CB_001.fastq.gz
```

After that, you are ready to throw `Undetermined_S0_R1_001.fastq.gz`, `Undetermined_S0_R2_001.fastq.gz` and `Undetermined_S0_CB_001.fastq.gz` into `chromap`. See the later section.

## Public Data

For the purpose of demonstration, we will use the __snATAC-seq__ data from the following paper:

```{eval-rst}
.. note::
  
  Fang R, Preissl S, Li Y, Hou X, Lucero J, Wang X, Motamedi A, Shiau AK, Zhou X, Xie F, Mukamel EA, Zhang K, Zhang Y, Behrens MM, Ecker JR, Ren B (2021) **Comprehensive analysis of single cell ATAC-seq data with SnapATAC.** *Nat Commun* 12:1337. https://doi.org/10.1038/s41467-021-21583-9

```

where the authors developed a computational tool called [SnapATAC](https://github.com/r3fang/SnapATAC/) for the analysis of scATAC-seq data. They used __snATAC-seq__ to generate chromatin accessibility profiles of tens of thousands of nuclei from mouse secondary motor cortex. The accession code is [__GSE126724__](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126724). There are quite a few samples there. I'm using the __snATAC MOs A rep1__ sample (__SRR8592624__). To get the reads, we use `fastq-dump` from the `sratools`:

```console
mkdir -p snATAC-seq/data
fastq-dump --split-files \
           --origfmt \
           --outdir snATAC-seq/data \
           --defline-seq '@rd.$si:$sg:$sn' SRR8592624
```

The reason of using the specific options above is a bit complicated. I wrote [__a post__](https://dbrg77.github.io/posts/2022-07-26-getting-index-reads-from-sra/) about this, and you can have a look to see the reason. Since the authors also uploaded the index reads, once the program finishes running, you will have four plain files: `SRR8592624_1.fastq`, `SRR8592624_2.fastq`, `SRR8592624_3.fastq` and `SRR8592624_4.fastq`. Unfortunately, it seems `fastq-dump` cannot generate gzipped file on the fly. I'm going to gzip them manually to save space:

```console
gzip snATAC-seq/data/SRR8592624_*.fastq
```

`SRR8592624_1.fastq` and `SRR8592624_2.fastq` are the genomic inserts. `SRR8592624_3.fastq` and `SRR8592624_4.fastq` are the index reads. Let's have a look at the first 5 reads from `SRR8592624_3.fastq` and `SRR8592624_4.fastq`:

```bash
# from SRR8592624_3.fastq
@rd.1::AGACGGAGTGCAGCTATGTCAGCTACTGCATA:7001113:948:HLFKYBCX2:1:1108:1237:1880
AGACGGAGTGCAGCTA
+AGACGGAGTGCAGCTATGTCAGCTACTGCATA:7001113:948:HLFKYBCX2:1:1108:1237:1880
DBDDDIHHIIEEFHGH
@rd.2::AGAGCAGTAACCAGGTACAAGGATAAGAGATG:7001113:948:HLFKYBCX2:1:1108:1027:1904
AGAGCAGTAACCAGGT
+AGAGCAGTAACCAGGTACAAGGATAAGAGATG:7001113:948:HLFKYBCX2:1:1108:1027:1904
ADDDDIIIIGGGHIIF
@rd.3::ACATTGGCGCTCATGAGGAAGACTTTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1056:1909
ACATTGGCGCTCATGA
+ACATTGGCGCTCATGAGGAAGACTTTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1056:1909
DDDDDIIIIIIIIIII
@rd.4::TGTTACCATCGACGTCCATTCAGTCCTAGAGT:7001113:948:HLFKYBCX2:1:1108:1119:1911
TGTTACCATCGACGTC
+TGTTACCATCGACGTCCATTCAGTCCTAGAGT:7001113:948:HLFKYBCX2:1:1108:1119:1911
DDD@DHHIIIIHIIHI
@rd.5::CTATTAGGCCTAAGACTGCTGGTATTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1028:1934
CTATTAGGCCTAAGAC
+CTATTAGGCCTAAGACTGCTGGTATTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1028:1934
CDCDDHIIIHIIIIGH

# from SRR8592624_4.fastq
@rd.1::AGACGGAGTGCAGCTATGTCAGCTACTGCATA:7001113:948:HLFKYBCX2:1:1108:1237:1880
TGTCAGCTACTGCATA
+AGACGGAGTGCAGCTATGTCAGCTACTGCATA:7001113:948:HLFKYBCX2:1:1108:1237:1880
IHFHIIFHGHFEEHHH
@rd.2::AGAGCAGTAACCAGGTACAAGGATAAGAGATG:7001113:948:HLFKYBCX2:1:1108:1027:1904
ACAAGGATAAGAGATG
+AGAGCAGTAACCAGGTACAAGGATAAGAGATG:7001113:948:HLFKYBCX2:1:1108:1027:1904
GHHHIIIEHHIIGHII
@rd.3::ACATTGGCGCTCATGAGGAAGACTTTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1056:1909
GGAAGACTTTCTAGCT
+ACATTGGCGCTCATGAGGAAGACTTTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1056:1909
IIIIIIIIIIIIIIII
@rd.4::TGTTACCATCGACGTCCATTCAGTCCTAGAGT:7001113:948:HLFKYBCX2:1:1108:1119:1911
CATTCAGTCCTAGAGT
+TGTTACCATCGACGTCCATTCAGTCCTAGAGT:7001113:948:HLFKYBCX2:1:1108:1119:1911
IIICHIIIHIHIHHIH
@rd.5::CTATTAGGCCTAAGACTGCTGGTATTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1028:1934
TGCTGGTATTCTAGCT
+CTATTAGGCCTAAGACTGCTGGTATTCTAGCT:7001113:948:HLFKYBCX2:1:1108:1028:1934
IIIIIIIIIEHHIIIH
```

Basically, `SRR8592624_3.fastq` is `I1`, and `SRR8592624_4.fastq` is `I2`. It seems the authors already trimmed off the linker sequences (27 bp in `I1` and 21 bp in `I2`) for us, so we could stich them together simply by:

```bash
paste <(zcat snATAC-seq/data/SRR8592624_3.fastq.gz) \
      <(zcat snATAC-seq/data/SRR8592624_4.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > snATAC-seq/data/SRR8592624_CB.fastq.gz
```

Then, we are ready to pass `SRR8592624_1.fastq.gz`,`1SRR8592624_2.fastq.gz` and `SRR8592624_CB.fastq.gz` to `chromap`.

## Prepare Whitelist

This is the most confusing part of the pipeline, but if you understand the procedures in the __snATAC-seq__ [GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/snATAC-seq.html), you should be able to understand.

```{eval-rst}
.. important::

  The data we are using in this documentation is from the `SnapATAC <https://www.nature.com/articles/s41467-021-21583-9>`_ paper. The protocol is the same as the `original snATAC-seq paper <https://www.nature.com/articles/s41593-018-0079-3>`_, but there are some minor modifications of the oligos. For the Tn5 barcode sequences, you need to check the `Supplementary Table S4 <https://teichlab.github.io/scg_lib_structs/data/41467_2021_21583_MOESM1_ESM.pdf>`_ from the **SnapATAC** paper. For the PCR primer (i5 and i7), you need to check the `Supplementary Table S5 <https://teichlab.github.io/scg_lib_structs/data/41593_2018_79_MOESM5_ESM.xlsx>`_ from the **original snATAC-seq paper**.

```

This is the Tn5-T7 oligo sequence: 5'- GTCTCGTGGGCTCGGCTGTCCCTGTCCNNNNNNNNCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG -3'. The 8 Ns are the barcode, and the barcode sequences (from 5' -> 3') are:

| Oligo | Barcode T7 |  Reverse complement  |
|-------|------------|:--------------------:|
| T7_12 |  CCTAATAG  |       CTATTAGG       |
| T7_13 |  CTGACATG  |       CATGTCAG       |
| T7_14 |  TGTATGAG  |       CTCATACA       |
| T7_15 |  GAAGATCT  |       AGATCTTC       |
| T7_16 |  ACTGCTCT  |       AGAGCAGT       |
| T7_17 |  CTCCGTCT  |       AGACGGAG       |
| T7_18 |  GTTGCACA  |       TGTGCAAC       |
| T7_19 |  GCCAATGT  |       ACATTGGC       |
| T7_20 |  TGGTAACA  |       TGTTACCA       |
| T7_21 |  CTTAGAGC  |       GCTCTAAG       |
| T7_22 |  AACCGATA  |       TATCGGTT       |
| T7_23 |  TGCAACTG  |       CAGTTGCA       |

This is the Tn5-T5 oligo sequence: 5'- TCGTCGGCAGCGTCTCCACGCNNNNNNNNGCGATCGAGGACGGCAGATGTGTATAAGAGACAG -3'. The 8 Ns are the barcode, and the barcode sequences (from 5' -> 3') are:

| Oligo | Barcode T5 |  Reverse complement  |
|-------|------------|:--------------------:|
| T5_9  |  TTCGGAAG  |       CTTCCGAA       |
| T5_10 |  AGGCTGGT  |       ACCAGCCT       |
| T5_11 |  GGAAGACT  |       AGTCTTCC       |
| T5_12 |  GAACGCAT  |       ATGCGTTC       |
| T5_13 |  TGCTGGTA  |       TACCAGCA       |
| T5_14 |  CATTCAGT  |       ACTGAATG       |
| T5_15 |  ACAAGGAT  |       ATCCTTGT       |
| T5_16 |  TGTCAGCT  |       AGCTGACA       |

This is the i7 primer sequence: 5'- CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGG -3'. The 8 Ns are the barcode, and the barcode sequences (from 5' -> 3') are:

| Primer | Barcode i7 | Reverse complement |
|--------|------------|:------------------:|
| N701   | TCGCCTTA   |      TAAGGCGA      |
| N702   | CTAGTACG   |      CGTACTAG      |
| N703   | TTCTGCCT   |      AGGCAGAA      |
| N704   | GCTCAGGA   |      TCCTGAGC      |
| N705   | AGGAGTCC   |      GGACTCCT      |
| N706   | CATGCCTA   |      TAGGCATG      |
| N707   | GTAGAGAG   |      CTCTCTAC      |
| N710   | CAGCCTCG   |      CGAGGCTG      |
| N711   | TGCCTCTT   |      AAGAGGCA      |
| N712   | TCCTCTAC   |      GTAGAGGA      |
| N714   | TCATGAGC   |      GCTCATGA      |
| N715   | CCTGAGAT   |      ATCTCAGG      |
| N716   | TAGCGAGT   |      ACTCGCTA      |
| N718   | GTAGCTCC   |      GGAGCTAC      |
| N719   | TACTACGC   |      GCGTAGTA      |
| N721   | GCAGCGTA   |      TACGCTGC      |
| N722   | CTGCGCAT   |      ATGCGCAG      |
| N723   | GAGCGCTA   |      TAGCGCTC      |
| N724   | CGCTCAGT   |      ACTGAGCG      |
| N726   | GTCTTAGG   |      CCTAAGAC      |
| N727   | ACTGATCG   |      CGATCAGT      |
| N728   | TAGCTGCA   |      TGCAGCTA      |
| N729   | GACGTCGA   |      TCGACGTC      |
| X730   | TACCAGAG   |      CTCTGGTA      |
| X731   | GGATGGAA   |      TTCCATCC      |
| X732   | ATTGAGGC   |      GCCTCAAT      |
| X733   | CGGATAGA   |      TCTATCCG      |
| X734   | TGGTAGAC   |      GTCTACCA      |
| X735   | ACCTGGTT   |      AACCAGGT      |
| X736   | CAGTTCTG   |      CAGAACTG      |
| X737   | TCGAACGT   |      ACGTTCGA      |
| X738   | CGTTGCTT   |      AAGCAACG      |
| X739   | TACCGTTC   |      GAACGGTA      |
| X740   | TAGGTTGC   |      GCAACCTA      |
| X741   | GAGGCTAA   |      TTAGCCTC      |
| X742   | CGACCATA   |      TATGGTCG      |
| X743   | AGGCAGTA   |      TACTGCCT      |
| X744   | ATCAAGCG   |      CGCTTGAT      |
| X745   | CATTGAAG   |      CTTCAATG      |
| X746   | CGACTTAT   |      ATAAGTCG      |
| X747   | TCTATACG   |      CGTATAGA      |
| X748   | AGCATTAG   |      CTAATGCT      |
| X749   | AATTGGCA   |      TGCCAATT      |
| X750   | AGATTCGT   |      ACGAATCT      |
| X751   | TTCATGAC   |      GTCATGAA      |
| X752   | TGAACTTG   |      CAAGTTCA      |
| X753   | ATGGCATA   |      TATGCCAT      |
| X754   | CGTAATTC   |      GAATTACG      |

This is the i5 primer sequence: 5'- AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTC -3'. The 8 Ns are the i5 barcode, and the barcode sequences (from 5' -> 3') are:

| Primer | Barcode i5 | Reverse complement |
|--------|------------|:------------------:|
| S502   | CTCTCTAT   |      ATAGAGAG      |
| S503   | TATCCTCT   |      AGAGGATA      |
| S505   | GTAAGGAG   |      CTCCTTAC      |
| S506   | ACTGCATA   |      TATGCAGT      |
| S507   | AAGGAGTA   |      TACTCCTT      |
| S508   | CTAAGCCT   |      AGGCTTAG      |
| S510   | CGTCTAAT   |      ATTAGACG      |
| S511   | TCTCTCCG   |      CGGAGAGA      |
| X512   | TCGACTAG   |      CTAGTCGA      |
| X513   | TTCTAGCT   |      AGCTAGAA      |
| X514   | CCTAGAGT   |      ACTCTAGG      |
| X515   | GCGTAAGA   |      TCTTACGC      |
| X516   | AAGGCTAT   |      ATAGCCTT      |
| X517   | GAGCCTTA   |      TAAGGCTC      |
| X518   | TTATGCGA   |      TCGCATAA      |
| X519   | ATCTGAGT   |      ACTCAGAT      |
| X520   | GGATACTA   |      TAGTATCC      |
| X521   | TAAGATCC   |      GGATCTTA      |
| X522   | AAGAGATG   |      CATCTCTT      |
| X523   | AATGACGT   |      ACGTCATT      |
| X524   | GAAGTATG   |      CATACTTC      |
| X525   | ATAGCCTT   |      AAGGCTAT      |
| X526   | TTGGAAGT   |      ACTTCCAA      |
| X527   | ATTCGTTG   |      CAACGAAT      |
| X528   | AGGATAAC   |      GTTATCCT      |
| X529   | TTCATCCA   |      TGGATGAA      |
| X530   | AACGAACG   |      CGTTCGTT      |
| X531   | TGCCTTAC   |      GTAAGGCA      |
| X532   | CGAATTCC   |      GGAATTCG      |
| X533   | GGTTAGAC   |      GTCTAACC      |
| X534   | TCCGGTAA   |      TTACCGGA      |
| X535   | TTACGACC   |      GGTCGTAA      |

The 32 bp cell barcode consists of the four parts listed above. To construct full whitelist, we need to get all possible combinations of `T7`, `i7`, `T5` and `i5`. In total, we should have a whitelist of __8 * 12 * 32 * 48 = 147456__ barcodes. However, the order of those four parts depends on Illumina machines. I have put those four tables into `csv` files and you can download them to have a look:

[scnATAC-T5.csv](https://teichlab.github.io/scg_lib_structs/data/snATAC-T5.csv)  
[scnATAC-T7.csv](https://teichlab.github.io/scg_lib_structs/data/snATAC-T7.csv)  
[scnATAC-i5.csv](https://teichlab.github.io/scg_lib_structs/data/snATAC-i5.csv)  
[scnATAC-i7.csv](https://teichlab.github.io/scg_lib_structs/data/snATAC-i7.csv)  

### Configuration 1

|  Position | Description                                                                                                                                                            |
|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|     1 - 8 | 8 bp Tn5 barcode at the i7 side. This is the `T7` barcode in Table S4 from the [SnapATAC](https://www.nature.com/articles/s41467-021-21583-9) paper                    |
|    9 - 16 | 8 bp `i7` index. This is the `i7` barcode in the `i7` primer in Table S5 from the [original snATAC-seq paper](https://www.nature.com/articles/s41593-018-0079-3) paper |
|   17 - 24 | 8 bp Tn5 barcode at the i5 side. This is the `T5` barcode in Table S4 from the [SnapATAC](https://www.nature.com/articles/s41467-021-21583-9) paper                    |
|   25 - 32 | 8 bp `i5` index. This is the `i5` barcode in the `i5` primer in Table S5 from the [original snATAC-seq paper](https://www.nature.com/articles/s41593-018-0079-3) paper |

In this configuration, `T7` and `i7` are sequenced using the bottom strand as the template, and `T5` and `i5` are sequenced using the top strand as the template. Therefore, we should take the reverse complementary (rc) in the fowllowing order:

`T7 rc` + `i7 rc` + `T5 rc` + `i5 rc`

Now let's generate the whitelist using the above order:

```bash
# download the tables
wget -P snATAC-seq/data/ \
    https://teichlab.github.io/scg_lib_structs/data/snATAC-T5.csv \
    https://teichlab.github.io/scg_lib_structs/data/snATAC-T7.csv \
    https://teichlab.github.io/scg_lib_structs/data/snATAC-i5.csv \
    https://teichlab.github.io/scg_lib_structs/data/snATAC-i7.csv

# generate combinations of them
for w in $(tail -n +2 snATAC-seq/data/snATAC-T7.csv | cut -f 3 -d,); do
    for x in $(tail -n +2 snATAC-seq/data/snATAC-i7.csv | cut -f 3 -d,); do
        for y in $(tail -n +2 snATAC-seq/data/snATAC-T5.csv | cut -f 3 -d,); do
            for z in $(tail -n +2 snATAC-seq/data/snATAC-i5.csv | cut -f 3 -d,); do
                echo "${w}${x}${y}${z}"
                done
            done
        done
    done > snATAC-seq/data/whitelist_config1.txt
```

### Configuration 2

|  Position | Description                                                                                                                                                            |
|-----------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|     1 - 8 | 8 bp Tn5 barcode at the i7 side. This is the `T7` barcode in Table S4 from the [SnapATAC](https://www.nature.com/articles/s41467-021-21583-9) paper                    |
|    9 - 16 | 8 bp `i7` index. This is the `i7` barcode in the `i7` primer in Table S5 from the [original snATAC-seq paper](https://www.nature.com/articles/s41593-018-0079-3) paper |
|   17 - 24 | 8 bp `i5` index. This is the `i5` barcode in the `i5` primer in Table S5 from the [original snATAC-seq paper](https://www.nature.com/articles/s41593-018-0079-3) paper |
|   25 - 32 | 8 bp Tn5 barcode at the i5 side. This is the `T5` barcode in Table S4 from the [SnapATAC](https://www.nature.com/articles/s41467-021-21583-9) paper                    |

In this configuration, `T7` and `i7` are sequenced using the bottom strand as the template, and we should take the reverse complementary (rc) of them. The `i5` and `T5` are sequenced using the bottom strand as the template as well, and we should take the sequence as they are. Therefore, we should take the index sequence in the fowllowing order:

`T7 rc` + `i7 rc` + `i5` + `T5`

Now let's generate the whitelist using the above order:

```bash
for w in $(tail -n +2 snATAC-seq/data/snATAC-T7.csv | cut -f 3 -d,); do
    for x in $(tail -n +2 snATAC-seq/data/snATAC-i7.csv | cut -f 3 -d,); do
        for y in $(tail -n +2 snATAC-seq/data/snATAC-i5.csv | cut -f 2 -d,); do
            for z in $(tail -n +2 snATAC-seq/data/snATAC-T5.csv | cut -f 2 -d,); do
                echo "${w}${x}${y}${z}"
                done
            done
        done
    done > snATAC-seq/data/whitelist_config2.txt
```

Normally, we are all set and ready to go now. However ... see below ...

### Some Extra work

According to the [SnapATAC](https://www.nature.com/articles/s41467-021-21583-9) paper, the data was generated by __HiSeq 2500__. Therefore, we should use `whitelist_config2.txt` as the whitelist. However, if you check carefully the content in `SRR8592624_4.fastq.gz`, which shold contain `i5` + `T5` barcodes, it seesm the order has been changed to `T5` + `i5`. Maybe the authors changed this before the submission, but I'm not sure. Therefore, for the sake of analysing this specific data set, you need to gerenate a new whitelist in the order of:

`T7 rc` + `i7 rc` + `T5` + `i5`

To this end, we do:

```bash
for w in $(tail -n +2 snATAC-seq/data/snATAC-T7.csv | cut -f 3 -d,); do
    for x in $(tail -n +2 snATAC-seq/data/snATAC-i7.csv | cut -f 3 -d,); do
        for y in $(tail -n +2 snATAC-seq/data/snATAC-T5.csv | cut -f 2 -d,); do
            for z in $(tail -n +2 snATAC-seq/data/snATAC-i5.csv | cut -f 2 -d,); do
                echo "${w}${x}${y}${z}"
                done
            done
        done
    done > snATAC-seq/data/whitelist_config2_GSE126724.txt
```

Now we are ready to go.

## From FastQ To Count Matrix

Now we are ready to map the reads to the genome using `chromap`:

```console
mkdir -p snATAC-seq/chromap_outs

# map and generate the fragment file

chromap -t 4 --preset atac \
        -x mm10/chromap_index/genome.index \
        -r mm10/mm10.fa \
        -1 snATAC-seq/data/SRR8592624_1.fastq.gz \
        -2 snATAC-seq/data/SRR8592624_2.fastq.gz \
        -b snATAC-seq/data/SRR8592624_CB.fastq.gz \
        --barcode-whitelist snATAC-seq/data/whitelist_config2_GSE126724.txt \
        --bc-error-threshold 2 \
        -o snATAC-seq/chromap_outs/fragments.tsv

# compress and index the fragment file

bgzip snATAC-seq/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed snATAC-seq/chromap_outs/fragments.tsv.gz
```

Two new files `fragments.tsv.gz` and `fragments.tsv.gz.tbi` are generated. They will be useful and sometimes required for other programs to perform downstream analysis.

### Explain chromap

If you understand the __snATAC-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/snATAC-seq.html) and the content in the previous sections, the command above should be straightforward to understand.

`-t 4`

>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`-x mm10/chromap_index/genome.index`

>>> The `chromap` index. The data we are analysing is from mouse.

`-r mm10/mm10.fa`

>>> Reference genome sequence in `fasta` format. This is basically the file which you used to create the `chromap` index file.

`-1`, `-2` and `-b`

>>> They are Read 1 (genomic), Read 2 (genomic) and cell barcode read, respectively. For ATAC-seq, the sequencing is usually done in pair-end mode. Therefore, you normally have two genomic reads for each genomic fragment: Read 1 and Read 2. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--barcode-whitelist snATAC-seq/data/whitelist_config2_GSE126724.txt`

>>> For the reason described previously, we are uisng the file `whitelist_config2_GSE126724.txt` specifically generated for this data set. However, in your own case, it should be either `whitelist_config1.txt` or `whitelist_config2.txt`, depending on your sequencing machines.

`--bc-error-threshold 2`

>>> Max Hamming distance allowed to correct a barcode. The default is 1. In this case, we have 32 bp cell barcodes, which is a bit long. Therefore, we would like to allow a bit more errors.

`-o snATAC-seq/chromap_outs/fragments.tsv`

>>> Direct the mapped fragments to a file. The format is described in the [10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments).

## Fragments To Reads

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

```text
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

zcat snATAC-seq/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin mm10/mm10.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > snATAC-seq/chromap_outs/reads.bed.gz
```

Note we also sort the output reads by `sort -k1,1 -k2,2n`. In this way, the order of chromosomes in the `reads.bed.gz` is the same as that in `mm10.chrom.sizes`, which makes downstream processes easier. The output `reads.bed.gz` are the reads in `bed` format, with the 4th column holding the cell barcodes.

## Peak Calling By MACS2

Now we can use the newly generated read file for the peak calling using `MACS2`:

```console
macs2 callpeak -t snATAC-seq/chromap_outs/reads.bed.gz \
               -g mm -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir snATAC-seq/chromap_outs \
               -n aggregate
```

### Explain MACS2

The reasons of choosing those specific parameters are a bit more complicated. I have dedicated a post for this a while ago. Please have a look at [__this post__](https://dbrg77.github.io/posts/2020-12-09-atac-seq-peak-calling-with-macs2/) if you are still confused. The following output files are particularly useful:

| File                       | Description                                                               |
|----------------------------|---------------------------------------------------------------------------|
| aggregate_peaks.narrowPeak | Open chromatin peak locations in the narrowPeak format                    |
| aggregate_peaks.xls        | More information about peaks                                              |
| aggregate_treat_pileup.bdg | Signal tracks. Can be used to generate the bigWig file for visualisation |

## Getting The Peak-By-Cell Count Matrix

Now that we have the peak and reads files, we can compute the number of reads in each peak for each cell. Then we could get the peak-by-cell count matrix. There are different ways of doing this. The following is the method I use.

### Find Reads In Peaks Per Cell

First, we use the `aggregate_peaks.narrowPeak` file. We only need the first 4 columns (chromosome, start, end, peak ID). You can also remove the peaks that overlap [the black list regions](https://www.nature.com/articles/s41598-019-45839-z). The black list is not available for every species and every build, so I'm not doing it here. We also need to sort the peak to make sure the order of the chromosomes in the peak file is the same as that in the `dm6.chrom.sizes` and `reads.bed.gz` files. Then we could find the overlap by `bedtools`. We need to do this in a specific way to get the number of reads in each peak from each cell:

```bash
# format and sort peaks

cut -f 1-4 snATAC-seq/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > snATAC-seq/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a snATAC-seq/chromap_outs/aggregate_peaks_sorted.bed \
    -b snATAC-seq/chromap_outs/reads.bed.gz \
    -wo -sorted -g mm10/mm10.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > snATAC-seq/chromap_outs/peak_read_ov.tsv.gz
```

#### Explain Finding Reads In Peaks Per Cell

We start with the command before the first pipe, that is, the intersection part. If you read the manual of the `bedtools intersect`, it should be straightforward to understand. The `-wo` option will output the records in both `-a` file and `-b` file. Since the `reads.bed.gz` file has the cell barcode information at the 4th column, we would get an output with both peak and cell information for the overlap. The `-sorted -g mm10/mm10.chrom.sizes` options make the program use very little memory. Here is an example (the first 5 lines) of the output of this part:

```text
chr1	3094857	3095503	aggregate_peak_1	chr1	3094812	3094862	CTATTAGGTGCAGCTATGCTGGTACCTAGAGT	.	+	5
chr1	3094857	3095503	aggregate_peak_1	chr1	3094819	3094869	AGAGCAGTATCTCAGGACAAGGATTTATGCGA	.	+	12
chr1	3094857	3095503	aggregate_peak_1	chr1	3094821	3094871	CATGTCAGACTGAGCGACAAGGATGTAAGGAG	.	+	14
chr1	3094857	3095503	aggregate_peak_1	chr1	3094840	3094890	GCTCTAAGCGTATAGATGTCAGCTAATGACGT	.	+	33
chr1	3094857	3095503	aggregate_peak_1	chr1	3094848	3094898	ACATTGGCGCTCATGAGGAAGACTCCTAGAGT	.	+	41
```

We see that the 8th column holds the cell barcode and we want to group them using `bedtools groupby`. Therefore, we need to sort by this column, that is the `sort -k8,8`. When we group by the 8th column, we are interested in how many times each peak appear per group, so we could gather the information of the peak ID (4th column). That is the `-g 8 -c 4 -o freqdesc`. The `-o freqdesc` option returns a `value:frequency` pair in descending order. Here are some records from `peak_read_ov.tsv.gz`:

```text
ACATTGGCAACCAGGTACAAGGATAAGAGATG        aggregate_peak_146375:1
ACATTGGCAACCAGGTACAAGGATCTCTCTAT        aggregate_peak_176781:2,aggregate_peak_44835:2,aggregate_peak_78312:2,aggregate_peak_87937:2
```

In a way, that is sort of a count matrix in an awkward format. For example:

- The first line means that in cell `ACATTGGCAACCAGGTACAAGGATAAGAGATG`, the peak `aggregate_peak_146375` has 1 count. All the rest peaks not mentioned here have 0 counts in this cell.
- The second line means that in cell `ACATTGGCAACCAGGTACAAGGATCTCTCTAT`, the peak `aggregate_peak_176781` has 2 counts, the peak `aggregate_peak_44835` has 2 counts, the peak `aggregate_peak_78312` has two peaks and the peak `aggregate_peak_87937` has 2 counts. All the rest peaks not mentioned here have 0 counts in this cell.

### Output The Peak-By-Cell Matrix

At this stage, we pretty much have all the things needed. Those two files `aggregate_peaks_sorted.bed` and `peak_read_ov.tsv.gz` contain all information for a peak-by-cell count matrix. We just need a final touch to make the output in a standard format: a [market exchange format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). Since most downstream software takes the input from the __10x Genomics Single Cell ATAC__ results, we are going to generate the MEX and the associated files similar to the output from 10x Genomics.

Here, I'm using a python script for this purpose. You don't have to do this. Choose whatever works for you. The point here is to just generate similar files as the __peak-barcode matrix__ described from [the 10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

First, let's make a directory to hold the output files and generate the `peaks.bed` and `barcodes.tsv` files, which are easy to do:

```bash
# create dirctory
mkdir -p snATAC-seq/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 snATAC-seq/chromap_outs/aggregate_peaks_sorted.bed > \
    snATAC-seq/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat snATAC-seq/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    snATAC-seq/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
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
    snATAC-seq/chromap_outs/aggregate_peaks_sorted.bed \
    snATAC-seq/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    snATAC-seq/chromap_outs/peak_read_ov.tsv.gz \
    snATAC-seq/chromap_outs/raw_peak_bc_matrix
```

After that, you should have the `matrix.mtx` in the `snATAC-seq/chromap_outs/raw_peak_bc_matrix` directory.

### Cell Calling (Filter Cell Barcodes)

Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, we could use `starolo`. It is used for scRNA-seq in general, but it does have a cell calling function that takes a directory containing raw mtx and associated files, and return the filtered ones. Since `starsolo` looks for the following three files in the input directory: `matrix.mtx`, `features.tsv` and `barcodes.tsv`. Those are the output from the 10x Genomics scRNA-seq workflow. In this case, we can use `peaks.bed` as our `features.tsv`:

``console
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed snATAC-seq/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     snATAC-seq/chromap_outs/raw_peak_bc_matrix \
     snATAC-seq/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv snATAC-seq/chromap_outs/filtered_peak_bc_matrix/peaks.bed
```

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/snATAC-seq/
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
└── data
    ├── snATAC-i5.csv
    ├── snATAC-i7.csv
    ├── snATAC-T5.csv
    ├── snATAC-T7.csv
    ├── SRR8592624_1.fastq.gz
    ├── SRR8592624_2.fastq.gz
    ├── SRR8592624_3.fastq.gz
    ├── SRR8592624_4.fastq.gz
    ├── SRR8592624_CB.fastq.gz
    ├── whitelist_config1.txt
    ├── whitelist_config2_GSE126724.txt
    └── whitelist_config2.txt

4 directories, 30 files
```