# sci-ATAC-seq

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-ATAC-seq.html) to see how __sci-ATAC-seq__ libraries are generated experimentally. This is a split-pool based combinatorial indexing strategy, where nuclei are transposed in mini-bulk with indexed transposase Tn5 (__C15__ + __D15__). Then all nuclei are pooled and a fraction of the nuclei is randomly distributed into new wells. Library preparation is performed using a PCR primers with different indices (__i5__ + __i7__). Single cells can be identified by the combination of __D15 + i7 + C15 + i5__. You might wonder what `C15` and `D15` are, those are the primer in the published paper. Keep reading and you will see.

```{eval-rst}
.. important::

   The library of the **sci-ATAC-seq** is kind of complicated, especially for beginners. First, make sure you are familiar with different sequencing modes from different Illumina machines by looking at `this page <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_. Then make sure you understand how the actual sequencing is done for **sci-ATAC-seq** by checking `this GitHub page <https://teichlab.github.io/scg_lib_structs/methods_html/sci-ATAC-seq.html>`_ and the associated publications.

```

## For Your Own Experiments

If you use this assay, you have to run the sequencing by yourself using a custom sequencing recipe or ask your core facility to do this for you. We don't have the proper setup to do this experiment here, so I haven't tried __sci-ATAC-seq__ by myself. This is the educational guess of the sequencing configuration. Using more recent machines and chemistries, like __iSeq 100__, __MiniSeq__, __NextSeq__, __HiSeq X__, __HiSeq 3000__, __HiSeq 4000__ and __NovaSeq 600 (v1.5)__, it should be (I call this __Configuration 1__):

| Order | Read             | Cycle                                 | Description                       |
|-------|------------------|---------------------------------------|-----------------------------------|
| 1     | Read 1           | >50                                   | `R1_001.fastq.gz`, Genomic insert |
| 2     | Index 1 (__i7__) | 8 normal + 27 dark + 10 normal cycles | `I1_001.fastq.gz`, D15 + i7 index |
| 3     | Index 2 (__i5__) | 8 normal + 21 dark + 10 normal cycles | `I2_001.fastq.gz`, C15 + i5 index |
| 4     | Read 2           | >50                                   | `R2_001.fastq.gz`, Genomic insert |

Using older machines and chemistries like __MiSeq__, __HiSeq 2000__, __HiSeq 2500__, __NovaSeq 6000 (v1.0)__, it should be (I call this __Configuration 2__):

| Order | Read             | Cycle                                 | Description                       |
|-------|------------------|---------------------------------------|-----------------------------------|
| 1     | Read 1           | >50                                   | `R1_001.fastq.gz`, Genomic insert |
| 2     | Index 1 (__i7__) | 8 normal + 27 dark + 10 normal cycles | `I1_001.fastq.gz`, D15 + i7 index |
| 3     | Index 2 (__i5__) | 10 normal + 21 dark + 8 normal cycles | `I2_001.fastq.gz`, i5 + C15 index |
| 4     | Read 2           | >50                                   | `R2_001.fastq.gz`, Genomic insert |

The reason for the dark cycles is to avoid sequencing those common adaptors with exact the same DNA sequences. Again, if you are not sure, check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-ATAC-seq.html) to see the reason.

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
Undetermined_S0_I1_001.fastq.gz  # 16 bp, D15 + i7
Undetermined_S0_I2_001.fastq.gz  # 16 bp, C15 + i5
Undetermined_S0_R1_001.fastq.gz  # >50 bp, genomic
Undetermined_S0_R2_001.fastq.gz  # >50 bp, genomic
```

Your cell barcodes will be `I1` + `I2`. We cloud stitch them together to get the cell barcodes to a single `fastq` file:

```bash
paste <(zcat Undetermined_S0_I1_001.fastq.gz) <(zcat Undetermined_S0_I2_001.fastq.gz) | \
    awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
    gzip > Undetermined_S0_CB_001.fastq.gz
```

After that, you are ready to throw `Undetermined_S0_R1_001.fastq.gz`, `Undetermined_S0_R2_001.fastq.gz` and `Undetermined_S0_CB_001.fastq.gz` into `chromap`. See the later section.

## Public Data

For the purpose of demonstration, we will use the __sci-ATAC-seq__ data from the following paper:

```{eval-rst}
.. note::

  Cusanovich DA, Reddington JP, Garfield DA, Daza RM, Aghamirzaie D, Marco-Ferreres R, Pliner HA, Christiansen L, Qiu X, Steemers FJ, Trapnell C, Shendure J, Furlong EEM (2018) **The cis-regulatory dynamics of embryonic development at single-cell resolution.** *Nature* 555:538–542. https://doi.org/10.1038/nature25981

```

where the authors used __sci-ATAC-seq__ to profile the chromatin accessibility of Drosophila embryos during three landmark embryonic stages. The data is in GEO under the accession code [GSE101581](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101581). The raw reads are in __SRR5837698__. To get the reads, we use `fastq-dump` from the `sratools`:

```console
mkdir -p sci-atac/data
fastq-dump --split-files \
           --origfmt \
           --outdir sci-atac/data \
           --defline-seq '@rd.$si:$sg:$sn' SRR5837698
```

The reason of using the specific options above is a bit complicated. I wrote [__a post__](https://dbrg77.github.io/posts/2022-07-26-getting-index-reads-from-sra/) about this, and you can have a look to see the reason. Once the program finishes running, you will have two plain files: `SRR5837698_1.fastq` and `SRR5837698_1.fastq`, each of which has a size of 86G. Unfortunately, it seems `fastq-dump` cannot generate gzipped file on the fly. I'm going to gzip them manually to save space:

```console
gzip sci-atac/data/SRR5837698_1.fastq sci-atac/data/SRR5837698_2.fastq
```

Anyway, let's have a look at the first 5 reads from Read 1:

```text
@rd.1::ATTCAGAACCGCTAAGAGNAAGATTATTAGATTCCG:1
GGCTTNTATTATGACCGCAATGAAGTCCGATCGCAGATAATCCGCAAAGGA
+ATTCAGAACCGCTAAGAGNAAGATTATTAGATTCCG:1
AAAAA#EEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@rd.2::CGCTCATTATACCTCGACNGGCTATAGATACCTAAG:2
CTTCTNTTCTCTTCTCTTCTCTTCTCTTCTCTTCTCTTCACTTCTATTCTA
+CGCTCATTATACCTCGACNGGCTATAGATACCTAAG:2
AAAAA#EEEAEEEEEEEEEEAEEAEEEEE/EEEAE</E/////E</////A
@rd.3::CTGAAGCTTCCTCTGCCGNGGCTATAAAGCCGGCTG:3
GTAGTNGAGTGTATGGCCCACGCCGTAGCACCAGAAAAAAATTGGGTCGTT
+CTGAAGCTTCCTCTGCCGNGGCTATAAAGCCGGCTG:3
AAAAA#EEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@rd.4::GAGATTCCGACTGGACCANCCTCTATATAAGAGGCG:4
ACCGANGAGGGTTGCGTCATTCGTCCTGTACGCTCCGCATTGTCCAGTAAA
+GAGATTCCGACTGGACCANCCTCTATATAAGAGGCG:4
AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
@rd.5::TAATGCGCGCGCGTTCATNTCAGTACGGTCCAGGAG:5
ATTTTNGATCAAATCCACATTCAAATGGGGATCTTTTTTCTCGGTCGATTA
+TAATGCGCGCGCGTTCATNTCAGTACGGTCCAGGAG:5
AAAAA#E/EEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEE
```

As you can see from the header, the __8 bp + 8 bp + 8 bp + 8 bp = 32 bp__ cell barcodes are in the 3rd field with `:` as the field separator. Unfortunately, the quality strings for them are lost, but we can just use them as if they are perfect. Now we need to generate a new `fastq` files holding these 32 bp cell barcodes with high Phred scores, say 40. The corresponding ASCII encoding is 40 + 33 = 73, which is `I`. You can use the header from either `SRR5837698_1.fastq.gz` or `SRR5837698_2.fastq.gz`, it does not matter:

```bash
zcat sci-atac/data/SRR5837698_1.fastq.gz | \
    awk 'NR%4==1' | \
    awk -F ':' '{print $0 "\n" $3 "\n+\n" "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"}' | \
    gzip > sci-atac/data/SRR5837698_CB.fastq.gz
```

Now we have the three `fastq` files `SRR5837698_1.fastq.gz`, `SRR5837698_2.fastq.gz` and `SRR5837698_CB.fastq.gz`. That's all we need for `chromap`.

## Prepare Whitelist

This is the most confusing part of the pipeline, but if you understand the procedures in the __sci-ATAC-seq__ [GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-ATAC-seq.html), you should be able to understand.

```{eval-rst}
.. important::

  For the exact primer sequences, check `Supplementary Table 12 <https://teichlab.github.io/scg_lib_structs/data/Cusanovich2018_Table_S12.xlsx>`_ from the Cusanovich 2018 *Nature* paper.

```

This is the D15 primer sequence: 5'- GTCTCGTGGGCTCGGCTGTCCCTGTCCNNNNNNNNCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG -3'. The 8 Ns are the barcode, and the barcode sequences (from 5' -> 3') are:

| Primer | Barcode D15 |  Reverse complement  |
|--------|-------------|:--------------------:|
| D15_1  | CGAGTAAT    |       ATTACTCG       |
| D15_2  | TCTCCGGA    |       TCCGGAGA       |
| D15_3  | AATGAGCG    |       CGCTCATT       |
| D15_4  | GGAATCTC    |       GAGATTCC       |
| D15_5  | TTCTGAAT    |       ATTCAGAA       |
| D15_6  | ACGAATTC    |       GAATTCGT       |
| D15_7  | AGCTTCAG    |       CTGAAGCT       |
| D15_8  | GCGCATTA    |       TAATGCGC       |
| D15_9  | CATAGCCG    |       CGGCTATG       |
| D15_10 | TTCGCGGA    |       TCCGCGAA       |
| D15_11 | GCGCGAGA    |       TCTCGCGC       |
| D15_12 | CTATCGCT    |       AGCGATAG       |

This is the C15 primer sequence: 5'- TCGTCGGCAGCGTCTCCACGCNNNNNNNNGCGATCGAGGACGGCAGATGTGTATAAGAGACAG -3'. The 8 Ns are the barcode, and the barcode sequences (from 5' -> 3') are:

| Primer | Barcode C15 |  Reverse complement  |
|--------|-------------|:--------------------:|
| C15_1  | TATAGCCT    |       AGGCTATA       |
| C15_2  | ATAGAGGC    |       GCCTCTAT       |
| C15_3  | CCTATCCT    |       AGGATAGG       |
| C15_4  | GGCTCTGA    |       TCAGAGCC       |
| C15_5  | AGGCGAAG    |       CTTCGCCT       |
| C15_6  | TAATCTTA    |       TAAGATTA       |
| C15_7  | CAGGACGT    |       ACGTCCTG       |
| C15_8  | GTACTGAC    |       GTCAGTAC       |

This is the P7 primer sequence: 5'- CAAGCAGAAGACGGCATACGAGATNNNNNNNNNNGTCTCGTGGGCTCGG -3'. The 10 Ns are the i7 barcode, and the barcode sequences (from 5' -> 3') are:

| Primer | Barcode i7 | Reverse complement |
|--------|------------|:------------------:|
| P7_1   | CCGAATCCGA |     TCGGATTCGG     |
| P7_2   | CCGCAGCCGC |     GCGGCTGCGG     |
| P7_3   | AACGTAATCT |     AGATTACGTT     |
| P7_4   | ACCTAGTTAG |     CTAACTAGGT     |
| P7_5   | GGTCGCTATG |     CATAGCGACC     |
| P7_6   | CTCTTAGCGG |     CCGCTAAGAG     |
| P7_7   | TTCGTTCCAT |     ATGGAACGAA     |
| P7_8   | AACGGAACGC |     GCGTTCCGTT     |
| P7_9   | TTCGATAACC |     GGTTATCGAA     |
| P7_10  | CATACGATGC |     GCATCGTATG     |
| P7_11  | TTATCGTATT |     AATACGATAA     |
| P7_12  | GTCGACGGAA |     TTCCGTCGAC     |
| P7_13  | ATAAGCCGGA |     TCCGGCTTAT     |
| P7_14  | TGCGCCTGGT |     ACCAGGCGCA     |
| P7_15  | ATTCTCCTCT |     AGAGGAGAAT     |
| P7_16  | ATAGGAGTAC |     GTACTCCTAT     |
| P7_17  | ATCCGTTAGC |     GCTAACGGAT     |
| P7_18  | TGATTCAACT |     AGTTGAATCA     |
| P7_19  | TACCTAATCA |     TGATTAGGTA     |
| P7_20  | GATGCTACGA |     TCGTAGCATC     |
| P7_21  | AACCTCAAGA |     TCTTGAGGTT     |
| P7_22  | AAGCTGACCT |     AGGTCAGCTT     |
| P7_23  | AAGTCTAATA |     TATTAGACTT     |
| P7_24  | ACTAATTGAG |     CTCAATTAGT     |
| P7_25  | CCGGCGGCGA |     TCGCCGCCGG     |
| P7_26  | AATCATACGG |     CCGTATGATT     |
| P7_27  | TCTGCGCGTT |     AACGCGCAGA     |
| P7_28  | CTACGACGAG |     CTCGTCGTAG     |
| P7_29  | TCGCAATTAG |     CTAATTGCGA     |
| P7_30  | TATGGCCGCG |     CGCGGCCATA     |
| P7_31  | AAGTAATATT |     AATATTACTT     |
| P7_32  | ATCTGCCAAT |     ATTGGCAGAT     |
| P7_33  | CAGGCGCCAT |     ATGGCGCCTG     |
| P7_34  | GAGTCCTTAT |     ATAAGGACTC     |
| P7_35  | CGGCTTACTA |     TAGTAAGCCG     |
| P7_36  | CTTGCATAAT |     ATTATGCAAG     |
| P7_37  | GGCTTGCCAA |     TTGGCAAGCC     |
| P7_38  | CGCCAATCAA |     TTGATTGGCG     |
| P7_39  | GCTCATATGC |     GCATATGAGC     |
| P7_40  | AGTCGAGTTC |     GAACTCGACT     |
| P7_41  | GGCTGGCTAG |     CTAGCCAGCC     |
| P7_42  | AGAGGTCGCA |     TGCGACCTCT     |
| P7_43  | AGCTAAGAAT |     ATTCTTAGCT     |
| P7_44  | ATCGTATCAA |     TTGATACGAT     |
| P7_45  | AACTATTATA |     TATAATAGTT     |
| P7_46  | CCTACGGCAA |     TTGCCGTAGG     |
| P7_47  | GATATGGTCT |     AGACCATATC     |
| P7_48  | TCCTTACCAA |     TTGGTAAGGA     |
| P7_49  | CCGCTAGCTG |     CAGCTAGCGG     |
| P7_50  | CAAGGCTTAG |     CTAAGCCTTG     |
| P7_51  | AGCGGTAACG |     CGTTACCGCT     |
| P7_52  | TGGTCCAGTC |     GACTGGACCA     |
| P7_53  | ACGGTCTTGC |     GCAAGACCGT     |
| P7_54  | AGGAGATTGA |     TCAATCTCCT     |
| P7_55  | GTCGAGGTAT |     ATACCTCGAC     |
| P7_56  | AACGCCTCTA |     TAGAGGCGTT     |
| P7_57  | AAGTTACCTA |     TAGGTAACTT     |
| P7_58  | AATATTCGAA |     TTCGAATATT     |
| P7_59  | TAGTCGTCCA |     TGGACGACTA     |
| P7_60  | TGCAGCCTAC |     GTAGGCTGCA     |
| P7_61  | CTTATCCTAC |     GTAGGATAAG     |
| P7_62  | GCGCTCGACG |     CGTCGAGCGC     |
| P7_63  | AATGAATAGT |     ACTATTCATT     |
| P7_64  | ATCTAAGCAA |     TTGCTTAGAT     |
| P7_65  | GCTCCATTCG |     CGAATGGAGC     |
| P7_66  | GGCTATATAG |     CTATATAGCC     |
| P7_67  | TTATTAGTAG |     CTACTAATAA     |
| P7_68  | ACGGCAACCA |     TGGTTGCCGT     |
| P7_69  | CGGCAGAGGA |     TCCTCTGCCG     |
| P7_70  | TTCAAGAATC |     GATTCTTGAA     |
| P7_71  | TAGCTGCTAC |     GTAGCAGCTA     |
| P7_72  | GGAGCTGAGG |     CCTCAGCTCC     |
| P7_73  | TGAGCTACTT |     AAGTAGCTCA     |
| P7_74  | TCCAGCAATA |     TATTGCTGGA     |
| P7_75  | CCGTATCTGG |     CCAGATACGG     |
| P7_76  | CGAATTCGTT |     AACGAATTCG     |
| P7_77  | ACGATAAGCG |     CGCTTATCGT     |
| P7_78  | TCGCGTACTT |     AAGTACGCGA     |
| P7_79  | TGCGAAGATC |     GATCTTCGCA     |
| P7_80  | CAGGCTAAGA |     TCTTAGCCTG     |
| P7_81  | GCCTCAATAA |     TTATTGAGGC     |
| P7_82  | ATGCTCGCAA |     TTGCGAGCAT     |
| P7_83  | CTCTTCAAGC |     GCTTGAAGAG     |
| P7_84  | GCAGCGGACT |     AGTCCGCTGC     |
| P7_85  | TCAGGACTTA |     TAAGTCCTGA     |
| P7_86  | CATGAGAACT |     AGTTCTCATG     |
| P7_87  | CCTTAGTCTG |     CAGACTAAGG     |
| P7_88  | CAGCGATAGA |     TCTATCGCTG     |
| P7_89  | ACCATAGCGC |     GCGCTATGGT     |
| P7_90  | AATAATAATG |     CATTATTATT     |
| P7_91  | AACTACGGCT |     AGCCGTAGTT     |
| P7_92  | CGCAATATCA |     TGATATTGCG     |
| P7_93  | TTAACGCCGT |     ACGGCGTTAA     |
| P7_94  | GGAGTAAGCC |     GGCTTACTCC     |
| P7_95  | ATGAACGCGC |     GCGCGTTCAT     |
| P7_96  | CATCGCGCTC |     GAGCGCGATG     |

This is the P5 primer sequence: 5'- AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNTCGTCGGCAGCGTC -3'. The 10 Ns are the i5 barcode, and the barcode sequences (from 5' -> 3') are:

| Primer | Barcode i5 | Reverse complement |
|--------|------------|:------------------:|
| P5_1   | CTCCATCGAG |     CTCGATGGAG     |
| P5_2   | TTGGTAGTCG |     CGACTACCAA     |
| P5_3   | GGCCGTCAAC |     GTTGACGGCC     |
| P5_4   | CCTAGACGAG |     CTCGTCTAGG     |
| P5_5   | TCGTTAGAGC |     GCTCTAACGA     |
| P5_6   | CGTTCTATCA |     TGATAGAACG     |
| P5_7   | CGGAATCTAA |     TTAGATTCCG     |
| P5_8   | ATGACTGATC |     GATCAGTCAT     |
| P5_9   | TCAATATCGA |     TCGATATTGA     |
| P5_10  | GTAGACCTGG |     CCAGGTCTAC     |
| P5_11  | TTATGACCAA |     TTGGTCATAA     |
| P5_12  | TTGGTCCGTT |     AACGGACCAA     |
| P5_13  | GGTACGTTAA |     TTAACGTACC     |
| P5_14  | CAATGAGTCC |     GGACTCATTG     |
| P5_15  | GATGCAGTTC |     GAACTGCATC     |
| P5_16  | CCATCGTTCC |     GGAACGATGG     |
| P5_17  | TTGAGAGAGT |     ACTCTCTCAA     |
| P5_18  | ACTGAGCGAC |     GTCGCTCAGT     |
| P5_19  | TGAGGAATCA |     TGATTCCTCA     |
| P5_20  | CCTCCGACGG |     CCGTCGGAGG     |
| P5_21  | CATTGACGCT |     AGCGTCAATG     |
| P5_22  | TCGTCCTTCG |     CGAAGGACGA     |
| P5_23  | TGATACTCAA |     TTGAGTATCA     |
| P5_24  | TTCTACCTCA |     TGAGGTAGAA     |
| P5_25  | TCGTCGGAAC |     GTTCCGACGA     |
| P5_26  | ATCGAGATGA |     TCATCTCGAT     |
| P5_27  | TAGACTAGTC |     GACTAGTCTA     |
| P5_28  | GTCGAAGCAG |     CTGCTTCGAC     |
| P5_29  | AGGCGCTAGG |     CCTAGCGCCT     |
| P5_30  | AGATGCAACT |     AGTTGCATCT     |
| P5_31  | AAGCCTACGA |     TCGTAGGCTT     |
| P5_32  | GTAGGCAATT |     AATTGCCTAC     |
| P5_33  | GGAGGCGGCG |     CGCCGCCTCC     |
| P5_34  | CCAGTACTTG |     CAAGTACTGG     |
| P5_35  | GGTCTCGCCG |     CGGCGAGACC     |
| P5_36  | GGCGGAGGTC |     GACCTCCGCC     |
| P5_37  | TAGTTCTAGA |     TCTAGAACTA     |
| P5_38  | TTGGAGTTAG |     CTAACTCCAA     |
| P5_39  | AGATCTTGGT |     ACCAAGATCT     |
| P5_40  | GTAATGATCG |     CGATCATTAC     |
| P5_41  | CAGAGAGGTC |     GACCTCTCTG     |
| P5_42  | TTAATTAGCC |     GGCTAATTAA     |
| P5_43  | CTCTAACTCG |     CGAGTTAGAG     |
| P5_44  | TACGATCATC |     GATGATCGTA     |
| P5_45  | AGGCGAGAGC |     GCTCTCGCCT     |
| P5_46  | TCAAGATAGT |     ACTATCTTGA     |
| P5_47  | TAATTGACCT |     AGGTCAATTA     |
| P5_48  | CAGCCGGCTT |     AAGCCGGCTG     |
| P5_49  | AGAACCGGAG |     CTCCGGTTCT     |
| P5_50  | GAGATGCATG |     CATGCATCTC     |
| P5_51  | GATTACCGGA |     TCCGGTAATC     |
| P5_52  | TCGTAACGGT |     ACCGTTACGA     |
| P5_53  | TGGCGACGGA |     TCCGTCGCCA     |
| P5_54  | AGTCATAGCC |     GGCTATGACT     |
| P5_55  | GTCAAGTCCA |     TGGACTTGAC     |
| P5_56  | ATTCGGAAGT |     ACTTCCGAAT     |
| P5_57  | GTCGGTAGTT |     AACTACCGAC     |
| P5_58  | AGGACGGACG |     CGTCCGTCCT     |
| P5_59  | CTCCTGGACC |     GGTCCAGGAG     |
| P5_60  | TAGCCTCGTT |     AACGAGGCTA     |
| P5_61  | GGTTGAACGT |     ACGTTCAACC     |
| P5_62  | AGGTCCTCGT |     ACGAGGACCT     |
| P5_63  | GGAAGTTATA |     TATAACTTCC     |
| P5_64  | TGGTAATCCT |     AGGATTACCA     |
| P5_65  | AAGCTAGGTT |     AACCTAGCTT     |
| P5_66  | TCCGCGGACT |     AGTCCGCGGA     |
| P5_67  | TGCGGATAGT |     ACTATCCGCA     |
| P5_68  | TGGCAGCTCG |     CGAGCTGCCA     |
| P5_69  | TGCTACGGTC |     GACCGTAGCA     |
| P5_70  | GCGCAATGAC |     GTCATTGCGC     |
| P5_71  | CTTAATCTTG |     CAAGATTAAG     |
| P5_72  | GGAGTTGCGT |     ACGCAACTCC     |
| P5_73  | ACTCGTATCA |     TGATACGAGT     |
| P5_74  | GGTAATAATG |     CATTATTACC     |
| P5_75  | TCCTTATAGA |     TCTATAAGGA     |
| P5_76  | CCGACTCCAA |     TTGGAGTCGG     |
| P5_77  | GCCAAGCTTG |     CAAGCTTGGC     |
| P5_78  | CATATCCTAT |     ATAGGATATG     |
| P5_79  | ACCTACGCCA |     TGGCGTAGGT     |
| P5_80  | GGAATTCAGT |     ACTGAATTCC     |
| P5_81  | TGGCGTAGAA |     TTCTACGCCA     |
| P5_82  | ATTGCGGCCA |     TGGCCGCAAT     |
| P5_83  | TTCAGCTTGG |     CCAAGCTGAA     |
| P5_84  | CCATCTGGCA |     TGCCAGATGG     |
| P5_85  | CTTATAAGTT |     AACTTATAAG     |
| P5_86  | GATTAGATGA |     TCATCTAATC     |
| P5_87  | TATAGGATCT |     AGATCCTATA     |
| P5_88  | AGCTTATAGG |     CCTATAAGCT     |
| P5_89  | GTCTGCAATC |     GATTGCAGAC     |
| P5_90  | CGCCTCTTAT |     ATAAGAGGCG     |
| P5_91  | GTTGGATCTT |     AAGATCCAAC     |
| P5_92  | GCGATTGCAG |     CTGCAATCGC     |
| P5_93  | TGCCAGTTGC |     GCAACTGGCA     |
| P5_94  | CTTAGGTATC |     GATACCTAAG     |
| P5_95  | GAGACCTACC |     GGTAGGTCTC     |
| P5_96  | ATTGACCGAG |     CTCGGTCAAT     |

The 36 bp cell barcode consists of the four parts listed above. To construct full whitelist, we need to get all possible combinations of `D15`, `i7`, `C15` and `i5`. In total, we should have a whitelist of __8 * 12 * 96 * 96 = 884736__ barcodes. However, the order of those four parts depends on Illumina machines. I have put those four tables into `csv` files and you can download them to have a look:

[sci-ATAC-C15.csv](https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-C15.csv)  
[sci-ATAC-D15.csv](https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-D15.csv)  
[sci-ATAC-P5.csv](https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-P5.csv)  
[sci-ATAC-P7.csv](https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-P7.csv)  

### Configuration 1

|  Position | Description                                                                                            |
|-----------|--------------------------------------------------------------------------------------------------------|
|     1 - 8 | 8 bp Tn5 barcode at the `i7` side. This is the barcode in the `D15` primer in Table S12 from the paper |
|    9 - 18 | 10 bp `i7` index. This is the barcode in the P7 primer in Table S12 from the paper                     |
|   19 - 26 | 8 bp Tn5 barcode at the `i5` side. This is the barcode in the `C15` primer in Table S12 from the paper |
|   27 - 36 | 10 bp `i5` index. This is the barcode in the P5 primer in Table S12 from the paper                     |

In this configuration, `D15` and `i7` are sequenced using the bottom strand as the template, and `C15` and `i5` are sequenced using the top strand as the template. Therefore, we should take the reverse complementary (rc) in the fowllowing order:

`D15 rc` + `i7 rc` + `C15 rc` + `i5 rc`

Now let's generate the whitelist using the above order:

```bash
# download the tables
wget -P sci-atac/data/ \
    https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-C15.csv \
    https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-D15.csv \
    https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-P5.csv \
    https://teichlab.github.io/scg_lib_structs/data/sci-ATAC-P7.csv

# generate combinations of them
for w in $(tail -n +2 sci-atac/data/sci-ATAC-D15.csv | cut -f 3 -d,); do
    for x in $(tail -n +2 sci-atac/data/sci-ATAC-P7.csv | cut -f 3 -d,); do
        for y in $(tail -n +2 sci-atac/data/sci-ATAC-C15.csv | cut -f 3 -d,); do
            for z in $(tail -n +2 sci-atac/data/sci-ATAC-P5.csv | cut -f 3 -d,); do
                echo "${w}${x}${y}${z}"
                done
            done
        done
    done > sci-atac/data/whitelist_config1.txt
```

### Configuration 2

|  Position | Description                                                                                            |
|-----------|--------------------------------------------------------------------------------------------------------|
|     1 - 8 | 8 bp Tn5 barcode at the `i7` side. This is the barcode in the `D15` primer in Table S12 from the paper |
|    9 - 18 | 10 bp `i7` index. This is the barcode in the P7 primer in Table S12 from the paper                     |
|   19 - 28 | 10 bp `i5` index. This is the barcode in the P5 primer in Table S12 from the paper                     |
|   29 - 36 | 8 bp Tn5 barcode at the `i5` side. This is the barcode in the `C15` primer in Table S12 from the paper |

In this configuration, `D15` and `i7` are sequenced using the bottom strand as the template, and we should take the reverse complementary (rc) of them. The `C15` and `i5` are sequenced using the bottom strand as the template as well, and we should take the sequence as they are. Therefore, we should take the index sequence in the fowllowing order:

`D15 rc` + `i7 rc` + `i5` + `C15`

Now let's generate the whitelist using the above order:

```bash
for w in $(tail -n +2 sci-atac/data/sci-ATAC-D15.csv | cut -f 3 -d,); do
    for x in $(tail -n +2 sci-atac/data/sci-ATAC-P7.csv | cut -f 3 -d,); do
        for y in $(tail -n +2 sci-atac/data/sci-ATAC-P5.csv | cut -f 2 -d,); do
            for z in $(tail -n +2 sci-atac/data/sci-ATAC-C15.csv | cut -f 2 -d,); do
                echo "${w}${x}${y}${z}"
                done
            done
        done
    done > sci-atac/data/whitelist_config2.txt
```

Normally, we are all set and ready to go now. However ... see below ...

### Some Extra Work

At this time of writing (31-July-2022), `chromap` only supports cell barcodes with `<32 bp` in length, based on [this thread](https://github.com/haowenz/chromap/issues/103). Therefore, we need a truncated version of the cell barcode and whitelist. We have a total of __884736__ valid cell barcodes with 36 bp. It turns out that the positions 3 - 33 (31 bp) are already unique. We could just use those 31 bp as our cell barcodes and whitelist. To this end, we do:

```bash
# truncate whitelist
cut -c 3-33 sci-atac/data/whitelist_config1.txt > sci-atac/data/whitelist_config1_base3-33.txt

# truncate cell barcode reads
zcat sci-atac/data/SRR5837698_CB.fastq.gz | \
    awk '{ if (NR%4==1||NR%4==3) {print $1} else {print substr($1, 3, 31)} }' | \
    gzip > sci-atac/data/SRR5837698_CB_base3-33.fastq.gz
```

Now we are ready to go.

## Index Genome

The data we are using here is from Drosophila embryos. We need to download the genome sequence and index it with `chromap`. In this case, we are using the `dm6` build from UCSC:

```console
mkdir -p dm6/chromap_index
wget -P dm6 -c https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
gunzip dm6/dm6.fa.gz
chromap -i -t 4 -r dm6/dm6.fa -o dm6/chromap_index/genome.index
```

## From FastQ To Count Matrix

Now we are ready to map the reads to the genome using `chromap`:

```console
mkdir -p sci-atac/chromap_outs

# map and generate the fragment file

chromap -t 4 --preset atac \
        -x dm6/chromap_index/genome.index \
        -r dm6/dm6.fa \
        -1 sci-atac/data/SRR5837698_1.fastq.gz \
        -2 sci-atac/data/SRR5837698_2.fastq.gz \
        -b sci-atac/data/SRR5837698_CB_base3-33.fastq.gz \
        --barcode-whitelist sci-atac/data/whitelist_config1_base3-33.txt \
        -o sci-atac/chromap_outs/fragments.tsv

# compress and index the fragment file

bgzip sci-atac/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed sci-atac/chromap_outs/fragments.tsv.gz
```

Two new files `fragments.tsv.gz` and `fragments.tsv.gz.tbi` are generated. They will be useful and sometimes required for other programs to perform downstream analysis.

### Explain chromap

If you understand the __sci-ATAC-seq__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/sci-ATAC-seq.html) and the content in the previous sections, the command above should be straightforward to understand.

`-t 4`

>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`-x dm6/chromap_index/genome.index`

>>> The `chromap` index file. We are dealing with fly samples in this case.

`-r dm6/dm6.fa`

>>> Reference genome sequence in `fasta` format. This is basically the file which you used to create the `chromap` index file.

`-1`, `-2` and `-b`

>>> They are Read 1 (genomic), Read 2 (genomic) and cell barcode read, respectively. For ATAC-seq, the sequencing is usually done in pair-end mode. Therefore, you normally have two genomic reads for each genomic fragment: Read 1 and Read 2. For the reason described previously, we are uisng the truncated version (31 bp) of the cell barcode reads here. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--barcode-whitelist sci-atac/data/whitelist_config1_base3-33.txt`

>>> The plain text file containing all possible valid cell barcodes. For the reason described in the previous sections, we are using the truncated version (31 bp) of the whitelist because the sequencing machines used to generate the data was __NextSeq__ based on the publication.

`-o sci-atac/chromap_outs/fragments.tsv`

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

faSize -detailed dm6/dm6.fa | \
    sort -k1,1 > dm6/dm6.chrom.sizes
```

This is the first 5 lines of `dm6/dm6.chrom.sizes`:

```
chr2L	23513712
chr2R	25286936
chr3L	28110227
chr3R	32079331
chr4	1348131
```

Now let's generate the reads from fragments:

```bash
# we use bedClip to remove reads outside the chromosome boundary
# we also remove reads mapped to the mitochondrial genome (chrM)

zcat sci-atac/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin dm6/dm6.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > sci-atac/chromap_outs/reads.bed.gz
```

Note we also sort the output reads by `sort -k1,1 -k2,2n`. In this way, the order of chromosomes in the `reads.bed.gz` is the same as that in `dm6.chrom.sizes`, which makes downstream processes easier. The output `reads.bed.gz` are the reads in `bed` format, with the 4th column holding the cell barcodes.

## Peak Calling By MACS2

Now we can use the newly generated read file for the peak calling using `MACS2`:

```console
macs2 callpeak -t sci-atac/chromap_outs/reads.bed.gz \
               -g dm -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir sci-atac/chromap_outs \
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

cut -f 1-4 sci-atac/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > sci-atac/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a sci-atac/chromap_outs/aggregate_peaks_sorted.bed \
    -b sci-atac/chromap_outs/reads.bed.gz \
    -wo -sorted -g dm6/dm6.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > sci-atac/chromap_outs/peak_read_ov.tsv.gz
```

#### Explain Finding Reads In Peaks Per Cell

We start with the command before the first pipe, that is, the intersection part. If you read the manual of the `bedtools intersect`, it should be straightforward to understand. The `-wo` option will output the records in both `-a` file and `-b` file. Since the `reads.bed.gz` file has the cell barcode information at the 4th column, we would get an output with both peak and cell information for the overlap. The `-sorted -g dm6/dm6.chrom.sizes` options make the program use very little memory. Here is an example (top 5 lines) of the output of this part:

```console
chr2L	5477	6081	aggregate_peak_1	chr2L	5430	5480	ATGCGCATTATGCAAGGCCTCTATAGCGTCA	.	-	3
chr2L	5477	6081	aggregate_peak_1	chr2L	5430	5480	ATGCGCTTATTGAGGCTAAGATTACCGTCGG	.	+	3
chr2L	5477	6081	aggregate_peak_1	chr2L	5430	5480	CGCGAACATAGCGACCGCCTCTATACGCAAC	.	-	3
chr2L	5477	6081	aggregate_peak_1	chr2L	5430	5480	CGCGAATTGGCAAGCCTCAGAGCCAGATCCT	.	+	3
chr2L	5477	6081	aggregate_peak_1	chr2L	5430	5480	CGGAGAAGCCGTAGTTAGGCTATAAACTACC	.	+	3
```

We see that the 8th column holds the cell barcode and we want to group them using `bedtools groupby`. Therefore, we need to sort by this column, that is the `sort -k8,8`. When we group by the 8th column, we are interested in how many times each peak appear per group, so we could gather the information of the peak ID (4th column). That is the `-g 8 -c 4 -o freqdesc`. The `-o freqdesc` option returns a `value:frequency` pair in descending order. Here are some records from `peak_read_ov.tsv.gz`:

```console
ATGCGCAACGAATTCGACGTCCTGAACCTAG	aggregate_peak_11321:2,aggregate_peak_12286:2,aggregate_peak_12350:2,aggregate_peak_15932:2,aggregate_peak_15948:2,aggregate_peak_17209:2,aggregate_peak_1747:2,aggregate_peak_17908:2,aggregate_peak_18231:2,aggregate_peak_287:2,aggregate_peak_3680:2,aggregate_peak_381:2,aggregate_peak_4107:2,aggregate_peak_4475:2,aggregate_peak_5726:2,aggregate_peak_9261:2,aggregate_peak_12324:1,aggregate_peak_7571:1
ATGCGCAACGAATTCGACGTCCTGAACGGAC	aggregate_peak_15000:2,aggregate_peak_3846:1
ATGCGCAACGAATTCGACGTCCTGAACTACC	aggregate_peak_7461:2,aggregate_peak_2739:1
```

In a way, that is sort of a count matrix in an awkward format. For example:

- The first line is too long.
- The second line means that in cell `ATGCGCAACGAATTCGACGTCCTGAACGGAC`, the peak `aggregate_peak_15000` has 2 counts, and the peak `aggregate_peak_3846` has 1 count. All the rest peaks not mentioned here have 0 counts in this cell.
- The third line means that in cell `ATGCGCAACGAATTCGACGTCCTGAACTACC`, the peak `aggregate_peak_7461` has 2 counts, and the peak `aggregate_peak_2739` has 1 count. All the rest peaks not mentioned here have 0 counts in this cell.

### Output The Peak-By-Cell Matrix

At this stage, we pretty much have all the things needed. Those two files `aggregate_peaks_sorted.bed` and `peak_read_ov.tsv.gz` contain all information for a peak-by-cell count matrix. We just need a final touch to make the output in a standard format: a [market exchange format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). Since most downstream software takes the input from the __10x Genomics Single Cell ATAC__ results, we are going to generate the MEX and the associated files similar to the output from 10x Genomics.

Here, I'm using a python script for this purpose. You don't have to do this. Choose whatever works for you. The point here is to just generate similar files as the __peak-barcode matrix__ described from [the 10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

First, let's make a directory to hold the output files and generate the `peaks.bed` and `barcodes.tsv` files, which are easy to do:

```bash
# create dirctory
mkdir -p sci-atac/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 sci-atac/chromap_outs/aggregate_peaks_sorted.bed > \
    sci-atac/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat sci-atac/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    sci-atac/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
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
    sci-atac/chromap_outs/aggregate_peaks_sorted.bed \
    sci-atac/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    sci-atac/chromap_outs/peak_read_ov.tsv.gz \
    sci-atac/chromap_outs/raw_peak_bc_matrix
```

After that, you should have the `matrix.mtx` in the `sci-atac/chromap_outs/raw_peak_bc_matrix` directory.

### Cell Calling (Filter Cell Barcodes)

Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, we could use `starolo`. It is used for scRNA-seq in general, but it does have a cell calling function that takes a directory containing raw mtx and associated files, and return the filtered ones. Since `starsolo` looks for the following three files in the input directory: `matrix.mtx`, `features.tsv` and `barcodes.tsv`. Those are the output from the 10x Genomics scRNA-seq workflow. In this case, we can use `peaks.bed` as our `features.tsv`:

```console
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed sci-atac/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     sci-atac/chromap_outs/raw_peak_bc_matrix \
     sci-atac/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv sci-atac/chromap_outs/filtered_peak_bc_matrix/peaks.bed
```

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/sci-atac/
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
    ├── sci-ATAC-C15.csv
    ├── sci-ATAC-D15.csv
    ├── sci-ATAC-P5.csv
    ├── sci-ATAC-P7.csv
    ├── SRR5837698_1.fastq.gz
    ├── SRR5837698_2.fastq.gz
    ├── SRR5837698_CB_base3-33.fastq.gz
    ├── SRR5837698_CB.fastq.gz
    ├── whitelist_config1_base3-33.txt
    ├── whitelist_config1.txt
    └── whitelist_config2.txt

4 directories, 29 files
```