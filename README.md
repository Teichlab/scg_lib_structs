# Single Cell RNA-seq (scRNA-seq) Library Structure
Collections of library structure and sequence of popular single cell genomic methods (mainly scRNA-seq).

## How to use?

Click the following links to veiw the methods. Notes:

1. the default alignment font is Monaco. Courier New font will be used if Monaco is not available.
2. In a dual-index library, how index2 (i5) is sequenced differs from machines to machines. Miseq and Hiseq2000/2500 use the bottom strand as template, which is why the index sequences are the same as the primer sequences in those machines. MiniSeq, NextSeq and Hiseq3000/4000 use the top strand as template, which is why the index sequences are reverse-complementary to the primer sequences in those machines. All methods listed below use Miseq, Hiseq2000/2500 as examples.

- [STRT-seq family (including STRT-seq, STRT-seq-C1 and STRT-seq-2i)](https://teichlab.github.io/scg_lib_structs/methods_html/STRT-seq_family.html)
- [SMART-seq family (including SMART-seq and SMART-seq2)](https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html)
- [Quartz-seq family (including Quartz-seq and Quartz-seq2)](https://teichlab.github.io/scg_lib_structs/methods_html/Quartz-seq_family.html)
- [CEL-seq family (including CEL-seq and CEL-seq2)](https://teichlab.github.io/scg_lib_structs/methods_html/CEL-seq_family.html)
- [10x Chromium Single Cell 3' V3 FeatureBarcoding](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html)
- [10x Chromium Single Cell 3' V2 and V3 GE](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html)
- [10x Chromium Single Cell 3' V1 GE](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html)
- [10x Chromium Single Cell 5' GE](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html)
- [10x Chromium Single Cell ATAC](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html)
- [SureCell 3' WTA for ddSEQ](https://teichlab.github.io/scg_lib_structs/methods_html/SureCell.html)
- [SCRB-seq / mcSCRB-seq](https://teichlab.github.io/scg_lib_structs/methods_html/SCRB-seq.html)
- [Drop-seq / Seq-Well](https://teichlab.github.io/scg_lib_structs/methods_html/Drop-seq.html)
- [Microwell-seq](https://teichlab.github.io/scg_lib_structs/methods_html/Microwell-seq.html)
- [sci-ATAC-seq](https://teichlab.github.io/scg_lib_structs/methods_html/sci-ATAC-seq.html)
- [sci-RNA-seq3](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq3.html)
- [sci-RNA-seq](https://teichlab.github.io/scg_lib_structs/methods_html/sci-RNA-seq.html)
- [Drop-ChIP](https://teichlab.github.io/scg_lib_structs/methods_html/Drop-ChIP.html)
- [MARS-seq](https://teichlab.github.io/scg_lib_structs/methods_html/MARS-seq.html)
- [SPLiT-seq](https://teichlab.github.io/scg_lib_structs/methods_html/SPLiT-seq.html)
- [inDrop](https://teichlab.github.io/scg_lib_structs/methods_html/inDrop.html)

## Technical comparisons (scRNA-seq only)

**The basic chemistry is very similar, the main differences among those scRNA-seq methods are summarised in the table below. For a detailed discussion, check the text boxes from our review: [From Tissues to Cell Types and Back: Single-Cell Gene Expression Analysis of Tissue Architecture](https://www.annualreviews.org/doi/10.1146/annurev-biodatasci-080917-013452)**

|                                  | Single cell isolation/capture |        2nd strand synthesis       | Full-length cDNA synthesis |                      Barcode addition                     | Pooling before library |  Library amplification | Gene coverage |
|:--------------------------------:|:-----------------------------:|:---------------------------------:|:--------------------------:|:---------------------------------------------------------:|:----------------------:|:----------------------:|:-------------:|
|    10x Chromium Single Cell 3'   |            Droplet            |                TSO                |             Yes            |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
|    10x Chromium Single Cell 5'   |            Droplet            |                TSO                |             Yes            |                    Barcoded TSO primers                   |           Yes          |           PCR          |       5'      |
|         CEL-seq/CEL-seq2         |              FACS             |       RNase H and DNA pol I       |             No             |                    Barcoded RT primers                    |           Yes          | In vitro transcription |       3'      |
|             Drop-seq             |            Droplet            |                TSO                |             Yes            |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
|             MARS-seq             |              FACS             |       RNase H and DNA pol I       |             No             |                    Barcoded RT primers                    |           Yes          | In vitro transcription |       3'      |
|           Microwell-seq          |           Nanowells           |                TSO                |             Yes            |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
|            Quartz-seq            |              FACS             | PolyA tailing and primer ligation |             Yes            |            Ligation of barcoded Truseq adapters           |           No           |           PCR          |       3'      |
|            Quartz-seq2           |              FACS             | PolyA tailing and primer ligation |             Yes            |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
|       SMART-seq/<br>SMART-seq2       |        FACS or Fluidigm C1       |                TSO                |             Yes            |             Library PCR with barcoded primers             |           No           |           PCR          |  full-length  |
|             SPLiT-seq            |           Not needed          |                TSO                |             Yes            |              Ligation of barcoded RT primers              |           Yes          |           PCR          |       3'      |
|             STRT-seq             |              FACS             |                TSO                |             Yes            |                    Barcoded TSO primers                   |           Yes          |           PCR          |       5'      |
|            STRT-seq-C1           |          Fluidigm C1          |                TSO                |             Yes            |                  Barcoded Tn5 transposase                 |           No           |           PCR          |       5'      |
|            STRT-seq-2i           |        FACS or dilution       |                TSO                |             Yes            |         Barcoded  PCR primers and Tn5 transposase         |           Yes          |           PCR          |       5'      |
|             Seq-Well             |           Nanowells           |                TSO                |             Yes            |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
| Illumina Bio-Rad SureCell 3' WTA |            Droplet            |       RNase H and DNA pol I       |             No             |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
|              inDrop              |            Droplet            |       RNase H and DNA pol I       |             No             |                    Barcoded RT primers                    |           Yes          | In vitro transcription |       3'      |
|        SCRB-seq/<br>mcSCRB-seq       |              FACS             |                TSO                |             Yes            |                    Barcoded RT primers                    |           Yes          |           PCR          |       3'      |
|            sci-RNA-seq           |           Not needed          |       RNase H and DNA pol I       |             No             | Barcoded RT primers and library PCR with barcoded primers |           Yes          |           PCR          |       3'      |

## Motivation

I was a little bit bombarded with all the single cell methods and got completely lost. To help myself understand all of them and future troubleshooting, I start to perform an on-paper library preparation whenever I see a new single cell method.

## Why bother?

Here I borrow from Feyman:

**What I cannot create, I do not understand.**

----

![](data/feyman.jpeg)

## TODO:

- Waiting for new methods ...

## Feedback

I would be very happy if you go through them and let me know what you think. If you spot some errors/mistakes, or I've missed some key methods. Feel free to contact me:

Xi Chen  
xc1@sanger.ac.uk
