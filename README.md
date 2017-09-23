# Single Cell RNA-seq (scRNA-seq) Library Structure
Collections of library structure and sequence of popular single cell genomic methods (mainly scRNA-seq).

## Motivation

I was a little bit bombarded with all the single cell methods and got completely lost. To help with myself understand all of them and future troubleshooting, I start to perform an on-paper library preparation whenever I see a new single cell method.

## Why bother?

Here I borrow from Feyman:

**What I cannot create on paper, I do not understand. Know how to re-construct every library that has been invented.**

----

![](data/feyman.jpeg)

## How to use?

Click the following links to veiw the methods. Notes:

1. the default alignment font is Monaco. Courier New font will be used if Monaco is not available.
2. In a dual-index library, how index2 (i5) is sequenced differs from machines to machines. Miseq and Hiseq2000/2500 use the bottom strand as template, which is why the index sequences are the same as the primer sequences in those machines. MiniSeq, NextSeq and Hiseq3000/4000 use the top strand as template, which is why the index sequences are reverse-complementary to the primer sequences in those machines. All methods listed below use Miseq, Hiseq/2000/2500 as examples.

- [STRT-seq family (including STRT-seq, STRT-seq-C1, STRT-seq-2i)](https://teichlab.github.io/scg_lib_structs/STRT-seq_family.html)
- [SMART-seq family (including SMART-seq, SMART-seq2)](https://teichlab.github.io/scg_lib_structs/SMART-seq_family.html)
- [CEL-seq family (including CEL-seq, CEL-seq2)](https://teichlab.github.io/scg_lib_structs/CEL-seq_family.html)
- [10x Chromium Single Cell 3' Solution v2](https://teichlab.github.io/scg_lib_structs/10xChromium.html)
- [Drop-seq / Seq-Well](https://teichlab.github.io/scg_lib_structs/Drop-seq.html)
- [inDrop](https://teichlab.github.io/scg_lib_structs/inDrop.html)
- [MARS-seq](https://teichlab.github.io/scg_lib_structs/MARS-seq.html)
- [sci-RNA-seq](https://teichlab.github.io/scg_lib_structs/sci-RNA-seq.html)
- [SPLiT-seq](https://teichlab.github.io/scg_lib_structs/SPLiT-seq.html)

## TODO:

- Quartz-seq/Quartz-seq2

## Feedback

I would be very happy if you go through them and let me know what you think. If you spot some errors/mistakes, or I've missed some key methods. Feel free to contact me:

Xi Chen  
xc1@sanger.ac.uk
