# 10x Chromium Single Cell Multiome ATAC + Gene Expression

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html) to see how __10x Chromium Single Cell Multiome ATAC + Gene Expression__ libraries are generated experimentally. This is a droplet-based method, where transposed nuclei are captured inside droplets. At the same time, gel beads with two types of barcoded oligos are also captured inside the droplet. One type of oligo contains oligo-dT to capture mRNA, and the other type of oligo contains an 8-bp linker sequence to capture ATAC fragments. See [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html) for more details. After droplet capture, reverse transcription and ATAC fragments capture are performed. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead. Then, all droplets from one sample is collected. After that, seven cycles of pre-amplification are performed. Then the reaction is split into two portions: one for gene expression library preparation and the other for ATAC library preparation.

## For Your Own Experiments

If you follow the user manual from 10x Genomics, you should have two libraries per sample: one for gene expression and the other for ATAC. Normally, you prepare them independently and sequence them separately. See below the sequencing read configurations.

### Gene Expression Read Configuration

Regardless of the sequencing machines, you always use the following configuration:

| Order | Read             | Cycle | Description                                           |
|-------|------------------|-------|-------------------------------------------------------|
| 1     | Read 1           | 28    | `R1_001.fastq.gz`, 16 bp cell barcodes + 12 bp UMI    |
| 2     | Index 1 (__i7__) | 10    | `I1_001.fastq.gz`, Sample index                       |
| 3     | Index 2 (__i5__) | 10    | `I2_001.fastq.gz`, Sample index (if using dual index) |
| 4     | Read 2           | >50   | `R2_001.fastq.gz`, cDNA reads                         |

If you sequence your data via your core facility or a company, you will need to provide the sample index sequence, which is the primer (__PN-1000215__) taken from the commercial kit from 10x Genomics (dual index in this case), to them and they will demultiplex for you. Tell them it is 10x Multiome Gene Expression library and they will know what to do.

If you sequence by yourself, you need to run `bcl2fastq` by yourself with a `SampleSheet.csv`. Here is an example of `SampleSheet.csv` of a NextSeq run with two different samples using the indexing primers from the A1 and B1 wells, respectively:

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
Sample01_GEX,,,,,,SI-TT-A1_i7,GTAACATGCG,SI-TT-A1_i5,AGGTAACACT,,
Sample02_GEX,,,,,,SI-TT-B1_i7,ACAGTAACTA,SI-TT-B1_i5,AACGAACTGT,,
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

After this, you will have `R1_001.fastq.gz` and `R2_001.fastq.gz` for each sample:

```bash
Sample01_GEX_S1_R1_001.fastq.gz # 28 bp: cell barcodes (16 bp) + UMI (12 bp)
Sample01_GEX_S1_R2_001.fastq.gz # >50 bp: cDNA reads
Sample02_GEX_S2_R1_001.fastq.gz # 28 bp: cell barcodes (16 bp) + UMI (12 bp)
Sample02_GEX_S2_R2_001.fastq.gz # >50 bp: cDNA reads
```

You are ready to go from here.

### ATAC Read Configuration

The ATAC library is slightly more complicated than the gene expression library. Make sure you understand how sequencing is done for this assay by checking [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html). Then it is relatively straightforward to see there are two configurations.

Using more recent machines and chemistries, like __iSeq 100__, __MiniSeq (Standard)__, __NextSeq__, __HiSeq X__, __HiSeq 3000__, __HiSeq 4000__ and __NovaSeq 600 (v1.5)__, it should be (I call this __Configuration 1__):

| Order | Read             | Cycle                     | Description                                                   |
|-------|------------------|---------------------------|---------------------------------------------------------------|
| 1     | Read 1           | >50                       | This normally yields `R1_001.fastq.gz`, Genomic insert        |
| 2     | Index 1 (__i7__) | 8                         | This normally yields `I1_001.fastq.gz`, Sample index          |
| 3     | Index 2 (__i5__) | 8 dark + 16 normal cycles | This normally yields `I2_001.fastq.gz`, Cell barcodes (16 bp) |
| 4     | Read 2           | >50                       | This normally yields `R2_001.fastq.gz`, Genomic insert        |

Using older machines and chemistries like __MiSeq__, __HiSeq 2000__, __HiSeq 2500__, __MiniSeq (Rapid)__ and __NovaSeq 6000 (v1.0)__, it should be (I call this __Configuration 2__):

| Order | Read             | Cycle | Description                                                   |
|-------|------------------|-------|---------------------------------------------------------------|
| 1     | Read 1           | >50   | This normally yields `R1_001.fastq.gz`, Genomic insert        |
| 2     | Index 1 (__i7__) | 8     | This normally yields `I1_001.fastq.gz`, Sample index          |
| 3     | Index 2 (__i5__) | 16    | This normally yields `I2_001.fastq.gz`, Cell barcodes (16 bp) |
| 4     | Read 2           | >50   | This normally yields `R2_001.fastq.gz`, Genomic insert        |

If you sequence your data via your core facility or a company, you will need to provide the sample index sequence, which is the primer (__PN-1000212__) taken from the commercial kit from 10x Genomics, to them and they will demultiplex for you. Tell them it is 10x Multiome ATAC library and they will know what to do.

If you sequence by yourself, you need to run `bcl2fastq` with a `SampleSheet.csv` in specific way. Here is an example of `SampleSheet.csv` of a NextSeq run with two different samples using the indexing primers from the A1 and B1 wells, respectively:

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
Sample01_ATAC,,,,,,SI-NA-A1_1,AAACGGCG,,,,
Sample01_ATAC,,,,,,SI-NA-A1_2,CCTACCAT,,,,
Sample01_ATAC,,,,,,SI-NA-A1_3,GGCGTTTC,,,,
Sample01_ATAC,,,,,,SI-NA-A1_4,TTGTAAGA,,,,
Sample02_ATAC,,,,,,SI-NA-B1_1,AGGCTACC,,,,
Sample02_ATAC,,,,,,SI-NA-B1_2,CTAGCTGT,,,,
Sample02_ATAC,,,,,,SI-NA-B1_3,GCCAACAA,,,,
Sample02_ATAC,,,,,,SI-NA-B1_4,TATTGGTG,,,,
```

You can see each sample actually has four different index sequences. This is because each well from the plate __PN-1000212__ actually contain four different indices for base balancing.

Now let's look at the order of the sequencing read configuration above, as you can see, the first (`R1`), the 3rd (`I2`) and the 4th (`R2`) reads are all important for us. Therefore, we would like to get all of them for each sample based on sample index, that is, the 2nd read (`I1`). To do this, you should run `bcl2fastq` in the following way:

```console
bcl2fastq --use-bases-mask=Y50,I8,Y16,Y50 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but I the important bit that needs explanation is `--use-bases-mask=Y50,I8,Y16,Y50`. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y50` at the first position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `I8` at the second position indicates "use the cycle as an index read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read.
3. `Y16` at the third position indicates "use the cycle as a real read", so you will get 16-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read, though it is the 3rd read overall.
4. `Y50` at the fourth position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

Therefore, you will get four fastq file per sample. Using the examples above, these are the files you should get:

```bash
# files for Sample01

Sample01_ATAC_S1_I1_001.fastq.gz # 8 bp: sample index
Sample01_ATAC_S1_R1_001.fastq.gz # 50 bp: genomic insert
Sample01_ATAC_S1_R2_001.fastq.gz # 16 bp: cell barcodes
Sample01_ATAC_S1_R3_001.fastq.gz # 50 bp: genomic insert 

# files for Sample02

Sample02_ATAC_S2_I1_001.fastq.gz # 8 bp: sample index
Sample02_ATAC_S2_R1_001.fastq.gz # 50 bp: genomic insert
Sample02_ATAC_S2_R2_001.fastq.gz # 16 bp: cell barcodes
Sample02_ATAC_S2_R3_001.fastq.gz # 50 bp: genomic insert
```

We can safely ignore the `I1` files, but the naming here is really different from our normal usage. The `R1` files are good. The `R2` files here actually means `I2` in our normal usage. The `R3` files here actually means `R2` in our normal usage. Anyway, __DO NOT get confused__. You are ready to go from here.

## Public Data

We will use the example data set from the 10x website, using the [Fresh Embryonic E18 Mouse Brain (5k)](https://www.10xgenomics.com/resources/datasets/fresh-embryonic-e-18-mouse-brain-5-k-1-standard-2-0-0) data.

```console
mkdir 10xMultiome
wget -P 10xMultiome -c https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/e18_mouse_brain_fresh_5k/e18_mouse_brain_fresh_5k_fastqs.tar
tar xf 10xMultiome/e18_mouse_brain_fresh_5k_fastqs.tar -C 10xMultiome/
```

These are the files after extraction:

```console
scg_prep_test/10xMultiome/
├── e18_mouse_brain_fresh_5k
│   ├── atac
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_I1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_R3_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L002_I1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz
│   │   └── e18_mouse_brain_fresh_5k_S1_L002_R3_001.fastq.gz
│   └── gex
│       ├── e18_mouse_brain_fresh_5k_S1_L001_I1_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L001_I2_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L002_I1_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L002_I2_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz
│       └── e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz
└── e18_mouse_brain_fresh_5k_fastqs.tar
```
Like said before, we could safely ignore all the `I1` reads. The data was also split based on lanes. For gene expression data, the `I2` reads can also be ignored.

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the __10x Chromium Single Cell Multiome ATAC + Gene Expression__ kit. If you have `cellranger-arc` in your computer, you will find two files, both of which are called `737K-arc-v1.txt.gz`. One is in the `lib/python/cellranger/barcodes` directory, which is the gene expression whitelist. The other is in the `lib/python/atac/barcodes` directory, which is the ATAC whitelist. If you don't have `cellranger-atac`, I have prepared and renamed the file to be a more informative for you:

```console
# download the GEX expression whitelist

wget -P 10xMultiome/ https://teichlab.github.io/scg_lib_structs/data/gex_737K-arc-v1.txt.gz
gunzip 10xMultiome/gex_737K-arc-v1.txt.gz

# download the atac whitelist

wget -P 10xMultiome/ https://teichlab.github.io/scg_lib_structs/data/atac_737K-arc-v1.txt.gz
gunzip 10xMultiome/atac_737K-arc-v1.txt.gz

# reverse complement the atac whitelist

cat 10xMultiome/atac_737K-arc-v1.txt | \
    rev | tr 'ACGT' 'TGCA' > \
    10xMultiome/atac_737K-arc-v1_rc.txt
```

### Explain Whitelist

The lengths of `gex_737K-arc-v1.txt` and `atac_737K-arc-v1.txt` are the same, both of which contain 736,320 barcodes. They are the barcodes in the capture oligos on the beads, one for mRNA and the other for ATAC. You can pair them to make a barcode translation file during the downstream analysis. The 1st barcode in the `gex_737K-arc-v1.txt` and the 1st barcode from the `atac_737K-arc-v1.txt` are from the same gel bead; The 2nd barcode in the `gex_737K-arc-v1.txt` and the 2nd barcode from the `atac_737K-arc-v1.txt` are from the same gel bead; ... The 736,320th barcode in the `gex_737K-arc-v1.txt` and the 736,320th barcode from the `atac_737K-arc-v1.txt` are from the same gel bead. Using this information, you are able to know what GEX reads and ATAC reads are from the same cell.

You may wonder what is the reverse complementary step about. The cell barcodes in the `atac_737K-arc-v1_rc.txt` are the sequences on the gel beads. These 16 bp cell barcodes are in the `i5` index location, that is, between Illumina P5 and the Nextera Read 1 sequence. It means they will be sequenced as `Index 2` (`I2`). How `i5` or `Index 2` is sequenced depends on the machine. Previously, __MiSeq__, __HiSeq 2000__, __HiSeq 2500__, __MiniSeq (Rapid)__ and __NovaSeq 6000 (v1.0)__ use the bottom strand as the template, so the index reads will be the same as the barcodes in the `737K-cratac-v1.txt`. However, more recent machines and chemistries, like __iSeq 100__, __MiniSeq (Standard)__, __NextSeq__, __HiSeq X__, __HiSeq 3000__, __HiSeq 4000__ and __NovaSeq 600 (v1.5)__, use the top strand as the template, so the index reads will be reverse complementary to the barcodes in the `atac_737K-arc-v1_rc.txt`. Therefore, we need to create a reverse complementary file as the whitelist for some data.

```{eval-rst}
.. tip::
  
  Whenever you are dealing with reads that come from ``index 2``, you should:

  1. Check the sequencing machine used to generate the data
  2. Make sure you are familiar with different sequencing modes from different Illumina machines by looking at `this page <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_.
  3. Extract some sequences from ``index 2``, compare them to the whitelist or the reverse complementary to the whitelist. 
```

## From FastQ To Count Matrices

Now we are ready to map the reads to the genome using `starsolo` for the RNA library and `chromap` for the ATAC library:

```console
mkdir -p 10xMultiome/star_outs
mkdir -p 10xMultiome/chromap_outs

# process the GEX library using starsolo

STAR --runThreadN 4 \
     --genomeDir mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix 10xMultiome/star_outs/ \
     --readFilesIn 10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz 10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/gex/e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 \
     --soloCBwhitelist 10xMultiome/gex_737K-arc-v1.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate

# process the ATAC library using chromap

## map and generate the fragment file

chromap -t 4 --preset atac \
        -x mm10/chromap_index/genome.index \
        -r mm10/mm10.fa \
        -1 10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz \
        -2 10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L001_R3_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L002_R3_001.fastq.gz \
        -b 10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz,10xMultiome/e18_mouse_brain_fresh_5k/atac/e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz \
        --barcode-whitelist 10xMultiome/atac_737K-arc-v1.txt \
        -o 10xMultiome/chromap_outs/fragments.tsv

## compress and index the fragment file

bgzip 10xMultiome/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed 10xMultiome/chromap_outs/fragments.tsv.gz
```

After this stage, we are done with the GEX library. The count matrix and other useful information can be found in the `star_outs` directory. For the ATAC library, two new files `fragments.tsv.gz` and `fragments.tsv.gz.tbi` are generated. They will be useful and sometimes required for other programs to perform downstream analysis. There are still some extra work.

### Explain star and chromap

If you understand the __10x Chromium Single Cell Multiome ATAC + Gene Expression__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html), the commands above should be straightforward to understand.

#### Explain star

`--runThreadN 4`
  
>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mm10/star_index`

>> Pointing to the directory of the star index. The public data we are analysing is from mouse brains.

`--readFilesCommand zcat`

>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix 10xMultiome/star_outs/`

>> We want to keep everything organised. This directs all output files inside the `10xMultiome/star_outs` directory.

`--readFilesIn`

>> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. In __10x Chromium Single Cell Multiome ATAC + Gene Expression__, cDNA reads come from Read 2, and the cell barcode and UMI come from Read 1. Check [the 10x Chromium Single Cell Multiome ATAC + Gene Expression GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html) if you are not sure. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--soloType CB_UMI_Simple`

>> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter.

`--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12`

>> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist 10xMultiome/gex_737K-arc-v1.txt`

>> The plain text file containing all possible valid cell barcodes, one per line. __10x Chromium Single Cell Multiome ATAC + Gene Expression__ is a commercial platform. The whitelist is taken from their commercial software `cellranger-arc`.

`--soloCellFilter EmptyDrops_CR`

>> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter.

`--soloStrand Forward`

>> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __10x Chromium Single Cell Multiome ATAC + Gene Expression__, the cDNA reads are from the Read 2 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing UMI and the Illumina Read 1 sequence. Therefore, Read 1 consists of cell barcodes and UMI comes from the first strand, complementary to the coding strand. Read 2 comes from the coding strand. Therefore, use `Forward` for __10x Chromium Single Cell Multiome ATAC + Gene Expression__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear. Check [the 10x Chromium Single Cell Multiome ATAC + Gene Expression GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html) if you are not sure.

`--outSAMattributes CB UB`

>> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

>> We want sorted `BAM` for easy handling by other programs.

#### Explain chromap

`-t 4`

>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`-x mm10/chromap_index/genome.index`

>> The `chromap` index file. The public data We are analysing is from mouse brains.

`-r mm10/mm10.fa`

>> Reference genome sequence in `fasta` format. This is basically the file which you used to create the `chromap` index file.

`-1`, `-2` and `-b`

>> They are Read 1 (genomic), Read 2 (genomic) and cell barcode read, respectively. For ATAC-seq, the sequencing is usually done in pair-end mode. Therefore, you normally have two genomic reads for each genomic fragment: Read 1 and Read 2. For the reason described previously, `R1` is the genomic Read 1 and should be passed to `-1`; `R3` is actually the genomic Read 2 and should be passed to `-2`; `R2` is the cell barcode read and should be passed to `-b`. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--barcode-whitelist 10xMultiome/atac_737K-arc-v1.txt`

>> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell ATAC__ is a commercial platform. The whitelist is taken from their commercial software `cellranger-arc`. In this example data, sequencing is done using __NovaSeq 6000 (v1.0)__. Therefore, we use the original whitelist. In other cases, you might want to use the reverse complementary version of the whitelist.

`-o 10xMultiome/chromap_outs/fragments.tsv`

>> Direct the mapped fragments to a file. The format is described in the [10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments).

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

faSize -detailed mm10/mm10.fa | \
    sort -k1,1 > mm10/mm10.chrom.sizes
```

This is the first 5 lines of `mm10/mm10.chrom.sizes`:

```
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

zcat 10xMultiome/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin mm10/mm10.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > 10xMultiome/chromap_outs/reads.bed.gz
```

Note we also sort the output reads by `sort -k1,1 -k2,2n`. In this way, the order of chromosomes in the `reads.bed.gz` is the same as that in `mm10.chrom.sizes`, which makes downstream processes easier. The output `reads.bed.gz` are the reads in `bed` format, with the 4th column holding the cell barcodes.

### Peak Calling By MACS2

Now we can use the newly generated read file for the peak calling using `MACS2`:

```console
macs2 callpeak -t 10xMultiome/chromap_outs/reads.bed.gz \
               -g mm -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir 10xMultiome/chromap_outs \
               -n aggregate
```

#### Explain MACS2

The reasons of choosing those specific parameters are a bit more complicated. I have dedicated a post for this a while ago. Please have a look at [__this post__](https://dbrg77.github.io/posts/2020-12-09-atac-seq-peak-calling-with-macs2/) if you are still confused. The following output files are particularly useful:

| File                       | Description                                                               |
|----------------------------|---------------------------------------------------------------------------|
| aggregate_peaks.narrowPeak | Open chromatin peak locations in the narrowPeak format                    |
| aggregate_peaks.xls        | More information about peaks                                              |
| aggregate_treat_pileup.bdg | Signal tracks. Can be used to generate the bigWig file for visualisation |

### Getting The Peak-By-Cell Count Matrix

Now that we have the peak and reads files, we can compute the number of reads in each peak for each cell. Then we could get the peak-by-cell count matrix. There are different ways of doing this. The following is the method I use.

#### Find Reads In Peaks Per Cell

First, we use the `aggregate_peaks.narrowPeak` file. We only need the first 4 columns (chromosome, start, end, peak ID). You can also remove the peaks that overlap [the black list regions](https://www.nature.com/articles/s41598-019-45839-z). The black list is not available for every species and every build, so I'm not doing it here. We also need to sort the peak to make sure the order of the chromosomes in the peak file is the same as that in the `mm10.chrom.sizes` and `reads.bed.gz` files. Then we could find the overlap by `bedtools`. We need to do this in a specific way to get the number of reads in each peak from each cell:

```bash
# format and sort peaks

cut -f 1-4 10xMultiome/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > 10xMultiome/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a 10xMultiome/chromap_outs/aggregate_peaks_sorted.bed \
    -b 10xMultiome/chromap_outs/reads.bed.gz \
    -wo -sorted -g mm10/mm10.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > 10xMultiome/chromap_outs/peak_read_ov.tsv.gz
```

##### Explain Finding Reads In Peaks Per Cell

We start with the command before the first pipe, that is, the intersection part. If you read the manual of the `bedtools intersect`, it should be straightforward to understand. The `-wo` option will output the records in both `-a` file and `-b` file. Since the `reads.bed.gz` file has the cell barcode information at the 4th column, we would get an output with both peak and cell information for the overlap. The `-sorted -g mm10/mm10.chrom.sizes` options make the program use very little memory. Here is an example (top 5 lines) of the output of this part:

```console
chr1	3060888	3061088	aggregate_peak_75	chr1	3060705	3060905	GGATTAGGTTAATCCA	.	+	17
chr1	3060888	3061088	aggregate_peak_75	chr1	3060728	3060928	GAATGACTCGGTCCAT	.	+	40
chr1	3060888	3061088	aggregate_peak_75	chr1	3060759	3060959	ATTCCGGGTACCGTTG	.	+	71
chr1	3060888	3061088	aggregate_peak_75	chr1	3060777	3060977	GCAATACAGTTACTAC	.	+	89
chr1	3060888	3061088	aggregate_peak_75	chr1	3060786	3060986	AATCATGAGTCATCAT	.	+	98
```

We see that the 8th column holds the cell barcode and we want to group them using `bedtools groupby`. Therefore, we need to sort by this column, that is the `sort -k8,8`. When we group by the 8th column, we are interested in how many times each peak appear per group, so we could gather the information of the peak ID (4th column). That is the `-g 8 -c 4 -o freqdesc`. The `-o freqdesc` option returns a `value:frequency` pair in descending order. Here are some records from `peak_read_ov.tsv.gz`:

```console
AAACAAGCAAAGAAGC	aggregate_peak_117806:2,aggregate_peak_35843:2,aggregate_peak_39908:2,aggregate_peak_42315:2
AAACAAGCAAAGGAAC	aggregate_peak_68626:2
AAACAAGCAAAGGCCT	aggregate_peak_116364:2,aggregate_peak_17597:2,aggregate_peak_48133:2
```

In a way, that is sort of a count matrix in an awkward format. For example:

- The first line means that in cell `AAACAAGCAAAGAAGC`, the peak `aggregate_peak_117806` has 2 counts, the peak `aggregate_peak_35843` has 2 counts, the peak `aggregate_peak_39908` has 2 counts and the peak `aggregate_peak_42315` has 2 counts. All the rest peaks not mentioned here have 0 counts in this cell.
- The second line means that in cell `AAACAAGCAAAGGAAC`, the peak `aggregate_peak_68626` has 2 counts. All the rest peaks not mentioned here have 0 counts in this cell.

#### Output The Peak-By-Cell Matrix

At this stage, we pretty much have all the things needed. Those two files `aggregate_peaks_sorted.bed` and `peak_read_ov.tsv.gz` contain all information for a peak-by-cell count matrix. We just need a final touch to make the output in a standard format: a [market exchange format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). Since most downstream software takes the input from the __10x Genomics Single Cell ATAC__ results, we are going to generate the MEX and the associated files similar to the output from 10x Genomics.

Here, I'm using a python script for this purpose. You don't have to do this. Choose whatever works for you. The point here is to just generate similar files as the __peak-barcode matrix__ described from [the 10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

First, let's make a directory to hold the output files and generate the `peaks.bed` and `barcodes.tsv` files, which are easy to do:

```bash
# create dirctory
mkdir -p 10xMultiome/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 10xMultiome/chromap_outs/aggregate_peaks_sorted.bed > \
    10xMultiome/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat 10xMultiome/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    10xMultiome/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
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
    10xMultiome/chromap_outs/aggregate_peaks_sorted.bed \
    10xMultiome/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    10xMultiome/chromap_outs/peak_read_ov.tsv.gz \
    10xMultiome/chromap_outs/raw_peak_bc_matrix
```

After that, you should have the `matrix.mtx` in the `10xMultiome/chromap_outs/raw_peak_bc_matrix` directory.

#### Cell Calling (Filter Cell Barcodes)

Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, we could use `starolo`. It is used for scRNA-seq in general, but it does have a cell calling function that takes a directory containing raw mtx and associated files, and return the filtered ones. Since `starsolo` looks for the following three files in the input directory: `matrix.mtx`, `features.tsv` and `barcodes.tsv`. Those are the output from the 10x Genomics scRNA-seq workflow. In this case, we can use `peaks.bed` as our `features.tsv`:

```console
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed 10xMultiome/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     10xMultiome/chromap_outs/raw_peak_bc_matrix \
     10xMultiome/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv 10xMultiome/chromap_outs/filtered_peak_bc_matrix/peaks.bed
```

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/10xMultiome/
├── atac_737K-arc-v1_rc.txt
├── atac_737K-arc-v1.txt
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
├── e18_mouse_brain_fresh_5k
│   ├── atac
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_I1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L001_R3_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L002_I1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz
│   │   ├── e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz
│   │   └── e18_mouse_brain_fresh_5k_S1_L002_R3_001.fastq.gz
│   └── gex
│       ├── e18_mouse_brain_fresh_5k_S1_L001_I1_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L001_I2_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L001_R1_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L001_R2_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L002_I1_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L002_I2_001.fastq.gz
│       ├── e18_mouse_brain_fresh_5k_S1_L002_R1_001.fastq.gz
│       └── e18_mouse_brain_fresh_5k_S1_L002_R2_001.fastq.gz
├── e18_mouse_brain_fresh_5k_fastqs.tar
├── gex_737K-arc-v1.txt
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

11 directories, 55 files
```