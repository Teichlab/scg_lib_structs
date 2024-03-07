# 10x Genomics Single Cell ATAC

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html) to see how __10x Chromium Single Cell ATAC__ libraries are generated experimentally. This is a droplet-based method, where transposed nuclei are captured inside droplets. At the same time, gel beads with barcoded primer containing Nextera Read 1 sequence are also captured inside the droplet. A few cycles of linear PCR (single primer PCR) are performed inside the droplet. Since the transposed nuclei have the Nextera Read 1 sequence in the open chromatin, the DNA fragment will be labelled with the barcodes during the linear PCR stage. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead. Then, all droplets from one sample is collected and a sample index is added during the library amplification stage.

## For Your Own Experiments

If you sequence your data via your core facility or a company, you will need to provide the sample index sequence, which is the primer (__PN-1000212/PN-1000084__) taken from the commercial kit from 10x Genomics, to them and they will demultiplex for you. Tell them it is 10x scATAC library and they will know what to do.

If you sequence by yourself, you need to run `bcl2fastq` by yourself with a `SampleSheet.csv` in a very specific way. Here is an example of `SampleSheet.csv` of a NextSeq run with two different samples using the indexing primers from the A1 and B1 wells, respectively:

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
Sample01,,,,,,SI-GA-A1_1,AAACGGCG,,,,
Sample01,,,,,,SI-GA-A1_2,CCTACCAT,,,,
Sample01,,,,,,SI-GA-A1_3,GGCGTTTC,,,,
Sample01,,,,,,SI-GA-A1_4,TTGTAAGA,,,,
Sample02,,,,,,SI-GA-B1_1,AGGCTACC,,,,
Sample02,,,,,,SI-GA-B1_2,CTAGCTGT,,,,
Sample02,,,,,,SI-GA-B1_3,GCCAACAA,,,,
Sample02,,,,,,SI-GA-B1_4,TATTGGTG,,,,
```

You can see each sample actually has four different index sequences. This is because each well from the plate __PN-1000212/PN-1000084__ actually contain four different indices for base balancing.

```{eval-rst}
.. important::
  
  Make sure you understand how sequencing is done for this assay by checking `this GitHub page <https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html>`_. There are a total of four reads in this order:

  =====  ================  =======  ========================================================
  Order  Read              Cycle    Description
  =====  ================  =======  ========================================================
    1    Read 1            > 50     This normally yields ``R1_001.fastq.gz``, Genomic insert
    2    Index 1 (**i7**)  8 or 10  This normally yields ``I1_001.fastq.gz``, Sample index
    3    Index 2 (**i5**)  16       This normally yields ``I2_001.fastq.gz``, Cell barcodes
    4    Read 2            > 50     This normally yields ``R2_001.fastq.gz``, Genomic insert
  =====  ================  =======  ========================================================
```

Look at the order of the sequencing, as you can see, the first (`R1`), the 3rd (`I2`) and the 4th (`R2`) reads are all important for us. Therefore, we would like to get all of them for each sample based on sample index, that is, the 2nd read (`I1`). To do this, you should run `bcl2fastq` in the following way:

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

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is `--use-bases-mask=Y50,I8,Y16,Y50`. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y50` at the first position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `I8` at the second position indicates "use the cycle as an index read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read.
3. `Y16` at the third position indicates "use the cycle as a real read", so you will get 16-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read, though it is the 3rd read overall.
4. `Y50` at the fourth position indicates "use the cycle as a real read", so you will get 50-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

Therefore, you will get four fastq file per sample. Using the examples above, these are the files you should get:

```bash
# files for Sample01

Sample01_S1_I1_001.fastq.gz
Sample01_S1_R1_001.fastq.gz
Sample01_S1_R2_001.fastq.gz
Sample01_S1_R3_001.fastq.gz

# files for Sample02

Sample02_S2_I1_001.fastq.gz
Sample02_S2_R1_001.fastq.gz
Sample02_S2_R2_001.fastq.gz
Sample02_S2_R3_001.fastq.gz
```

We can safely ignore the `I1` files, but the naming here is really different from our normal usage. The `R1` files are good. The `R2` files here actually mean `I2` in our normal usage. The `R3` files here actually mean `R2` in our normal usage. Anyway, __DO NOT get confused__.

## Public Data

We will use the example data set from the 10x website, using the [500 Peripheral Blood Mononuclear Cells (PBMCs) from a Healthy Donor (Next GEM v1.1)](https://www.10xgenomics.com/resources/datasets/500-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0) data.

```console
mkdir -p 10xscATAC/pbmc500
wget -P 10xscATAC/pbmc500 -c https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fastqs.tar
tar xf 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs.tar -C 10xscATAC/pbmc500/
```

These are the files after extraction:

```console
scg_prep_test/10xscATAC/
└── pbmc500
    ├── atac_pbmc_500_nextgem_fastqs
    │   ├── atac_pbmc_500_nextgem_S1_L001_I1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L001_R1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L001_R2_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L001_R3_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L002_I1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L002_R1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L002_R2_001.fastq.gz
    │   └── atac_pbmc_500_nextgem_S1_L002_R3_001.fastq.gz
    └── atac_pbmc_500_nextgem_fastqs.tar
```
Like said before, we could safely ignore all the `I1` reads.

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the Single Cell ATAC kit. If you have `cellranger-atac` in your computer, you will find a file called `737K-cratac-v1.txt.gz` in the `lib/python/atac/barcodes` directory. If you don't have `cellranger-atac`, I have prepared the file for you:

```console
# download the whitelist

wget -P 10xscATAC/pbmc500/ https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-cratac-v1.txt.gz
gunzip 10xscATAC/pbmc500/737K-cratac-v1.txt.gz

# reverse complement the whitelist

cat 10xscATAC/pbmc500/737K-cratac-v1.txt | \
    rev | tr 'ACGT' 'TGCA' > \
    10xscATAC/pbmc500/737K-cratac-v1_rc.txt
```

### Explain Whitelist

You may wonder what is the reverse complementary step about. The cell barcodes in the `737K-cratac-v1.txt` are the sequences on the gel beads. These 16 bp cell barcodes are in the `i5` index location, that is, between Illumina P5 and the Nextera Read 1 sequence. It means they will be sequenced as `Index 2` (`I2`). How `i5` or `Index 2` is sequenced depends on the machine. Previously, __MiSeq__, __HiSeq 2000__, __HiSeq 2500__, __MiniSeq (Rapid)__ and __NovaSeq 6000 (v1.0)__ use the bottom strand as the template, so the index reads will be the same as the barcodes in the `737K-cratac-v1.txt`. However, more recent machines and chemistries, like __iSeq 100__, __MiniSeq (Standard)__, __NextSeq__, __HiSeq X__, __HiSeq 3000__, __HiSeq 4000__ and __NovaSeq 600 (v1.5)__, use the top strand as the template, so the index reads will be reverse complementary to the barcodes in the `737K-cratac-v1.txt`. Therefore, we need to create a reverse complementary file as the whitelist for some data.

```{eval-rst}
.. tip::
  
  Whenever you are dealing with reads that come from ``index 2``, you should:

  1. Check the sequencing machine used to generate the data
  2. Make sure you are familiar with different sequencing modes from different Illumina machines by looking at `this page <https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html>`_.
  3. Extract some sequences from ``index 2``, compare them to the whitelist or the reverse complementary to the whitelist. 
```

## From FastQ To Fragments

Now we are ready to map the reads to the genome using `chromap`:

```console
mkdir -p 10xscATAC/chromap_outs

# map and generate the fragment file

chromap -t 4 --preset atac \
        -x hg38/chromap_index/genome.index \
        -r hg38/hg38.analysisSet.fa \
        -1 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L001_R1_001.fastq.gz,10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L002_R1_001.fastq.gz \
        -2 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L001_R3_001.fastq.gz,10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L002_R3_001.fastq.gz \
        -b 10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L001_R2_001.fastq.gz,10xscATAC/pbmc500/atac_pbmc_500_nextgem_fastqs/atac_pbmc_500_nextgem_S1_L002_R2_001.fastq.gz \
        --barcode-whitelist 10xscATAC/pbmc500/737K-cratac-v1.txt \
        -o 10xscATAC/chromap_outs/fragments.tsv

# compress and index the fragment file

bgzip 10xscATAC/chromap_outs/fragments.tsv
tabix -s 1 -b 2 -e 3 -p bed 10xscATAC/chromap_outs/fragments.tsv.gz
```

Two new files `fragments.tsv.gz` and `fragments.tsv.gz.tbi` are generated. They will be useful and sometimes required for other programs to perform downstream analysis.

### Explain chromap

If you understand the __10x Genomics Single Cell ATAC__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html), the command above should be straightforward to understand.

`-t 4`

>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`-x hg38/chromap_index/genome.index`

>> The `chromap` index file. We are dealing with human samples in this case.

`-r hg38/hg38.analysisSet.fa`

>> Reference genome sequence in `fasta` format. This is basically the file which you used to create the `chromap` index file.

`-1`, `-2` and `-b`

>> They are Read 1 (genomic), Read 2 (genomic) and cell barcode read, respectively. For ATAC-seq, the sequencing is usually done in pair-end mode. Therefore, you normally have two genomic reads for each genomic fragment: Read 1 and Read 2. For the reason described previously, `R1` is the genomic Read 1 and should be passed to `-1`; `R3` is actually the genomic Read 2 and should be passed to `-2`; `R2` is the cell barcode read and should be passed to `-b`. Multiple input files are supported and they can be listed in a comma-separated manner. In that case, they must be in the same order.

`--barcode-whitelist 10xscATAC/pbmc500/737K-cratac-v1.txt`

>> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell ATAC__ is a commercial platform. The whitelist is taken from their commercial software `cellranger-atac`. In this example data, sequencing is done using __NovaSeq 6000 (v1.0)__. Therefore, we use the original whitelist. In other cases, you might want to use the reverse complementary version of the whitelist.

`-o 10xscATAC/chromap_outs/fragments.tsv`

>> Direct the mapped fragments to a file. The format is described in the [10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments).

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

faSize -detailed hg38/hg38.analysisSet.fa | \
    sort -k1,1 > hg38/hg38.chrom.sizes
```

This is the first 5 lines of `hg38/hg38.chrom.sizes`:

```
chr1    248956422
chr10   133797422
chr11   135086622
chr11_KI270721v1_random 100316
chr12   133275309
```

Now let's generate the reads from fragments:

```bash
# we use bedClip to remove reads outside the chromosome boundary
# we also remove reads mapped to the mitochondrial genome (chrM)

zcat 10xscATAC/chromap_outs/fragments.tsv.gz | \
    awk 'BEGIN{OFS="\t"}{print $1, $2, $2+50, $4, ".", "+" "\n" $1, $3-50, $3, $4, ".", "-"}' | \
    sed '/chrM/d' | \
    bedClip stdin hg38/hg38.chrom.sizes stdout | \
    sort -k1,1 -k2,2n | \
    gzip > 10xscATAC/chromap_outs/reads.bed.gz
```

Note we also sort the output reads by `sort -k1,1 -k2,2n`. In this way, the order of chromosomes in the `reads.bed.gz` is the same as that in `hg38.chrom.sizes`, which makes downstream processes easier. The output `reads.bed.gz` are the reads in `bed` format, with the 4th column holding the cell barcodes.

## Peak Calling By MACS2

Now we can use the newly generated read file for the peak calling using `MACS2`:

```console
macs2 callpeak -t 10xscATAC/chromap_outs/reads.bed.gz \
               -g hs -f BED -q 0.01 \
               --nomodel --shift -100 --extsize 200 \
               --keep-dup all \
               -B --SPMR \
               --outdir 10xscATAC/chromap_outs \
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

First, we use the `aggregate_peaks.narrowPeak` file. We only need the first 4 columns (chromosome, start, end, peak ID). You can also remove the peaks that overlap [the black list regions](https://www.nature.com/articles/s41598-019-45839-z). The black list is not available for every species and every build, so I'm not doing it here. We also need to sort the peak to make sure the order of the chromosomes in the peak file is the same as that in the `hg38.chrom.sizes` and `reads.bed.gz` files. Then we could find the overlap by `bedtools`. We need to do this in a specific way to get the number of reads in each peak from each cell:

```bash
# format and sort peaks

cut -f 1-4 10xscATAC/chromap_outs/aggregate_peaks.narrowPeak | \
    sort -k1,1 -k2,2n > 10xscATAC/chromap_outs/aggregate_peaks_sorted.bed

# prepare the overlap

bedtools intersect \
    -a 10xscATAC/chromap_outs/aggregate_peaks_sorted.bed \
    -b 10xscATAC/chromap_outs/reads.bed.gz \
    -wo -sorted -g hg38/hg38.chrom.sizes | \
    sort -k8,8 | \
    bedtools groupby -g 8 -c 4 -o freqdesc | \
    gzip > 10xscATAC/chromap_outs/peak_read_ov.tsv.gz
```

#### Explain Finding Reads In Peaks Per Cell

We start with the command before the first pipe, that is, the intersection part. If you read the manual of the `bedtools intersect`, it should be straightforward to understand. The `-wo` option will output the records in both `-a` file and `-b` file. Since the `reads.bed.gz` file has the cell barcode information at the 4th column, we would get an output with both peak and cell information for the overlap. The `-sorted -g hg38/hg38.chrom.sizes` options make the program use very little memory. Here is an example (top 5 lines) of the output of this part:

```console
chr1	181365	181579	aggregate_peak_1	chr1	181336	181386	GCACGCACAGGTGGTA	.	+	21
chr1	181365	181579	aggregate_peak_1	chr1	181352	181402	ACTACCCGTCCTCAGG	.	-	37
chr1	181365	181579	aggregate_peak_1	chr1	181375	181425	ACTACCCGTCCTCAGG	.	+	50
chr1	181365	181579	aggregate_peak_1	chr1	181418	181468	GCACGCACAGGTGGTA	.	-	50
chr1	181365	181579	aggregate_peak_1	chr1	181423	181473	TTATGTCCATTACTTC	.	+	50
```

We see that the 8th column holds the cell barcode and we want to group them using `bedtools groupby`. Therefore, we need to sort by this column, that is the `sort -k8,8`. When we group by the 8th column, we are interested in how many times each peak appear per group, so we could gather the information of the peak ID (4th column). That is the `-g 8 -c 4 -o freqdesc`. The `-o freqdesc` option returns a `value:frequency` pair in descending order. Here are the first 5 records from `peak_read_ov.tsv.gz`:

```console
AAACGAAAGAATCAAC	aggregate_peak_32073:1
AAACGAAAGACACGGT	aggregate_peak_18915:2,aggregate_peak_41671:2,aggregate_peak_54745:2
AAACGAAAGACACTTC	aggregate_peak_14784:2,aggregate_peak_16349:2,aggregate_peak_37921:2
AAACGAAAGACCATAA	aggregate_peak_51252:2
AAACGAAAGACGACTG	aggregate_peak_32356:2
```

In a way, that is sort of a count matrix in an awkward format. For example:

- The first line means that in cell `AAACGAAAGAATCAAC`, the peak `aggregate_peak_32073` has 1 count. All the rest peaks not mentioned here have 0 counts in this cell.
- The second line means that in cell `AAACGAAAGACACGGT`, the peak `aggregate_peak_18915` has 2 counts, the peak `aggregate_peak_41671` has 2 counts and the peak `aggregate_peak_54745` has 2 counts. All the rest peaks not mentioned here have 0 counts in this cell.

### Output The Peak-By-Cell Matrix

At this stage, we pretty much have all the things needed. Those two files `aggregate_peaks_sorted.bed` and `peak_read_ov.tsv.gz` contain all information for a peak-by-cell count matrix. We just need a final touch to make the output in a standard format: a [market exchange format (MEX)](https://math.nist.gov/MatrixMarket/formats.html). Since most downstream software takes the input from the __10x Genomics Single Cell ATAC__ results, we are going to generate the MEX and the associated files similar to the output from 10x Genomics.

Here, I'm using a python script for this purpose. You don't have to do this. Choose whatever works for you. The point here is to just generate similar files as the __peak-barcode matrix__ described from [the 10x Genomics website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

First, let's make a directory to hold the output files and generate the `peaks.bed` and `barcodes.tsv` files, which are easy to do:

```bash
# create dirctory
mkdir -p 10xscATAC/chromap_outs/raw_peak_bc_matrix

# The 10x Genomics peaks.bed is a 3-column bed file, so we do
cut -f 1-3 10xscATAC/chromap_outs/aggregate_peaks_sorted.bed > \
    10xscATAC/chromap_outs/raw_peak_bc_matrix/peaks.bed

# The barcode is basically the first column of the file peak_read_ov.tsv.gz
zcat 10xscATAC/chromap_outs/peak_read_ov.tsv.gz | \
    cut -f 1 > \
    10xscATAC/chromap_outs/raw_peak_bc_matrix/barcodes.tsv
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
    10xscATAC/chromap_outs/aggregate_peaks_sorted.bed \
    10xscATAC/chromap_outs/raw_peak_bc_matrix/barcodes.tsv \
    10xscATAC/chromap_outs/peak_read_ov.tsv.gz \
    10xscATAC/chromap_outs/raw_peak_bc_matrix
```

After that, you should have the `matrix.mtx` in the `10xscATAC/chromap_outs/raw_peak_bc_matrix` directory.

### Cell Calling (Filter Cell Barcodes)

Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, we could use `starolo`. It is used for scRNA-seq in general, but it does have a cell calling function that takes a directory containing raw mtx and associated files, and return the filtered ones. Since `starsolo` looks for the following three files in the input directory: `matrix.mtx`, `features.tsv` and `barcodes.tsv`. Those are the output from the 10x Genomics scRNA-seq workflow. In this case, we can use `peaks.bed` as our `features.tsv`:

```console
# trick starsolo to use peaks.bed as features.tsv by creating symlink

ln -s peaks.bed 10xscATAC/chromap_outs/raw_peak_bc_matrix/features.tsv

# filter cells using starsolo

STAR --runMode soloCellFiltering \
     10xscATAC/chromap_outs/raw_peak_bc_matrix \
     10xscATAC/chromap_outs/filtered_peak_bc_matrix/ \
     --soloCellFilter EmptyDrops_CR

# rename the new feature.tsv to peaks.bed or just create symlink
ln -s features.tsv 10xscATAC/chromap_outs/filtered_peak_bc_matrix/peaks.bed
```

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/10xscATAC/
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
└── pbmc500
    ├── 737K-cratac-v1_rc.txt
    ├── 737K-cratac-v1.txt
    ├── atac_pbmc_500_nextgem_fastqs
    │   ├── atac_pbmc_500_nextgem_S1_L001_I1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L001_R1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L001_R2_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L001_R3_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L002_I1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L002_R1_001.fastq.gz
    │   ├── atac_pbmc_500_nextgem_S1_L002_R2_001.fastq.gz
    │   └── atac_pbmc_500_nextgem_S1_L002_R3_001.fastq.gz
    └── atac_pbmc_500_nextgem_fastqs.tar

5 directories, 29 files
```

```{eval-rst}
.. note:: How does this workflow compared to cellranger-atac?

  It is certainly much much faster, but what about the results? We can compare them with the ``cellranger-atac`` generated matrix from the `example page <https://www.10xgenomics.com/resources/datasets/500-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0>`_. Just to do a very brief comparison:

  ===========  ===============  ===============  =======
  Comparison   chromap + macs2  cellranger-atac  overlap
  ===========  ===============  ===============  =======
  Peak number      71,000           65,908       61,042
  Cell number       643              484           482
  ===========  ===============  ===============  =======
```
