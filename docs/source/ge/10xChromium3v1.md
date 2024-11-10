# 10x Genomics Single Cell 3' V1

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html) to see how __10x Genomics Single Cell 3' V1__ libraries are generated experimentally. This is a droplet-based method, where cells are captured inside droplets. At the same time, gel beads with barcoded oligo-dT primer containing UMIs are also captured inside the droplet. Reverse transcription happens inside the droplet. The cells and gel beads are loaded on the microfluidic device at certain concentrations, such that a fraction of droplets contain only one cell __AND__ one bead.

## For Your Own Experiments

The `V1` chemistry is already obsolete, but I'm still providing the preprocessing pipeline for the sake of keeping a record. Although it is highly unlikely that you will do this on your own in future, but just in case, this is the configuration:

| Order | Read             | Cycle | Description   | Comment                                                                            |
|-------|------------------|-------|---------------|------------------------------------------------------------------------------------|
| 1     | Read 1           | >50   | cDNA reads    | Normally, this yields `R1_001.fastq.gz`.                                           |
| 2     | Index 1 (__i7__) | 14    | Cell barcodes | Normally, this yields `I1_001.fastq.gz`, but here, it is called `R2_001.fastq.gz`. |
| 3     | Index 2 (__i5__) | 8     | Sample index  | Normally, this yields `I2_001.fastq.gz`, but here, it is called `I1_001.fastq.gz`. |
| 4     | Read 2           | 10    | UMI           | Normally, this yields `R2_001.fastq.gz`, but here, it is called `R3_001.fastq.gz`. |

Look at the order of the sequencing, as you can see, the first (`R1`), the 2nd (`I1`) and the 4th (`R2`) reads are all important for us. Therefore, you would like to get all of them for each sample based on sample index, that is, the 3rd read (`I2`). You could prepare a `SampleSheet.csv` with the sample index information. Here is an example of `SampleSheet.csv` of a NextSeq run with a sample using standard `i5` indexing primers:

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
75,,,,,,,,,,,
10,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Sample01,,,,,,,,SI-GA-A1_1,AGGCTGGT,,
Sample01,,,,,,,,SI-GA-A1_2,CACAACTA,,
Sample01,,,,,,,,SI-GA-A1_3,GTTGGTCC,,
Sample01,,,,,,,,SI-GA-A1_4,TTGTAAGA,,
```

You can see each sample actually has four different index sequences. This is because each well from the index plate actually contains four different indices for base balancing. To get the reads you need, you should run `bcl2fastq` in the following way:

```console
bcl2fastq --use-bases-mask=Y75,Y14,I8,Y10 \
          --create-fastq-for-index-reads \
          --no-lane-splitting \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r 4 -w 4 -p 4
```

You can check the [bcl2fastq manual](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software/documentation.html) for more information, but the important bit that needs explanation is `--use-bases-mask=Y75,Y14,I8,Y10`. We have four reads, and that parameter specify how we treat each read in the stated order:

1. `Y75` at the first position indicates "use the cycle as a real read", so you will get 75-nt sequences, output as `R1_001.fastq.gz`, because this is the 1st real read.
2. `Y14` at the second position indicates "use the cycle as a real read", so you will get 14-nt sequences, output as `R2_001.fastq.gz`, because this is the 2nd real read.
3. `I8` at the third position indicates "use the cycle as an index read", so you will get 8-nt sequences, output as `I1_001.fastq.gz`, because this is the 1st index read, though it is the 3rd read overall.
4. `Y10` at the fourth position indicates "use the cycle as a real read", so you will get 10-nt sequences, output as `R3_001.fastq.gz`, because this is the 3rd real read, though it is the 4th read overall.

Therefore, you will get four fastq file per sample. Using the examples above, these are the files you should get:

```bash
Sample01_S1_I1_001.fastq.gz # 8 bp: sample index, this is actually I2 in a conventional sense
Sample01_S1_R1_001.fastq.gz # 75 bp: cDNA reads
Sample01_S1_R2_001.fastq.gz # 14 bp: cell barcodes, this is actually I1 in a conventional sense
Sample01_S1_R3_001.fastq.gz # 10 bp: UMI, this is actually R2 in a conventional sense
```

We can safely ignore the `I1` files, but the naming here is really different from our normal usage. The `R1` files are good. The `R2` files here actually mean `I1` in our normal usage. The `R3` files here actually mean `R2` in our normal usage. Anyway, __DO NOT get confused__.

To run `starsolo`, we need to get the cell barcodes and the UMI into the same fastq file. This can be simply achieved by stitching `R2` and `R3` together:

```bash
paste <(zcat Sample01_S1_R2_001.fastq.gz) \
      <(zcat Sample01_S1_R3_001.fastq.gz) | \
      awk -F '\t' '{ if(NR%4==1||NR%4==3) {print $1} else {print $1 $2} }' | \
      gzip > Sample01_S1_CB_UMI.fastq.gz
```

After that, you are ready to go.

## Public Data

The data is from the following paper:

```{eval-rst}
.. note::
  Pijuan-Sala B, Griffiths JA, Guibentif C, Hiscock TW, Jawaid W, Calero-Nieto FJ, Mulas C, Ibarra-Soria X, Tyser RCV, Ho DLL, Reik W, Srinivas S, Simons BD, Nichols J, Marioni JC, Göttgens B (2019) **A single-cell molecular map of mouse gastrulation and early organogenesis.** *Nature* 566:490–495. https://doi.org/10.1038/s41586-019-0933-9
```

where the authors profiled the gene expression levels of single cells from mouse embryos during gastrulation stages. This is a milestone paper that provides a molecular map of cell fate specifications spanning the time window of gastrulation when three germ layers are formed. The data contains over 100,000 cells from early-, mid- and late-gastrulation, including nine time points from E6.5 to E8.5. Note that the 10x Genomics V1 chemistry has already become obsolete. This mouse gastrulation project was conceived in the very early stage of single-cell sequencing era, so the V1 chemistry was used. The data can be found under the accession code [E-MTAB-6967](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6967) from ArrayExpress. Since the whole data set is huge, we are just going to use the first sample `Sample 1` which corresponds to `embryo pool 1` as the demonstration.

This sample was sequenced on the Illumina HiSeq 2500 machine across many different runs and lanes. On top of that, 10x Genomics uses different sample indices even though there is only one sample. Therefore, there are many files associated with `embryo pool 1`, a total of 96 files for this sample alone. We organise them into subdirectories based on their run number. The trick is to take the file [E-MTAB-6967.sdrf.txt](https://www.ebi.ac.uk/biostudies/files/E-MTAB-6967/E-MTAB-6967.sdrf.txt) located at the right-hand side of the ArrayExpress page and use scripts to parse out the URLs to download them systematically. Anyway, I just put the whole `wget` commands here, and you can see that I'm ignoring the `I1` files:

```console
mkdir -p pijuan-sala2019/data/22089 pijuan-sala2019/data/22108 pijuan-sala2019/data/22109 pijuan-sala2019/data/22111
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_AAACGGCG_S1_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_AAACGGCG_S1_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_AAACGGCG_S1_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_CCTACCAT_S2_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_CCTACCAT_S2_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_CCTACCAT_S2_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_AAACGGCG_S45_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_AAACGGCG_S45_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_AAACGGCG_S45_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_CCTACCAT_S46_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_CCTACCAT_S46_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_CCTACCAT_S46_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22089 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22089_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_AAACGGCG_S1_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_AAACGGCG_S1_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_AAACGGCG_S1_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_CCTACCAT_S2_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_CCTACCAT_S2_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_CCTACCAT_S2_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_AAACGGCG_S45_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_AAACGGCG_S45_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_AAACGGCG_S45_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_CCTACCAT_S46_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_CCTACCAT_S46_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_CCTACCAT_S46_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22108 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22108_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_AAACGGCG_S1_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_AAACGGCG_S1_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_AAACGGCG_S1_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_CCTACCAT_S2_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_CCTACCAT_S2_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_CCTACCAT_S2_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_AAACGGCG_S45_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_AAACGGCG_S45_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_AAACGGCG_S45_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_CCTACCAT_S46_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_CCTACCAT_S46_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_CCTACCAT_S46_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22109 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22109_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_AAACGGCG_S1_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_AAACGGCG_S1_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_AAACGGCG_S1_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_CCTACCAT_S2_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_CCTACCAT_S2_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_CCTACCAT_S2_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_AAACGGCG_S45_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_AAACGGCG_S45_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_AAACGGCG_S45_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_CCTACCAT_S46_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_CCTACCAT_S46_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_CCTACCAT_S46_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
wget -P pijuan-sala2019/data/22111 -c ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-6967/22111_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
```

Since we have so many downloads, we should check the `md5` to make sure the download is successful. However, I cannot find the `md5` information, so let's just proceed. You should see the following files:

```console
scg_prep_test/pijuan-sala2019/data/
├── 22089
│   ├── 22089_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   ├── 22089_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   ├── 22089_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   ├── 22089_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   ├── 22089_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   ├── 22089_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   ├── 22089_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   ├── 22089_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   ├── 22089_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   ├── 22089_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   ├── 22089_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   ├── 22089_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   ├── 22089_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   ├── 22089_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   ├── 22089_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   ├── 22089_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   ├── 22089_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   ├── 22089_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   ├── 22089_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   ├── 22089_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   ├── 22089_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   ├── 22089_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   ├── 22089_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   └── 22089_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
├── 22108
│   ├── 22108_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   ├── 22108_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   ├── 22108_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   ├── 22108_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   ├── 22108_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   ├── 22108_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   ├── 22108_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   ├── 22108_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   ├── 22108_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   ├── 22108_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   ├── 22108_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   ├── 22108_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   ├── 22108_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   ├── 22108_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   ├── 22108_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   ├── 22108_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   ├── 22108_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   ├── 22108_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   ├── 22108_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   ├── 22108_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   ├── 22108_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   ├── 22108_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   ├── 22108_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   └── 22108_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
├── 22109
│   ├── 22109_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   ├── 22109_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   ├── 22109_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   ├── 22109_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   ├── 22109_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   ├── 22109_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   ├── 22109_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   ├── 22109_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   ├── 22109_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   ├── 22109_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   ├── 22109_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   ├── 22109_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   ├── 22109_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   ├── 22109_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   ├── 22109_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   ├── 22109_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   ├── 22109_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   ├── 22109_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   ├── 22109_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   ├── 22109_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   ├── 22109_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   ├── 22109_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   ├── 22109_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   └── 22109_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
├── 22111
│   ├── 22111_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   ├── 22111_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   ├── 22111_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   ├── 22111_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   ├── 22111_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   ├── 22111_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   ├── 22111_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   ├── 22111_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   ├── 22111_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   ├── 22111_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   ├── 22111_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   ├── 22111_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   ├── 22111_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   ├── 22111_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   ├── 22111_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   ├── 22111_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   ├── 22111_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   ├── 22111_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   ├── 22111_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   ├── 22111_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   ├── 22111_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   ├── 22111_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   ├── 22111_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   └── 22111_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
└── E-MTAB-6967.sdrf.txt

4 directories, 97 files
```

## Reformat FastQ Files

To use `starsolo`, we need to prepare `FASTQ` files into a file containing cDNA reads and a file with cell barcode + UMI. The cDNA reads are just all those files that match `*_R1_001.fastq.gz`. To get the CB+UMI reads, we need to stitch `R2` and `R3` together. That is, append the 10-bp UMI in `R3` to the 14-bp cell barcodes in `R2`. Since there are so many files, we also combine the same type of files for simplicity, even though it takes more space:

```bash
# for cDNA reads, combine them
cat pijuan-sala2019/data/*/*_R1_001.fastq.gz > pijuan-sala2019/data/cDNA_reads.fastq.gz

# stitch R2 and R3
paste <(zcat pijuan-sala2019/data/*/*_R2_001.fastq.gz) <(zcat pijuan-sala2019/data/*/*_R3_001.fastq.gz) | \
    awk -F "\t" '{if (NR%4==1||NR%4==3) {print $1} else {print $1 $2}}' | \
    gzip > pijuan-sala2019/data/CB_UMI.fastq.gz
```

The resulting files `cDNA_reads.fastq.gz` and `CB_UMI_reads.fastq.gz` are just what we need.

## Prepare Whitelist

The barcodes on the gel beads of the 10x Genomics platform are well defined. We need the information for the `V1` chemistry. If you have `cellranger` in your computer, you will find a file called `737K-april-2014_rc.txt` in the `lib/python/cellranger/barcodes/` directory. If you don't have `cellranger`, I have prepared the file for you:

```console
# download the whitelist 
wget -P pijuan-sala2019/data/ https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-april-2014_rc.txt.gz
gunzip pijuan-sala2019/data/737K-april-2014_rc.txt.gz
```

## From FastQ To Count Matrix

Now we could start the preprocessing by simply doing:

```console
STAR --runThreadN 4 \
     --genomeDir mm10/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix pijuan-sala2019/star_outs/ \
     --readFilesIn pijuan-sala2019/data/cDNA_reads.fastq.gz pijuan-sala2019/data/CB_UMI_reads.fastq.gz \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10 \
     --soloCBwhitelist pijuan-sala2019/data/737K-april-2014_rc.txt \
     --soloCellFilter EmptyDrops_CR \
     --soloStrand Forward \
     --outSAMattributes CB UB \
     --outSAMtype BAM SortedByCoordinate
```

## Explanation

If you understand the __10x Genomics Single Cell 3' V1__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3v1.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir mm10/star_index`

> Pointing to the directory of the star index. The public data from the paper is investigating mouse gastrulation. Therefore, we need to use the mouse reference genome.

`--readFilesCommand zcat`

> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix pijuan-sala2019/star_outs/`

> We want to keep everything organised. This directs all output files inside the `pijuan-sala2019/star_outs` directory.

`--readFilesIn pijuan-sala2019/data/cDNA_reads.fastq.gz pijuan-sala2019/data/CB_UMI_reads.fastq.gz`

> If you check the manual, we should put two files here. The first file is the reads that come from cDNA, and the second the file should contain cell barcode and UMI. We have gone through all the trouble to generate those files using the procedures described above.

`--soloType CB_UMI_Simple`

> Most of the time, you should use this option, and specify the configuration of cell barcodes and UMI in the command line (see immediately below). Sometimes, it is actually easier to prepare the cell barcode and UMI file upfront so that we could use this parameter. That is why went through those procedures to reformat the `fastq` files.

`--soloCBstart 1 --soloCBlen 14 --soloUMIstart 15 --soloUMIlen 10`

> The name of the parameter is pretty much self-explanatory. If using `--soloType CB_UMI_Simple`, we can specify where the cell barcode and UMI start and how long they are in the reads from the first file passed to `--readFilesIn`. Note the position is 1-based (the first base of the read is 1, NOT 0).

`--soloCBwhitelist pijuan-sala2019/data/737K-april-2014_rc.txt`

> The plain text file containing all possible valid cell barcodes, one per line. __10x Genomics Single Cell 3' V1__ is a commercial platform. The whitelist is taken from their commercial software `cellranger`.

`--soloCellFilter EmptyDrops_CR`

> Experiments are never perfect. Even for droplets that do not contain any cell, you may still get some reads. In general, the number of reads from those droplets should be much smaller, often orders of magnitude smaller, than those droplets with cells. In order to identify true cells from the background, you can apply different algorithms. Check the `star` manual for more information. We use `EmptyDrops_CR` which is the most frequently used parameter. 

`--soloStrand Forward`

> The choice of this parameter depends on where the cDNA reads come from, i.e. the reads from the first file passed to `--readFilesIn`. You need to check the experimental protocol. If the cDNA reads are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __10x Genomics Single Cell 3' V1__, the cDNA reads are from the Read 1 file. During the experiment, the mRNA molecules are captured by barcoded oligo-dT primer containing the Illumina Read 2 sequence. Therefore, Read 2 comes from the first strand, complementary to the coding strand. Read 1 comes from the coding strand. Therefore, use `Forward` for __10x Genomics Single Cell 3' V1__ data. This `Forward` parameter is the default, because many protocols generate data like this, but I still specified it here to make it clear.

`--outSAMattributes CB UB`

> We want the cell barcode and UMI sequences in the `CB` and `UB` attributes of the output, respectively. The information will be very helpful for downstream analysis. 

`--outSAMtype BAM SortedByCoordinate`

> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/pijuan-sala2019/
├── data
│   ├── 22089
│   │   ├── 22089_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   │   ├── 22089_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   │   ├── 22089_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   │   ├── 22089_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   │   ├── 22089_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   │   ├── 22089_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   │   ├── 22089_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   │   ├── 22089_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   │   ├── 22089_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   │   ├── 22089_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   │   ├── 22089_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   │   ├── 22089_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   │   ├── 22089_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   │   ├── 22089_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   │   ├── 22089_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   │   ├── 22089_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   │   ├── 22089_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   │   ├── 22089_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   │   ├── 22089_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   │   ├── 22089_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   │   ├── 22089_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   │   ├── 22089_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   │   ├── 22089_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   │   └── 22089_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
│   ├── 22108
│   │   ├── 22108_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   │   ├── 22108_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   │   ├── 22108_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   │   ├── 22108_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   │   ├── 22108_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   │   ├── 22108_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   │   ├── 22108_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   │   ├── 22108_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   │   ├── 22108_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   │   ├── 22108_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   │   ├── 22108_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   │   ├── 22108_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   │   ├── 22108_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   │   ├── 22108_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   │   ├── 22108_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   │   ├── 22108_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   │   ├── 22108_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   │   ├── 22108_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   │   ├── 22108_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   │   ├── 22108_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   │   ├── 22108_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   │   ├── 22108_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   │   ├── 22108_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   │   └── 22108_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
│   ├── 22109
│   │   ├── 22109_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   │   ├── 22109_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   │   ├── 22109_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   │   ├── 22109_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   │   ├── 22109_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   │   ├── 22109_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   │   ├── 22109_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   │   ├── 22109_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   │   ├── 22109_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   │   ├── 22109_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   │   ├── 22109_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   │   ├── 22109_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   │   ├── 22109_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   │   ├── 22109_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   │   ├── 22109_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   │   ├── 22109_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   │   ├── 22109_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   │   ├── 22109_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   │   ├── 22109_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   │   ├── 22109_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   │   ├── 22109_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   │   ├── 22109_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   │   ├── 22109_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   │   └── 22109_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
│   ├── 22111
│   │   ├── 22111_1_AAACGGCG_S1_L001_R1_001.fastq.gz
│   │   ├── 22111_1_AAACGGCG_S1_L001_R2_001.fastq.gz
│   │   ├── 22111_1_AAACGGCG_S1_L001_R3_001.fastq.gz
│   │   ├── 22111_1_CCTACCAT_S2_L001_R1_001.fastq.gz
│   │   ├── 22111_1_CCTACCAT_S2_L001_R2_001.fastq.gz
│   │   ├── 22111_1_CCTACCAT_S2_L001_R3_001.fastq.gz
│   │   ├── 22111_1_GGCGTTTC_S3_L001_R1_001.fastq.gz
│   │   ├── 22111_1_GGCGTTTC_S3_L001_R2_001.fastq.gz
│   │   ├── 22111_1_GGCGTTTC_S3_L001_R3_001.fastq.gz
│   │   ├── 22111_1_TTGTAAGA_S4_L001_R1_001.fastq.gz
│   │   ├── 22111_1_TTGTAAGA_S4_L001_R2_001.fastq.gz
│   │   ├── 22111_1_TTGTAAGA_S4_L001_R3_001.fastq.gz
│   │   ├── 22111_2_AAACGGCG_S45_L002_R1_001.fastq.gz
│   │   ├── 22111_2_AAACGGCG_S45_L002_R2_001.fastq.gz
│   │   ├── 22111_2_AAACGGCG_S45_L002_R3_001.fastq.gz
│   │   ├── 22111_2_CCTACCAT_S46_L002_R1_001.fastq.gz
│   │   ├── 22111_2_CCTACCAT_S46_L002_R2_001.fastq.gz
│   │   ├── 22111_2_CCTACCAT_S46_L002_R3_001.fastq.gz
│   │   ├── 22111_2_GGCGTTTC_S47_L002_R1_001.fastq.gz
│   │   ├── 22111_2_GGCGTTTC_S47_L002_R2_001.fastq.gz
│   │   ├── 22111_2_GGCGTTTC_S47_L002_R3_001.fastq.gz
│   │   ├── 22111_2_TTGTAAGA_S48_L002_R1_001.fastq.gz
│   │   ├── 22111_2_TTGTAAGA_S48_L002_R2_001.fastq.gz
│   │   └── 22111_2_TTGTAAGA_S48_L002_R3_001.fastq.gz
│   ├── 737K-april-2014_rc.txt
│   ├── CB_UMI.fastq.gz
│   └── cDNA_reads.fastq.gz
├── E-MTAB-6967.sdrf.txt
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

13 directories, 127 files
```