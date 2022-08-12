# SMART-seq/SMART-seq2

Check [this GitHub page](https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html) to see how SMART-seq/SMART-seq2 libraries are generated experimentally. Those are plate-based methods, where single cells are sorted into 96- or 384-well plates in a one-cell-per-well manner. The library from each well is prepared separately and indexed using Illumina Nextera Index primers.

## For Your Own Experiments

Your sequencing read configuration is like this:

| Order | Read             | Cycle   | Description                      |
|-------|------------------|---------|----------------------------------|
| 1     | Read 1           | >50     | `R1_001.fastq.gz`, cDNA reads    |
| 2     | Index 1 (__i7__) | 8 or 10 | `I1_001.fastq.gz`, Cell barcodes |
| 3     | Index 2 (__i5__) | 8 or 10 | `I2_001.fastq.gz`, Cell barcodes |
| 4     | Read 2           | >50     | `R2_001.fastq.gz`, cDNA reads    |

If you sequence your data via your core facility or a company, you will need to provide the index sequence to them and they will demultiplex for you. You will get one (single-end) or two (pair-end) `fastq` files per cell.

If you sequence by yourself, you probably need to run `bcl2fatq` by yourself. Prepare a `SampleSheet.csv` with index information for each well (cell), and get the `fastq` for each cell by running `bcl2fastq`. Here is an example of `SampleSheet.csv` of a NextSeq run (only showing 5 random cells):

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
75,,,,,,,,,,,
,,,,,,,,,,,
[Settings],,,,,,,,,,,
,,,,,,,,,,,
[Data],,,,,,,,,,,
Sample_ID,Sample_Name,Sample_Plate,Sample_Well,Index_Plate,Index_Plate_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description
Meso_D3_A01,,,,,,N701,TAAGGCGA,S502,ATAGAGAG,,
Meso_D3_B02,,,,,,N702,CGTACTAG,S503,AGAGGATA,,
Meso_D3_C03,,,,,,N703,AGGCAGAA,S505,CTCCTTAC,,
Meso_D3_D04,,,,,,N704,TCCTGAGC,S506,TATGCAGT,,
Meso_D3_E05,,,,,,N705,GGACTCCT,S507,TACTCCTT,,
```

## Public Data

For the purpose of demonstration, we will use the data from the following paper:

```{eval-rst}
.. note::
  Hagai T, Chen X, Miragaia RJ, Rostom R, Gomes T, Kunowska N, Henriksson J, Park J- E, Proserpio V, Donati G, Bossini-Castillo L, Braga FAV, Naamati G, Fletcher J, Stephenson E, Vegh P, Trynka G, Kondova I, Dennis M, Haniffa M, Nourmohammad A, Lässig M, Teichmann SA (2018) **Gene expression variability across cells and species shapes innate immunity.** *Nature* 563:197–202. https://doi.org/10.1038/s41586-018-0657-2
```

where __SMART-seq2__ was used to investigate the cellular responses of dermal fibroblasts from various species after polyI:C transfection which mimics virus attack.

The whole data can be accessed from ArrayExpress under the accession code [E-MTAB-5920](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5920/). There are a lot of cells, and we only focused on the unstimulated human cells that passed quality control. To get the `URL` of the `fastq` files, we could do:

```console
curl https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-5920/E-MTAB-5920.sdrf.txt | \
    grep HS_unst | \
    grep -w pass | \
    cut -f 38 > fq_url.txt
```

The output `fq_url.txt` contain the links to the `fastq` files, two per cell. There are 114 files (57 cells) there. To download them, do

```console
mkdir -p hagai2018/data
wget -P hagai2018/data -c -i fq_url.txt
```

## From FastQ To Count Matrix

After downloading the data, we are now ready for the data preprocessing. To use `starsolo` with __SMART-seq2__ data, we need to create a manifest, which is a tab-delimited file mapping the `fastq` files to `cell IDs`. To prepare the manifest, you can run the following command in `bash`:

```bash
for i in $(ls -l hagai2018/data/*.gz | cut -f 2 -d/ | cut -f 1 -d '_' | sort -u); \
    do echo -e "hagai2018/data/${i}_1.fastq.gz\thagai2018/data/${i}_2.fastq.gz\t${i}"; \
    done > hagai2018_manifest.tsv
```

These are the first 5 lines of `hagai2018_manifest.tsv`:

```text
hagai2018/data/ERR2060664_1.fastq.gz    hagai2018/data/ERR2060664_2.fastq.gz    ERR2060664
hagai2018/data/ERR2060665_1.fastq.gz    hagai2018/data/ERR2060665_2.fastq.gz    ERR2060665
hagai2018/data/ERR2060666_1.fastq.gz    hagai2018/data/ERR2060666_2.fastq.gz    ERR2060666
hagai2018/data/ERR2060667_1.fastq.gz    hagai2018/data/ERR2060667_2.fastq.gz    ERR2060667
hagai2018/data/ERR2060668_1.fastq.gz    hagai2018/data/ERR2060668_2.fastq.gz    ERR2060668
```

Once we have the manifest, we could start the preprocessing by simply doing:

```console
STAR --runThreadN 4 \
     --genomeDir hg38/star_index \
     --readFilesCommand zcat \
     --outFileNamePrefix hagai2018/star_outs/ \
     --soloType SmartSeq \
     --readFilesManifest hagai2018_manifest.tsv \
     --soloUMIdedup Exact NoDedup \
     --soloStrand Unstranded \
     --outSAMtype BAM SortedByCoordinate
```

It should run successfully.

## Explanation

If you understand the __SMART-seq/SMART-seq2__ experimental procedures described in [this GitHub Page](https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html), the command above should be straightforward to understand.

`--runThreadN 4`
  
>>> Use 4 cores for the preprocessing. Change accordingly if using more or less cores.

`--genomeDir hg38/star_index`

>>> Pointing to the directory of the star index.

`--readFilesCommand zcat`

>>> Since the `fastq` files are in `.gz` format, we need the `zcat` command to extract them on the fly.

`--outFileNamePrefix hagai2018/star_outs/`

>>> We want to keep everything organised. This directs all output files inside the `hagai2018/star_outs` directory.

`--soloType SmartSeq`

>>> The technology used to generate the data. The data generated by __SMART-seq__ and __SMART-seq2__ are exactly the same. The differences are in the experimental conditions, not the sequencing data.

`--readFilesManifest hagai2018_manifest.tsv`

>>> Apparently, the manifest file we just created.

`--soloUMIdedup Exact NoDedup`

>>> The __SMART-seq2__ data do not have UMI in the reads. `Exact` means perform the deduplication using the genomic coordinates, that is, fragments with the exact same starts and ends will be treated as duplicates. `NoDedup` means do not perform deduplication. In ChIP-seq, deduplication is standard. In non-UMI RNA-seq, it seems deduplication is not always enforced (I might be wrong). I'm not sure if this makes a huge difference. Putting both options here will generated two versions of count matrices, one with and one without deduplication.

`--soloStrand Unstranded`

>>> The choice of this parameter depends on where the cDNA reads come from. You need to check the experimental protocol. If they are from the same strand as the mRNA (the coding strand), this parameter will be `Forward` (this is the default). If they are from the opposite strand as the mRNA, which is often called the first strand, this parameter will be `Reverse`. In the case of __SMART-seq2__, we have pair-end mode, and both Read 1 and Read 2 come from cDNA, so we need to ask: which strand does Read 1 come from? The library is generated from full-length double-stranded cDNA. Very often, it is constructed using the Illumina Nextera Library Preparation Kit. Alternatively, Fragmentase + adaptor ligation can also be used. In both cases, the strand information is lost during the procedure: the reads can originate from either strand of the gene. Therefore, use `Unstranded` for __SMART-seq__ and __SMART-seq2__ data.

`--outSAMtype BAM SortedByCoordinate`

>>> We want sorted `BAM` for easy handling by other programs.

If everything goes well, your directory should look the same as the following:

```console
scg_prep_test/hagai2018
├── data
│   ├── ERR2060664_1.fastq.gz
│   ├── ERR2060664_2.fastq.gz
│   ├── ERR2060665_1.fastq.gz
│   ├── ERR2060665_2.fastq.gz
│   ├── ERR2060666_1.fastq.gz
│   ├── ERR2060666_2.fastq.gz
│   ├── ERR2060667_1.fastq.gz
│   ├── ERR2060667_2.fastq.gz
│   ├── ERR2060668_1.fastq.gz
│   ├── ERR2060668_2.fastq.gz
│   ├── ERR2060669_1.fastq.gz
│   ├── ERR2060669_2.fastq.gz
│   ├── ERR2060670_1.fastq.gz
│   ├── ERR2060670_2.fastq.gz
│   ├── ERR2060671_1.fastq.gz
│   ├── ERR2060671_2.fastq.gz
│   ├── ERR2060672_1.fastq.gz
│   ├── ERR2060672_2.fastq.gz
│   ├── ERR2060673_1.fastq.gz
│   ├── ERR2060673_2.fastq.gz
│   ├── ERR2060674_1.fastq.gz
│   ├── ERR2060674_2.fastq.gz
│   ├── ERR2060675_1.fastq.gz
│   ├── ERR2060675_2.fastq.gz
│   ├── ERR2060676_1.fastq.gz
│   ├── ERR2060676_2.fastq.gz
│   ├── ERR2060677_1.fastq.gz
│   ├── ERR2060677_2.fastq.gz
│   ├── ERR2060678_1.fastq.gz
│   ├── ERR2060678_2.fastq.gz
│   ├── ERR2060679_1.fastq.gz
│   ├── ERR2060679_2.fastq.gz
│   ├── ERR2060680_1.fastq.gz
│   ├── ERR2060680_2.fastq.gz
│   ├── ERR2060681_1.fastq.gz
│   ├── ERR2060681_2.fastq.gz
│   ├── ERR2060682_1.fastq.gz
│   ├── ERR2060682_2.fastq.gz
│   ├── ERR2060683_1.fastq.gz
│   ├── ERR2060683_2.fastq.gz
│   ├── ERR2060684_1.fastq.gz
│   ├── ERR2060684_2.fastq.gz
│   ├── ERR2060685_1.fastq.gz
│   ├── ERR2060685_2.fastq.gz
│   ├── ERR2060686_1.fastq.gz
│   ├── ERR2060686_2.fastq.gz
│   ├── ERR2060687_1.fastq.gz
│   ├── ERR2060687_2.fastq.gz
│   ├── ERR2060688_1.fastq.gz
│   ├── ERR2060688_2.fastq.gz
│   ├── ERR2060689_1.fastq.gz
│   ├── ERR2060689_2.fastq.gz
│   ├── ERR2060690_1.fastq.gz
│   ├── ERR2060690_2.fastq.gz
│   ├── ERR2060691_1.fastq.gz
│   ├── ERR2060691_2.fastq.gz
│   ├── ERR2060692_1.fastq.gz
│   ├── ERR2060692_2.fastq.gz
│   ├── ERR2060693_1.fastq.gz
│   ├── ERR2060693_2.fastq.gz
│   ├── ERR2060694_1.fastq.gz
│   ├── ERR2060694_2.fastq.gz
│   ├── ERR2060695_1.fastq.gz
│   ├── ERR2060695_2.fastq.gz
│   ├── ERR2060696_1.fastq.gz
│   ├── ERR2060696_2.fastq.gz
│   ├── ERR2060697_1.fastq.gz
│   ├── ERR2060697_2.fastq.gz
│   ├── ERR2060698_1.fastq.gz
│   ├── ERR2060698_2.fastq.gz
│   ├── ERR2060699_1.fastq.gz
│   ├── ERR2060699_2.fastq.gz
│   ├── ERR2060700_1.fastq.gz
│   ├── ERR2060700_2.fastq.gz
│   ├── ERR2060701_1.fastq.gz
│   ├── ERR2060701_2.fastq.gz
│   ├── ERR2060702_1.fastq.gz
│   ├── ERR2060702_2.fastq.gz
│   ├── ERR2060703_1.fastq.gz
│   ├── ERR2060703_2.fastq.gz
│   ├── ERR2060704_1.fastq.gz
│   ├── ERR2060704_2.fastq.gz
│   ├── ERR2060705_1.fastq.gz
│   ├── ERR2060705_2.fastq.gz
│   ├── ERR2060706_1.fastq.gz
│   ├── ERR2060706_2.fastq.gz
│   ├── ERR2060707_1.fastq.gz
│   ├── ERR2060707_2.fastq.gz
│   ├── ERR2060708_1.fastq.gz
│   ├── ERR2060708_2.fastq.gz
│   ├── ERR2060709_1.fastq.gz
│   ├── ERR2060709_2.fastq.gz
│   ├── ERR2060710_1.fastq.gz
│   ├── ERR2060710_2.fastq.gz
│   ├── ERR2060711_1.fastq.gz
│   ├── ERR2060711_2.fastq.gz
│   ├── ERR2060712_1.fastq.gz
│   ├── ERR2060712_2.fastq.gz
│   ├── ERR2060713_1.fastq.gz
│   ├── ERR2060713_2.fastq.gz
│   ├── ERR2060714_1.fastq.gz
│   ├── ERR2060714_2.fastq.gz
│   ├── ERR2060715_1.fastq.gz
│   ├── ERR2060715_2.fastq.gz
│   ├── ERR2060716_1.fastq.gz
│   ├── ERR2060716_2.fastq.gz
│   ├── ERR2060717_1.fastq.gz
│   ├── ERR2060717_2.fastq.gz
│   ├── ERR2060718_1.fastq.gz
│   ├── ERR2060718_2.fastq.gz
│   ├── ERR2060719_1.fastq.gz
│   ├── ERR2060719_2.fastq.gz
│   ├── ERR2060720_1.fastq.gz
│   └── ERR2060720_2.fastq.gz
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

6 directories, 130 files
```