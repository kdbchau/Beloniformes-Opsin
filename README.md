# Multiple Ecological Axes Drive Cone Opsin Evolution in Beloniformes
This pipeline focuses on the extraction of cone opsin sequences from whole exome sequence data of fishes. It will also focus on multiple sequence alignment using both reference and makeshift reference sequences, cleaning up the multiple sequence alignments, and testing for changes in selection patterns with CODEML from the PAML program. Additional scripts and software are shown for phylogeny reconstruction and ancestral amino acid reconstruction.

This project was conducted at the [University of Toronto Scarborough](https://www.utsc.utoronto.ca/home/) in the [Lovejoy Lab](http://www.utsc.utoronto.ca/~lovejoy/).

Authors:
* Katherine D. Chau (me)
* [Frances E. Hauser](https://fehauser.wordpress.com/)
* [Jacob M. Daane](https://www.daanelab.org/)
* [Matthew P. Harris](http://www.fishbonelab.org/harris/Home.html)
* [Belinda S. W. Chang](https://chang.eeb.utoronto.ca/)
* [Nathan R. Lovejoy](http://www.utsc.utoronto.ca/~lovejoy/)


# Table of Contents
1. [Install software](#1-install-software)
    * [Additional Software](#11-additional-software)
3. [Quality Control](#2-quality-control)
4. [Read Mapping](#3-read-mapping)
    * [Round 1](#31-round-1)
    * [Round 2+](#32-round-2+)
5. [Cleaning Multiple Sequence Alignments](#4-cleaning-multiple-sequence-alignments)
6. [Phylogeny Reconstruction](#5-phylogeny-reconstruction)
7. [Ancestral Habitat and Diet Reconstruction](#6-ancestral-habitat-and-diet-reconstruction)


# 1. Install Software
Software with a * next to the name were already available in the Niagara cluster from Compute Canada.

1. [__FastQC*__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.8)
2. [__Trimmomatic*__](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.38)
3. [__BWA*__](http://bio-bwa.sourceforge.net/) (version 0.7.17)
4. [__SAMtools*__](http://www.htslib.org/) (version 1.9)
5. [__IQTree__](http://www.iqtree.org/) (version 1.6.0)
6. [__MrBayes__](http://nbisweden.github.io/MrBayes/) (version 3.2.6)
7. [__BEAST__](https://beast.community/) (version 1.8.4)
8. [__PAML__](http://abacus.gene.ucl.ac.uk/software/paml.html) (version 4)

## 1.1. Additional Software
This includes software that was used but not installed onto a cluster, rather, they are available online.

A fast, unconstrained Bayesian approximation for inferring selection ([FUBAR](http://www.datamonkey.org/fubar)) was conducted using the [Datamonkey Server](http://www.datamonkey.org/).

# 2. Quality Control

Whole exome sequence data for 36 beloniform species was extracted by [Dr. Jake Daane](https://www.daanelab.org/) and [Dr. Matthew Harris Laboratory](http://www.fishbonelab.org/harris/Home.html). This project is in collaboration with the Harris Lab from Harvard Medical School.

For tissue and DNA extraction method of exome sequences, see this [paper](https://www.cell.com/current-biology/fulltext/S0960-9822(21)01190-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982221011908%3Fshowall%3Dtrue) and the _Targeted sequence capture design_ and _Specimen tissue collection and sequencing library preparation_ section in this [paper](https://www.biorxiv.org/content/10.1101/2021.03.05.434157v1.full).

FastQC was run on each fastq file for each beloniform species using default settings. I used this [script](https://github.com/kdbchau/Beloniformes/tree/main/Scripts/fastqc.sh)

Next, Trimmomatic was run on each file separately.

```
module load java
module load StdEnv/2020
module load trimmomatic

for file in *fastq; do java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 $file "${file%%_*}_trimmed.fastq" ILLUMINACLIP:TruSeq3-SE:2:30:10 SLIDINGWINDOW:4:20 HEADCROP:5
```

Then rerun FastQC on the trimmed files to ensure adapters were removed and sequences look good.

# 3. Read Mapping
The trimmed beloniform exome sequences were then used for mapping against _Oryzias latipes_ (medaka) sequences to extract opsins.

Medaka are beloniforms that have their whole genome sequenced and can be used as a good reference for protein coding sequences.

Medaka have 8 cone opsin sequences readily available on Genbank:

| Opsin Type | Opsin ID | Genbank ID |
| --- | --- | --- |
| red opsin 1 | LWSA | [AB223051](https://www.ncbi.nlm.nih.gov/nuccore/AB223051) |
| red opsin 2 | LWSB | [AB223052](https://www.ncbi.nlm.nih.gov/nuccore/AB223052) |
| green opsin 1 | RH2A | [AB223053](https://www.ncbi.nlm.nih.gov/nuccore/AB223053) |
| green opsin 2 | RH2B | [AB223054](https://www.ncbi.nlm.nih.gov/nuccore/AB223054) |
| green opsin 3 | RH2C | [AB223055](https://www.ncbi.nlm.nih.gov/nuccore/AB223055) |
| blue opsin 1 | SWS2A | [AB223056](https://www.ncbi.nlm.nih.gov/nuccore/AB223056) |
| blue opsin 2 | SWS2B | [AB223057](https://www.ncbi.nlm.nih.gov/nuccore/AB223057) |
| ultraviolet opsin | SWS1 | [AB223058](https://www.ncbi.nlm.nih.gov/nuccore/AB223058) |


## 3.1. Round 1
## 3.2. Round 2+
# 4. Cleaning Multiple Sequence Alignments
# 5. Phylogeny Reconstruction
# 6. Ancestral Habitat and Diet Reconstruction
