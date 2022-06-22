# Multiple Ecological Axes Drive Cone Opsin Evolution in Beloniformes

This project looks at extracting cone opsin sequences from fishes belonging to the order Beloniformes, and testing for changes in selection patterns within the cone opsins. This project uses bioinformatics and PAML to analyze selection patterns and the evolution of vision in fishes that transitioned from marine to freshwater habitats, as well as transitioned from a generalized diet to more specialized diets.

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
All software used in the pipeline were implemented in the Niagara cluster from Compute Canada.

1. [__FastQC__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.8)
2. [__Trimmomatic__](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.38)
3. [__BWA__](http://bio-bwa.sourceforge.net/) (version 0.7.17)
4. [__SAMtools__](http://www.htslib.org/) (version 1.9)
5. [__IQTree__](http://www.iqtree.org/) (version 1.6.0)
6. [__MrBayes__](http://nbisweden.github.io/MrBayes/) (version 3.2.6)
7. [__BEAST__](https://beast.community/) (version 1.8.4)
8. [__PAML__](http://abacus.gene.ucl.ac.uk/software/paml.html) (version 4)

## 1.1. Additional Software
This includes software that was used but not installed onto a cluster, rather, they are available online.

A fast, unconstrained Bayesian approximation for inferring selection ([FUBAR](http://www.datamonkey.org/fubar)) was conducted using the [Datamonkey Server](http://www.datamonkey.org/).

# 2. Quality Control

Whole exome sequence data for 36 beloniform species was extracted by [Dr. Jake Daane](https://www.daanelab.org/) and [Dr. Matthew Harris Laboratory](http://www.fishbonelab.org/harris/Home.html). This project is in collaboration with the Harris Lab from Harvard Medical School.

For tissue and DNA extraction method of exome sequences, see this [paper](https://www.cell.com/current-biology/fulltext/S0960-9822(21)01190-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982221011908%3Fshowall%3Dtrue) and the **Targeted sequence capture design** and **Specimen tissue collection and sequencing library preparation** section in this [paper](https://www.biorxiv.org/content/10.1101/2021.03.05.434157v1.full).

# 3. Read Mapping
## 3.1. Round 1
## 3.2. Round 2+
# 4. Cleaning Multiple Sequence Alignments
# 5. Phylogeny Reconstruction
# 6. Ancestral Habitat and Diet Reconstruction
