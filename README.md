# Multiple Ecological Axes Drive Cone Opsin Evolution in Beloniformes
This pipeline focuses on the extraction of cone opsin sequences from whole exome sequence data of fishes. It will also focus on multiple sequence alignment using both reference and makeshift reference sequences, cleaning up the multiple sequence alignments, and testing for changes in selection patterns with CODEML from the PAML program. Additional scripts and software are shown for phylogeny reconstruction and ancestral amino acid reconstruction.

Authors:
* Katherine D. Chau (me)
* [Frances E. Hauser](https://fehauser.wordpress.com/)
* [Jacob M. Daane](https://www.daanelab.org/)
* [Matthew P. Harris](http://www.fishbonelab.org/harris/Home.html)
* [Belinda S. W. Chang](https://chang.eeb.utoronto.ca/lab-members/belinda/)
* [Nathan R. Lovejoy](https://scholar.google.ca/citations?user=kAzTYTsAAAAJ&hl=en)

Acknowledgements:
* [Erik Spence](https://www.linkedin.com/in/erik-spence-25229a31/?originalSubdomain=ca): for help with creation of scripts and optimizing performance on Niagara cluster
* [SciNet](https://www.scinethpc.ca/) - Computations were performed on the [Niagara]() supercomputer at the SciNet HPC Consortium.
    * [Marcelo Ponce et. al (2019). Deploying a Top-100 Supercomputer for Large Parallel Workloads: The Niagara Supercomputer PEARC'19 Proceedings](https://dl.acm.org/doi/10.1145/3332186.3332195)
    * [Chris Loken et al. (2010). SciNet: Lessons learned from Building a Power-efficient Top-20 System and Data Centre. J. Phys. Conf. Ser.](https://iopscience.iop.org/article/10.1088/1742-6596/256/1/012026)
* [Lovejoy Lab](http://www.utsc.utoronto.ca/~lovejoy/) - For help with additional computations on the Lovejoy server
* [Chang Lab](https://chang.eeb.utoronto.ca/) - Editing and help with analysis of codeml output and opsin evolution

# Table of Contents
1. [Install software](#1-install-software)
    * [1.1. Additional Software](#11-additional-software)
2. [Quality Control](#2-quality-control)
3. [Read Mapping](#3-read-mapping)
    * [3.1. Round 1](#31-round-1)
    * [3.2. Round 2+](#32-round-2)
4. [Phylogeny Reconstruction (Gene Trees)](#4-phylogeny-reconstruction-gene-trees)
    * [4.1. MrBayes](#41-mrbayes)
    * [4.2. IQ-TREE](#42-iq-tree)
5. [Ancestral Habitat and Diet Reconstruction [BEAST]](#5-ancestral-habitat-and-diet-reconstruction-beast)
6. [Species Tree](#6-species-tree)
7. [Cleaning Multiple Sequence Alignments](#7-cleaning-multiple-sequence-alignments)
8. [Molecular Evolutionary Analyses](#8-molecular-evolutionary-analyses)
    * [8.1. Random Sites and FUBAR](#81-random-sites-and-fubar)
    * [8.2. Clade and Branch Models](#82-clade-and-branch-models)
    * [8.3. Ancestral Amino Acid Reconstruction](#83-ancestral-amino-acid-reconstruction)

# 1. Install Software
Software with a * next to the name were already available in the Niagara cluster from Compute Canada. PAML was run on the Lovejoy server.

1. [__FastQC*__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.8)
2. [__Trimmomatic*__](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.38)
3. [__BWA*__](http://bio-bwa.sourceforge.net/) (version 0.7.17)
4. [__SAMtools*__](http://www.htslib.org/) (version 1.9)
5. [__IQTree__](http://www.iqtree.org/) (version 1.6.0)
6. [__MrBayes__](http://nbisweden.github.io/MrBayes/) (version 3.2.6)
7. [__BEAST__](https://beast.community/) (version 1.8.4)
8. [__Tracer__](https://beast.community/tracer) (version 1.5)
9. [__TreeAnnotator__](https://beast.community/treeannotator) (version 2.5.2 - part of BEAST package)
10. [__Mesquite__](https://www.mesquiteproject.org/Installation.html) (version 3.51)
11. [__PAML__](http://abacus.gene.ucl.ac.uk/software/paml.html) (version 4)
12. [__FigTree__](http://tree.bio.ed.ac.uk/software/figtree/)

## 1.1. Additional Software
This includes software that was used but not installed onto a cluster or personal computer, rather, they are available online.

A fast, unconstrained Bayesian approximation for inferring selection ([FUBAR](http://www.datamonkey.org/fubar)) was conducted using the [Datamonkey Server](http://www.datamonkey.org/).

# 2. Quality Control

Whole exome sequence data for 36 beloniform species was extracted by [Jacob M. Daane](https://www.daanelab.org/) and [Dr. Matthew Harris Laboratory](http://www.fishbonelab.org/harris/Home.html). This project is in collaboration with the Harris Lab from Harvard Medical School.

For tissue and DNA extraction method of exome sequences, see this paper by [Daane et al., 2021 in _Current Biology_](https://www.cell.com/current-biology/fulltext/S0960-9822(21)01190-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0960982221011908%3Fshowall%3Dtrue) and the _Targeted sequence capture design_ and _Specimen tissue collection and sequencing library preparation_ section in [Daane et al., 2021 in _bioRxiv_](https://www.biorxiv.org/content/10.1101/2021.03.05.434157v1.full).

FastQC was run on each fastq file for each beloniform species using default settings. I used the [fastqc.sh](https://github.com/kdbchau/Beloniformes/tree/main/Scripts/fastqc.sh) script.

Next, Trimmomatic was run on each file separately.

```
module load java
module load StdEnv/2020
module load trimmomatic

for file in *fastq; do java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar SE -phred33 $file "${file%%_*}_trimmed.fastq" ILLUMINACLIP:TruSeq3-SE:2:30:10 SLIDINGWINDOW:4:20 HEADCROP:5
```
Then rerun FastQC on the trimmed files to ensure adapters were removed and sequences look good. Although I didn't use it, I recommend [MultiQC](https://multiqc.info/) for easier visualization of quality cleanup in section 2. I just checked the fastqc output for the trimmed files to ensure adapters were removed and base quality was high, which it was. There were some warnings such as sequenceduplication in a couple species and some species having variable bases at the beginning (since I didn't use HEADCROP in trimmomatic) but I didn't want to trim because I wanted to keep as much information as possible for protein-coding gene assembly. 

# 3. Read Mapping
The trimmed beloniform exome sequences were then used for mapping against _Oryzias latipes_ (freshwater medaka) sequences to extract opsins.

Freshwater medaka are beloniforms that have their whole genome sequenced and can be used as a good reference for protein coding sequences.

Freshwater medaka have 8 cone opsin sequences readily available on [Genbank](https://www.ncbi.nlm.nih.gov/genbank/):

### Freshwater Medaka (_Oryzias latipes_) Opsins
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

The opsin sequence names and wavelengths discussed in the paper were obtained from [Matsumoto et al. 2006 _Gene_ paper](https://www.sciencedirect.com/science/article/abs/pii/S0378111905007456?via%3Dihub). See image below.

![medaka opsins](https://github.com/kdbchau/Beloniformes/blob/main/Images/medaka_wavelengths.png)

We will also later integrate the marine medaka (_Oryzias melastigma_) sequences which are available at [Ensembl](https://useast.ensembl.org/index.html) and [Genbank](https://www.ncbi.nlm.nih.gov/genbank/). Note: it is called "Indian medaka" on Ensembl.

### Marine Medaka (_Oryzias melastigma_) Opsins
| Opsin Type | Opsin ID | Genbank ID |
| --- | --- | --- |
| red opsin 1 | LWSA | [XM_024269265.2](https://www.ncbi.nlm.nih.gov/nucleotide/XM_024269265.2?report=genbank&log$=nucltop&blast_rank=1&RID=B9DZXEVN016) |
| red opsin 2 | LWSB | [XM_024269264.2](https://www.ncbi.nlm.nih.gov/nucleotide/XM_024269264.2?report=genbank&log$=nucltop&blast_rank=1&RID=B9E6Y2JF016) |
| green opsin 1 | RH2A | [ENSOMET00000006525.1](https://useast.ensembl.org/Oryzias_melastigma/Transcript/Summary?db=core;g=ENSOMEG00000007059;r=NVQA01000009.1:3707619-3711640;t=ENSOMET00000006525) |
| green opsin 2 | RH2B | [ENSOMET00000006471.1](https://useast.ensembl.org/Oryzias_melastigma/Transcript/Summary?db=core;g=ENSOMEG00000007091;r=NVQA01000009.1:3694457-3696234;t=ENSOMET00000006471) |
| green opsin 3 | RH2C | [ENSOMET00000006376.1](https://useast.ensembl.org/Oryzias_melastigma/Transcript/Summary?db=core;g=ENSOMEG00000007148;r=NVQA01000009.1:3683498-3685163;t=ENSOMET00000006376) |
| blue opsin 1 | SWS2A | [XM_024270161.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_024270161.1) |
| blue opsin 2 | SWS2B | [ENSOMET00000036756.1](https://useast.ensembl.org/Oryzias_melastigma/Transcript/Summary?db=core;g=ENSOMEG00000023704;r=NVQA01000009.1:29646285-29652531;t=ENSOMET00000036756) |
| ultraviolet opsin | SWS1 | [ENSOMET00000035584.1](https://useast.ensembl.org/Oryzias_melastigma/Transcript/Summary?db=core;g=ENSOMEG00000021044;r=NVQA01000003.1:19506721-19509232;t=ENSOMET00000035584) |

## 3.1. Round 1

In the very first round of read mapping, use the protein-coding opsin sequence from freshwater medaka. I used [mapping.sh](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/mapping.sh). This was run on the Niagara cluster from ComputeCanada. BWA is used as the mapping software.
See BWA paper by [Heng Li, 2013 _arXiv_](https://arxiv.org/pdf/1303.3997.pdf). 

However, upon trial and error it was found that because of high sequence similarity in the LWS opsins (see this [paper](https://www.nature.com/articles/s41598-019-39978-6)), and in the RH2 opsins for the medakas, I used zebrafish (_Danio rerio_) LWS and guppy RH2 (_Poecilia reticulata_) opsins in the first round of read mapping.

```
# Ran the script as follows on the cluster:

# Sbatch to set up the script on the job queue, then gave it a job name "lwsa", then called the reference sequence (a fasta file with just the freshwater medaka opsin sequence), and an output filename.
sbatch --job-name="lwsa" mapping.sh sws2a_medaka.fa sws2a_medaka_output
```
The above script will use a python script called [consensus.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/consensus.py) which will then take all the mpileup files generated after the mapping, and create a multiple sequence alignment (MSA) which we can alter edit.


## 3.2. Round 2+
Sometimes, using just freshwater medaka as a reference is not good enough. Because medaka are quite divergent from the rest of the beloniformes, we can use information ffrom the first round mapping to create "makeshift" or "chimeric" reference sequences that should be helpful to fill in more gaps. See this [figure](https://github.com/kdbchau/Beloniformes/blob/main/Images/readassembly.pdf) for a visual on the process.

The second and subsequent rounds use these makeshift sequences that are filled in with nucleotides from the top most complete beloniform (after medaka which would be 100% complete because it is a direct match to itself).

To identify the next top beloniform, we used the [ideal_species.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/ideal_species.py) script. But this might not come out sorted properly so we run it the following way in terminal:

```
# name is the name of the file running without its extension, for example if it was rh2a_medaka_MSA.fa
name="rh2a_medaka_MSA" && python ideal_species.py ${name}.fa > ${name}.txt && cat ${name}.txt | sort -nr -k2 -o ${name}.txt
```

A textfile will be created with the species in the first column, sorted by the number of codons in the second column. Medakas are usually always the first two, followed by the beloniform with the next most complete sequence. Then using that beloniform's sequence, manually fill in its gaps with either the medaka sequence or a different nucleotide/amino acid from another beloniform. Typically this is the best approach, as refilling in the gaps with medaka bases doesn't pull out more information. Rather, fill it with bases from the seccond most complete beloniform, or in some cases, if a paritcular family (e.g. all the flyingfishes) are missing codons in an area, fill the gaps with their bases to try and pull out more of their sequences in the next round mapping.

In my study, this is how many round mappings were needed to obtain opsin MSAs that were fairly complete:

| Opsin | # of rounds | Round 1 reference species|
| --- | --- | --- |
| LWSA | 3 | zebrafish LWS1 [KT008400.1](https://www.ncbi.nlm.nih.gov/nuccore/KT008400.1)|
| LWSB | 2 | zebrafish LWS2 [KT008401.1](https://www.ncbi.nlm.nih.gov/nuccore/KT008401.1)|
| RH2A | 1 | guppy RH2A [ENSPRET00000003992.1 ](https://www.ensembl.org/Poecilia_reticulata/Transcript/Summary?db=core;g=ENSPREG00000002751;r=LG5:3901638-3906113;t=ENSPRET00000003992)|
| RH2B | 1 | medaka RH2B [AB223054](https://www.ncbi.nlm.nih.gov/nuccore/AB223054)|
| RH2C | 1 | guppy RH2B/C [ENSPRET00000004125.1](https://www.ensembl.org/Poecilia_reticulata/Transcript/Summary?db=core;g=ENSPREG00000002827;r=LG5:3918624-3921999;t=ENSPRET00000004125)|
| SWS2A | 3 | medaka SWS2A [AB223056](https://www.ncbi.nlm.nih.gov/nuccore/AB223056)|
| SWS2B | 3 | medaka SWS2B [AB223057](https://www.ncbi.nlm.nih.gov/nuccore/AB223057)|
| SWS1 | 3 | medaka SWS1 [ENSOMET00000035584.1](https://useast.ensembl.org/Oryzias_melastigma/Transcript/Summary?db=core;g=ENSOMEG00000021044;r=NVQA01000003.1:19506721-19509232;t=ENSOMET00000035584)|

After completing the round mapping, we will notice that some parts of the alignment are more complete during the first round and other parts are more complete during the second part (becuase different starting references can cause this discrepancy). In this case, we shall merge the MSAs to obtain the most complete MSA for each opsin, then conduct further alignment cleaning prior to codeml analysis.

Merging is done with [merge_seqs.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/merge_seqs.py) which will take one MSA, and fill any gaps it comes across with bases if they are present in the second called MSA. Also, sometimes very minor base differences can occur between the MSAs, (like TTT or TTC) - this script will keep the base from the most up-to-date MSA which is called second.

```
python merge_seqs.py round1_MSA round3_MSA # where the second MSA called is the most up-to-date OR most complete (fewest gaps)
```

After this, we will have our opsin MSAs to work with. But despite all this read mapping and refining, there are still regions in the MSA that are very gapped or have premature stop codons and need to be removed for effective codeml analysis.

# 4. Phylogeny Reconstruction [Gene Trees]

### _For species tree construction, see [section 6](#6-species-tree)._

Using the full MSAs (prior to cleaning), we can construct our gene tree phylogenies. I used IQTREE for a maximum likelihood reconstruction and MrBayes for a Bayesion reconstruction.

In both cases, all cone opsin MSAs were combined into one file and aligned as translated amino acids using MUSCLE alignment (use any program for this - I used [AliView](https://ormbunkar.se/aliview/)), and I included a few outgroup species for better root determination. I included the zebrafish and guppy opsin reference sequences I used for round 1 mapping.

For the aligned fasta file of cone opsins combined, download file [here](https://github.com/kdbchau/Beloniformes/blob/main/Downloadables/combined_cone_opsins.fa).

Programs like [AliView](https://ormbunkar.se/aliview/)) was used to convert the aligned MSA into Phylip foramt, which is needed for MrBayes and IQ-TREE (but IQ-TREE can use a fasta format as well).

## 4.1. MrBayes
To run MrBayes, a nexus file must be created. Templates can be found [here](https://github.com/NBISweden/MrBayes/tree/v3.2.7a/examples). 

For my data, I set dimensions to:
* ntax=310 (because I combined 8 cone opsins each with 38 beloniforms, 8\*38=304, plus guppy RH2BC and RH2A and zebrafish LWS1 and LWS2; totally 310 species).
* nchar=1104 (after aligning everything, the total length of the MSA was 1104 bases).
* format gap = - (hyphens represent gaps)
* datatype = dna (I used dna sequences, not amino acid although everythign was in proper reaading frame)
* missing = ? (for any missing info - my alignments just had gaps, no question marks).

Then using the phylip file version of the MSA, input it into the nexus file where it says to. Afterwards there is a MrBayes block with additional parameters that need to be set. These are my settings:
* outgroup Danio_rerio_LWS2;
* lset nst=6 rates=gamma;
* set autoclose=yes nowarn=yes;
* mcmc ngen=1000000 printfreq=100 samplefreq=100 nruns=2 nchains=4 savebrlens=yes;
* sump; (leaving this blank uses default of 25% burn-in)
* sumt nruns=2 ntrees=1 contype=halfcompat conformat=simple;

Then double click the MrBayes.exe and a terminal window should open up.
![](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes.png)

Once executed, it will run for a while and then complete. The run was successful if met with convergence, noted by a plot produced at the end that looks randomized.
In my case, it did have an upward trend but my standard deviation values approached 0...

![](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes_StandardDeviations.png) 

...and my PSRF values were around 1.

![](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes_PSRF.png)

This is okay although this can be improved using a different burn-in, more runs, etc. The tree ends up looking the same as IQ-TREE so not a huge problem.

## 4.2. IQ-TREE

Ultrafast bootstrap IQ-TREE is simpler, once installed it can be called from a terminal.

Simply run it as:
```
iqtree -s combined_cone_opsins.fa -bb 1000  # -bb 1000 to get node values
```

The output ```.contree``` will have the node bootstrap values.

Both the Bayesian and maximum likelihood tree can be visualized using any tree visualizing program like FigTree. Here are the figures for the Beloniformes [Bayesian cone opsin tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes_AllConeOpsins.nexus.con.tre.pdf) 
and the [maximum likelihood cone opsin tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/IQTREE_AllOpsinsCombined_noRH1.fa.contree.pdf).


# 5. Ancestral Habitat and Diet Reconstruction [BEAST]

To best determine when marine to freshwater transitions occurred within the Beloniformes, using a larger tree is best for ancestral habitat reconstruction.
[Bloom & Lovejoy 2017](https://devinbloom.files.wordpress.com/2017/08/bloom-lovejoy-2017-jbiogeography.pdf) had constructed a phylogeny for 104 beloniforms and 7 outgroup species using RAG1, RAG2, TMO-4C4, and CYTB genes. We obtained these sequences and I included TMO-4C4 for _Hemiramphodon pogonognathus_, and all four genes for _Rhynchorhamphus georgii_ and _Cheilopogon papilio_ which were missing. Our tree totalled 120 beloniforms and 7 outgroups.

I followed the procedure for a maximum credibility phylogeny reconstruction using BEAST according to Bloom & Lovejoy 2017 paper, using the descriptions given in their supporting information.

In brief:
* A relaxed log-normal tree prior was used which allows rates to differ across brancches
* A birth-death prior for variable speciation rates
* I partitioned each gene using a general time reversible with gamma distribution (GTR + G)
* Ran everything three times with MCMC length of 60,000,000 and 10% burn-in
* I ran this with random starting seeds to allow random initial parameters.

After the three runs, I used [Tracer 1.5](https://beast.community/tracer) to mix all three runs and verify convergence. Following this, [TreeAnnotator](https://beast.community/treeannotator) was used to generate a maximum credibility (MC) tree. This tree was used for ancestral habitat reconstruction.

Then, using [Mesquite](https://www.mesquiteproject.org/home.html), I added habitat characters (i.e. freshwater or marine) for all 127 species in the MC tree. Habitat reconstruction was run using two methods:
1. Maximum likelihood - MK1 model
2. Maximum parsimony

While habitat reconstruction was identical for majority of the species, it slightly differed for needlefishes. I opted to use only the maximum likelihood method because it corroborated past research that showed multiple individual transitions from marine to freshwater habitats. See [Lovejoy & Collette 2001 in _Copeia_](https://www.jstor.org/stable/1447878).

Here is how the large Beloniformes tree looks under maximum likelihood habitat reconstruction, followed by a pruned version of the tree down to my 38 species.

![large tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/LargeTreeBeloniformHabitat.png)

For a higher resolution download image [here](https://github.com/kdbchau/Beloniformes/blob/main/Images/ML_bigtree_habitatreconstruction.pdf)

![pruned tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/Pruned_BelTree_ML.jpg)

Using the pruned tree and literature research, reconstructing diet was easier as we based it off each family. Needlefishes are primarily piscivorous (except for _Belonion dibranchodon_ and _Potamorrhaphis guianensis_), flyingfishes and sauries feed on zooplankton, halfbeaks are herbivorous feeding on plants and algae. Medakas, viviparous halfbeaks, and the two non-piscivorous needlefishes feed on insects, algae or other small taxa like zooplankton, so they are classified as generalists. 

Below is the figure showing both habitat and diet transitions for the 38 Beloniformes.
![diet and habitat beloniformes](https://github.com/kdbchau/Beloniformes/blob/main/Images/Partitions_Clean_3Diet.jpg)

# 6. Species Tree
To generate a species tree for the 38 beloniforms, we had obtained an alignment from the whole exome sequencing data that encompassed single-copy exons >100bp, with at least 85% coverage in all species, concatenated. The total length of the alignment was 1,579,692 bases encompassing 8,768 exons. This alignment included the freshwater medaka but did not include the marine medaka. Because marine medaka is most closely related to the freshwater medaka, this was easy to manually add in as an outgroup.

Using this tree, I ran MrBayes and IQ-TREE on the alignment (same process as above) to obtain my tree. The only difference here is MrBayes has a limit of 99,999 bases, so I trimmed the alignment down to 99,999 bases. Both methods produced identical topologies - resulting in one coherent species tree. ![species tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/Bel_0.85_species_tree.jpg)

Branch colouring is based of the ancestral habitat reconstruction using maximum likelihood. Nodes on branches show bootstrap/Bayesian posterior probabilities. Scale represents number of substitutions per nucleotide.

# 7. Cleaning Multiple Sequence Alignments
Premature stop codons or highly gapped regions will hinder calculations in codeml from the PAML program. It is advised to remove all stop codons (the very end of a protein-coding sequence; TAA, TAG, or TGA codons) or manually delete positions that are highly gapped.

To be as accurate as possible and keep things consistent between the 8 cone opsin MSAs, use [alignment_editor.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/alignment_editor.py) which will do three things:
1. Delete columns with a certain % of stop codons (premature stop codons from erroneous sequencing)
2. Columns with premature stop codons that are not removed, will convert the stop codon into a gap ("---") (because codeml will not run otherwise).
3. Delete columns with gaps exceeding a certain %.

For this work, we delete codon positions that have 10%+ stop codons present, and 30%+ gaps (including incomplete codons such as AT-, A--).

To run this script:

```
python alignment_editor.py input_msa.fa 0.1 0.3 output_msa.fa # if using 10% stops and 30% gaps removal
```

Eight files are moved into a new folder with the opsin name as the folder name and a #### attached to indicate the stop-gaps % used (i.e. rh2a_msa_1030):

| # | File | Example | Description |
| --- | --- | --- | --- |
|1 |__input_msa.fa__ |rh2a_msa.fa | original, unedited msa|
|2 | __output_msa\_####.fa__| rh2a_msa_1030.fa| cleaned msa; value indicates cleanup thresholds use (i.e. 1030 = 10% stop codon removal and 30 for 30% gap removal)|
|3 |__output_msa\_####.phy__ | rh2a_msa_1030.phy | Phylip version of the output_msa.fa|
|4 |__output_msa\_####\_RemovedAttributes.xlsx__ | rh2a_msa_1030_RemovedAttributed.xlsx| Collection of the textfiles that show the codons removed/edited|
|5 |__output_msa\_####\_RemovedChange.txt__ |rh2a_msa_1030_RemovedChange.txt | Textfile of codon positions where stops that were present but not removed changed to a gap|
|6 |__output_msa\_####\_RemovedGaps.txt__ | rh2a_msa_1030_RemovedGaps.txt|Textfile of codon positions with removed gaps |
|7 |__output_msa\_####\_RemovedStops.txt__ |rh2a_msa_1030_RemovedStops.txt | Textfile of codon positions with removed stops|
|8 | __output_msa\_####\_STATS.txt__|rh2a_msa_1030_STATS.txt | Overall stats of the cleanup, shows how many total bases removed, most common codon in remaining msa, etc.|

# 8. Molecular Evolutionary Analyses
In this section, we implement CODEML from the PAML program to determine patterns of selection in each of the cone opsins and how they may have evolved with respect to changes in habitat and diet.

In order to run codeml, we need to set up codeml.ctl files for each opsin and create a newick tree file for the program to call with each alignment. Labeling of branches is done in [TreeView](https://treeview.co.uk/download-file/?v=2) but any editor program will work, including Notepad.

Note: All CODEML analyses were run three times using different starting kappa and omega values. Out of these three, the models with the lowest lnL values were used for drawing conclusions.

* kappa = 1, omega = 0.5
* kappa = 2, omega = 1
* kappa = 3, omega = 2

## 8.1. Random-Sites and FUBAR
Random-sites models (M0, M1a, M2a, M3, M7, M8, and M8a) was used to assess positive selection for each group. Our data has five groups being assessed. Species belonging to each group were extracted from each cone opsin and run individually. For example, out of 38  beloniforms, 10 are freshwater and 28 are marine - so each cone opsin was split into two MSAs (each with their own mini newick tree file) for freshwater and marine. 

There was an alignment and tree (generated in Mesquite which allows one to remove taxa and recreate the tree) for: freshwater, marine, piscivores, zooplanktivores, and herbivores.

Positive selection is determined when M2a vs. M1a is significant, M8 vs. M7 was significant, as well as M8 vs. M8a was significant. For a full explanation on testing positive selection using these models, see [Yang 2007](https://academic.oup.com/mbe/article/24/8/1586/1103731). Model M3 vs. M0 just assesses variation for three site classes.

The same alignments and trees were used for FUBAR analyses on the [Datamonkey Server](http://www.datamonkey.org/fubar). FUBAR allows for synonymous substitutions (dS) to be estimated independently, whereas random-sites models in CODEML does not account for selection on synonymous sites.

An example of my codeml.ctl file is shown [here](https://github.com/kdbchau/Beloniformes/blob/main/Downloadables/codeml.ctl) for the freshwater alignment for LWSA using kappa 1 and omega 0.5. An example of the treefile is shown [here](https://github.com/kdbchau/Beloniformes/blob/main/Downloadables/bel_freshwater.phy).

## 8.2. Clade and Branch Models
Clade and branch models are run in the same manner, using a codeml.ctl file and a treefile. However, this time the full MSA is used (no separating groups) and instead the branches are labeled to indicate the "foreground" partition. For more explanation on PAML and partitions, see [Yang's PAML Doc](http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf).

In this case, five different labelings were used. The fifth one includes multiple foreground partitions (which cannot be used for branch-site models).
1. Foreground = freshwater, background = marine
2. Foreground = piscivores, background = remaining beloniforms
3. Foreground = zooplanktivores, background = remaining beloniforms
4. Foreground = herbivores, background = remaining beloniforms
5. Foreground 1 = piscivores, foreground 2 = herbivores, foreground 3 = zooplanktivores, background = generalists.

![partitions](https://github.com/kdbchau/Beloniformes/blob/main/Images/Partitions_Clean_3Diet.jpg)

Because we are comparing many different scenerios, we also implemented Akaike Information Criterion ([AIC](https://ieeexplore.ieee.org/document/1100705)). Models with ΔAIC < 2 are deemed the most likely to explain the observed patterns of selection for each cone opsin.

To calculate ΔAIC, the following was used from parameters obtained in CODEML output files (log file names are specified in codeml.ctl).

$AIC = (-2 * lnL) + (2 * np)$

For example looking at one cone opsin, once AIC is calculated for each of the 5 models, the model with the lowest AIC is identified as the best fit model, and we use its AIC value to see if any other models have ΔAIC < 2. Simply take the AIC value for a model and subtract it from the model with the lowest AIC (this means the best fit model will automatically have ΔAIC = 0).

Here is an example of my CmC model output and values for SWS1 (image below).

![CmC for SWS1](https://github.com/kdbchau/Beloniformes/blob/main/Images/CmCAIC.png).

Here, setting the zooplanktivores as the foreground only was the best fit model for the ultraviolet-sensitive opsin SWS1. However, using habitat where freshwater beloniforms were the foreground can also be considered a good model to explain variation in SWS1 codons, since its ΔAIC was less than 2. 


## 8.3. Ancestral Amino Acid Reconstruction
The last step was reconstructing ancestral amino acids for each opsin. First, we used AIC from [ModelFinder]([http://www.iqtree.org/ModelFinder/](http://www.iqtree.org/doc/Tutorial)) in IQ-TREE to determine the best model to explain amino acid evolution for each opsin.

```
# Example of how to get the best model for a cone opsin
iqtree -s lwsa_msa_1030.fa -st NT2AA -m TEST -AIC -nt 1  #
```

This will produce a log file output that will show testing for all the models. At the bottom there will be several output for different algorithms deciding the best amino acid model. We use AIC.

Here were the best models for each opsin. Based on this, we need to specifiy the approprate dat file in the codeml.ctl file for ancestral amino acid reconstruction. Here is an example of a codeml.ctl file for [amino acid reconstruction]() and changing ```aaRatefile``` variable to call the appropriate amino acid model saved in the codeml bin folder.

|Opsin|AIC Model| aaRatefile .dat file|
|---|---|---|
| LWSA | LG+F+I+G4  | lg.dat |
| LWSB | LG+F+I+G4  |lg.dat |
| RH2A |mtZOA+F+I+G4  | MtZoa.dat|
| RH2B | LG+F+I+G4  |lg.dat |
| RH2C | LG+F+I+G4 |lg.dat |
| SWS2A |LG+F+I+G4  |lg.dat |
| SWS2B |JTT+F+I+G4  |jones.dat |
| SWS1 |JTTDCMut+F+I+G4  |jones-dcmut.dat |

The outputs will produce several files, and every amino acid reconstructions will be shown in the ```rst``` file, for every branch in the tree. A snippet of how that looks:

![](https://github.com/kdbchau/Beloniformes/blob/main/Images/AAreconstructions.png)

Here, there will be shown the  branch the reconstructions are off (based on the input tree, this file will also give an output tree at the top with branch labels for easier viewing in a program like TreeView or FigTree). Below the branch label are the amino acid positions, followed by the amino acid, its posterior probability and conservion to the next amino acid. If it is a terminal branch, the amino acid will not have a posterior probability because it is not "guessing" that amino acid (since it is the actual sequence).

In the figure above, we see that branch 34 (nodes 55 to 56) have several amino acid reconstructions that have posterior probablities all above 80%. However, we are only interested in the codon positions that were identified to be under positive selection according to M8, FUBAR, and CmC models.

In CODEML M8 output, this can be found in the output file with Bayes' Empircal Bayes (BEB) values shown after scrolling to the M8 model. It will look like this:

![](https://github.com/kdbchau/Beloniformes/blob/main/Images/BEB_m8.png)

Where the first column is the codon position according to the first species sequence (which was _2-Ablennes hians_). The second column is the codon, third is posterior probablitity that that codon is under positive selection. Next column is the omega (dN/dS) ratio +- the standard deviation. Those with ** or * are codons under significant positive selection.

This is compared against FUBAR output which is downloaded off the server. Finally, CmC also outputs BEB values for every codon position in the ```rst``` file, with similar posterior probabilties and omega values. If there is a codon that is identified as under significant positive selection from all three methods, it is very likely to be experiecing positive selection.

Those sites are focused on the ancestral amino acid reconstructions. For example if I know that codon at position 217 is under positive selection, in the figure above I see that codon 217 (branch 34) transitioned from V to A and I can mark that on a tree. I can find all instances where amino acid 217 shows up in this output and see how the amino acids evolved along the branches.

### Thank you!
For the full analyis and results of this project see our [paper]().
