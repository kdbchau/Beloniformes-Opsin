# Multiple Ecological Axes Drive Cone Opsin Evolution in Beloniformes
This pipeline focuses on the extraction of cone opsin sequences from whole exome sequence data of fishes. It will also focus on multiple sequence alignment using both reference and makeshift reference sequences, cleaning up the multiple sequence alignments, and testing for changes in selection patterns with CODEML from the PAML program. Additional scripts and software are shown for phylogeny reconstruction and ancestral amino acid reconstruction.

This project was conducted at the [University of Toronto Scarborough](https://www.utsc.utoronto.ca/home/) in the [Lovejoy Lab](http://www.utsc.utoronto.ca/~lovejoy/).

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
    * [Additional Software](#11-additional-software)
2. [Quality Control](#2-quality-control)
3. [Read Mapping](#3-read-mapping)
    * [Round 1](#31-round-1)
    * [Round 2+](#32-round-2+)
4. [Phylogeny Reconstruction](#4-phylogeny-reconstruction)
    * [MrBayes](#41-mrbayes)
    * [IQ-TREE](#42-iq-tree)
5. [Ancestral Habitat and Diet Reconstruction](#6-ancestral-habitat-and-diet-reconstruction)
6. [Cleaning Multiple Sequence Alignments](#5-cleaning-multiple-sequence-alignments)



# 1. Install Software
Software with a * next to the name were already available in the Niagara cluster from Compute Canada. PAML was run on the Lovejoy server.

1. [__FastQC*__](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.8)
2. [__Trimmomatic*__](http://www.usadellab.org/cms/?page=trimmomatic) (version 0.38)
3. [__BWA*__](http://bio-bwa.sourceforge.net/) (version 0.7.17)
4. [__SAMtools*__](http://www.htslib.org/) (version 1.9)
5. [__IQTree__](http://www.iqtree.org/) (version 1.6.0)
6. [__MrBayes__](http://nbisweden.github.io/MrBayes/) (version 3.2.6)
7. [__BEAST__](https://beast.community/) (version 1.8.4)
8. [__PAML__](http://abacus.gene.ucl.ac.uk/software/paml.html) (version 4)
9. [__FigTree__](http://tree.bio.ed.ac.uk/software/figtree/)

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

However, upon trial and error it was found that because of high sequence similarity in the LWS opsins (see this [paper](https://www.nature.com/articles/s41598-019-39978-6)), and in the RH2 opsins for the medakas, I used zebrafish (_Danio rerio_) LWS and guppy RH2 (_Poecilia reticulata_) opsins in the first round of read mapping.

```
# Ran the script as follows on the cluster:

# Sbatch to set up the script on the job queue, then gave it a job name "lwsa", then called the reference sequence (a fasta file with just the freshwater medaka opsin sequence), and an output filename.
sbatch --job-name="lwsa" mapping.sh sws2a_medaka.fa sws2a_medaka_output
```
The above script will use a python script called [consensus.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/consensus.py) which will then take all the mpileup files generated after the mapping, and create a multiple sequence alignment (MSA) which we can alter edit.


## 3.2. Round 2+
Sometimes, using just freshwater medaka as a reference is not good enough. Because medaka are quite divergent from the rest of the beloniformes, we can use information ffrom the first round mapping to create "makeshift" or "chimeric" reference sequences that should be helpful to fill in more gaps. See the figure below.

![](https://github.com/kdbchau/Beloniformes/blob/main/Images/Screenshot%202022-06-23%20151643.png)

The second and subsequent rounds use these makeshift sequences that are filled in with nucleotides from the top most complete beloniform (after medaka which would be 100% complete because it is a direct match to itself).

To identify the next top beloniform, we used the [ideal_species.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/ideal_species.py) script. But this might not come out sorted properly so we run it the following way in terminal:

```
# name is the name of the file you are running without its extension, for example if it was rh2a_medaka_MSA.fa
name="rh2a_medaka_MSA" && python ideal_species.py ${name}.fa > ${name}.txt && cat ${name}.txt | sort -nr -k2 -o ${name}.txt
```

A textfile will be created with the species in the first column, sorted by the number of codons in the second column. Medakas are usually always the first two, followed by the beloniform with the next most complete sequence. Then using that beloniform's sequence, you can then fill in its gaps with either the medaka sequence or a different nucleotide/amino acid from another beloniform. Typically this is the best approach, as refilling in the gaps with medaka bases doesn't pull out more information. Rather, fill it with bases from the seccond most complete beloniform, or in some cases, if a paritcular family (e.g. all the flyingfishes) are missing codons in an area, fill the gaps with their bases to try and pull out more of their sequences in the next round mapping.

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

# 4. Phylogeny Reconstruction

Using the full MSAs (prior to cleaning), we can construct our phylogenys. I used IQTREE for a maximum likelihood reconstruction and MrBayes for a Bayesion reconstruction.

In both cases, all cone opsin MSAs were combined into one file and aligned as translated amino acids using MUSCLE alignemnt (use any program for this - I used [AliView](https://ormbunkar.se/aliview/)), and include some outgroup species for better root determination. I included the zebrafish and guppy opsin reference sequences I used for round 1 mapping.

You can also use programs like [AliView](https://ormbunkar.se/aliview/)) to convert the aligned MSA into Phylip foramt, which is needed for MrBayes and IQ-TREE (but IQ-TREE can use a fasta format as well).

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
In my case, it did have an upward trend but my [standard deviation values approached 0](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes_StandardDeviations.png) and my [PSRF values were around 1](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes_PSRF.png) - this is okay although this can be improved using a different burn-in, more runs, etc. The tree ends up looking the same as IQ-TREE so not a huge problem.

## 4.2. IQ-TREE

Ultrafast bootstrap IQ-TREE is simpler, once installed it can be called from a terminal.

Simply run it as:
```
iqtree -s combined_cone_opsins.fa -bb 1000  # -bb 1000 to get node values
```

The output ```.contree``` will have the node bootstrap values.

Both the Bayesian and maximum likelihood tree can be visualized using any tree visualizing program like FigTree. Here are the figures for the [Beloniformes Bayesian Cone Opsins Tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/MrBayes_AllConeOpsins.nexus.con.tre.pdf) and the [Beloniformes Maximum Likelihood Cone Opsins Tree](https://github.com/kdbchau/Beloniformes/blob/main/Images/IQTREE_AllOpsinsCombined_noRH1.fa.contree.pdf).

# 5. Ancestral Habitat and Diet Reconstruction

# 6. Cleaning Multiple Sequence Alignments
Premature stop codons or highly gapped regions will hinder calculations in codeml from the PAML program. It is advised to remove all stop codons (the very end of a protein-coding sequence; TAA, TAG, or TGA codons) or you can manually delete positions that are highly gapped.

To be as accurate as possible and keep things consistent between the 8 cone opsin MSAs, I use [alignment_editor.py](https://github.com/kdbchau/Beloniformes/blob/main/Scripts/alignment_editor.py) which will do three things:
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

