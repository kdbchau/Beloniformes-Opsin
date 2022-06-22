#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks=96
#SBATCH --time=05:00:00
#SBATCH --output=mpi_output_%j.txt
#SBATCH --mail-type=FAIL

cd $SLURM_SUBMIT_DIR

module load openmpi/2.2.1
module load gcc
module load samtools
module load bwa
module load gnu-parallel

REF=$1
GENE=$(basename $1)
path=$PWD

## INDEX REFERENCE FILES
if [ -f $REF.bwt ]
then
	echo "$REF already indexed, skipping"
else
	echo "Creating index for $REF"
	bwa index $REF
fi

## BWA MEM FOR EACH SPECIES IN TRIMMED FOLDER
for files in ../trimmed/*_trimmed.fastq
do
	if [ -f $(basename ${files%%_trimmed.fastq}_${GENE}_aln.sam) ]
	then
		echo "$(basename ${files%%_trimmed.fastq}_${GENE}_aln.sam) already exists, skipping."
	elif [ -f $(basename ${files%%_trimmed.fastq}_${GENE}_aln_sorted.bam) ]
	then
		echo "Sorted bam exists for ${files}, skipping."
	else
		echo "Running BWA MEM for $(basename ${files})"
		bwa mem -B 2 -M -t 20 $REF ${files} > ${files%%_trimmed.fastq}_${GENE}_aln.sam && mv ${files%%_trimmed.fastq}_${GENE}_aln.sam $path
		sleep 10
	fi
done


## CONVERT SAM TO BAM AND REMOVE SAM
for files2 in *_${GENE}_aln.sam
do
	if [ -f ${files2%%.sam}.bam ]
	then
		echo "${files2%%.sam}.bam already exists, skipping"
	else
		echo "Converting ${files2} to BAM and removing SAM file"
		samtools view -Sb -@ 48 ${files2} > ${files2%%.sam}.bam && echo "Removing ${files2}" | rm -f ${files2}
	fi
done

## SORT FILES AND REMOVE UNSORTED FILES
for files3 in *_${GENE}_aln.bam
do
	if [ -f ${files3%%.bam}_sorted.bam ]
	then
		echo "${files3} already has SORTED BAM, skipping."
	else
		echo "Sorting ${files3} and removing unsorted file"
		samtools sort -@ 48 ${files3} > ${files3%%.bam}_sorted.bam && echo "Removing ${files3}" | rm -f ${files3}
	fi
done

## INDEX SORTED FILES
for files4 in *${GENE}_aln_sorted.bam
do
	if [ -f ${files4}.bai ]
	then
		echo "Indexed ${files4} already exists, skipping."
	else
		echo "Indexing ${files4}"
		samtools index ${files4}
	fi
done

## MPILEUP
for files5 in *${GENE}_aln_sorted.bam
do
	if [ -f ${files5%%_sorted.bam}.mileup ]
	then
		echo "MPILEUP already exists for ${files5}, skipping."
	else
		echo "Creating mpileup for ${files5}"
		samtools mpileup -s -a -f $REF ${files5} -o ${files5%%_aln_sorted.bam}.mpileup
	fi
done

## RUN PYTHON SCRIPT TO CREATE UNAMBIGUOUS CONSENSUS FASTA FILES
for files in *${GENE}.mpileup
do
        if [ -f ${files%%.mpileup}_consensus.fa ]
        then
                echo "Conensus FASTA already exists for ${files}, skipping"
        else
                echo "Creating unambiguous consensus seq for ${files}"
                python ../consensus.py ${files}
        fi
done

## CONCATENATE FASTA FILES
for files2 in *${GENE}_consensus.fa
do
        extension="${GENE}*"
        NAME=${files2%%_$extension}
        echo "Concatenating the FASTA files"
        echo "$NAME"
        echo ">$NAME" >> ${GENE}_MSA.fa
        tail -n +2 ${files2} | tr '[:lower:]' '[:upper:]' | tr 'N/!/?' '-' >> ${GENE}_MSA.fa
        echo >> ${GENE}_MSA.fa
done

## REMOVE EXTRA FILES
for files in ./*{.mpileup,_consensus.fa,.fa.,.bam,.bai}
do
	echo "Removing consensus, mpileup, bam, bai, and indexed files"
	rm -f *.fa.* *.mpileup *.ba{i,m} *_consensus.fa
done
