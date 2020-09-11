#!/bin/bash

# this is the header for SEDNA
#SBATCH --job-name=Het_angsd    ## job name
#SBATCH -e Het_angsd_%j.e.txt    ## error message name
#SBATCH -o Het_angsd.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#  SBATCH -c 20    ## <number of cores to ask for> (cpu cores per task?)
#SBATCH --mem=20G
#SBATCH -t 10:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 
#SBATCH --array=15-21
#SBATCH --ntasks=1                   ## Run a single task
#SBATCH -D /scratch/pmorin/temp		## <folder to change to when starting the job>


### Sliding window heterozygosity based on Angsd variant detection
##########################################################################################
# set parameters and folders
##########################################################################################
module load bio/samtools/1.10
module load bio/angsd/0.931 
module load bio/bedtools/2.29.2 

# directories
REFDIR=~/Ref_genomes/Psin/mPhoSin1.pri_ref_genome
BAMDIR=~/Ref_genomes/Psin/mPhoSin1_cur_10x_assembly_080620
OUTDIR=~/projects/Psin/ROH/angsd/test
TEMPDIR=/scratch/pmorin/temp

# files
BAMFILE=mPhoSin1_10x2cur_align_merged180620.bam
REF=GCF_008692025.1_mPhoSin1.pri_genomic.fna
OUTFILE=Psin_10x2mPhoSin1_pri
SCAFFOLDLIST=${REF}_1MB_scaffold.lengths.txt # see below for how to generate list
WINDOWSLIST=${REF}_1MB_windows.txt

# variables
THREADS=5
MINDEPTH=20   # 1/3x average coverage
MAXDEPTH=123	# 2x average coverage
MBQ=20  # minimum base quality filter
MAPQ=30  # minimum map quality filter
################################################################
# Index your fasta file
################################################################
# comment out if this step is already done
# samtools faidx ~/Ref_genomes/Psin/mPhoSin1.pri_ref_genome/GCF_008692025.1_mPhoSin1.pri_genomic.fna

# samtools index ${BAMDIR}/${BAMFILE}

################################################################
# Get scaffold lengths for those that are over 1MB
################################################################
# do this once before running script, so that the scaffold lengths list is already in the OUTDIR

# awk '$2 > 1000000 {print $1"\t"$2}' ${REFDIR}/${REF}.fai > ${OUTDIR}/${REF}_1MB_scaffold.lengths.txt

################################################################
# Get sliding windows of 1,000,000 bp
################################################################
# do this once before running script, so that the scaffold lengths list is already in the OUTDIR

# bedtools makewindows -g ${OUTDIR}/${REF}_1MB_scaffold.lengths.txt -w 1000000 | awk '$3 ~ "000000" {print$1":"$2"-"$3}' > ${OUTDIR}/${REF}_1MB_windows.txt

################################################################
# Get working scaffold based on array number
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})

# generate windows list for each scaffold
CHR=$(head -n ${NUM} ${OUTDIR}/${SCAFFOLDLIST} | tail -n 1)
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
CHR=$(echo ${CHR} | awk -F " " '{$0=$0":"}{print}')
echo ${CHR}

grep ${CHR} ${OUTDIR}/${WINDOWSLIST} > ${OUTDIR}/${REF}_${NUM}_windows.txt

################################################################
# Run angsd
################################################################
# .saf files are put into temp directory (scratch) to speed up processing and reduce clutter in output directory. See http://www.popgen.dk/angsd/index.php/Input#Genotype_Likelihood_Files for explanation of options.

angsd -GL 2 -setMinDepth ${MINDEPTH} -setMaxDepth ${MAXDEPTH} -minmapq ${MAPQ} -minq ${MBQ} -uniqueonly 1 -remove_bads 1 -only_proper_pairs 1 -docounts 1 -i ${BAMDIR}/${BAMFILE} -ref ${REFDIR}/${REF} -P ${THREADS} -out ${TEMPDIR}/${OUTFILE}_${NUM} -doSaf 1 -anc ${REFDIR}/${REF} -rf ${OUTDIR}/${REF}_${NUM}_windows.txt -baq 2 -fold 1

################################################################
# Run realSFS
################################################################
while read -r line; do realSFS -r $line ${TEMPDIR}/${OUTFILE}_${NUM}.saf.idx  -P ${THREADS} -tole 1e-8 2> log >> ${OUTDIR}/${OUTFILE}_1MB_combined_sfs_${NUM}.txt; done < ${OUTDIR}/${REF}_${NUM}_windows.txt
# -P = number of threads

################################################################
# Calculate heterozygosity (ends up being the value in column 4)
################################################################
awk '{print $1,$2,$3=$1+$2,$4=$2/$3}' ${OUTDIR}/${OUTFILE}_1MB_combined_sfs_${NUM}.txt >> ${OUTDIR}/${OUTFILE}_1MB_ROH_${NUM}_summary.het

########### Run in R script to plot and summarize outfile_ROH_summary.het #############