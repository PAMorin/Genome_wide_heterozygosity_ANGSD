# Genome_wide_heterozygosity_ANGSD

ANGSD_SlidingWindow_heterozygosity
Phil Morin, Sept. 11, 2020

Scripts to detect variants in a genome assembly and count them across non-overlapping windows of specified size (e.g., 1MB).

These are adapted from methods originally described in (Westbury, M.V., Hartmann, S., Barlow, A., Wiesel, I., Leo, V., Welch, R., Parker, D.M., Sicks, F., Ludwig, A., Dalen, L., Hofreiter, M., 2018. Extended and continuous decline in effective population size results in low genomic diversity in the world's rarest hyena species, the brown hyena. Molecular Biology & Evolution 35, 1225-1237. doi:10.1093/molbev/msy037)

PROGRAM/version

###############

samtools/1.10

angsd/0.931 

bedtools/2.29.2 


SCRIPTS

###############

ROH_angsd_1MB_byChrom_SEDNA.sh

#########################################################################################

# Required files (prior to running sliding window script):

## Index your fasta and bam files

comment out if this step is already done

samtools faidx ${REFDIR}/${REF}

samtools index ${BAMDIR}/${BAMFILE}

## Get scaffold lengths for those that are over 1MB

do this once before running script, so that the scaffold lengths list is already in the OUTDIR

awk '$2 > 1000000 {print $1"\t"$2}' ${REFDIR}/${REF}.fai > ${OUTDIR}/${REF}_1MB_scaffold.lengths.txt

################################################################

## Get sliding windows of 1,000,000 bp

do this once before running script, so that the scaffold lengths list is already in the OUTDIR

bedtools makewindows -g ${OUTDIR}/${REF}_1MB_scaffold.lengths.txt -w 1000000 | awk '$3 ~ "000000" {print$1":"$2"-"$3}' > ${OUTDIR}/${REF}_1MB_windows.tx


# Variables

###############

THREADS=5

MINDEPTH=20   # 1/3x average coverage

MAXDEPTH=123	# 2x average coverage

MBQ=20  # minimum base quality filter
	

####################################################################

# Running script
Set the paths to files and scripts in the shell script ROH_angsd_1MB_byChrom_SEDNA.sh

Determine which scaffolds to scan for SNPs and set them as the array numbers (e.g., in the sbatch header: #SBATCH --array=1-20,22,24

##########################

# OUTPUT files

1) The .saf files ${TEMPDIR}/${OUTFILE}_${NUM}.saf. contains the site frequency spectrum data as binary files, stored in TEMPDIR, which also speeds the script if writing those files is on the scratch directory. 

2) The .sfs and .het files (stored in OUTDIR), contain the heterozygote frequency data. The .het files are used to plot heterozygosity across 1MB windows and summarize heterozygosity data.

##########################

# Summarizing and plotting data
R script "ROH_angsd_1MB_plot_100920.R". 

The R script loads the files based on filenames and filters based on the specified scaffold number (based on order in scaffolds list). 
OUTPUT files = 5 pdf files



