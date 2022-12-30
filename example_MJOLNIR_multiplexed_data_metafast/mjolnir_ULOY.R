# Load MJOLNIR silently
suppressPackageStartupMessages(library(mjolnir))

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c("ULO1_R1.fastq.gz","ULO2_R1.fastq.gz","ULO3_R1.fastq.gz","ULO4_R1.fastq.gz")

# Input identifiers for the individual libraries to be used. It should be a 4-character name, matching the information in the ngsfilter files
lib_prefixes <- c("ULO1","ULO2","ULO3","ULO4")

# Input name for the final combined library (should be a 4-character name)
lib <- "ULOY"


cores <- 3; setwd('~/Nextcloud/2_PROJECTES/MJOLNIR/adriantich_tests/');
obipath="/home/adriantich/obi3-env/bin/"; lib_prefixes <- c("ULO1"); lib <- "ULOY"

####################
# MJOLNIR pipeline #
####################

# Enter number of cores to be used in parallel for RAN and FREYJA. For best performance, the number of libraries to process x the number of cores should be less or equal than the total cores available in the system.
# e.g.: Here we will use 3 cores x 4 libraries = 12 cores. This can be run in any system with 12 cores or more.
cores <- 3

# RAN will distribute R1 & R2 fastq files into equal-sized pieces for parallel computing
mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1.",R2_motif="_R2.")

# FREYJA will do the paired-end alignment, demultiplexing & length filtering. It will give individual filtered sample files as an output.
mjolnir2_FREYJA(lib_prefix = lib_prefixes,lib = lib,cores = cores,Lmin=299,Lmax=320,obipath = obipath)

# Now you can enter the total number of cores available in the system, for full computing power.
# cores <- 16
cores <- 3

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib,cores)

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib,cores,d=13)

# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib,cores,tax_dir="~/taxo_NCBI",ref_db="DUFA_Leray_20200610.fasta",taxo_db="taxo_NCBI")

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)

# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib,min_id=.84)

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib,"ULOY_metadata.tsv")

