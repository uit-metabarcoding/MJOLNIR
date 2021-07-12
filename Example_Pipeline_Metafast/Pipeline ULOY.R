suppressPackageStartupMessages(library(mjolnir))

# Define input fastq files (only names of R1 files are needed)
R1_filenames <-c("ULO1_R1.fastq","ULO2_R1.fastq","ULO3_R1.fastq","ULO4_R1.fastq")

# Define number of cores to be used in parallel. For best performance, the number of libraries to process x the number of cores should be less or equal 
# than the total cores available in the computing server. e.g.: here we will use 30 cores x 4 libraries = 120 cores
cores <- 120

# Input names of the individual libraries to be used. It should be a 4-character name, matching the information of the ngsfilter files
lib_prefixes <- c("ULO1","ULO2","ULO3","ULO4")

# Input name for the final combined library (should be a 4-character name)
lib <- "ULOY"

####################
# MJOLNIR pipeline #
####################

# RAN will distribute R1 & R2 fastq files into equal-sized pieces for parallel computing
mjolnir1_RAN(R1_filenames,cores,lib_prefixes,R1_motif="_R1.",R2_motif="_R2.") 

# FREYJA will do the paired-end alignment, demultiplexing & length filtering
mjolnir2_FREYJA(lib_prefixes,cores,Lmin=299,Lmax=320)

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib_prefixes,lib,cores)

# Here we change the number of cores to full computing power. All libraries will be treated as a combined one from here on
cores <- 120

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib,cores,d=13)

# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib,cores,tax_dir="~/taxo_NCBI",ref_db="DUFA_Leray_20200610.fasta",taxo_db="taxo_NCBI") 

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)

# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib,min_id=.84)

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib,"ULOY_metadata.csv")
