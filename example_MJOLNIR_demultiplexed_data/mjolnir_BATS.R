# Load MJOLNIR silently
suppressPackageStartupMessages(library(mjolnir))

# Define number of cores to be used in parallel.
cores <- 16

# Input name for the final combined library (should be a 4-character name)
lib <- "BATS"

####################
# MJOLNIR pipeline #
####################

# FREYJA will do the paired-end alignment, demultiplexing & length filtering. This will be done for each sample file separately.
mjolnir2_FREYJA("", cores, Lmin=299, Lmax=320,lib,demultiplexed=T,R1_motif="_R1.",R2_motif="_R2.")

# HELA will remove chimaeric sequences in a sample-by-sample basis, will change identifiers of remaining unique sequences & will generate a table of their abundances in each sample & a fasta file with unique sequences and their total abundance for ODIN
mjolnir3_HELA(lib, cores)

# ODIN will do the clustering & will generate a table with the abundances of each MOTU in each sample
mjolnir4_ODIN(lib, cores, d=13, generate_ESV=FALSE)

# THOR will asign the taxonomy to the representative sequence of each MOTU
mjolnir5_THOR(lib, cores, tax_dir="~/taxo_NCBI", ref_db="DUFA_COLR_20210913.fasta", taxo_db="taxo_NCBI_20210720", run_ecotag=T)

# FRIGGA will integrate the information of MOTU abundances and taxonomy assignment from ODIN & THOR in a single table
mjolnir6_FRIGGA(lib)

# LOKI kill remove the pseudogenes and will keep track of the taxonomic information of the removed MOTUs
mjolnir7_LOKI(lib, min_id=.84)

# RAGNAROC will change the names of the samples to recover the original names and will remove unnecessary columns
mjolnir8_RAGNAROC(lib, "BATS_metadata.csv", "BATS_final_dataset.csv", sort_MOTUs="taxonomy", remove_bacteria=T, remove_contamination=F, min_reads=3)
