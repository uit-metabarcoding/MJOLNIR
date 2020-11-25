# FRIGGA: Final Recount and Integration of Generated Genealogies and Abundances
mjolnir_FRIGGA <- function(lib){
message("Producing the combined file.")
owi_combine(infile=paste0(lib,".ecotag.fasta.annotated.csv"),abundances=paste0(lib,".SWARM_output.counts.csv"),outfile=paste0(lib,".All_MOTUs.csv"))
}
