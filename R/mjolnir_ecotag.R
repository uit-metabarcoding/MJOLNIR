mjolnir_ecotag <- function(lib,cores,tax_dir,ref_db,taxo_db){
  message("Splitting seeds file in ",cores," fragments.")
  system(paste0("obidistribute -n ",cores," -p ",lib,".seeds ",lib,".seeds_nonsingleton.fasta"),intern=T,wait=T)
  message("Taxonomic assignment with ecotag")
  library(parallel)
  no_cores <- cores
  clust <- makeCluster(no_cores)
  X <- NULL
  for (i in 1:cores) X <- c(X,paste0("ecotag -d ",tax_dir,"/",taxo_db," -R ",tax_dir,"/",ref_db," ",lib,".seeds_",sprintf("%02d",i),".fasta > ",lib,".seeds.ecotag_",sprintf("%02d",i),".fasta"))
  clusterExport(clust, "X",envir = environment())
  parLapply(clust,X, function(x) system(x,intern=T,wait=T))
  stopCluster(clust)
  message("Adding higher taxonomy ranks.")
  system(paste0("cat ",lib,".seeds.ecotag_??.fasta > ",lib,".ecotag.fasta"),intern=T,wait=T)
  owi_add_taxonomy(paste0(lib,".ecotag.fasta"),tax_dir)
  message("Producing the combined file.")
  owi_combine(infile=paste0(lib,".ecotag.fasta.annotated.csv"),abundances=paste0(lib,".SWARM_output.counts.csv"),outfile=paste0(lib,".All_MOTUs.csv"))
}
